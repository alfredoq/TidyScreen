import os
import csv
import sqlite3
import pandas as pd
import re
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import numpy as np
import time

# Try to import tqdm for progress bars
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# Add the parent directory to path to import our local tidyscreen module
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir)

# Import from our local tidyscreen module
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject

def _process_chunk_worker(chunk_df: pd.DataFrame, smiles_column: str, 
                         name_column: Optional[str], flag_column: Optional[str],
                         name_available: bool, flag_available: bool, 
                         start_idx: int) -> List[Tuple[str, str, str]]:
    """
    Worker function to process a chunk of DataFrame in parallel.
    This function must be at module level to be pickleable for multiprocessing.
    
    Args:
        chunk_df (pd.DataFrame): Chunk of DataFrame to process
        smiles_column (str): Name of the SMILES column
        name_column (Optional[str]): Name of the name column
        flag_column (Optional[str]): Name of the flag column
        name_available (bool): Whether name column is available
        flag_available (bool): Whether flag column is available
        start_idx (int): Starting index for compound naming
        
    Returns:
        List[Tuple[str, str, str]]: List of tuples (smiles, name, flag)
    """
    compounds_data = []
    
    for local_idx, (_, row) in enumerate(chunk_df.iterrows()):
        global_idx = start_idx + local_idx
        smiles = str(row[smiles_column]).strip()
        
        # Handle name column - use provided column or generate from index
        if name_available:
            name = str(row[name_column]).strip()
            if not name or name.lower() in ['nan', 'none', '']:
                name = f"compound_{global_idx + 1}"  # Generate name if empty
        else:
            name = f"compound_{global_idx + 1}"  # Generate name if column doesn't exist
        
        # Handle flag column - use provided column or default to "nd"
        if flag_available:
            flag = str(row[flag_column]).strip()
            if not flag or flag.lower() in ['nan', 'none', '']:
                flag = "nd"  # Default flag if empty
        else:
            flag = "nd"  # Default flag if column doesn't exist
        
        # Skip empty SMILES rows (only SMILES is mandatory)
        if not smiles or smiles.lower() in ['nan', 'none', '']:
            continue
        
        compounds_data.append((smiles, name, flag))
    
    return compounds_data

def _compute_inchi_keys_worker(chunk_data: List[Tuple[int, str, str]]) -> List[Tuple[int, str, str]]:
    """
    Worker function to compute InChI keys for a chunk of SMILES strings in parallel.
    This function must be at module level to be pickleable for multiprocessing.
    
    Args:
        chunk_data (List[Tuple[int, str, str]]): List of tuples (row_index, smiles, compound_name)
        
    Returns:
        List[Tuple[int, str, str]]: List of tuples (row_index, inchi_key, compound_name)
    """
    try:
        # Import RDKit inside worker to avoid import issues
        from rdkit import Chem
        
        # Suppress RDKit warnings for cleaner parallel output
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        
    except ImportError:
        # If RDKit is not available, return error for all compounds in chunk
        return [(idx, 'ERROR', name) for idx, smiles, name in chunk_data]
    
    results = []
    
    for idx, smiles, name in chunk_data:
        try:
            # Skip empty or None SMILES
            if not smiles or smiles.strip() == '':
                results.append((idx, 'INVALID_SMILES', name))
                continue
            
            # Parse SMILES with RDKit
            mol = Chem.MolFromSmiles(str(smiles).strip())
            
            if mol is not None:
                # Compute InChI key
                inchi_key = Chem.MolToInchiKey(mol)
                if inchi_key:  # Ensure InChI key is not empty
                    results.append((idx, inchi_key, name))
                else:
                    results.append((idx, 'ERROR', name))
            else:
                # Invalid SMILES
                results.append((idx, 'INVALID_SMILES', name))
                
        except Exception:
            # Error computing InChI key
            results.append((idx, 'ERROR', name))
    
    return results

class ChemSpace:
    """
    ChemSpace class for managing chemical compound data within a project.
    Uses an ActivateProject object to access project information and database functionality.
    
    Handles reading CSV files containing SMILES strings, compound names, and flags,
    and stores them in a dedicated SQLite database within the project directory.
    """
    
    def __init__(self, project_obj: ActivateProject):
        """
        Initialize ChemSpace with an ActivateProject object.
        
        Args:
            project_obj (ActivateProject): An instantiated ActivateProject object
        """
        # Validate that we received a proper ActivateProject object
        if not isinstance(project_obj, ActivateProject):
            raise TypeError("ChemSpace requires an ActivateProject object")
        
        # Check if project exists
        if not project_obj.project_exists():
            raise ValueError(f"Project '{project_obj.name}' not found. Please create the project first.")
        
        # Store the project object and its attributes
        self.project = project_obj
        self.name = project_obj.name
        self.path = project_obj.path
        self.description = getattr(project_obj, 'description', None)
        self.id = getattr(project_obj, 'id', None)
        self.created_date = getattr(project_obj, 'created_date', None)
        
        # Set up chemspace database path within the project directory
        self.__chemspace_db = os.path.join(self.path, 'chemspace/processed_data', 'chemspace.db')
        
        # Ensure data directory exists
        data_dir = os.path.dirname(self.__chemspace_db)
        os.makedirs(data_dir, exist_ok=True)
        
        # Initialize the chemspace database
        self._initialize_chemspace_database()
        
        # Track loaded compounds
        self._compounds_loaded = False
        self._total_compounds = 0
    
    def _initialize_chemspace_database(self) -> bool:
        """
        Initialize the chemspace SQLite database file.
        
        Note: Tables are created dynamically based on CSV filenames during load_csv_file operations.
        This method only ensures the database file can be created/accessed.
        
        Returns:
            bool: True if database was initialized successfully
        """
        try:
            # Simply connect to create the database file if it doesn't exist
            conn = sqlite3.connect(self.__chemspace_db)
            conn.close()
            
            print(f"✅ ChemSpace database initialized: {self.__chemspace_db}")
            return True
            
        except Exception as e:
            print(f"❌ Error initializing ChemSpace database: {e}")
            return False
    
    def _sanitize_table_name(self, table_name: str) -> str:
        """
        Sanitize table name for SQLite compatibility.
        
        Args:
            table_name (str): Raw table name from CSV filename
            
        Returns:
            str: Sanitized table name safe for SQLite
        """
        import re
        
        # Remove special characters and replace with underscores
        sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', table_name)
        
        # Ensure it starts with a letter or underscore
        if sanitized and not sanitized[0].isalpha() and sanitized[0] != '_':
            sanitized = 't_' + sanitized
        
        # Ensure it's not empty
        if not sanitized:
            sanitized = 'compounds_table'
        
        # Limit length to reasonable size
        if len(sanitized) > 50:
            sanitized = sanitized[:50]
        
        return sanitized
    
    def _create_compounds_table(self, table_name: str) -> bool:
        """
        Create a compounds table with the specified name.
        
        Args:
            table_name (str): Name of the table to create
            
        Returns:
            bool: True if table was created successfully
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create compounds table
            create_table_query = f"""
            CREATE TABLE IF NOT EXISTS {table_name} (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT NOT NULL,
                name TEXT NOT NULL,
                flag TEXT,
                inchi_key TEXT,
                UNIQUE(name, smiles)
            )
            """
            
            cursor.execute(create_table_query)
            
            # Create indexes for faster searches
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_name ON {table_name}(name)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_smiles ON {table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_flag ON {table_name}(flag)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_inchi_key ON {table_name}(inchi_key)")
            
            conn.commit()
            conn.close()
            
            return True
            
        except Exception as e:
            print(f"❌ Error creating table '{table_name}': {e}")
            return False
    
    def drop_table(self, table_name: str, confirm: bool = True) -> bool:
        """
        Drop a specific table from the database.
        
        Args:
            table_name (str): Name of the table to drop
            confirm (bool): Whether to ask for confirmation
            
        Returns:
            bool: True if table was dropped successfully
        """
        try:
            # Check if table exists
            tables = self.get_all_tables()
            if table_name not in tables:
                print(f"⚠️  Table '{table_name}' does not exist")
                return False
            
            if confirm:
                response = input(f"⚠️  Are you sure you want to drop table '{table_name}'? This cannot be undone! (yes/no): ").strip().lower()
                if response not in ['yes', 'y']:
                    print("Table drop operation cancelled.")
                    return False
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Drop the table
            cursor.execute(f"DROP TABLE {table_name}")
            
            # Also drop associated indexes
            try:
                cursor.execute(f"DROP INDEX IF EXISTS idx_{table_name}_name")
                cursor.execute(f"DROP INDEX IF EXISTS idx_{table_name}_smiles")
                cursor.execute(f"DROP INDEX IF EXISTS idx_{table_name}_flag")
                cursor.execute(f"DROP INDEX IF EXISTS idx_{table_name}_inchi_key")
            except:
                pass  # Indexes might not exist
            
            conn.commit()
            conn.close()
            
            print(f"✅ Table '{table_name}' dropped successfully")
            return True
            
        except Exception as e:
            print(f"❌ Error dropping table '{table_name}': {e}")
            return False
    
    def get_all_tables(self) -> List[str]:
        """
        Get list of all compound tables in the database.
        
        Returns:
            List[str]: List of table names
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")
            tables = [row[0] for row in cursor.fetchall()]
            
            conn.close()
            return tables
            
        except Exception as e:
            print(f"❌ Error retrieving table list: {e}")
            return []
    
    def rename_table(self, old_name: str, new_name: str) -> bool:
        """
        Rename a table in the database.
        
        Args:
            old_name (str): Current name of the table
            new_name (str): New name for the table
            
        Returns:
            bool: True if table was renamed successfully
        """
        try:
            # Check if old table exists
            tables = self.get_all_tables()
            if old_name not in tables:
                print(f"⚠️  Table '{old_name}' does not exist")
                return False
            
            # Sanitize new table name
            new_name = self._sanitize_table_name(new_name)
            
            # Check if new table name already exists
            if new_name in tables:
                print(f"⚠️  Table '{new_name}' already exists")
                return False
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # SQLite doesn't support RENAME TABLE directly, so we need to:
            # 1. Create new table with same structure
            # 2. Copy data
            # 3. Drop old table
            
            # Create new table
            cursor.execute(f"""
            CREATE TABLE {new_name} AS 
            SELECT * FROM {old_name}
            """)
            
            # Recreate indexes for new table
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_name}_name ON {new_name}(name)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_name}_smiles ON {new_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_name}_flag ON {new_name}(flag)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_name}_inchi_key ON {new_name}(inchi_key)")
            
            # Drop old table and its indexes
            cursor.execute(f"DROP INDEX IF EXISTS idx_{old_name}_name")
            cursor.execute(f"DROP INDEX IF EXISTS idx_{old_name}_smiles")
            cursor.execute(f"DROP INDEX IF EXISTS idx_{old_name}_flag")
            cursor.execute(f"DROP INDEX IF EXISTS idx_{old_name}_inchi_key")
            cursor.execute(f"DROP TABLE {old_name}")
            
            conn.commit()
            conn.close()
            
            print(f"✅ Table '{old_name}' renamed to '{new_name}' successfully")
            return True
            
        except Exception as e:
            print(f"❌ Error renaming table '{old_name}' to '{new_name}': {e}")
            return False
    
    def load_csv_file(self, csv_file_path: str, 
                     smiles_column: str = 'smiles', 
                     name_column: Optional[str] = 'name', 
                     flag_column: Optional[str] = 'flag',
                     skip_duplicates: bool = True,
                     compute_inchi: bool = True,
                     parallel_threshold: int = 10000,
                     max_workers: Optional[int] = None,
                     chunk_size: Optional[int] = None) -> Dict[str, Any]:
        """
        Load compounds from a CSV file into the chemspace database.
        Creates a table named after the CSV file prefix.
        Automatically computes InChI keys for loaded SMILES.
        Uses parallel processing for large files to improve performance.
        
        Args:
            csv_file_path (str): Path to the CSV file
            smiles_column (str): Name of the column containing SMILES strings
            name_column (Optional[str]): Name of the column containing compound names. If None or not found, fills with "nd"
            flag_column (Optional[str]): Name of the column containing flags. If None or not found, fills with "nd"
            skip_duplicates (bool): Whether to skip duplicate compounds
            compute_inchi (bool): Whether to automatically compute InChI keys after loading
            parallel_threshold (int): Number of rows above which parallel processing is used (default: 10000)
            max_workers (Optional[int]): Maximum number of parallel workers. If None, uses cpu_count()
            chunk_size (Optional[int]): Size of each chunk for parallel processing. If None, automatically calculated
            
        Returns:
            dict: Results containing success status, counts, and messages
        """
        try:
            # Validate file exists
            if not os.path.exists(csv_file_path):
                return {
                    'success': False,
                    'message': f"CSV file not found: {csv_file_path}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1
                }
            
            # Extract table name from CSV file prefix
            csv_filename = os.path.basename(csv_file_path)
            table_name = os.path.splitext(csv_filename)[0]  # Remove .csv extension
            
            # Sanitize table name for SQLite (remove special characters, ensure it starts with letter)
            table_name = self._sanitize_table_name(table_name)
            
            # Create table for this CSV file
            if not self._create_compounds_table(table_name):
                return {
                    'success': False,
                    'message': f"Failed to create table '{table_name}' for CSV file",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1
                }
            
            # Read CSV file
            print(f"📖 Reading CSV file: {csv_file_path}")
            print(f"📋 Target table: {table_name}")
            df = pd.read_csv(csv_file_path)
            
            # Validate required columns (only SMILES is mandatory)
            if smiles_column not in df.columns:
                return {
                    'success': False,
                    'message': f"Missing required SMILES column: {smiles_column}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1
                }
            
            # Check if optional columns exist and warn if they don't
            name_available = name_column and name_column in df.columns
            flag_available = flag_column and flag_column in df.columns
            
            if not name_available:
                print(f"⚠️  Name column '{name_column}' not found. Will use 'nd' as default.")
            if not flag_available:
                print(f"⚠️  Flag column '{flag_column}' not found. Will use 'nd' as default.")
            
            # Determine if parallel processing should be used
            use_parallel = len(df) > parallel_threshold
            
            if use_parallel:
                print(f"🚀 Large dataset detected ({len(df)} rows). Using parallel processing...")
                
                # Set up parallel processing parameters
                if max_workers is None:
                    max_workers = min(cpu_count(), 8)  # Limit to reasonable number
                
                if chunk_size is None:
                    chunk_size = max(1000, len(df) // (max_workers * 4))  # Ensure reasonable chunk size
                
                print(f"   👥 Workers: {max_workers}")
                print(f"   📦 Chunk size: {chunk_size}")
                
                # Process in parallel chunks
                result = self._process_csv_parallel(
                    df, table_name, smiles_column, name_column, flag_column,
                    name_available, flag_available, skip_duplicates,
                    max_workers, chunk_size
                )
            else:
                print(f"📝 Processing {len(df)} rows sequentially...")
                
                # Sequential processing for smaller files
                compounds_data = []
                
                for idx, row in df.iterrows():
                    smiles = str(row[smiles_column]).strip()
                    
                    # Handle name column - use provided column or generate from index
                    if name_available:
                        name = str(row[name_column]).strip()
                        if not name or name.lower() in ['nan', 'none', '']:
                            name = f"compound_{idx + 1}"  # Generate name if empty
                    else:
                        name = f"compound_{idx + 1}"  # Generate name if column doesn't exist
                    
                    # Handle flag column - use provided column or default to "nd"
                    if flag_available:
                        flag = str(row[flag_column]).strip()
                        if not flag or flag.lower() in ['nan', 'none', '']:
                            flag = "nd"  # Default flag if empty
                    else:
                        flag = "nd"  # Default flag if column doesn't exist
                    
                    # Skip empty SMILES rows (only SMILES is mandatory)
                    if not smiles or smiles.lower() in ['nan', 'none', '']:
                        continue
                    
                    compounds_data.append((smiles, name, flag))
                
                # Insert compounds into database
                result = self._insert_compounds(compounds_data, table_name, skip_duplicates)
            
            if result['success']:
                self._compounds_loaded = True
                table_count = self.get_compound_count(table_name=table_name)
                
                print(f"✅ Successfully loaded compounds from {csv_filename}")
                print(f"   📋 Table: {table_name}")
                print(f"   📊 Compounds added: {result['compounds_added']}")
                print(f"   🔄 Duplicates skipped: {result['duplicates_skipped']}")
                print(f"   ❌ Errors: {result['errors']}")
                print(f"   📈 Total compounds in table '{table_name}': {table_count}")
                
                # Automatically compute InChI keys if requested
                if compute_inchi and result['compounds_added'] > 0:
                    print(f"\n🧪 Computing InChI keys for loaded compounds...")
                    try:
                        inchi_df = self.compute_inchi_keys(table_name, update_database=True)
                        if not inchi_df.empty:
                            # Count successful InChI computations
                            valid_inchi_count = len(inchi_df[
                                (inchi_df['inchi_key'].notna()) & 
                                (~inchi_df['inchi_key'].isin(['INVALID_SMILES', 'ERROR']))
                            ])
                            result['inchi_keys_computed'] = valid_inchi_count
                            print(f"   🔬 InChI keys computed: {valid_inchi_count}")
                        else:
                            result['inchi_keys_computed'] = 0
                            print(f"   ⚠️  No InChI keys computed")
                    except Exception as inchi_error:
                        print(f"   ⚠️  Warning: InChI key computation failed: {inchi_error}")
                        result['inchi_keys_computed'] = 0
                else:
                    result['inchi_keys_computed'] = 0
            
            result['table_name'] = table_name
            return result
            
        except Exception as e:
            return {
                'success': False,
                'message': f"Error loading CSV file: {e}",
                'compounds_added': 0,
                'duplicates_skipped': 0,
                'errors': 1
            }
    
    def _insert_compounds(self, compounds_data: List[Tuple], table_name: str, skip_duplicates: bool = True) -> Dict[str, Any]:
        """
        Insert compounds data into the database.
        
        Args:
            compounds_data (List[Tuple]): List of tuples containing compound data
            table_name (str): Name of the table to insert into
            skip_duplicates (bool): Whether to skip duplicate compounds
            
        Returns:
            dict: Results containing counts and status
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            compounds_added = 0
            duplicates_skipped = 0
            errors = 0
            
            insert_query = f"""
            INSERT INTO {table_name} (smiles, name, flag)
            VALUES (?, ?, ?)
            """
            
            for compound_data in compounds_data:
                try:
                    cursor.execute(insert_query, compound_data)
                    compounds_added += 1
                except sqlite3.IntegrityError:
                    if skip_duplicates:
                        duplicates_skipped += 1
                    else:
                        errors += 1
                except Exception as e:
                    print(f"⚠️  Error inserting compound {compound_data[1]}: {e}")
                    errors += 1
            
            conn.commit()
            conn.close()
            
            return {
                'success': True,
                'message': f"Inserted {compounds_added} compounds successfully",
                'compounds_added': compounds_added,
                'duplicates_skipped': duplicates_skipped,
                'errors': errors
            }
            
        except Exception as e:
            return {
                'success': False,
                'message': f"Database error: {e}",
                'compounds_added': 0,
                'duplicates_skipped': 0,
                'errors': len(compounds_data)
            }
    
    def get_compounds(self, table_name: Optional[str] = None,
                     limit: Optional[int] = None, 
                     flag_filter: Optional[str] = None,
                     name_pattern: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Retrieve compounds from the database with optional filtering.
        
        Args:
            table_name (Optional[str]): Name of the table to query from. If None, uses 'compounds'
            limit (Optional[int]): Maximum number of compounds to return
            flag_filter (Optional[str]): Filter by specific flag value
            name_pattern (Optional[str]): Filter by name pattern (SQL LIKE)
            
        Returns:
            List[Dict]: List of compound dictionaries
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Use default table if none specified
            if table_name is None:
                table_name = 'compounds'
            
            # Build query with filters
            query = f"SELECT * FROM {table_name} WHERE 1=1"
            params = []
            
            if flag_filter:
                query += " AND flag = ?"
                params.append(flag_filter)
            
            if name_pattern:
                query += " AND name LIKE ?"
                params.append(f"%{name_pattern}%")
            
            query += " ORDER BY id ASC"
            
            if limit:
                query += " LIMIT ?"
                params.append(limit)
            
            cursor.execute(query, params)
            rows = cursor.fetchall()
            
            # Get column names
            column_names = [description[0] for description in cursor.description]
            
            # Convert to dictionaries
            compounds = []
            for row in rows:
                compound_dict = {}
                for i, column_name in enumerate(column_names):
                    compound_dict[column_name] = row[i]
                compounds.append(compound_dict)
            
            conn.close()
            return compounds
            
        except Exception as e:
            print(f"❌ Error retrieving compounds from table '{table_name}': {e}")
            return []
    
    def get_compound_count(self, table_name: Optional[str] = None, flag_filter: Optional[str] = None) -> int:
        """
        Get the total number of compounds in the database.
        
        Args:
            table_name (Optional[str]): Name of the table to count from. If None, uses 'compounds'
            flag_filter (Optional[str]): Filter by specific flag value
            
        Returns:
            int: Number of compounds
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Use default table if none specified
            if table_name is None:
                table_name = 'compounds'
            
            if flag_filter:
                cursor.execute(f"SELECT COUNT(*) FROM {table_name} WHERE flag = ?", (flag_filter,))
            else:
                cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            
            count = cursor.fetchone()[0]
            conn.close()
            return count
            
        except Exception as e:
            print(f"❌ Error counting compounds in table '{table_name}': {e}")
            return 0
    
    def export_to_csv(self, output_path: str, 
                     table_name: Optional[str] = None,
                     flag_filter: Optional[str] = None) -> bool:
        """
        Export compounds to a CSV file.
        
        Args:
            output_path (str): Path for the output CSV file
            table_name (Optional[str]): Name of the table to export from. If None, exports all tables
            flag_filter (Optional[str]): Filter by specific flag value
            
        Returns:
            bool: True if export was successful
        """
        try:
            compounds = self.get_compounds(table_name=table_name, flag_filter=flag_filter)
            
            if not compounds:
                table_desc = f" from table '{table_name}'" if table_name else " from database"
                print(f"⚠️  No compounds found to export{table_desc}")
                return False
            
            # Convert to DataFrame and save
            df = pd.DataFrame(compounds)
            df.to_csv(output_path, index=False)
            
            table_desc = f" from table '{table_name}'" if table_name else ""
            print(f"✅ Exported {len(compounds)} compounds{table_desc} to: {output_path}")
            return True
            
        except Exception as e:
            table_desc = f" from table '{table_name}'" if table_name else ""
            print(f"❌ Error exporting compounds{table_desc}: {e}")
            return False
    
    def filter_using_smarts(self, table_name: str, id: int):
        """
        Filter compounds in a table using a SMARTS pattern retrieved by ID from the projects database.
        
        Args:
            table_name (str): Name of the table to filter
            id (int): ID of the SMARTS filter in the projects database
        """
    
        try:
            # Check if the table exists
            tables = self.get_all_tables()
            if table_name not in tables:
                print(f"❌ Table '{table_name}' does not exist in chemspace database")
                return

            # Get the projects database path
            import site
            projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"
            
            if not os.path.exists(projects_db):
                print(f"❌ Projects database not found: {projects_db}")
                return
            
            # Connect to projects database and retrieve the SMARTS pattern by ID
            projects_conn = sqlite3.connect(projects_db)
            projects_cursor = projects_conn.cursor()
            
            # Query to get the SMARTS pattern and filter name by ID
            smarts_query = "SELECT filter_name, smarts FROM chem_filters WHERE id = ?"
            projects_cursor.execute(smarts_query, (id,))
            smarts_result = projects_cursor.fetchone()
            
            projects_conn.close()
            
            if not smarts_result:
                print(f"❌ No SMARTS filter found with ID {id} in projects database")
                return
            
            filter_name, smarts_pattern = smarts_result
            print(f"🔍 Found SMARTS filter: '{filter_name}' (ID: {id})")
            print(f"🧪 SMARTS pattern: {smarts_pattern}")

            # Connect to chemspace database and retrieve the table
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()

            # Query to retrieve all compounds from the specified table
            query = f"SELECT id, smiles, name, flag FROM {table_name}"
            cursor.execute(query)
            compounds = cursor.fetchall()
        
            conn.close()

            if not compounds:
                print(f"⚠️  No compounds found in table '{table_name}'")
                return
        
            print(f"📊 Retrieved {len(compounds)} compounds from table '{table_name}' for filtering")
            print(f"🔬 Starting SMARTS pattern matching...")
            
            # Import RDKit for SMARTS pattern matching
            try:
                from rdkit import Chem
                from rdkit import RDLogger
                # Suppress RDKit warnings for cleaner output
                RDLogger.DisableLog('rdApp.*')
            except ImportError:
                print("❌ RDKit not installed. Please install RDKit to use SMARTS filtering:")
                print("   conda install -c conda-forge rdkit")
                print("   or")
                print("   pip install rdkit")
                return
            
            # Parse the SMARTS pattern
            smarts_mol = Chem.MolFromSmarts(smarts_pattern)
            if smarts_mol is None:
                print(f"❌ Invalid SMARTS pattern: {smarts_pattern}")
                return
                
            # Initialize counters
            matching_compounds = []
            invalid_smiles_count = 0
            
            # Process compounds and apply SMARTS filter
            progress_bar = None
            if TQDM_AVAILABLE:
                progress_bar = tqdm(total=len(compounds), desc="Filtering compounds", unit="compound")
            
            for compound in compounds:
                compound_id, smiles, name, flag = compound
                
                # Parse SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    invalid_smiles_count += 1
                    if progress_bar:
                        progress_bar.update(1)
                    continue
                
                # Check if compound matches the SMARTS pattern
                if mol.HasSubstructMatch(smarts_mol):
                    matching_compounds.append({
                        'id': compound_id,
                        'smiles': smiles,
                        'name': name,
                        'flag': flag
                    })
                
                if progress_bar:
                    progress_bar.update(1)
            
            if progress_bar:
                progress_bar.close()
            
            # Print results summary
            print(f"✅ SMARTS filtering completed!")
            print(f"📋 Filter applied: '{filter_name}' (ID: {id})")
            print(f"🔍 SMARTS pattern: {smarts_pattern}")
            print(f"📊 Results:")
            print(f"   • Total compounds processed: {len(compounds)}")
            print(f"   • Matching compounds: {len(matching_compounds)}")
            print(f"   • Non-matching compounds: {len(compounds) - len(matching_compounds) - invalid_smiles_count}")
            print(f"   • Invalid SMILES: {invalid_smiles_count}")
            print(f"   • Match percentage: {len(matching_compounds)/len(compounds)*100:.1f}%")
            
            # Optionally, save filtered results to a new table or file
            if matching_compounds:
                print(f"\n💾 Found {len(matching_compounds)} compounds matching the filter.")
                save_option = input("Do you want to save the filtered compounds to a new table? (y/n): ").strip().lower()
                if save_option in ['y', 'yes']:
                    new_table_name = input(f"Enter new table name (default: {table_name}_filtered_{filter_name.lower().replace(' ', '_').replace('-', '_')}): ").strip()
                    if not new_table_name:
                        new_table_name = f"{table_name}_filtered_{filter_name.lower().replace(' ', '_').replace('-', '_')}"
                    
                    # Save matching compounds to new table
                    self._save_filtered_compounds(matching_compounds, new_table_name)
            else:
                print("⚠️  No compounds matched the filter criteria.")
            
            
        except Exception as e:
            print(f"❌ Error retrieving SMARTS pattern or compounds for filtering: {e}")
            return []   

    def _save_filtered_compounds(self, compounds: List[Dict], table_name: str) -> bool:
        """
        Save filtered compounds to a new table in the database.
        Removes duplicates based on SMILES strings before saving.
        
        Args:
            compounds (List[Dict]): List of compound dictionaries with keys: id, smiles, name, flag
            table_name (str): Name of the new table to create
            
        Returns:
            bool: True if compounds were saved successfully, False otherwise
        """
        try:
            # Remove duplicates based on SMILES strings
            unique_compounds = []
            seen_smiles = set()
            duplicates_count = 0
            
            print(f"📋 Checking for duplicates among {len(compounds)} filtered compounds...")
            
            for compound in compounds:
                smiles = compound['smiles']
                if smiles not in seen_smiles:
                    unique_compounds.append(compound)
                    seen_smiles.add(smiles)
                else:
                    duplicates_count += 1
            
            print(f"🔍 Found {duplicates_count} duplicate compounds (removed)")
            print(f"✨ Retaining {len(unique_compounds)} unique compounds for saving")
            
            if not unique_compounds:
                print("⚠️  No unique compounds to save after duplicate removal")
                return True
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create the new table with the same schema as other compound tables
            # Add UNIQUE constraint on SMILES to prevent future duplicates
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {table_name} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT NOT NULL UNIQUE,
                    name TEXT,
                    flag INTEGER DEFAULT 1
                )
            ''')
            
            # Create indexes for better performance
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_smiles ON {table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_name ON {table_name}(name)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_flag ON {table_name}(flag)")
            
            # Insert the unique filtered compounds
            inserted_count = 0
            database_duplicates = 0
            
            for compound in unique_compounds:
                try:
                    cursor.execute(f'''
                        INSERT INTO {table_name} (smiles, name, flag)
                        VALUES (?, ?, ?)
                    ''', (compound['smiles'], compound['name'], compound['flag']))
                    inserted_count += 1
                except sqlite3.IntegrityError:
                    # Handle case where compound already exists in database
                    database_duplicates += 1
            
            conn.commit()
            conn.close()
            
            print(f"✅ Successfully saved {inserted_count} unique filtered compounds to table '{table_name}'")
            if database_duplicates > 0:
                print(f"⚠️  Skipped {database_duplicates} compounds that were already in the database")
            
            return True
            
        except Exception as e:
            print(f"❌ Error saving filtered compounds to table '{table_name}': {e}")
            return False

    def clear_database(self, table_name: Optional[str] = None, confirm: bool = True) -> bool:
        """
        Clear compounds from the database.
        
        Args:
            table_name (Optional[str]): Name of the table to clear. If None, clears all tables
            confirm (bool): Whether to ask for confirmation
            
        Returns:
            bool: True if database was cleared successfully
        """
        if confirm:
            if table_name:
                response = input(f"⚠️  Are you sure you want to clear all compounds from table '{table_name}'? This cannot be undone! (yes/no): ").strip().lower()
            else:
                response = input("⚠️  Are you sure you want to clear all compounds from all tables? This cannot be undone! (yes/no): ").strip().lower()
            
            if response not in ['yes', 'y']:
                print("Database clear operation cancelled.")
                return False
        
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            if table_name:
                # Clear specific table
                cursor.execute(f"DELETE FROM {table_name}")
                print(f"✅ All compounds cleared from table '{table_name}'")
            else:
                # Clear all tables
                tables = self.get_all_tables()
                for table in tables:
                    cursor.execute(f"DELETE FROM {table}")
                print(f"✅ All compounds cleared from {len(tables)} tables")
            
            conn.commit()
            conn.close()
            
            self._compounds_loaded = False
            self._total_compounds = 0
            
            return True
            
        except Exception as e:
            table_desc = f" from table '{table_name}'" if table_name else ""
            print(f"❌ Error clearing compounds{table_desc}: {e}")
            return False
    
    def _get_table_as_dataframe(self, table_name) -> pd.DataFrame:
        """
        Return a table from the database as a pandas DataFrame.
        
        Args:
            table_name: Name of the table to retrieve
            
        Returns:
            pd.DataFrame: DataFrame containing the table data, empty DataFrame if error occurs
        """
        try:
            # Check if table exists
            if table_name is not None:
                tables = self.get_all_tables()
                if table_name not in tables:
                    print(f"⚠️  Table '{table_name}' does not exist")
                    return pd.DataFrame()
            
            conn = sqlite3.connect(self.__chemspace_db)
            
            # Build query with filters
            query = f"SELECT * FROM {table_name}"
            
            df = pd.read_sql_query(query, conn)
            
            conn.close()
            
            print(f"✅ Retrieved table '{table_name}' as DataFrame with {len(df)} rows")
            return df
            
        except Exception as e:
            print(f"❌ Error retrieving table '{table_name}' as DataFrame: {e}")
            return pd.DataFrame()
    
    def compute_inchi_keys(self, table_name: str, 
                          update_database: bool = True,
                          parallel_threshold: int = 1000,
                          max_workers: Optional[int] = None,
                          chunk_size: Optional[int] = None) -> pd.DataFrame:
        """
        Retrieve a table as DataFrame and compute InChI keys for each SMILES string using parallel processing.
        
        Args:
            table_name (str): Name of the table to retrieve
            update_database (bool): Whether to add inchi_key column to the database table
            parallel_threshold (int): Number of compounds above which parallel processing is used (default: 1000)
            max_workers (Optional[int]): Maximum number of parallel workers. If None, uses cpu_count()
            chunk_size (Optional[int]): Size of each chunk for parallel processing. If None, automatically calculated
            
        Returns:
            pd.DataFrame: DataFrame with added 'inchi_key' column, empty DataFrame if error occurs
        """
        try:
            # Check RDKit availability early
            try:
                from rdkit import Chem
            except ImportError:
                print("❌ RDKit not installed. Please install RDKit to compute InChI keys:")
                print("   conda install -c conda-forge rdkit")
                print("   or")
                print("   pip install rdkit")
                return pd.DataFrame()
            
            # Get the DataFrame using existing method
            print(f"📊 Retrieving table '{table_name}' for InChI key computation...")
            df = self._get_table_as_dataframe(table_name)
            
            if df.empty:
                print(f"⚠️  No data retrieved from table '{table_name}'")
                return pd.DataFrame()
            
            # Check if SMILES column exists
            if 'smiles' not in df.columns:
                print("❌ No 'smiles' column found in the DataFrame")
                return pd.DataFrame()
            
            # Initialize InChI key column
            df['inchi_key'] = None
            
            num_compounds = len(df)
            print(f"🔬 Computing InChI keys for {num_compounds} compounds...")
            
            # Determine processing method
            use_parallel = num_compounds > parallel_threshold
            
            if use_parallel:
                print(f"🚀 Large dataset detected ({num_compounds} compounds). Using parallel processing...")
                
                # Set up parallel processing parameters
                if max_workers is None:
                    max_workers = min(cpu_count(), 8)  # Limit to reasonable number
                
                if chunk_size is None:
                    chunk_size = max(100, num_compounds // (max_workers * 4))  # Ensure reasonable chunk size
                
                print(f"   👥 Workers: {max_workers}")
                print(f"   📦 Chunk size: {chunk_size}")
                
                # Process in parallel
                successful_computations, failed_computations = self._compute_inchi_keys_parallel(
                    df, max_workers, chunk_size
                )
            else:
                print(f"📝 Processing {num_compounds} compounds sequentially...")
                # Sequential processing for smaller datasets
                successful_computations, failed_computations = self._compute_inchi_keys_sequential(df)
            
            # Print summary
            print(f"✅ InChI key computation completed:")
            print(f"   🎯 Successful: {successful_computations}")
            print(f"   ❌ Failed: {failed_computations}")
            
            # Optionally update the database table
            if update_database and successful_computations > 0:
                print(f"💾 Updating database table '{table_name}' with InChI keys...")
                success = self._update_table_with_inchi_keys(table_name, df)
                if success:
                    print(f"✅ Database table '{table_name}' updated with InChI keys")
                else:
                    print(f"❌ Failed to update database table '{table_name}'")
            
            return df
            
        except Exception as e:
            print(f"❌ Error computing InChI keys for table '{table_name}': {e}")
            return pd.DataFrame()
    
    def _compute_inchi_keys_parallel(self, df: pd.DataFrame, max_workers: int, chunk_size: int) -> Tuple[int, int]:
        """
        Compute InChI keys using parallel processing.
        
        Args:
            df (pd.DataFrame): DataFrame containing SMILES data
            max_workers (int): Maximum number of parallel workers
            chunk_size (int): Size of each chunk for parallel processing
            
        Returns:
            Tuple[int, int]: (successful_computations, failed_computations)
        """
        try:
            start_time = time.time()
            
            # Prepare data for parallel processing
            print(f"   📦 Preparing data for parallel processing...")
            chunk_data_list = []
            
            # Split data into chunks
            for i in range(0, len(df), chunk_size):
                chunk_df = df.iloc[i:i + chunk_size]
                chunk_data = []
                
                for idx, row in chunk_df.iterrows():
                    smiles = row['smiles']
                    name = row.get('name', 'unknown')
                    chunk_data.append((idx, smiles, name))
                
                if chunk_data:  # Only add non-empty chunks
                    chunk_data_list.append(chunk_data)
            
            num_chunks = len(chunk_data_list)
            print(f"   📊 Split {len(df)} compounds into {num_chunks} chunks")
            
            # Initialize progress tracking
            successful_computations = 0
            failed_computations = 0
            processed_chunks = 0
            
            # Initialize progress bar if tqdm is available
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(df),
                    desc="Computing InChI keys",
                    unit="compounds",
                    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
                )
            else:
                progress_bar = None
                print(f"   🚀 Starting parallel processing of {num_chunks} chunks...")
            
            # Process chunks in parallel
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all chunk processing jobs
                future_to_chunk = {}
                chunk_submission_start = time.time()
                
                for i, chunk_data in enumerate(chunk_data_list):
                    future = executor.submit(_compute_inchi_keys_worker, chunk_data)
                    future_to_chunk[future] = i
                
                chunk_submission_time = time.time() - chunk_submission_start
                print(f"   ⏱️  All {num_chunks} chunks submitted in {chunk_submission_time:.2f}s")
                
                # Collect results
                failed_chunks = 0
                
                for future in as_completed(future_to_chunk):
                    chunk_idx = future_to_chunk[future]
                    
                    try:
                        chunk_results = future.result()
                        
                        # Update DataFrame with results
                        chunk_successful = 0
                        chunk_failed = 0
                        
                        for idx, inchi_key, name in chunk_results:
                            df.at[idx, 'inchi_key'] = inchi_key
                            
                            if inchi_key in ['INVALID_SMILES', 'ERROR']:
                                chunk_failed += 1
                                # Show only first few errors
                                if failed_computations + chunk_failed <= 5:
                                    error_type = "Invalid SMILES" if inchi_key == 'INVALID_SMILES' else "Computation error"
                                    print(f"⚠️  {error_type} for compound '{name}' (chunk {chunk_idx})")
                            else:
                                chunk_successful += 1
                        
                        successful_computations += chunk_successful
                        failed_computations += chunk_failed
                        processed_chunks += 1
                        
                        # Update progress
                        chunk_size_actual = len(chunk_results)
                        if progress_bar:
                            progress_bar.update(chunk_size_actual)
                        else:
                            # Print progress every 5 chunks if no tqdm
                            if processed_chunks % 5 == 0 or processed_chunks == num_chunks:
                                progress_pct = (processed_chunks / num_chunks) * 100
                                print(f"   📈 Progress: {processed_chunks}/{num_chunks} chunks ({progress_pct:.1f}%)")
                    
                    except Exception as e:
                        failed_chunks += 1
                        chunk_size_est = len(chunk_data_list[chunk_idx]) if chunk_idx < len(chunk_data_list) else chunk_size
                        failed_computations += chunk_size_est
                        processed_chunks += 1
                        
                        print(f"   ❌ Error processing chunk {chunk_idx}: {e}")
                        
                        # Update progress even for failed chunks
                        if progress_bar:
                            progress_bar.update(chunk_size_est)
            
            # Close progress bar
            if progress_bar:
                progress_bar.close()
            
            processing_time = time.time() - start_time
            print(f"   ⏱️  Parallel processing completed in {processing_time:.2f}s")
            
            if failed_chunks > 0:
                print(f"   ⚠️  {failed_chunks} chunks failed during processing")
            
            return successful_computations, failed_computations
            
        except Exception as e:
            print(f"❌ Error in parallel InChI key computation: {e}")
            return 0, len(df)
    
    def _compute_inchi_keys_sequential(self, df: pd.DataFrame) -> Tuple[int, int]:
        """
        Compute InChI keys using sequential processing.
        
        Args:
            df (pd.DataFrame): DataFrame containing SMILES data
            
        Returns:
            Tuple[int, int]: (successful_computations, failed_computations)
        """
        try:
            from rdkit import Chem
        except ImportError:
            print("❌ RDKit not available for sequential processing")
            return 0, len(df)
        
        successful_computations = 0
        failed_computations = 0
        
        # Initialize progress bar if tqdm is available
        if TQDM_AVAILABLE:
            progress_bar = tqdm(
                total=len(df),
                desc="Computing InChI keys",
                unit="compounds"
            )
        else:
            progress_bar = None
        
        for idx, row in df.iterrows():
            smiles = row['smiles']
            name = row.get('name', 'unknown')
            
            try:
                # Parse SMILES with RDKit
                mol = Chem.MolFromSmiles(smiles)
                
                if mol is not None:
                    # Compute InChI key
                    inchi_key = Chem.MolToInchiKey(mol)
                    df.at[idx, 'inchi_key'] = inchi_key
                    successful_computations += 1
                else:
                    # Invalid SMILES
                    df.at[idx, 'inchi_key'] = 'INVALID_SMILES'
                    failed_computations += 1
                    if failed_computations <= 5:  # Show only first 5 errors
                        print(f"⚠️  Invalid SMILES for compound '{name}': {smiles}")
            
            except Exception as e:
                df.at[idx, 'inchi_key'] = 'ERROR'
                failed_computations += 1
                if failed_computations <= 5:  # Show only first 5 errors
                    print(f"⚠️  Error computing InChI key for '{name}': {e}")
            
            # Update progress
            if progress_bar:
                progress_bar.update(1)
        
        # Close progress bar
        if progress_bar:
            progress_bar.close()
        
        return successful_computations, failed_computations

    def _update_table_with_inchi_keys(self, table_name: str, df: pd.DataFrame) -> bool:
        """
        Update the database table with computed InChI key values.
        Handles both new tables (with inchi_key column) and existing tables (adds column if needed).
        
        Args:
            table_name (str): Name of the table to update
            df (pd.DataFrame): DataFrame containing the computed InChI keys
            
        Returns:
            bool: True if update was successful
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Check if inchi_key column exists, add if it doesn't (for backward compatibility)
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [row[1] for row in cursor.fetchall()]
            
            if 'inchi_key' not in columns:
                try:
                    cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN inchi_key TEXT")
                    cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_inchi_key ON {table_name}(inchi_key)")
                    print(f"   📋 Added 'inchi_key' column to existing table '{table_name}'")
                except sqlite3.OperationalError as e:
                    print(f"   ⚠️  Warning: Could not add inchi_key column: {e}")
                    return False
            
            # Update each row with computed InChI key
            update_query = f"UPDATE {table_name} SET inchi_key = ? WHERE id = ?"
            
            updates_made = 0
            for _, row in df.iterrows():
                if pd.notna(row.get('inchi_key')) and row['inchi_key'] not in ['INVALID_SMILES', 'ERROR']:
                    cursor.execute(update_query, (row['inchi_key'], row['id']))
                    updates_made += 1
            
            conn.commit()
            conn.close()
            
            print(f"   📊 Updated {updates_made} rows with InChI keys")
            return True
            
        except Exception as e:
            print(f"❌ Error updating database table '{table_name}' with InChI keys: {e}")
            return False
    
    def _create_sql_registers_table(self) -> bool:
        """
        Create the sql_registers table to track SQL import operations.
        
        Returns:
            bool: True if table was created successfully
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create sql_registers table
            create_table_query = """
            CREATE TABLE IF NOT EXISTS sql_registers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sql_query TEXT NOT NULL,
                source_db_path TEXT NOT NULL,
                target_table_name TEXT NOT NULL,
                creation_date TEXT NOT NULL,
                compounds_added INTEGER DEFAULT 0,
                duplicates_skipped INTEGER DEFAULT 0,
                errors INTEGER DEFAULT 0,
                inchi_keys_computed INTEGER DEFAULT 0,
                success BOOLEAN DEFAULT 0
            )
            """
            
            cursor.execute(create_table_query)
            
            # Create indexes for faster searches
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_sql_registers_target_table ON sql_registers(target_table_name)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_sql_registers_creation_date ON sql_registers(creation_date)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_sql_registers_source_db ON sql_registers(source_db_path)")
            
            conn.commit()
            conn.close()
            
            return True
            
        except Exception as e:
            print(f"❌ Error creating sql_registers table: {e}")
            return False
    
    def _register_sql_operation(self, sql_query: str, source_db_path: str, 
                               target_table_name: str, result: Dict[str, Any]) -> bool:
        """
        Register an SQL import operation in the sql_registers table.
        
        Args:
            sql_query (str): The SQL query that was executed
            source_db_path (str): Path to the source database
            target_table_name (str): Name of the target table created
            result (Dict[str, Any]): Result dictionary from the import operation
            
        Returns:
            bool: True if registration was successful
        """
        try:
            # Ensure sql_registers table exists
            if not self._create_sql_registers_table():
                return False
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get current timestamp
            creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # Insert register entry
            insert_query = """
            INSERT INTO sql_registers 
            (sql_query, source_db_path, target_table_name, creation_date, 
             compounds_added, duplicates_skipped, errors, inchi_keys_computed, success)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            
            cursor.execute(insert_query, (
                sql_query,
                source_db_path,
                target_table_name,
                creation_date,
                result.get('compounds_added', 0),
                result.get('duplicates_skipped', 0),
                result.get('errors', 0),
                result.get('inchi_keys_computed', 0),
                result.get('success', False)
            ))
            
            conn.commit()
            conn.close()
            
            print(f"📝 SQL operation registered in sql_registers table")
            return True
            
        except Exception as e:
            print(f"⚠️  Warning: Failed to register SQL operation: {e}")
            return False
    
    def get_sql_registers(self, limit: Optional[int] = None, 
                         target_table_filter: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Retrieve SQL import operation registers from the database.
        
        Args:
            limit (Optional[int]): Maximum number of registers to return
            target_table_filter (Optional[str]): Filter by target table name
            
        Returns:
            List[Dict]: List of register dictionaries
        """
        try:
            # Ensure sql_registers table exists
            if not self._create_sql_registers_table():
                return []
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Build query with filters
            query = "SELECT * FROM sql_registers WHERE 1=1"
            params = []
            
            if target_table_filter:
                query += " AND target_table_name = ?"
                params.append(target_table_filter)
            
            query += " ORDER BY creation_date DESC"
            
            if limit:
                query += " LIMIT ?"
                params.append(limit)
            
            cursor.execute(query, params)
            rows = cursor.fetchall()
            
            # Get column names
            column_names = [description[0] for description in cursor.description]
            
            # Convert to dictionaries
            registers = []
            for row in rows:
                register_dict = {}
                for i, column_name in enumerate(column_names):
                    register_dict[column_name] = row[i]
                registers.append(register_dict)
            
            conn.close()
            return registers
            
        except Exception as e:
            print(f"❌ Error retrieving SQL registers: {e}")
            return []
    
    def print_sql_registers_summary(self, limit: int = 10) -> None:
        """
        Print a formatted summary of SQL import operation registers.
        
        Args:
            limit (int): Maximum number of recent registers to show
        """
        registers = self.get_sql_registers(limit=limit)
        
        if not registers:
            print("📝 No SQL import operations registered yet")
            return
        
        print("\n" + "="*80)
        print(f"SQL IMPORT OPERATIONS REGISTER - Project: {self.name}")
        print("="*80)
        
        for i, reg in enumerate(registers, 1):
            success_icon = "✅" if reg.get('success') else "❌"
            print(f"\n{success_icon} Operation {i} - {reg.get('creation_date', 'Unknown date')}")
            print(f"   📋 Target Table: {reg.get('target_table_name', 'N/A')}")
            print(f"   🗄️  Source Database: {reg.get('source_db_path', 'N/A')}")
            print(f"   📊 Results: {reg.get('compounds_added', 0)} added, "
                  f"{reg.get('duplicates_skipped', 0)} duplicates, "
                  f"{reg.get('errors', 0)} errors")
            print(f"   🔬 InChI Keys: {reg.get('inchi_keys_computed', 0)} computed")
            
            # Show truncated SQL query
            sql_query = reg.get('sql_query', '')
            if len(sql_query) > 100:
                sql_display = sql_query[:97] + "..."
            else:
                sql_display = sql_query
            print(f"   🔍 Query: {sql_display}")
        
        print("="*80)
    
    def create_table_from_sql(self, sql_query: str, 
                             source_db_path: str,
                             target_table_name: str,
                             smiles_column: str = 'smiles',
                             name_column: Optional[str] = 'name', 
                             flag_column: Optional[str] = 'flag',
                             compute_inchi: bool = True) -> Dict[str, Any]:
        """
        Create a molecules table in the chemspace database from an SQL query executed on an external database.
        
        Args:
            sql_query (str): SQL SELECT query to execute on the source database
            source_db_path (str): Path to the source SQLite database file
            target_table_name (str): Name for the new table in the chemspace database
            smiles_column (str): Name of the column containing SMILES strings
            name_column (Optional[str]): Name of the column containing compound names. If None or not found, fills with "nd"
            flag_column (Optional[str]): Name of the column containing flags. If None or not found, fills with "nd"
            compute_inchi (bool): Whether to automatically compute InChI keys after loading
            
        Returns:
            dict: Results containing success status, counts, and messages
        """
        try:
            # Validate source database exists
            if not os.path.exists(source_db_path):
                return {
                    'success': False,
                    'message': f"Source database not found: {source_db_path}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
            # Sanitize target table name
            target_table_name = self._sanitize_table_name(target_table_name)
            
            print(f"📊 Executing SQL query on source database: {source_db_path}")
            print(f"📋 Target table: {target_table_name}")
            print(f"🔍 Query: {sql_query}")
            
            # Connect to source database and execute query
            source_conn = sqlite3.connect(source_db_path)
            
            try:
                df = pd.read_sql_query(sql_query, source_conn)
                source_conn.close()
            except Exception as e:
                source_conn.close()
                return {
                    'success': False,
                    'message': f"Error executing SQL query: {e}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
            if df.empty:
                return {
                    'success': False,
                    'message': "SQL query returned no results",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
            print(f"📊 Query returned {len(df)} rows")
            
            # Validate required columns exist in query result (only SMILES is mandatory)
            if smiles_column not in df.columns:
                return {
                    'success': False,
                    'message': f"Query result missing required SMILES column: {smiles_column}. Available columns: {list(df.columns)}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
            # Check if optional columns exist and warn if they don't
            name_available = name_column and name_column in df.columns
            flag_available = flag_column and flag_column in df.columns
            
            if not name_available:
                print(f"⚠️  Name column '{name_column}' not found in query result. Will use 'nd' as default.")
            if not flag_available:
                print(f"⚠️  Flag column '{flag_column}' not found in query result. Will use 'nd' as default.")
            
            # Create target table in chemspace database
            if not self._create_compounds_table(target_table_name):
                return {
                    'success': False,
                    'message': f"Failed to create target table '{target_table_name}'",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
            # Prepare data for insertion
            compounds_data = []
            skipped_rows = 0
            
            for _, row in df.iterrows():
                smiles = str(row[smiles_column]).strip()
                name = str(row[name_column]).strip()
                flag = str(row.get(flag_column, '')).strip() if flag_column in df.columns else ''
                
                # Skip empty or invalid rows
                if not smiles or not name or smiles.lower() in ['nan', 'none', '']:
                    skipped_rows += 1
                    continue
                
                compounds_data.append((smiles, name, flag))
            
            if skipped_rows > 0:
                print(f"⚠️  Skipped {skipped_rows} rows with missing/invalid data")
            
            if not compounds_data:
                return {
                    'success': False,
                    'message': "No valid compound data found in query results",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
            # Insert compounds into chemspace database
            result = self._insert_compounds(compounds_data, target_table_name, skip_duplicates=True)
            
            if result['success']:
                self._compounds_loaded = True
                table_count = self.get_compound_count(table_name=target_table_name)
                
                print(f"✅ Successfully created table '{target_table_name}' from SQL query")
                print(f"   📊 Compounds added: {result['compounds_added']}")
                print(f"   🔄 Duplicates skipped: {result['duplicates_skipped']}")
                print(f"   ❌ Errors: {result['errors']}")
                print(f"   📈 Total compounds in table: {table_count}")
                
                # Automatically compute InChI keys if requested
                if compute_inchi and result['compounds_added'] > 0:
                    print(f"\n🧪 Computing InChI keys for imported compounds...")
                    try:
                        inchi_df = self.compute_inchi_keys(target_table_name, update_database=True)
                        if not inchi_df.empty:
                            # Count successful InChI computations
                            valid_inchi_count = len(inchi_df[
                                (inchi_df['inchi_key'].notna()) & 
                                (~inchi_df['inchi_key'].isin(['INVALID_SMILES', 'ERROR']))
                            ])
                            result['inchi_keys_computed'] = valid_inchi_count
                            print(f"   🔬 InChI keys computed: {valid_inchi_count}")
                        else:
                            result['inchi_keys_computed'] = 0
                            print(f"   ⚠️  No InChI keys computed")
                    except Exception as inchi_error:
                        print(f"   ⚠️  Warning: InChI key computation failed: {inchi_error}")
                        result['inchi_keys_computed'] = 0
                else:
                    result['inchi_keys_computed'] = 0
            
            result['table_name'] = target_table_name
            result['source_database'] = source_db_path
            result['rows_processed'] = len(df)
            result['rows_skipped'] = skipped_rows
            
            # Register the SQL operation
            self._register_sql_operation(sql_query, source_db_path, target_table_name, result)
            
            return result
            
        except Exception as e:
            error_result = {
                'success': False,
                'message': f"Error creating table from SQL query: {e}",
                'compounds_added': 0,
                'duplicates_skipped': 0,
                'errors': 1,
                'inchi_keys_computed': 0
            }
            
            # Register the failed SQL operation
            try:
                self._register_sql_operation(sql_query, source_db_path, target_table_name, error_result)
            except:
                pass  # Don't fail if registration fails
            
            return error_result
    
    def execute_custom_query(self, sql_query: str, 
                           source_db_path: Optional[str] = None) -> pd.DataFrame:
        """
        Execute a custom SQL query on either the chemspace database or an external database.
        
        Args:
            sql_query (str): SQL query to execute
            source_db_path (Optional[str]): Path to external database. If None, uses chemspace database
            
        Returns:
            pd.DataFrame: Query results as DataFrame, empty DataFrame if error occurs
        """
        try:
            # Determine which database to use
            db_path = source_db_path if source_db_path else self.__chemspace_db
            
            if not os.path.exists(db_path):
                print(f"❌ Database not found: {db_path}")
                return pd.DataFrame()
            
            print(f"📊 Executing custom query on: {db_path}")
            print(f"🔍 Query: {sql_query}")
            
            # Connect to database and execute query
            conn = sqlite3.connect(db_path)
            
            try:
                df = pd.read_sql_query(sql_query, conn)
                conn.close()
                
                print(f"✅ Query executed successfully, returned {len(df)} rows")
                return df
                
            except Exception as e:
                conn.close()
                print(f"❌ Error executing query: {e}")
                return pd.DataFrame()
            
        except Exception as e:
            print(f"❌ Error with custom query execution: {e}")
            return pd.DataFrame()
    
    def _process_csv_parallel(self, df: pd.DataFrame, table_name: str,
                             smiles_column: str, name_column: Optional[str], flag_column: Optional[str],
                             name_available: bool, flag_available: bool, skip_duplicates: bool,
                             max_workers: int, chunk_size: int) -> Dict[str, Any]:
        """
        Process CSV data using parallel processing with chunks and comprehensive progress tracking.
        
        Args:
            df (pd.DataFrame): DataFrame containing the CSV data
            table_name (str): Name of the target table
            smiles_column (str): Name of the SMILES column
            name_column (Optional[str]): Name of the name column
            flag_column (Optional[str]): Name of the flag column
            name_available (bool): Whether name column is available
            flag_available (bool): Whether flag column is available
            skip_duplicates (bool): Whether to skip duplicate compounds
            max_workers (int): Maximum number of parallel workers
            chunk_size (int): Size of each chunk
            
        Returns:
            dict: Results containing counts and status
        """
        try:
            # Record start time for performance tracking
            start_time = time.time()
            
            # Split DataFrame into chunks
            print(f"   📦 Splitting {len(df)} rows into chunks...")
            chunks = self._split_dataframe_into_chunks(df, chunk_size)
            num_chunks = len(chunks)
            
            print(f"   📊 Parallel Processing Setup:")
            print(f"      📦 Total chunks: {num_chunks}")
            print(f"      📏 Chunk size: {chunk_size}")
            print(f"      👥 Workers: {max_workers}")
            
            # Initialize progress tracking
            total_compounds_added = 0
            total_duplicates_skipped = 0
            total_errors = 0
            processed_rows = 0
            
            # Initialize progress bar if tqdm is available
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(df),
                    desc="Processing compounds",
                    unit="compounds",
                    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
                )
            else:
                progress_bar = None
                print(f"   🚀 Starting parallel processing of {num_chunks} chunks...")
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all chunk processing jobs with timing
                future_to_chunk = {}
                chunk_submission_start = time.time()
                
                for i, chunk in enumerate(chunks):
                    future = executor.submit(
                        _process_chunk_worker,
                        chunk, smiles_column, name_column, flag_column,
                        name_available, flag_available, i * chunk_size
                    )
                    future_to_chunk[future] = {
                        'index': i,
                        'size': len(chunk),
                        'submitted_at': time.time()
                    }
                
                chunk_submission_time = time.time() - chunk_submission_start
                print(f"   ⏱️  All {num_chunks} chunks submitted in {chunk_submission_time:.2f}s")
                
                # Collect results and insert into database with detailed progress
                processed_chunks = 0
                failed_chunks = 0
                
                for future in as_completed(future_to_chunk):
                    chunk_info = future_to_chunk[future]
                    chunk_idx = chunk_info['index']
                    chunk_size_actual = chunk_info['size']
                    
                    try:
                        # Process chunk result
                        chunk_start_time = time.time()
                        compounds_data = future.result()
                        chunk_process_time = time.time() - chunk_start_time
                        
                        if compounds_data:
                            # Insert this chunk's data into database
                            insert_start_time = time.time()
                            chunk_result = self._insert_compounds(compounds_data, table_name, skip_duplicates)
                            insert_time = time.time() - insert_start_time
                            
                            if chunk_result['success']:
                                total_compounds_added += chunk_result['compounds_added']
                                total_duplicates_skipped += chunk_result['duplicates_skipped']
                                total_errors += chunk_result['errors']
                                
                                # Detailed chunk statistics (only show every 5th chunk to avoid spam)
                                if not TQDM_AVAILABLE and (processed_chunks + 1) % 5 == 0:
                                    print(f"   ✅ Chunk {chunk_idx + 1}/{num_chunks}: "
                                          f"{chunk_result['compounds_added']} added, "
                                          f"{chunk_result['duplicates_skipped']} duplicates, "
                                          f"{chunk_result['errors']} errors "
                                          f"(Process: {chunk_process_time:.2f}s, Insert: {insert_time:.2f}s)")
                        
                        processed_rows += chunk_size_actual
                        processed_chunks += 1
                        
                        # Update progress bar or print progress
                        if progress_bar:
                            progress_bar.update(chunk_size_actual)
                            progress_bar.set_postfix({
                                'chunks': f"{processed_chunks}/{num_chunks}",
                                'added': total_compounds_added,
                                'duplicates': total_duplicates_skipped,
                                'errors': total_errors
                            })
                        elif processed_chunks % 10 == 0 or processed_chunks == num_chunks:
                            elapsed_time = time.time() - start_time
                            remaining_chunks = num_chunks - processed_chunks
                            estimated_remaining = (elapsed_time / processed_chunks) * remaining_chunks if processed_chunks > 0 else 0
                            
                            print(f"   ⏳ Progress: {processed_chunks}/{num_chunks} chunks "
                                  f"({processed_rows:,}/{len(df):,} rows) "
                                  f"[{elapsed_time:.1f}s elapsed, ~{estimated_remaining:.1f}s remaining]")
                            print(f"      📊 Current totals: {total_compounds_added:,} added, "
                                  f"{total_duplicates_skipped:,} duplicates, {total_errors:,} errors")
                        
                    except Exception as e:
                        failed_chunks += 1
                        total_errors += chunk_size_actual  # Estimate errors for failed chunk
                        processed_chunks += 1
                        processed_rows += chunk_size_actual
                        
                        print(f"   ❌ Chunk {chunk_idx + 1} failed: {e}")
                        
                        # Update progress even for failed chunks
                        if progress_bar:
                            progress_bar.update(chunk_size_actual)
                            progress_bar.set_postfix({
                                'chunks': f"{processed_chunks}/{num_chunks}",
                                'failed': failed_chunks,
                                'errors': total_errors
                            })
                
                # Close progress bar
                if progress_bar:
                    progress_bar.close()
            
            # Final timing and statistics
            total_time = time.time() - start_time
            rows_per_second = len(df) / total_time if total_time > 0 else 0
            
            print(f"\n   🏁 Parallel Processing Complete!")
            print(f"      ⏱️  Total time: {total_time:.2f}s")
            print(f"      📈 Processing rate: {rows_per_second:,.0f} rows/second")
            print(f"      ✅ Successful chunks: {processed_chunks - failed_chunks}/{num_chunks}")
            if failed_chunks > 0:
                print(f"      ❌ Failed chunks: {failed_chunks}/{num_chunks}")
            print(f"      📊 Final counts: {total_compounds_added:,} added, "
                  f"{total_duplicates_skipped:,} duplicates, {total_errors:,} errors")
            
            return {
                'success': True,
                'message': f"Parallel processing completed: {total_compounds_added} compounds added",
                'compounds_added': total_compounds_added,
                'duplicates_skipped': total_duplicates_skipped,
                'errors': total_errors
            }
            
        except Exception as e:
            return {
                'success': False,
                'message': f"Parallel processing error: {e}",
                'compounds_added': 0,
                'duplicates_skipped': 0,
                'errors': len(df)
            }
    
    def _split_dataframe_into_chunks(self, df: pd.DataFrame, chunk_size: int) -> List[pd.DataFrame]:
        """
        Split a DataFrame into chunks of specified size.
        
        Args:
            df (pd.DataFrame): DataFrame to split
            chunk_size (int): Size of each chunk
            
        Returns:
            List[pd.DataFrame]: List of DataFrame chunks
        """
        chunks = []
        num_chunks = len(df) // chunk_size + (1 if len(df) % chunk_size != 0 else 0)
        
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, len(df))
            chunk = df.iloc[start_idx:end_idx].copy()
            chunks.append(chunk)
        
        return chunks


