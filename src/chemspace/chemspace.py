import os
import csv
import sqlite3
import pandas as pd
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any
import sys

# Add the parent directory to path to import our local tidyscreen module
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir)

# Import from our local tidyscreen module
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject

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
    
    def describe_table(self, table_name: str) -> Dict[str, Any]:
        """
        Get detailed information about a specific table.
        
        Args:
            table_name (str): Name of the table to describe
            
        Returns:
            Dict: Table information including schema, stats, and sample data
        """
        try:
            tables = self.get_all_tables()
            if table_name not in tables:
                print(f"⚠️  Table '{table_name}' does not exist")
                return {}
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get table schema
            cursor.execute(f"PRAGMA table_info({table_name})")
            schema = cursor.fetchall()
            
            # Get sample data
            cursor.execute(f"SELECT * FROM {table_name} LIMIT 3")
            sample_rows = cursor.fetchall()
            column_names = [description[0] for description in cursor.description]
            
            # Get statistics
            stats = self.get_database_stats(table_name)
            
            conn.close()
            
            # Format schema information
            schema_info = []
            for col_info in schema:
                schema_info.append({
                    'column_id': col_info[0],
                    'name': col_info[1],
                    'type': col_info[2],
                    'not_null': bool(col_info[3]),
                    'default_value': col_info[4],
                    'primary_key': bool(col_info[5])
                })
            
            # Format sample data
            sample_data = []
            for row in sample_rows:
                row_dict = {}
                for i, col_name in enumerate(column_names):
                    row_dict[col_name] = row[i]
                sample_data.append(row_dict)
            
            return {
                'table_name': table_name,
                'schema': schema_info,
                'statistics': stats,
                'sample_data': sample_data
            }
            
        except Exception as e:
            print(f"❌ Error describing table '{table_name}': {e}")
            return {}
    
    def print_table_info(self, table_name: str) -> None:
        """
        Print formatted information about a specific table.
        
        Args:
            table_name (str): Name of the table to describe
        """
        info = self.describe_table(table_name)
        
        if not info:
            return
        
        print("\n" + "="*60)
        print(f"TABLE INFORMATION: {table_name}")
        print("="*60)
        
        # Schema
        print("📋 Schema:")
        for col in info['schema']:
            pk_marker = " (PK)" if col['primary_key'] else ""
            null_marker = " NOT NULL" if col['not_null'] else ""
            default = f" DEFAULT {col['default_value']}" if col['default_value'] else ""
            print(f"   {col['name']}: {col['type']}{pk_marker}{null_marker}{default}")
        
        # Statistics
        if info['statistics']:
            stats = info['statistics']
            print(f"\n📊 Statistics:")
            print(f"   Total compounds: {stats.get('total_compounds', 0)}")
            if stats.get('flag_counts'):
                print("   Compounds by flag:")
                for flag, count in stats['flag_counts'].items():
                    flag_display = flag if flag else '(no flag)'
                    print(f"      {flag_display}: {count}")
        
        # Sample data
        if info['sample_data']:
            print(f"\n🔍 Sample data (first 3 rows):")
            for i, row in enumerate(info['sample_data'], 1):
                print(f"   Row {i}:")
                for col, value in row.items():
                    if col == 'smiles':
                        # Truncate long SMILES for display
                        display_value = value[:50] + "..." if len(str(value)) > 50 else value
                    else:
                        display_value = value
                    print(f"      {col}: {display_value}")
        
        print("="*60)
    
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
                     name_column: str = 'name', 
                     flag_column: str = 'flag',
                     skip_duplicates: bool = True,
                     compute_inchi: bool = True) -> Dict[str, Any]:
        """
        Load compounds from a CSV file into the chemspace database.
        Creates a table named after the CSV file prefix.
        Automatically computes InChI keys for loaded SMILES.
        
        Args:
            csv_file_path (str): Path to the CSV file
            smiles_column (str): Name of the column containing SMILES strings
            name_column (str): Name of the column containing compound names
            flag_column (str): Name of the column containing flags
            skip_duplicates (bool): Whether to skip duplicate compounds
            compute_inchi (bool): Whether to automatically compute InChI keys after loading
            
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
            
            # Validate required columns
            missing_columns = []
            if smiles_column not in df.columns:
                missing_columns.append(smiles_column)
            if name_column not in df.columns:
                missing_columns.append(name_column)
            if flag_column not in df.columns and flag_column != 'flag':
                missing_columns.append(flag_column)
            
            if missing_columns:
                return {
                    'success': False,
                    'message': f"Missing required columns: {missing_columns}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1
                }
            
            # Prepare data for insertion
            compounds_data = []
            
            for _, row in df.iterrows():
                smiles = str(row[smiles_column]).strip()
                name = str(row[name_column]).strip()
                flag = str(row.get(flag_column, '')).strip() if flag_column in df.columns else ''
                
                # Skip empty rows
                if not smiles or not name or smiles.lower() in ['nan', 'none', '']:
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
    
    def search_compounds(self, search_term: str, search_in: str = 'name', table_name: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search for compounds by name, SMILES, flag, or InChI key pattern.
        
        Args:
            search_term (str): Term to search for
            search_in (str): Field to search in ('name', 'smiles', 'flag', 'inchi_key')
            table_name (Optional[str]): Name of the table to search in. If None, searches all tables
            
        Returns:
            List[Dict]: List of matching compounds
        """
        valid_fields = ['name', 'smiles', 'flag', 'inchi_key']
        if search_in not in valid_fields:
            print(f"❌ Invalid search field. Use one of: {valid_fields}")
            return []
        
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            if table_name:
                # Search in specific table
                query = f"SELECT *, '{table_name}' as table_name FROM {table_name} WHERE {search_in} LIKE ? ORDER BY id DESC"
                cursor.execute(query, (f"%{search_term}%",))
                rows = cursor.fetchall()
            else:
                # Search in all tables
                tables = self.get_all_tables()
                all_rows = []
                
                for table in tables:
                    try:
                        query = f"SELECT *, '{table}' as table_name FROM {table} WHERE {search_in} LIKE ? ORDER BY id DESC"
                        cursor.execute(query, (f"%{search_term}%",))
                        table_rows = cursor.fetchall()
                        all_rows.extend(table_rows)
                    except Exception as e:
                        print(f"⚠️  Warning: Error searching table '{table}': {e}")
                
                rows = all_rows
            
            # Get column names
            if rows:
                cursor.execute(f"SELECT *, 'dummy' as table_name FROM {tables[0] if not table_name else table_name} LIMIT 1")
                column_names = [description[0] for description in cursor.description]
            else:
                conn.close()
                return []
            
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
            print(f"❌ Error searching compounds: {e}")
            return []
    
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
    
    def get_database_stats(self, table_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Get statistics about the chemspace database.
        
        Args:
            table_name (Optional[str]): Specific table to analyze. If None, analyzes all tables.
        
        Returns:
            Dict: Database statistics including counts by flag, total compounds, etc.
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            if table_name:
                # Stats for specific table
                cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
                total_compounds = cursor.fetchone()[0]
                
                cursor.execute(f"SELECT flag, COUNT(*) FROM {table_name} GROUP BY flag ORDER BY COUNT(*) DESC")
                flag_counts = dict(cursor.fetchall())
                
                return {
                    'table_name': table_name,
                    'total_compounds': total_compounds,
                    'flag_counts': flag_counts,
                    'database_path': self.__chemspace_db
                }
            else:
                # Stats for all tables
                tables = self.get_all_tables()
                all_stats = {}
                total_all_compounds = 0
                
                for table in tables:
                    try:
                        cursor.execute(f"SELECT COUNT(*) FROM {table}")
                        table_count = cursor.fetchone()[0]
                        total_all_compounds += table_count
                        
                        cursor.execute(f"SELECT flag, COUNT(*) FROM {table} GROUP BY flag ORDER BY COUNT(*) DESC")
                        table_flag_counts = dict(cursor.fetchall())
                        
                        all_stats[table] = {
                            'total_compounds': table_count,
                            'flag_counts': table_flag_counts
                        }
                    except Exception as e:
                        print(f"⚠️  Warning: Error analyzing table '{table}': {e}")
                
                conn.close()
                
                return {
                    'total_tables': len(tables),
                    'total_compounds_all_tables': total_all_compounds,
                    'tables': all_stats,
                    'database_path': self.__chemspace_db
                }
            
            conn.close()
            
        except Exception as e:
            print(f"❌ Error getting database stats: {e}")
            return {}
    
    def print_database_summary(self, table_name: Optional[str] = None) -> None:
        """
        Print a formatted summary of the chemspace database.
        
        Args:
            table_name (Optional[str]): Specific table to summarize. If None, summarizes all tables.
        """
        stats = self.get_database_stats(table_name)
        
        if not stats:
            print("❌ Unable to retrieve database statistics")
            return
        
        print("\n" + "="*60)
        print(f"CHEMSPACE DATABASE SUMMARY - Project: {self.name}")
        print("="*60)
        
        if table_name:
            # Summary for specific table
            print(f"📋 Table: {stats['table_name']}")
            print(f"📊 Total Compounds: {stats['total_compounds']}")
            
            if stats['flag_counts']:
                print("\n🏷️  Compounds by Flag:")
                for flag, count in stats['flag_counts'].items():
                    flag_display = flag if flag else '(no flag)'
                    print(f"   {flag_display}: {count}")
        else:
            # Summary for all tables
            print(f"📊 Total Tables: {stats['total_tables']}")
            print(f"📈 Total Compounds (All Tables): {stats['total_compounds_all_tables']}")
            
            if stats['tables']:
                print("\n📋 Tables:")
                for table, table_stats in stats['tables'].items():
                    print(f"   {table}: {table_stats['total_compounds']} compounds")
                    if table_stats['flag_counts']:
                        for flag, count in table_stats['flag_counts'].items():
                            flag_display = flag if flag else '(no flag)'
                            print(f"      └─ {flag_display}: {count}")
        
        print(f"\n🗄️  Database Path: {stats['database_path']}")
        print("="*60)
    
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
                          update_database: bool = True) -> pd.DataFrame:
        """
        Retrieve a table as DataFrame and compute InChI keys for each SMILES string.
        
        Args:
            table_name (str): Name of the table to retrieve
            flag_filter (Optional[str]): Filter by specific flag value
            name_pattern (Optional[str]): Filter by name pattern (SQL LIKE)
            update_database (bool): Whether to add inchi_key column to the database table
            
        Returns:
            pd.DataFrame: DataFrame with added 'inchi_key' column, empty DataFrame if error occurs
        """
        try:
            # Import RDKit
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
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
            
            print(f"🔬 Computing InChI keys for {len(df)} compounds...")
            
            # Compute InChI keys for each SMILES
            successful_computations = 0
            failed_computations = 0
            
            for idx, row in df.iterrows():
                smiles = row['smiles']
                
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
                            print(f"⚠️  Invalid SMILES for compound '{row.get('name', 'unknown')}': {smiles}")
                
                except Exception as e:
                    df.at[idx, 'inchi_key'] = 'ERROR'
                    failed_computations += 1
                    if failed_computations <= 5:  # Show only first 5 errors
                        print(f"⚠️  Error computing InChI key for '{row.get('name', 'unknown')}': {e}")
            
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
                             name_column: str = 'name', 
                             flag_column: str = 'flag',
                             compute_inchi: bool = True) -> Dict[str, Any]:
        """
        Create a molecules table in the chemspace database from an SQL query executed on an external database.
        
        Args:
            sql_query (str): SQL SELECT query to execute on the source database
            source_db_path (str): Path to the source SQLite database file
            target_table_name (str): Name for the new table in the chemspace database
            smiles_column (str): Name of the column containing SMILES strings
            name_column (str): Name of the column containing compound names
            flag_column (str): Name of the column containing flags (optional)
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
            
            # Validate required columns exist in query result
            missing_columns = []
            if smiles_column not in df.columns:
                missing_columns.append(smiles_column)
            if name_column not in df.columns:
                missing_columns.append(name_column)
            # flag_column is optional
            
            if missing_columns:
                return {
                    'success': False,
                    'message': f"Query result missing required columns: {missing_columns}. Available columns: {list(df.columns)}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1,
                    'inchi_keys_computed': 0
                }
            
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


