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
import json
from tidyscreen import tidyscreen


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
#from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject

## Module level helper functions for multiprocessing workers

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

def _filter_chunk_worker_by_instances(chunk_data: List[Tuple], filters_list: List[Tuple]) -> List[Dict]:
    """
    Memory-optimized version that only returns essential data.
    """
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        return [{"error": "RDKit not available"}]
    
    results = []
    sorted_filters = sorted(filters_list, key=lambda x: (x[1] == 0, x[1]))
    
    for compound_id, smiles, name, flag, inchi_key in chunk_data:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue  # Skip invalid compounds entirely
            
            passes_all_filters = True
            
            # Apply filters with early termination but minimal metadata
            for smarts_pattern, required_instances in sorted_filters:
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        passes_all_filters = False
                        break
                    
                    matches = mol.GetSubstructMatches(pattern)
                    match_count = len(matches)
                    
                    if match_count != required_instances:
                        passes_all_filters = False
                        break  # Early exit - don't process more filters
                        
                except Exception:
                    passes_all_filters = False
                    break
            
            # Only store compounds that pass ALL filters
            if passes_all_filters:
                results.append({
                    'id': compound_id,
                    'smiles': smiles,
                    'name': name or 'unknown',
                    'flag': flag or 'nd',
                    'inchi_key': inchi_key
                    # Remove all the heavy metadata
                })
                
        except Exception:
            continue  # Skip problematic compounds
    
    return results

def _process_bimolecular_chunk_worker(chunk_data: List[Tuple], reaction_smarts: str, 
                                    reaction_name: str, workflow_name: str) -> List[Dict[str, Any]]:
    """
    Worker function to process a chunk of bimolecular reaction combinations.
    This function must be at module level to be pickleable for multiprocessing.
    
    Args:
        chunk_data (List[Tuple]): List of (primary_compound, secondary_compound) tuples
        reaction_smarts (str): SMARTS pattern for the reaction
        reaction_name (str): Name of the reaction
        workflow_name (str): Name of the workflow
        
    Returns:
        List[Dict]: List of reaction products
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        
        products = []
        
        # Parse reaction
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return []
        
        for primary_compound, secondary_compound in chunk_data:
            try:
                # Parse molecules
                primary_mol = Chem.MolFromSmiles(primary_compound['smiles'])
                secondary_mol = Chem.MolFromSmiles(secondary_compound['smiles'])
                
                if primary_mol is None or secondary_mol is None:
                    continue
                
                # Run reaction
                reaction_results = rxn.RunReactants((primary_mol, secondary_mol))
                
                # Process products
                for product_set_idx, product_set in enumerate(reaction_results):
                    for product_idx, product_mol in enumerate(product_set):
                        try:
                            Chem.SanitizeMol(product_mol)
                            product_smiles = Chem.MolToSmiles(product_mol)
                            
                            # Generate product name
                            primary_name = primary_compound.get('name', f"cpd_{primary_compound.get('id', 'unk')}")
                            secondary_name = secondary_compound.get('name', f"cpd_{secondary_compound.get('id', 'unk')}")
                            product_name = f"{primary_name}+{secondary_name}_{reaction_name}_{product_set_idx}_{product_idx}"
                            
                            products.append({
                                'smiles': product_smiles,
                                'name': product_name,
                                'flag': 'parallel_bimolecular_product',
                                'reactant1_name': primary_name,
                                'reactant1_smiles': primary_compound['smiles'],
                                'reactant2_name': secondary_name,
                                'reactant2_smiles': secondary_compound['smiles'],
                                'reaction_name': reaction_name,
                                'workflow': workflow_name
                            })
                            
                        except Exception:
                            continue
            
            except Exception:
                continue
        
        return products
        
    except Exception:
        return []

def _process_unimolecular_chunk_worker(chunk_data: List[Dict], reaction_smarts: str, 
                                     reaction_name: str, workflow_name: str, 
                                     name_prefix: str) -> List[Dict[str, Any]]:
    """
    Worker function to process a chunk of unimolecular reactions.
    This function must be at module level to be pickleable for multiprocessing.
    
    Args:
        chunk_data (List[Dict]): List of compound dictionaries
        reaction_smarts (str): SMARTS pattern for the reaction
        reaction_name (str): Name of the reaction
        workflow_name (str): Name of the workflow
        name_prefix (str): Prefix for product names
        
    Returns:
        List[Dict]: List of reaction products
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        
        products = []
        
        # Parse reaction
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return []
        
        for compound in chunk_data:
            try:
                # Parse molecule
                reactant_mol = Chem.MolFromSmiles(compound['smiles'])
                if reactant_mol is None:
                    continue
                
                # Run reaction
                reaction_results = rxn.RunReactants((reactant_mol,))
                
                # Process products
                for product_set_idx, product_set in enumerate(reaction_results):
                    for product_idx, product_mol in enumerate(product_set):
                        try:
                            Chem.SanitizeMol(product_mol)
                            product_smiles = Chem.MolToSmiles(product_mol)
                            
                            # Generate product name
                            original_name = compound.get('name', f"cpd_{compound.get('id', 'unk')}")
                            product_name = f"{name_prefix}{original_name}_{reaction_name}_{product_set_idx}_{product_idx}"
                            
                            products.append({
                                'smiles': product_smiles,
                                'name': product_name,
                                'flag': 'parallel_unimolecular_product',
                                'reactant_name': original_name,
                                'reactant_smiles': compound['smiles'],
                                'reaction_name': reaction_name,
                                'workflow': workflow_name
                            })
                            
                        except Exception:
                            continue
                            
            except Exception:
                continue
        
        return products
        
    except Exception:
        return []

def _process_bimolecular_chunk_to_file_worker(chunk_data: List[Tuple], reaction_smarts: str, 
                                            reaction_name: str, workflow_name: str, 
                                            output_file_path: str) -> int:
    """
    Worker function to process a bimolecular chunk and write results directly to a CSV file.
    
    Args:
        chunk_data (List[Tuple]): List of (primary_compound, secondary_compound) tuples
        reaction_smarts (str): SMARTS pattern for the reaction
        reaction_name (str): Name of the reaction
        workflow_name (str): Name of the workflow
        output_file_path (str): Path to output CSV file
        
    Returns:
        int: Number of products written to file
    """
    try:
        import csv
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        
        products_count = 0
        
        # Parse reaction
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return 0
        
        # Open file for writing
        with open(output_file_path, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['smiles', 'name', 'flag', 'reactant1_name', 'reactant1_smiles',
                         'reactant2_name', 'reactant2_smiles', 'reaction_name', 'workflow']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for primary_compound, secondary_compound in chunk_data:
                try:
                    # Parse molecules
                    primary_mol = Chem.MolFromSmiles(primary_compound['smiles'])
                    secondary_mol = Chem.MolFromSmiles(secondary_compound['smiles'])
                    
                    if primary_mol is None or secondary_mol is None:
                        continue
                    
                    # Run reaction
                    reaction_results = rxn.RunReactants((primary_mol, secondary_mol))
                    
                    # Process products and write immediately
                    for product_set_idx, product_set in enumerate(reaction_results):
                        for product_idx, product_mol in enumerate(product_set):
                            try:
                                Chem.SanitizeMol(product_mol)
                                product_smiles = Chem.MolToSmiles(product_mol)
                                
                                # Generate product name
                                primary_name = primary_compound.get('name', f"cpd_{primary_compound.get('id', 'unk')}")
                                secondary_name = secondary_compound.get('name', f"cpd_{secondary_compound.get('id', 'unk')}")
                                product_name = f"{primary_name}+{secondary_name}_{reaction_name}_{product_set_idx}_{product_idx}"
                                
                                # Write product directly to file
                                writer.writerow({
                                    'smiles': product_smiles,
                                    'name': product_name,
                                    'flag': 'stream_bimolecular_product',
                                    'reactant1_name': primary_name,
                                    'reactant1_smiles': primary_compound['smiles'],
                                    'reactant2_name': secondary_name,
                                    'reactant2_smiles': secondary_compound['smiles'],
                                    'reaction_name': reaction_name,
                                    'workflow': workflow_name
                                })
                                products_count += 1
                                
                            except Exception:
                                continue
                                
                except Exception:
                    continue
        
        return products_count
        
    except Exception:
        return 0

def _process_unimolecular_chunk_to_file_worker(chunk_data: List[Dict], reaction_smarts: str, 
                                             reaction_name: str, workflow_name: str,
                                             name_prefix: str, output_file_path: str) -> int:
    """
    Worker function to process a unimolecular chunk and write results directly to a CSV file.
    
    Args:
        chunk_data (List[Dict]): List of compound dictionaries
        reaction_smarts (str): SMARTS pattern for the reaction
        reaction_name (str): Name of the reaction
        workflow_name (str): Name of the workflow
        name_prefix (str): Prefix for product names
        output_file_path (str): Path to output CSV file
        
    Returns:
        int: Number of products written to file
    """
    try:
        import csv
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        
        products_count = 0
        
        # Parse reaction
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return 0
        
        # Open file for writing
        with open(output_file_path, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['smiles', 'name', 'flag', 'reactant_name', 
                         'reactant_smiles', 'reaction_name', 'workflow']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for compound in chunk_data:
                try:
                    # Parse molecule
                    reactant_mol = Chem.MolFromSmiles(compound['smiles'])
                    if reactant_mol is None:
                        continue
                    
                    # Run reaction
                    reaction_results = rxn.RunReactants((reactant_mol,))
                    
                    # Process products and write immediately
                    for product_set_idx, product_set in enumerate(reaction_results):
                        for product_idx, product_mol in enumerate(product_set):
                            try:
                                Chem.SanitizeMol(product_mol)
                                product_smiles = Chem.MolToSmiles(product_mol)
                                
                                # Generate product name
                                original_name = compound.get('name', f"cpd_{compound.get('id', 'unk')}")
                                product_name = f"{name_prefix}{original_name}_{reaction_name}_{product_set_idx}_{product_idx}"
                                
                                # Write product directly to file
                                writer.writerow({
                                    'smiles': product_smiles,
                                    'name': product_name,
                                    'flag': 'stream_unimolecular_product',
                                    'reactant_name': original_name,
                                    'reactant_smiles': compound['smiles'],
                                    'reaction_name': reaction_name,
                                    'workflow': workflow_name
                                })
                                products_count += 1
                                
                            except Exception:
                                continue
                                
                except Exception:
                    continue
        
        return products_count
        
    except Exception:
        return 0

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
                flag_description TEXT,
                inchi_key TEXT,
                UNIQUE(name, smiles)
            )
            """
            
            cursor.execute(create_table_query)
            
            conn.commit()
            conn.close()
            
            return True
            
        except Exception as e:
            print(f"❌ Error creating table '{table_name}': {e}")
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
    
    def load_csv_file(self, csv_file_path: str = None, 
                     smiles_column: str = 'smiles', 
                     name_column: str = 'name', 
                     flag_column: str = 'flag',
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
            # Prompt for file path if not provided
            if not csv_file_path:
                csv_file_path = input("Enter the path to the CSV file: ").strip()
                while not csv_file_path:
                    print("❌ CSV file path cannot be empty.")
                    csv_file_path = input("Enter the path to the CSV file: ").strip()
            # Validate file exists
            if not os.path.exists(csv_file_path):
                return {
                    'success': False,
                    'message': f"CSV file not found: {csv_file_path}",
                    'compounds_added': 0,
                    'duplicates_skipped': 0,
                    'errors': 1
                }
            
            ## Prepare table name asking the user the intended behavior
            csv_filename = os.path.basename(csv_file_path)
            default_table_name = os.path.splitext(csv_filename)[0]  # Remove .csv extension
            print(f"Default table name (from file): '{default_table_name}'")
            use_custom = input("Do you want to input a custom table name? (y/N): ").strip().lower()
            if use_custom in ['y', 'yes']:
                custom_name = ''
                while not custom_name:
                    custom_name = input("Enter custom table name: ").strip()
                    if not custom_name:
                        table_name = default_table_name
                        print(f"Using custom table name: {default_table_name}.")
                table_name = custom_name
            else:
                table_name = default_table_name
            
            # Sanitize table name for SQLite (remove special characters, ensure it starts with letter)
            table_name = self._sanitize_table_name(table_name)

            # Check if table already exists
            existing_tables = self.get_all_tables()
            if table_name in existing_tables:
                print(f"⚠️  Table '{table_name}' already exists in the database.")
                replace = input(f"Do you want to replace the existing table '{table_name}'? (y/N): ").strip().lower()
                if replace in ['y', 'yes']:
                    # Drop the existing table
                    try:
                        conn = sqlite3.connect(self.__chemspace_db)
                        cursor = conn.cursor()
                        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
                        conn.commit()
                        conn.close()
                        print(f"✅ Table '{table_name}' is to be replaced.")
                    except Exception as e:
                        return {
                            'success': False,
                            'message': f"Failed to replace table '{table_name}': {e}",
                            'compounds_added': 0,
                            'duplicates_skipped': 0,
                            'errors': 1
                        }
                else:
                    print("Quiting .csv file addition")
                    
                    return {
                        'success': False,
                        'message': f"Table '{table_name}' already exists and was not replaced.",
                        'compounds_added': 0,
                        'duplicates_skipped': 0,
                        'errors': 1
                    }

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
            
            # Read the .csv file to be inputed
            df = pd.read_csv(csv_file_path)
            
            # Parse the df to check columns and information
            df = self._parse_df_from_csv_file(df, smiles_column, name_column, flag_column)
            
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
            flag_description_available = 'flag_description' in df.columns

            if not name_available:
                print(f"⚠️  Name column '{name_column}' not found. Will use 'nd' as default.")
            if not flag_available:
                print(f"⚠️  Flag column '{flag_column}' not found. Will use 'nd' as default.")
            if not flag_description_available:
                print(f"⚠️  Flag description column 'flag_description' not found. Will use 'nd' as default.")
            
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
                    name_available, flag_available, #flag_description_available, 
                    skip_duplicates,
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
                            name = f"compound_{idx + 1}"
                    else:
                        name = f"compound_{idx + 1}"
                    # Handle flag column - use provided column or default to "nd"
                    if flag_available:
                        flag = str(row[flag_column]).strip()
                        if not flag or flag.lower() in ['nan', 'none', '']:
                            flag = "nd"
                    else:
                        flag = "nd"
                    # Handle flag_description column - use provided column or default to "nd"
                    if flag_description_available:
                        flag_description = str(row['flag_description']) if 'flag_description' in row else 'nd'
                        if not flag_description or flag_description.lower() in ['nan', 'none', '']:
                            flag_description = "nd"
                    else:
                        flag_description = "nd"
                    # Skip empty SMILES rows (only SMILES is mandatory)
                    if not smiles or smiles.lower() in ['nan', 'none', '']:
                        continue
                    compounds_data.append((smiles, name, flag, flag_description))
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
            error_message = f"Error loading CSV file: {e}"
            print(f"❌ {error_message}")
            return {
                'success': False,
                'message': error_message,
                'compounds_added': 0,
                'duplicates_skipped': 0,
                'errors': 1
            }
    
    def load_single_smiles(self) -> Dict[str, Any]:
        """
        Will prompt the user for a smiles, a name, a flag, and a table name. 
        Then the smiles will be loaded to the chemspace.db using the approach followed by the load_csv_file method
        
        Returns:
            Dict[str, Any]: Dictionary containing operation results with keys:
                - success (bool): Whether the operation was successful
                - message (str): Descriptive message about the operation
                - compounds_added (int): Number of compounds added
                - table_name (str): Name of the table where compound was added
                - smiles (str): The SMILES string that was added
        """
        try:
            print("\n" + "="*70)
            print("LOAD SINGLE SMILES")
            print("="*70)
            
            # Prompt for SMILES string
            print("\n📋 Enter compound information:")
            smiles = input("🧬 Enter SMILES string: ").strip()
            
            if not smiles:
                print("❌ SMILES string cannot be empty.")
                return {
                    'success': False,
                    'message': 'SMILES string cannot be empty',
                    'compounds_added': 0,
                    'table_name': None,
                    'smiles': None
                }
            
            # Validate SMILES using RDKit
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    print(f"❌ Invalid SMILES string: {smiles}")
                    return {
                        'success': False,
                        'message': f'Invalid SMILES string: {smiles}',
                        'compounds_added': 0,
                        'table_name': None,
                        'smiles': smiles
                    }
            except ImportError:
                print("⚠️  Warning: RDKit not available. Skipping SMILES validation.")
            
            # Prompt for compound name
            name = input("📝 Enter compound name: ").strip()
            if not name:
                print("❌ Compound name cannot be empty.")
                return {
                    'success': False,
                    'message': 'Compound name cannot be empty',
                    'compounds_added': 0,
                    'table_name': None,
                    'smiles': smiles
                }
            
            # Prompt for flag
            flag = input("🏷️  Enter flag (optional, press Enter to skip): ").strip()
            if not flag:
                flag = None
            
            # Prompt for flag description
            flag_description = input("📌 Enter flag description (optional): ").strip()
            if not flag_description:
                flag_description = None
            
            # Prompt for table name
            existing_tables = self.get_all_tables()
            
            print("\n📊 Available tables:")
            if existing_tables:
                for idx, table in enumerate(existing_tables, 1):
                    print(f"  [{idx}] {table}")
            else:
                print("  No existing tables found.")
            
            print(f"  [N] Create new table")
            
            table_choice = input("\nSelect table (enter number, table name, or 'N' for new): ").strip()
            
            if table_choice.upper() == 'N':
                # Create new table
                table_name = input("Enter new table name: ").strip()
                if not table_name:
                    print("❌ Table name cannot be empty.")
                    return {
                        'success': False,
                        'message': 'Table name cannot be empty',
                        'compounds_added': 0,
                        'table_name': None,
                        'smiles': smiles
                    }
                
                table_name = self._sanitize_table_name(table_name)
                
                # Create the table
                if not self._create_compounds_table(table_name):
                    print(f"❌ Failed to create table: {table_name}")
                    return {
                        'success': False,
                        'message': f'Failed to create table: {table_name}',
                        'compounds_added': 0,
                        'table_name': table_name,
                        'smiles': smiles
                    }
                
                print(f"✅ Table created: {table_name}")
            
            elif table_choice.isdigit():
                # Select existing table by number
                idx = int(table_choice) - 1
                if idx < 0 or idx >= len(existing_tables):
                    print("❌ Invalid table selection.")
                    return {
                        'success': False,
                        'message': 'Invalid table selection',
                        'compounds_added': 0,
                        'table_name': None,
                        'smiles': smiles
                    }
                table_name = existing_tables[idx]
            
            else:
                # Use the provided table name directly
                table_name = self._sanitize_table_name(table_choice)
                if table_name not in existing_tables:
                    print(f"⚠️  Table '{table_name}' does not exist. Creating new table...")
                    if not self._create_compounds_table(table_name):
                        print(f"❌ Failed to create table: {table_name}")
                        return {
                            'success': False,
                            'message': f'Failed to create table: {table_name}',
                            'compounds_added': 0,
                            'table_name': table_name,
                            'smiles': smiles
                        }
            
            # Show summary
            print("\n" + "="*70)
            print("COMPOUND SUMMARY")
            print("="*70)
            print(f"📊 Table: {table_name}")
            print(f"🧬 SMILES: {smiles}")
            print(f"📝 Name: {name}")
            print(f"🏷️  Flag: {flag or '(none)'}")
            print(f"📌 Flag Description: {flag_description or '(none)'}")
            print("="*70)
            
            # Confirm insertion
            confirm = input("\nProceed with insertion? (y/n): ").strip().lower()
            if confirm not in ['y', 'yes']:
                print("Operation cancelled.")
                return {
                    'success': False,
                    'message': 'Operation cancelled by user',
                    'compounds_added': 0,
                    'table_name': table_name,
                    'smiles': smiles
                }
            
            # Insert the compound
            compound_data = [(smiles, name, flag, flag_description)]
            result = self._insert_compounds(compound_data, table_name, skip_duplicates=True)
            
            if result['success']:
                print(f"\n✅ Compound loaded successfully to table '{table_name}'")
                
                # Compute InChI key for the inserted compound
                print(f"\n🧪 Computing InChI key for the loaded compound...")
                try:
                    inchi_df = self.compute_inchi_keys(table_name, update_database=True)
                    if not inchi_df.empty and 'inchi_key' in inchi_df.columns:
                        match_df = inchi_df[inchi_df['smiles'] == smiles]
                        if not match_df.empty:
                            inchi_key_value = match_df.iloc[0].get('inchi_key')
                            if inchi_key_value and inchi_key_value not in ['INVALID_SMILES', 'ERROR']:
                                print(f"   🔬 InChI key: {inchi_key_value}")
                            else:
                                print("   ⚠️  InChI key not computed for the inserted compound")
                    else:
                        print("   ⚠️  No InChI keys computed")
                except Exception as inchi_error:
                    print(f"   ⚠️  Warning: InChI key computation failed: {inchi_error}")

                return {
                    'success': True,
                    'message': f"Compound loaded to table '{table_name}'",
                    'compounds_added': result['compounds_added'],
                    'table_name': table_name,
                    'smiles': smiles
                }

            

            else:
                print(f"\n❌ Failed to load compound: {result['message']}")
                return {
                    'success': False,
                    'message': result['message'],
                    'compounds_added': 0,
                    'table_name': table_name,
                    'smiles': smiles
                }
        
        except KeyboardInterrupt:
            print("\n\n⚠️  Operation cancelled by user.")
            return {
                'success': False,
                'message': 'Operation cancelled by user',
                'compounds_added': 0,
                'table_name': None,
                'smiles': None
            }
        
        except Exception as e:
            print(f"\n❌ Unexpected error: {e}")
            return {
                'success': False,
                'message': f'Unexpected error: {e}',
                'compounds_added': 0,
                'table_name': None,
                'smiles': None
            }    
    
    def _parse_df_from_csv_file(self, df, smiles_column, name_column, flag_column):
        
        print("Parsing dataframe from .csv file to check columns and information...")
        
        # Parsing the smiles column
        if smiles_column not in df.columns:
            print(f"❌ The specified SMILES column '{smiles_column}' was not found in the CSV file.")
            print("Available columns:")
            for idx, col in enumerate(df.columns):
                print(f"  [{idx}] {col}")
            selected_idx = None
            while selected_idx is None:
                try:
                    user_input = input("Please enter the number of the column to use as SMILES: ").strip()
                    selected_idx = int(user_input)
                    if selected_idx < 0 or selected_idx >= len(df.columns):
                        print("Invalid selection. Please try again.")
                        selected_idx = None
                except ValueError:
                    print("Invalid input. Please enter a valid number.")
            selected_col = df.columns[selected_idx]
            print(f"Renaming column '{selected_col}' to '{smiles_column}'.")
            df.rename(columns={selected_col: smiles_column}, inplace=True)
        
        # Parsing the name column        
        if name_column not in df.columns:
            print(f"❌ The specified NAMES column '{name_column}' was not found in the CSV file.")
            print("Available columns:")
            for idx, col in enumerate(df.columns):
                print(f"  [{idx}] {col}")
            selected_idx = None
            while selected_idx is None:
                try:
                    user_input = input("Please enter the number of the column to use as NAMES: ").strip()
                    selected_idx = int(user_input)
                    if selected_idx < 0 or selected_idx >= len(df.columns):
                        print("Invalid selection. Please try again.")
                        selected_idx = None
                except ValueError:
                    print("Invalid input. Please enter a valid number.")
            selected_col = df.columns[selected_idx]
            print(f"Renaming column '{selected_col}' to '{name_column}'.")
            df.rename(columns={selected_col: name_column}, inplace=True)
        
        # Parsing the flag column
        if flag_column not in df.columns:
            print(f"❌ The specified FLAG column '{flag_column}' was not found in the CSV file.")
            print("Available columns:")
            for idx, col in enumerate(df.columns):
                print(f"  [{idx}] {col}")
            selected_idx = None
            while selected_idx is None:
                try:
                    user_input = input("Please enter the number of the column to use as FLAGS: ").strip()
                    selected_idx = int(user_input)
                    if selected_idx < 0 or selected_idx >= len(df.columns):
                        print("Invalid selection. Please try again.")
                        selected_idx = None
                except ValueError:
                    print("Invalid input. Please enter a valid number.")
            selected_col = df.columns[selected_idx]
            print(f"Renaming column '{selected_col}' to '{flag_column}'.")
            df.rename(columns={selected_col: flag_column}, inplace=True)
        
        # Add 'flag_description' column to the dataframe
        flag_description = input("Enter a description for the 'flag' column to be added to all rows: ").strip()
        df['flag_description'] = flag_description

        # Only keep the required columns
        required_columns = [smiles_column, name_column, flag_column, 'flag_description']
        # If any required column is missing, add it with default value 'nd' (except flag_description)
        for col in [smiles_column, name_column, flag_column]:
            if col not in df.columns:
                df[col] = 'nd'
        df = df[required_columns]
        
        return df
    
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
            INSERT INTO {table_name} (smiles, name, flag, flag_description)
            VALUES (?, ?, ?, ?)
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
    
    def get_sql_registers(self, limit: Optional[int] = None, target_table_filter: Optional[str] = None) -> List[Dict[str, Any]]:
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
    
    def create_table_from_sql(self, sql_query: Optional[str] = None, 
                         source_db_path: Optional[str] = None,
                         clean_files: bool = True,
                         dry_run: bool = False):
        """
        Create a molecules table in the chemspace database from an SQL query executed on an external database.
        
        Args:
            sql_query (Optional[str]): SQL SELECT query to execute on the source database.
                If not provided, will be prompted from user input.
            source_db_path (Optional[str]): Path to the source SQLite database file.
                If not provided, will be prompted from user input.
            clean_files (bool): Whether to remove temporary files after loading. Default is True.
            
        Returns:
            dict: Results containing success status, counts, and messages
        """
    
        try:
            
            # Prompt for sql_query if not provided
            if sql_query is None:
                print("📝 Please enter your SQL SELECT query:")
                print("💡 Example: SELECT smiles, name, flag FROM compounds LIMIT 100")
                print("   (Press Enter twice when done)")
                lines = []
                while True:
                    line = input()
                    if line:
                        lines.append(line)
                    else:
                        if lines:
                            break
                sql_query = " ".join(lines).strip()
                
                if not sql_query:
                    return {
                        'success': False,
                        'message': "No SQL query provided",
                        'compounds_added': 0,
                        'duplicates_skipped': 0,
                        'errors': 1,
                        'inchi_keys_computed': 0
                    }
            
            # Prompt for source_db_path if not provided
            if source_db_path is None:
                user_input = input("🗄️  Enter the path to the source database file (press Enter for default: chemspace.db): ").strip()
                
                if user_input:
                    # User provided a path
                    source_db_path = user_input
                else:
                    # User pressed Enter without input, use default path
                    source_db_path = os.path.join(self.path, 'chemspace', 'processed_data', 'chemspace.db')
                    print(f"📂 Using default source database: {source_db_path}")
            
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
            
            print(f"📊 Executing SQL query on source database: {source_db_path}")
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
            
            
            ## If dry_run is false, store the query the db
            if not dry_run:
            
                ## Write the df to a tempory dataframe for further loading
                temp_csv_file = os.path.join(self.path, 'chemspace/processed_data', 'temp_csv.db')
                
                df.to_csv(temp_csv_file, header=True, index=False)
                
                # Use the standard method to load the generated .csv file
                result = self.load_csv_file(temp_csv_file)
                
                # Register the SQL query in the chemspace.db
                target_table_name = result['table_name']
                self._register_sql_operation(sql_query, source_db_path, target_table_name, result)                
                
                if clean_files: 
                    # After loading the .csv file remove it if requested
                    os.remove(temp_csv_file)
                    
            else:
                # If a dry_run is requested, just print the query df
                print(df)
            
        except Exception as e:
            print(f"Error loading table from sql: {e}")
            
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
                target_table_name = result['table_name']
                self._register_sql_operation(sql_query, source_db_path, target_table_name, error_result)
            except Exception as e:
                print(f"Error registering the table creation process: {e}")
            
            return error_result
    
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

    def filter_using_workflow(self, table_name: Optional[str] = None, workflow_name: Optional[str] = None, 
                    save_results: Optional[bool] = None, 
                    result_table_name: Optional[str] = None,
                    parallel_threshold: int = 10000,
                    max_workers: Optional[int] = None,
                    chunk_size: Optional[int] = None) -> pd.DataFrame:
        """
        Apply a chemical filtering workflow to compounds in a table with parallel processing support.
        Shows available tables and workflows for interactive selection if not provided.
        
        Args:
            table_name (Optional[str]): Name of the table containing compounds to filter. If None, shows selection.
            workflow_name (Optional[str]): Name of the workflow to apply. If None, shows selection.
            save_results (Optional[bool]): Whether to save filtered results to database (prompts if None)
            result_table_name (Optional[str]): Name for the result table (prompts if save_results=True and None)
            parallel_threshold (int): Minimum number of compounds to trigger parallel processing
            max_workers (Optional[int]): Maximum number of worker processes (default: min(cpu_count(), 8))
            chunk_size (Optional[int]): Size of chunks for parallel processing (auto-calculated if None)
            
        Returns:
            pd.DataFrame: DataFrame containing filtered compounds with match information
        """
        try:
            print(f"\n🔬 Starting workflow filtering...")
            
            # Interactive table selection if not provided
            if table_name is None:
                table_name = self._select_table_for_filtering()
                if not table_name:
                    print("❌ No table selected for filtering")
                    return pd.DataFrame()
            
            # Interactive workflow selection if not provided
            if workflow_name is None:
                workflow_name = self._select_workflow_for_filtering()
                if not workflow_name:
                    print("❌ No workflow selected for filtering")
                    return pd.DataFrame()
            
            print(f"   📋 Table: '{table_name}'")
            print(f"   🧪 Workflow: '{workflow_name}'")
            
            # Load workflow filters with progress
            print(f"   🔍 Loading workflow filters...")
            workflow_filters, filters_dict = self._load_workflow_filters(workflow_name)
            
            if not workflow_filters:
                print(f"❌ No filters found for workflow '{workflow_name}'")
                return pd.DataFrame()
            
            print(f"   ✅ Loaded {len(workflow_filters)} filters")
            
            # Print filters showing their corresponding key from filters_dict
            keys = list(filters_dict.keys())
            for i, (filter_smarts, _) in enumerate(workflow_filters, 1):
                key = keys[i-1] if i-1 < len(keys) else f"unknown_key_{i}"
                print(f"      {i}. {key} -> {filter_smarts}")
            
            # Get compounds from table with progress
            print(f"   📊 Loading compounds from table '{table_name}'...")
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return pd.DataFrame()
            
            total_compounds = len(compounds_df)
            print(f"   ✅ Loaded {total_compounds:,} compounds to filter")
            
            # Determine processing method
            use_parallel = total_compounds >= parallel_threshold
            
            if use_parallel:
                print(f"   🚀 Using parallel processing (threshold: {parallel_threshold:,})")
                
                # Set default parameters for parallel processing
                if max_workers is None:
                    max_workers = min(os.cpu_count() or 4, 8)
                
                if chunk_size is None:
                    chunk_size = max(1000, total_compounds // (max_workers * 4))
                
                print(f"      👥 Workers: {max_workers}")
                print(f"      📦 Chunk size: {chunk_size:,}")
            
                filtered_df = self._apply_filters_parallel(
                    compounds_df, workflow_filters, max_workers, chunk_size
                )
            else:
                print(f"   🔄 Using sequential processing")
                filtered_df = self._apply_filters_sequential(compounds_df, workflow_filters)
            
            if filtered_df.empty:
                print("❌ No compounds passed the filtering workflow")
                return pd.DataFrame()
            
            # Calculate statistics
            compounds_removed = total_compounds - len(filtered_df)
            retention_rate = (len(filtered_df) / total_compounds) * 100
            
            print(f"\n📊 Filtering Results:")
            print(f"   ✅ Compounds passed: {len(filtered_df):,}")
            print(f"   ❌ Compounds removed: {compounds_removed:,}")
            print(f"   📈 Retention rate: {retention_rate:.2f}%")
            
            # Prompt for save_results if not provided
            if save_results is None:
                save_choice = input("\n💾 Do you want to save the filtered results to a new table? (y/n): ").strip().lower()
                save_results = save_choice in ['y', 'yes']
            
            # Save results if requested with progress
            if save_results:
                # Prompt for table name if not provided
                if result_table_name is None:
                    default_name = f"{table_name}_filtered_{workflow_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                    user_table_name = input(f"📝 Enter table name for results (default: {default_name}): ").strip()
                    result_table_name = user_table_name if user_table_name else default_name
                
                print(f"   💾 Saving filtered results to '{result_table_name}'...")
                
                # Create filter results summary for metadata
                filter_results = {
                    'workflow_name': workflow_name,
                    'initial_compounds': total_compounds,
                    'final_compounds': len(filtered_df),
                    'compounds_removed': compounds_removed,
                    'retention_rate': retention_rate,
                    'processing_method': 'parallel' if use_parallel else 'sequential'
                }
                
                success = self._save_workflow_filtered_compounds(
                    filtered_df, result_table_name, workflow_name, filter_results
                )
                
                if success:
                    print(f"   ✅ Results saved to table: '{result_table_name}'")
                else:
                    print(f"   ❌ Failed to save results to database")
            else:
                print("   📄 Results not saved to database")
            
            return filtered_df
            
        except Exception as e:
            print(f"❌ Error in filter_using_workflow: {e}")
            return pd.DataFrame()

    def _select_table_for_filtering(self) -> Optional[str]:
        """
        Interactive selection of table for filtering.
        
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        try:
            available_tables = self.get_all_tables()
            
            if not available_tables:
                print("❌ No tables available for filtering")
                return None
            
            print(f"\n📋 SELECT TABLE FOR FILTERING")
            print("=" * 60)
            print(f"Available tables ({len(available_tables)} total):")
            print("-" * 60)
            
            # Display tables without compound counts for better performance
            table_info = []
            for i, table_name in enumerate(available_tables, 1):
                try:
                    table_type = self._classify_table_type(table_name)
                    type_icon = self._get_table_type_icon(table_type)
                    
                    print(f"{i:3d}. {type_icon} {table_name:<35} [{table_type}]")
                    table_info.append({
                        'name': table_name,
                        'type': table_type
                    })
                    
                except Exception as e:
                    print(f"{i:3d}. ❓ {table_name:<35} [Error: {e}]")
                    table_info.append({
                        'name': table_name,
                        'type': 'unknown'
                    })
            
            print("-" * 60)
            print("Commands: Enter table number, table name, or 'cancel' to abort")
            
            while True:
                try:
                    selection = input(f"\n🔍 Select table for filtering: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    
                    # Try as number first
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(available_tables):
                            selected_table = available_tables[table_idx]
                            selected_info = table_info[table_idx]
                            
                            print(f"\n✅ Selected table: '{selected_table}'")
                            print(f"   🏷️  Type: {selected_info['type']}")
                            return selected_table
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(available_tables)}")
                            continue
                    except ValueError:
                        # Try as table name
                        matching_tables = [t for t in available_tables if t.lower() == selection.lower()]
                        if matching_tables:
                            selected_table = matching_tables[0]
                            # Find table info
                            selected_info = next((info for info in table_info if info['name'] == selected_table), 
                                            {'type': 'unknown'})
                            
                            print(f"\n✅ Selected table: '{selected_table}'")
                            print(f"   🏷️  Type: {selected_info['type']}")
                            return selected_table
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Table selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error selecting table for filtering: {e}")
            return None

    def _select_workflow_for_filtering(self) -> Optional[str]:
        """
        Interactive selection of workflow for filtering.
        
        Returns:
            Optional[str]: Selected workflow name or None if cancelled
        """
        try:
            # Get available filtering workflows
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Check if filtering_workflows table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='filtering_workflows'")
            if not cursor.fetchone():
                print("❌ No filtering workflows table found")
                print("   Create workflows first using create_filtering_workflow()")
                conn.close()
                return None
            
            # Get all workflows
            cursor.execute("""
                SELECT workflow_name, creation_date, description, filter_count, total_instances
                FROM filtering_workflows 
                ORDER BY creation_date DESC
            """)
            workflows = cursor.fetchall()
            conn.close()
            
            if not workflows:
                print("❌ No filtering workflows found")
                print("   Create workflows first using create_filtering_workflow()")
                return None
            
            print(f"\n🧪 SELECT FILTERING WORKFLOW")
            print("=" * 80)
            print(f"Available workflows ({len(workflows)} total):")
            print("-" * 80)
            print(f"{'#':<3} {'Workflow Name':<25} {'Filters':<8} {'Instances':<10} {'Created':<12} {'Description':<50}")
            print("-" * 80)
            
            workflow_info = []
            for i, (name, date, desc, filter_count, total_instances) in enumerate(workflows, 1):
                date_short = date[:10] if date else "Unknown"
                name_display = name[:24] if len(name) <= 24 else name[:21] + "..."
                desc_display = (desc or "No description")[:50]
                if len(desc_display) > 50:
                    desc_display = desc_display[:50] + "..."
                
                print(f"{i:<3} {name_display:<25} {filter_count or 0:<8} {total_instances or 0:<10} "
                    f"{date_short:<12} {desc_display:<25}")
                workflow_info.append({
                    'name': name,
                    'filter_count': filter_count or 0,
                    'total_instances': total_instances or 0,
                    'description': desc or "No description",
                    'date': date
                })
            
            print("-" * 80)
            print("Commands: Enter workflow number, workflow name, 'list' for details, or 'cancel' to abort")
            
            while True:
                try:
                    selection = input(f"\n🔍 Select filtering workflow: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    elif selection.lower() == 'list':
                        self._show_detailed_workflow_info(workflow_info)
                        continue
                    
                    # Try as number first
                    try:
                        workflow_idx = int(selection) - 1
                        if 0 <= workflow_idx < len(workflows):
                            selected_workflow = workflows[workflow_idx][0]  # workflow name
                            selected_info = workflow_info[workflow_idx]
                            
                            print(f"\n✅ Selected workflow: '{selected_workflow}'")
                            print(f"   🔍 Filters: {selected_info['filter_count']}")
                            print(f"   🔢 Total instances: {selected_info['total_instances']}")
                            print(f"   📄 Description: {selected_info['description']}")
                            return selected_workflow
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(workflows)}")
                            continue
                    except ValueError:
                        # Try as workflow name
                        matching_workflows = [w[0] for w in workflows if w[0].lower() == selection.lower()]
                        if matching_workflows:
                            selected_workflow = matching_workflows[0]
                            # Find workflow info
                            selected_info = next((info for info in workflow_info if info['name'] == selected_workflow), 
                                            {'filter_count': 0, 'total_instances': 0, 'description': 'Unknown'})
                            
                            print(f"\n✅ Selected workflow: '{selected_workflow}'")
                            print(f"   🔍 Filters: {selected_info['filter_count']}")
                            print(f"   🔢 Total instances: {selected_info['total_instances']}")
                            print(f"   📄 Description: {selected_info['description']}")
                            return selected_workflow
                        else:
                            print(f"❌ Workflow '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Workflow selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error selecting workflow for filtering: {e}")
            return None

    def _show_detailed_workflow_info(self, workflow_info: List[Dict]) -> None:
        """
        Show detailed information about filtering workflows.
        
        Args:
            workflow_info (List[Dict]): List of workflow information dictionaries
        """
        try:
            print(f"\n📋 DETAILED WORKFLOW INFORMATION")
            print("=" * 80)
            
            for i, info in enumerate(workflow_info, 1):
                print(f"\n{i}. {info['name']}")
                print(f"   📅 Created: {info['date'] or 'Unknown'}")
                print(f"   🔍 Filters: {info['filter_count']}")
                print(f"   🔢 Total instances: {info['total_instances']}")
                print(f"   📄 Description: {info['description']}")
            
            print("=" * 80)
            
        except Exception as e:
            print(f"❌ Error showing detailed workflow information: {e}")

    def _load_workflow_filters(self, workflow_name: str) -> List[Tuple[str, str]]:
        """
        Load workflow filters from the database and retrieve their SMARTS patterns.
        
        Args:
            workflow_name (str): Name of the workflow to load
            
        Returns:
            List[Tuple[str, str]]: List of (filter_name, smarts_pattern) tuples, 
                                empty list if workflow not found or error occurs
        """
        try:
            # Check if the table exists in chemspace database
            tables = self.get_all_tables()
            if 'filtering_workflows' not in [table.lower() for table in tables]:
                print("❌ No filtering workflows table found in chemspace database")
                print("   Create workflows first using create_filtering_workflow()")
                return []
            
            # Connect to chemspace database to retrieve the workflow
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Retrieve the workflow by name
            cursor.execute(
                "SELECT filters_dict FROM filtering_workflows WHERE workflow_name = ?", 
                (workflow_name,)
            )
            workflow_result = cursor.fetchone()
            conn.close()
            
            if not workflow_result:
                print(f"❌ Workflow '{workflow_name}' not found in filtering_workflows table")
                print("\n📋 Available workflows:")
                self._show_available_workflows()
                return []
            
            # Parse the filters dictionary from JSON string
            try:
                filters_dict = json.loads(workflow_result[0])
            
            except json.JSONDecodeError as e:
                print(f"❌ Error parsing filters dictionary for workflow '{workflow_name}': {e}")
                return []
            
            if not filters_dict:
                print(f"⚠️  Workflow '{workflow_name}' contains no filters")
                return []
            
            # Get the projects database path to retrieve SMARTS patterns
            import site
            projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"
            
            if not os.path.exists(projects_db):
                print(f"❌ Projects database not found: {projects_db}")
                return []
            
            # Connect to projects database and retrieve SMARTS patterns
            projects_conn = sqlite3.connect(projects_db)
            projects_cursor = projects_conn.cursor()
            
            workflow_filters = []
            missing_filters = []
            
            for filter_name, required_instances in filters_dict.items():
                # Get SMARTS pattern for this filter name
                projects_cursor.execute(
                    "SELECT smarts FROM chem_filters WHERE filter_name = ?", 
                    (filter_name,)
                )
                smarts_result = projects_cursor.fetchone()
                
                if smarts_result:
                    smarts_pattern = smarts_result[0]
                    #workflow_filters.append((filter_name, smarts_pattern))
                    workflow_filters.append((smarts_pattern, required_instances))
                else:
                    missing_filters.append(filter_name)
            
            projects_conn.close()
            
            # Report any missing filters
            if missing_filters:
                print(f"⚠️  Warning: {len(missing_filters)} filter(s) not found in projects database:")
                for filter_name in missing_filters:
                    print(f"   - '{filter_name}'")
                print("   These filters will be skipped during workflow execution")
            
            if not workflow_filters:
                print(f"❌ No valid filters found for workflow '{workflow_name}'")
                return []
            
            print(f"✅ Loaded {len(workflow_filters)} valid filters for workflow '{workflow_name}'")
            return workflow_filters, filters_dict
            
        except Exception as e:
            print(f"❌ Error loading workflow filters for '{workflow_name}': {e}")
            return []

    def _apply_filters_parallel(self, compounds_df: pd.DataFrame, 
                                    workflow_filters: List[Tuple[str, str]], 
                                    max_workers: int, chunk_size: int) -> pd.DataFrame:
        """
        Ultra-memory-efficient version that streams results to temporary files with progress tracking.
        """
        import tempfile
        import os
        
        # Import tqdm locally to avoid multiprocessing issues
        try:
            from tqdm import tqdm
            tqdm_available = True
        except ImportError:
            tqdm_available = False
        
        try:
            # Create temporary file for results
            temp_dir = tempfile.mkdtemp(prefix='chemspace_filter_')
            temp_file = os.path.join(temp_dir, 'filtered_results.csv')
            
            # Write header
            with open(temp_file, 'w') as f:
                f.write('id,smiles,name,flag,inchi_key\n')
            
            # Process in batches and append to file
            compound_data = [(row.get('id', 0), row['smiles'], row.get('name', 'unknown'),
                            row.get('flag', 'nd'), row.get('inchi_key', None))
                            for _, row in compounds_df.iterrows()]
            
            chunks = [compound_data[i:i + chunk_size] 
                    for i in range(0, len(compound_data), chunk_size)]
            
            total_passed = 0
            processed_chunks = 0
            
            print(f"   📦 Processing {len(chunks)} chunks with {max_workers} workers...")
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Process chunks and write results immediately
                future_to_chunk = {
                    executor.submit(_filter_chunk_worker_by_instances, chunk, workflow_filters): i 
                    for i, chunk in enumerate(chunks[:max_workers])  # Start with first batch
                }
                
                remaining_chunks = chunks[max_workers:]
                
                # Initialize progress bar with local tqdm check
                if tqdm_available:
                    progress_bar = tqdm(
                        total=len(chunks),
                        desc="Filtering chunks",
                        unit="chunks",
                        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
                    )
                else:
                    progress_bar = None
                    print(f"   🔄 Started processing {len(chunks)} chunks...")
                
                with open(temp_file, 'a') as f:
                    while future_to_chunk:
                        for future in as_completed(list(future_to_chunk.keys())):
                            chunk_idx = future_to_chunk.pop(future)
                            
                            try:
                                chunk_results = future.result()
                                
                                # Write results immediately to disk
                                chunk_passed = 0
                                for result in chunk_results:
                                    f.write(f"{result['id']},{result['smiles']},{result['name']},"
                                        f"{result['flag']},{result.get('inchi_key', '')}\n")
                                    chunk_passed += 1
                                
                                total_passed += chunk_passed
                                processed_chunks += 1
                                
                                # Update progress
                                if progress_bar:
                                    progress_bar.update(1)
                                    progress_bar.set_postfix({
                                        'passed': total_passed,
                                        'chunks': f"{processed_chunks}/{len(chunks)}"
                                    })
                                else:
                                    # Manual progress for non-tqdm
                                    if processed_chunks % max(1, len(chunks) // 10) == 0:
                                        progress_pct = (processed_chunks / len(chunks)) * 100
                                        print(f"      📊 Progress: {progress_pct:.1f}% ({processed_chunks}/{len(chunks)} chunks, {total_passed} passed)")
                                
                                # Submit next chunk if available
                                if remaining_chunks:
                                    next_chunk = remaining_chunks.pop(0)
                                    new_future = executor.submit(
                                        _filter_chunk_worker_by_instances, next_chunk, workflow_filters
                                    )
                                    future_to_chunk[new_future] = len(chunks) - len(remaining_chunks) - 1
                                    
                            except Exception as e:
                                processed_chunks += 1
                                if progress_bar:
                                    progress_bar.update(1)
                                print(f"   ❌ Chunk {chunk_idx} error: {e}")
                
                if progress_bar:
                    progress_bar.close()
            
            # Read results back from file with progress
            print(f"   📊 Reading {total_passed} filtered compounds from disk...")
            try:
                if tqdm_available:
                    # Use tqdm for reading large files with local import
                    from tqdm import tqdm
                    import pandas as pd
                    
                    # Get file size for progress bar
                    file_size = os.path.getsize(temp_file)
                    
                    with tqdm(total=file_size, unit='B', unit_scale=True, desc="Reading results") as pbar:
                        result_df = pd.read_csv(temp_file)
                        pbar.update(file_size)  # Update to completion
                        
                else:
                    result_df = pd.read_csv(temp_file)
                
                print(f"   ✅ Successfully loaded {len(result_df)} filtered compounds")
                return result_df
            finally:
                # Clean up temporary files
                try:
                    os.unlink(temp_file)
                    os.rmdir(temp_dir)
                except:
                    pass
                    
        except Exception as e:
            print(f"❌ Error in streaming parallel filtering: {e}")
            return pd.DataFrame()
   
    def _apply_filters_sequential(self, compounds_df: pd.DataFrame, 
                                workflow_filters: List[Tuple[str, str]]) -> pd.DataFrame:
        """
        Apply SMARTS filters using sequential processing with enhanced progress tracking.
        """
        # Import tqdm locally
        try:
            from tqdm import tqdm
            tqdm_available = True
        except ImportError:
            tqdm_available = False
        
        try:
            from rdkit import Chem
            from rdkit import RDLogger
            RDLogger.DisableLog('rdApp.*')  # Suppress RDKit warnings
            
            filtered_compounds = []
            total_compounds = len(compounds_df)
            processed_count = 0
            compounds_passed = 0
            
            print(f"   🔄 Processing {total_compounds:,} compounds sequentially...")
            
            # Initialize progress tracking with enhanced information
            if tqdm_available:
                progress_bar = tqdm(
                    compounds_df.iterrows(),
                    total=total_compounds,
                    desc="Filtering compounds",
                    unit="compounds",
                    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] Passed: {postfix}",
                    postfix="0"
                )
            else:
                progress_bar = compounds_df.iterrows()
                print_interval = max(1, total_compounds // 20)  # Print progress 20 times
            
            for _, compound in progress_bar:
                processed_count += 1
                passes_all_filters = True
                match_counts = {}
                
                try:
                    mol = Chem.MolFromSmiles(compound['smiles'])
                    if mol is None:
                        continue
                    
                    # Apply each filter
                    for filter_name, smarts_pattern in workflow_filters:
                        try:
                            pattern = Chem.MolFromSmarts(smarts_pattern)
                            if pattern is None:
                                match_counts[f"{filter_name}_matches"] = 0
                                continue
                            
                            matches = mol.GetSubstructMatches(pattern)
                            match_count = len(matches)
                            match_counts[f"{filter_name}_matches"] = match_count
                            
                            # If any filter has matches, compound fails (exclusion filters)
                            if match_count > 0:
                                passes_all_filters = False
                                break
                        
                        except Exception:
                            match_counts[f"{filter_name}_matches"] = 0
                    
                    if passes_all_filters:
                        compound_data = {
                            'smiles': compound['smiles'],
                            'name': compound.get('name', 'unknown'),
                            'flag': compound.get('flag', 'nd'),
                            'inchi_key': compound.get('inchi_key', None)
                        }
                        compound_data.update(match_counts)
                        filtered_compounds.append(compound_data)
                        compounds_passed += 1
                
                except Exception:
                    continue
                
                # Update progress display
                if tqdm_available:
                    progress_bar.set_postfix(str(compounds_passed))
                else:
                    if processed_count % print_interval == 0 or processed_count == total_compounds:
                        progress_pct = (processed_count / total_compounds) * 100
                        retention_rate = (compounds_passed / processed_count) * 100 if processed_count > 0 else 0
                        print(f"      📊 Progress: {progress_pct:.1f}% ({processed_count:,}/{total_compounds:,}) | "
                            f"Passed: {compounds_passed:,} ({retention_rate:.1f}%)")
            
            if tqdm_available and hasattr(progress_bar, 'close'):
                progress_bar.close()
            
            final_retention_rate = (len(filtered_compounds) / total_compounds) * 100 if total_compounds > 0 else 0
            print(f"   ✅ Sequential processing completed")
            print(f"   📊 Final results: {len(filtered_compounds):,}/{total_compounds:,} compounds passed ({final_retention_rate:.2f}%)")
            
            return pd.DataFrame(filtered_compounds)
            
        except Exception as e:
            print(f"❌ Error in sequential filtering: {e}")
            return pd.DataFrame()

    def _save_workflow_filtered_compounds(self, compounds_df: pd.DataFrame, 
                                    new_table_name: str, 
                                    workflow_name: str, 
                                    filter_results: Dict[str, Dict]) -> bool:
        """
        Save filtered compounds from workflow to a new table in chemspace.db with metadata and progress tracking.
        
        Args:
            compounds_df (pd.DataFrame): DataFrame containing filtered compounds
            new_table_name (str): Name of the new table to create
            workflow_name (str): Name of the workflow that was applied
            filter_results (Dict[str, Dict]): Results from each filter step
            
        Returns:
            bool: True if compounds were saved successfully, False otherwise
        """
        # Import tqdm locally to avoid scope issues
        try:
            from tqdm import tqdm
            tqdm_available = True
        except ImportError:
            tqdm_available = False
        
        try:
            if compounds_df.empty:
                print("⚠️  No compounds to save from the filtering workflow.")
                return False
            
            print(f"   🔍 Checking for duplicate SMILES among {len(compounds_df)} compounds...")
        
            # Remove duplicates based on SMILES strings
            original_count = len(compounds_df)
            compounds_df_unique = compounds_df.drop_duplicates(subset=['smiles'], keep='first')
            duplicates_removed = original_count - len(compounds_df_unique)
            
            if duplicates_removed > 0:
                print(f"   🔄 Removed {duplicates_removed} duplicate SMILES from filtered results")
                print(f"   ✨ Proceeding with {len(compounds_df_unique)} unique compounds")
            else:
                print(f"   ✅ No duplicate SMILES found in filtered results")
            
            # Use the deduplicated DataFrame for the rest of the process
            compounds_df = compounds_df_unique
            
            if compounds_df.empty:
                print("   ⚠️  No unique compounds remaining after duplicate removal")
                return False
            
            # Sanitize table name
            new_table_name = self._sanitize_table_name(new_table_name)
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create the new table with extended schema including workflow metadata
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {new_table_name} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT NOT NULL,
                    name TEXT,
                    flag TEXT,
                    inchi_key TEXT,
                    UNIQUE(smiles)
                )
            ''')
            
            # Create indexes for better performance
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_table_name}_smiles ON {new_table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_table_name}_name ON {new_table_name}(name)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{new_table_name}_inchi_key ON {new_table_name}(inchi_key)")
            
            # Get current timestamp
            filter_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # Insert the filtered compounds with progress tracking
            inserted_count = 0
            database_duplicates = 0
            errors = 0
            
            insert_query = f'''
                INSERT INTO {new_table_name} 
                (smiles, name, flag, inchi_key)
                VALUES (?, ?, ?, ?)
            '''
            
            print(f"   💾 Saving {len(compounds_df)} filtered compounds to table '{new_table_name}'...")
            
            # Initialize progress bar for saving with proper error handling
            use_tqdm = tqdm_available and len(compounds_df) > 1000
            
            if use_tqdm:
                try:
                    progress_bar = tqdm(
                        compounds_df.iterrows(),
                        total=len(compounds_df),
                        desc="Saving compounds",
                        unit="compounds",
                        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] Saved: {postfix}",
                        postfix="0"
                    )
                except Exception as e:
                    print(f"   ⚠️  Progress bar initialization failed, using fallback: {e}")
                    use_tqdm = False
                    progress_bar = compounds_df.iterrows()
            else:
                progress_bar = compounds_df.iterrows()
            
            # Set up manual progress tracking for non-tqdm case
            if not use_tqdm:
                save_interval = max(1, len(compounds_df) // 10)
                saved_count = 0
            
            # Process compounds
            try:
                for row_data in progress_bar:
                    # Handle the unpacking based on whether we're using tqdm or not
                    if use_tqdm:
                        # tqdm returns (index, row) tuples
                        try:
                            _, compound = row_data
                        except (ValueError, TypeError):
                            # Fallback if unpacking fails
                            compound = row_data
                    else:
                        # Direct iterrows() returns (index, row) tuples
                        try:
                            _, compound = row_data
                        except (ValueError, TypeError):
                            # Fallback if unpacking fails
                            compound = row_data
                    
                    try:
                        cursor.execute(insert_query, (
                            compound.get('smiles', ''),
                            compound.get('name', 'unknown'),
                            compound.get('flag', 'nd'),
                            compound.get('inchi_key', None)
                        ))
                        inserted_count += 1
                        
                        if use_tqdm:
                            try:
                                progress_bar.set_postfix(str(inserted_count))
                            except:
                                pass  # Ignore progress bar update errors
                        else:
                            saved_count += 1
                            if saved_count % save_interval == 0:
                                progress_pct = (saved_count / len(compounds_df)) * 100
                                print(f"      💾 Saving progress: {progress_pct:.1f}% ({saved_count}/{len(compounds_df)})")
                                
                    except sqlite3.IntegrityError:
                        # Handle duplicate SMILES
                        database_duplicates += 1
                    except Exception as e:
                        errors += 1
                        if errors <= 5:  # Show only first few errors
                            print(f"      ⚠️  Error inserting compound '{compound.get('name', 'unknown')}': {e}")
            
            finally:
                # Close progress bar if it was created successfully
                if use_tqdm and hasattr(progress_bar, 'close'):
                    try:
                        progress_bar.close()
                    except:
                        pass  # Ignore close errors
            
            conn.commit()
            conn.close()
            
            # Print summary
            print(f"   ✅ Successfully saved filtered compounds!")
            print(f"      📋 Table name: '{new_table_name}'")
            print(f"      📊 Compounds inserted: {inserted_count}")
            print(f"      🔄 Duplicates skipped: {database_duplicates}")
            print(f"      ❌ Errors: {errors}")
            print(f"      🏷️  Workflow: '{workflow_name}'")
            print(f"      📅 Filter date: {filter_date}")
                    
            return True
            
        except Exception as e:
            print(f"❌ Error saving workflow filtered compounds to table '{new_table_name}': {e}")
            return False

    def depict_table(self, table_name: Optional[str] = None, 
                    output_dir: Optional[str] = None,
                    image_size: Tuple[int, int] = (300, 300),
                    molecules_per_image: int = 25,
                    max_molecules: Optional[int] = None,
                    depict_random: bool = False,
                    highlight_substructures: bool = False,
                    smarts_pattern: Optional[str] = None) -> bool:
        """
        Generate chemical structure depictions from SMILES strings in a table and save as PNG images.
        
        Args:
            table_name (Optional[str]): Name of the table containing compounds with SMILES. If None, prompts user.
            output_dir (Optional[str]): Directory to save images. If None, creates directory named after table
            image_size (Tuple[int, int]): Size of each molecule image (width, height)
            molecules_per_image (int): Number of molecules per image (1 for individual, >1 for grids)
            max_molecules (Optional[int]): Maximum number of molecules to depict. If None, prompts user.
            depict_random (bool): If True and max_molecules is set, randomly select molecules instead of taking first N
            highlight_substructures (bool): Whether to highlight substructures matching a SMARTS pattern
            smarts_pattern (Optional[str]): SMARTS pattern for substructure highlighting
            
        Returns:
            bool: True if depiction was successful
        """
        try:
            # Import RDKit for molecule depiction
            try:
                from rdkit import Chem
                from rdkit.Chem import Draw
                from rdkit.Chem.Draw import rdMolDraw2D
                import PIL
                from PIL import Image, ImageDraw, ImageFont
            except ImportError:
                print("❌ RDKit and/or PIL not installed. Please install them:")
                print("   conda install -c conda-forge rdkit pillow")
                print("   or")
                print("   pip install rdkit pillow")
                return False
            
            # Interactive table selection if not provided
            if table_name is None:
                table_name = self._select_table_for_depiction()
                if not table_name:
                    print("❌ No table selected for depiction")
                    return False
            
            # Check if table exists
            tables = self.get_all_tables()
            if table_name not in tables:
                print(f"❌ Table '{table_name}' does not exist in chemspace database")
                return False
            
            # Get compounds from table
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return False
            
            # Check if SMILES column exists
            if 'smiles' not in compounds_df.columns:
                print("❌ No 'smiles' column found in the table")
                return False
            
            total_available = len(compounds_df)
            print(f"📊 Found {total_available:,} molecules in table '{table_name}'")
            
            # Interactive molecule count selection if not provided
            if max_molecules is None:
                max_molecules = self._prompt_for_molecule_count(total_available)
                if max_molecules is None:
                    print("❌ Depiction cancelled by user")
                    return False
            
            # Handle special case: -1 means all molecules
            if max_molecules == -1:
                max_molecules = total_available
                print(f"🎨 Depicting all {total_available:,} molecules in the table")
            elif max_molecules > 0:
                if depict_random:
                    compounds_df = compounds_df.sample(n=min(max_molecules, total_available))
                    print(f"🔍 Randomly selected {min(max_molecules, total_available)} molecules for depiction")
                else:
                    compounds_df = compounds_df.head(max_molecules)
                    print(f"🔍 Limited to first {min(max_molecules, total_available)} molecules")
            
            total_molecules = len(compounds_df)
            print(f"🎨 Generating depictions for {total_molecules:,} molecules from table '{table_name}'")
            
            # Set up output directory
            if output_dir is None:
                output_dir = os.path.join(self.path, 'chemspace', 'depictions', table_name)
            
            os.makedirs(output_dir, exist_ok=True)
            print(f"📁 Output directory: {output_dir}")
            
            # Parse SMARTS pattern for highlighting if provided
            highlight_mol = None
            if highlight_substructures and smarts_pattern:
                try:
                    highlight_mol = Chem.MolFromSmarts(smarts_pattern)
                    if highlight_mol:
                        print(f"🎯 Highlighting substructures matching: {smarts_pattern}")
                    else:
                        print(f"⚠️  Invalid SMARTS pattern: {smarts_pattern}")
                        highlight_substructures = False
                except:
                    print(f"⚠️  Error parsing SMARTS pattern: {smarts_pattern}")
                    highlight_substructures = False
            
            # Generate depictions
            successful_depictions = 0
            failed_depictions = 0
            invalid_smiles = 0
            
            if molecules_per_image == 1:
                # Individual molecule images
                successful_depictions, failed_depictions, invalid_smiles = self._depict_individual_molecules(
                    compounds_df, output_dir, image_size, highlight_mol, highlight_substructures
                )
            else:
                # Grid images with multiple molecules
                successful_depictions, failed_depictions, invalid_smiles = self._depict_molecule_grids(
                    compounds_df, output_dir, image_size, molecules_per_image, highlight_mol, highlight_substructures
                )
            
            # Generate summary report
            self._generate_depiction_report(
                table_name, output_dir, total_molecules, successful_depictions, 
                failed_depictions, invalid_smiles, image_size, molecules_per_image
            )
            
            print(f"\n✅ Depiction completed!")
            print(f"   🎨 Successful depictions: {successful_depictions}")
            print(f"   ❌ Failed depictions: {failed_depictions}")
            print(f"   ⚠️  Invalid SMILES: {invalid_smiles}")
            print(f"   📁 Images saved in: {output_dir}")
            
            return successful_depictions > 0
            
        except Exception as e:
            print(f"❌ Error in depict_table: {e}")
            return False

    def _select_table_for_depiction(self) -> Optional[str]:
        """
        Interactive selection of table for molecular depiction.
        
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        try:
            available_tables = self.get_all_tables()
            
            if not available_tables:
                print("❌ No tables available for depiction")
                return None
            
            print(f"\n🎨 SELECT TABLE FOR MOLECULAR DEPICTION")
            print("=" * 70)
            print(f"Available tables ({len(available_tables)} total):")
            print("-" * 70)
            
            # Display tables with compound counts and types
            table_info = []
            for i, table_name in enumerate(available_tables, 1):
                try:
                    compound_count = self.get_compound_count(table_name=table_name)
                    table_type = self._classify_table_type(table_name)
                    type_icon = self._get_table_type_icon(table_type)
                    
                    print(f"{i:3d}. {type_icon} {table_name:<30} ({compound_count:>6,} compounds) [{table_type}]")
                    table_info.append({
                        'name': table_name,
                        'type': table_type,
                        'count': compound_count
                    })
                    
                except Exception as e:
                    print(f"{i:3d}. ❓ {table_name:<30} (Error: {e}) [Unknown]")
                    table_info.append({
                        'name': table_name,
                        'type': 'unknown',
                        'count': 0
                    })
            
            print("-" * 70)
            print("Commands: Enter table number, table name, or 'cancel' to abort")
            
            while True:
                try:
                    selection = input(f"\n🔍 Select table for depiction: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    
                    # Try as number first
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(available_tables):
                            selected_table = available_tables[table_idx]
                            selected_info = table_info[table_idx]
                            
                            print(f"\n✅ Selected table: '{selected_table}'")
                            print(f"   🏷️  Type: {selected_info['type']}")
                            print(f"   📊 Compounds: {selected_info['count']:,}")
                            return selected_table
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(available_tables)}")
                            continue
                    except ValueError:
                        # Try as table name
                        matching_tables = [t for t in available_tables if t.lower() == selection.lower()]
                        if matching_tables:
                            selected_table = matching_tables[0]
                            # Find table info
                            selected_info = next((info for info in table_info if info['name'] == selected_table), 
                                            {'type': 'unknown', 'count': 0})
                            
                            print(f"\n✅ Selected table: '{selected_table}'")
                            print(f"   🏷️  Type: {selected_info['type']}")
                            print(f"   📊 Compounds: {selected_info['count']:,}")
                            return selected_table
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Table selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error selecting table for depiction: {e}")
            return None

    def _prompt_for_molecule_count(self, total_available: int) -> Optional[int]:
        """
        Interactive prompt for the number of molecules to depict.
        
        Args:
            total_available (int): Total number of molecules available in the table
            
        Returns:
            Optional[int]: Number of molecules to depict, or None if cancelled
        """
        try:
            print(f"\n📊 MOLECULE COUNT SELECTION")
            print("=" * 40)
            print(f"📋 Total molecules available: {total_available:,}")
            print("\nOptions:")
            print(f"   • Enter a number (1-{total_available:,}) to limit molecules")
            print("   • Enter -1 to depict ALL molecules")
            print("   • Enter 0 or 'cancel' to abort")
            
            # Suggest reasonable defaults based on table size
            if total_available <= 25:
                suggested = total_available
                print(f"\n💡 Suggestion: {suggested} (all molecules)")
            elif total_available <= 100:
                suggested = 50
                print(f"\n💡 Suggestion: {suggested} (good sample size)")
            elif total_available <= 1000:
                suggested = 100
                print(f"\n💡 Suggestion: {suggested} (manageable sample)")
            else:
                suggested = 250
                print(f"\n💡 Suggestion: {suggested} (representative sample)")
            
            print("=" * 40)
            
            while True:
                try:
                    user_input = input(f"🔢 Number of molecules to depict (default: {suggested}): ").strip()
                    
                    # Handle empty input (use default)
                    if not user_input:
                        print(f"✅ Using default: {suggested} molecules")
                        return suggested
                    
                    # Handle cancel
                    if user_input.lower() in ['cancel', 'quit', 'exit', 'abort']:
                        return None
                    
                    # Handle numeric input
                    try:
                        count = int(user_input)
                        
                        # Handle special cases
                        if count == 0:
                            print("❌ Depiction cancelled (count = 0)")
                            return None
                        elif count == -1:
                            print(f"✅ Will depict ALL {total_available:,} molecules")
                            return -1
                        elif count < 0:
                            print("❌ Invalid input. Use -1 for all molecules or positive numbers")
                            continue
                        elif count > total_available:
                            print(f"⚠️  Requested {count:,} molecules, but only {total_available:,} available")
                            use_all = input(f"   Use all {total_available:,} molecules instead? (y/n): ").strip().lower()
                            if use_all in ['y', 'yes']:
                                print(f"✅ Will depict all {total_available:,} molecules")
                                return total_available
                            else:
                                continue
                        else:
                            print(f"✅ Will depict {count:,} molecules")
                            
                            # Ask about random sampling for large requests
                            if count > 25 and count < total_available:
                                random_choice = input(f"   Select molecules randomly? (y/n, default: n): ").strip().lower()
                                if random_choice in ['y', 'yes']:
                                    print("   🎲 Will use random sampling")
                                else:
                                    print("   📋 Will use first molecules in table")
                            
                            return count
                            
                    except ValueError:
                        print("❌ Invalid input. Please enter a number, -1, or 'cancel'")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ Molecule count selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error prompting for molecule count: {e}")
            return None

    def _depict_individual_molecules(self, compounds_df: pd.DataFrame, output_dir: str, 
                                   image_size: Tuple[int, int], highlight_mol, 
                                   highlight_substructures: bool) -> Tuple[int, int, int]:
        """
        Generate individual PNG images for each molecule.
        
        Args:
            compounds_df (pd.DataFrame): DataFrame containing compounds
            output_dir (str): Directory to save images
            image_size (Tuple[int, int]): Size of each image
            highlight_mol: RDKit molecule object for highlighting
            highlight_substructures (bool): Whether to highlight substructures
            
        Returns:
            Tuple[int, int, int]: (successful_depictions, failed_depictions, invalid_smiles)
        """
        from rdkit import Chem
        from rdkit.Chem import Draw
        
        successful_depictions = 0
        failed_depictions = 0
        invalid_smiles = 0
        
        # Initialize progress tracking
        if TQDM_AVAILABLE:
            progress_bar = tqdm(
                compounds_df.iterrows(),
                total=len(compounds_df),
                desc="Generating images",
                unit="molecules"
            )
        else:
            progress_bar = compounds_df.iterrows()
            processed = 0
        
        for idx, compound in progress_bar:
            try:
                smiles = compound['smiles']
                name = compound.get('name', f'compound_{idx}')
                compound_id = compound.get('id', idx)
                
                # Parse SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    invalid_smiles += 1
                    continue
                
                # Prepare highlighting if requested
                highlight_atoms = []
                if highlight_substructures and highlight_mol:
                    try:
                        matches = mol.GetSubstructMatches(highlight_mol)
                        if matches:
                            # Flatten all match tuples to get highlighted atoms
                            highlight_atoms = [atom for match in matches for atom in match]
                    except:
                        pass
                
                # Generate image
                if highlight_atoms:
                    img = Draw.MolToImage(
                        mol, 
                        size=image_size, 
                        highlightAtoms=highlight_atoms,
                        highlightColor=(1, 0.8, 0.8)  # Light red highlighting
                    )
                else:
                    img = Draw.MolToImage(mol, size=image_size)
                
                # Create filename
                safe_name = self._sanitize_filename(name)
                filename = f"{compound_id:06d}_{safe_name}.png"
                filepath = os.path.join(output_dir, filename)
                
                # Save image
                img.save(filepath)
                successful_depictions += 1
                
            except Exception as e:
                failed_depictions += 1
                if failed_depictions <= 5:  # Show first few errors
                    print(f"⚠️  Error depicting molecule '{name}': {e}")
            
            # Update progress for non-tqdm case
            if not TQDM_AVAILABLE:
                processed += 1
                if processed % max(1, len(compounds_df) // 10) == 0:
                    print(f"   📊 Processed {processed}/{len(compounds_df)} molecules...")
        
        return successful_depictions, failed_depictions, invalid_smiles
    
    def _depict_molecule_grids(self, compounds_df: pd.DataFrame, output_dir: str, 
                             image_size: Tuple[int, int], molecules_per_image: int,
                             highlight_mol, highlight_substructures: bool) -> Tuple[int, int, int]:
        """
        Generate grid images with multiple molecules per image.
        
        Args:
            compounds_df (pd.DataFrame): DataFrame containing compounds
            output_dir (str): Directory to save images
            image_size (Tuple[int, int]): Size of each molecule in the grid
            molecules_per_image (int): Number of molecules per grid image
            highlight_mol: RDKit molecule object for highlighting
            highlight_substructures (bool): Whether to highlight substructures
            
        Returns:
            Tuple[int, int, int]: (successful_depictions, failed_depictions, invalid_smiles)
        """
        from rdkit import Chem
        from rdkit.Chem import Draw
        
        successful_depictions = 0
        failed_depictions = 0
        invalid_smiles = 0
        
        # Calculate grid dimensions
        grid_cols = int(molecules_per_image ** 0.5)
        grid_rows = (molecules_per_image + grid_cols - 1) // grid_cols
        
        print(f"   📐 Grid layout: {grid_rows}x{grid_cols} ({molecules_per_image} molecules per image)")
        
        # Process molecules in batches
        total_molecules = len(compounds_df)
        num_grids = (total_molecules + molecules_per_image - 1) // molecules_per_image
        
        if TQDM_AVAILABLE:
            progress_bar = tqdm(range(num_grids), desc="Generating grid images", unit="grids")
        else:
            progress_bar = range(num_grids)
        
        for grid_idx in progress_bar:
            try:
                start_idx = grid_idx * molecules_per_image
                end_idx = min(start_idx + molecules_per_image, total_molecules)
                batch_df = compounds_df.iloc[start_idx:end_idx]
                
                mols = []
                legends = []
                highlight_atoms_list = []
                
                # Process each molecule in the batch
                for _, compound in batch_df.iterrows():
                    try:
                        smiles = compound['smiles']
                        name = compound.get('name', f'compound_{compound.get("id", "")}')
                        id = compound.get('id', f'compound_{compound.get("id", "")}')
                        
                        mol = Chem.MolFromSmiles(smiles)
                        if mol is None:
                            invalid_smiles += 1
                            continue
                        
                        mols.append(mol)
                        legends.append(f"{name[:20]} \n id: {str(id)}")  # Truncate long names
                        
                        # Prepare highlighting
                        highlight_atoms = []
                        if highlight_substructures and highlight_mol:
                            try:
                                matches = mol.GetSubstructMatches(highlight_mol)
                                if matches:
                                    highlight_atoms = [atom for match in matches for atom in match]
                            except:
                                pass
                        
                        highlight_atoms_list.append(highlight_atoms)
                        
                    except Exception:
                        failed_depictions += 1
                        continue
                
                if mols:
                    # Generate grid image
                    if highlight_substructures and any(highlight_atoms_list):
                        img = Draw.MolsToGridImage(
                            mols,
                            molsPerRow=grid_cols,
                            subImgSize=image_size,
                            legends=legends,
                            highlightAtomLists=highlight_atoms_list
                        )
                    else:
                        img = Draw.MolsToGridImage(
                            mols,
                            molsPerRow=grid_cols,
                            subImgSize=image_size,
                            legends=legends
                        )
                    
                    # Save grid image
                    filename = f"grid_{grid_idx + 1:04d}_{len(mols)}molecules.png"
                    filepath = os.path.join(output_dir, filename)
                    img.save(filepath)
                    
                    successful_depictions += len(mols)
                
            except Exception as e:
                failed_depictions += molecules_per_image
                print(f"⚠️  Error generating grid {grid_idx + 1}: {e}")
            
            # Progress update for non-tqdm case
            if not TQDM_AVAILABLE and (grid_idx + 1) % max(1, num_grids // 10) == 0:
                print(f"   📊 Generated {grid_idx + 1}/{num_grids} grid images...")
        
        return successful_depictions, failed_depictions, invalid_smiles
    
    def _sanitize_filename(self, filename: str) -> str:
        """
        Sanitize a string to be safe for use as a filename.
        
        Args:
            filename (str): Original filename
            
        Returns:
            str: Sanitized filename
        """
        import re
        
        # Replace problematic characters with underscores
        sanitized = re.sub(r'[<>:"/\\|?*\n\r\t]', '_', filename)
        
        # Remove multiple consecutive underscores
        sanitized = re.sub(r'_+', '_', sanitized)
        
        # Remove leading/trailing underscores and spaces
        sanitized = sanitized.strip('_ ')
        
        # Limit length
        if len(sanitized) > 50:
            sanitized = sanitized[:50]
        
        # Ensure it's not empty
        if not sanitized:
            sanitized = 'compound'
        
        return sanitized
    
    def _generate_depiction_report(self, table_name: str, output_dir: str, 
                                 total_molecules: int, successful_depictions: int,
                                 failed_depictions: int, invalid_smiles: int,
                                 image_size: Tuple[int, int], molecules_per_image: int) -> None:
        """
        Generate a summary report for the depiction process.
        
        Args:
            table_name (str): Name of the source table
            output_dir (str): Directory where images were saved
            total_molecules (int): Total number of molecules processed
            successful_depictions (int): Number of successful depictions
            failed_depictions (int): Number of failed depictions
            invalid_smiles (int): Number of invalid SMILES
            image_size (Tuple[int, int]): Size of each image
            molecules_per_image (int): Molecules per image
        """
        try:
            report_path = os.path.join(output_dir, 'depiction_report.txt')
            
            with open(report_path, 'w') as f:
                f.write("CHEMICAL STRUCTURE DEPICTION REPORT\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Project: {self.name}\n")
                f.write(f"Source Table: {table_name}\n")
                f.write(f"Generation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Output Directory: {output_dir}\n\n")
                
                f.write("PROCESSING SUMMARY\n")
                f.write("-" * 20 + "\n")
                f.write(f"Total molecules processed: {total_molecules}\n")
                f.write(f"Successful depictions: {successful_depictions}\n")
                f.write(f"Failed depictions: {failed_depictions}\n")
                f.write(f"Invalid SMILES: {invalid_smiles}\n")
                f.write(f"Success rate: {(successful_depictions/total_molecules)*100:.1f}%\n\n")
                
                f.write("IMAGE SETTINGS\n")
                f.write("-" * 15 + "\n")
                f.write(f"Image size: {image_size[0]}x{image_size[1]} pixels\n")
                f.write(f"Molecules per image: {molecules_per_image}\n")
                
                if molecules_per_image == 1:
                    f.write("Format: Individual molecule images\n")
                else:
                    f.write("Format: Grid images\n")
                
                f.write(f"\nImage format: PNG\n")
                
                # List generated files
                f.write(f"\nGENERATED FILES\n")
                f.write("-" * 15 + "\n")
                
                try:
                    files = sorted([f for f in os.listdir(output_dir) if f.endswith('.png')])
                    f.write(f"Total PNG files: {len(files)}\n\n")
                    
                    if len(files) <= 20:  # List all files if not too many
                        for file in files:
                            f.write(f"  {file}\n")
                    else:  # List first 10 and last 10
                        f.write("First 10 files:\n")
                        for file in files[:10]:
                            f.write(f"  {file}\n")
                        f.write(f"  ... ({len(files)-20} files omitted) ...\n")
                        f.write("Last 10 files:\n")
                        for file in files[-10:]:
                            f.write(f"  {file}\n")
                except:
                    f.write("Could not list generated files\n")
            
            print(f"📋 Depiction report saved: {report_path}")
            
        except Exception as e:
            print(f"⚠️  Warning: Could not generate depiction report: {e}")
    
    def depict_filtered_compounds(self, table_name: str, workflow_name: str,
                                save_images: bool = True, **kwargs) -> pd.DataFrame:
        """
        Apply a filtering workflow and then generate depictions of the filtered compounds.
        
        Args:
            table_name (str): Name of the table to filter and depict
            workflow_name (str): Name of the filtering workflow to apply
            save_images (bool): Whether to save structure images
            **kwargs: Additional arguments passed to depict_table()
            
        Returns:
            pd.DataFrame: Filtered compounds DataFrame
        """
        try:
            print(f"🔬 Filtering and depicting compounds from table '{table_name}'")
            
            # Apply filtering workflow
            filtered_df = self.filter_using_workflow(
                table_name, workflow_name, 
                save_results=False,  # Don't save to database, just get DataFrame
                **{k: v for k, v in kwargs.items() if k in ['parallel_threshold', 'max_workers', 'chunk_size']}
            )
            
            if filtered_df.empty:
                print("❌ No compounds passed the filtering workflow")
                return pd.DataFrame()
            
            if save_images:
                # Create temporary table for depiction
                temp_table_name = f"temp_filtered_{workflow_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                
                # Save filtered compounds to temporary table
                filter_results = {'workflow_name': workflow_name}
                success = self._save_workflow_filtered_compounds(
                    filtered_df, temp_table_name, workflow_name, filter_results
                )
                
                if success:
                    # Generate depictions
                    depiction_kwargs = {k: v for k, v in kwargs.items() 
                                      if k in ['output_dir', 'image_size', 'molecules_per_image', 
                                             'max_molecules', 'highlight_substructures', 'smarts_pattern']}
                    
                    # Set default output directory if not specified
                    if 'output_dir' not in depiction_kwargs:
                        depiction_kwargs['output_dir'] = os.path.join(
                            self.path, 'chemspace', 'depictions', f"{table_name}_filtered_{workflow_name}"
                        )
                    
                    depiction_success = self.depict_table(temp_table_name, **depiction_kwargs)
                    
                    # Clean up temporary table
                    self.drop_table(temp_table_name, confirm=False)
                    
                    if depiction_success:
                        print(f"✅ Filtered compounds depicted successfully")
                    else:
                        print(f"❌ Failed to generate depictions for filtered compounds")
                else:
                    print(f"❌ Failed to save filtered compounds for depiction")
            
            return filtered_df
            
        except Exception as e:
            print(f"❌ Error in depict_filtered_compounds: {e}")
            return pd.DataFrame()
        
    def create_reaction_workflow(self):
        """
        Create a reaction workflow dictionary by selecting reactions stored in the chem_reactions table in the projects database. A series of id are requested to the user to build the reaction workflow, until a -1 is provided. id are used to as key to retrieve the smarts value from chem_reactions table in the main projects database.
        
        Returns:
            dictionary of reactions to be applied in the reaction workflow
        """
        try:
            # Get the projects database path
            import site
            projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"
            
            if not os.path.exists(projects_db):
                print(f"❌ Projects database not found: {projects_db}")
                return {}
            
            print("🧪 Creating Reaction Workflow")
            print("=" * 80)
            print("Enter reaction IDs to build your reaction workflow.")
            print("Enter -1 when you're done adding reactions.")
            print("=" * 80)
            
            workflow_reactions = {}
            reaction_counter = 1
            
            while True:
                try:
                    # Request reaction ID from user
                    user_input = input(f"\n🔬 Enter reaction ID #{reaction_counter} (or -1 to finish): ").strip()
                    
                    # Check if user wants to finish
                    if user_input == "-1":
                        break
                    
                    # Validate input is a number
                    try:
                        reaction_id = int(user_input)
                    except ValueError:
                        print("❌ Please enter a valid number or -1 to finish")
                        continue
                    
                    # Skip if already added
                    if reaction_id in workflow_reactions:
                        print(f"⚠️  Reaction ID {reaction_id} is already in the workflow")
                        continue
                    
                    # Retrieve reaction from projects database
                    reaction_info = self._get_reaction_by_id(projects_db, reaction_id)
                    
                    if reaction_info:
                        reaction_name, reaction_smarts = reaction_info
                        workflow_reactions[reaction_id] = {
                            'name': reaction_name,
                            'smarts': reaction_smarts,
                            'order': reaction_counter
                        }
                        
                        print(f"✅ Added reaction: '{reaction_name}' (ID: {reaction_id})")
                        print(f"   🧪 Reaction SMARTS: {reaction_smarts}")
                        reaction_counter += 1
                    else:
                        print(f"❌ No reaction found with ID {reaction_id}")
                        
                except KeyboardInterrupt:
                    print("\n\n⏹️  Reaction workflow creation cancelled by user")
                    return {}
                except Exception as e:
                    print(f"❌ Error processing input: {e}")
                    continue
            
            # Display final workflow summary
            if workflow_reactions:
                print(f"\n✅ Reaction Workflow Created!")
                print("=" * 80)
                print(f"🧪 Total reactions: {len(workflow_reactions)}")
                print("\n📋 Workflow Summary:")
                
                for reaction_id, info in workflow_reactions.items():
                    print(f"   {info['order']}. '{info['name']}' (ID: {reaction_id})")
                    print(f"      🧪 SMARTS: {info['smarts']}")
                
                print("=" * 80)
                
                # Ask if user wants to save the workflow
                save_workflow = input("\n💾 Do you want to save this reaction workflow? (y/n): ").strip().lower()
                if save_workflow in ['y', 'yes']:
                    self._save_reaction_workflow(workflow_reactions)
                    
            else:
                print("\n⚠️  No reactions were added to the workflow")
            
            return workflow_reactions
            
        except Exception as e:
            print(f"❌ Error creating reaction workflow: {e}")
            return {}
    
    def _get_reaction_by_id(self, projects_db_path: str, reaction_id: int) -> Optional[Tuple[str, str]]:
        """
        Retrieve a specific reaction by ID from the projects database.
        
        Args:
            projects_db_path (str): Path to the projects database
            reaction_id (int): ID of the reaction to retrieve
            
        Returns:
            Optional[Tuple[str, str]]: (reaction_name, reaction_smarts) or None if not found
        """
        try:
            conn = sqlite3.connect(projects_db_path)
            cursor = conn.cursor()
            
            cursor.execute("SELECT reaction_name, smarts FROM chem_reactions WHERE id = ?", (reaction_id,))
            result = cursor.fetchone()
            conn.close()
            
            return result if result else None
            
        except Exception as e:
            print(f"❌ Error retrieving reaction ID {reaction_id}: {e}")
            return None
    
    def _save_reaction_workflow(self, workflow_reactions: Dict[int, Dict[str, Any]]) -> None:
        """
        Save the reaction workflow to a table named 'reaction_workflows' within chemspace.db.
        
        Args:
            workflow_reactions (Dict[int, Dict]): Dictionary of reaction workflow data
        """
        try:
            if not workflow_reactions:
                print("⚠️  No workflow reactions to save")
                return
            
            # Get workflow name from user
            workflow_name = input("📝 Enter a name for this reaction workflow: ").strip()
            if not workflow_name:
                workflow_name = f"reaction_workflow_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                print(f"📋 Using default name: {workflow_name}")
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create reaction_workflows table if it doesn't exist
            create_table_query = """
            CREATE TABLE IF NOT EXISTS reaction_workflows (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                workflow_name TEXT NOT NULL UNIQUE,
                reactions_dict TEXT NOT NULL,
                creation_date TEXT NOT NULL,
                description TEXT
            )
            """
            
            cursor.execute(create_table_query)
            
            # Create index for faster searches
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_reaction_workflows_name ON reaction_workflows(workflow_name)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_reaction_workflows_date ON reaction_workflows(creation_date)")
            
            # Get current timestamp
            creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # Get optional description from user
            description = input("📄 Enter a description for this reaction workflow (optional): ").strip()
            if not description:
                description = f"Reaction workflow with {len(workflow_reactions)} reactions"
            
            # Check if workflow name already exists
            cursor.execute("SELECT COUNT(*) FROM reaction_workflows WHERE workflow_name = ?", (workflow_name,))
            if cursor.fetchone()[0] > 0:
                overwrite = input(f"⚠️  Reaction workflow '{workflow_name}' already exists. Overwrite? (y/n): ").strip().lower()
                if overwrite in ['y', 'yes']:
                    # Delete existing workflow
                    cursor.execute("DELETE FROM reaction_workflows WHERE workflow_name = ?", (workflow_name,))
                    print(f"🔄 Overwriting existing reaction workflow '{workflow_name}'")
                else:
                    # Generate unique name
                    counter = 1
                    original_name = workflow_name
                    while True:
                        new_name = f"{original_name}_{counter}"
                        cursor.execute("SELECT COUNT(*) FROM reaction_workflows WHERE workflow_name = ?", (new_name,))
                        if cursor.fetchone()[0] == 0:
                            workflow_name = new_name
                            break
                        counter += 1
                    print(f"📝 Using name: {workflow_name}")
            
            # Convert workflow to JSON string for storage
            reactions_json = json.dumps(workflow_reactions, sort_keys=True)
            
            # Insert workflow
            insert_query = """
            INSERT INTO reaction_workflows 
            (workflow_name, reactions_dict, creation_date, description)
            VALUES (?, ?, ?, ?)
            """
            
            cursor.execute(insert_query, (workflow_name, reactions_json, creation_date, description))
            
            conn.commit()
            conn.close()
            
            print(f"✅ Successfully saved reaction workflow!")
            print(f"   📋 Workflow name: '{workflow_name}'")
            print(f"   🧪 Reactions saved: {len(workflow_reactions)}")
            print(f"   📅 Created: {creation_date}")
            print(f"   📄 Description: {description}")
            
        except sqlite3.IntegrityError as e:
            print(f"❌ Database integrity error saving reaction workflow: {e}")
        except Exception as e:
            print(f"❌ Error saving reaction workflow: {e}")

    def list_reaction_workflows(self) -> None:
        """
        Display all saved reaction workflows in a formatted table with IDs.
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Check if reaction_workflows table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='reaction_workflows'")
            if not cursor.fetchone():
                print("📝 No reaction workflows table found")
                conn.close()
                return
            
            # Get all workflows with summary information including ID
            cursor.execute("""
            SELECT id, workflow_name, creation_date, description, reactions_dict
            FROM reaction_workflows 
            ORDER BY creation_date DESC
            """)
            
            workflows = cursor.fetchall()
            conn.close()
            
            if not workflows:
                print("📝 No saved reaction workflows found")
                return
            
            print("\n" + "="*100)
            print(f"SAVED REACTION WORKFLOWS - Project: {self.name}")
            print("="*100)
            
            for i, (workflow_id, name, date, desc, reactions_dict_str) in enumerate(workflows, 1):
                try:
                    reactions_dict = json.loads(reactions_dict_str)
                    reaction_count = len(reactions_dict)
                    
                    # Get reaction names for summary
                    reaction_names = [info['name'] for info in reactions_dict.values()]
                    reactions_summary = ', '.join(reaction_names[:3])
                    if len(reaction_names) > 3:
                        reactions_summary += f" and {len(reaction_names) - 3} more..."
                    
                except json.JSONDecodeError:
                    reaction_count = 0
                    reactions_summary = "Error parsing reactions"
                
                print(f"\n🧪 Workflow name: '{name}' (ID: {workflow_id})")
                print(f"   📅 Created: {date}")
                print(f"   🔬 Reactions: {reaction_count}")
                print(f"   📄 Description: {desc}")
                print(f"   🧪 Reactions: {reactions_summary}")
            
            print("="*100)
            
        except Exception as e:
            print(f"❌ Error listing reaction workflows: {e}")

    def _load_reaction_workflow_by_id(self, workflow_id: int) -> Optional[Dict[str, Any]]:
        """
        Load a reaction workflow from the database by its ID.
        
        Args:
            workflow_id (int): ID of the reaction workflow to load
            
        Returns:
            Optional[Dict]: Workflow data dictionary or None if not found
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Check if reaction_workflows table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='reaction_workflows'")
            if not cursor.fetchone():
                print("❌ No reaction workflows table found")
                conn.close()
                return None
            
            # Load workflow by ID
            cursor.execute("""
            SELECT workflow_name, reactions_dict, creation_date, description
            FROM reaction_workflows 
            WHERE id = ?
            """, (workflow_id,))
            
            result = cursor.fetchone()
            conn.close()
            
            if not result:
                # Show available workflows if ID not found
                print(f"❌ No workflow found with ID {workflow_id}")
                self._show_available_reaction_workflows()
                return None
            
            workflow_name, reactions_dict_str, creation_date, description = result
            
            # Parse reactions dictionary
            try:
                reactions_dict = json.loads(reactions_dict_str)
            except json.JSONDecodeError as e:
                print(f"❌ Error parsing workflow data: {e}")
                return None
            
            return {
                'workflow_name': workflow_name,
                'reactions_dict': reactions_dict,
                'creation_date': creation_date,
                'description': description
            }
            
        except Exception as e:
            print(f"❌ Error loading reaction workflow {workflow_id}: {e}")
            return None
    
    def _show_available_reaction_workflows(self) -> None:
        """
        Display available reaction workflows with their IDs.
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("""
            SELECT id, workflow_name, creation_date, description
            FROM reaction_workflows 
            ORDER BY creation_date DESC
            """)
            
            workflows = cursor.fetchall()
            conn.close()
            
            if workflows:
                print("\n📋 Available Reaction Workflows:")
                print("-" * 70)
                print(f"{'ID':<4} {'Name':<25} {'Created':<12} {'Description':<25}")
                print("-" * 70)
                
                for wf_id, name, date, desc in workflows:
                    date_short = date[:10] if date else "Unknown"
                    desc_short = desc[:22] + "..." if len(desc) > 25 else desc
                    print(f"{wf_id:<4} {name[:24]:<25} {date_short:<12} {desc_short:<25}")
                
                print("-" * 70)
            else:
                print("📝 No reaction workflows found.")
                
        except Exception as e:
            print(f"❌ Error showing available workflows: {e}")
    
    def _select_tables_for_reaction(self, available_tables: List[str]) -> List[str]:
        """
        Let user select tables for reaction processing.
        
        Args:
            available_tables (List[str]): List of available table names
            
        Returns:
            List[str]: List of selected table names
        """
        try:
            print(f"\n📋 Table Selection for Reaction Workflow")
            print("-" * 50)
            print("Enter table numbers or names (comma-separated).")
            
            user_input = input("Selection: ").strip()
            
            if not user_input:
                print("No tables selected. Stopping.")
                sys.exit(0)
            
            selected_tables = []
            selections = [s.strip() for s in user_input.split(',')]
            
            for selection in selections:
                try:
                    # Try as table number first
                    table_idx = int(selection) - 1
                    if 0 <= table_idx < len(available_tables):
                        table_name = available_tables[table_idx]
                        if table_name not in selected_tables:
                            selected_tables.append(table_name)
                    else:
                        print(f"⚠️  Invalid table number: {selection}")
                except ValueError:
                    # Treat as table name
                    if selection in available_tables:
                        if selection not in selected_tables:
                            selected_tables.append(selection)
                    else:
                        print(f"⚠️  Table '{selection}' not found")
            
            return selected_tables
            
        except Exception as e:
            print(f"❌ Error selecting tables: {e}")
            return []
    
    def _apply_reactions_to_table(self, compounds_df: pd.DataFrame, 
                                reactions_dict: Dict[int, Dict[str, Any]], 
                                workflow_name: str, 
                                source_table_name: str) -> Dict[str, int]:
        """
        Apply reaction workflow to compounds in a table.
        
        Args:
            compounds_df (pd.DataFrame): DataFrame containing compounds
            reactions_dict (Dict): Dictionary of reactions to apply
            workflow_name (str): Name of the workflow
            source_table_name (str): Name of the source table
            
        Returns:
            Dict[str, int]: Results dictionary with counts
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Sort reactions by order
            sorted_reactions = sorted(reactions_dict.items(), 
                                    key=lambda x: x[1]['order'])
            
            # Initialize tracking
            all_products = []
            compounds_processed = 0
            products_generated = 0
            failed_reactions = 0
            

            # Process each compound through the reaction sequence
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    compounds_df.iterrows(),
                    total=len(compounds_df),
                    desc=f"Processing {source_table_name}",
                    unit="compounds"
                )
            else:
                progress_bar = compounds_df.iterrows()
                processed_count = 0
            
            for _, compound in progress_bar:
                compounds_processed += 1
                try:
                    # Start with original compound
                    reactant_smiles = compound['smiles']
                    compound_name = compound.get('name', f'compound_{compound.get("id", compounds_processed)}')
                    
                    # Parse reactant molecule
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol is None:
                        failed_reactions += 1
                        continue
                    
                    # Apply each reaction in sequence
                    current_products = [(reactant_mol, compound_name, 0)]  # (mol, name, reaction_step)
                    
                    for reaction_id, reaction_info in sorted_reactions:
                        new_products = []
                        reaction_smarts = reaction_info['smarts']
                        reaction_name = reaction_info['name']
                        step = reaction_info['order']
                        
                        try:
                            # Parse reaction SMARTS
                            rxn = AllChem.ReactionFromSmarts(reaction_smarts)
                            if rxn is None:
                                print(f"   ⚠️  Invalid reaction SMARTS: {reaction_smarts}")
                                continue
                            
                            # Apply reaction to all current products
                            for product_mol, product_name, prev_step in current_products:
                                try:
                                    # Run reaction
                                    reaction_products = rxn.RunReactants((product_mol,))
                                    
                                    

                                    for product_set_idx, product_set in enumerate(reaction_products):
                                        for product_idx, product_mol_new in enumerate(product_set):
                                            # Sanitize product
                                            try:
                                                Chem.SanitizeMol(product_mol_new)
                                                
                                                # Generate new product name
                                                new_product_name = f"{product_name}_{reaction_name}_{step}_{product_set_idx}_{product_idx}"
                                                new_products.append((product_mol_new, new_product_name, step))
                                                
                                            except:
                                                # Skip invalid products
                                                continue
                                                
                                except Exception:
                                    # Skip failed reactions
                                    continue
                        
                        except Exception:
                            print(f"   ⚠️  Error with reaction '{reaction_name}': skipping")
                            continue
                        
                        # Update current products for next reaction
                        if new_products:
                            current_products = new_products
                        # If no new products, keep current products for next reaction
                    
                    # Add final products to results
                    for product_mol, product_name, step in current_products:
                        if step > 0:  # Only count products from actual
                            try:
                                product_smiles = Chem.MolToSmiles(product_mol)
                                all_products.append({
                                    'smiles': product_smiles,
                                    'name': product_name,
                                    'flag': 'reaction_product',
                                    'source_compound': compound_name,
                                    'source_table': source_table_name,
                                    'workflow': workflow_name,
                                    'final_step': step
                                })
                                products_generated += 1
                            except:
                                failed_reactions += 1
                    
                except Exception:
                    failed_reactions += 1
                
                # Progress update for non-tqdm case
                if not TQDM_AVAILABLE:
                    processed_count += 1
                    if processed_count % max(1, len(compounds_df) // 10) == 0:
                        print(f"      📊 Processed {processed_count}/{len(compounds_df)} compounds...")
            
            # Save products to new table
            output_table_name = None
            if all_products:

                print(f"\n💾 Found {products_generated} reaction products from table '{source_table_name}'")

                save_choice = input("Do you want to save these products to a new table? (y/n): ").strip().lower()

                if save_choice in ['y', 'yes']:

                    output_table_name = self._save_reaction_products(
                        all_products, workflow_name, source_table_name
                    )
            
            return {
                'compounds_processed': compounds_processed,
                'products_generated': products_generated,
                'failed_reactions': failed_reactions,
                'output_table': output_table_name
            }
            
        except Exception as e:
            print(f"   ❌ Error processing table reactions: {e}")
            return {
                'compounds_processed': 0,
                'products_generated': 0,
                'failed_reactions': len(compounds_df),
                'output_table': None
            }
    
    def _save_reaction_products(self, products: List[Dict[str, Any]], 
                              workflow_name: str, 
                              source_table_name: str) -> Optional[str]:
        """
        Save reaction products to a new table.
        
        Args:
            products (List[Dict]): List of product dictionaries
            workflow_name (str): Name of the workflow
            source_table_name (str): Name of the source table
            
        Returns:
            Optional[str]: Name of the created table or None if failed
        """
        try:
            if not products:
                return None
            
            # Generate output table name
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_table_name = f"{source_table_name}_products_{workflow_name}_{timestamp}"
            output_table_name = self._sanitize_table_name(output_table_name)
            
            # Create table
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create products table with extended schema
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {output_table_name} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT NOT NULL,
                    name TEXT,
                    flag TEXT,
                    source_compound TEXT,
                    source_table TEXT,
                    workflow TEXT,
                    final_step INTEGER,
                    creation_date TEXT,
                    UNIQUE(smiles, name)
                )
            ''')
            
            # Create indexes
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_smiles ON {output_table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_source ON {output_table_name}(source_compound)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_workflow ON {output_table_name}(workflow)")
            
            # Insert products
            insert_query = f'''
                INSERT OR IGNORE INTO {output_table_name} 
                (smiles, name, flag, source_compound, source_table, workflow, final_step, creation_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            '''
            
            creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            inserted_count = 0
            
            for product in products:
                try:
                    cursor.execute(insert_query, (
                        product['smiles'],
                        product['name'],
                        product.get('flag', 'reaction_product'),
                        product.get('source_compound', 'unknown'),
                        product.get('source_table', source_table_name),
                        product.get('workflow', workflow_name),
                        product.get('final_step', 1),
                        creation_date
                    ))
                    inserted_count += 1
                except sqlite3.IntegrityError:
                    # Skip duplicates
                    continue
            
            conn.commit()
            conn.close()
            
            print(f"      💾 Created table '{output_table_name}' with {inserted_count} products")
            return output_table_name
            
        except Exception as e:
            print(f"      ❌ Error saving reaction products: {e}")
            return None
    
    def list_reaction_products_tables(self) -> None:
        """
        List all tables containing reaction products.
        """
        try:
            all_tables = self.get_all_tables()
            product_tables = [table for table in all_tables if '_products_' in table]
            
            if not product_tables:
                print("📝 No reaction product tables found")
                return
            
            print(f"\n📋 Reaction Product Tables ({len(product_tables)} found):")
            print("=" * 80)
            
            for i, table_name in enumerate(product_tables, 1):
                compound_count = self.get_compound_count(table_name=table_name)
                
                # Try to extract workflow info from table name
                parts = table_name.split('_products_')
                source_table = parts[0] if len(parts) > 1 else "unknown"
                workflow_info = parts[1] if len(parts) > 1 else "unknown"
                
                print(f"{i:2d}. {table_name}")
                print(f"     📊 Products: {compound_count}")
                print(f"     📋 Source: {source_table}")
                print(f"     🧪 Workflow: {workflow_info}")
                print()
            
        except Exception as e:
            print(f"❌ Error listing product tables: {e}")    

    def create_filtering_workflow(self):
        """
        Create a filtering workflow dictionary by selecting filters stored in the chem_filters table in the projects database. A series of id are requested to the user to build the filtering workflow, until a -1 is provided. id are used to as key to retrieve the smarts value from chem_filters table in the main projects database.
        
        Returns:
            dict: Dictionary of filters to be applied in the filtering workflow
        """
        try:
            # Get the projects database path
            import site
            projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"
            
            if not os.path.exists(projects_db):
                print(f"❌ Projects database not found: {projects_db}")
                return {}
            
            # Show available filters with enhanced display
            #self._display_available_filters_enhanced(projects_db)
            tidyscreen.list_chemical_filters()
            
            workflow_filters = {}
            filter_counter = 1
            
            print(f"\n📋 Building Filtering Workflow")
            print("=" * 80)
            print("Commands:")
            print("   • Enter filter ID to add a filter")
            print("   • Type 'list' to show available filters again")
            print("   • Type 'search <term>' to search filters by name")
            print("   • Type 'preview <id>' to preview a filter")
            print("   • Type 'batch <id1,id2,id3>' to add multiple filters")
            print("   • Type 'remove <filter_name>' to remove a filter")
            print("   • Type 'show' to display current workflow")
            print("   • Type '-1' or 'done' to finish")
            print("=" * 80)
            
            while True:
                try:
                    # Enhanced command input
                    user_input = input(f"\n🔍 Command or Filter ID #{filter_counter}: ").strip().lower()
                    
                    # Handle different commands
                    if user_input in ['-1', 'done', 'finish', 'exit']:
                        break
                    elif user_input == 'list':
                        self._display_available_filters_enhanced(projects_db)
                        continue
                    elif user_input.startswith('search '):
                        search_term = user_input[7:].strip()
                        self._search_filters(projects_db, search_term)
                        continue
                    elif user_input.startswith('preview '):
                        try:
                            preview_id = int(user_input[8:].strip())
                            self._preview_filter(projects_db, preview_id)
                        except ValueError:
                            print("❌ Invalid filter ID for preview")
                        continue
                    elif user_input.startswith('batch '):
                        batch_ids = user_input[6:].strip()
                        self._add_batch_filters(projects_db, batch_ids, workflow_filters)
                        filter_counter += len([id for id in batch_ids.split(',') if id.strip().isdigit()])
                        continue
                    elif user_input.startswith('remove '):
                        filter_name = user_input[7:].strip()
                        if filter_name in workflow_filters:
                            del workflow_filters[filter_name]
                            print(f"✅ Removed filter: '{filter_name}'")
                        else:
                            print(f"❌ Filter '{filter_name}' not found in workflow")
                        continue
                    elif user_input == 'show':
                        self._show_current_workflow(workflow_filters)
                        continue
                    elif user_input == 'help':
                        self._show_workflow_help()
                        continue
                    
                    # Handle numeric filter ID input
                    try:
                        filter_id = int(user_input)
                    except ValueError:
                        print("❌ Invalid command. Type 'help' for available commands or enter a filter ID")
                        continue
                    
                    # Skip if already added
                    filter_exists = any(fid == filter_id for fid in [self._get_filter_id_by_name(name, projects_db) 
                                        for name in workflow_filters.keys()])
                    
                    if filter_exists:
                        print(f"⚠️  Filter ID {filter_id} is already in the workflow")
                        continue
                    
                    # Get required instances with validation
                    instances = self._get_required_instances()
                    if instances is None:
                        continue
                    
                    # Retrieve and validate filter from projects database
                    filter_info = self._get_filter_by_id_enhanced(projects_db, filter_id)
                    
                    if filter_info:
                        filter_name, smarts_pattern = filter_info
                        
                        # Validate SMARTS pattern
                        if self._validate_smarts_pattern(smarts_pattern):
                            workflow_filters[filter_name] = {
                                'instances': instances,
                                'smarts': smarts_pattern,
                                'filter_id': filter_id
                            }
                            
                            print(f"✅ Added filter #{filter_counter}: '{filter_name}' (ID: {filter_id})")
                            print(f"   🧪 SMARTS: {smarts_pattern}")
                            print(f"   🔢 Required instances: {instances}")
                            filter_counter += 1
                        else:
                            print(f"❌ Invalid SMARTS pattern for filter ID {filter_id}")
                    else:
                        print(f"❌ No filter found with ID {filter_id}")
                        
                except KeyboardInterrupt:
                    print("\n\n⏹️  Workflow creation cancelled by user")
                    return {}
                except Exception as e:
                    print(f"❌ Error processing input: {e}")
                    continue
            
            # Display final workflow summary
            if workflow_filters:
                print(f"\n✅ Enhanced Filtering Workflow Created!")
                print("=" * 70)
                print(f"📊 Total filters: {len(workflow_filters)}")
                print(f"🧪 Total required instances: {sum(f['instances'] for f in workflow_filters.values())}")
                
                self._display_filtering_workflow_summary(workflow_filters)
                
                # Enhanced workflow options
                print("\n📋 Workflow Options:")
                print("=" * 30)
                
                # Option to save workflow
                save_workflow = input("💾 Save this workflow? (y/n): ").strip().lower()
                if save_workflow in ['y', 'yes']:
                    # Convert to simple format for compatibility
                    simple_workflow = {name: info['instances'] for name, info in workflow_filters.items()}
                    self._save_filtering_workflow(simple_workflow)
                    
            else:
                print("\n⚠️  No filters were added to the workflow")
            
            print("\n🏁 Filtering workflow creation completed.")
            
            # Return simple format for compatibility with existing methods
            return {name: info['instances'] for name, info in workflow_filters.items()}

        except Exception as e:
            print(f"❌ Error creating filtering workflow: {e}")
            return {}

    def _get_filter_id_by_name(self, filter_name: str, projects_db_path: str = None) -> Optional[int]:
        """
        Get filter ID from the projects database by filter name.
        
        Args:
            filter_name (str): Name of the filter to search for
            projects_db_path (str, optional): Path to projects database. If None, uses default path
            
        Returns:
            Optional[int]: Filter ID if found, None otherwise
        """
        try:
            # Get projects database path if not provided
            if projects_db_path is None:
                import site
                projects_db_path = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"
            
            # Check if database exists
            if not os.path.exists(projects_db_path):
                print(f"❌ Projects database not found: {projects_db_path}")
                return None
            
            # Connect to database and search for filter
            conn = sqlite3.connect(projects_db_path)
            cursor = conn.cursor()
            
            # Query for exact filter name match
            cursor.execute(
                "SELECT id FROM chem_filters WHERE filter_name = ? LIMIT 1",
                (filter_name,)
            )
            
            result = cursor.fetchone()
            conn.close()
            
            if result:
                filter_id = result[0]
                return filter_id
            else:
                return None
                
        except Exception as e:
            print(f"❌ Error retrieving filter ID for '{filter_name}': {e}")
            return None

    def _display_available_filters_enhanced(self, projects_db_path: str) -> None:
        """
        Display available filters with enhanced formatting and pagination.
        
        Args:
            projects_db_path (str): Path to the projects database
        """
        try:
            conn = sqlite3.connect(projects_db_path)
            cursor = conn.cursor()
            
            # Get all available filters with additional info
            cursor.execute("""
                SELECT id, filter_name, smarts 
                FROM chem_filters 
                ORDER BY id
            """)
            filters = cursor.fetchall()
            conn.close()
            
            if not filters:
                print("⚠️  No chemical filters found in projects database")
                return
            
            print(f"\n📋 Available Chemical Filters ({len(filters)} total)")
            print("=" * 90)
            print(f"{'ID':<4} {'Filter Name':<25} {'SMARTS Preview':<30}")
            print("=" * 90)
            
            # Display filters with pagination
            page_size = 15
            total_pages = (len(filters) + page_size - 1) // page_size
            
            for page in range(total_pages):
                start_idx = page * page_size
                end_idx = min(start_idx + page_size, len(filters))
                
                for filter_id, filter_name, smarts in filters[start_idx:end_idx]:
                    # Truncate long strings for display
                    name_display = filter_name[:24] if len(filter_name) > 24 else filter_name
                    smarts_display = smarts[:29] + "..." if len(smarts) > 30 else smarts
                    
                    
                    print(f"{filter_id:<4} {name_display:<25} {smarts_display:<30}")
                
                # Show pagination info
                if total_pages > 1:
                    print(f"\n📄 Page {page + 1}/{total_pages} - Showing filters {start_idx + 1}-{end_idx}")
                    if page < total_pages - 1:
                        continue_display = input("Press Enter to see more filters, or 'q' to stop: ").strip().lower()
                        if continue_display == 'q':
                            break
            
            print("=" * 90)
            
        except Exception as e:
            print(f"❌ Error displaying available filters: {e}")
    
    def _search_filters(self, projects_db_path: str, search_term: str) -> None:
        """
        Search filters by name or description.
        
        Args:
            projects_db_path (str): Path to the projects database
            search_term (str): Term to search for
        """
        try:
            conn = sqlite3.connect(projects_db_path)
            cursor = conn.cursor()
            
            # Search in filter name and description
            cursor.execute("""
                SELECT id, filter_name, smarts, description 
                FROM chem_filters 
                WHERE filter_name LIKE ? OR description LIKE ?
                ORDER BY filter_name
            """, (f"%{search_term}%", f"%{search_term}%"))
            
            results = cursor.fetchall()
            conn.close()
            
            if results:
                print(f"\n🔍 Search Results for '{search_term}' ({len(results)} found):")
                print("-" * 80)
                print(f"{'ID':<4} {'Filter Name':<25} {'SMARTS Preview':<30} {'Description':<25}")
                print("-" * 80)
                
                for filter_id, filter_name, smarts, description in results:
                    name_display = filter_name[:24] if len(filter_name) > 24 else filter_name
                    smarts_display = smarts[:29] + "..." if len(smarts) > 30 else smarts
                    desc_display = (description or "No description")[:24]
                    if len(desc_display) > 24:
                        desc_display = desc_display[:21] + "..."
                    
                    print(f"{filter_id:<4} {name_display:<25} {smarts_display:<30} {desc_display:<25}")
                print("-" * 80)
            else:
                print(f"🔍 No filters found matching '{search_term}'")
                
        except Exception as e:
            print(f"❌ Error searching filters: {e}")
    
    def _preview_filter(self, projects_db_path: str, filter_id: int) -> None:
        """
        Display detailed preview of a specific filter.
        
        Args:
            projects_db_path (str): Path to the projects database
            filter_id (int): ID of the filter to preview
        """
        try:
            filter_info = self._get_filter_by_id_enhanced(projects_db_path, filter_id)
            
            if filter_info:
                filter_name, smarts_pattern, description = filter_info
                
                print(f"\n🔍 Filter Preview - ID: {filter_id}")
                print("=" * 50)
                print(f"📋 Name: {filter_name}")
                print(f"🧪 SMARTS: {smarts_pattern}")
                print(f"📝 Description: {description or 'No description available'}")
                
                # Validate SMARTS pattern
                is_valid = self._validate_smarts_pattern(smarts_pattern)
                validity_icon = "✅" if is_valid else "❌"
                print(f"{validity_icon} SMARTS Validity: {'Valid' if is_valid else 'Invalid'}")
                
                # Show pattern complexity
                complexity = len(smarts_pattern)
                if complexity < 20:
                    complexity_level = "Simple"
                elif complexity < 50:
                    complexity_level = "Moderate"
                else:
                    complexity_level = "Complex"
                
                print(f"📏 Pattern Complexity: {complexity_level} ({complexity} characters)")
                print("=" * 50)
            else:
                print(f"❌ No filter found with ID {filter_id}")
                
        except Exception as e:
            print(f"❌ Error previewing filter: {e}")
    
    def _add_batch_filters(self, projects_db_path: str, batch_ids: str, workflow_filters: Dict) -> None:
        """
        Add multiple filters at once using batch input.
        
        Args:
            projects_db_path (str): Path to the projects database
            batch_ids (str): Comma-separated filter IDs
            workflow_filters (Dict): Current workflow filters dictionary
        """
        try:
            ids = [id_str.strip() for id_str in batch_ids.split(',') if id_str.strip()]
            valid_ids = []
            
            # Validate all IDs first
            for id_str in ids:
                try:
                    filter_id = int(id_str)
                    if self._get_filter_by_id_enhanced(projects_db_path, filter_id):
                        valid_ids.append(filter_id)
                    else:
                        print(f"⚠️  Invalid filter ID: {filter_id}")
                except ValueError:
                    print(f"⚠️  Invalid ID format: {id_str}")
            
            if not valid_ids:
                print("❌ No valid filter IDs provided")
                return
            
            # Get default instances for batch
            print(f"\n📦 Adding {len(valid_ids)} filters in batch mode")
            default_instances = self._get_required_instances("Enter default instances for all filters")
            if default_instances is None:
                return
            
            added_count = 0
            for filter_id in valid_ids:
                filter_info = self._get_filter_by_id_enhanced(projects_db_path, filter_id)
                if filter_info:
                    filter_name, smarts_pattern, description = filter_info
                    
                    # Skip if already exists
                    if filter_name in workflow_filters:
                        print(f"⚠️  Filter '{filter_name}' already in workflow, skipping")
                        continue
                    
                    # Validate SMARTS
                    if self._validate_smarts_pattern(smarts_pattern):
                        workflow_filters[filter_name] = {
                            'instances': default_instances,
                            'smarts': smarts_pattern,
                            'description': description or 'No description available',
                            'filter_id': filter_id
                        }
                        added_count += 1
                        print(f"✅ Added: '{filter_name}' (ID: {filter_id})")
                    else:
                        print(f"❌ Invalid SMARTS for filter ID {filter_id}, skipped")
            
            print(f"📦 Batch operation completed: {added_count} filters added")
            
        except Exception as e:
            print(f"❌ Error in batch filter addition: {e}")
    
    def _get_required_instances(self, prompt: str = "🔢 Enter required instances") -> Optional[int]:
        """
        Get required instances with validation.
        
        Args:
            prompt (str): Custom prompt text
            
        Returns:
            Optional[int]: Number of required instances or None if invalid
        """
        try:
            while True:
                user_input = input(f"{prompt}: ").strip()
                
                if user_input.lower() in ['cancel', 'skip', 'back']:
                    return None
                
                try:
                    instances = int(user_input)
                    if instances < 0:
                        print("❌ Instances must be non-negative")
                        continue
                    return instances
                except ValueError:
                    print("❌ Please enter a valid number")
                    continue
        except KeyboardInterrupt:
            return None
    
    def _get_filter_by_id_enhanced(self, projects_db_path: str, filter_id: int) -> Optional[Tuple[str, str, str]]:
        """
        Enhanced version of filter retrieval with description.
        
        Args:
            projects_db_path (str): Path to the projects database
            filter_id (int): ID of the filter to retrieve
            
        Returns:
            Optional[Tuple[str, str, str]]: (filter_name, smarts_pattern, description) or None if not found
        """
        try:
            conn = sqlite3.connect(projects_db_path)
            cursor = conn.cursor()
            
            cursor.execute("""
                SELECT filter_name, smarts 
                FROM chem_filters 
                WHERE id = ?
            """, (filter_id,))
            result = cursor.fetchone()
            conn.close()
            
            return result if result else None
            
        except Exception as e:
            print(f"❌ Error retrieving filter ID {filter_id}: {e}")
            return None
    
    def _validate_smarts_pattern(self, smarts_pattern: str) -> bool:
        """
        Validate a SMARTS pattern using RDKit.
        
        Args:
            smarts_pattern (str): SMARTS pattern to validate
            
        Returns:
            bool: True if valid, False otherwise
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmarts(smarts_pattern)
            return mol is not None
        except:
            return False
    
    def _show_current_workflow(self, workflow_filters: Dict) -> None:
        """
        Display the current workflow being built.
        
        Args:
            workflow_filters (Dict): Current workflow filters
        """
        if not workflow_filters:
            print("📋 Current workflow is empty")
            return
        
        print(f"\n📋 Current Filtering Workflow ({len(workflow_filters)} filters)")
        print("=" * 80)
        print(f"{'#':<3} {'Filter Name':<25} {'Instances':<10} {'SMARTS Preview':<35}")
        print("=" * 80)
        
        for i, (name, info) in enumerate(workflow_filters.items(), 1):
            smarts_preview = info['smarts'][:34] + "..." if len(info['smarts']) > 35 else info['smarts']
            print(f"{i:<3} {name[:24]:<25} {info['instances']:<10} {smarts_preview:<35}")
        
        print("=" * 80)
        total_instances = sum(info['instances'] for info in workflow_filters.values())
        print(f"📊 Total filters: {len(workflow_filters)} | Total instances: {total_instances}")
    
    def _show_workflow_help(self) -> None:
        """
        Display help information for workflow creation commands.
        """
        print("\n📖 Enhanced Filtering Workflow Help")
        print("=" * 50)
        print("Available Commands:")
        print("   list                    - Show all available filters")
        print("   search <term>          - Search filters by name/description")
        print("   preview <id>           - Preview detailed filter information")
        print("   batch <id1,id2,id3>    - Add multiple filters at once")
        print("   remove <filter_name>   - Remove a filter from workflow")
        print("   show                   - Display current workflow")
        print("   help                   - Show this help message")
        print("   done / -1              - Finish workflow creation")
        print("\nFilter Management:")
        print("   • Enter a number to add a filter by ID")
        print("   • Each filter requires a number of instances")
        print("   • SMARTS patterns are validated automatically")
        print("   • Duplicate filters are automatically detected")
        print("=" * 50)
    
    def _display_filtering_workflow_summary(self, workflow_filters: Dict) -> None:
        """
        Display a comprehensive summary of the created workflow.
        
        Args:
            workflow_filters (Dict): Complete workflow filters with metadata
        """
        print("\n📋 Workflow Summary:")
        print("-" * 60)
        
        for i, (name, info) in enumerate(workflow_filters.items(), 1):
            print(f"\n{i}. {name} (ID: {info['filter_id']})")
            print(f"   🔢 Required instances: {info['instances']}")
            print(f"   🧪 SMARTS: {info['smarts']}")
        
        print("-" * 60)
        
        # Calculate workflow statistics
        total_instances = sum(info['instances'] for info in workflow_filters.values())
        avg_instances = total_instances / len(workflow_filters) if workflow_filters else 0
        
        print(f"\n📊 Workflow Statistics:")
        print(f"   • Total filters: {len(workflow_filters)}")
        print(f"   • Total instances: {total_instances}")
        print(f"   • Average instances per filter: {avg_instances:.1f}")
        
        # Pattern complexity analysis
        simple_patterns = sum(1 for info in workflow_filters.values() if len(info['smarts']) < 20)
        moderate_patterns = sum(1 for info in workflow_filters.values() if 20 <= len(info['smarts']) < 50)
        complex_patterns = sum(1 for info in workflow_filters.values() if len(info['smarts']) >= 50)
        
        print(f"   • Simple patterns: {simple_patterns}")
        print(f"   • Moderate patterns: {moderate_patterns}")
        print(f"   • Complex patterns: {complex_patterns}")
    
    def _save_filtering_workflow(self, workflow_filters: Dict[str, int]) -> None:
        """
        Save a filtering workflow to the chemspace database.
        
        Args:
            workflow_filters (Dict[str, int]): Dictionary mapping filter names to required instances
        """
        try:
            if not workflow_filters:
                print("⚠️  No workflow filters to save")
                return
            
            # Get workflow name from user
            workflow_name = input("📝 Enter a name for this filtering workflow: ").strip()
            if not workflow_name:
                workflow_name = f"filtering_workflow_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                print(f"📋 Using default name: {workflow_name}")
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create filtering_workflows table if it doesn't exist
            create_table_query = """
            CREATE TABLE IF NOT EXISTS filtering_workflows (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                workflow_name TEXT NOT NULL UNIQUE,
                filters_dict TEXT NOT NULL,
                creation_date TEXT NOT NULL,
                description TEXT,
                filter_count INTEGER DEFAULT 0,
                total_instances INTEGER DEFAULT 0
            )
            """
            
            cursor.execute(create_table_query)
            
            # Create indexes for faster searches
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_filtering_workflows_name ON filtering_workflows(workflow_name)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_filtering_workflows_date ON filtering_workflows(creation_date)")
            
            # Get current timestamp
            creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            # Get optional description from user
            description = input("📄 Enter a description for this filtering workflow (optional): ").strip()
            if not description:
                filter_count = len(workflow_filters)
                total_instances = sum(workflow_filters.values())
                description = f"Filtering workflow with {filter_count} filters and {total_instances} total instances"
            
            # Check if workflow name already exists
            cursor.execute("SELECT COUNT(*) FROM filtering_workflows WHERE workflow_name = ?", (workflow_name,))
            if cursor.fetchone()[0] > 0:
                overwrite = input(f"⚠️  Filtering workflow '{workflow_name}' already exists. Overwrite? (y/n): ").strip().lower()
                if overwrite in ['y', 'yes']:
                    # Delete existing workflow
                    cursor.execute("DELETE FROM filtering_workflows WHERE workflow_name = ?", (workflow_name,))
                    print(f"🔄 Overwriting existing filtering workflow '{workflow_name}'")
                else:
                    # Generate unique name
                    counter = 1
                    original_name = workflow_name
                    while True:
                        new_name = f"{original_name}_{counter}"
                        cursor.execute("SELECT COUNT(*) FROM filtering_workflows WHERE workflow_name = ?", (new_name,))
                        if cursor.fetchone()[0] == 0:
                            workflow_name = new_name
                            break
                        counter += 1
                    print(f"📝 Using name: {workflow_name}")
            
            # Calculate statistics
            filter_count = len(workflow_filters)
            total_instances = sum(workflow_filters.values())
            
            # Convert workflow to JSON string for storage
            filters_json = json.dumps(workflow_filters, sort_keys=True)
            
            # Insert workflow
            insert_query = """
            INSERT INTO filtering_workflows 
            (workflow_name, filters_dict, creation_date, description, filter_count, total_instances)
            VALUES (?, ?, ?, ?, ?, ?)
            """
            
            cursor.execute(insert_query, (
                workflow_name, 
                filters_json, 
                creation_date, 
                description,
                filter_count,
                total_instances
            ))
            
            conn.commit()
            conn.close()
            
            print(f"✅ Successfully saved filtering workflow!")
            print(f"   📋 Workflow name: '{workflow_name}'")
            print(f"   🔍 Filters saved: {filter_count}")
            print(f"   🔢 Total instances: {total_instances}")
            print(f"   📅 Created: {creation_date}")
            print(f"   📄 Description: {description}")
            
        except sqlite3.IntegrityError as e:
            print(f"❌ Database integrity error saving filtering workflow: {e}")
        except Exception as e:
            print(f"❌ Error saving filtering workflow: {e}")
    
    def list_filtering_workflows(self) -> None:
        """
        Display all saved filtering workflows in a formatted table.
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Check if filtering_workflows table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='filtering_workflows'")
            if not cursor.fetchone():
                print("📝 No filtering workflows table found")
                conn.close()
                return
            
            # Get all workflows with summary information
            cursor.execute("""
            SELECT workflow_name, creation_date, description, filters_dict, filter_count, total_instances
            FROM filtering_workflows 
            ORDER BY creation_date DESC
            """)
            
            workflows = cursor.fetchall()
            conn.close()
            
            if not workflows:
                print("📝 No saved filtering workflows found")
                return
            
            print("\n" + "="*100)
            print(f"SAVED FILTERING WORKFLOWS - Project: {self.name}")
            print("="*100)
            
            for i, (name, date, desc, filters_dict_str, filter_count, total_instances) in enumerate(workflows, 1):
                try:
                    filters_dict = json.loads(filters_dict_str)
                    actual_filter_count = len(filters_dict) if filters_dict else filter_count
                    
                    # Get filter names for summary
                    filter_names = list(filters_dict.keys()) if filters_dict else []
                    filters_summary = ', '.join(filter_names[:3])
                    if len(filter_names) > 3:
                        filters_summary += f" and {len(filter_names) - 3} more..."
                    elif not filters_summary:
                        filters_summary = "No filters available"
                    
                except json.JSONDecodeError:
                    actual_filter_count = filter_count or 0
                    filters_summary = "Error parsing filters"
                
                print(f"\n🔍 Workflow {i}: '{name}'")
                print(f"   📅 Created: {date}")
                print(f"   🔬 Filters: {actual_filter_count}")
                print(f"   🔢 Total instances: {total_instances or 0}")
                print(f"   📄 Description: {desc}")
                print(f"   🔍 Filters: {filters_summary}")
                print(f"   🔍 Filters dictionary: {filters_dict}")
            
            print("="*100)
            
        except Exception as e:
            print(f"❌ Error listing filtering workflows: {e}")

    def check_duplicates(self, table_name: str, 
                    duplicate_by: str = 'smiles',
                    show_duplicates: bool = True,
                    remove_duplicates: bool = False,
                    keep_first: bool = True,
                    export_duplicates: bool = False,
                    export_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Check for duplicate molecules in a table based on specified criteria.
        
        Args:
            table_name (str): Name of the table to check for duplicates
            duplicate_by (str): Column to check for duplicates ('smiles', 'inchi_key', 'name', or 'all')
            show_duplicates (bool): Whether to display duplicate entries
            remove_duplicates (bool): Whether to remove duplicate entries from the table
            keep_first (bool): If removing duplicates, whether to keep first occurrence (True) or last (False)
            export_duplicates (bool): Whether to export duplicate entries to CSV
            export_path (Optional[str]): Path for export file. If None, uses default naming
            
        Returns:
            Dict[str, Any]: Results dictionary with duplicate statistics and details
        """
        try:
            # Check if table exists
            tables = self.get_all_tables()
            if table_name not in tables:
                print(f"❌ Table '{table_name}' does not exist in chemspace database")
                return {'success': False, 'message': f"Table '{table_name}' not found"}
            
            # Get table data
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return {'success': False, 'message': f"No compounds in table '{table_name}'"}
            
            print(f"🔍 Checking for duplicates in table '{table_name}'")
            print(f"   📊 Total compounds: {len(compounds_df):,}")
            print(f"   🎯 Duplicate criteria: {duplicate_by}")
            
            # Initialize results
            results = {
                'success': True,
                'table_name': table_name,
                'total_compounds': len(compounds_df),
                'duplicate_by': duplicate_by,
                'duplicates_found': {},
                'unique_compounds': 0,
                'removed_duplicates': 0,
                'export_path': None
            }
            
            # Check duplicates based on specified criteria
            if duplicate_by == 'all':
                # Check all available columns for duplicates
                duplicate_results = self._check_all_duplicate_types(compounds_df)
                results['duplicates_found'] = duplicate_results
                
                # Display comprehensive results
                self._display_comprehensive_duplicate_report(duplicate_results, len(compounds_df))
                
            else:
                # Check specific column for duplicates
                duplicate_info = self._check_column_duplicates(compounds_df, duplicate_by)
                results['duplicates_found'][duplicate_by] = duplicate_info
                
                # Display results for specific column
                self._display_column_duplicate_report(duplicate_by, duplicate_info, len(compounds_df))
                
                # Show duplicate entries if requested
                if show_duplicates and duplicate_info['duplicate_count'] > 0:
                    self._show_duplicate_entries(compounds_df, duplicate_by, duplicate_info)
                
                # Export duplicates if requested
                if export_duplicates and duplicate_info['duplicate_count'] > 0:
                    export_success, export_file = self._export_duplicates(
                        compounds_df, duplicate_by, duplicate_info, export_path, table_name
                    )
                    if export_success:
                        results['export_path'] = export_file
                
                # Remove duplicates if requested
                if remove_duplicates and duplicate_info['duplicate_count'] > 0:
                    removed_count = self._remove_duplicates_from_table(
                        table_name, duplicate_by, keep_first
                    )
                    results['removed_duplicates'] = removed_count
                    results['unique_compounds'] = len(compounds_df) - removed_count
            
            return results
            
        except Exception as e:
            print(f"❌ Error checking duplicates in table '{table_name}': {e}")
            return {'success': False, 'message': f"Error checking duplicates: {e}"}

    def _check_all_duplicate_types(self, df: pd.DataFrame) -> Dict[str, Dict]:
        """
        Check for duplicates in all available columns.
        
        Args:
            df (pd.DataFrame): DataFrame to check
            
        Returns:
            Dict[str, Dict]: Dictionary of duplicate information for each column
        """
        duplicate_results = {}
        
        # Define columns to check and their priorities
        columns_to_check = {
            'smiles': 'SMILES strings',
            'inchi_key': 'InChI keys', 
            'name': 'Compound names'
        }
        
        for column, description in columns_to_check.items():
            if column in df.columns:
                print(f"   🔍 Checking {description}...")
                duplicate_info = self._check_column_duplicates(df, column)
                duplicate_results[column] = duplicate_info
                duplicate_results[column]['description'] = description
        
        return duplicate_results

    def _check_column_duplicates(self, df: pd.DataFrame, column: str) -> Dict[str, Any]:
        """
        Check for duplicates in a specific column.
        
        Args:
            df (pd.DataFrame): DataFrame to check
            column (str): Column name to check for duplicates
            
        Returns:
            Dict[str, Any]: Dictionary with duplicate information
        """
        if column not in df.columns:
            return {
                'column_exists': False,
                'duplicate_count': 0,
                'unique_count': 0,
                'duplicate_groups': [],
                'duplicate_values': []
            }
        
        # Remove null/empty values for duplicate checking
        valid_data = df[df[column].notna() & (df[column] != '') & (df[column] != 'nan')]
        
        if valid_data.empty:
            return {
                'column_exists': True,
                'duplicate_count': 0,
                'unique_count': 0,
                'duplicate_groups': [],
                'duplicate_values': [],
                'null_count': len(df) - len(valid_data)
            }
        
        # Find duplicates
        value_counts = valid_data[column].value_counts()
        duplicated_values = value_counts[value_counts > 1]
        
        # Get duplicate groups with details
        duplicate_groups = []
        for value, count in duplicated_values.items():
            duplicate_rows = valid_data[valid_data[column] == value]
            group_info = {
                'value': value,
                'count': count,
                'indices': duplicate_rows.index.tolist(),
                'compounds': []
            }
            
            # Add compound details for each duplicate
            for _, row in duplicate_rows.iterrows():
                compound_info = {
                    'id': row.get('id', 'N/A'),
                    'name': row.get('name', 'N/A'),
                    'smiles': row.get('smiles', 'N/A'),
                    'inchi_key': row.get('inchi_key', 'N/A')
                }
                group_info['compounds'].append(compound_info)
            
            duplicate_groups.append(group_info)
        
        # Sort duplicate groups by count (highest first)
        duplicate_groups.sort(key=lambda x: x['count'], reverse=True)
        
        return {
            'column_exists': True,
            'duplicate_count': sum(duplicated_values) - len(duplicated_values),  # Total duplicate entries
            'unique_count': len(value_counts),
            'duplicate_groups_count': len(duplicated_values),
            'duplicate_groups': duplicate_groups,
            'duplicate_values': duplicated_values.index.tolist(),
            'null_count': len(df) - len(valid_data),
            'valid_entries': len(valid_data)
        }

    def _display_comprehensive_duplicate_report(self, duplicate_results: Dict, total_compounds: int) -> None:
        """
        Display comprehensive duplicate report for all columns.
        
        Args:
            duplicate_results (Dict): Results from checking all columns
            total_compounds (int): Total number of compounds
        """
        print(f"\n📊 COMPREHENSIVE DUPLICATE ANALYSIS")
        print("=" * 70)
        
        for column, info in duplicate_results.items():
            if not info.get('column_exists', False):
                print(f"\n❌ {info.get('description', column)}: Column not found")
                continue
            
            print(f"\n🔍 {info.get('description', column)}:")
            print(f"   📊 Valid entries: {info.get('valid_entries', 0):,}")
            print(f"   🔄 Duplicate entries: {info.get('duplicate_count', 0):,}")
            print(f"   📋 Unique values: {info.get('unique_count', 0):,}")
            print(f"   🗂️  Duplicate groups: {info.get('duplicate_groups_count', 0):,}")
            
            if info.get('null_count', 0) > 0:
                print(f"   ⚠️  Null/empty entries: {info.get('null_count', 0):,}")
            
            if info.get('duplicate_count', 0) > 0:
                duplicate_rate = (info['duplicate_count'] / info['valid_entries']) * 100
                print(f"   📈 Duplication rate: {duplicate_rate:.2f}%")
        
        print("=" * 70)

    def _display_column_duplicate_report(self, column: str, duplicate_info: Dict, total_compounds: int) -> None:
        """
        Display duplicate report for a specific column.
        
        Args:
            column (str): Column name
            duplicate_info (Dict): Duplicate information for the column
            total_compounds (int): Total number of compounds
        """
        if not duplicate_info.get('column_exists', False):
            print(f"❌ Column '{column}' not found in table")
            return
        
        print(f"\n📊 DUPLICATE ANALYSIS - {column.upper()}")
        print("=" * 50)
        print(f"📊 Total compounds: {total_compounds:,}")
        print(f"📋 Valid entries: {duplicate_info.get('valid_entries', 0):,}")
        print(f"🔄 Duplicate entries: {duplicate_info.get('duplicate_count', 0):,}")
        print(f"📝 Unique values: {duplicate_info.get('unique_count', 0):,}")
        print(f"🗂️  Duplicate groups: {duplicate_info.get('duplicate_groups_count', 0):,}")
        
        if duplicate_info.get('null_count', 0) > 0:
            print(f"⚠️  Null/empty entries: {duplicate_info.get('null_count', 0):,}")
        
        if duplicate_info.get('duplicate_count', 0) > 0:
            duplicate_rate = (duplicate_info['duplicate_count'] / duplicate_info['valid_entries']) * 100
            print(f"📈 Duplication rate: {duplicate_rate:.2f}%")
            
            # Show top duplicate groups
            top_groups = duplicate_info.get('duplicate_groups', [])[:5]
            if top_groups:
                print(f"\n🔝 Top duplicate groups:")
                for i, group in enumerate(top_groups, 1):
                    value_display = group['value'][:50] + "..." if len(str(group['value'])) > 50 else group['value']
                    print(f"   {i}. '{value_display}' - {group['count']} occurrences")
        else:
            print("✅ No duplicates found!")
        
        print("=" * 50)

    def _show_duplicate_entries(self, df: pd.DataFrame, column: str, duplicate_info: Dict) -> None:
        """
        Display detailed information about duplicate entries.
        
        Args:
            df (pd.DataFrame): DataFrame containing the data
            column (str): Column being checked
            duplicate_info (Dict): Duplicate information
        """
        duplicate_groups = duplicate_info.get('duplicate_groups', [])
        
        if not duplicate_groups:
            return
        
        print(f"\n📋 DUPLICATE ENTRIES DETAILS - {column.upper()}")
        print("=" * 80)
        
        # Show first 10 duplicate groups
        max_groups_to_show = min(10, len(duplicate_groups))
        
        for i, group in enumerate(duplicate_groups[:max_groups_to_show], 1):
            value_display = str(group['value'])[:60] + "..." if len(str(group['value'])) > 60 else str(group['value'])
            print(f"\n🔄 Group {i}: '{value_display}' ({group['count']} occurrences)")
            print("-" * 60)
            
            for j, compound in enumerate(group['compounds'], 1):
                print(f"   {j}. ID: {compound.get('id', 'N/A')} | "
                    f"Name: {str(compound.get('name', 'N/A'))[:20]} | "
                    f"SMILES: {str(compound.get('smiles', 'N/A'))[:30]}")
        
        if len(duplicate_groups) > max_groups_to_show:
            remaining = len(duplicate_groups) - max_groups_to_show
            print(f"\n... and {remaining} more duplicate groups")
        
        print("=" * 80)

    def _export_duplicates(self, df: pd.DataFrame, column: str, duplicate_info: Dict, 
                        export_path: Optional[str], table_name: str) -> Tuple[bool, Optional[str]]:
        """
        Export duplicate entries to a CSV file.
        
        Args:
            df (pd.DataFrame): DataFrame containing the data
            column (str): Column being checked
            duplicate_info (Dict): Duplicate information
            export_path (Optional[str]): Custom export path
            table_name (str): Name of the source table
            
        Returns:
            Tuple[bool, Optional[str]]: (success, export_file_path)
        """
        try:
            if duplicate_info.get('duplicate_count', 0) == 0:
                print("⚠️  No duplicates to export")
                return False, None
            
            # Generate export filename
            if export_path is None:
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                export_filename = f"{table_name}_duplicates_{column}_{timestamp}.csv"
                export_path = os.path.join(self.path, 'chemspace', 'reports', export_filename)
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(export_path), exist_ok=True)
            
            # Collect all duplicate entries
            duplicate_rows = []
            
            for group in duplicate_info.get('duplicate_groups', []):
                for idx in group['indices']:
                    row_data = df.loc[idx].to_dict()
                    row_data['duplicate_group'] = group['value']
                    row_data['duplicate_count'] = group['count']
                    duplicate_rows.append(row_data)
            
            # Create DataFrame and export
            duplicates_df = pd.DataFrame(duplicate_rows)
            duplicates_df.to_csv(export_path, index=False)
            
            print(f"📤 Duplicates exported to: {export_path}")
            print(f"   📊 Exported {len(duplicate_rows)} duplicate entries")
            
            return True, export_path
            
        except Exception as e:
            print(f"❌ Error exporting duplicates: {e}")
            return False, None

    def _remove_duplicates_from_table(self, table_name: str, column: str, keep_first: bool = True) -> int:
        """
        Remove duplicate entries from the database table.
        
        Args:
            table_name (str): Name of the table
            column (str): Column to use for duplicate detection
            keep_first (bool): Whether to keep first occurrence (True) or last (False)
            
        Returns:
            int: Number of duplicate entries removed
        """
        try:
            print(f"\n🗑️  Removing duplicates from table '{table_name}' based on '{column}'...")
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get duplicate entries to remove
            if keep_first:
                # Keep first occurrence, remove subsequent ones
                delete_query = f"""
                DELETE FROM {table_name} 
                WHERE id NOT IN (
                    SELECT MIN(id) 
                    FROM {table_name} 
                    WHERE {column} IS NOT NULL AND {column} != '' AND {column} != 'nan'
                    GROUP BY {column}
                )
                AND {column} IS NOT NULL AND {column} != '' AND {column} != 'nan'
                """
            else:
                # Keep last occurrence, remove previous ones
                delete_query = f"""
                DELETE FROM {table_name} 
                WHERE id NOT IN (
                    SELECT MAX(id) 
                    FROM {table_name} 
                    WHERE {column} IS NOT NULL AND {column} != '' AND {column} != 'nan'
                    GROUP BY {column}
                )
                AND {column} IS NOT NULL AND {column} != '' AND {column} != 'nan'
                """
            
            # Execute deletion
            cursor.execute(delete_query)
            removed_count = cursor.rowcount
            
            conn.commit()
            conn.close()
            
            print(f"✅ Removed {removed_count} duplicate entries")
            print(f"   📋 Strategy: Keep {'first' if keep_first else 'last'} occurrence")
            
            return removed_count
            
        except Exception as e:
            print(f"❌ Error removing duplicates: {e}")
            return 0

    def find_similar_molecules(self, table_name: str, 
                            similarity_threshold: float = 0.8,
                            fingerprint_type: str = 'morgan',
                            max_comparisons: Optional[int] = None) -> pd.DataFrame:
        """
        Find chemically similar molecules in a table using molecular fingerprints.
        
        Args:
            table_name (str): Name of the table to analyze
            similarity_threshold (float): Minimum Tanimoto similarity (0.0 to 1.0)
            fingerprint_type (str): Type of fingerprint ('morgan', 'rdkit', 'maccs')
            max_comparisons (Optional[int]): Maximum number of molecules to compare (for performance)
            
        Returns:
            pd.DataFrame: DataFrame with similar molecule pairs
        """
        try:
            # Import required libraries
            try:
                from rdkit import Chem
                from rdkit.Chem import rdMolDescriptors, DataStructs
                from rdkit import RDLogger
                RDLogger.DisableLog('rdApp.*')
            except ImportError:
                print("❌ RDKit not installed. Cannot perform similarity analysis.")
                return pd.DataFrame()
            
            # Get compounds from table
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return pd.DataFrame()
            
            # Limit for performance if specified
            if max_comparisons and len(compounds_df) > max_comparisons:
                compounds_df = compounds_df.head(max_comparisons)
                print(f"🔍 Limited analysis to first {max_comparisons} molecules for performance")
            
            print(f"🧬 Analyzing molecular similarity in table '{table_name}'")
            print(f"   📊 Compounds: {len(compounds_df):,}")
            print(f"   🎯 Similarity threshold: {similarity_threshold:.2f}")
            print(f"   🔬 Fingerprint type: {fingerprint_type}")
            
            # Generate molecular fingerprints
            print("🔬 Computing molecular fingerprints...")
            fingerprints = []
            valid_molecules = []
            
            for idx, row in compounds_df.iterrows():
                try:
                    mol = Chem.MolFromSmiles(row['smiles'])
                    if mol is not None:
                        # Generate fingerprint based on type
                        if fingerprint_type == 'morgan':
                            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                        elif fingerprint_type == 'rdkit':
                            fp = Chem.RDKFingerprint(mol)
                        elif fingerprint_type == 'maccs':
                            fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
                        else:
                            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                        
                        fingerprints.append(fp)
                        valid_molecules.append(row)
                except Exception:
                    continue
            
            print(f"   ✅ Generated {len(fingerprints)} valid fingerprints")
            
            # Find similar pairs
            print("🔍 Finding similar molecule pairs...")
            similar_pairs = []
            
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(fingerprints) * (len(fingerprints) - 1) // 2,
                    desc="Comparing molecules",
                    unit="comparisons"
                )
            else:
                progress_bar = None
                comparisons_done = 0
                total_comparisons = len(fingerprints) * (len(fingerprints) - 1) // 2
            
            for i in range(len(fingerprints)):
                for j in range(i + 1, len(fingerprints)):
                    # Calculate Tanimoto similarity
                    similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                    
                    if similarity >= similarity_threshold:
                        mol1 = valid_molecules[i]
                        mol2 = valid_molecules[j]
                        
                        similar_pairs.append({
                            'molecule1_id': mol1.get('id', 'N/A'),
                            'molecule1_name': mol1.get('name', 'N/A'),
                            'molecule1_smiles': mol1['smiles'],
                            'molecule2_id': mol2.get('id', 'N/A'),
                            'molecule2_name': mol2.get('name', 'N/A'),
                            'molecule2_smiles': mol2['smiles'],
                            'similarity': similarity,
                            'fingerprint_type': fingerprint_type
                        })
                    
                    if progress_bar:
                        progress_bar.update(1)
                    else:
                        comparisons_done += 1
                        if comparisons_done % max(1, total_comparisons // 20) == 0:
                            progress = (comparisons_done / total_comparisons) * 100
                            print(f"   📊 Progress: {progress:.1f}%")
            
            if progress_bar:
                progress_bar.close()
            
            # Create results DataFrame
            results_df = pd.DataFrame(similar_pairs)
            
            if not results_df.empty:
                # Sort by similarity (highest first)
                results_df = results_df.sort_values('similarity', ascending=False)
                
                print(f"\n✅ Similarity analysis completed!")
                print(f"   🎯 Similar pairs found: {len(results_df):,}")
                print(f"   📈 Highest similarity: {results_df['similarity'].max():.4f}")
                print(f"   📉 Lowest similarity: {results_df['similarity'].min():.4f}")
                
                # Show top 5 most similar pairs
                if len(results_df) > 0:
                    print(f"\n🔝 Top 5 most similar pairs:")
                    for i, (_, pair) in enumerate(results_df.head(5).iterrows(), 1):
                        print(f"   {i}. {pair['molecule1_name']} ↔ {pair['molecule2_name']}")
                        print(f"      Similarity: {pair['similarity']:.4f}")
            else:
                print(f"\n📊 No similar pairs found above threshold {similarity_threshold:.2f}")
            
            return results_df
            
        except Exception as e:
            print(f"❌ Error in similarity analysis: {e}")
            return pd.DataFrame()

    def generate_duplicate_report(self, table_name: str, output_path: Optional[str] = None) -> bool:
        """
        Generate a comprehensive duplicate analysis report.
        
        Args:
            table_name (str): Name of the table to analyze
            output_path (Optional[str]): Path for the report file
            
        Returns:
            bool: True if report was generated successfully
        """
        try:
            # Run comprehensive duplicate check
            results = self.check_duplicates(table_name, duplicate_by='all', show_duplicates=False)
            
            if not results['success']:
                return False
            
            # Generate report filename
            if output_path is None:
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                report_filename = f"{table_name}_duplicate_report_{timestamp}.txt"
                output_path = os.path.join(self.path, 'chemspace', 'reports', report_filename)
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # Write comprehensive report
            with open(output_path, 'w') as f:
                f.write("MOLECULAR DUPLICATE ANALYSIS REPORT\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Project: {self.name}\n")
                f.write(f"Table: {table_name}\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total Compounds: {results['total_compounds']:,}\n\n")
                
                # Summary section
                f.write("DUPLICATE SUMMARY\n")
                f.write("-" * 20 + "\n")
                
                for column, info in results['duplicates_found'].items():
                    if info.get('column_exists', False):
                        f.write(f"\n{info.get('description', column)}:\n")
                        f.write(f"  Valid entries: {info.get('valid_entries', 0):,}\n")
                        f.write(f"  Duplicate entries: {info.get('duplicate_count', 0):,}\n")
                        f.write(f"  Unique values: {info.get('unique_count', 0):,}\n")
                        f.write(f"  Duplicate groups: {info.get('duplicate_groups_count', 0):,}\n")
                        
                        if info.get('duplicate_count', 0) > 0:
                            duplicate_rate = (info['duplicate_count'] / info['valid_entries']) * 100
                            f.write(f"  Duplication rate: {duplicate_rate:.2f}%\n")
                
                # Detailed duplicate groups
                f.write(f"\nDETAILED DUPLICATE GROUPS\n")
                f.write("-" * 30 + "\n")
                
                for column, info in results['duplicates_found'].items():
                    if info.get('column_exists', False) and info.get('duplicate_count', 0) > 0:
                        f.write(f"\n{column.upper()} DUPLICATES:\n")
                        
                        for i, group in enumerate(info.get('duplicate_groups', [])[:20], 1):
                            f.write(f"\nGroup {i}: '{group['value']}' ({group['count']} occurrences)\n")
                            for j, compound in enumerate(group['compounds'], 1):
                                f.write(f"  {j}. ID: {compound.get('id', 'N/A')} | "
                                    f"Name: {compound.get('name', 'N/A')} | "
                                    f"SMILES: {compound.get('smiles', 'N/A')}\n")
                        
                        if len(info.get('duplicate_groups', [])) > 20:
                            remaining = len(info['duplicate_groups']) - 20
                            f.write(f"\n... and {remaining} more duplicate groups\n")
            
            print(f"📋 Duplicate analysis report generated: {output_path}")
            return True
            
        except Exception as e:
            print(f"❌ Error generating duplicate report: {e}")
            return False
        
    def delete_reaction_workflow(self, workflow_identifier: Optional[str] = None) -> bool:
        """
        Delete a reaction workflow from the chemspace database.
        Can delete by workflow name or ID, or show interactive selection.
        
        Args:
            workflow_identifier (Optional[str]): Workflow name or ID to delete. If None, shows interactive selection.
            
        Returns:
            bool: True if workflow was deleted successfully, False otherwise
        """
        try:
            # Check if reaction_workflows table exists
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='reaction_workflows'")
            if not cursor.fetchone():
                print("📝 No reaction workflows table found - no workflows to delete")
                conn.close()
                return False
            
            # Get all available workflows
            cursor.execute("""
                SELECT id, workflow_name, creation_date, description, reactions_dict
                FROM reaction_workflows 
                ORDER BY creation_date DESC
            """)
            workflows = cursor.fetchall()
            
            if not workflows:
                print("📝 No reaction workflows found to delete")
                conn.close()
                return False
            
            # If no identifier provided, show interactive selection
            if workflow_identifier is None:
                workflow_to_delete = self._select_workflow_for_deletion(workflows)
                if workflow_to_delete is None:
                    conn.close()
                    return False
            else:
                # Find workflow by name or ID
                workflow_to_delete = self._find_workflow_by_identifier(workflows, workflow_identifier)
                if workflow_to_delete is None:
                    print(f"❌ No workflow found with identifier '{workflow_identifier}'")
                    self._show_available_reaction_workflows()
                    conn.close()
                    return False
            
            # Extract workflow information
            workflow_id, workflow_name, creation_date, description, reactions_dict_str = workflow_to_delete
            
            # Parse reactions for display
            try:
                reactions_dict = json.loads(reactions_dict_str)
                reaction_count = len(reactions_dict)
            except json.JSONDecodeError:
                reactions_dict = {}
                reaction_count = 0
            
            # Show workflow details and confirm deletion
            print(f"\n⚠️  WORKFLOW DELETION CONFIRMATION")
            print("=" * 60)
            print(f"🆔 ID: {workflow_id}")
            print(f"📋 Name: '{workflow_name}'")
            print(f"📅 Created: {creation_date}")
            print(f"🧪 Reactions: {reaction_count}")
            print(f"📄 Description: {description or 'No description'}")
            
            if reaction_count > 0:
                print(f"\n🧪 Reactions in this workflow:")
                for reaction_id, reaction_info in list(reactions_dict.items())[:5]:  # Show first 5
                    print(f"   • {reaction_info.get('name', 'Unknown')} (ID: {reaction_id})")
                if reaction_count > 5:
                    print(f"   ... and {reaction_count - 5} more reactions")
            
            print("=" * 60)
            
            # Double confirmation for safety
            confirm1 = input(f"❓ Are you sure you want to delete workflow '{workflow_name}'? (yes/no): ").strip().lower()
            if confirm1 not in ['yes', 'y']:
                print("❌ Workflow deletion cancelled")
                conn.close()
                return False
            
            confirm2 = input(f"❓ Type the workflow name '{workflow_name}' to confirm deletion: ").strip()
            if confirm2 != workflow_name:
                print("❌ Workflow name doesn't match. Deletion cancelled for safety.")
                conn.close()
                return False
            
            # Perform deletion
            print(f"🗑️  Deleting workflow '{workflow_name}'...")
            
            cursor.execute("DELETE FROM reaction_workflows WHERE id = ?", (workflow_id,))
            deleted_rows = cursor.rowcount
            
            conn.commit()
            conn.close()
            
            if deleted_rows > 0:
                print(f"✅ Successfully deleted reaction workflow!")
                print(f"   📋 Workflow: '{workflow_name}' (ID: {workflow_id})")
                print(f"   🧪 Reactions removed: {reaction_count}")
                print(f"   📅 Original creation date: {creation_date}")
                
                # Show remaining workflows
                remaining_count = self._count_remaining_workflows()
                if remaining_count > 0:
                    print(f"   📊 Remaining workflows: {remaining_count}")
                else:
                    print("   📝 No reaction workflows remaining in database")
                
                return True
            else:
                print(f"❌ Failed to delete workflow (no rows affected)")
                return False
            
        except Exception as e:
            print(f"❌ Error deleting reaction workflow: {e}")
            return False

    def _select_workflow_for_deletion(self, workflows: List[Tuple]) -> Optional[Tuple]:
        """
        Interactive workflow selection for deletion.
        
        Args:
            workflows (List[Tuple]): List of workflow tuples from database
            
        Returns:
            Optional[Tuple]: Selected workflow tuple or None if cancelled
        """
        try:
            print(f"\n🗑️  SELECT WORKFLOW FOR DELETION")
            print("=" * 60)
            print(f"📋 Available Reaction Workflows ({len(workflows)} total):")
            print("-" * 60)
            print(f"{'#':<3} {'Name':<25} {'Created':<12} {'Reactions':<10}")
            print("-" * 60)
            
            workflow_map = {}
            for i, (wf_id, name, date, desc, reactions_dict_str) in enumerate(workflows, 1):
                try:
                    reactions_dict = json.loads(reactions_dict_str)
                    reaction_count = len(reactions_dict)
                except json.JSONDecodeError:
                    reaction_count = 0
                
                date_short = date[:10] if date else "Unknown"
                name_display = name[:24] if len(name) <= 24 else name[:21] + "..."
                
                print(f"{i:<3} {name_display:<25} {date_short:<12} {reaction_count:<10}")
                workflow_map[str(i)] = workflows[i-1]
            
            print("-" * 60)
            print("Commands: Enter workflow number, 'cancel' to abort, 'list' to show details")
            
            while True:
                selection = input("\n🔍 Select workflow to delete: ").strip().lower()
                
                if selection in ['cancel', 'abort', 'exit', 'quit']:
                    print("❌ Workflow deletion cancelled by user")
                    return None
                elif selection == 'list':
                    self._show_detailed_workflow_list(workflows)
                    continue
                elif selection in workflow_map:
                    return workflow_map[selection]
                else:
                    try:
                        selection_int = int(selection)
                        if 1 <= selection_int <= len(workflows):
                            return workflows[selection_int - 1]
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(workflows)}")
                    except ValueError:
                        print("❌ Invalid input. Please enter a number, 'list', or 'cancel'")
            
        except KeyboardInterrupt:
            print("\n❌ Selection cancelled by user")
            return None
        except Exception as e:
            print(f"❌ Error in workflow selection: {e}")
            return None

    def _find_workflow_by_identifier(self, workflows: List[Tuple], identifier: str) -> Optional[Tuple]:
        """
        Find workflow by name or ID.
        
        Args:
            workflows (List[Tuple]): List of workflow tuples
            identifier (str): Workflow name or ID to search for
            
        Returns:
            Optional[Tuple]: Matching workflow tuple or None if not found
        """
        try:
            # Try to find by ID first
            try:
                search_id = int(identifier)
                for workflow in workflows:
                    if workflow[0] == search_id:  # workflow[0] is the ID
                        return workflow
            except ValueError:
                pass  # Not a number, continue to name search
            
            # Search by name (case-insensitive)
            identifier_lower = identifier.lower()
            for workflow in workflows:
                workflow_name = workflow[1].lower()  # workflow[1] is the name
                if workflow_name == identifier_lower:
                    return workflow
            
            # Partial name match if exact match not found
            for workflow in workflows:
                workflow_name = workflow[1].lower()
                if identifier_lower in workflow_name:
                    return workflow
            
            return None
            
        except Exception as e:
            print(f"❌ Error searching for workflow: {e}")
            return None

    def _show_detailed_workflow_list(self, workflows: List[Tuple]) -> None:
        """
        Show detailed information about available workflows.
        
        Args:
            workflows (List[Tuple]): List of workflow tuples
        """
        try:
            print(f"\n📋 DETAILED WORKFLOW INFORMATION")
            print("=" * 80)
            
            for i, (wf_id, name, date, desc, reactions_dict_str) in enumerate(workflows, 1):
                try:
                    reactions_dict = json.loads(reactions_dict_str)
                    reaction_count = len(reactions_dict)
                    
                    # Get reaction names for preview
                    reaction_names = [info.get('name', 'Unknown') for info in reactions_dict.values()]
                    reactions_preview = ', '.join(reaction_names[:3])
                    if len(reaction_names) > 3:
                        reactions_preview += f" and {len(reaction_names) - 3} more..."
                except json.JSONDecodeError:
                    reaction_count = 0
                    reactions_preview = "Error parsing reactions"
                
                print(f"\n{i}. {name} (ID: {wf_id})")
                print(f"   📅 Created: {date}")
                print(f"   🧪 Reactions: {reaction_count}")
                print(f"   📄 Description: {desc or 'No description'}")
                print(f"   🔬 Reactions: {reactions_preview}")
            
            print("=" * 80)
            
        except Exception as e:
            print(f"❌ Error showing detailed workflow list: {e}")

    def _count_remaining_workflows(self) -> int:
        """
        Count remaining workflows after deletion.
        
        Returns:
            int: Number of remaining workflows
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT COUNT(*) FROM reaction_workflows")
            count = cursor.fetchone()[0]
            
            conn.close()
            return count
            
        except Exception:
            return 0

    def delete_multiple_reaction_workflows(self, workflow_identifiers: List[str]) -> Dict[str, bool]:
        """
        Delete multiple reaction workflows in batch.
        
        Args:
            workflow_identifiers (List[str]): List of workflow names or IDs to delete
            
        Returns:
            Dict[str, bool]: Dictionary mapping workflow identifiers to deletion success status
        """
        try:
            print(f"🗑️  BATCH WORKFLOW DELETION")
            print("=" * 50)
            print(f"📋 Workflows to delete: {len(workflow_identifiers)}")
            
            results = {}
            successful_deletions = 0
            failed_deletions = 0
            
            for identifier in workflow_identifiers:
                print(f"\n🔍 Processing workflow: '{identifier}'")
                success = self.delete_reaction_workflow(identifier)
                results[identifier] = success
                
                if success:
                    successful_deletions += 1
                else:
                    failed_deletions += 1
            
            # Summary
            print(f"\n📊 BATCH DELETION SUMMARY")
            print("=" * 30)
            print(f"✅ Successful deletions: {successful_deletions}")
            print(f"❌ Failed deletions: {failed_deletions}")
            print(f"📋 Total processed: {len(workflow_identifiers)}")
            
            return results
            
        except Exception as e:
            print(f"❌ Error in batch workflow deletion: {e}")
            return {identifier: False for identifier in workflow_identifiers}

    def clear_all_reaction_workflows(self, confirm_with_count: bool = True) -> bool:
        """
        Delete all reaction workflows from the database.
        
        Args:
            confirm_with_count (bool): Require user to enter the exact count for confirmation
            
        Returns:
            bool: True if all workflows were cleared successfully
        """
        try:
            # Check if table exists and get workflow count
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='reaction_workflows'")
            if not cursor.fetchone():
                print("📝 No reaction workflows table found")
                conn.close()
                return True
            
            cursor.execute("SELECT COUNT(*) FROM reaction_workflows")
            workflow_count = cursor.fetchone()[0]
            
            if workflow_count == 0:
                print("📝 No reaction workflows found to clear")
                conn.close()
                return True
            
            # Show warning and get confirmation
            print(f"⚠️  CLEAR ALL REACTION WORKFLOWS")
            print("=" * 50)
            print(f"🚨 WARNING: This will permanently delete ALL {workflow_count} reaction workflows!")
            print("🚨 This action cannot be undone!")
            print("=" * 50)
            
            # First confirmation
            confirm1 = input(f"❓ Are you sure you want to delete all {workflow_count} workflows? (yes/no): ").strip().lower()
            if confirm1 not in ['yes', 'y']:
                print("❌ Operation cancelled")
                conn.close()
                return False
            
            # Second confirmation with count
            if confirm_with_count:
                confirm2 = input(f"❓ Type the number {workflow_count} to confirm deletion of all workflows: ").strip()
                if confirm2 != str(workflow_count):
                    print("❌ Count doesn't match. Operation cancelled for safety.")
                    conn.close()
                    return False
            
            # Third confirmation
            confirm3 = input("❓ Type 'DELETE ALL WORKFLOWS' in capitals to proceed: ").strip()
            if confirm3 != 'DELETE ALL WORKFLOWS':
                print("❌ Confirmation phrase doesn't match. Operation cancelled.")
                conn.close()
                return False
            
            # Perform deletion
            print(f"🗑️  Clearing all reaction workflows...")
            cursor.execute("DELETE FROM reaction_workflows")
            deleted_count = cursor.rowcount
            
            conn.commit()
            conn.close()
            
            print(f"✅ Successfully cleared all reaction workflows!")
            print(f"   🗑️  Workflows deleted: {deleted_count}")
            print(f"   📊 Database is now empty of reaction workflows")
            
            return True
            
        except Exception as e:
            print(f"❌ Error clearing all reaction workflows: {e}")
            return False
        
    def _analyze_reaction_types(self, reactions_dict: Dict[int, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Analyze the reaction workflow to identify bimolecular and unimolecular reactions.
        
        Args:
            reactions_dict (Dict): Dictionary of reactions in the workflow
            
        Returns:
            Dict[str, Any]: Analysis results including reaction types and details
        """
        try:
            analysis = {
                'has_bimolecular': False,
                'has_unimolecular': False,
                'bimolecular_reactions': [],
                'unimolecular_reactions': [],
                'invalid_reactions': [],
                'total_reactions': len(reactions_dict)
            }
            
            for reaction_id, reaction_info in reactions_dict.items():
                smarts = reaction_info.get('smarts', '')
                name = reaction_info.get('name', f'reaction_{reaction_id}')
                
                # Check if SMARTS contains reactant separator (dot)
                if '>>' in smarts:
                    reactants_part = smarts.split('>>')[0]
                    
                    # Count dots in reactants part to determine if bimolecular
                    dot_count = reactants_part.count('.')
                    
                    if dot_count >= 1:
                        # Bimolecular or multi-molecular reaction
                        analysis['has_bimolecular'] = True
                        analysis['bimolecular_reactions'].append({
                            'id': reaction_id,
                            'name': name,
                            'smarts': smarts,
                            'reactant_count': dot_count + 1,
                            'order': reaction_info.get('order', 0)
                        })
                    else:
                        # Unimolecular reaction
                        analysis['has_unimolecular'] = True
                        analysis['unimolecular_reactions'].append({
                            'id': reaction_id,
                            'name': name,
                            'smarts': smarts,
                            'order': reaction_info.get('order', 0)
                        })
                else:
                    # Invalid SMARTS format
                    analysis['invalid_reactions'].append({
                        'id': reaction_id,
                        'name': name,
                        'smarts': smarts,
                        'error': 'Missing reaction arrow (>>)'
                    })
            
            return analysis
            
        except Exception as e:
            print(f"❌ Error analyzing reaction types: {e}")
            return {
                'has_bimolecular': False,
                'has_unimolecular': False,
                'bimolecular_reactions': [],
                'unimolecular_reactions': [],
                'invalid_reactions': [],
                'total_reactions': 0
            }

    def _display_reaction_analysis(self, analysis: Dict[str, Any]) -> None:
        """
        Display the analysis results of reaction types in the workflow.
        
        Args:
            analysis (Dict): Analysis results from _analyze_reaction_types
        """
        print(f"\n🔬 REACTION WORKFLOW ANALYSIS")
        print("=" * 60)
        print(f"📊 Total reactions: {analysis['total_reactions']}")
        print(f"🧪 Unimolecular reactions: {len(analysis['unimolecular_reactions'])}")
        print(f"🧬 Bimolecular reactions: {len(analysis['bimolecular_reactions'])}")
        
        if analysis['invalid_reactions']:
            print(f"❌ Invalid reactions: {len(analysis['invalid_reactions'])}")
        
        # Display bimolecular reactions details
        if analysis['bimolecular_reactions']:
            print(f"\n🧬 BIMOLECULAR REACTIONS DETECTED:")
            print("-" * 40)
            for reaction in analysis['bimolecular_reactions']:
                print(f"   {reaction['order']}. {reaction['name']} (ID: {reaction['id']})")
                print(f"      🔗 Requires {reaction['reactant_count']} reactants")
                print(f"      🧪 SMARTS: {reaction['smarts']}")
        
        # Display unimolecular reactions
        if analysis['unimolecular_reactions']:
            print(f"\n🧪 UNIMOLECULAR REACTIONS:")
            print("-" * 25)
            for reaction in analysis['unimolecular_reactions']:
                print(f"   {reaction['order']}. {reaction['name']} (ID: {reaction['id']})")
        
        # Display invalid reactions
        if analysis['invalid_reactions']:
            print(f"\n❌ INVALID REACTIONS:")
            print("-" * 20)
            for reaction in analysis['invalid_reactions']:
                print(f"   ⚠️  {reaction['name']} (ID: {reaction['id']})")
                print(f"      Error: {reaction['error']}")
                print(f"      SMARTS: {reaction['smarts']}")
        
        print("=" * 60)

    def _select_tables_for_bimolecular_reactions(self, available_tables: List[str], 
                                            analysis: Dict[str, Any]) -> Dict[str, Any]:
        """
        Select tables for workflows containing bimolecular reactions.
        
        Args:
            available_tables (List[str]): List of available table names
            analysis (Dict): Reaction analysis results
            
        Returns:
            Dict[str, Any]: Selected tables configuration for bimolecular reactions
        """
        try:
            print(f"\n🧬 BIMOLECULAR REACTION WORKFLOW SETUP")
            print("=" * 50)
            print("This workflow contains bimolecular reactions that require two reactant tables.")
            print("You need to specify:")
            print("   • Primary reactant table (Table A)")
            print("   • Secondary reactant table (Table B)")
            print("   • How to handle unimolecular reactions (if present)")
            
            # Select primary reactant table
            print(f"\n📋 Select PRIMARY reactant table (Table A):")
            primary_table = self._select_single_table(available_tables, "primary reactant")
            if not primary_table:
                return {}
            
            # Select secondary reactant table
            print(f"\n📋 Select SECONDARY reactant table (Table B):")
            remaining_tables = [t for t in available_tables if t != primary_table]
            secondary_table = self._select_single_table(remaining_tables, "secondary reactant")
            if not secondary_table:
                return {}
            
            # Handle unimolecular reactions if present
            unimolecular_strategy = 'primary'  # Default
            if analysis['has_unimolecular']:
                print(f"\n🧪 This workflow also contains {len(analysis['unimolecular_reactions'])} unimolecular reaction(s).")
                print("Which table should be used for unimolecular reactions?")
                print("   1. Primary table (Table A)")
                print("   2. Secondary table (Table B)")
                print("   3. Both tables separately")
                
                choice = input("Enter choice (1-3): ").strip()
                if choice == '2':
                    unimolecular_strategy = 'secondary'
                elif choice == '3':
                    unimolecular_strategy = 'both'
                else:
                    unimolecular_strategy = 'primary'
            
            # Validate table selections
            primary_count = self.get_compound_count(table_name=primary_table)
            secondary_count = self.get_compound_count(table_name=secondary_table)
            
            print(f"\n✅ Table Selection Summary:")
            print(f"   🅰️  Primary table: '{primary_table}' ({primary_count} compounds)")
            print(f"   🅱️  Secondary table: '{secondary_table}' ({secondary_count} compounds)")
            print(f"   🧪 Unimolecular strategy: {unimolecular_strategy}")
            
            # Estimate computational complexity
            if analysis['has_bimolecular']:
                combinations = primary_count * secondary_count
                print(f"   🔢 Potential combinations: {combinations:,}")
                
                if combinations > 1000000:
                    print(f"   ⚠️  WARNING: Large number of combinations may take significant time!")
                    proceed = input("   Continue anyway? (y/n): ").strip().lower()
                    if proceed not in ['y', 'yes']:
                        return {}
            
            return {
                'type': 'bimolecular',
                'primary_table': primary_table,
                'secondary_table': secondary_table,
                'unimolecular_strategy': unimolecular_strategy,
                'primary_count': primary_count,
                'secondary_count': secondary_count
            }
            
        except Exception as e:
            print(f"❌ Error selecting tables for bimolecular reactions: {e}")
            return {}

    def _select_tables_for_unimolecular_reactions(self, available_tables: List[str]) -> Dict[str, Any]:
        """
        Select tables for workflows containing only unimolecular reactions.
        
        Args:
            available_tables (List[str]): List of available table names
            
        Returns:
            Dict[str, Any]: Selected tables configuration for unimolecular reactions
        """
        try:
            print(f"\n🧪 UNIMOLECULAR REACTION WORKFLOW SETUP")
            print("=" * 50)
            print("This workflow contains only unimolecular reactions.")
            print("Select one or more tables to process:")
            
            selected_tables = self._select_multiple_tables(available_tables)
            
            if selected_tables:
                total_compounds = sum(self.get_compound_count(table_name=table) for table in selected_tables)
                print(f"\n✅ Selected {len(selected_tables)} table(s) ({total_compounds} total compounds)")
                
                return {
                    'type': 'unimolecular',
                    'tables': selected_tables,
                    'total_compounds': total_compounds
                }
            else:
                return {}
            
        except Exception as e:
            print(f"❌ Error selecting tables for unimolecular reactions: {e}")
            return {}

    def _select_single_table(self, available_tables: List[str], table_type: str) -> Optional[str]:
        """
        Select a single table from available options.
        
        Args:
            available_tables (List[str]): List of available table names
            table_type (str): Description of table type for user prompt
            
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        try:
            if not available_tables:
                print(f"❌ No tables available for {table_type} selection")
                return None
            
            print(f"\nAvailable tables for {table_type}:")
            for i, table in enumerate(available_tables, 1):
                compound_count = self.get_compound_count(table_name=table)
                print(f"   {i}. {table} ({compound_count} compounds)")
            
            while True:
                try:
                    selection = input(f"\nSelect {table_type} table (1-{len(available_tables)} or table name): ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    
                    # Try as number first
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(available_tables):
                            selected_table = available_tables[table_idx]
                            print(f"✅ Selected {table_type} table: '{selected_table}'")
                            return selected_table
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(available_tables)}")
                            continue
                    except ValueError:
                        # Treat as table name
                        if selection in available_tables:
                            print(f"✅ Selected {table_type} table: '{selection}'")
                            return selection
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print(f"\n❌ {table_type} table selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error selecting {table_type} table: {e}")
            return None

    def _select_multiple_tables(self, available_tables: List[str]) -> List[str]:
        """
        Select multiple tables from available options.
        
        Args:
            available_tables (List[str]): List of available table names
            
        Returns:
            List[str]: List of selected table names
        """
        try:
            print("Enter table numbers or names (comma-separated), or 'all' for all tables:")
            
            user_input = input("Selection: ").strip().lower()
            
            if not user_input or user_input in ['cancel', 'quit', 'exit']:
                return []
            
            if user_input == 'all':
                print(f"✅ Selected all {len(available_tables)} tables")
                return available_tables
            
            selected_tables = []
            selections = [s.strip() for s in user_input.split(',')]
            
            for selection in selections:
                try:
                    # Try as table number first
                    table_idx = int(selection) - 1
                    if 0 <= table_idx < len(available_tables):
                        table_name = available_tables[table_idx]
                        if table_name not in selected_tables:
                            selected_tables.append(table_name)
                    else:
                        print(f"⚠️  Invalid table number: {selection}")
                except ValueError:
                    # Treat as table name
                    if selection in available_tables:
                        if selection not in selected_tables:
                            selected_tables.append(selection)
                    else:
                        print(f"⚠️  Table '{selection}' not found")
            
            return selected_tables
            
        except Exception as e:
            print(f"❌ Error selecting multiple tables: {e}")
            return []

    def _process_bimolecular_workflow(self, table_config: Dict[str, Any], 
                                    reactions_dict: Dict[int, Dict[str, Any]], 
                                    workflow_name: str, 
                                    analysis: Dict[str, Any]) -> None:
        """
        Process a workflow containing bimolecular reactions.
        
        Args:
            table_config (Dict): Table configuration for bimolecular reactions
            reactions_dict (Dict): Dictionary of reactions
            workflow_name (str): Name of the workflow
            analysis (Dict): Reaction analysis results
        """
        try:
            print(f"\n🧬 PROCESSING BIMOLECULAR REACTION WORKFLOW")
            print("=" * 50)
            
            primary_table = table_config['primary_table']
            secondary_table = table_config['secondary_table']
            
            # Get compound data
            primary_df = self._get_table_as_dataframe(primary_table)
            secondary_df = self._get_table_as_dataframe(secondary_table)
            
            print(f"📊 Primary reactants: {len(primary_df)} compounds")
            print(f"📊 Secondary reactants: {len(secondary_df)} compounds")
            
            # Process bimolecular reactions
            bimolecular_products = []
            for reaction_info in analysis['bimolecular_reactions']:
                print(f"\n🔬 Processing bimolecular reaction: {reaction_info['name']}")
                
                reaction_products = self._apply_bimolecular_reaction(
                    primary_df, secondary_df, reaction_info, workflow_name
                )
                bimolecular_products.extend(reaction_products)
            
            # Process unimolecular reactions if present
            unimolecular_products = []
            if analysis['has_unimolecular']:
                strategy = table_config['unimolecular_strategy']
                
                for reaction_info in analysis['unimolecular_reactions']:
                    print(f"\n🧪 Processing unimolecular reaction: {reaction_info['name']}")
                    
                    if strategy == 'both':
                        # Apply to both tables separately
                        products_a = self._apply_unimolecular_reaction(
                            primary_df, reaction_info, workflow_name, f"{primary_table}_"
                        )
                        products_b = self._apply_unimolecular_reaction(
                            secondary_df, reaction_info, workflow_name, f"{secondary_table}_"
                        )
                        unimolecular_products.extend(products_a)
                        unimolecular_products.extend(products_b)
                    elif strategy == 'secondary':
                        products = self._apply_unimolecular_reaction(
                            secondary_df, reaction_info, workflow_name, f"{secondary_table}_"
                        )
                        unimolecular_products.extend(products)
                    else:  # primary
                        products = self._apply_unimolecular_reaction(
                            primary_df, reaction_info, workflow_name, f"{primary_table}_"
                        )
                        unimolecular_products.extend(products)
            
            # Combine and save all products
            all_products = bimolecular_products + unimolecular_products
            
            print(f"\n📊 REACTION WORKFLOW RESULTS:")
            print(f"   🧬 Bimolecular products: {len(bimolecular_products)}")
            print(f"   🧪 Unimolecular products: {len(unimolecular_products)}")
            print(f"   🎯 Total products: {len(all_products)}")
            
            if all_products:
                save_choice = input("\n💾 Save all products to a new table? (y/n): ").strip().lower()
                if save_choice in ['y', 'yes']:
                    self._save_reaction_products(
                        all_products, workflow_name, f"{primary_table}_and_{secondary_table}"
                    )
            
        except Exception as e:
            print(f"❌ Error processing bimolecular workflow: {e}")

    def _process_unimolecular_workflow(self, table_config: Dict[str, Any], 
                                    reactions_dict: Dict[int, Dict[str, Any]], 
                                    workflow_name: str) -> None:
        """
        Process a workflow containing only unimolecular reactions.
        
        Args:
            table_config (Dict): Table configuration for unimolecular reactions
            reactions_dict (Dict): Dictionary of reactions
            workflow_name (str): Name of the workflow
        """
        try:
            print(f"\n🧪 PROCESSING UNIMOLECULAR REACTION WORKFLOW")
            print("=" * 50)
            
            tables = table_config['tables']
            all_products = []
            
            # Process each table
            for table_name in tables:
                print(f"\n🔬 Processing table: '{table_name}'")
                compounds_df = self._get_table_as_dataframe(table_name)
                
                if compounds_df.empty:
                    print(f"   ⚠️  No compounds in table '{table_name}', skipping...")
                    continue
                
                print(f"   📊 Processing {len(compounds_df)} compounds...")
                
                # Apply each unimolecular reaction
                table_products = []
                sorted_reactions = sorted(reactions_dict.items(), key=lambda x: x[1]['order'])
                
                for reaction_id, reaction_info in sorted_reactions:
                    print(f"   🧪 Applying reaction: {reaction_info['name']}")
                    
                    reaction_products = self._apply_unimolecular_reaction(
                        compounds_df, reaction_info, workflow_name, f"{table_name}_"
                    )
                    table_products.extend(reaction_products)
                
                all_products.extend(table_products)
                print(f"   ✅ Generated {len(table_products)} products from '{table_name}'")
            
            # Final results
            print(f"\n📊 UNIMOLECULAR WORKFLOW RESULTS:")
            print(f"   📋 Tables processed: {len(tables)}")
            print(f"   🎯 Total products: {len(all_products)}")
            
            if all_products:
                save_choice = input("\n💾 Save all products to a new table? (y/n): ").strip().lower()
                if save_choice in ['y', 'yes']:
                    combined_source = "_".join(tables[:3])  # Limit length
                    if len(tables) > 3:
                        combined_source += f"_and_{len(tables)-3}_more"
                    
                    self._save_reaction_products(all_products, workflow_name, combined_source)
            
        except Exception as e:
            print(f"❌ Error processing unimolecular workflow: {e}")

    def _apply_bimolecular_reaction(self, primary_df: pd.DataFrame, secondary_df: pd.DataFrame,
                                reaction_info: Dict[str, Any], workflow_name: str) -> List[Dict[str, Any]]:
        """
        Apply a bimolecular reaction between compounds from two tables.
        
        Args:
            primary_df (pd.DataFrame): Primary reactants dataframe
            secondary_df (pd.DataFrame): Secondary reactants dataframe  
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            
        Returns:
            List[Dict]: List of reaction products
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            reaction_smarts = reaction_info['smarts']
            reaction_name = reaction_info['name']
            products = []
            
            # Parse reaction
            rxn = AllChem.ReactionFromSmarts(reaction_smarts)
            if rxn is None:
                print(f"   ❌ Invalid reaction SMARTS: {reaction_smarts}")
                return []
            
            # Get expected number of reactants
            expected_reactants = rxn.GetNumReactantTemplates()
            print(f"   🔗 Reaction expects {expected_reactants} reactants")
            
            # Process combinations with progress tracking
            total_combinations = len(primary_df) * len(secondary_df)
            print(f"   📊 Processing {total_combinations:,} reactant combinations...")
            
            processed = 0
            successful_reactions = 0
            
            # Use progress bar if available
            if TQDM_AVAILABLE and total_combinations > 1000:
                progress_bar = tqdm(
                    total=total_combinations,
                    desc=f"Bimolecular {reaction_name}",
                    unit="combinations"
                )
            else:
                progress_bar = None
            
            for _, primary_compound in primary_df.iterrows():
                primary_mol = Chem.MolFromSmiles(primary_compound['smiles'])
                if primary_mol is None:
                    if progress_bar:
                        progress_bar.update(len(secondary_df))
                    processed += len(secondary_df)
                    continue
                
                for _, secondary_compound in secondary_df.iterrows():
                    secondary_mol = Chem.MolFromSmiles(secondary_compound['smiles'])
                    if secondary_mol is None:
                        if progress_bar:
                            progress_bar.update(1)
                        processed += 1
                        continue
                    
                    try:
                        # Run bimolecular reaction
                        if expected_reactants >= 2:
                            reaction_results = rxn.RunReactants((primary_mol, secondary_mol))
                        else:
                            # Fallback for reactions that might work with single reactant
                            reaction_results = rxn.RunReactants((primary_mol,))
                        
                        # Process products
                        for product_set_idx, product_set in enumerate(reaction_results):
                            for product_idx, product_mol in enumerate(product_set):
                                try:
                                    Chem.SanitizeMol(product_mol)
                                    product_smiles = Chem.MolToSmiles(product_mol)
                                    
                                    # Generate product name
                                    primary_name = primary_compound.get('name', f"cpd_{primary_compound.get('id', 'unk')}")
                                    secondary_name = secondary_compound.get('name', f"cpd_{secondary_compound.get('id', 'unk')}")
                                    product_name = f"{primary_name}+{secondary_name}_{reaction_name}_{product_set_idx}_{product_idx}"
                                    
                                    products.append({
                                        'smiles': product_smiles,
                                        'name': product_name,
                                        'flag': 'bimolecular_product',
                                        'reactant1_name': primary_name,
                                        'reactant1_smiles': primary_compound['smiles'],
                                        'reactant2_name': secondary_name,
                                        'reactant2_smiles': secondary_compound['smiles'],
                                        'reaction_name': reaction_name,
                                        'workflow': workflow_name
                                    })
                                    successful_reactions += 1
                                    
                                except Exception:
                                    continue  # Skip invalid products
                                    
                    except Exception:
                        pass  # Skip failed reactions
                    
                    processed += 1
                    if progress_bar:
                        progress_bar.update(1)
                    elif processed % max(1, total_combinations // 20) == 0:
                        progress_pct = (processed / total_combinations) * 100
                        print(f"      📊 Progress: {progress_pct:.1f}% ({successful_reactions} products so far)")
            
            if progress_bar:
                progress_bar.close()
            
            print(f"   ✅ Generated {len(products)} products from {successful_reactions} successful reactions")
            return products
            
        except Exception as e:
            print(f"   ❌ Error in bimolecular reaction {reaction_info['name']}: {e}")
            return []

    def _apply_unimolecular_reaction(self, compounds_df: pd.DataFrame, reaction_info: Dict[str, Any],
                                workflow_name: str, name_prefix: str = "") -> List[Dict[str, Any]]:
        """
        Apply a unimolecular reaction to compounds from a table.
        
        Args:
            compounds_df (pd.DataFrame): Compounds dataframe
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            name_prefix (str): Prefix for product names
            
        Returns:
            List[Dict]: List of reaction products
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            reaction_smarts = reaction_info['smarts']
            reaction_name = reaction_info['name']
            products = []
            
            # Parse reaction
            rxn = AllChem.ReactionFromSmarts(reaction_smarts)
            if rxn is None:
                print(f"   ❌ Invalid reaction SMARTS: {reaction_smarts}")
                return []
            
            print(f"   📊 Processing {len(compounds_df)} compounds...")
            
            successful_reactions = 0
            
            # Use progress bar if available
            if TQDM_AVAILABLE and len(compounds_df) > 100:
                progress_bar = tqdm(
                    compounds_df.iterrows(),
                    total=len(compounds_df),
                    desc=f"Unimolecular {reaction_name}",
                    unit="compounds"
                )
            else:
                progress_bar = compounds_df.iterrows()
                processed_count = 0
            
            for _, compound in progress_bar:
                try:
                    reactant_mol = Chem.MolFromSmiles(compound['smiles'])
                    if reactant_mol is None:
                        continue
                    
                    # Run unimolecular reaction
                    reaction_results = rxn.RunReactants((reactant_mol,))
                    
                    # Process products
                    for product_set_idx, product_set in enumerate(reaction_results):
                        for product_idx, product_mol in enumerate(product_set):
                            try:
                                Chem.SanitizeMol(product_mol)
                                product_smiles = Chem.MolToSmiles(product_mol)
                                
                                # Generate product name
                                original_name = compound.get('name', f"cpd_{compound.get('id', 'unk')}")
                                product_name = f"{name_prefix}{original_name}_{reaction_name}_{product_set_idx}_{product_idx}"
                                
                                products.append({
                                    'smiles': product_smiles,
                                    'name': product_name,
                                    'flag': 'unimolecular_product',
                                    'reactant_name': original_name,
                                    'reactant_smiles': compound['smiles'],
                                    'reaction_name': reaction_name,
                                    'workflow': workflow_name
                                })
                                successful_reactions += 1
                                
                            except Exception:
                                continue  # Skip invalid products
                                
                except Exception:
                    continue  # Skip failed reactions
                
                # Progress update for non-tqdm case
                if not TQDM_AVAILABLE:
                    processed_count += 1
                    if processed_count % max(1, len(compounds_df) // 10) == 0:
                        print(f"      📊 Processed {processed_count}/{len(compounds_df)} compounds ({successful_reactions} products)")
            
            print(f"   ✅ Generated {len(products)} products from {successful_reactions} successful reactions")
            return products
            
        except Exception as e:
            print(f"   ❌ Error in unimolecular reaction {reaction_info['name']}: {e}")
            return []

    def apply_reaction_workflow(self, workflow_id: Optional[int] = None, 
                                        max_workers: Optional[int] = None,
                                        chunk_size: Optional[int] = None,
                                        parallel_threshold: int = 1000,
                                        memory_limit_mb: int = 500,
                                        stream_threshold: int = 10000) -> bool:
        """
        Apply a saved reaction workflow to compounds with streaming to disk to avoid memory exhaustion.
        Uses existing helper methods for consistency while adding disk streaming for large datasets.
        
        Args:
            workflow_id (Optional[int]): ID of the reaction workflow to apply. If None, prompts user for selection.
            max_workers (Optional[int]): Maximum number of parallel workers. If None, uses cpu_count()
            chunk_size (Optional[int]): Size of chunks for parallel processing. If None, automatically calculated
            parallel_threshold (int): Minimum number of compounds to trigger parallel processing
            memory_limit_mb (int): Memory limit in MB to trigger streaming mode
            stream_threshold (int): Minimum number of products to trigger streaming to disk
            
        Returns:
            bool: True if workflow was applied successfully
        """
        try:
            print(f"🚀 Starting streaming reaction workflow application...")
            
            # Reuse existing workflow selection logic
            if workflow_id is None:
                workflow_id = self._select_reaction_workflow_for_application()
                if workflow_id is None:
                    print("❌ No workflow selected for application")
                    return False
            
            print(f"   🆔 Workflow ID: {workflow_id}")
            
            # Reuse existing workflow loading logic
            workflow_data = self._load_reaction_workflow_by_id(workflow_id)
            if not workflow_data:
                print(f"❌ No reaction workflow found with ID {workflow_id}")
                return False
            
            workflow_name = workflow_data['workflow_name']
            reactions_dict = workflow_data['reactions_dict']
            
            print(f"   📋 Workflow: '{workflow_name}'")
            print(f"   🧪 Total reactions: {len(reactions_dict)}")
            print(f"   💾 Memory limit: {memory_limit_mb}MB")
            print(f"   📦 Stream threshold: {stream_threshold:,} products")
            
            # Set up parallel processing parameters
            if max_workers is None:
                max_workers = min(cpu_count() or 4, 8)
            
            print(f"   👥 Max workers: {max_workers}")
            print(f"   📦 Parallel threshold: {parallel_threshold:,} compounds")
            
            # Sort reactions by order for sequential processing
            sorted_reactions = sorted(reactions_dict.items(), key=lambda x: x[1]['order'])
            
            print(f"\n🔄 STREAMING REACTION WORKFLOW ANALYSIS")
            print("=" * 60)
            for i, (reaction_id, reaction_info) in enumerate(sorted_reactions, 1):
                reaction_type = self._analyze_single_reaction_type(reaction_info['smarts'])
                print(f"   Step {i}: {reaction_info['name']} ({reaction_type})")
            print("=" * 60)
            
            # Confirm execution
            input("Press Enter to continue with streaming workflow execution...")
            
            # Import required libraries
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
                from rdkit import RDLogger
                RDLogger.DisableLog('rdApp.*')  # Suppress RDKit warnings
            except ImportError:
                print("❌ RDKit not installed. Please install RDKit to use reaction workflows:")
                print("   conda install -c conda-forge rdkit")
                return False
            
            # Initialize streaming workflow state
            workflow_state = {
                'step_results': {},
                'available_tables': self.get_all_tables(),
                'step_products': {},
                'streaming_stats': {
                    'total_chunks_streamed': 0,
                    'total_temp_files': 0,
                    'total_disk_writes': 0,
                    'peak_memory_usage': 0
                }
            }
            
            if not workflow_state['available_tables']:
                print("❌ No tables available in chemspace database")
                return False
            
            # Process each reaction step with streaming support
            for step_num, (reaction_id, reaction_info) in enumerate(sorted_reactions, 1):
                print(f"\n{'='*80}")
                print(f"🌊 STREAMING REACTION STEP {step_num}/{len(sorted_reactions)}: {reaction_info['name']}")
                print(f"{'='*80}")
                
                step_result = self._process_streaming_reaction_step(
                    step_num, reaction_id, reaction_info, workflow_name, workflow_state,
                    max_workers, chunk_size, parallel_threshold, memory_limit_mb, stream_threshold
                )
                
                workflow_state['step_results'][step_num] = step_result
                
                if not step_result['success']:
                    print(f"❌ Step {step_num} failed. Stopping streaming workflow execution.")
                    break
                
                # Store step products for potential use in next steps
                if step_result['products_generated'] > 0:
                    workflow_state['step_products'][step_num] = {
                        'table_name': step_result.get('output_table'),
                        'product_count': step_result['products_generated'],
                        'reaction_name': reaction_info['name']
                    }
            
            # Reuse existing cleanup and summary methods
            self._query_and_delete_prev_step_tables(workflow_state)
            self._display_streaming_workflow_summary(workflow_state, workflow_name)
            
            return True
            
        except Exception as e:
            print(f"❌ Error applying streaming reaction workflow: {e}")
            return False

    def _process_streaming_reaction_step(self, step_num: int, reaction_id: int, 
                                    reaction_info: Dict[str, Any], workflow_name: str,
                                    workflow_state: Dict[str, Any], max_workers: int,
                                    chunk_size: Optional[int], parallel_threshold: int,
                                    memory_limit_mb: int, stream_threshold: int) -> Dict[str, Any]:
        """
        Process a single reaction step with streaming to disk support.
        Reuses existing helper methods where possible.
        
        Args:
            step_num (int): Current step number
            reaction_id (int): ID of the reaction
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            workflow_state (Dict): Current workflow state
            max_workers (int): Maximum number of parallel workers
            chunk_size (Optional[int]): Size of chunks for parallel processing
            parallel_threshold (int): Minimum compounds to trigger parallel processing
            memory_limit_mb (int): Memory limit in MB
            stream_threshold (int): Threshold to trigger streaming
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            reaction_name = reaction_info['name']
            reaction_smarts = reaction_info['smarts']
            
            # Reuse existing reaction analysis
            reaction_type = self._analyze_single_reaction_type(reaction_smarts)
            
            print(f"🔍 Streaming Reaction Analysis:")
            print(f"   📋 Name: {reaction_name}")
            print(f"   🧪 Type: {reaction_type}")
            print(f"   ⚗️  SMARTS: {reaction_smarts}")
            
            # Reuse existing input source detection
            input_sources = self._get_available_input_sources(step_num, workflow_state)
            
            if not input_sources:
                return {
                    'success': False,
                    'message': 'No input sources available',
                    'products_generated': 0,
                    'streaming_used': False
                }
            
            # Reuse existing table selection logic
            if reaction_type == 'bimolecular':
                table_config = self._select_tables_for_step_bimolecular(
                    input_sources, step_num, reaction_name
                )
            else:  # unimolecular
                table_config = self._select_tables_for_step_unimolecular(
                    input_sources, step_num, reaction_name
                )
            
            if not table_config:
                return {
                    'success': False,
                    'message': 'No tables selected for reaction step',
                    'products_generated': 0,
                    'streaming_used': False
                }
            
            # Apply the reaction with streaming support
            if reaction_type == 'bimolecular':
                step_result = self._apply_streaming_step_bimolecular_reaction(
                    table_config, reaction_info, workflow_name, step_num,
                    max_workers, chunk_size, parallel_threshold, memory_limit_mb, stream_threshold
                )
            else:  # unimolecular
                step_result = self._apply_streaming_step_unimolecular_reaction(
                    table_config, reaction_info, workflow_name, step_num,
                    max_workers, chunk_size, parallel_threshold, memory_limit_mb, stream_threshold
                )
            
            # Update streaming statistics
            if step_result.get('streaming_used', False):
                stats = workflow_state['streaming_stats']
                stats['total_chunks_streamed'] += step_result.get('chunks_streamed', 0)
                stats['total_temp_files'] += step_result.get('temp_files_used', 0)
                stats['total_disk_writes'] += step_result.get('disk_writes', 0)
                stats['peak_memory_usage'] = max(
                    stats['peak_memory_usage'],
                    step_result.get('peak_memory_mb', 0)
                )
            
            return step_result
            
        except Exception as e:
            print(f"❌ Error processing streaming reaction step {step_num}: {e}")
            return {
                'success': False,
                'message': f'Error in streaming step {step_num}: {e}',
                'products_generated': 0,
                'streaming_used': False
            }

    def _apply_streaming_step_bimolecular_reaction(self, table_config: Dict[str, Any], 
                                                reaction_info: Dict[str, Any], 
                                                workflow_name: str, step_num: int,
                                                max_workers: int, chunk_size: Optional[int],
                                                parallel_threshold: int, memory_limit_mb: int,
                                                stream_threshold: int) -> Dict[str, Any]:
        """
        Apply a bimolecular reaction with streaming to disk support.
        
        Args:
            table_config (Dict): Table configuration
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            step_num (int): Current step number
            max_workers (int): Maximum parallel workers
            chunk_size (Optional[int]): Chunk size for parallel processing
            parallel_threshold (int): Threshold to trigger parallel processing
            memory_limit_mb (int): Memory limit in MB
            stream_threshold (int): Threshold to trigger streaming
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            import tempfile
            import os
            import time
            
            print(f"\n🌊 Executing Streaming Bimolecular Reaction - Step {step_num}")
            print("-" * 60)
            
            # Get compound data
            primary_df = self._get_table_as_dataframe(table_config['primary_table'])
            secondary_df = self._get_table_as_dataframe(table_config['secondary_table'])
            
            primary_count = len(primary_df)
            secondary_count = len(secondary_df)
            total_combinations = primary_count * secondary_count
            
            print(f"📊 Primary reactants: {primary_count:,} compounds")
            print(f"📊 Secondary reactants: {secondary_count:,} compounds")
            print(f"🔢 Total combinations: {total_combinations:,}")
            
            # Estimate memory usage and determine if streaming is needed
            estimated_products = total_combinations * 0.1  # Conservative estimate
            use_streaming = estimated_products > stream_threshold
            
            if use_streaming:
                
                # Create temporary directory for streaming
                temp_dir = tempfile.mkdtemp(prefix=f'chemspace_stream_step{step_num}_')
                temp_files_created = []
                
                try:
                    start_time = time.time()
                    total_products = self._stream_bimolecular_reaction_to_disk(
                        primary_df, secondary_df, reaction_info, workflow_name,
                        temp_dir, max_workers, chunk_size, temp_files_created
                    )
                    processing_time = time.time() - start_time
                    
                    # Consolidate temporary files into final table
                    output_table_name = None
                    if total_products > 0:
                        print(f"\n💾 Consolidating {total_products:,} products from {len(temp_files_created)} files...")
                        
                        save_choice = input("Save consolidated products to a new table? (y/n): ").strip().lower()
                        if save_choice in ['y', 'yes']:
                            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                            default_name = f"stream_step{step_num}_{reaction_info['name']}_{timestamp}"
                            user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                            chosen_name = user_table_name if user_table_name else default_name
                            output_table_name = self._sanitize_table_name(chosen_name)
                            
                            success = self._consolidate_temp_files_to_table(
                                temp_files_created, output_table_name, workflow_name, step_num
                            )
                            
                            if not success:
                                output_table_name = None
                    
                    result = {
                        'success': True,
                        'products_generated': total_products,
                        'reaction_type': 'bimolecular',
                        'step_num': step_num,
                        'streaming_used': True,
                        'temp_files_used': len(temp_files_created),
                        'chunks_streamed': len(temp_files_created),
                        'disk_writes': len(temp_files_created),
                        'processing_time': processing_time,
                        'output_table': output_table_name,
                        'combinations_processed': total_combinations
                    }
                    
                finally:
                    # Clean up temporary files
                    self._cleanup_temp_files(temp_files_created, temp_dir)
            
            else:
                print(f"📊 Using in-memory processing (below streaming threshold)")
                
                # Use existing parallel/sequential method
                use_parallel = total_combinations >= parallel_threshold
                
                if use_parallel:
                    # Set default chunk size if not provided
                    if chunk_size is None:
                        chunk_size = max(100, total_combinations // (max_workers * 4))
                    
                    products = self._apply_bimolecular_reaction_parallel(
                        primary_df, secondary_df, reaction_info, workflow_name,
                        max_workers, chunk_size
                    )
                else:
                    # Use existing sequential method
                    products = self._apply_bimolecular_reaction(
                        primary_df, secondary_df, reaction_info, workflow_name
                    )
                
                # Save products using existing logic
                output_table_name = None
                if products:
                    print(f"\n💾 Generated {len(products)} products from bimolecular reaction")
                    
                    save_choice = input("Save products to a new table? (y/n): ").strip().lower()
                    if save_choice in ['y', 'yes']:
                        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                        default_name = f"step{step_num}_{reaction_info['name']}_{timestamp}"
                        user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                        chosen_name = user_table_name if user_table_name else default_name
                        output_table_name = self._sanitize_table_name(chosen_name)
                        
                        success = self._save_step_products(
                            products, output_table_name, workflow_name, step_num
                        )
                        
                        if not success:
                            output_table_name = None
                
                result = {
                    'success': True,
                    'products_generated': len(products),
                    'reaction_type': 'bimolecular',
                    'step_num': step_num,
                    'streaming_used': False,
                    'combinations_processed': total_combinations,
                    'output_table': output_table_name
                }
            
            return result
            
        except Exception as e:
            print(f"❌ Error applying streaming bimolecular reaction: {e}")
            return {
                'success': False,
                'products_generated': 0,
                'error': str(e),
                'streaming_used': False
            }

    def _apply_streaming_step_unimolecular_reaction(self, table_config: Dict[str, Any], 
                                                reaction_info: Dict[str, Any], 
                                                workflow_name: str, step_num: int,
                                                max_workers: int, chunk_size: Optional[int],
                                                parallel_threshold: int, memory_limit_mb: int,
                                                stream_threshold: int) -> Dict[str, Any]:
        """
        Apply a unimolecular reaction with streaming to disk support.
        
        Args:
            table_config (Dict): Table configuration
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            step_num (int): Current step number
            max_workers (int): Maximum parallel workers
            chunk_size (Optional[int]): Chunk size for parallel processing
            parallel_threshold (int): Threshold to trigger parallel processing
            memory_limit_mb (int): Memory limit in MB
            stream_threshold (int): Threshold to trigger streaming
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            import tempfile
            import os
            import time
            
            print(f"\n🌊 Executing Streaming Unimolecular Reaction - Step {step_num}")
            print("-" * 60)
            
            total_compounds = table_config.get('total_compounds', 0)
            
            # Estimate if streaming is needed
            estimated_products = total_compounds * 2  # Conservative estimate for unimolecular
            use_streaming = estimated_products > stream_threshold
            
            if use_streaming:
                print(f"🌊 Using streaming mode (estimated products: {estimated_products:,.0f})")
                
                # Create temporary directory for streaming
                temp_dir = tempfile.mkdtemp(prefix=f'chemspace_stream_step{step_num}_')
                temp_files_created = []
                
                try:
                    start_time = time.time()
                    total_products = 0
                    
                    # Process each source with streaming
                    for source in table_config['sources']:
                        print(f"\n🔬 Streaming processing source: '{source['name']}'")
                        compounds_df = self._get_table_as_dataframe(source['name'])
                        
                        if compounds_df.empty:
                            print(f"   ⚠️  No compounds found, skipping...")
                            continue
                        
                        source_prefix = f"stream_step{step_num}_{source['name']}_"
                        source_products = self._stream_unimolecular_reaction_to_disk(
                            compounds_df, reaction_info, workflow_name, source_prefix,
                            temp_dir, max_workers, chunk_size, temp_files_created
                        )
                        
                        total_products += source_products
                        print(f"   ✅ Streamed {source_products} products to disk")
                    
                    processing_time = time.time() - start_time
                    
                    # Consolidate temporary files into final table
                    output_table_name = None
                    if total_products > 0:
                        print(f"\n💾 Consolidating {total_products:,} products from {len(temp_files_created)} files...")
                        
                        save_choice = input("Save consolidated products to a new table? (y/n): ").strip().lower()
                        if save_choice in ['y', 'yes']:
                            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                            default_name = f"stream_step{step_num}_{reaction_info['name']}_{timestamp}"
                            user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                            chosen_name = user_table_name if user_table_name else default_name
                            output_table_name = self._sanitize_table_name(chosen_name)
                            
                            success = self._consolidate_temp_files_to_table(
                                temp_files_created, output_table_name, workflow_name, step_num
                            )
                            
                            if not success:
                                output_table_name = None
                    
                    result = {
                        'success': True,
                        'products_generated': total_products,
                        'output_table': output_table_name,
                        'reaction_type': 'unimolecular',
                        'step_num': step_num,
                        'streaming_used': True,
                        'temp_files_used': len(temp_files_created),
                        'chunks_streamed': len(temp_files_created),
                        'disk_writes': len(temp_files_created),
                        'processing_time': processing_time,
                        'compounds_processed': total_compounds
                    }
                    
                finally:
                    # Clean up temporary files
                    self._cleanup_temp_files(temp_files_created, temp_dir)
            
            else:
                print(f"📊 Using in-memory processing (below streaming threshold)")
                
                # Use existing sequential/parallel method
                all_products = []
                
                for source in table_config['sources']:
                    print(f"\n🔬 Processing source: '{source['name']}'")
                    compounds_df = self._get_table_as_dataframe(source['name'])
                    
                    if compounds_df.empty:
                        print(f"   ⚠️  No compounds found, skipping...")
                        continue
                    
                    compound_count = len(compounds_df)
                    use_parallel = compound_count >= parallel_threshold
                    
                    if use_parallel:
                        # Set default chunk size if not provided
                        if chunk_size is None:
                            source_chunk_size = max(100, compound_count // (max_workers * 4))
                        else:
                            source_chunk_size = chunk_size
                        
                        source_prefix = f"step{step_num}_{source['name']}_"
                        products = self._apply_unimolecular_reaction_parallel(
                            compounds_df, reaction_info, workflow_name, source_prefix,
                            max_workers, source_chunk_size
                        )
                    else:
                        # Use existing sequential method
                        source_prefix = f"step{step_num}_{source['name']}_"
                        products = self._apply_unimolecular_reaction(
                            compounds_df, reaction_info, workflow_name, source_prefix
                        )
                    
                    all_products.extend(products)
                    print(f"   ✅ Generated {len(products)} products")
                
                # Save products using existing logic
                output_table_name = None
                if all_products:
                    print(f"\n💾 Total products generated: {len(all_products)}")
                    
                    save_choice = input("Save products to a new table? (y/n): ").strip().lower()
                    if save_choice in ['y', 'yes']:
                        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                        default_name = f"step{step_num}_{reaction_info['name']}_{timestamp}"
                        user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                        chosen_name = user_table_name if user_table_name else default_name
                        output_table_name = self._sanitize_table_name(chosen_name)
                        
                        success = self._save_step_products(
                            all_products, output_table_name, workflow_name, step_num
                        )
                        
                        if not success:
                            output_table_name = None
                
                result = {
                    'success': True,
                    'products_generated': len(all_products),
                    'output_table': output_table_name,
                    'reaction_type': 'unimolecular',
                    'step_num': step_num,
                    'streaming_used': False,
                    'compounds_processed': total_compounds
                }
            
            return result
            
        except Exception as e:
            print(f"❌ Error applying streaming unimolecular reaction: {e}")
            return {
                'success': False,
                'products_generated': 0,
                'error': str(e),
                'streaming_used': False
            }

    def _stream_bimolecular_reaction_to_disk(self, primary_df: pd.DataFrame, secondary_df: pd.DataFrame,
                                        reaction_info: Dict[str, Any], workflow_name: str,
                                        temp_dir: str, max_workers: int, chunk_size: Optional[int],
                                        temp_files_created: List[str]) -> int:
        """
        Stream bimolecular reaction results directly to disk files.
        
        Args:
            primary_df (pd.DataFrame): Primary reactants
            secondary_df (pd.DataFrame): Secondary reactants
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            temp_dir (str): Temporary directory for files
            max_workers (int): Maximum workers
            chunk_size (Optional[int]): Chunk size
            temp_files_created (List[str]): List to track created files
            
        Returns:
            int: Total number of products streamed
        """
        try:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            import csv
            import os
            
            reaction_smarts = reaction_info['smarts']
            reaction_name = reaction_info['name']
            
            # Set default chunk size if not provided
            total_combinations = len(primary_df) * len(secondary_df)
            if chunk_size is None:
                chunk_size = max(1000, total_combinations // (max_workers * 8))  # Smaller chunks for streaming
            
            # Create chunks for streaming processing
            chunks = self._create_bimolecular_chunks(primary_df, secondary_df, chunk_size)
            
            print(f"   📦 Created {len(chunks)} chunks for streaming processing")
            print(f"   🌊 Each chunk will be streamed directly to disk")
            
            total_products = 0
            processed_chunks = 0
            
            # Initialize progress tracking
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(chunks),
                    desc=f"Streaming {reaction_name}",
                    unit="chunks"
                )
            else:
                progress_bar = None
            
            # Process chunks and stream to disk
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit streaming jobs
                future_to_chunk = {}
                for i, chunk in enumerate(chunks):
                    temp_file_path = os.path.join(temp_dir, f'bimolecular_chunk_{i:06d}.csv')
                    temp_files_created.append(temp_file_path)
                    
                    future = executor.submit(
                        _process_bimolecular_chunk_to_file_worker,
                        chunk, reaction_smarts, reaction_name, workflow_name, temp_file_path
                    )
                    future_to_chunk[future] = i
                
                # Collect results
                for future in as_completed(future_to_chunk):
                    chunk_idx = future_to_chunk[future]
                    
                    try:
                        chunk_products_count = future.result()
                        total_products += chunk_products_count
                        processed_chunks += 1
                        
                        if progress_bar:
                            progress_bar.update(1)
                            progress_bar.set_postfix({
                                'products': total_products,
                                'chunks': f"{processed_chunks}/{len(chunks)}"
                            })
                        else:
                            if processed_chunks % max(1, len(chunks) // 10) == 0:
                                progress = (processed_chunks / len(chunks)) * 100
                                print(f"      📊 Progress: {progress:.1f}% ({total_products:,} products streamed)")
                    
                    except Exception as e:
                        print(f"   ❌ Error processing chunk {chunk_idx}: {e}")
                        processed_chunks += 1
                        if progress_bar:
                            progress_bar.update(1)
            
            if progress_bar:
                progress_bar.close()
            
            print(f"   ✅ Streaming completed: {total_products:,} products in {len(temp_files_created)} files")
            return total_products
            
        except Exception as e:
            print(f"   ❌ Error in streaming bimolecular reaction: {e}")
            return 0

    def _stream_unimolecular_reaction_to_disk(self, compounds_df: pd.DataFrame, 
                                            reaction_info: Dict[str, Any], workflow_name: str,
                                            name_prefix: str, temp_dir: str, max_workers: int,
                                            chunk_size: Optional[int], temp_files_created: List[str]) -> int:
        """
        Stream unimolecular reaction results directly to disk files.
        
        Args:
            compounds_df (pd.DataFrame): Compounds dataframe
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            name_prefix (str): Prefix for product names
            temp_dir (str): Temporary directory for files
            max_workers (int): Maximum workers
            chunk_size (Optional[int]): Chunk size
            temp_files_created (List[str]): List to track created files
            
        Returns:
            int: Total number of products streamed
        """
        try:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            import os
            
            reaction_smarts = reaction_info['smarts']
            reaction_name = reaction_info['name']
            
            # Set default chunk size if not provided
            if chunk_size is None:
                chunk_size = max(500, len(compounds_df) // (max_workers * 8))  # Smaller chunks for streaming
            
            # Create chunks for streaming processing
            chunks = self._create_unimolecular_chunks(compounds_df, chunk_size)
            
            print(f"      📦 Created {len(chunks)} chunks for streaming processing")
            
            total_products = 0
            processed_chunks = 0
            
            # Initialize progress tracking
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(chunks),
                    desc=f"Streaming {reaction_name}",
                    unit="chunks"
                )
            else:
                progress_bar = None
            
            # Process chunks and stream to disk
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit streaming jobs
                future_to_chunk = {}
                for i, chunk in enumerate(chunks):
                    temp_file_path = os.path.join(temp_dir, f'unimolecular_chunk_{i:06d}.csv')
                    temp_files_created.append(temp_file_path)
                    
                    future = executor.submit(
                        _process_unimolecular_chunk_to_file_worker,
                        chunk, reaction_smarts, reaction_name, workflow_name, name_prefix, temp_file_path
                    )
                    future_to_chunk[future] = i
                
                # Collect results
                for future in as_completed(future_to_chunk):
                    chunk_idx = future_to_chunk[future]
                    
                    try:
                        chunk_products_count = future.result()
                        total_products += chunk_products_count
                        processed_chunks += 1
                        
                        if progress_bar:
                            progress_bar.update(1)
                            progress_bar.set_postfix({
                                'products': total_products
                            })
                    
                    except Exception as e:
                        print(f"      ❌ Error processing chunk {chunk_idx}: {e}")
                        processed_chunks += 1
                        if progress_bar:
                            progress_bar.update(1)
            
            if progress_bar:
                progress_bar.close()
            
            print(f"      ✅ Streaming completed: {total_products:,} products")
            return total_products
            
        except Exception as e:
            print(f"      ❌ Error in streaming unimolecular reaction: {e}")
            return 0

    def _consolidate_temp_files_to_table(self, temp_files: List[str], output_table_name: str,
                                        workflow_name: str, step_num: int) -> bool:
        """
        Consolidate temporary CSV files into a database table with streaming to avoid memory issues.
        
        Args:
            temp_files (List[str]): List of temporary file paths
            output_table_name (str): Name for the output table
            workflow_name (str): Workflow name
            step_num (int): Step number
            
        Returns:
            bool: True if consolidation was successful
        """
        try:
            import csv
            
            if not temp_files:
                return False
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create products table
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {output_table_name} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT NOT NULL,
                    name TEXT,
                    flag TEXT,
                    workflow_step INTEGER,
                    workflow_name TEXT,
                    creation_date TEXT,
                    UNIQUE(smiles, name)
                )
            ''')
            
            # Create indexes
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_smiles ON {output_table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_step ON {output_table_name}(workflow_step)")
            
            # Insert products from temp files with streaming
            insert_query = f'''
                INSERT OR IGNORE INTO {output_table_name} 
                (smiles, name, flag, workflow_step, workflow_name, creation_date)
                VALUES (?, ?, ?, ?, ?, ?)
            '''
            
            creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            total_inserted = 0
            
            # Process files in batches to avoid memory issues
            batch_size = 1000
            
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    temp_files,
                    desc="Consolidating files",
                    unit="files"
                )
            else:
                progress_bar = temp_files
            
            for temp_file in progress_bar:
                try:
                    with open(temp_file, 'r', newline='', encoding='utf-8') as f:
                        reader = csv.DictReader(f)
                        batch = []
                        
                        for row in reader:
                            batch.append((
                                row.get('smiles', ''),
                                row.get('name', ''),
                                row.get('flag', 'stream_product'),
                                step_num,
                                workflow_name,
                                creation_date
                            ))
                            
                            if len(batch) >= batch_size:
                                cursor.executemany(insert_query, batch)
                                total_inserted += len(batch)
                                batch = []
                        
                        # Insert remaining items in batch
                        if batch:
                            cursor.executemany(insert_query, batch)
                            total_inserted += len(batch)
                            
                except Exception as e:
                    print(f"   ⚠️  Error processing file {temp_file}: {e}")
                    continue
            
            conn.commit()
            conn.close()
            
            print(f"   💾 Consolidated {total_inserted:,} products to table '{output_table_name}'")
            return True
            
        except Exception as e:
            print(f"   ❌ Error consolidating temp files: {e}")
            return False

    def _cleanup_temp_files(self, temp_files: List[str], temp_dir: str) -> None:
        """
        Clean up temporary files and directory.
        
        Args:
            temp_files (List[str]): List of temporary file paths
            temp_dir (str): Temporary directory path
        """
        try:
            import os
            import shutil
            
            # Remove individual files
            for temp_file in temp_files:
                try:
                    if os.path.exists(temp_file):
                        os.unlink(temp_file)
                except Exception:
                    pass  # Continue cleanup even if individual file fails
            
            # Remove temporary directory
            try:
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
            except Exception:
                pass  # Continue even if directory removal fails
            
            print(f"   🧹 Cleaned up {len(temp_files)} temporary files")
            
        except Exception as e:
            print(f"   ⚠️  Warning: Error during cleanup: {e}")

    def _display_streaming_workflow_summary(self, workflow_state: Dict[str, Any], workflow_name: str) -> None:
        """
        Display a comprehensive summary of the streaming workflow execution.
        
        Args:
            workflow_state (Dict): Final workflow state
            workflow_name (str): Name of the workflow
        """
        try:
            # Reuse existing summary display
            self._display_reaction_workflow_summary(workflow_state, workflow_name)
            
            # Add streaming statistics
            streaming_stats = workflow_state.get('streaming_stats', {})
            
            if streaming_stats.get('total_chunks_streamed', 0) > 0:
                print(f"\n🌊 STREAMING PROCESSING STATISTICS:")
                print(f"   📦 Total chunks streamed: {streaming_stats.get('total_chunks_streamed', 0)}")
                print(f"   📄 Total temp files created: {streaming_stats.get('total_temp_files', 0)}")
                print(f"   💾 Total disk writes: {streaming_stats.get('total_disk_writes', 0)}")
                print(f"   📊 Peak memory usage: {streaming_stats.get('peak_memory_usage', 0):.1f}MB")
                
                print(f"{'='*80}")
            
        except Exception as e:
            print(f"❌ Error displaying streaming workflow summary: {e}")

    def _process_parallel_reaction_step(self, step_num: int, reaction_id: int, 
                                    reaction_info: Dict[str, Any], workflow_name: str,
                                    workflow_state: Dict[str, Any], max_workers: int,
                                    chunk_size: Optional[int], parallel_threshold: int) -> Dict[str, Any]:
        """
        Process a single reaction step with parallel processing support.
        Reuses existing helper methods where possible.
        
        Args:
            step_num (int): Current step number
            reaction_id (int): ID of the reaction
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            workflow_state (Dict): Current workflow state
            max_workers (int): Maximum number of parallel workers
            chunk_size (Optional[int]): Size of chunks for parallel processing
            parallel_threshold (int): Minimum compounds to trigger parallel processing
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            reaction_name = reaction_info['name']
            reaction_smarts = reaction_info['smarts']
            
            # Reuse existing reaction analysis
            reaction_type = self._analyze_single_reaction_type(reaction_smarts)
            
            print(f"🔍 Parallel Reaction Analysis:")
            print(f"   📋 Name: {reaction_name}")
            print(f"   🧪 Type: {reaction_type}")
            print(f"   ⚗️  SMARTS: {reaction_smarts}")
            
            # Reuse existing input source detection
            input_sources = self._get_available_input_sources(step_num, workflow_state)
            
            if not input_sources:
                return {
                    'success': False,
                    'message': 'No input sources available',
                    'products_generated': 0,
                    'parallel_used': False
                }
            
            # Reuse existing table selection logic
            if reaction_type == 'bimolecular':
                table_config = self._select_tables_for_step_bimolecular(
                    input_sources, step_num, reaction_name
                )
            else:  # unimolecular
                table_config = self._select_tables_for_step_unimolecular(
                    input_sources, step_num, reaction_name
                )
            
            if not table_config:
                return {
                    'success': False,
                    'message': 'No tables selected for reaction step',
                    'products_generated': 0,
                    'parallel_used': False
                }
            
            # Apply the reaction with parallel support
            if reaction_type == 'bimolecular':
                step_result = self._apply_parallel_step_bimolecular_reaction(
                    table_config, reaction_info, workflow_name, step_num,
                    max_workers, chunk_size, parallel_threshold
                )
            else:  # unimolecular
                step_result = self._apply_parallel_step_unimolecular_reaction(
                    table_config, reaction_info, workflow_name, step_num,
                    max_workers, chunk_size, parallel_threshold
                )
            
            # Update parallel statistics
            if step_result.get('parallel_used', False):
                stats = workflow_state['parallel_stats']
                stats['total_workers_used'] += step_result.get('workers_used', 0)
                stats['total_chunks_processed'] += step_result.get('chunks_processed', 0)
                stats['total_processing_time'] += step_result.get('processing_time', 0.0)
            
            return step_result
            
        except Exception as e:
            print(f"❌ Error processing parallel reaction step {step_num}: {e}")
            return {
                'success': False,
                'message': f'Error in parallel step {step_num}: {e}',
                'products_generated': 0,
                'parallel_used': False
            }

    def _apply_parallel_step_bimolecular_reaction(self, table_config: Dict[str, Any], 
                                                reaction_info: Dict[str, Any], 
                                                workflow_name: str, step_num: int,
                                                max_workers: int, chunk_size: Optional[int],
                                                parallel_threshold: int) -> Dict[str, Any]:
        """
        Apply a bimolecular reaction with parallel processing support.
        
        Args:
            table_config (Dict): Table configuration
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            step_num (int): Current step number
            max_workers (int): Maximum parallel workers
            chunk_size (Optional[int]): Chunk size for parallel processing
            parallel_threshold (int): Threshold to trigger parallel processing
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            import time
            
            print(f"\n🚀 Executing Parallel Bimolecular Reaction - Step {step_num}")
            print("-" * 60)
            
            # Get compound data
            primary_df = self._get_table_as_dataframe(table_config['primary_table'])
            secondary_df = self._get_table_as_dataframe(table_config['secondary_table'])
            
            primary_count = len(primary_df)
            secondary_count = len(secondary_df)
            total_combinations = primary_count * secondary_count
            
            print(f"📊 Primary reactants: {primary_count:,} compounds")
            print(f"📊 Secondary reactants: {secondary_count:,} compounds")
            print(f"🔢 Total combinations: {total_combinations:,}")
            
            # Determine if parallel processing should be used
            use_parallel = total_combinations >= parallel_threshold
            
            if use_parallel:
                print(f"🚀 Using parallel processing (threshold: {parallel_threshold:,})")
                
                # Set default chunk size if not provided
                if chunk_size is None:
                    chunk_size = max(100, total_combinations // (max_workers * 4))
                
                print(f"   👥 Workers: {max_workers}")
                print(f"   📦 Chunk size: {chunk_size:,}")
                
                # Apply parallel bimolecular reaction
                start_time = time.time()
                products = self._apply_bimolecular_reaction_parallel(
                    primary_df, secondary_df, reaction_info, workflow_name,
                    max_workers, chunk_size
                )
                processing_time = time.time() - start_time
                
                chunks_processed = (total_combinations + chunk_size - 1) // chunk_size
                
                result = {
                    'success': True,
                    'products_generated': len(products),
                    'reaction_type': 'bimolecular',
                    'step_num': step_num,
                    'parallel_used': True,
                    'workers_used': max_workers,
                    'chunks_processed': chunks_processed,
                    'processing_time': processing_time,
                    'combinations_processed': total_combinations
                }
            else:
                print(f"🔄 Using sequential processing (below threshold)")
                
                # Use existing sequential method
                products = self._apply_bimolecular_reaction(
                    primary_df, secondary_df, reaction_info, workflow_name
                )
                
                result = {
                    'success': True,
                    'products_generated': len(products),
                    'reaction_type': 'bimolecular',
                    'step_num': step_num,
                    'parallel_used': False,
                    'combinations_processed': total_combinations
                }
            
            # Save products using existing logic
            output_table_name = None
            if products:
                print(f"\n💾 Generated {len(products)} products from bimolecular reaction")
                
                save_choice = input("Save products to a new table? (y/n): ").strip().lower()
                if save_choice in ['y', 'yes']:
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    default_name = f"parallel_step{step_num}_{reaction_info['name']}_{timestamp}"
                    user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                    chosen_name = user_table_name if user_table_name else default_name
                    output_table_name = self._sanitize_table_name(chosen_name)
                    
                    # Reuse existing save method
                    success = self._save_step_products(
                        products, output_table_name, workflow_name, step_num
                    )
                    
                    if not success:
                        output_table_name = None
            
            result['output_table'] = output_table_name
            return result
            
        except Exception as e:
            print(f"❌ Error applying parallel bimolecular reaction: {e}")
            return {
                'success': False,
                'products_generated': 0,
                'error': str(e),
                'parallel_used': False
            }

    def _apply_parallel_step_unimolecular_reaction(self, table_config: Dict[str, Any], 
                                                reaction_info: Dict[str, Any], 
                                                workflow_name: str, step_num: int,
                                                max_workers: int, chunk_size: Optional[int],
                                                parallel_threshold: int) -> Dict[str, Any]:
        """
        Apply a unimolecular reaction with parallel processing support.
        
        Args:
            table_config (Dict): Table configuration
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            step_num (int): Current step number
            max_workers (int): Maximum parallel workers
            chunk_size (Optional[int]): Chunk size for parallel processing
            parallel_threshold (int): Threshold to trigger parallel processing
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            import time
            
            print(f"\n🚀 Executing Parallel Unimolecular Reaction - Step {step_num}")
            print("-" * 60)
            
            all_products = []
            total_compounds = 0
            parallel_used = False
            total_processing_time = 0.0
            total_chunks_processed = 0
            
            # Process each selected source
            for source in table_config['sources']:
                print(f"\n🔬 Processing source: '{source['name']}'")
                compounds_df = self._get_table_as_dataframe(source['name'])
                
                if compounds_df.empty:
                    print(f"   ⚠️  No compounds found, skipping...")
                    continue
                
                compound_count = len(compounds_df)
                total_compounds += compound_count
                print(f"   📊 Processing {compound_count:,} compounds...")
                
                # Determine if parallel processing should be used for this source
                use_parallel_for_source = compound_count >= parallel_threshold
                
                if use_parallel_for_source:
                    print(f"   🚀 Using parallel processing for this source")
                    parallel_used = True
                    
                    # Set default chunk size if not provided
                    if chunk_size is None:
                        source_chunk_size = max(100, compound_count // (max_workers * 4))
                    else:
                        source_chunk_size = chunk_size
                    
                    print(f"      👥 Workers: {max_workers}")
                    print(f"      📦 Chunk size: {source_chunk_size:,}")
                    
                    # Apply parallel unimolecular reaction
                    start_time = time.time()
                    source_prefix = f"parallel_step{step_num}_{source['name']}_"
                    products = self._apply_unimolecular_reaction_parallel(
                        compounds_df, reaction_info, workflow_name, source_prefix,
                        max_workers, source_chunk_size
                    )
                    processing_time = time.time() - start_time
                    total_processing_time += processing_time
                    
                    chunks_for_source = (compound_count + source_chunk_size - 1) // source_chunk_size
                    total_chunks_processed += chunks_for_source
                    
                    print(f"   ⚡ Parallel processing completed in {processing_time:.2f}s")
                    print(f"   📦 Chunks processed: {chunks_for_source}")
                else:
                    print(f"   🔄 Using sequential processing (below threshold)")
                    
                    # Use existing sequential method
                    source_prefix = f"step{step_num}_{source['name']}_"
                    products = self._apply_unimolecular_reaction(
                        compounds_df, reaction_info, workflow_name, source_prefix
                    )
                
                all_products.extend(products)
                print(f"   ✅ Generated {len(products)} products")
            
            # Save products using existing logic
            output_table_name = None
            if all_products:
                print(f"\n💾 Total products generated: {len(all_products)}")
                
                save_choice = input("Save products to a new table? (y/n): ").strip().lower()
                if save_choice in ['y', 'yes']:
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    default_name = f"parallel_step{step_num}_{reaction_info['name']}_{timestamp}"
                    user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                    chosen_name = user_table_name if user_table_name else default_name
                    output_table_name = self._sanitize_table_name(chosen_name)
                    
                    # Reuse existing save method
                    success = self._save_step_products(
                        all_products, output_table_name, workflow_name, step_num
                    )
                    
                    if not success:
                        output_table_name = None
            
            return {
                'success': True,
                'products_generated': len(all_products),
                'output_table': output_table_name,
                'reaction_type': 'unimolecular',
                'step_num': step_num,
                'parallel_used': parallel_used,
                'workers_used': max_workers if parallel_used else 0,
                'chunks_processed': total_chunks_processed,
                'processing_time': total_processing_time,
                'compounds_processed': total_compounds
            }
            
        except Exception as e:
            print(f"❌ Error applying parallel unimolecular reaction: {e}")
            return {
                'success': False,
                'products_generated': 0,
                'error': str(e),
                'parallel_used': False
            }

    def _apply_bimolecular_reaction_parallel(self, primary_df: pd.DataFrame, secondary_df: pd.DataFrame,
                                        reaction_info: Dict[str, Any], workflow_name: str,
                                        max_workers: int, chunk_size: int) -> List[Dict[str, Any]]:
        """
        Apply a bimolecular reaction using parallel processing.
        
        Args:
            primary_df (pd.DataFrame): Primary reactants dataframe
            secondary_df (pd.DataFrame): Secondary reactants dataframe
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            max_workers (int): Maximum number of parallel workers
            chunk_size (int): Size of each chunk for parallel processing
            
        Returns:
            List[Dict]: List of reaction products
        """
        try:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            
            reaction_smarts = reaction_info['smarts']
            reaction_name = reaction_info['name']
            
            print(f"   🔗 Reaction: {reaction_name}")
            print(f"   ⚗️  SMARTS: {reaction_smarts}")
            
            # Prepare data chunks for parallel processing
            total_combinations = len(primary_df) * len(secondary_df)
            chunks = self._create_bimolecular_chunks(primary_df, secondary_df, chunk_size)
            
            print(f"   📦 Created {len(chunks)} chunks for parallel processing")
            
            all_products = []
            processed_chunks = 0
            successful_reactions = 0
            
            # Initialize progress tracking
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(chunks),
                    desc=f"Parallel {reaction_name}",
                    unit="chunks"
                )
            else:
                progress_bar = None
            
            # Process chunks in parallel
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all chunks
                future_to_chunk = {
                    executor.submit(
                        _process_bimolecular_chunk_worker,
                        chunk, reaction_smarts, reaction_name, workflow_name
                    ): i for i, chunk in enumerate(chunks)
                }
                
                # Collect results
                for future in as_completed(future_to_chunk):
                    chunk_idx = future_to_chunk[future]
                    
                    try:
                        chunk_products = future.result()
                        all_products.extend(chunk_products)
                        successful_reactions += len(chunk_products)
                        processed_chunks += 1
                        
                        if progress_bar:
                            progress_bar.update(1)
                            progress_bar.set_postfix({
                                'products': len(all_products),
                                'chunks': f"{processed_chunks}/{len(chunks)}"
                            })
                        else:
                            if processed_chunks % max(1, len(chunks) // 10) == 0:
                                progress = (processed_chunks / len(chunks)) * 100
                                print(f"      📊 Progress: {progress:.1f}% ({len(all_products)} products)")
                    
                    except Exception as e:
                        print(f"   ❌ Error processing chunk {chunk_idx}: {e}")
                        processed_chunks += 1
                        if progress_bar:
                            progress_bar.update(1)
            
            if progress_bar:
                progress_bar.close()
            
            print(f"   ✅ Parallel processing completed: {len(all_products)} products from {successful_reactions} reactions")
            return all_products
            
        except Exception as e:
            print(f"   ❌ Error in parallel bimolecular reaction: {e}")
            return []

    def _apply_unimolecular_reaction_parallel(self, compounds_df: pd.DataFrame, 
                                            reaction_info: Dict[str, Any], workflow_name: str,
                                            name_prefix: str, max_workers: int, 
                                            chunk_size: int) -> List[Dict[str, Any]]:
        """
        Apply a unimolecular reaction using parallel processing.
        
        Args:
            compounds_df (pd.DataFrame): Compounds dataframe
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            name_prefix (str): Prefix for product names
            max_workers (int): Maximum number of parallel workers
            chunk_size (int): Size of each chunk for parallel processing
            
        Returns:
            List[Dict]: List of reaction products
        """
        try:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            
            reaction_smarts = reaction_info['smarts']
            reaction_name = reaction_info['name']
            
            # Create data chunks for parallel processing
            chunks = self._create_unimolecular_chunks(compounds_df, chunk_size)
            
            print(f"      📦 Created {len(chunks)} chunks for parallel processing")
            
            all_products = []
            processed_chunks = 0
            successful_reactions = 0
            
            # Initialize progress tracking
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    total=len(chunks),
                    desc=f"Parallel {reaction_name}",
                    unit="chunks"
                )
            else:
                progress_bar = None
            
            # Process chunks in parallel
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all chunks
                future_to_chunk = {
                    executor.submit(
                        _process_unimolecular_chunk_worker,
                        chunk, reaction_smarts, reaction_name, workflow_name, name_prefix
                    ): i for i, chunk in enumerate(chunks)
                }
                
                # Collect results
                for future in as_completed(future_to_chunk):
                    chunk_idx = future_to_chunk[future]
                    
                    try:
                        chunk_products = future.result()
                        all_products.extend(chunk_products)
                        successful_reactions += len(chunk_products)
                        processed_chunks += 1
                        
                        if progress_bar:
                            progress_bar.update(1)
                            progress_bar.set_postfix({
                                'products': len(all_products)
                            })
                        
                    except Exception as e:
                        print(f"      ❌ Error processing chunk {chunk_idx}: {e}")
                        processed_chunks += 1
                        if progress_bar:
                            progress_bar.update(1)
            
            if progress_bar:
                progress_bar.close()
            
            print(f"      ✅ Parallel processing completed: {len(all_products)} products")
            return all_products
            
        except Exception as e:
            print(f"      ❌ Error in parallel unimolecular reaction: {e}")
            return []

    def _create_bimolecular_chunks(self, primary_df: pd.DataFrame, secondary_df: pd.DataFrame, 
                                chunk_size: int) -> List[List[Tuple]]:
        """
        Create chunks for parallel bimolecular reaction processing.
        
        Args:
            primary_df (pd.DataFrame): Primary reactants
            secondary_df (pd.DataFrame): Secondary reactants
            chunk_size (int): Size of each chunk
            
        Returns:
            List[List[Tuple]]: List of chunks, each containing compound pairs
        """
        chunks = []
        current_chunk = []
        
        for _, primary_compound in primary_df.iterrows():
            for _, secondary_compound in secondary_df.iterrows():
                current_chunk.append((
                    primary_compound.to_dict(),
                    secondary_compound.to_dict()
                ))
                
                if len(current_chunk) >= chunk_size:
                    chunks.append(current_chunk)
                    current_chunk = []
        
        # Add remaining compounds
        if current_chunk:
            chunks.append(current_chunk)
        
        return chunks

    def _create_unimolecular_chunks(self, compounds_df: pd.DataFrame, 
                                chunk_size: int) -> List[List[Dict]]:
        """
        Create chunks for parallel unimolecular reaction processing.
        
        Args:
            compounds_df (pd.DataFrame): Compounds dataframe
            chunk_size (int): Size of each chunk
            
        Returns:
            List[List[Dict]]: List of chunks, each containing compound dictionaries
        """
        chunks = []
        current_chunk = []
        
        for _, compound in compounds_df.iterrows():
            current_chunk.append(compound.to_dict())
            
            if len(current_chunk) >= chunk_size:
                chunks.append(current_chunk)
                current_chunk = []
        
        # Add remaining compounds
        if current_chunk:
            chunks.append(current_chunk)
        
        return chunks

    def _display_parallel_workflow_summary(self, workflow_state: Dict[str, Any], workflow_name: str) -> None:
        """
        Display a comprehensive summary of the parallel workflow execution.
        Extends the existing summary with parallel processing statistics.
        
        Args:
            workflow_state (Dict): Final workflow state
            workflow_name (str): Name of the workflow
        """
        try:
            # Reuse existing summary display
            self._display_reaction_workflow_summary(workflow_state, workflow_name)
            
            # Add parallel processing statistics
            parallel_stats = workflow_state.get('parallel_stats', {})
            
            if parallel_stats.get('total_processing_time', 0) > 0:
                print(f"\n🚀 PARALLEL PROCESSING STATISTICS:")
                print(f"   👥 Total workers used: {parallel_stats.get('total_workers_used', 0)}")
                print(f"   📦 Total chunks processed: {parallel_stats.get('total_chunks_processed', 0)}")
                print(f"   ⏱️  Total processing time: {parallel_stats.get('total_processing_time', 0):.2f}s")
                
                # Calculate processing efficiency
                total_time = parallel_stats.get('total_processing_time', 0)
                total_workers = parallel_stats.get('total_workers_used', 1)
                if total_time > 0 and total_workers > 0:
                    efficiency = (total_time / total_workers)
                    print(f"   ⚡ Average time per worker: {efficiency:.2f}s")
                
                print(f"{'='*80}")
            
        except Exception as e:
            print(f"❌ Error displaying parallel workflow summary: {e}")

    def _select_reaction_workflow_for_application(self) -> Optional[int]:
        """
        Interactive selection of reaction workflow for application.
        Reuses existing helper methods for consistency.
        
        Returns:
            Optional[int]: Selected workflow ID or None if cancelled
        """
        try:
            # Check if reaction_workflows table exists
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='reaction_workflows'")
            if not cursor.fetchone():
                print("❌ No reaction workflows table found")
                print("   Create workflows first using create_reaction_workflow()")
                conn.close()
                return None
            
            # Get all available workflows (reuse existing query pattern)
            cursor.execute("""
                SELECT id, workflow_name, creation_date, description, reactions_dict
                FROM reaction_workflows 
                ORDER BY creation_date DESC
            """)
            workflows = cursor.fetchall()
            conn.close()
            
            if not workflows:
                print("❌ No reaction workflows found")
                print("   Create workflows first using create_reaction_workflow()")
                return None
            
            print(f"\n🔬 SELECT REACTION WORKFLOW FOR APPLICATION")
            print("=" * 70)
            print(f"Available workflows ({len(workflows)} total):")
            print("-" * 70)
            print(f"{'#':<3} {'Workflow Name':<25} {'Reactions':<10} {'Created':<12} {'Description':<20}")
            print("-" * 70)
            
            workflow_info = []
            for i, (wf_id, name, date, desc, reactions_dict_str) in enumerate(workflows, 1):
                try:
                    reactions_dict = json.loads(reactions_dict_str)
                    reaction_count = len(reactions_dict)
                except json.JSONDecodeError:
                    reaction_count = 0
                
                date_short = date[:10] if date else "Unknown"
                name_display = name[:24] if len(name) <= 24 else name[:21] + "..."
                desc_display = (desc or "No description")[:19]
                if len(desc_display) > 19:
                    desc_display = desc_display[:16] + "..."
                
                print(f"{i:<3} {name_display:<25} {reaction_count:<10} {date_short:<12} {desc_display:<20}")
                workflow_info.append({
                    'id': wf_id,
                    'name': name,
                    'reaction_count': reaction_count,
                    'description': desc or "No description",
                    'date': date
                })
            
            print("-" * 70)
            print("Commands: Enter workflow number, workflow name, 'list' for details, or 'cancel' to abort")
            
            while True:
                try:
                    selection = input(f"\n🔍 Select workflow for application: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    elif selection.lower() == 'list':
                        # Reuse existing detailed workflow display (adapted for application context)
                        self._show_detailed_workflow_list(workflows)
                        continue
                    
                    # Try as number first
                    try:
                        workflow_idx = int(selection) - 1
                        if 0 <= workflow_idx < len(workflows):
                            selected_workflow_id = workflows[workflow_idx][0]  # workflow ID
                            selected_info = workflow_info[workflow_idx]
                            
                            print(f"\n✅ Selected workflow: '{selected_info['name']}' (ID: {selected_workflow_id})")
                            print(f"   🧪 Reactions: {selected_info['reaction_count']}")
                            print(f"   📄 Description: {selected_info['description']}")
                            print(f"   📅 Created: {selected_info['date']}")
                            return selected_workflow_id
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(workflows)}")
                            continue
                    except ValueError:
                        # Try as workflow name (reuse existing search logic)
                        matching_workflows = [w for w in workflows if w[1].lower() == selection.lower()]
                        if matching_workflows:
                            selected_workflow = matching_workflows[0]
                            selected_workflow_id = selected_workflow[0]
                            # Find workflow info
                            selected_info = next((info for info in workflow_info if info['id'] == selected_workflow_id), 
                                            {'name': 'Unknown', 'reaction_count': 0, 'description': 'Unknown'})
                            
                            print(f"\n✅ Selected workflow: '{selected_info['name']}' (ID: {selected_workflow_id})")
                            print(f"   🧪 Reactions: {selected_info['reaction_count']}")
                            print(f"   📄 Description: {selected_info['description']}")
                            return selected_workflow_id
                        else:
                            print(f"❌ Workflow '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Workflow selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error selecting workflow for application: {e}")
            return None

    def _process_sequential_reaction_step(self, step_num: int, reaction_id: int, 
                                        reaction_info: Dict[str, Any], workflow_name: str,
                                        workflow_state: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a single reaction step in the sequential workflow.
        
        Args:
            step_num (int): Current step number
            reaction_id (int): ID of the reaction
            reaction_info (Dict): Reaction information
            workflow_name (str): Name of the workflow
            workflow_state (Dict): Current workflow state
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            reaction_name = reaction_info['name']
            reaction_smarts = reaction_info['smarts']
            
            # Analyze this specific reaction
            reaction_type = self._analyze_single_reaction_type(reaction_smarts)
            
            print(f"🔍 Reaction Analysis:")
            print(f"   📋 Name: {reaction_name}")
            print(f"   🧪 Type: {reaction_type}")
            print(f"   ⚗️  SMARTS: {reaction_smarts}")
            
            # Get available input sources
            input_sources = self._get_available_input_sources(step_num, workflow_state)
            
            if not input_sources:
                return {
                    'success': False,
                    'message': 'No input sources available',
                    'products_generated': 0
                }
            
            # Select input tables based on reaction type and available sources
            if reaction_type == 'bimolecular':
                table_config = self._select_tables_for_step_bimolecular(
                    input_sources, step_num, reaction_name
                )
            else:  # unimolecular
                table_config = self._select_tables_for_step_unimolecular(
                    input_sources, step_num, reaction_name
                )
            
            if not table_config:
                return {
                    'success': False,
                    'message': 'No tables selected for reaction step',
                    'products_generated': 0
                }
            
            # Apply the reaction
            if reaction_type == 'bimolecular':
                step_result = self._apply_step_bimolecular_reaction(
                    table_config, reaction_info, workflow_name, step_num
                )
            else:  # unimolecular
                step_result = self._apply_step_unimolecular_reaction(
                    table_config, reaction_info, workflow_name, step_num
                )
            
            return step_result
            
        except Exception as e:
            print(f"❌ Error processing reaction step {step_num}: {e}")
            return {
                'success': False,
                'message': f'Error in step {step_num}: {e}',
                'products_generated': 0
            }

    def _analyze_single_reaction_type(self, smarts: str) -> str:
        """
        Analyze a single reaction SMARTS to determine if it's unimolecular or bimolecular.
        
        Args:
            smarts (str): SMARTS pattern of the reaction
            
        Returns:
            str: 'unimolecular', 'bimolecular', or 'invalid'
        """
        try:
            if '>>' not in smarts:
                return 'invalid'
            
            reactants_part = smarts.split('>>')[0]
            dot_count = reactants_part.count('.')
            
            if dot_count >= 1:
                return 'bimolecular'
            else:
                return 'unimolecular'
                
        except Exception:
            return 'invalid'

    def _get_available_input_sources(self, step_num: int, workflow_state: Dict[str, Any]) -> Dict[str, List[str]]:
        """
        Get available input sources for the current reaction step.
        
        Args:
            step_num (int): Current step number
            workflow_state (Dict): Current workflow state
            
        Returns:
            Dict[str, List[str]]: Dictionary of available input sources
        """
        input_sources = {
            'original_tables': workflow_state['available_tables'],
            'previous_products': []
        }
        
        # Add previous step products as potential inputs
        for prev_step, product_info in workflow_state['step_products'].items():
            if prev_step < step_num and product_info.get('table_name'):
                input_sources['previous_products'].append({
                    'table_name': product_info['table_name'],
                    'step': prev_step,
                    'reaction_name': product_info['reaction_name'],
                    'count': product_info['product_count']
                })
        
        return input_sources

    def _select_tables_for_step_bimolecular(self, input_sources: Dict[str, Any], 
                                        step_num: int, reaction_name: str) -> Dict[str, Any]:
        """
        Select tables for a bimolecular reaction step.
        
        Args:
            input_sources (Dict): Available input sources
            step_num (int): Current step number
            reaction_name (str): Name of the reaction
            
        Returns:
            Dict[str, Any]: Table configuration
        """
        try:
            print(f"\n🧬 BIMOLECULAR REACTION STEP {step_num} SETUP")
            print("=" * 60)
            print(f"Reaction: {reaction_name}")
            print("This reaction requires two reactant sources.")
            print("Note: The same table can be used for both primary and secondary reactants.")
            
            # Display available sources
            print(f"\n📋 Available Input Sources:")
            all_sources = []
            
            print(f"\n🗂️  Original Tables:")
            for i, table in enumerate(input_sources['original_tables'], 1):
                compound_count = self.get_compound_count(table_name=table)
                print(f"   {len(all_sources)+1}. {table} ({compound_count} compounds) [Original]")
                all_sources.append({
                    'name': table,
                    'type': 'original',
                    'count': compound_count
                })
            
            if input_sources['previous_products']:
                print(f"\n🔬 Previous Step Products:")
                for product_info in input_sources['previous_products']:
                    print(f"   {len(all_sources)+1}. {product_info['table_name']} "
                        f"({product_info['count']} products from Step {product_info['step']}: "
                        f"{product_info['reaction_name']}) [Products]")
                    all_sources.append({
                        'name': product_info['table_name'],
                        'type': 'products',
                        'count': product_info['count'],
                        'step': product_info['step']
                    })
            
            # Select primary reactant source
            print(f"\n📋 Select PRIMARY reactant source:")
            primary_idx = self._select_source_by_index(all_sources, "primary reactant")
            if primary_idx is None:
                return {}
            
            primary_source = all_sources[primary_idx]
            
            # Select secondary reactant source (ALLOW same table as primary)
            print(f"\n📋 Select SECONDARY reactant source:")
            print(f"   💡 You can select the same table as primary for intra-molecular reactions")
            secondary_idx = self._select_source_by_index(all_sources, "secondary reactant")
            if secondary_idx is None:
                return {}
            
            secondary_source = all_sources[secondary_idx]
            
            print(f"\n✅ Bimolecular Reaction Configuration:")
            print(f"   🅰️  Primary: '{primary_source['name']}' ({primary_source['count']} compounds)")
            print(f"   🅱️  Secondary: '{secondary_source['name']}' ({secondary_source['count']} compounds)")
            
            # Special handling when same table is used for both reactants
            if primary_source['name'] == secondary_source['name']:
                print(f"   🔄 Using same table for both reactants (intra-molecular reactions)")
                # For same table, combinations = n * (n-1) for non-self reactions, or n * n for self reactions
                same_table_combinations = primary_source['count'] * secondary_source['count']
                print(f"   🔢 Potential combinations: {same_table_combinations:,} (including self-reactions)")
                
                # Ask if user wants to exclude self-reactions
                exclude_self = input("   ❓ Exclude self-reactions (same compound with itself)? (y/n): ").strip().lower()
                exclude_self_reactions = exclude_self in ['y', 'yes']
                
                if exclude_self_reactions:
                    unique_combinations = primary_source['count'] * (primary_source['count'] - 1)
                    print(f"   🔢 Unique combinations (excluding self): {unique_combinations:,}")
                    combinations = unique_combinations
                else:
                    combinations = same_table_combinations
            else:
                # Estimate computational load for different tables
                combinations = primary_source['count'] * secondary_source['count']
                exclude_self_reactions = False  # Not applicable for different tables
            
            print(f"   🔢 Total combinations to process: {combinations:,}")
            
            if combinations > 1000000:
                print(f"   ⚠️  WARNING: Large number of combinations may take significant time!")
                proceed = input("   Continue anyway? (y/n): ").strip().lower()
                if proceed not in ['y', 'yes']:
                    return {}
            
            return {
                'type': 'bimolecular',
                'primary_table': primary_source['name'],
                'secondary_table': secondary_source['name'],
                'primary_source_type': primary_source['type'],
                'secondary_source_type': secondary_source['type'],
                'step_num': step_num,
                'same_table': primary_source['name'] == secondary_source['name'],
                'exclude_self_reactions': exclude_self_reactions
            }
            
        except Exception as e:
            print(f"❌ Error selecting tables for bimolecular step: {e}")
            return {}

    def _select_tables_for_step_unimolecular(self, input_sources: Dict[str, Any], 
                                            step_num: int, reaction_name: str) -> Dict[str, Any]:
        """
        Select tables for a unimolecular reaction step.
        
        Args:
            input_sources (Dict): Available input sources
            step_num (int): Current step number
            reaction_name (str): Name of the reaction
            
        Returns:
            Dict[str, Any]: Table configuration
        """
        try:
            print(f"\n🧪 UNIMOLECULAR REACTION STEP {step_num} SETUP")
            print("=" * 60)
            print(f"Reaction: {reaction_name}")
            print("This reaction processes compounds individually.")
            
            # Display available sources
            print(f"\n📋 Available Input Sources:")
            all_sources = []
            
            print(f"\n🗂️  Original Tables:")
            for table in input_sources['original_tables']:
                compound_count = self.get_compound_count(table_name=table)
                print(f"   {len(all_sources)+1}. {table} ({compound_count} compounds) [Original]")
                all_sources.append({
                    'name': table,
                    'type': 'original',
                    'count': compound_count
                })
            
            if input_sources['previous_products']:
                print(f"\n🔬 Previous Step Products:")
                for product_info in input_sources['previous_products']:
                    print(f"   {len(all_sources)+1}. {product_info['table_name']} "
                        f"({product_info['count']} products from Step {product_info['step']}: "
                        f"{product_info['reaction_name']}) [Products]")
                    all_sources.append({
                        'name': product_info['table_name'],
                        'type': 'products',
                        'count': product_info['count'],
                        'step': product_info['step']
                    })
            
            # Select input sources
            selected_sources = []
            
            if step_num == 1:
                # First step: allow multiple original tables
                print(f"\nSelect input sources (comma-separated numbers or 'all'):")
                selection = input("Selection: ").strip().lower()
                
                if selection == 'all':
                    selected_sources = all_sources
                else:
                    try:
                        indices = [int(x.strip()) - 1 for x in selection.split(',')]
                        selected_sources = [all_sources[i] for i in indices if 0 <= i < len(all_sources)]
                    except:
                        print("❌ Invalid selection")
                        return {}
            else:
                # Later steps: offer choice between original tables and previous products
                print(f"\nChoose input strategy:")
                print(f"   1. Use original tables")
                print(f"   2. Use products from previous steps")
                print(f"   3. Use both original tables and previous products")
                
                strategy = input("Enter choice (1-3): ").strip()
                
                if strategy == '1':
                    selected_sources = [s for s in all_sources if s['type'] == 'original']
                elif strategy == '2':
                    selected_sources = [s for s in all_sources if s['type'] == 'products']
                elif strategy == '3':
                    selected_sources = all_sources
                else:
                    print("❌ Invalid strategy selection")
                    return {}
            
            if not selected_sources:
                print("❌ No sources selected")
                return {}
            
            print(f"\n✅ Selected {len(selected_sources)} input source(s):")
            total_compounds = 0
            for source in selected_sources:
                source_type_label = "Products" if source['type'] == 'products' else "Original"
                print(f"   • {source['name']} ({source['count']} compounds) [{source_type_label}]")
                total_compounds += source['count']
            
            print(f"   📊 Total compounds to process: {total_compounds:,}")
            
            return {
                'type': 'unimolecular',
                'sources': selected_sources,
                'total_compounds': total_compounds,
                'step_num': step_num
            }
            
        except Exception as e:
            print(f"❌ Error selecting tables for unimolecular step: {e}")
            return {}

    def _select_source_by_index(self, sources: List[Dict], source_type: str) -> Optional[int]:
        """
        Select a source by index from available options.
        
        Args:
            sources (List[Dict]): List of available sources
            source_type (str): Description of source type
            
        Returns:
            Optional[int]: Selected index or None if cancelled
        """
        try:
            while True:
                selection = input(f"Select {source_type} (1-{len(sources)} or 'cancel'): ").strip()
                
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    return None
                
                try:
                    idx = int(selection) - 1
                    if 0 <= idx < len(sources):
                        return idx
                    else:
                        print(f"❌ Invalid selection. Please enter 1-{len(sources)}")
                except ValueError:
                    print("❌ Please enter a valid number or 'cancel'")
                    
        except KeyboardInterrupt:
            return None

    def _apply_step_bimolecular_reaction(self, table_config: Dict[str, Any], 
                                    reaction_info: Dict[str, Any], 
                                    workflow_name: str, step_num: int) -> Dict[str, Any]:
        """
        Apply a bimolecular reaction for a specific step.
        
        Args:
            table_config (Dict): Table configuration
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            step_num (int): Current step number
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            print(f"\n🧬 Executing Bimolecular Reaction - Step {step_num}")
            print("-" * 50)
            
            # Get compound data
            primary_df = self._get_table_as_dataframe(table_config['primary_table'])
            secondary_df = self._get_table_as_dataframe(table_config['secondary_table'])
            
            print(f"📊 Primary reactants: {len(primary_df)} compounds")
            print(f"📊 Secondary reactants: {len(secondary_df)} compounds")
            
            # Apply the bimolecular reaction
            products = self._apply_bimolecular_reaction(
                primary_df, secondary_df, reaction_info, workflow_name
            )
            
            # Save products if any were generated
            output_table_name = None
            if products:
                print(f"\n💾 Generated {len(products)} products from bimolecular reaction")
                
                while True:
                    save_choice = input("Save products to a new table? (y/n): ").strip().lower()
                    if save_choice in ['y', 'yes', 'n', 'no']:
                        break
                    print("❗ Please answer 'y'/'yes' or 'n'/'no'.")
                
                if save_choice in ['y', 'yes']:
                    # Ask user for a table name, use default if empty
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    default_name = f"step{step_num}_{reaction_info['name']}_{timestamp}"
                    user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                    chosen_name = user_table_name if user_table_name else default_name
                    output_table_name = self._sanitize_table_name(chosen_name)
                    
                    success = self._save_step_products(
                        products, output_table_name, workflow_name, step_num
                    )
                    
                    if not success:
                        output_table_name = None
                    
                    # Return now to avoid running the default-save block below
                    return {
                        'success': True,
                        'products_generated': len(products),
                        'output_table': output_table_name,
                        'reaction_type': 'bimolecular',
                        'step_num': step_num
                    }
                else:
                    # User chose not to save; return without creating a table
                    return {
                        'success': True,
                        'products_generated': len(products),
                        'output_table': None,
                        'reaction_type': 'bimolecular',
                        'step_num': step_num
                    }
            
        except Exception as e:
            print(f"❌ Error applying bimolecular reaction: {e}")
            return {
                'success': False,
                'products_generated': 0,
                'error': str(e)
            }

    def _apply_step_unimolecular_reaction(self, table_config: Dict[str, Any], 
                                        reaction_info: Dict[str, Any], 
                                        workflow_name: str, step_num: int) -> Dict[str, Any]:
        """
        Apply a unimolecular reaction for a specific step.
        
        Args:
            table_config (Dict): Table configuration
            reaction_info (Dict): Reaction information
            workflow_name (str): Workflow name
            step_num (int): Current step number
            
        Returns:
            Dict[str, Any]: Step execution results
        """
        try:
            print(f"\n🧪 Executing Unimolecular Reaction - Step {step_num}")
            print("-" * 50)
            
            all_products = []
            
            # Process each selected source
            for source in table_config['sources']:
                print(f"\n🔬 Processing source: '{source['name']}'")
                compounds_df = self._get_table_as_dataframe(source['name'])
                
                if compounds_df.empty:
                    print(f"   ⚠️  No compounds found, skipping...")
                    continue
                
                print(f"   📊 Processing {len(compounds_df)} compounds...")
                
                # Apply unimolecular reaction
                source_prefix = f"step{step_num}_{source['name']}_"
                products = self._apply_unimolecular_reaction(
                    compounds_df, reaction_info, workflow_name, source_prefix
                )
                
                all_products.extend(products)
                print(f"   ✅ Generated {len(products)} products")
            
            # Save products if any were generated
            output_table_name = None
            if all_products:
                print(f"\n💾 Total products generated: {len(all_products)}")
                
                while True:
                    save_choice = input("Save products to a new table? (y/n): ").strip().lower()
                    if save_choice in ['y', 'yes', 'n', 'no']:
                        break
                    print("❗ Please answer 'y'/'yes' or 'n'/'no'.")
                
                if save_choice in ['y', 'yes']:
                    # Ask user for a table name, use default if empty
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    default_name = f"step{step_num}_{reaction_info['name']}_{timestamp}"
                    user_table_name = input(f"Enter table name (default: {default_name}): ").strip()
                    chosen_name = user_table_name if user_table_name else default_name
                    output_table_name = self._sanitize_table_name(chosen_name)
                    
                    success = self._save_step_products(
                        all_products, output_table_name, workflow_name, step_num
                    )
                    
                    if not success:
                        output_table_name = None
            
            return {
                'success': True,
                'products_generated': len(all_products),
                'output_table': output_table_name,
                'reaction_type': 'unimolecular',
                'step_num': step_num
            }
            
        except Exception as e:
            print(f"❌ Error applying unimolecular reaction: {e}")
            return {
                'success': False,
                'products_generated': 0,
                'error': str(e)
            }

    def _save_step_products(self, products: List[Dict[str, Any]], table_name: str, 
                        workflow_name: str, step_num: int) -> bool:
        """
        Save products from a reaction step to a new table.
        
        Args:
            products (List[Dict]): List of product dictionaries
            table_name (str): Name for the new table
            workflow_name (str): Name of the workflow
            step_num (int): Current step number
            
        Returns:
            bool: True if saved successfully
        """
        try:
            if not products:
                return False
            
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create products table with step metadata
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {table_name} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT NOT NULL,
                    name TEXT,
                    flag TEXT,
                    workflow_step INTEGER,
                    workflow_name TEXT,
                    creation_date TEXT,
                    UNIQUE(smiles, name)
                )
            ''')
            
            # Create indexes
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_smiles ON {table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_step ON {table_name}(workflow_step)")
            
            # Insert products
            insert_query = f'''
                INSERT OR IGNORE INTO {table_name} 
                (smiles, name, flag, workflow_step, workflow_name, creation_date)
                VALUES (?, ?, ?, ?, ?, ?)
            '''
            
            creation_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            inserted_count = 0
            
            for product in products:
                try:
                    cursor.execute(insert_query, (
                        product['smiles'],
                        product['name'],
                        product.get('flag', 'step_product'),
                        step_num,
                        workflow_name,
                        creation_date
                    ))
                    inserted_count += 1
                except sqlite3.IntegrityError:
                    continue
            
            conn.commit()
            conn.close()
            
            print(f"   💾 Saved {inserted_count} products to table '{table_name}'")
            return True
            
        except Exception as e:
            print(f"   ❌ Error saving step products: {e}")
            return False

    def _display_reaction_workflow_summary(self, workflow_state: Dict[str, Any], workflow_name: str) -> None:
        """
        Display a comprehensive summary of the sequential workflow execution.
        
        Args:
            workflow_state (Dict): Final workflow state
            workflow_name (str): Name of the workflow
        """
        try:
            print(f"\n{'='*80}")
            print(f"🎯 SEQUENTIAL WORKFLOW SUMMARY: '{workflow_name}'")
            print(f"{'='*80}")
            
            total_products = 0
            successful_steps = 0
            
            for step_num, result in workflow_state['step_results'].items():
                status_icon = "✅" if result['success'] else "❌"
                products = result.get('products_generated', 0)
                total_products += products
                
                if result['success']:
                    successful_steps += 1
                
                print(f"\n{status_icon} Step {step_num}: {products:,} products generated")
                
                if result.get('output_table'):
                    print(f"   💾 Saved to: '{result['output_table']}'")
                
                if not result['success']:
                    print(f"   ❌ Error: {result.get('error', 'Unknown error')}")
            
            print(f"\n📊 FINAL RESULTS:")
            print(f"   🔬 Total steps: {len(workflow_state['step_results'])}")
            print(f"   ✅ Successful steps: {successful_steps}")
            print(f"   🧪 Total products generated: {total_products:,}")
            
            if workflow_state['step_products']:
                print(f"\n📋 Product Tables Created:")
                for step, product_info in workflow_state['step_products'].items():
                    print(f"   Step {step}: {product_info['table_name']} "
                        f"({product_info['product_count']} products)")
            
            print(f"{'='*80}")
            
        except Exception as e:
            print(f"❌ Error displaying workflow summary: {e}")
            
    def _query_and_delete_prev_step_tables(self, state: Dict[str, Any]) -> None:
                try:
                    step_products = state.get('step_products', {})
                    if not step_products:
                        print("\nℹ️  No step product tables recorded; nothing to delete.")
                        return

                    # Collect valid tables excluding the last step since it is considered to be the final products
                    tables = []
                    step_items = []
                    for step, info in step_products.items():
                        try:
                            step_int = int(step)
                        except Exception:
                            try:
                                step_int = int(str(step))
                            except Exception:
                                continue
                        tbl = info.get('table_name')
                        try:
                            cnt = int(info.get('product_count', 0) or 0)
                        except Exception:
                            cnt = 0
                        step_items.append((step_int, tbl, cnt))

                    if step_items:
                        # determine the highest step number and skip it
                        max_step = max(s for s, _, _ in step_items)
                        for s, tbl, cnt in sorted(step_items):
                            if s == max_step:
                                # skip the table corresponding to the last step
                                continue
                            if tbl:
                                tables.append((s, tbl, cnt))

                    if not tables:
                        print("\nℹ️  No named product tables found for previous steps.")
                        return

                    print("\n🗂️  Intermediate product tables from previous steps:")
                    for step, tbl, cnt in tables:
                        print(f"   • Step {step}: {tbl} ({cnt} products)")

                    
                    while True:
                        choice = input("\n❓ Delete intermediate product tables? (y = all / n = none / s = select): ").strip().lower()
                        if choice in ['y', 'yes', 'n', 'no', 's', 'select']:
                            break
                        print("❗ Please answer 'y'/'yes' , 'n'/'no' or 's'/'select'.")
                    
                    
                    if choice in ['', 'n', 'no']:
                        print("✅ Keeping all intermediate tables.")
                        return

                    to_delete = []
                    if choice in ['y', 'yes', 'all']:
                        to_delete = [tbl for _, tbl, _ in tables]
                    elif choice in ['s', 'select']:
                        sel = input("Enter step numbers or table names (comma-separated): ").strip()
                        if not sel:
                            print("⚠️  No selection provided. Aborting deletion.")
                            return
                        items = [s.strip() for s in sel.split(',') if s.strip()]
                        # map selection to table names
                        step_map = {str(step): tbl for step, tbl, _ in tables}
                        name_set = {tbl for _, tbl, _ in tables}
                        for it in items:
                            if it in step_map:
                                to_delete.append(step_map[it])
                            elif it in name_set:
                                to_delete.append(it)
                            else:
                                print(f"⚠️  Ignoring unknown selection: {it}")
                        to_delete = list(dict.fromkeys(to_delete))  # preserve order, remove duplicates
                    else:
                        print("⚠️  Unknown choice, aborting deletion.")
                        return

                    if not to_delete:
                        print("ℹ️  No tables selected for deletion.")
                        return

                    # Final confirmation
                    print("\n⚠️  Tables to be deleted:")
                    for tbl in to_delete:
                        print(f"   - {tbl}")
                    confirm = input(f"❗ Confirm deletion of {len(to_delete)} table(s)? This cannot be undone (yes/no): ").strip().lower()
                    if confirm not in ['yes', 'y']:
                        print("❌ Deletion cancelled by user.")
                        return

                    deleted = []
                    for tbl in to_delete:
                        try:
                            # Use confirm=False to avoid nested prompts; drop_table already prints status
                            success = self.drop_table(tbl, confirm=False)
                            if success:
                                deleted.append(tbl)
                        except Exception as e:
                            print(f"❌ Failed to delete table '{tbl}': {e}")

                    # Record deletions in workflow state
                    state.setdefault('deleted_step_tables', []).extend(deleted)

                    print(f"🗑️  Deleted {len(deleted)}/{len(to_delete)} requested tables.")
                except Exception as e:
                    print(f"⚠️  Error during deletion query: {e}")

    def drop_tables(self, confirm_each: bool = True) -> Dict[str, bool]:
        """
        Interactive method to list and delete multiple tables from the chemspace database.
        
        Args:
            confirm_each (bool): Whether to confirm each individual table deletion
            
        Returns:
            Dict[str, bool]: Dictionary mapping table names to deletion success status
        """
        try:
            # Get all available tables
            all_tables = self.get_all_tables()
            
            if not all_tables:
                print("📝 No tables found in chemspace database")
                return {}
            
            print(f"\n🗂️  CHEMSPACE DATABASE TABLES")
            print("=" * 80)
            print(f"📊 Total tables: {len(all_tables)}")
            print("=" * 80)
            
            # Display tables with details
            table_info = []
            for i, table_name in enumerate(all_tables, 1):
                try:
                    compound_count = self.get_compound_count(table_name=table_name)
                    table_type = self._classify_table_type(table_name)
                    table_info.append({
                        'index': i,
                        'name': table_name,
                        'count': compound_count,
                        'type': table_type
                    })
                    
                    type_icon = self._get_table_type_icon(table_type)
                    print(f"{i:3d}. {type_icon} {table_name:<35} ({compound_count:>6,} compounds) [{table_type}]")
                    
                except Exception as e:
                    print(f"{i:3d}. ❓ {table_name:<35} (Error: {e}) [Unknown]")
                    table_info.append({
                        'index': i,
                        'name': table_name,
                        'count': 0,
                        'type': 'unknown',
                        'error': str(e)
                    })
            
            print("=" * 80)
            
            # Interactive table selection
            print(f"\n🗑️  TABLE DELETION INTERFACE")
            print("-" * 50)
            print("Commands:")
            print("   • Enter table numbers (comma-separated): 1,3,5")
            print("   • Enter table names (comma-separated): table1,table2")
            print("   • Type 'all' to delete all tables")
            print("   • Type 'filter <type>' to delete by type: filter products")
            print("   • Type 'range <start>-<end>' to delete range: range 1-5")
            print("   • Type 'show' to display table list again")
            print("   • Type 'cancel' to abort")
            print("-" * 50)
            
            while True:
                try:
                    user_input = input("\n🔍 Select tables to delete: ").strip()
                    
                    if not user_input or user_input.lower() == 'cancel':
                        print("❌ Table deletion cancelled by user")
                        return {}
                    elif user_input.lower() == 'show':
                        self._redisplay_table_list(table_info)
                        continue
                    
                    # Parse user selection
                    selected_tables = self._parse_table_selection(user_input, table_info)
                    
                    if not selected_tables:
                        print("❌ No valid tables selected. Try again or type 'cancel'")
                        continue
                    
                    # Confirm selection
                    if not self._confirm_table_deletion(selected_tables, all_tables):
                        continue
                    
                    # Perform deletions
                    return self._execute_table_deletions(selected_tables, confirm_each)
                    
                except KeyboardInterrupt:
                    print("\n❌ Table deletion cancelled by user")
                    return {}
                except Exception as e:
                    print(f"❌ Error in table selection: {e}")
                    continue
        
        except Exception as e:
            print(f"❌ Error in drop_tables method: {e}")
            return {}

    def _classify_table_type(self, table_name: str) -> str:
        """
        Classify table type based on naming patterns.
        
        Args:
            table_name (str): Name of the table
            
        Returns:
            str: Table type classification
        """
        name_lower = table_name.lower()
        
        if '_products_' in name_lower:
            return 'reaction_products'
        elif '_filtered_' in name_lower:
            return 'filtered_compounds'
        elif 'step' in name_lower and any(x in name_lower for x in ['reaction', 'product']):
            return 'step_products'
        elif name_lower.startswith('temp_') or name_lower.startswith('tmp_'):
            return 'temporary'
        elif any(word in name_lower for word in ['workflow', 'batch', 'processed']):
            return 'processed'
        else:
            return 'original'

    def _get_table_type_icon(self, table_type: str) -> str:
        """
        Get emoji icon for table type.
        
        Args:
            table_type (str): Table type
            
        Returns:
            str: Emoji icon
        """
        icons = {
            'original': '📋',
            'reaction_products': '🧪',
            'filtered_compounds': '🔍',
            'step_products': '⚗️',
            'processed': '🔄',
            'temporary': '📄',
            'unknown': '❓'
        }
        return icons.get(table_type, '📊')

    def _redisplay_table_list(self, table_info: List[Dict]) -> None:
        """
        Redisplay the table list.
        
        Args:
            table_info (List[Dict]): Table information list
        """
        print(f"\n📋 TABLE LIST:")
        print("-" * 80)
        for info in table_info:
            type_icon = self._get_table_type_icon(info['type'])
            if 'error' in info:
                print(f"{info['index']:3d}. ❓ {info['name']:<35} (Error) [Unknown]")
            else:
                print(f"{info['index']:3d}. {type_icon} {info['name']:<35} ({info['count']:>6,} compounds) [{info['type']}]")
        print("-" * 80)

    def _parse_table_selection(self, user_input: str, table_info: List[Dict]) -> List[str]:
        """
        Parse user input to determine selected tables.
        
        Args:
            user_input (str): User input string
            table_info (List[Dict]): Table information list
            
        Returns:
            List[str]: List of selected table names
        """
        try:
            user_input = user_input.strip().lower()
            selected_tables = []
            
            if user_input == 'all':
                # Select all tables
                selected_tables = [info['name'] for info in table_info]
                print(f"✅ Selected all {len(selected_tables)} tables")
                
            elif user_input.startswith('filter '):
                # Filter by table type
                table_type = user_input[7:].strip()
                type_mapping = {
                    'original': 'original',
                    'products': 'reaction_products',
                    'reaction': 'reaction_products',
                    'filtered': 'filtered_compounds',
                    'step': 'step_products',
                    'processed': 'processed',
                    'temp': 'temporary',
                    'temporary': 'temporary',
                    'unknown': 'unknown'
                }
                
                target_type = type_mapping.get(table_type)
                if target_type:
                    selected_tables = [info['name'] for info in table_info if info['type'] == target_type]
                    print(f"✅ Selected {len(selected_tables)} tables of type '{table_type}'")
                else:
                    print(f"❌ Unknown table type: '{table_type}'")
                    print(f"   Available types: {', '.join(type_mapping.keys())}")
                    
            elif user_input.startswith('range '):
                # Range selection
                range_part = user_input[6:].strip()
                if '-' in range_part:
                    try:
                        start_str, end_str = range_part.split('-', 1)
                        start_idx = int(start_str.strip())
                        end_idx = int(end_str.strip())
                        
                        if 1 <= start_idx <= len(table_info) and 1 <= end_idx <= len(table_info):
                            if start_idx <= end_idx:
                                selected_tables = [table_info[i-1]['name'] for i in range(start_idx, end_idx + 1)]
                                print(f"✅ Selected range {start_idx}-{end_idx}: {len(selected_tables)} tables")
                            else:
                                print("❌ Invalid range: start must be <= end")
                        else:
                            print(f"❌ Range out of bounds. Valid range: 1-{len(table_info)}")
                    except ValueError:
                        print("❌ Invalid range format. Use: range 1-5")
                else:
                    print("❌ Invalid range format. Use: range 1-5")
                    
            else:
                # Parse comma-separated values (numbers or names)
                items = [item.strip() for item in user_input.split(',') if item.strip()]
                
                for item in items:
                    try:
                        # Try as table number first
                        table_idx = int(item) - 1
                        if 0 <= table_idx < len(table_info):
                            table_name = table_info[table_idx]['name']
                            if table_name not in selected_tables:
                                selected_tables.append(table_name)
                        else:
                            print(f"⚠️  Invalid table number: {item} (valid: 1-{len(table_info)})")
                    except ValueError:
                        # Treat as table name
                        matching_tables = [info['name'] for info in table_info if info['name'].lower() == item.lower()]
                        if matching_tables:
                            for table_name in matching_tables:
                                if table_name not in selected_tables:
                                    selected_tables.append(table_name)
                        else:
                            print(f"⚠️  Table not found: '{item}'")
            
            return selected_tables
            
        except Exception as e:
            print(f"❌ Error parsing table selection: {e}")
            return []

    def _confirm_table_deletion(self, selected_tables: List[str], all_tables: List[str]) -> bool:
        """
        Confirm table deletion with user.
        
        Args:
            selected_tables (List[str]): Tables selected for deletion
            all_tables (List[str]): All available tables
            
        Returns:
            bool: True if confirmed, False otherwise
        """
        try:
            print(f"\n⚠️  DELETION CONFIRMATION")
            print("=" * 50)
            print(f"🗑️  Tables to delete: {len(selected_tables)}")
            print(f"📊 Remaining tables: {len(all_tables) - len(selected_tables)}")
            print("=" * 50)
            
            # Show tables to be deleted with their compound counts
            total_compounds = 0
            for i, table_name in enumerate(selected_tables, 1):
                try:
                    compound_count = self.get_compound_count(table_name=table_name)
                    table_type = self._classify_table_type(table_name)
                    type_icon = self._get_table_type_icon(table_type)
                    total_compounds += compound_count
                    
                    print(f"   {i:2d}. {type_icon} {table_name} ({compound_count:,} compounds)")
                except Exception:
                    print(f"   {i:2d}. ❓ {table_name} (unknown count)")
            
            print("=" * 50)
            print(f"📊 Total compounds to be deleted: {total_compounds:,}")
            print("🚨 WARNING: This action cannot be undone!")
            print("=" * 50)
            
            # Multiple confirmation steps
            confirm1 = input(f"❓ Delete {len(selected_tables)} tables? (yes/no): ").strip().lower()
            if confirm1 not in ['yes', 'y']:
                print("❌ Deletion cancelled")
                return False
            
            # Show critical warnings for important table types
            critical_tables = []
            for table_name in selected_tables:
                table_type = self._classify_table_type(table_name)
                if table_type == 'original':
                    critical_tables.append(table_name)
            
            if critical_tables:
                print(f"\n🚨 CRITICAL WARNING:")
                print(f"   You are about to delete {len(critical_tables)} ORIGINAL compound table(s)!")
                print(f"   This will permanently remove your source data!")
                
                for table in critical_tables:
                    print(f"   📋 {table}")
                
                confirm2 = input(f"\n❗ Type 'DELETE ORIGINAL TABLES' to confirm: ").strip()
                if confirm2 != 'DELETE ORIGINAL TABLES':
                    print("❌ Critical table deletion cancelled")
                    return False
            
            # Final confirmation with count
            confirm3 = input(f"❓ Type the number {len(selected_tables)} to proceed: ").strip()
            if confirm3 != str(len(selected_tables)):
                print("❌ Count confirmation failed. Deletion cancelled.")
                return False
            
            return True
            
        except KeyboardInterrupt:
            print("\n❌ Deletion cancelled by user")
            return False
        except Exception as e:
            print(f"❌ Error in deletion confirmation: {e}")
            return False

    def _execute_table_deletions(self, selected_tables: List[str], confirm_each: bool) -> Dict[str, bool]:
        """
        Execute the deletion of selected tables.
        
        Args:
            selected_tables (List[str]): Tables to delete
            confirm_each (bool): Whether to confirm each deletion
            
        Returns:
            Dict[str, bool]: Deletion results
        """
        try:
            results = {}
            successful_deletions = 0
            failed_deletions = 0
            
            print(f"\n🗑️  EXECUTING TABLE DELETIONS")
            print("=" * 50)
            
            # Use progress bar if many tables
            if TQDM_AVAILABLE and len(selected_tables) > 5:
                progress_bar = tqdm(
                    selected_tables,
                    desc="Deleting tables",
                    unit="tables"
                )
            else:
                progress_bar = selected_tables
            
            for table_name in progress_bar:
                try:
                    print(f"\n🗑️  Deleting table: '{table_name}'")
                    
                    # Use existing drop_table method with appropriate confirmation
                    success = self.drop_table(table_name, confirm=confirm_each)
                    results[table_name] = success
                    
                    if success:
                        successful_deletions += 1
                        if not TQDM_AVAILABLE:
                            print(f"   ✅ Successfully deleted '{table_name}'")
                    else:
                        failed_deletions += 1
                        if not TQDM_AVAILABLE:
                            print(f"   ❌ Failed to delete '{table_name}'")
                    
                except Exception as e:
                    results[table_name] = False
                    failed_deletions += 1
                    print(f"   ❌ Error deleting '{table_name}': {e}")
            
            # Display summary
            print(f"\n📊 DELETION SUMMARY")
            print("=" * 30)
            print(f"✅ Successfully deleted: {successful_deletions}")
            print(f"❌ Failed deletions: {failed_deletions}")
            print(f"📊 Total processed: {len(selected_tables)}")
            
            if successful_deletions > 0:
                print(f"🗑️  {successful_deletions} table(s) have been permanently removed")
            
            if failed_deletions > 0:
                print(f"\n❌ Failed deletions:")
                for table_name, success in results.items():
                    if not success:
                        print(f"   • {table_name}")
            
            print("=" * 30)
            
            return results
            
        except Exception as e:
            print(f"❌ Error executing table deletions: {e}")
            return {table: False for table in selected_tables}

    def drop_table(self, table_name: Optional[str] = None, confirm: bool = True) -> bool:
        """
        Drop a single table from the chemspace database.
        
        Args:
            table_name (Optional[str]): Name of the table to drop. If None, prompts user to select one.
            confirm (bool): Whether to ask for confirmation before deletion
            
        Returns:
            bool: True if table was dropped successfully, False otherwise
        """
        try:
            # If no table name provided, prompt user to select one
            if not table_name or not str(table_name).strip():
                tables = self.get_all_tables()
                if not tables:
                    print("📝 No tables found in chemspace database")
                    return False

                print("\n📊 Available tables:")
                # Show tables with type icon for clarity
                for idx, t in enumerate(tables, start=1):
                    try:
                        t_type = self._classify_table_type(t)
                        t_icon = self._get_table_type_icon(t_type)
                    except Exception:
                        t_type, t_icon = "unknown", "❓"
                    print(f"  [{idx}] {t_icon} {t} ({t_type})")

                while True:
                    choice = input("\nSelect table to drop (number/name) or 'q' to cancel: ").strip()
                    if choice.lower() in ['q', 'quit', 'cancel', 'n']:
                        print("❌ Table deletion cancelled.")
                        return False
                    if choice.isdigit():
                        idx = int(choice) - 1
                        if 0 <= idx < len(tables):
                            table_name = tables[idx]
                            break
                        else:
                            print("❗ Invalid selection. Please enter a valid number or name.")
                            continue
                    else:
                        # Name provided: validate it exists
                        if choice in tables:
                            table_name = choice
                            break
                        else:
                            print(f"❗ Table '{choice}' not found. Please enter a valid number or name.")
                            continue

            # Check if table exists
            tables = self.get_all_tables()
            if table_name not in tables:
                print(f"❌ Table '{table_name}' does not exist in chemspace database")
                return False
            
            # Get table information for confirmation
            try:
                compound_count = self.get_compound_count(table_name=table_name)
                table_type = self._classify_table_type(table_name)
                type_icon = self._get_table_type_icon(table_type)
            except Exception:
                compound_count = 0
                table_type = 'unknown'
                type_icon = '❓'
            
            # Show confirmation if requested
            if confirm:
                print(f"\n🗑️  TABLE DELETION CONFIRMATION")
                print("-" * 40)
                print(f"{type_icon} Table: '{table_name}'")
                print(f"📊 Compounds: {compound_count:,}")
                print(f"🏷️  Type: {table_type}")
                print(f"⚠️  This action cannot be undone!")
                print("-" * 40)
                
                # Get user confirmation
                while True:
                    user_confirm = input(f"❓ Delete table '{table_name}'? (y/n): ").strip().lower()
                    if user_confirm in ['y', 'yes']:
                        break
                    elif user_confirm in ['n', 'no']:
                        print(f"❌ Table deletion cancelled for '{table_name}'")
                        return False
                    else:
                        print("❗ Please answer 'y'/'yes' or 'n'/'no'")
            
            # Connect to database and drop table
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Drop associated indexes first (if any)
            index_prefixes = [f"idx_{table_name}_", f"index_{table_name}_"]
            for prefix in index_prefixes:
                try:
                    # Get all indexes for this table
                    cursor.execute("SELECT name FROM sqlite_master WHERE type='index' AND name LIKE ?", (f"{prefix}%",))
                    indexes = cursor.fetchall()
                    
                    for (index_name,) in indexes:
                        try:
                            cursor.execute(f"DROP INDEX IF EXISTS {index_name}")
                        except Exception:
                            pass  # Continue even if index drop fails
                except Exception:
                    pass  # Continue if index query fails
            
            # Drop the table
            cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
            
            # Verify table was dropped
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name = ?", (table_name,))
            table_still_exists = cursor.fetchone() is not None
            
            conn.commit()
            conn.close()
            
            if not table_still_exists:
                print(f"✅ Successfully dropped table '{table_name}' ({compound_count:,} compounds)")
                return True
            else:
                print(f"❌ Failed to drop table '{table_name}' (table still exists)")
                return False
            
        except sqlite3.Error as e:
            print(f"❌ Database error dropping table '{table_name}': {e}")
            return False
        except Exception as e:
            print(f"❌ Error dropping table '{table_name}': {e}")
            return False

    def clear_all_tables(self, confirm_with_project_name: bool = True) -> bool:
        """
        Delete ALL tables from the chemspace database.
        
        Args:
            confirm_with_project_name (bool): Require typing project name for confirmation
            
        Returns:
            bool: True if all tables were cleared successfully
        """
        try:
            all_tables = self.get_all_tables()
            
            if not all_tables:
                print("📝 No tables found in chemspace database")
                return True
            
            # Calculate total compounds
            total_compounds = 0
            for table_name in all_tables:
                try:
                    total_compounds += self.get_compound_count(table_name=table_name)
                except:
                    pass
            
            print(f"🚨 CLEAR ALL CHEMSPACE TABLES")
            print("=" * 60)
            print(f"⚠️  WARNING: This will delete ALL {len(all_tables)} tables!")
            print(f"📊 Total compounds to be lost: {total_compounds:,}")
            print(f"🗂️  Project: {self.name}")
            print(f"🚨 THIS ACTION CANNOT BE UNDONE!")
            print("=" * 60)
            
            # Multiple confirmation steps
            confirm1 = input(f"❓ Delete ALL {len(all_tables)} tables? (yes/no): ").strip().lower()
            if confirm1 not in ['yes', 'y']:
                print("❌ Operation cancelled")
                return False
            
            if confirm_with_project_name:
                confirm2 = input(f"❓ Type the project name '{self.name}' to confirm: ").strip()
                if confirm2 != self.name:
                    print("❌ Project name doesn't match. Operation cancelled.")
                    return False
            
            confirm3 = input("❓ Type 'CLEAR ALL CHEMSPACE TABLES' in capitals to proceed: ").strip()
            if confirm3 != 'CLEAR ALL CHEMSPACE TABLES':
                print("❌ Confirmation phrase doesn't match. Operation cancelled.")
                return False
            
            # Execute mass deletion
            print(f"\n🗑️  Clearing all chemspace tables...")
            
            deleted_count = 0
            failed_count = 0
            
            for table_name in all_tables:
                try:
                    success = self.drop_table(table_name, confirm=False)
                    if success:
                        deleted_count += 1
                    else:
                        failed_count += 1
                except Exception as e:
                    print(f"❌ Error deleting '{table_name}': {e}")
                    failed_count += 1
            
            print(f"\n✅ Chemspace database cleared!")
            print(f"   🗑️  Tables deleted: {deleted_count}")
            print(f"   ❌ Deletion failures: {failed_count}")
            print(f"   📊 Compounds removed: {total_compounds:,}")
            
            return failed_count == 0
            
        except Exception as e:
            print(f"❌ Error clearing all tables: {e}")
            return False

    def delete_filtering_workflow(self, workflow_identifier: Optional[str] = None) -> bool:
        """
        Delete a filtering workflow from the chemspace database.
        Can delete by workflow name or ID, or show interactive selection.
        
        Args:
            workflow_identifier (Optional[str]): Workflow name or ID to delete. If None, shows interactive selection.
            
        Returns:
            bool: True if workflow was deleted successfully, False otherwise
        """
        try:
            # Check if filtering_workflows table exists
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='filtering_workflows'")
            if not cursor.fetchone():
                print("📝 No filtering workflows table found - no workflows to delete")
                conn.close()
                return False
            
            # Get all available workflows
            cursor.execute("""
                SELECT id, workflow_name, creation_date, description, filters_dict, filter_count, total_instances
                FROM filtering_workflows 
                ORDER BY creation_date DESC
            """)
            workflows = cursor.fetchall()
            
            if not workflows:
                print("📝 No filtering workflows found to delete")
                conn.close()
                return False
            
            # If no identifier provided, show interactive selection
            if workflow_identifier is None:
                workflow_to_delete = self._select_filtering_workflow_for_deletion(workflows)
                if workflow_to_delete is None:
                    conn.close()
                    return False
            else:
                # Find workflow by name or ID using existing method (adapted for filtering workflows)
                workflow_to_delete = self._find_filtering_workflow_by_identifier(workflows, workflow_identifier)
                if workflow_to_delete is None:
                    print(f"❌ No filtering workflow found with identifier '{workflow_identifier}'")
                    self._show_available_filtering_workflows()
                    conn.close()
                    return False
            
            # Extract workflow information
            workflow_id, workflow_name, creation_date, description, filters_dict_str, filter_count, total_instances = workflow_to_delete
            
            # Parse filters for display
            try:
                filters_dict = json.loads(filters_dict_str)
                actual_filter_count = len(filters_dict)
            except json.JSONDecodeError:
                filters_dict = {}
                actual_filter_count = filter_count or 0
            
            # Show workflow details and confirm deletion
            print(f"\n⚠️  FILTERING WORKFLOW DELETION CONFIRMATION")
            print("=" * 70)
            print(f"🆔 ID: {workflow_id}")
            print(f"📋 Name: '{workflow_name}'")
            print(f"📅 Created: {creation_date}")
            print(f"🔍 Filters: {actual_filter_count}")
            print(f"🔢 Total instances: {total_instances or 0}")
            print(f"📄 Description: {description or 'No description'}")
            
            if filters_dict:
                print(f"\n🔍 Filters in this workflow:")
                for i, (filter_name, instances) in enumerate(list(filters_dict.items())[:5], 1):  # Show first 5
                    print(f"   {i}. {filter_name} (instances: {instances})")
                if len(filters_dict) > 5:
                    print(f"   ... and {len(filters_dict) - 5} more filters")
            
            print("=" * 70)
            
            # Double confirmation for safety
            confirm1 = input(f"❓ Are you sure you want to delete filtering workflow '{workflow_name}'? (yes/no): ").strip().lower()
            if confirm1 not in ['yes', 'y']:
                print("❌ Filtering workflow deletion cancelled")
                conn.close()
                return False
            
            confirm2 = input(f"❓ Type the workflow name '{workflow_name}' to confirm deletion: ").strip()
            if confirm2 != workflow_name:
                print("❌ Workflow name doesn't match. Deletion cancelled for safety.")
                conn.close()
                return False
            
            # Perform deletion
            print(f"🗑️  Deleting filtering workflow '{workflow_name}'...")
            
            cursor.execute("DELETE FROM filtering_workflows WHERE id = ?", (workflow_id,))
            deleted_rows = cursor.rowcount
            
            conn.commit()
            conn.close()
            
            if deleted_rows > 0:
                print(f"✅ Successfully deleted filtering workflow!")
                print(f"   📋 Workflow: '{workflow_name}' (ID: {workflow_id})")
                print(f"   🔍 Filters removed: {actual_filter_count}")
                print(f"   🔢 Total instances removed: {total_instances or 0}")
                print(f"   📅 Original creation date: {creation_date}")
                
                # Show remaining workflows using existing method
                remaining_count = self._count_remaining_filtering_workflows()
                if remaining_count > 0:
                    print(f"   📊 Remaining workflows: {remaining_count}")
                else:
                    print("   📝 No filtering workflows remaining in database")
                
                return True
            else:
                print(f"❌ Failed to delete workflow (no rows affected)")
                return False
            
        except Exception as e:
            print(f"❌ Error deleting filtering workflow: {e}")
            return False

    def _select_filtering_workflow_for_deletion(self, workflows: List[Tuple]) -> Optional[Tuple]:
        """
        Interactive filtering workflow selection for deletion.
        Uses similar pattern to existing _select_workflow_for_deletion but adapted for filtering workflows.
        
        Args:
            workflows (List[Tuple]): List of workflow tuples from database
            
        Returns:
            Optional[Tuple]: Selected workflow tuple or None if cancelled
        """
        try:
            print(f"\n🗑️  SELECT FILTERING WORKFLOW FOR DELETION")
            print("=" * 70)
            print(f"📋 Available Filtering Workflows ({len(workflows)} total):")
            print("-" * 70)
            print(f"{'#':<3} {'Name':<25} {'Created':<12} {'Filters':<8} {'Instances':<10}")
            print("-" * 70)
            
            workflow_map = {}
            for i, (wf_id, name, date, desc, filters_dict_str, filter_count, total_instances) in enumerate(workflows, 1):
                date_short = date[:10] if date else "Unknown"
                name_display = name[:24] if len(name) <= 24 else name[:21] + "..."
                
                # Try to get actual filter count from JSON
                try:
                    filters_dict = json.loads(filters_dict_str)
                    actual_filter_count = len(filters_dict)
                except:
                    actual_filter_count = filter_count or 0
                
                print(f"{i:<3} {name_display:<25} {date_short:<12} {actual_filter_count:<8} {total_instances or 0:<10}")
                workflow_map[str(i)] = workflows[i-1]
            
            print("-" * 70)
            print("Commands: Enter workflow number, 'cancel' to abort, 'list' for details")
            
            while True:
                selection = input("\n🔍 Select filtering workflow to delete: ").strip().lower()
                
                if selection in ['cancel', 'abort', 'exit', 'quit']:
                    print("❌ Filtering workflow deletion cancelled by user")
                    return None
                elif selection == 'list':
                    # Reuse existing detailed workflow display method
                    self._show_detailed_filtering_workflow_list(workflows)
                    continue
                elif selection in workflow_map:
                    return workflow_map[selection]
                else:
                    try:
                        selection_int = int(selection)
                        if 1 <= selection_int <= len(workflows):
                            return workflows[selection_int - 1]
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(workflows)}")
                    except ValueError:
                        print("❌ Invalid input. Please enter a number, 'list', or 'cancel'")
            
        except KeyboardInterrupt:
            print("\n❌ Selection cancelled by user")
            return None
        except Exception as e:
            print(f"❌ Error in filtering workflow selection: {e}")
            return None

    def _find_filtering_workflow_by_identifier(self, workflows: List[Tuple], identifier: str) -> Optional[Tuple]:
        """
        Find filtering workflow by name or ID. 
        Reuses the same logic as existing _find_workflow_by_identifier but adapted for filtering workflows.
        
        Args:
            workflows (List[Tuple]): List of workflow tuples
            identifier (str): Workflow name or ID to search for
            
        Returns:
            Optional[Tuple]: Matching workflow tuple or None if not found
        """
        try:
            # Try to find by ID first
            try:
                search_id = int(identifier)
                for workflow in workflows:
                    if workflow[0] == search_id:  # workflow[0] is the ID
                        return workflow
            except ValueError:
                pass  # Not a number, continue to name search
            
            # Search by name (case-insensitive)
            identifier_lower = identifier.lower()
            for workflow in workflows:
                workflow_name = workflow[1].lower()  # workflow[1] is the name
                if workflow_name == identifier_lower:
                    return workflow
            
            # Partial name match if exact match not found
            for workflow in workflows:
                workflow_name = workflow[1].lower()
                if identifier_lower in workflow_name:
                    return workflow
            
            return None
            
        except Exception as e:
            print(f"❌ Error searching for filtering workflow: {e}")
            return None

    def _show_detailed_filtering_workflow_list(self, workflows: List[Tuple]) -> None:
        """
        Show detailed information about available filtering workflows.
        Adapted from existing _show_detailed_workflow_list method.
        
        Args:
            workflows (List[Tuple]): List of workflow tuples
        """
        try:
            print(f"\n📋 DETAILED FILTERING WORKFLOW INFORMATION")
            print("=" * 90)
            
            for i, (wf_id, name, date, desc, filters_dict_str, filter_count, total_instances) in enumerate(workflows, 1):
                try:
                    filters_dict = json.loads(filters_dict_str)
                    actual_filter_count = len(filters_dict)
                    
                    # Get filter names for preview
                    filter_names = list(filters_dict.keys())
                    filters_preview = ', '.join(filter_names[:3])
                    if len(filter_names) > 3:
                        filters_preview += f" and {len(filter_names) - 3} more..."
                except json.JSONDecodeError:
                    actual_filter_count = filter_count or 0
                    filters_preview = "Error parsing filters"
                
                print(f"\n{i}. {name} (ID: {wf_id})")
                print(f"   📅 Created: {date}")
                print(f"   🔍 Filters: {actual_filter_count}")
                print(f"   🔢 Total instances: {total_instances or 0}")
                print(f"   📄 Description: {desc or 'No description'}")
                print(f"   🔍 Filter names: {filters_preview}")
            
            print("=" * 90)
            
        except Exception as e:
            print(f"❌ Error showing detailed filtering workflow list: {e}")

    def _show_available_filtering_workflows(self) -> None:
        """
        Display available filtering workflows with their IDs.
        Reuses the pattern from existing _show_available_reaction_workflows.
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("""
            SELECT id, workflow_name, creation_date, description, filter_count, total_instances
            FROM filtering_workflows 
            ORDER BY creation_date DESC
            """)
            
            workflows = cursor.fetchall()
            conn.close()
            
            if workflows:
                print("\n📋 Available Filtering Workflows:")
                print("-" * 80)
                print(f"{'ID':<4} {'Name':<25} {'Created':<12} {'Filters':<8} {'Instances':<10} {'Description':<25}")
                print("-" * 80)
                
                for wf_id, name, date, desc, filter_count, total_instances in workflows:
                    date_short = date[:10] if date else "Unknown"
                    desc_short = (desc or "No description")[:22] + "..." if len(str(desc or "")) > 25 else (desc or "No description")
                    print(f"{wf_id:<4} {name[:24]:<25} {date_short:<12} {filter_count or 0:<8} {total_instances or 0:<10} {desc_short:<25}")
                
                print("-" * 80)
            else:
                print("📝 No filtering workflows found.")
                
        except Exception as e:
            print(f"❌ Error showing available filtering workflows: {e}")

    def _count_remaining_filtering_workflows(self) -> int:
        """
        Count remaining filtering workflows after deletion.
        Reuses the pattern from existing _count_remaining_workflows.
        
        Returns:
            int: Number of remaining workflows
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT COUNT(*) FROM filtering_workflows")
            count = cursor.fetchone()[0]
            
            conn.close()
            return count
            
        except Exception:
            return 0

    def delete_multiple_filtering_workflows(self, workflow_identifiers: List[str]) -> Dict[str, bool]:
        """
        Delete multiple filtering workflows in batch.
        Reuses the exact pattern from existing delete_multiple_reaction_workflows.
        
        Args:
            workflow_identifiers (List[str]): List of workflow names or IDs to delete
            
        Returns:
            Dict[str, bool]: Dictionary mapping workflow identifiers to deletion success status
        """
        try:
            print(f"🗑️  BATCH FILTERING WORKFLOW DELETION")
            print("=" * 60)
            print(f"📋 Workflows to delete: {len(workflow_identifiers)}")
            
            results = {}
            successful_deletions = 0
            failed_deletions = 0
            
            for identifier in workflow_identifiers:
                print(f"\n🔍 Processing filtering workflow: '{identifier}'")
                success = self.delete_filtering_workflow(identifier)
                results[identifier] = success
                
                if success:
                    successful_deletions += 1
                else:
                    failed_deletions += 1
            
            # Summary
            print(f"\n📊 BATCH DELETION SUMMARY")
            print("=" * 40)
            print(f"✅ Successful deletions: {successful_deletions}")
            print(f"❌ Failed deletions: {failed_deletions}")
            print(f"📋 Total processed: {len(workflow_identifiers)}")
            
            return results
            
        except Exception as e:
            print(f"❌ Error in batch filtering workflow deletion: {e}")
            return {identifier: False for identifier in workflow_identifiers}

    def clear_all_filtering_workflows(self, confirm_with_count: bool = True) -> bool:
        """
        Delete all filtering workflows from the database.
        Reuses the exact pattern from existing clear_all_reaction_workflows.
        
        Args:
            confirm_with_count (bool): Require user to enter the exact count for confirmation
            
        Returns:
            bool: True if all workflows were cleared successfully
        """
        try:
            # Check if table exists and get workflow count
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='filtering_workflows'")
            if not cursor.fetchone():
                print("📝 No filtering workflows table found")
                conn.close()
                return True
            
            cursor.execute("SELECT COUNT(*) FROM filtering_workflows")
            workflow_count = cursor.fetchone()[0]
            
            if workflow_count == 0:
                print("📝 No filtering workflows found to clear")
                conn.close()
                return True
            
            # Show warning and get confirmation
            print(f"⚠️  CLEAR ALL FILTERING WORKFLOWS")
            print("=" * 60)
            print(f"🚨 WARNING: This will permanently delete ALL {workflow_count} filtering workflows!")
            print(f"🚨 This action cannot be undone!")
            print("=" * 60)
            
            # First confirmation
            confirm1 = input(f"❓ Are you sure you want to delete all {workflow_count} filtering workflows? (yes/no): ").strip().lower()
            if confirm1 not in ['yes', 'y']:
                print("❌ Operation cancelled")
                conn.close()
                return False
            
            # Second confirmation with count
            if confirm_with_count:
                confirm2 = input(f"❓ Type the number {workflow_count} to confirm deletion of all filtering workflows: ").strip()
                if confirm2 != str(workflow_count):
                    print("❌ Count doesn't match. Operation cancelled for safety.")
                    conn.close()
                    return False
            
            # Third confirmation
            confirm3 = input("❓ Type 'DELETE ALL FILTERING WORKFLOWS' in capitals to proceed: ").strip()
            if confirm3 != 'DELETE ALL FILTERING WORKFLOWS':
                print("❌ Confirmation phrase doesn't match. Operation cancelled.")
                conn.close()
                return False
            
            # Perform deletion
            print(f"🗑️  Clearing all filtering workflows...")
            cursor.execute("DELETE FROM filtering_workflows")
            deleted_count = cursor.rowcount
            
            conn.commit()
            conn.close()
            
            print(f"✅ Successfully cleared all filtering workflows!")
            print(f"   🗑️  Workflows deleted: {deleted_count}")
            print(f"   📊 Database is now empty of filtering workflows")
            
            return True
            
        except Exception as e:
            print(f"❌ Error clearing all filtering workflows: {e}")
            return False

    def apply_ersilia_model(self) -> Optional[List[str]]:
        """
        Apply an Ersilia model by selecting a table and returning SMILES strings.
        Prompts user for table selection, model name, and returns SMILES as a list.
        
        Returns:
            Optional[List[str]]: List of SMILES strings from selected table, or None if cancelled
        """
        try:
            # Check if Ersilia is available
            try:
                from ersilia.api import Model
            except ImportError:
                print("❌ Ersilia not installed. Please install Ersilia:")
                print("   pip install ersilia")
                return None
            
            # Get available tables
            available_tables = self.get_all_tables()
            
            if not available_tables:
                print("❌ No tables found in chemspace database")
                return None
            
            print(f"\n🧬 ERSILIA MODEL APPLICATION")
            print("=" * 50)
            print("Select a table to process with Ersilia models:")
            print("=" * 50)
            
            # Display available tables
            table_info = []
            for i, table_name in enumerate(available_tables, 1):
                try:
                    compound_count = self.get_compound_count(table_name=table_name)
                    table_type = self._classify_table_type(table_name)
                    type_icon = self._get_table_type_icon(table_type)
                    
                    print(f"{i:3d}. {type_icon} {table_name:<35} ({compound_count:>6,} compounds) [{table_type}]")
                    table_info.append({
                        'name': table_name,
                        'count': compound_count,
                        'type': table_type
                    })
                    
                except Exception as e:
                    print(f"{i:3d}. ❓ {table_name:<35} (Error: {e}) [Unknown]")
                    table_info.append({
                        'name': table_name,
                        'count': 0,
                        'type': 'unknown'
                    })
            
            print("=" * 50)
            print("Commands: Enter table number, table name, or 'cancel' to abort")
            
            # Interactive table selection
            while True:
                try:
                    selection = input(f"\n🔍 Select table for Ersilia processing: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Ersilia model application cancelled")
                        return None
                    
                    # Try as number first
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(available_tables):
                            selected_table = available_tables[table_idx]
                            selected_info = table_info[table_idx]
                            break
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(available_tables)}")
                            continue
                    except ValueError:
                        # Try as table name
                        matching_tables = [t for t in available_tables if t.lower() == selection.lower()]
                        if matching_tables:
                            selected_table = matching_tables[0]
                            # Find table info
                            selected_info = next((info for info in table_info if info['name'] == selected_table), 
                                            {'count': 0, 'type': 'unknown'})
                            break
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Table selection cancelled")
                    return None
            
            print(f"\n✅ Selected table: '{selected_table}'")
            print(f"   🏷️  Type: {selected_info['type']}")
            print(f"   📊 Compounds: {selected_info['count']:,}")
            
            # Prompt for Ersilia model name
            print(f"\n🤖 ERSILIA MODEL SELECTION")
            print("-" * 30)
            print("Enter the Ersilia model identifier to use for predictions.")
            print("Examples:")
            print("   • eos3b5e (Molecular weight)")
            print("   • eos4e40 (Broad spectrum antibiotic activity)")
            print("   • eos7a45 (Molecule price prediction)")
            print("   • eos1amr (Blood-brain barrier penetration)")
            print("-" * 30)
            
            while True:
                try:
                    model_name = input("🧬 Enter Ersilia model name (or 'cancel' to abort): ").strip()
                    
                    if model_name.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Ersilia model application cancelled")
                        return None
                    
                    if not model_name:
                        print("❌ Model name cannot be empty. Please enter a valid model identifier.")
                        continue
                    
                    # Validate model name format (Ersilia models typically start with 'eos')
                    if not model_name.startswith('eos') and len(model_name) < 6:
                        print("⚠️  Warning: Model name doesn't follow typical Ersilia format (e.g., 'eos3b5e')")
                        proceed = input("   Continue anyway? (y/n): ").strip().lower()
                        if proceed not in ['y', 'yes']:
                            continue
                    
                    print(f"✅ Selected Ersilia model: '{model_name}'")
                    break
                    
                except KeyboardInterrupt:
                    print("\n❌ Model selection cancelled")
                    return None
            
            # Get compounds from selected table
            print(f"\n📊 Retrieving SMILES from table '{selected_table}'...")
            compounds_df = self._get_table_as_dataframe(selected_table)
            
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{selected_table}'")
                return None
            
            # Check if SMILES column exists
            if 'smiles' not in compounds_df.columns:
                print(f"❌ No 'smiles' column found in table '{selected_table}'")
                print(f"   Available columns: {list(compounds_df.columns)}")
                return None
            
            # Extract SMILES strings and filter out invalid ones
            print(f"🧪 Extracting and validating SMILES strings...")
            smiles_list = []
            invalid_count = 0
            
            for _, row in compounds_df.iterrows():
                smiles = row['smiles']
                
                # Skip invalid/empty SMILES
                if pd.isna(smiles) or str(smiles).strip() == '' or str(smiles).lower() in ['nan', 'none']:
                    invalid_count += 1
                    continue
                
                smiles_clean = str(smiles).strip()
                if smiles_clean:
                    smiles_list.append(smiles_clean)
                else:
                    invalid_count += 1
            
            # Display extraction results
            print(f"\n📊 SMILES EXTRACTION RESULTS:")
            print(f"   ✅ Valid SMILES: {len(smiles_list):,}")
            print(f"   ❌ Invalid/empty SMILES: {invalid_count:,}")
            print(f"   📈 Success rate: {(len(smiles_list)/(len(smiles_list)+invalid_count))*100:.1f}%")
            
            if not smiles_list:
                print("❌ No valid SMILES found in the selected table")
                return None
            
            # Show sample of SMILES
            print(f"\n🔬 Sample SMILES (first 5):")
            for i, smiles in enumerate(smiles_list[:5], 1):
                display_smiles = smiles[:50] + "..." if len(smiles) > 50 else smiles
                print(f"   {i}. {display_smiles}")
            
            if len(smiles_list) > 5:
                print(f"   ... and {len(smiles_list) - 5} more SMILES")
            
            
            # Process with Ersilia and display results
            results_df = self._process_smiles_for_ersilia(smiles_list, model_name)
                
            if not results_df.empty:
                print(f"\n🎯 ERSILIA MODEL RESULTS:")
                print(f"   📊 Results shape: {results_df.shape}")
                print(f"   📋 Columns: {list(results_df.columns)}")
                    
                # Show sample results
                print(f"\n📋 Sample results (first 5 rows):")
                print(results_df.head())
                    
                # Ask if user wants to merge results back to original table
                print(f"\n🔗 MERGE OPTIONS:")
                print("1. Merge results to new table (recommended)")
                print("2. Update original table with predictions")
                print("3. Save results as separate CSV file")
                print("4. Just display results (no saving)")
                    
                merge_choice = input("Select option (1-4): ").strip()
                    
                if merge_choice == '1':
                    # Merge to new table
                    merged_table = self._merge_ersilia_results_to_table(
                        selected_table, results_df, model_name, save_to_new_table=True
                    )
                    if merged_table:
                        print(f"🎉 Merged results saved to new table: '{merged_table}'")
                    
                elif merge_choice == '2':
                    # Update original table
                    confirm = input(f"⚠️  This will modify the original table '{selected_table}'. Continue? (y/n): ").strip().lower()
                    if confirm in ['y', 'yes']:
                        merged_table = self._merge_ersilia_results_to_table(
                            selected_table, results_df, model_name, save_to_new_table=False
                        )
                        if merged_table:
                            print(f"🎉 Original table '{selected_table}' updated with Ersilia predictions")
                    
                elif merge_choice == '3':
                    # Save as CSV
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    csv_filename = f"ersilia_{model_name}_{selected_table}_{timestamp}.csv"
                    csv_path = os.path.join(self.path, 'chemspace', 'ersilia_results', csv_filename)
                    
                    # Ensure directory exists
                    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
                    
                    results_df.to_csv(csv_path, index=False)
                    print(f"📁 Results saved to CSV: {csv_path}")
                
                elif merge_choice == '4':
                    print("📊 Results displayed only (not saved)")
                
                else:
                    print("❌ Invalid option selected")
                        
        except Exception as e:
            print(f"❌ Error in apply_ersilia_model: {e}")
            return None    
        
    def _process_smiles_for_ersilia(self, smiles_list: List[str], model_name: str) -> pd.DataFrame:
        """
        Process a SMILES list with a given Ersilia model.
        
        Args:
            smiles_list (List[str]): List of SMILES strings to process
            model_name (str): Name of the Ersilia model to use
            
        Returns:
            pd.DataFrame: Results from the Ersilia model with SMILES column preserved
        """
        try:
            # Import Ersilia Model here since it's needed in this method
            from ersilia.api import Model
            
            print(f"\n🧬 Processing {len(smiles_list)} SMILES with Ersilia model '{model_name}'...")
            
            # Initialize and set up the model
            print(f"   📥 Fetching model '{model_name}'...")
            mdl = Model(model_name)
            
            print(f"   🔄 Loading model...")
            mdl.fetch()
            
            print(f"   🚀 Starting model server...")
            mdl.serve()
            
            print(f"   🧪 Running predictions...")
            results_df = mdl.run(smiles_list)
            
            print(f"   🛑 Closing model server...")
            mdl.close()
            
            print(f"   ✅ Model processing completed!")
            print(f"   📊 Results shape: {results_df.shape}")
            
            # Add SMILES column to results if it's missing
            if 'smiles' not in results_df.columns:
                print(f"   🔗 Adding SMILES column to results...")
                
                # Ensure we have the same number of SMILES as results
                if len(results_df) == len(smiles_list):
                    results_df.insert(0, 'smiles', smiles_list)
                    print(f"   ✅ SMILES column added successfully")
                else:
                    print(f"   ⚠️  Warning: Mismatch between input SMILES ({len(smiles_list)}) and results ({len(results_df)})")
                    # Try to match by creating a truncated or padded SMILES list
                    if len(results_df) < len(smiles_list):
                        # Results are fewer than input - take first N SMILES
                        matched_smiles = smiles_list[:len(results_df)]
                        results_df.insert(0, 'smiles', matched_smiles)
                        print(f"   🔗 Added SMILES for first {len(results_df)} results")
                    else:
                        # Results are more than input - this shouldn't happen, but handle it
                        extended_smiles = smiles_list + ['unknown'] * (len(results_df) - len(smiles_list))
                        results_df.insert(0, 'smiles', extended_smiles)
                        print(f"   🔗 Extended SMILES list to match {len(results_df)} results")
            else:
                print(f"   ✅ SMILES column already present in results")
            
            # Verify final structure
            print(f"   📋 Final columns: {list(results_df.columns)}")
            
            return results_df
            
        except ImportError:
            print("❌ Ersilia not installed or not available")
            return pd.DataFrame()
        except Exception as e:
            print(f"❌ Error processing SMILES with Ersilia model: {e}")
            # Return empty DataFrame with at least the smiles column for error handling
            return pd.DataFrame({'smiles': []})
    
    def _merge_ersilia_results_to_table(self, selected_table: str, results_df: pd.DataFrame, 
                                    model_name: str, save_to_new_table: bool = True) -> Optional[str]:
        """
        Merge Ersilia model results back to the original table by matching SMILES.
        
        Args:
            selected_table (str): Name of the original table
            results_df (pd.DataFrame): Results from Ersilia model
            model_name (str): Name of the Ersilia model used
            save_to_new_table (bool): Whether to save to a new table or update existing
            
        Returns:
            Optional[str]: Name of the table with merged results, or None if failed
        """
        try:
            if results_df.empty:
                print("❌ No Ersilia results to merge")
                return None
            
            print(f"\n🔄 MERGING ERSILIA RESULTS TO ORIGINAL TABLE")
            print("=" * 60)
            
            # Get the original table data
            original_df = self._get_table_as_dataframe(selected_table)
            
            if original_df.empty:
                print(f"❌ Original table '{selected_table}' is empty")
                return None
            
            print(f"📊 Original table: {len(original_df)} compounds")
            print(f"📊 Ersilia results: {len(results_df)} predictions")
            
            # Check if both have SMILES column
            if 'smiles' not in original_df.columns:
                print("❌ Original table missing 'smiles' column")
                return None
            
            if 'smiles' not in results_df.columns:
                print("❌ Ersilia results missing 'smiles' column")
                return None
            
            # Get prediction columns (exclude specified columns and 'smiles')
            excluded_columns = {'smiles', 'key', 'input', 'ersilia_model', 'prediction_date'}
            prediction_columns = [col for col in results_df.columns if col not in excluded_columns]
            
            if not prediction_columns:
                print("❌ No prediction columns found in Ersilia results after filtering")
                print(f"   Available columns: {list(results_df.columns)}")
                print(f"   Excluded columns: {excluded_columns}")
                return None
            
            print(f"🧪 Prediction columns to add: {prediction_columns}")
            
            # Merge DataFrames on SMILES
            print(f"🔗 Merging tables on SMILES column...")
            
            # Use left join to keep all original compounds
            merged_df = original_df.merge(
                results_df, 
                on='smiles', 
                how='left', 
                suffixes=('', f'_{model_name}')
            )
            
            # Count successful matches
            successful_matches = 0
            missing_predictions = 0
            
            for col in prediction_columns:
                non_null_count = merged_df[col].notna().sum()
                successful_matches = max(successful_matches, non_null_count)
            
            missing_predictions = len(merged_df) - successful_matches
            
            print(f"📊 MERGE RESULTS:")
            print(f"   ✅ Successful matches: {successful_matches}")
            print(f"   ❌ Missing predictions: {missing_predictions}")
            print(f"   📈 Match rate: {(successful_matches/len(merged_df))*100:.1f}%")
            
            # Show sample of merged data
            print(f"\n📋 Sample merged data (first 3 rows):")
            print("=" * 80)
            
            # Display key columns
            display_cols = ['smiles', 'name'] + prediction_columns[:3]  # Show first 3 prediction cols
            available_display_cols = [col for col in display_cols if col in merged_df.columns]
            
            for i, (_, row) in enumerate(merged_df.head(3).iterrows(), 1):
                print(f"Row {i}:")
                for col in available_display_cols:
                    value = row[col]
                    if pd.isna(value):
                        value_str = "NaN"
                    else:
                        value_str = str(value)[:30] + "..." if len(str(value)) > 30 else str(value)
                    print(f"   {col}: {value_str}")
                print()
            
            print("=" * 80)
            
            # Determine output table name
            if save_to_new_table:
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                default_name = f"{selected_table}_{model_name}_{timestamp}"
                
                user_input = input(f"\n💾 Enter table name for merged results (default: {default_name}): ").strip()
                output_table_name = user_input if user_input else default_name
                output_table_name = self._sanitize_table_name(output_table_name)
            else:
                output_table_name = selected_table
                print(f"\n💾 Updating existing table '{selected_table}' with Ersilia results...")
            
            # Save merged results to database
            success = self._save_merged_ersilia_results(merged_df, output_table_name, model_name, prediction_columns)
            
            if success:
                print(f"✅ Successfully saved merged results to table '{output_table_name}'")
                print(f"   📊 Total compounds: {len(merged_df)}")
                print(f"   🧪 Prediction columns: {len(prediction_columns)}")
                print(f"   🎯 Model used: {model_name}")
                return output_table_name
            else:
                print(f"❌ Failed to save merged results")
                return None
            
        except Exception as e:
            print(f"❌ Error merging Ersilia results: {e}")
            return None

    def _save_merged_ersilia_results(self, merged_df: pd.DataFrame, table_name: str, 
                                    model_name: str, prediction_columns: List[str]) -> bool:
        """
        Save merged Ersilia results to database table.
        
        Args:
            merged_df (pd.DataFrame): Merged DataFrame with original data and predictions
            table_name (str): Name of the table to create/update
            model_name (str): Name of the Ersilia model used
            prediction_columns (List[str]): List of prediction column names
            
        Returns:
            bool: True if saved successfully
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create table with extended schema including prediction columns (no metadata columns)
            base_schema = """
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT NOT NULL,
                name TEXT,
                flag TEXT,
                inchi_key TEXT
            """
            
            # Add prediction columns to schema with model name prefix
            prediction_schema_parts = []
            prefixed_column_mapping = {}
            
            for col in prediction_columns:
                # Create prefixed column name: modelname_property
                prefixed_col_name = f"{model_name}_{col}"
                
                # Sanitize column name for SQL
                safe_col_name = self._sanitize_column_name(prefixed_col_name)
                prediction_schema_parts.append(f"{safe_col_name} REAL")
                
                # Keep mapping of original column to prefixed column for data insertion
                prefixed_column_mapping[col] = safe_col_name
            
            # Add unique constraint (no metadata columns)
            unique_constraint = "UNIQUE(smiles, name)"
            
            full_schema = base_schema
            if prediction_schema_parts:
                full_schema += ",\n            " + ",\n            ".join(prediction_schema_parts)
            full_schema += ",\n            " + unique_constraint
            
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {table_name} (
                    {full_schema}
                )
            ''')
            
            # Create indexes for better performance (no model index)
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_smiles ON {table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_name ON {table_name}(name)")
            
            # Prepare insert query with prefixed column names (no metadata columns)
            columns = ['smiles', 'name', 'flag', 'inchi_key'] + list(prefixed_column_mapping.values())
            placeholders = ', '.join(['?'] * len(columns))
            
            insert_query = f'''
                INSERT OR REPLACE INTO {table_name} 
                ({', '.join(columns)})
                VALUES ({placeholders})
            '''
            
            # Insert data (no metadata)
            inserted_count = 0
            
            for _, row in merged_df.iterrows():
                try:
                    # Prepare row data (no metadata columns)
                    row_data = [
                        row.get('smiles', ''),
                        row.get('name', 'unknown'),
                        row.get('flag', 'nd'),
                        row.get('inchi_key', None)
                    ]
                    
                    # Add prediction values only with original column names
                    for original_col in prediction_columns:
                        value = row.get(original_col)
                        # Convert to float if possible, otherwise None
                        if pd.isna(value):
                            row_data.append(None)
                        else:
                            try:
                                row_data.append(float(value))
                            except (ValueError, TypeError):
                                row_data.append(None)
                    
                    # No metadata added here
                    cursor.execute(insert_query, row_data)
                    inserted_count += 1
                    
                except Exception as e:
                    print(f"⚠️  Error inserting row: {e}")
                    continue
            
            conn.commit()
            conn.close()
            
            # Show the prefixed column names that were created
            print(f"   💾 Inserted {inserted_count} compounds with Ersilia predictions")
            print(f"   🏷️  Created columns with model prefix:")
            for original_col, prefixed_col in prefixed_column_mapping.items():
                print(f"      • {original_col} → {prefixed_col}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error saving merged Ersilia results: {e}")
            return False

    def _sanitize_column_name(self, column_name: str) -> str:
        """
        Sanitize column name for SQLite compatibility.
        
        Args:
            column_name (str): Original column name
            
        Returns:
            str: Sanitized column name
        """
        import re
        
        # Replace special characters with underscores
        sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', str(column_name))
        
        # Ensure it starts with a letter or underscore
        if sanitized and not sanitized[0].isalpha() and sanitized[0] != '_':
            sanitized = 'col_' + sanitized
        
        # Ensure it's not empty
        if not sanitized:
            sanitized = 'prediction_col'
        
        # Limit length
        if len(sanitized) > 50:
            sanitized = sanitized[:50]
        
        return sanitized.lower()

    def enumerate_stereoisomers(self, table_name: Optional[str] = None,
                            max_stereoisomers: int = 20,
                            save_results: bool = True,
                            result_table_name: Optional[str] = None,
                            include_original: bool = False,
                            filter_duplicates: bool = True) -> pd.DataFrame:
        """
        Enumerate stereoisomers for molecules in a table using RDKit.
        Computes stereoisomeric configuration and total count for each compound.
        
        Args:
            table_name (Optional[str]): Name of the table to process. If None, prompts user.
            max_stereoisomers (int): Maximum number of stereoisomers to generate per molecule
            save_results (bool): Whether to save results to a new table
            result_table_name (Optional[str]): Name for the result table
            include_original (bool): Whether to include the original molecule in results
            filter_duplicates (bool): Whether to filter out duplicate stereoisomers
            
        Returns:
            pd.DataFrame: DataFrame containing original molecules and their stereoisomers
        """
        try:
            # Import required libraries
            try:
                from rdkit import Chem
                from rdkit.Chem import EnumerateStereoisomers
                from rdkit import RDLogger
                RDLogger.DisableLog('rdApp.*')
            except ImportError:
                print("❌ RDKit not installed. Please install RDKit to enumerate stereoisomers:")
                print("   conda install -c conda-forge rdkit")
                return pd.DataFrame()
            
            # Reuse existing table selection method
            if table_name is None:
                table_name = self._select_table_interactive("SELECT TABLE FOR STEREOISOMER ENUMERATION")
                if not table_name:
                    print("❌ No table selected for stereoisomer enumeration")
                    return pd.DataFrame()
            
            # Get compounds from table using existing method
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return pd.DataFrame()
            
            # Check if SMILES column exists (reuse existing pattern)
            if 'smiles' not in compounds_df.columns:
                print("❌ No 'smiles' column found in the table")
                return pd.DataFrame()
            
            # Analyze SMILES for existing stereochemistry
            stereochemistry_analysis = self._analyze_existing_stereochemistry(compounds_df)
            
            # Show stereochemistry analysis and get user preferences
            enumeration_strategy = self._get_stereochemistry_enumeration_strategy(stereochemistry_analysis)
            
            if enumeration_strategy is None:
                print("❌ Stereoisomer enumeration cancelled by user")
                return pd.DataFrame()
            
            print(f"🧬 Enumerating stereoisomers for table '{table_name}'")
            print(f"   📊 Input molecules: {len(compounds_df):,}")
            print(f"   🔄 Max stereoisomers per molecule: {max_stereoisomers}")
            print(f"   📋 Include original molecules: {'Yes' if include_original else 'No'}")
            print(f"   🔍 Filter duplicates: {'Yes' if filter_duplicates else 'No'}")
            print(f"   ⚗️  Stereochemistry strategy: {enumeration_strategy['strategy_name']}")
            
            # Process molecules and enumerate stereoisomers
            print("\n🔬 Processing molecules for stereoisomer enumeration...")
            all_stereoisomers = []
            processing_stats = {
                'processed': 0,
                'with_stereocenters': 0,
                'total_stereoisomers': 0,
                'invalid_smiles': 0,
                'failed_enumeration': 0,
                'inchi_computed': 0,
                'inchi_failed': 0,
                'predefined_kept': 0,
                'predefined_enumerated': 0
            }
            
            # Use progress bar if available (reuse existing pattern)
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    compounds_df.iterrows(),
                    total=len(compounds_df),
                    desc="Enumerating stereoisomers",
                    unit="molecules"
                )
            else:
                progress_bar = compounds_df.iterrows()
                processed_count = 0
            
            for _, compound in progress_bar:
                processing_stats['processed'] += 1
                
                try:
                    original_smiles = compound['smiles']
                    compound_name = compound.get('name', f"compound_{compound.get('id', processing_stats['processed'])}")
                    
                    # Parse original molecule
                    mol = Chem.MolFromSmiles(original_smiles)
                    if mol is None:
                        processing_stats['invalid_smiles'] += 1
                        continue
                    
                    # Check if this molecule has predefined stereochemistry
                    has_predefined_stereo = self._has_stereochemistry_definition(original_smiles)
                    
                    # Decide whether to enumerate based on user strategy
                    should_enumerate = self._should_enumerate_compound(
                        has_predefined_stereo, enumeration_strategy
                    )
                    
                    if has_predefined_stereo and not should_enumerate:
                        # Keep the predefined stereochemistry as-is
                        processing_stats['predefined_kept'] += 1
                        
                        # Analyze stereoisomeric properties for metadata
                        stereo_info = self._analyze_compound_stereochemistry(mol, original_smiles)
                        
                        # Compute InChI key
                        original_inchi_key = self._compute_single_inchi_key(original_smiles)
                        if original_inchi_key and original_inchi_key not in ['INVALID_SMILES', 'ERROR']:
                            processing_stats['inchi_computed'] += 1
                        else:
                            processing_stats['inchi_failed'] += 1
                        
                        # Add the predefined stereoisomer to results
                        all_stereoisomers.append({
                            'parent_name': compound_name,
                            'parent_smiles': original_smiles,
                            'stereoisomer_smiles': original_smiles,  # ← Keep original SMILES
                            'stereoisomer_name': compound_name,  # ← Keep original name
                            'stereoisomer_type': 'predefined',  # ← Correct type
                            'stereoisomer_index': 0,  # ← Index 0 for predefined
                            'flag': compound.get('flag', 'nd'),
                            'inchi_key': original_inchi_key,  # ← Use the computed InChI key
                            'total_stereoisomers': stereo_info['total_stereoisomers'],
                            'stereochemical_config': stereo_info['original_config'],  # ← Use original config
                            'num_stereocenters': stereo_info['num_stereocenters'],
                            'stereocenter_atoms': ','.join(map(str, stereo_info['stereocenter_atoms']))
                        })
                        
                        
                        processing_stats['total_stereoisomers'] += 1
                        continue
                    
                    # If we reach here, we should enumerate stereoisomers
                    if has_predefined_stereo:
                        processing_stats['predefined_enumerated'] += 1
                    
                    # Analyze stereoisomeric properties using helper method
                    stereo_info = self._analyze_compound_stereochemistry(mol, original_smiles)
                    
                    # Include original molecule if requested
                    if include_original:
                        # Compute InChI key for original molecule
                        original_inchi_key = self._compute_single_inchi_key(original_smiles)
                        if original_inchi_key and original_inchi_key not in ['INVALID_SMILES', 'ERROR']:
                            processing_stats['inchi_computed'] += 1
                        else:
                            processing_stats['inchi_failed'] += 1
                        
                        all_stereoisomers.append({
                            'parent_name': compound_name,
                            'parent_smiles': original_smiles,
                            'stereoisomer_smiles': original_smiles,
                            'stereoisomer_name': f"{compound_name}_original",
                            'stereoisomer_type': 'original',
                            'stereoisomer_index': 0,
                            'flag': compound.get('flag', 'nd'),
                            'inchi_key': original_inchi_key,
                            'total_stereoisomers': stereo_info['total_stereoisomers'],
                            'stereochemical_config': stereo_info['original_config'],
                            'num_stereocenters': stereo_info['num_stereocenters'],
                            'stereocenter_atoms': ','.join(map(str, stereo_info['stereocenter_atoms']))
                        })
                        processing_stats['total_stereoisomers'] += 1
                    
                    # Check for stereocenters and enumerate
                    if stereo_info['total_stereoisomers'] > 1:
                        processing_stats['with_stereocenters'] += 1
                        
                        # Enumerate stereoisomers
                        stereoisomer_options = EnumerateStereoisomers.EnumerateStereoisomers(mol)
                        generated_smiles = set()
                        stereoisomer_count = 0
                        
                        for i, stereoisomer in enumerate(stereoisomer_options):
                            if stereoisomer_count >= max_stereoisomers:
                                break
                            
                            try:
                                # Generate SMILES for stereoisomer
                                stereo_smiles = Chem.MolToSmiles(stereoisomer, isomericSmiles=True)
                                
                                # Skip if duplicate and filtering is enabled
                                if filter_duplicates:
                                    if stereo_smiles in generated_smiles:
                                        continue
                                    generated_smiles.add(stereo_smiles)
                                
                                # Skip original molecule if we already included it
                                if include_original and stereo_smiles == original_smiles:
                                    continue
                                
                                # Compute InChI key for stereoisomer using existing helper
                                stereo_inchi_key = self._compute_single_inchi_key(stereo_smiles)
                                if stereo_inchi_key and stereo_inchi_key not in ['INVALID_SMILES', 'ERROR']:
                                    processing_stats['inchi_computed'] += 1
                                else:
                                    processing_stats['inchi_failed'] += 1
                                
                                # Analyze stereochemical configuration for this stereoisomer
                                stereo_config = self._determine_stereoisomer_configuration(stereoisomer)
                                
                                # Add stereoisomer to results
                                stereoisomer_name = f"{compound_name}_stereo_{stereoisomer_count + 1}"
                                
                                all_stereoisomers.append({
                                    'parent_name': compound_name,
                                    'parent_smiles': original_smiles,
                                    'stereoisomer_smiles': stereo_smiles,
                                    'stereoisomer_name': stereoisomer_name,
                                    'stereoisomer_type': 'stereoisomer',
                                    'stereoisomer_index': stereoisomer_count + 1,
                                    'flag': compound.get('flag', 'nd'),  # ← Preserve original flag
                                    'inchi_key': stereo_inchi_key,
                                    'total_stereoisomers': stereo_info['total_stereoisomers'],
                                    'stereochemical_config': stereo_config,
                                    'num_stereocenters': stereo_info['num_stereocenters'],
                                    'stereocenter_atoms': ','.join(map(str, stereo_info['stereocenter_atoms']))
                                })
                                
                                stereoisomer_count += 1
                                processing_stats['total_stereoisomers'] += 1
                                
                            except Exception as e:
                                processing_stats['failed_enumeration'] += 1
                                continue
                
                except Exception as e:
                    processing_stats['failed_enumeration'] += 1
                    continue
                
                # Progress update for non-tqdm case (reuse existing pattern)
                if not TQDM_AVAILABLE:
                    processed_count += 1
                    if processed_count % max(1, len(compounds_df) // 10) == 0:
                        print(f"   📊 Processed {processed_count}/{len(compounds_df)} molecules...")
            
            if TQDM_AVAILABLE and hasattr(progress_bar, 'close'):
                progress_bar.close()
            
            # Create results DataFrame
            results_df = pd.DataFrame(all_stereoisomers)
            
            # Display processing statistics (enhanced with stereochemistry info)
            print(f"\n✅ Stereoisomer enumeration completed!")
            print(f"   📊 Molecules processed: {processing_stats['processed']:,}")
            print(f"   🎯 Molecules with stereocenters: {processing_stats['with_stereocenters']:,}")
            print(f"   🧪 Total stereoisomers generated: {processing_stats['total_stereoisomers']:,}")
            print(f"   ❌ Invalid SMILES: {processing_stats['invalid_smiles']:,}")
            print(f"   ⚠️  Failed enumerations: {processing_stats['failed_enumeration']:,}")
            print(f"   🔬 InChI keys computed: {processing_stats['inchi_computed']:,}")
            print(f"   ❌ InChI computation failures: {processing_stats['inchi_failed']:,}")
            
            # Show stereochemistry handling statistics
            if processing_stats['predefined_kept'] > 0 or processing_stats['predefined_enumerated'] > 0:
                print(f"\n⚗️  Stereochemistry Handling:")
                print(f"   ✅ Predefined stereochemistry kept: {processing_stats['predefined_kept']:,}")
                print(f"   🔄 Predefined stereochemistry enumerated: {processing_stats['predefined_enumerated']:,}")
            
            if not results_df.empty:
                expansion_factor = len(results_df) / len(compounds_df)
                inchi_success_rate = (processing_stats['inchi_computed'] / processing_stats['total_stereoisomers']) * 100 if processing_stats['total_stereoisomers'] > 0 else 0
                
                print(f"   📈 Expansion factor: {expansion_factor:.2f}x")
                print(f"   🔑 InChI success rate: {inchi_success_rate:.1f}%")
                
                # Show sample results with stereochemical information
                print(f"\n📋 Sample stereoisomers with stereochemical information (first 5):")
                for i, (_, row) in enumerate(results_df.head(5).iterrows(), 1):
                    if row['stereoisomer_type'] == 'original':
                        stereo_type = "Original"
                    elif row['stereoisomer_type'] == 'predefined':
                        stereo_type = "Predefined"
                    else:
                        stereo_type = f"Stereo {row['stereoisomer_index']}"
                    
                    inchi_status = "✅" if row.get('inchi_key') and row['inchi_key'] not in ['INVALID_SMILES', 'ERROR'] else "❌"
                    
                    # Clean the configuration for display
                    clean_config = self._clean_stereochemical_config(row['stereochemical_config'])
                    
                    print(f"   {i}. {row['parent_name']} → {stereo_type} {inchi_status}")
                    print(f"      SMILES: {row['stereoisomer_smiles'][:40]}{'...' if len(row['stereoisomer_smiles']) > 40 else ''}")
                    print(f"      Total stereoisomers: {row['total_stereoisomers']}")
                    print(f"      Stereocenters: {row['num_stereocenters']}")
                    print(f"      Configuration: {clean_config}")
                    if row.get('inchi_key') and row['inchi_key'] not in ['INVALID_SMILES', 'ERROR']:
                        print(f"      InChI Key: {row['inchi_key'][:27]}...")
            
            # Save results if requested (reuse existing pattern)
            if save_results and not results_df.empty:
                if result_table_name is None:
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    default_name = f"{table_name}_stereoisomers_{timestamp}"
                    user_table_name = input(f"\n💾 Enter table name for results (default: {default_name}): ").strip()
                    result_table_name = user_table_name if user_table_name else default_name
                
                result_table_name = self._sanitize_table_name(result_table_name)
                
                print(f"   💾 Saving stereoisomers to table '{result_table_name}'...")
                success = self._save_stereoisomer_results(results_df, result_table_name, table_name, processing_stats)
                
                if success:
                    print(f"   ✅ Results saved to table: '{result_table_name}'")
                else:
                    print(f"   ❌ Failed to save results to database")
            
            return results_df
            
        except Exception as e:
            print(f"❌ Error enumerating stereoisomers: {e}")
            return pd.DataFrame()

    def _analyze_existing_stereochemistry(self, compounds_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Analyze the compounds to detect existing stereochemistry definitions.
        
        Args:
            compounds_df (pd.DataFrame): DataFrame containing compounds
            
        Returns:
            Dict[str, Any]: Analysis results about existing stereochemistry
        """
        try:
            total_compounds = len(compounds_df)
            compounds_with_stereo = 0
            compounds_without_stereo = 0
            sample_with_stereo = []
            sample_without_stereo = []
            
            for _, row in compounds_df.iterrows():
                smiles = row.get('smiles', '')
                name = row.get('name', 'unknown')
                
                has_stereo = self._has_stereochemistry_definition(smiles)
                
                if has_stereo:
                    compounds_with_stereo += 1
                    if len(sample_with_stereo) < 5:
                        sample_with_stereo.append({
                            'name': name,
                            'smiles': smiles[:50] + ('...' if len(smiles) > 50 else '')
                        })
                else:
                    compounds_without_stereo += 1
                    if len(sample_without_stereo) < 5:
                        sample_without_stereo.append({
                            'name': name,
                            'smiles': smiles[:50] + ('...' if len(smiles) > 50 else '')
                        })
            
            return {
                'total_compounds': total_compounds,
                'with_stereochemistry': compounds_with_stereo,
                'without_stereochemistry': compounds_without_stereo,
                'sample_with_stereo': sample_with_stereo,
                'sample_without_stereo': sample_without_stereo,
                'percentage_with_stereo': (compounds_with_stereo / total_compounds * 100) if total_compounds > 0 else 0
            }
            
        except Exception as e:
            print(f"❌ Error analyzing existing stereochemistry: {e}")
            return {
                'total_compounds': 0,
                'with_stereochemistry': 0,
                'without_stereochemistry': 0,
                'sample_with_stereo': [],
                'sample_without_stereo': [],
                'percentage_with_stereo': 0
            }

    def _has_stereochemistry_definition(self, smiles: str) -> bool:
        """
        Check if a SMILES string contains stereochemistry definitions.
        
        Args:
            smiles (str): SMILES string to check
            
        Returns:
            bool: True if stereochemistry is defined, False otherwise
        """
        try:
            if not smiles or pd.isna(smiles):
                return False
            
            smiles_str = str(smiles)
            
            # Check for common stereochemistry markers
            stereo_markers = ['@', '/', '\\']
            
            return any(marker in smiles_str for marker in stereo_markers)
            
        except Exception:
            return False

    def _get_stereochemistry_enumeration_strategy(self, analysis: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Get user preferences for handling existing stereochemistry during enumeration.
        
        Args:
            analysis (Dict): Results from stereochemistry analysis
            
        Returns:
            Optional[Dict]: Strategy configuration or None if cancelled
        """
        try:
            print(f"\n⚗️  EXISTING STEREOCHEMISTRY ANALYSIS")
            print("=" * 60)
            print(f"📊 Total compounds: {analysis['total_compounds']:,}")
            print(f"🧬 With stereochemistry: {analysis['with_stereochemistry']:,} ({analysis['percentage_with_stereo']:.1f}%)")
            print(f"🔄 Without stereochemistry: {analysis['without_stereochemistry']:,}")
            
            if analysis['with_stereochemistry'] > 0:
                print(f"\n📋 Sample compounds with stereochemistry (first 5):")
                for i, sample in enumerate(analysis['sample_with_stereo'], 1):
                    print(f"   {i}. {sample['name']}: {sample['smiles']}")
            
            if analysis['without_stereochemistry'] > 0:
                print(f"\n📋 Sample compounds without stereochemistry (first 5):")
                for i, sample in enumerate(analysis['sample_without_stereo'], 1):
                    print(f"   {i}. {sample['name']}: {sample['smiles']}")
            
            print("=" * 60)
            
            if analysis['with_stereochemistry'] == 0:
                # No stereochemistry found, proceed normally
                return {
                    'strategy': 'enumerate_all',
                    'strategy_name': 'Enumerate all (no predefined stereochemistry found)'
                }
            
            # Ask user what to do with predefined stereochemistry
            print(f"\n❓ STEREOCHEMISTRY HANDLING OPTIONS")
            print("-" * 40)
            print("How should compounds with predefined stereochemistry be handled?")
            print()
            print("1. Keep predefined stereochemistry as-is (recommended)")
            print("   → Compounds with @ / \\ symbols will be preserved")
            print("   → Only compounds without stereochemistry will be enumerated")
            print()
            print("2. Enumerate all compounds (ignore predefined stereochemistry)")
            print("   → All compounds will be enumerated regardless of existing stereochemistry")
            print("   → May generate duplicates of already defined stereoisomers")
            print()
            print("3. Only enumerate compounds without stereochemistry")
            print("   → Compounds with predefined stereochemistry will be kept as-is")
            print("   → Compounds without stereochemistry will be enumerated")
            print()
            print("4. Cancel enumeration")
            print("-" * 40)
            
            while True:
                try:
                    choice = input("Select option (1-4): ").strip()
                    
                    if choice == '1':
                        return {
                            'strategy': 'keep_predefined',
                            'strategy_name': 'Keep predefined stereochemistry, enumerate undefined',
                            'enumerate_predefined': False,
                            'enumerate_undefined': True
                        }
                    elif choice == '2':
                        print("⚠️  Warning: This will enumerate all compounds, potentially creating duplicates")
                        print("   of stereoisomers that are already defined in the SMILES.")
                        confirm = input("Continue anyway? (y/n): ").strip().lower()
                        if confirm in ['y', 'yes']:
                            return {
                                'strategy': 'enumerate_all',
                                'strategy_name': 'Enumerate all compounds (ignore predefined stereochemistry)',
                                'enumerate_predefined': True,
                                'enumerate_undefined': True
                            }
                        else:
                            continue
                    elif choice == '3':
                        return {
                            'strategy': 'enumerate_undefined_only',
                            'strategy_name': 'Keep predefined, enumerate only undefined stereochemistry',
                            'enumerate_predefined': False,
                            'enumerate_undefined': True
                        }
                    elif choice == '4':
                        return None
                    else:
                        print("❌ Invalid choice. Please enter 1, 2, 3, or 4.")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ Enumeration cancelled by user")
                    return None
        
        except Exception as e:
            print(f"❌ Error getting stereochemistry strategy: {e}")
            return None

    def _should_enumerate_compound(self, has_predefined_stereo: bool, strategy: Dict[str, Any]) -> bool:
        """
        Determine whether a compound should be enumerated based on its stereochemistry status and user strategy.
        
        Args:
            has_predefined_stereo (bool): Whether the compound has predefined stereochemistry
            strategy (Dict): User's enumeration strategy
            
        Returns:
            bool: True if the compound should be enumerated, False if it should be kept as-is
        """
        try:
            strategy_type = strategy.get('strategy', 'enumerate_all')
            
            if strategy_type == 'enumerate_all':
                return True
            elif strategy_type == 'keep_predefined':
                # Enumerate compounds without predefined stereochemistry
                return not has_predefined_stereo
            elif strategy_type == 'enumerate_undefined_only':
                # Only enumerate compounds without predefined stereochemistry
                return not has_predefined_stereo
            else:
                # Default: enumerate everything
                return True
                
        except Exception:
            return True

    def _analyze_compound_stereochemistry(self, mol, smiles: str) -> Dict[str, Any]:
        """
        Analyze stereochemical properties of a compound.
        
        Args:
            mol: RDKit molecule object
            smiles (str): SMILES string
            
        Returns:
            Dict[str, Any]: Dictionary containing stereochemical analysis
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import EnumerateStereoisomers
            
            # Find stereocenters
            stereo_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            num_stereocenters = len(stereo_info)
            stereocenter_atoms = [info[0] for info in stereo_info]
            
            # Count total possible stereoisomers
            stereoisomer_options = list(EnumerateStereoisomers.EnumerateStereoisomers(mol))
            total_stereoisomers = len(stereoisomer_options)
            
            # Determine original configuration
            original_config = self._determine_stereoisomer_configuration(mol)
            
            return {
                'num_stereocenters': num_stereocenters,
                'stereocenter_atoms': stereocenter_atoms,
                'total_stereoisomers': total_stereoisomers,
                'original_config': original_config
            }
            
        except Exception as e:
            return {
                'num_stereocenters': 0,
                'stereocenter_atoms': [],
                'total_stereoisomers': 1,
                'original_config': 'unknown'
            }

    def _determine_stereoisomer_configuration(self, mol) -> str:
        """
        Determine the stereochemical configuration of a molecule in simplified R/S format.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            str: String describing the stereochemical configuration (e.g., 'R,S', 'S,S')
        """
        try:
            from rdkit import Chem
            
            # Find chiral centers with their configurations
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            if not chiral_centers:
                return 'achiral'
            
            # Extract just the R/S designations
            rs_designations = []
            for atom_idx, chirality in chiral_centers:
                if chirality in ['R', 'S']:
                    rs_designations.append(chirality)
                else:
                    rs_designations.append('?')
            
            if not rs_designations:
                return 'unassigned'
            
            # Sort for consistent ordering (optional - you can remove this if you want to maintain atom order)
            rs_designations.sort()
            return ','.join(rs_designations)
            
        except Exception:
            return 'unknown'

    def _save_stereoisomer_results(self, results_df: pd.DataFrame, table_name: str, 
                                source_table: str, processing_stats: Dict[str, int]) -> bool:
        """
        Save stereoisomer enumeration results to database table with stereochemical information.
        Reuses existing table creation patterns for consistency.
        
        Args:
            results_df (pd.DataFrame): DataFrame containing stereoisomer results
            table_name (str): Name of the table to create
            source_table (str): Name of the source table
            processing_stats (Dict): Processing statistics
            
        Returns:
            bool: True if saved successfully
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Create stereoisomers table with extended schema (excluding stereocenter_atoms column)
            cursor.execute(f'''
                CREATE TABLE IF NOT EXISTS {table_name} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT NOT NULL,
                    name TEXT,
                    flag TEXT,
                    inchi_key TEXT,
                    total_stereoisomers INTEGER,
                    stereochemical_config TEXT,
                    num_stereocenters INTEGER,
                    UNIQUE(smiles, name)
                )
            ''')
            
            # Create indexes for better performance (excluding stereocenter_atoms)
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_smiles ON {table_name}(smiles)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_name ON {table_name}(name)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_flag ON {table_name}(flag)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_inchi_key ON {table_name}(inchi_key)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_stereocenters ON {table_name}(num_stereocenters)")
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_total_stereo ON {table_name}(total_stereoisomers)")
            
            # Prepare insert query with stereochemical columns (excluding stereocenter_atoms)
            insert_query = f'''
                INSERT OR IGNORE INTO {table_name} 
                (smiles, name, flag, inchi_key, total_stereoisomers, stereochemical_config, 
                num_stereocenters)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            '''
            
            # Insert stereoisomer data including stereochemical information
            inserted_count = 0
            
            for _, row in results_df.iterrows():
                try:
                    # Clean up the stereochemical configuration format
                    config = row.get('stereochemical_config', 'unknown')
                    clean_config = self._clean_stereochemical_config(config)
                    
                    cursor.execute(insert_query, (
                        row['stereoisomer_smiles'],
                        row['stereoisomer_name'],
                        row.get('flag', 'nd'),  # ← Use original flag or default to 'nd'
                        row.get('inchi_key'),
                        row.get('total_stereoisomers', 1),
                        clean_config,
                        row.get('num_stereocenters', 0)
                    ))
                    inserted_count += 1
                except sqlite3.IntegrityError:
                    continue  # Skip duplicates
            
            conn.commit()
            conn.close()
            
            # Display summary (reuse existing pattern)
            print(f"   💾 Inserted {inserted_count} stereoisomers into table '{table_name}'")
            print(f"   📊 Processing summary:")
            print(f"      • Molecules processed: {processing_stats['processed']}")
            print(f"      • With stereocenters: {processing_stats['with_stereocenters']}")
            print(f"      • Total stereoisomers: {processing_stats['total_stereoisomers']}")
            print(f"   🧬 Table columns:")
            print(f"      • total_stereoisomers: Total possible stereoisomers for parent")
            print(f"      • stereochemical_config: R/S configuration (e.g., 'R,S', 'S,S')")
            print(f"      • num_stereocenters: Number of stereogenic centers")
            
            return True
            
        except Exception as e:
            print(f"❌ Error saving stereoisomer results: {e}")
            return False

    def _clean_stereochemical_config(self, config: str) -> str:
        """
        Clean stereochemical configuration to simple R/S format.
        
        Args:
            config (str): Original configuration string (e.g., "5R,10S" or "achiral")
            
        Returns:
            str: Cleaned configuration (e.g., "R,S" or "achiral")
        """
        try:
            if not config or config in ['unknown', 'achiral', 'unassigned']:
                return config
            
            # Extract R/S designations from atom-specific format
            import re
            rs_matches = re.findall(r'[RS]', config)
            
            if rs_matches:
                # Sort to ensure consistent ordering (optional)
                rs_matches.sort()
                return ','.join(rs_matches)
            else:
                return config
                
        except Exception:
            return 'unknown'
    
    def _compute_single_inchi_key(self, smiles: str) -> Optional[str]:
        """
        Compute InChI key for a single SMILES string.
        Reuses existing RDKit patterns from the codebase.
        
        Args:
            smiles (str): SMILES string to process
            
        Returns:
            Optional[str]: InChI key or error indicator
        """
        try:
            from rdkit import Chem
            
            # Skip empty or None SMILES (reuse existing validation pattern)
            if not smiles or smiles.strip() == '':
                return 'INVALID_SMILES'
            
            # Parse SMILES with RDKit (reuse existing pattern)
            mol = Chem.MolFromSmiles(str(smiles).strip())
            
            if mol is not None:
                # Compute InChI key (reuse existing pattern)
                inchi_key = Chem.MolToInchiKey(mol)
                if inchi_key:  # Ensure InChI key is not empty
                    return inchi_key
                else:
                    return 'ERROR'
            else:
                # Invalid SMILES (reuse existing pattern)
                return 'INVALID_SMILES'
                
        except Exception:
            # Error computing InChI key (reuse existing pattern)
            return 'ERROR' 

    def enumerate_stereoisomers_for_specific_molecules(self, table_name: str, 
                                                    molecule_names: List[str],
                                                    max_stereoisomers: int = 20) -> pd.DataFrame:
        """
        Enumerate stereoisomers for specific molecules by name.
        Reuses the main enumeration method for consistency.
        
        Args:
            table_name (str): Name of the table containing molecules
            molecule_names (List[str]): List of molecule names to process
            max_stereoisomers (int): Maximum stereoisomers per molecule
            
        Returns:
            pd.DataFrame: DataFrame containing stereoisomers for specified molecules
        """
        try:
            print(f"🎯 Enumerating stereoisomers for specific molecules in table '{table_name}'")
            
            # Get the full table using existing method
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return pd.DataFrame()
            
            # Filter for specific molecules (reuse existing validation pattern)
            if 'name' not in compounds_df.columns:
                print("❌ No 'name' column found in table for molecule selection")
                return pd.DataFrame()
            
            filtered_df = compounds_df[compounds_df['name'].isin(molecule_names)]
            
            if filtered_df.empty:
                print(f"❌ None of the specified molecules found in table '{table_name}'")
                print(f"   Requested: {molecule_names}")
                return pd.DataFrame()
            
            found_molecules = filtered_df['name'].tolist()
            missing_molecules = [name for name in molecule_names if name not in found_molecules]
            
            print(f"   ✅ Found {len(found_molecules)} molecules")
            if missing_molecules:
                print(f"   ⚠️  Missing {len(missing_molecules)} molecules: {missing_molecules}")
            
            # Use the main enumeration method on the filtered data
            # Temporarily replace the original dataframe data for processing
            original_method = self._get_table_as_dataframe
            
            def mock_get_table(table):
                if table == table_name:
                    return filtered_df
                return original_method(table)
            
            # Temporarily replace the method
            self._get_table_as_dataframe = mock_get_table
            
            try:
                results_df = self.enumerate_stereoisomers(
                    table_name=table_name,
                    max_stereoisomers=max_stereoisomers,
                    save_results=False,
                    include_original=True,
                    filter_duplicates=True
                )
            finally:
                # Restore original method
                self._get_table_as_dataframe = original_method
            
            return results_df
            
        except Exception as e:
            print(f"❌ Error enumerating stereoisomers for specific molecules: {e}")
            return pd.DataFrame()

    def get_stereocenter_info(self, table_name: str) -> pd.DataFrame:
        """
        Analyze molecules in a table to identify stereocenters and stereoisomer potential.
        Reuses existing table access patterns for consistency.
        
        Args:
            table_name (str): Name of the table to analyze
            
        Returns:
            pd.DataFrame: DataFrame with stereocenter analysis for each molecule
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors, EnumerateStereoisomers
            
            print(f"🔬 Analyzing stereocenters in table '{table_name}'")
            
            # Use existing method to get table data
            compounds_df = self._get_table_as_dataframe(table_name)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{table_name}'")
                return pd.DataFrame()
            
            stereocenter_analysis = []
            
            # Use existing progress bar pattern
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    compounds_df.iterrows(),
                    total=len(compounds_df),
                    desc="Analyzing stereocenters",
                    unit="molecules"
                )
            else:
                progress_bar = compounds_df.iterrows()
            
            for _, compound in progress_bar:
                try:
                    smiles = compound['smiles']
                    name = compound.get('name', 'unknown')
                    
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        # Use consistent error handling pattern
                        stereocenter_analysis.append({
                            'name': name,
                            'smiles': smiles,
                            'valid_smiles': False,
                            'num_stereocenters': 0,
                            'num_undefined_stereocenters': 0,
                            'potential_stereoisomers': 0,
                            'has_stereochemistry': False,
                            'stereocenter_atoms': []
                        })
                        continue
                    
                    # Find stereocenters
                    stereo_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                    num_stereocenters = len(stereo_info)
                    
                    # Count undefined stereocenters
                    undefined_stereocenters = len([info for info in stereo_info if info[1] == '?'])
                    
                    # Estimate potential stereoisomers
                    stereoisomer_options = list(EnumerateStereoisomers.EnumerateStereoisomers(mol))
                    potential_stereoisomers = len(stereoisomer_options)
                    
                    # Check if molecule already has defined stereochemistry
                    has_stereochemistry = '@' in smiles
                    
                    # Get stereocenter atom indices
                    stereocenter_atoms = [info[0] for info in stereo_info]
                    
                    stereocenter_analysis.append({
                        'name': name,
                        'smiles': smiles,
                        'valid_smiles': True,
                        'num_stereocenters': num_stereocenters,
                        'num_undefined_stereocenters': undefined_stereocenters,
                        'potential_stereoisomers': potential_stereoisomers,
                        'has_stereochemistry': has_stereochemistry,
                        'stereocenter_atoms': stereocenter_atoms
                    })
                    
                except Exception:
                    # Use consistent error handling
                    stereocenter_analysis.append({
                        'name': name,
                        'smiles': smiles,
                        'valid_smiles': False,
                        'num_stereocenters': 0,
                        'num_undefined_stereocenters': 0,
                        'potential_stereoisomers': 0,
                        'has_stereochemistry': False,
                        'stereocenter_atoms': []
                    })
            
            if TQDM_AVAILABLE and hasattr(progress_bar, 'close'):
                progress_bar.close()
            
            analysis_df = pd.DataFrame(stereocenter_analysis)
            
            # Display summary statistics (reuse existing display patterns)
            if not analysis_df.empty:
                valid_molecules = analysis_df[analysis_df['valid_smiles']].shape[0]
                with_stereocenters = analysis_df[analysis_df['num_stereocenters'] > 0].shape[0]
                with_stereochemistry = analysis_df[analysis_df['has_stereochemistry']].shape[0]
                
                print(f"\n📊 Stereocenter Analysis Summary:")
                print(f"   📋 Total molecules: {len(analysis_df):,}")
                print(f"   ✅ Valid SMILES: {valid_molecules:,}")
                print(f"   🎯 With stereocenters: {with_stereocenters:,} ({(with_stereocenters/valid_molecules)*100:.1f}%)")
                print(f"   🧬 With defined stereochemistry: {with_stereochemistry:,} ({(with_stereochemistry/valid_molecules)*100:.1f}%)")
                
                if with_stereocenters > 0:
                    max_stereocenters = analysis_df['num_stereocenters'].max()
                    avg_stereocenters = analysis_df[analysis_df['num_stereocenters'] > 0]['num_stereocenters'].mean()
                    total_potential = analysis_df['potential_stereoisomers'].sum()
                    
                    print(f"   📈 Max stereocenters in one molecule: {max_stereocenters}")
                    print(f"   📊 Average stereocenters (for molecules with stereocenters): {avg_stereocenters:.1f}")
                    print(f"   🚀 Total potential stereoisomers: {total_potential:,}")
            
            return analysis_df
            
        except Exception as e:
            print(f"❌ Error analyzing stereocenters: {e}")
            return pd.DataFrame()

    def _select_table_interactive(self, prompt_title: str = "SELECT TABLE", 
                                filter_type: Optional[str] = None,
                                show_compound_count: bool = True) -> Optional[str]:
        """
        Unified interactive table selection method.
        
        Args:
            prompt_title (str): Title for the selection prompt
            filter_type (Optional[str]): Optional filter for table types
            show_compound_count (bool): Whether to show compound counts
            
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        try:
            available_tables = self.get_all_tables()
            
            if not available_tables:
                print("❌ No tables available for selection")
                return None
            
            # Filter tables if type specified
            if filter_type:
                filtered_tables = []
                for table in available_tables:
                    table_type = self._classify_table_type(table)
                    if filter_type.lower() in table_type.lower():
                        filtered_tables.append(table)
                available_tables = filtered_tables
                
                if not available_tables:
                    print(f"❌ No tables found matching filter: {filter_type}")
                    return None
            
            print(f"\n🔍 {prompt_title}")
            print("=" * 70)
            print(f"Available tables ({len(available_tables)} total):")
            print("-" * 70)
            
            # Display tables with optional compound counts
            table_info = []
            for i, table_name in enumerate(available_tables, 1):
                try:
                    if show_compound_count:
                        compound_count = self.get_compound_count(table_name=table_name)
                        table_type = self._classify_table_type(table_name)
                        type_icon = self._get_table_type_icon(table_type)
                        print(f"{i:3d}. {type_icon} {table_name:<30} ({compound_count:>6,} compounds) [{table_type}]")
                        table_info.append({
                            'name': table_name,
                            'type': table_type,
                            'count': compound_count
                        })
                    else:
                        print(f"{i:3d}. 📊 {table_name}")
                        table_info.append({
                            'name': table_name,
                            'type': 'unknown',
                            'count': 0
                        })
                except Exception as e:
                    print(f"{i:3d}. ❓ {table_name:<30} (Error: {e}) [Unknown]")
                    table_info.append({
                        'name': table_name,
                        'type': 'unknown',
                        'count': 0
                    })
            
            print("-" * 70)
            print("Commands: Enter table number, table name, or 'cancel' to abort")
            
            return self._get_user_table_selection(available_tables, table_info)
            
        except Exception as e:
            print(f"❌ Error in table selection: {e}")
            return None

    def _get_user_table_selection(self, available_tables: List[str], 
                                table_info: List[Dict]) -> Optional[str]:
        """
        Helper method for getting user table selection input.
        
        Args:
            available_tables (List[str]): List of available table names
            table_info (List[Dict]): List of table information dictionaries
            
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        while True:
            try:
                selection = input(f"\n🔍 Select table: ").strip()
                
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    return None
                
                # Try as number first
                try:
                    table_idx = int(selection) - 1
                    if 0 <= table_idx < len(available_tables):
                        selected_table = available_tables[table_idx]
                        selected_info = table_info[table_idx]
                        
                        print(f"\n✅ Selected table: '{selected_table}'")
                        if selected_info.get('type') != 'unknown':
                            print(f"   🏷️  Type: {selected_info['type']}")
                            print(f"   📊 Compounds: {selected_info['count']:,}")
                        
                        return selected_table
                    else:
                        print(f"❌ Invalid selection. Please enter 1-{len(available_tables)}")
                        continue
                except ValueError:
                    # Try as table name
                    matching_tables = [t for t in available_tables if t.lower() == selection.lower()]
                    if matching_tables:
                        selected_table = matching_tables[0]
                        # Find table info
                        selected_info = next((info for info in table_info if info['name'] == selected_table), 
                                        {'type': 'unknown', 'count': 0})
                        
                        print(f"\n✅ Selected table: '{selected_table}'")
                        if selected_info.get('type') != 'unknown':
                            print(f"   🏷️  Type: {selected_info['type']}")
                            print(f"   📊 Compounds: {selected_info['count']:,}")
                        
                        return selected_table
                    else:
                        print(f"❌ Table '{selection}' not found")
                        continue
                        
            except KeyboardInterrupt:
                print("\n❌ Table selection cancelled")
                return None

    def generate_mols_in_table(self, max_molecules: Optional[int] = None,
                            generate_conformers: bool = True,
                            nbr_confs: int = 10,
                            conf_percentile: float = 0.25) -> Optional[List[Dict[str, Any]]]:
        """
        Generate RDKit molecule objects from SMILES in a selected table.
        Uses existing helper methods for table selection and follows established patterns.
        
        Args:
            max_molecules (Optional[int]): Maximum number of molecules to process. If None, prompts user.
            generate_conformers (bool): Whether to generate 3D conformers for molecules
            conformer_count (int): Number of conformers to generate per molecule
            
        Returns:
            Optional[List[Dict]]: List of molecule dictionaries with RDKit objects and metadata
        """
        try:
            print(f"🧬 GENERATE MOLECULE OBJECTS FROM TABLE")
            print("=" * 60)
            
            # Use existing interactive table selection method
            selected_table = self._select_table_interactive(
                prompt_title="SELECT TABLE FOR MOLECULE GENERATION",
                show_compound_count=True
            )
            
            if not selected_table:
                print("❌ No table selected for molecule generation")
                return None
            
            # Get table data using existing method
            compounds_df = self._get_table_as_dataframe(selected_table)
            if compounds_df.empty:
                print(f"❌ No compounds found in table '{selected_table}'")
                return None
            
            # Validate SMILES column exists (reuse existing validation pattern)
            if 'smiles' not in compounds_df.columns:
                print("❌ No 'smiles' column found in the selected table")
                print(f"   Available columns: {list(compounds_df.columns)}")
                return None
            
            total_available = len(compounds_df)
            print(f"\n📊 Table Analysis:")
            print(f"   📋 Selected table: '{selected_table}'")
            print(f"   🧪 Available compounds: {total_available:,}")
            
            # Interactive molecule count selection if not provided
            if max_molecules is None:
                max_molecules = self._prompt_for_processing_count(total_available, "molecules to process")
                if max_molecules is None:
                    print("❌ Molecule generation cancelled by user")
                    return None
            
            # Handle special case: -1 means all molecules
            if max_molecules == -1:
                max_molecules = total_available
                print(f"🧪 Processing all {total_available:,} molecules in the table")
            elif max_molecules > 0:
                compounds_df = compounds_df.head(max_molecules)
                print(f"🔍 Limited to first {min(max_molecules, total_available):,} molecules")
            
            processing_count = len(compounds_df)
            print(f"🧬 Generating molecule objects for {processing_count:,} compounds...")
            
            # Show processing configuration
            print(f"\n⚙️  Processing Configuration:")
            print(f"   🧪 Generate conformers: {'Yes' if generate_conformers else 'No'}")
            if generate_conformers:
                print(f"   🔄 Conformers per molecule: {nbr_confs}")
            
            # Check RDKit availability (reuse existing pattern)
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
                from rdkit import RDLogger
                RDLogger.DisableLog('rdApp.*')  # Suppress RDKit warnings
            except ImportError:
                print("❌ RDKit not installed. Please install RDKit:")
                print("   conda install -c conda-forge rdkit")
                return None
            
            # Process molecules with progress tracking
            print(f"\n🔬 Processing molecules...")
            generated_molecules = []
            processing_stats = {
                'processed': 0,
                'successful': 0,
                'invalid_smiles': 0,
                'conformer_failures': 0,
                'property_calculations': 0
            }
            
            # Use existing progress bar pattern
            if TQDM_AVAILABLE:
                progress_bar = tqdm(
                    compounds_df.iterrows(),
                    total=len(compounds_df),
                    desc="Generating molecules",
                    unit="compounds"
                )
            else:
                progress_bar = compounds_df.iterrows()
                processed_count = 0
            
            self._check_sdf_blob_column(selected_table)
            
            for _, compound in progress_bar:
                processing_stats['processed'] += 1
                
                try:
                    smiles = compound['smiles']
                    compound_name = compound.get('name', f'compound_{compound.get("id", processing_stats["processed"])}')
                    compound_id = compound.get('id', processing_stats['processed'])
                    
                    # Generate RDKit molecule object
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        processing_stats['invalid_smiles'] += 1
                        continue
                    
                    # Add hydrogens for 3D conformer generation
                    mol_hs = Chem.AddHs(mol) if generate_conformers else mol
                    
                    # # Calculate molecular properties (reuse existing patterns)
                    # mol_properties = self._calculate_molecular_properties(mol, smiles)
                    # processing_stats['property_calculations'] += 1
                    
                    # Generate conformers if requested
                    conformers_data = []
                    conformer_success = True
                    
                    # Create molecule data dictionary
                    molecule_data = {
                        'compound_id': compound_id,
                        'compound_name': compound_name,
                        'smiles': smiles,
                        'rdkit_mol':mol,
                        'source_table': selected_table,
                        'original_flag': compound.get('flag', 'nd'),
                        'inchi_key': compound.get('inchi_key', None)
                    }
                    
                    if generate_conformers:
                        try:

                            mol, molecule_datamol = self._generate_molecule_conformer(
                                mol_hs, nbr_confs, compound_name, conf_percentile
                            )
                        
                            self._store_mol(mol, molecule_data)
                        
                        except Exception as e:
                            conformer_success = False
                            processing_stats['conformer_failures'] += 1
                    
                    generated_molecules.append(molecule_data)
                    processing_stats['successful'] += 1
                    
                except Exception as e:
                    # Log error but continue processing
                    continue
                
                # Progress update for non-tqdm case (reuse existing pattern)
                if not TQDM_AVAILABLE:
                    processed_count += 1
                    if processed_count % max(1, len(compounds_df) // 10) == 0:
                        print(f"   📊 Processed {processed_count}/{len(compounds_df)} compounds...")
            
            if TQDM_AVAILABLE and hasattr(progress_bar, 'close'):
                progress_bar.close()
            
            # Display processing results
            print(f"\n✅ Molecule generation completed!")
            print(f"   📊 Processing Statistics:")
            print(f"      • Total processed: {processing_stats['processed']:,}")
            print(f"      • Successful generations: {processing_stats['successful']:,}")
            print(f"      • Invalid SMILES: {processing_stats['invalid_smiles']:,}")
            print(f"      • Property calculations: {processing_stats['property_calculations']:,}")
            
            if generate_conformers:
                print(f"      • Conformer failures: {processing_stats['conformer_failures']:,}")
                successful_conformers = processing_stats['successful'] - processing_stats['conformer_failures']
                print(f"      • Successful conformer generations: {successful_conformers:,}")
            
            if processing_stats['successful'] > 0:
                success_rate = (processing_stats['successful'] / processing_stats['processed']) * 100
                print(f"   📈 Success rate: {success_rate:.1f}%")
            
            return generated_molecules if generated_molecules else None
            
        except Exception as e:
            print(f"❌ Error in generate_mols_in_table: {e}")
            return None

    def _prompt_for_processing_count(self, total_available: int, item_type: str = "items") -> Optional[int]:
        """
        Interactive prompt for the number of items to process.
        Reuses existing prompt patterns for consistency.
        
        Args:
            total_available (int): Total number of items available
            item_type (str): Description of items being processed
            
        Returns:
            Optional[int]: Number of items to process, or None if cancelled
        """
        try:
            print(f"\n📊 {item_type.upper()} COUNT SELECTION")
            print("=" * 50)
            print(f"📋 Total {item_type} available: {total_available:,}")
            print("\nOptions:")
            print(f"   • Enter a number (1-{total_available:,}) to limit {item_type}")
            print(f"   • Enter -1 to process ALL {item_type}")
            print(f"   • Enter 0 or 'cancel' to abort")
            
            # Suggest reasonable defaults based on size (reuse existing logic)
            if total_available <= 50:
                suggested = total_available
                print(f"\n💡 Suggestion: {suggested} (all {item_type})")
            elif total_available <= 200:
                suggested = 100
                print(f"\n💡 Suggestion: {suggested} (good sample size)")
            elif total_available <= 1000:
                suggested = 250
                print(f"\n💡 Suggestion: {suggested} (manageable sample)")
            else:
                suggested = 500
                print(f"\n💡 Suggestion: {suggested} (representative sample)")
            
            print("=" * 50)
            
            while True:
                try:
                    user_input = input(f"🔢 Number of {item_type} to process (default: {suggested}): ").strip()
                    
                    # Handle empty input (use default)
                    if not user_input:
                        print(f"✅ Using default: {suggested} {item_type}")
                        return suggested
                    
                    # Handle cancel (reuse existing pattern)
                    if user_input.lower() in ['cancel', 'quit', 'exit', 'abort']:
                        print(f"❌ {item_type} processing cancelled")
                        return None
                    
                    # Handle numeric input (reuse existing logic)
                    try:
                        count = int(user_input)
                        
                        if count == 0:
                            print(f"❌ Processing cancelled (0 {item_type} selected)")
                            return None
                        elif count == -1:
                            print(f"✅ Processing ALL {total_available:,} {item_type}")
                            return -1
                        elif 1 <= count <= total_available:
                            print(f"✅ Processing {count:,} {item_type}")
                            return count
                        else:
                            print(f"❌ Invalid count. Please enter 1-{total_available:,}, -1 for all, or 0 to cancel")
                            continue
                            
                    except ValueError:
                        print("❌ Please enter a valid number")
                        continue
                        
                except KeyboardInterrupt:
                    print(f"\n❌ {item_type} count selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error prompting for {item_type} count: {e}")
            return None

    def _generate_molecule_conformer(self, mol_hs, nbr_confs: int, 
                                    compound_name: str, conf_percentile) -> List[Dict[str, Any]]:
        """
        Generate 3D conformers for a molecule.
        Reuses existing conformer generation patterns from _generate_conformers.
        
        Args:
            mol_with_hs: RDKit molecule object with hydrogens
            conformer_count (int): Number of conformers to generate
            compound_name (str): Name of the compound
            
        Returns:
            List[Dict]: List of conformer data dictionaries
        """
        try:
            from rdkit.Chem import AllChem
            
            # Use existing conformer generation method
            mol_hs, confs, ps = self._generate_conformers(
                mol_hs, nbr_confs=nbr_confs, mmff='MMFF94s'
            )

            if not confs:
                return []
            
            # Create conformer data
            conformers_data = []
            for i, conf_id in enumerate(confs):
                try:
                    # Get conformer energy if MMFF properties (ps) available
                    energy = None
                    if ps is not None:
                        try:
                            # Use the molecule with hydrogens returned by _generate_conformers (mol_hs)
                            ff = AllChem.MMFFGetMoleculeForceField(mol_hs, ps, confId=int(conf_id))
                            if ff is not None:
                                energy = float(ff.CalcEnergy())
                        except Exception:
                            energy = None
                    
                    conformer_info = {
                        'conformer_id': conf_id,
                        'conformer_index': i,
                        'energy': energy,
                        'molecule_object': mol_hs,  # Reference to the molecule with this conformer
                        'generation_successful': True
                    }
                    
                    conformers_data.append(conformer_info)
                    
                except Exception as e:
                    # Add failed conformer entry
                    conformers_data.append({
                        'conformer_id': conf_id,
                        'conformer_index': i,
                        'energy': None,
                        'molecule_object': None,
                        'generation_successful': False,
                        'error': str(e)
                    })
            
            # Determine lowest-energy conformer object (may be None or a mol object)
            target_conformer = self._sort_conf_by_energy(conformers_data, mol_hs, ps, conf_percentile)
            
            target_conformer_renumbered = self._renumber_mol_object(target_conformer)

            # Also return the associated data dict for the lowest-energy conformer.
            # _sort_conf_by_energy marks the record with 'is_lowest_energy' when possible
            molecule_data = None
            if conformers_data:
                molecule_data = next((r for r in conformers_data if r.get('is_lowest_energy')), conformers_data[0])
            
            return target_conformer_renumbered, molecule_data
            
        except Exception as e:
            return []

    def _generate_conformers(self, mol_hs, nbr_confs, mmff):
        """
        Generate 3D conformers for a SMILES string with robust handling and faster defaults.

        Returns:
            (mol_hs, conf_ids, ps)
            - mol_hs: RDKit Mol with hydrogens (or None on failure)
            - conf_ids: list of conformer ids (could be empty)
            - ps: MMFF property object (or None)
        """
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            import multiprocessing

            
            # Prefer ETKDGv3 when available, fallback to ETKDG
            try:
                params = AllChem.ETKDGv3()
            except Exception:
                try:
                    params = AllChem.ETKDG()
                except Exception:
                    params = None

            # Configure embedding parameters if present
            prune_rms = 0.25
            max_attempts = 100
            threads = max(1, (multiprocessing.cpu_count() or 1) - 1)

            if params is not None:
                # defensive attribute setting (different RDKit versions expose different names)
                if hasattr(params, "pruneRmsThresh"):
                    params.pruneRmsThresh = float(prune_rms)
                if hasattr(params, "useRandomCoords"):
                    params.useRandomCoords = False  # deterministic & typically faster
                if hasattr(params, "numThreads"):
                    try:
                        params.numThreads = int(threads)
                    except Exception:
                        pass
                if hasattr(params, "maxAttempts"):
                    try:
                        params.maxAttempts = int(max_attempts)
                    except Exception:
                        pass
                # encourage better small-ring geometry handling when available
                if hasattr(params, "useSmallRingTorsions"):
                    try:
                        params.useSmallRingTorsions = True
                    except Exception:
                        pass

            # Embed multiple conformers (defensive calling to support different RDKit signatures)
            try:
                if params is not None:
                    conf_ids = list(AllChem.EmbedMultipleConfs(mol_hs, numConfs=nbr_confs, params=params))
                else:
                    conf_ids = list(AllChem.EmbedMultipleConfs(mol_hs, nbr_confs))
            except TypeError:
                # older signature
                try:
                    conf_ids = list(AllChem.EmbedMultipleConfs(mol_hs, nbr_confs, params))
                except Exception:
                    conf_ids = []
            except Exception:
                conf_ids = []

            if not conf_ids:
                return mol_hs, [], None

            # Prepare MMFF properties
            try:
                ps = AllChem.MMFFGetMoleculeProperties(mol_hs, mmffVariant='MMFF94s')
            except Exception:
                try:
                    ps = AllChem.MMFFGetMoleculeProperties(mol_hs)
                except Exception:
                    ps = None

            # Try bulk optimization; fall back to per-conformer optimize if signature differs
            try:
                # modern RDKit: MMFFOptimizeMoleculeConfs returns list of (confId, status)
                try:
                    AllChem.MMFFOptimizeMoleculeConfs(mol_hs, confIds=conf_ids, maxIters=100, mmffVariant='MMFF94s')
                except TypeError:
                    # older signature without keyword args
                    AllChem.MMFFOptimizeMoleculeConfs(mol_hs, conf_ids, 100)
            except Exception:
                # best-effort per-conformer optimization
                for cid in conf_ids:
                    try:
                        ff = AllChem.MMFFGetMoleculeForceField(mol_hs, ps, confId=int(cid))
                        if ff is not None:
                            ff.Minimize(maxIts=100)
                    except Exception:
                        continue

            return mol_hs, conf_ids, ps

        except Exception as error:
            
            print(f"Exiting conformer generation with failure... Error: {error}")
            
            return None, [], None

    def _sort_conf_by_energy(self, conformers_data, mol_hs, ps, conf_percentile):
        """
        Compute/collect energies for conformers_data, mark the lowest-energy conformer
        with 'is_lowest_energy' = True and return that conformer dict.
        The conformers_data list is sorted in-place by energy (lowest first, None last).
        Returns:
            dict or None: lowest energy conformer dict or None if none available
        """
        try:
            # Defensive defaults
            if not conformers_data:
                return None

            # Try to compute energies if missing, using RDKit MMFF if available
            try:
                rdkit_available = True
            except Exception:
                rdkit_available = False

            for rec in conformers_data:
                # Ensure conformer_id present
                cid = rec.get('conformer_id')
                # Normalize existing energy values to float or None
                energy = rec.get('energy', None)
                if energy is not None:
                    try:
                        rec['energy'] = float(energy)
                    except Exception:
                        rec['energy'] = None
                    continue

            # Collect records with valid energy
            valid_records = [r for r in conformers_data if r.get('energy') is not None]

            if valid_records:
                # Find minimum energy
                min_energy = min(r['energy'] for r in valid_records)
                # Mark lowest; if multiple equal energies, first one kept as lowest
                lowest = None
                for r in conformers_data:
                    is_low = (r.get('energy') is not None and abs(r['energy'] - min_energy) <= 1e-9)
                    r['is_lowest_energy'] = bool(is_low)
                    if is_low and lowest is None:
                        lowest = r

                # Sort conformers_data in-place: lowest energies first, None energies last
                conformers_data.sort(key=lambda r: (r.get('energy') is None, r.get('energy') if r.get('energy') is not None else float('inf')))

                # If all valid energies equal, fallback to lowest-energy behavior
                if not valid_records:
                    return None

                min_energy = min(r['energy'] for r in valid_records)
                max_energy = max(r['energy'] for r in valid_records)

                # Determine target energy along the range [min_energy, max_energy]
                if max_energy == min_energy:
                    target_energy = min_energy
                else:
                    target_energy = min_energy + conf_percentile * (max_energy - min_energy)

                # Find record whose energy is closest to the target_energy
                selected = min(valid_records, key=lambda r: abs((r.get('energy') or 0.0) - target_energy))

                # Mark selection flags for clarity
                for r in conformers_data:
                    r['is_selected_percentile'] = (r is selected)
                    # keep existing is_lowest_energy marking
                    if r.get('energy') is None:
                        r['is_lowest_energy'] = False

                # Return the molecule object for the selected conformer, fallback to lowest if neededprint(mol)
                # fallback to lowest-energy conformer object if selected has no molecule_object
                return lowest.get('molecule_object') if lowest else None

            else:
                # No energies available: annotate all as not lowest, pick first as fallback
                for i, r in enumerate(conformers_data):
                    r['energy'] = None
                    r['is_lowest_energy'] = (i == 0)
                return conformers_data[0] if conformers_data else None

        except Exception:
            # On any unexpected error, fail gracefully
            try:
                for i, r in enumerate(conformers_data or []):
                    r.setdefault('energy', None)
                    r['is_lowest_energy'] = (i == 0)
                return (conformers_data[0] if conformers_data else None)
            except Exception:
                return None
    
    def _renumber_mol_object(self, mol):
        """
        Efficiently annotate atoms in an RDKit Mol with human-readable labels.

        For each atom sets:
          - _Name (string) -> e.g. "C1", "O2"
          - atomLabel (string) -> same as _Name
          - atom map number (integer) -> idx+1 (uses SetAtomMapNum when available)

        Returns the same mol object (modified in place). Safe if mol is None.
        """
        
        if mol is None:
            return None

        # Localize method lookups for small speed gain
        get_atoms = mol.GetAtoms

        for atom in get_atoms():
            # Use GetIdx once, GetSymbol once
            idx = atom.GetIdx()
            symbol = atom.GetSymbol() or ""
            label = f"{symbol}{idx + 1}"

            # Prefer integer atom map for fast numeric access; ignore errors if not supported
            try:
                atom.SetAtomMapNum(idx + 1)
            except Exception:
                # Some RDKit contexts may not allow or need atom map numbers — ignore safely
                pass

            # Set textual properties; use try/except to avoid unexpected failures
            try:
                atom.SetProp("_Name", label)
            except Exception:
                pass

            try:
                atom.SetProp("atomLabel", label)
            except Exception:
                pass

        return mol

    def _check_sdf_blob_column(self, source_table) -> bool:
        """
        Checks for existing sdf_blob column and prompts user for deletion if needed.
        
        """

        try:
            if not source_table:
                print("❌ Error: No 'source_table' specified in molecule_data")
                return
            
            # Check if table exists
            available_tables = self.get_all_tables()
            if source_table not in available_tables:
                print(f"❌ Error: Table '{source_table}' does not exist in chemspace database")
                return
            
            # Check if sdf_blob column exists
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            cursor.execute(f"PRAGMA table_info({source_table})")
            columns = [row[1] for row in cursor.fetchall()]
            
            sdf_blob_exists = 'sdf_blob' in columns
            conn.close()
            
            if sdf_blob_exists:
                print(f"⚠️  Column 'sdf_blob' already exists in table '{source_table}'")
                
                # Count existing blobs
                conn = sqlite3.connect(self.__chemspace_db)
                cursor = conn.cursor()
                cursor.execute(f"SELECT COUNT(*) FROM {source_table} WHERE sdf_blob IS NOT NULL")
                existing_blobs = cursor.fetchone()[0]
                conn.close()
                
                if existing_blobs > 0:
                    print(f"   📊 Found {existing_blobs} existing SDF blob(s) in the table")
                    
                    while True:
                        user_choice = input(f"❓ Do you want to delete the existing 'sdf_blob' column? (y/n): ").strip().lower()
                        
                        if user_choice in ['y', 'yes']:
                            print(f"🗑️  Removing 'sdf_blob' column from table '{source_table}'...")
                            success = self._remove_sdf_blob_column(source_table)
                            if success:
                                print(f"   ✅ Successfully removed 'sdf_blob' column")
                                break
                            else:
                                print(f"   ❌ Failed to remove 'sdf_blob' column")
                                return
                        elif user_choice in ['n', 'no']:
                            print(f"   ⚠️  Keeping existing 'sdf_blob' column - new blobs will be added/updated")
                            break
                        else:
                            print("❗ Please answer 'y'/'yes' or 'n'/'no'")
                else:
                    print(f"   ℹ️  Column exists but contains no data - proceeding with molecule storage")
                    
        except Exception as e:
            print(f"❌ Error checking sdf_blob column: {e}")

    def _store_mol(self, mol, molecule_data):
        """
        Store an RDKit molecule as an SDF blob in the database.
        
        
        Args:
            mol: RDKit molecule object
            molecule_data: Dictionary containing molecule metadata
        """
        try:
            # Proceed with storing the molecule
            destination_filename = self._create_sdf_file(mol, molecule_data)
            self._save_mol_to_table(destination_filename, molecule_data)
            
        except Exception as e:
            print(f"❌ Error in _store_mol: {e}")

    def _remove_sdf_blob_column(self, table_name: str) -> bool:
        """
        Remove the sdf_blob column from a table by recreating the table without it.
        SQLite doesn't support DROP COLUMN directly, so we use the standard approach.
        
        Args:
            table_name (str): Name of the table to modify
            
        Returns:
            bool: True if column was removed successfully
        """
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get current table schema
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns_info = cursor.fetchall()
            
            # Filter out the sdf_blob column
            filtered_columns = []
            for col_info in columns_info:
                col_name = col_info[1]
                col_type = col_info[2]
                col_notnull = col_info[3]
                col_default = col_info[4]
                col_pk = col_info[5]
                
                if col_name.lower() != 'sdf_blob':
                    # Reconstruct column definition
                    col_def = f"{col_name} {col_type}"
                    
                    if col_notnull:
                        col_def += " NOT NULL"
                    
                    if col_default is not None:
                        col_def += f" DEFAULT {col_default}"
                    
                    if col_pk:
                        col_def += " PRIMARY KEY"
                    
                    filtered_columns.append(col_def)
            
            if not filtered_columns:
                print(f"❌ Error: Cannot remove sdf_blob column - no other columns found")
                conn.close()
                return False
            
            # Create new table without sdf_blob column
            temp_table_name = f"{table_name}_temp_no_sdf"
            
            # Add unique constraint if it was in original table
            unique_constraint = ""
            cursor.execute(f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}'")
            original_sql = cursor.fetchone()
            if original_sql and 'UNIQUE' in original_sql[0]:
                # Try to preserve UNIQUE constraints (simplified approach)
                if 'UNIQUE(smiles, name)' in original_sql[0] or 'UNIQUE(name, smiles)' in original_sql[0]:
                    unique_constraint = ", UNIQUE(smiles, name)"
            
            create_query = f"""
            CREATE TABLE {temp_table_name} (
                {', '.join(filtered_columns)}{unique_constraint}
            )
            """
            
            cursor.execute(create_query)
            
            # Copy data from original table (excluding sdf_blob column)
            column_names = [col.split()[0] for col in filtered_columns]  # Extract just column names
            columns_str = ', '.join(column_names)
            
            cursor.execute(f"""
            INSERT INTO {temp_table_name} ({columns_str})
            SELECT {columns_str} FROM {table_name}
            """)
            
            # Drop original table
            cursor.execute(f"DROP TABLE {table_name}")
            
            # Rename temp table to original name
            cursor.execute(f"ALTER TABLE {temp_table_name} RENAME TO {table_name}")
            
            # Recreate indexes (basic ones)
            try:
                cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_smiles ON {table_name}(smiles)")
                cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_name ON {table_name}(name)")
                if 'flag' in [col.split()[0].lower() for col in filtered_columns]:
                    cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_flag ON {table_name}(flag)")
                if 'inchi_key' in [col.split()[0].lower() for col in filtered_columns]:
                    cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_inchi_key ON {table_name}(inchi_key)")
            except Exception as e:
                # Continue even if index creation fails
                print(f"   ⚠️  Warning: Some indexes could not be recreated: {e}")
            
            conn.commit()
            conn.close()
            
            print(f"   🔄 Table '{table_name}' recreated without 'sdf_blob' column")
            return True
            
        except Exception as e:
            print(f"❌ Error removing sdf_blob column from table '{table_name}': {e}")
            try:
                conn.close()
            except:
                pass
            return False
        
    def _create_sdf_file(self, mol, molecule_data):
        
        """Write single RDKit molecule to SDF file with atom names preserved"""
        from rdkit.Chem import SDWriter
        
        # Ensure atom names are stored as molecule properties before writing
        atom_names = []
        
        for atom in mol.GetAtoms():
            if atom.HasProp("_Name"):
                atom_names.append(atom.GetProp("_Name"))
            else:
                atom_names.append(f"{atom.GetSymbol()}{atom.GetIdx() + 1}")
        
        # Store atom names as a property in the molecule
        mol.SetProp("AtomNames", "|".join(atom_names))
        
        # Also store individual atom names as properties
        for i, name in enumerate(atom_names):
            mol.SetProp(f"Atom_{i}_Name", name)
        
        inchi_key = str(molecule_data.get('inchi_key') or '')
        if not inchi_key:
            inchi_key = 'unknown_inchi'
        
        destination_filename = f"/tmp/{inchi_key}.sdf"
        
        writer = SDWriter(destination_filename)
        writer.write(mol)
        writer.close()
        
        return destination_filename
    
    def _save_mol_to_table(self, filename: str, molecule_data: Dict[str, Any]) -> bool:
        """
        Will receive a filename and molecule_data dictionary, and will store the filename as a blob object 
        to the table whose name is in molecule_data['source_table']. The row in which the file is stored 
        is selected matching the molecule_data['inchi_key'] with the corresponding column in the table.
        
        Args:
            filename (str): Path to the SDF file to store as blob
            molecule_data (Dict[str, Any]): Dictionary containing molecule metadata including 'source_table' and 'inchi_key'
            
        Returns:
            bool: True if the file was successfully stored, False otherwise
        """
        try:
            # Extract required information from molecule_data
            source_table = molecule_data.get('source_table')
            inchi_key = molecule_data.get('inchi_key')
            
            if not source_table:
                print("❌ Error: No 'source_table' specified in molecule_data")
                return False
                
            if not inchi_key:
                print("❌ Error: No 'inchi_key' specified in molecule_data")
                return False
            
            # Check if the file exists
            if not os.path.exists(filename):
                print(f"❌ Error: SDF file not found: {filename}")
                return False
            
            # Check if the table exists
            available_tables = self.get_all_tables()
            if source_table not in available_tables:
                print(f"❌ Error: Table '{source_table}' does not exist in chemspace database")
                return False
            
            #print(f"💾 Storing SDF file as blob in table '{source_table}' for InChI key: {inchi_key}")
            
            # Read the SDF file as binary data
            try:
                with open(filename, 'rb') as sdf_file:
                    sdf_blob = sdf_file.read()
            except Exception as e:
                print(f"❌ Error reading SDF file '{filename}': {e}")
                return False
            
            # Connect to database and update the table
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # First, check if the table has a 'sdf_blob' column, add it if it doesn't exist
            cursor.execute(f"PRAGMA table_info({source_table})")
            columns = [row[1] for row in cursor.fetchall()]
            
            if 'sdf_blob' not in columns:
                #print(f"   🔧 Adding 'sdf_blob' column to table '{source_table}'")
                try:
                    cursor.execute(f"ALTER TABLE {source_table} ADD COLUMN sdf_blob BLOB")
                    conn.commit()
                except Exception as e:
                    print(f"❌ Error adding sdf_blob column: {e}")
                    conn.close()
                    return False
            
            # Check if a row with the matching inchi_key exists
            cursor.execute(f"SELECT COUNT(*) FROM {source_table} WHERE inchi_key = ?", (inchi_key,))
            row_count = cursor.fetchone()[0]
            
            if row_count == 0:
                print(f"⚠️  Warning: No row found with inchi_key '{inchi_key}' in table '{source_table}'")
                conn.close()
                return False
            
            if row_count > 1:
                print(f"⚠️  Warning: Multiple rows found with inchi_key '{inchi_key}' in table '{source_table}'. Updating the first match.")
            
            # Update the row with the SDF blob data
            update_query = f"UPDATE {source_table} SET sdf_blob = ? WHERE inchi_key = ?"
            cursor.execute(update_query, (sdf_blob, inchi_key))
            
            # Check if the update was successful
            updated_rows = cursor.rowcount
            
            conn.commit()
            conn.close()
            
            if updated_rows > 0:
                #print(f"   ✅ Successfully stored SDF blob ({len(sdf_blob)} bytes) for InChI key: {inchi_key}")
                
                # Clean up temporary file after successful storage
                try:
                    os.unlink(filename)
                    #print(f"   🧹 Cleaned up temporary SDF file: {filename}")
                except Exception as e:
                    print(f"   ⚠️  Warning: Could not delete temporary file '{filename}': {e}")
                
                return True
            else:
                print(f"❌ Error: No rows were updated in table '{source_table}'")
                return False
                
        except sqlite3.Error as e:
            print(f"❌ Database error storing SDF blob: {e}")
            return False
        except Exception as e:
            print(f"❌ Error storing SDF file to table: {e}")
            return False
        

    def subset_table(self):
        """
        Create a subset of a database table based on user-defined filter criteria.
        This method provides an interactive workflow for users to:
        1. Select a table from the chemspace database
        2. Choose a column to filter on
        3. Specify filter values for that column
        4. Choose whether to INCLUDE or EXCLUDE rows matching the filter values
        5. Create a new table containing only the filtered rows
        The process includes validation at each step and displays relevant metadata
        (e.g., row counts, unique value counts) to help users make informed selections.
        The filtered subset is saved as a new table in the database with the same schema
        as the original table, including indexes on 'smiles' and 'name' columns.
        Workflow:
            - Lists all available tables with compound counts
            - Displays column information (type, unique values) for the selected table
            - Accepts user input for filter values with data type parsing
            - Prompts for INCLUSION or EXCLUSION mode
            - Shows filtering statistics (original rows, filtered rows, retention rate)
            - Prompts for output table name (with timestamp-based default)
            - Creates new table and inserts filtered data, skipping duplicates
            - Displays final operation summary
        Returns:
            None. Creates a new table in the database as a side effect.
        Raises:
            Handled internally - prints error messages on failure rather than raising exceptions.
            KeyboardInterrupt: Gracefully cancels the operation at any interactive prompt.
        Notes:
            - User can cancel at any selection prompt (table, column, filter values, or mode)
            - Filter values are case-insensitive for string matching
            - Numeric values are parsed according to the column's data type
            - User must choose between INCLUSION (keep matching rows) or EXCLUSION (remove matching rows)
            - Duplicate rows (determined by UNIQUE constraints) are skipped during insertion
            - Table names are sanitized before use
            - The process requires an active database connection
            - If the output table already exists, user is prompted to replace it
        """

        
        try:
            # Get available tables
            available_tables = self.get_all_tables()
            
            if not available_tables:
                print("❌ No tables available in chemspace database")
                return
            
            # Interactive table selection
            print(f"\n🔍 TABLE SUBSETTING")
            print("=" * 70)
            print("Select a table to subset:")
            print("-" * 70)
            
            table_info = []
            for i, table_name in enumerate(available_tables, 1):
                try:
                    compound_count = self.get_compound_count(table_name=table_name)
                    print(f"{i:3d}. {table_name:<35} ({compound_count:>6,} compounds)")
                    table_info.append({
                        'index': i,
                        'name': table_name,
                        'count': compound_count
                    })
                except Exception as e:
                    print(f"{i:3d}. {table_name:<35} (Error: {e})")
                    table_info.append({
                        'index': i,
                        'name': table_name,
                        'count': 0,
                        'error': str(e)
                    })
            
            print("-" * 70)
            
            # Get table selection from user
            while True:
                try:
                    selection = input("Select table (number or name, or 'cancel'): ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Table subsetting cancelled")
                        return
                    
                    # Try as number
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(available_tables):
                            selected_table = available_tables[table_idx]
                            break
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(available_tables)}")
                            continue
                    except ValueError:
                        # Try as table name
                        matching_tables = [t for t in available_tables if t.lower() == selection.lower()]
                        if matching_tables:
                            selected_table = matching_tables[0]
                            break
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                except KeyboardInterrupt:
                    print("\n❌ Table selection cancelled")
                    return
            
            print(f"\n✅ Selected table: '{selected_table}'")
            
            # Get table data
            table_df = self._get_table_as_dataframe(selected_table)
            
            if table_df.empty:
                print(f"❌ No data found in table '{selected_table}'")
                return
            
            # Display available columns
            print(f"\n📋 COLUMN SELECTION")
            print("=" * 70)
            print(f"Available columns in table '{selected_table}':")
            print("-" * 70)
            
            available_columns = list(table_df.columns)
            for i, col_name in enumerate(available_columns, 1):
                unique_count = table_df[col_name].nunique()
                col_type = str(table_df[col_name].dtype)
                print(f"{i:2d}. {col_name:<20} ({col_type:<10}) - {unique_count:,} unique values")
            
            print("-" * 70)
            
            # Get column selection from user
            while True:
                try:
                    col_selection = input("Select column for filtering (number or name): ").strip()
                    
                    # Try as number
                    try:
                        col_idx = int(col_selection) - 1
                        if 0 <= col_idx < len(available_columns):
                            selected_column = available_columns[col_idx]
                            break
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(available_columns)}")
                            continue
                    except ValueError:
                        # Try as column name
                        if col_selection in available_columns:
                            selected_column = col_selection
                            break
                        else:
                            print(f"❌ Column '{col_selection}' not found")
                            continue
                except KeyboardInterrupt:
                    print("\n❌ Column selection cancelled")
                    return
            
            print(f"\n✅ Selected column: '{selected_column}'")
            
            # Show unique values in the column
            unique_values = sorted(table_df[selected_column].unique())
            print(f"\n📊 Unique values in column '{selected_column}':")
            print("=" * 70)
            
            # Display unique values (limit to first 50 to avoid clutter)
            max_display = min(50, len(unique_values))
            for i, value in enumerate(unique_values[:max_display], 1):
                value_str = str(value)[:50] + "..." if len(str(value)) > 50 else str(value)
                print(f"   {i:3d}. {value_str}")
            
            if len(unique_values) > max_display:
                print(f"   ... and {len(unique_values) - max_display} more values")
            
            print("-" * 70)
            
            # Collect filter values from user
            filter_values = []
            filter_counter = 1
            
            print(f"\n🔍 FILTERING VALUES INPUT")
            print("=" * 70)
            print(f"Enter values to filter by in column '{selected_column}'")
            print("Enter -1 when finished adding values")
            print("=" * 70)
            
            while True:
                try:
                    user_input = input(f"\n🔢 Filter value {filter_counter} (or -1 to finish): ").strip()
                    
                    if user_input == "-1":
                        if not filter_values:
                            print("❌ You must enter at least one filter value")
                            continue
                        print(f"\n✅ Finished entering filter values")
                        break
                    
                    if not user_input:
                        print("⚠️  Please enter a value or -1 to finish")
                        continue
                    
                    # Try to match the value in the column
                    # Handle different data types
                    matching_values = []
                    for unique_val in unique_values:
                        if str(unique_val).lower() == user_input.lower():
                            matching_values.append(unique_val)
                    
                    if matching_values:
                        selected_value = matching_values[0]
                        filter_values.append(selected_value)
                        print(f"   ✅ Added: {selected_value}")
                        filter_counter += 1
                    else:
                        # Try to parse as the column's data type
                        try:
                            col_dtype = table_df[selected_column].dtype
                            
                            if col_dtype in ['int64', 'int32', 'int16', 'int8']:
                                parsed_value = int(user_input)
                            elif col_dtype in ['float64', 'float32']:
                                parsed_value = float(user_input)
                            else:
                                parsed_value = user_input
                            
                            filter_values.append(parsed_value)
                            print(f"   ✅ Added: {parsed_value}")
                            filter_counter += 1
                            
                        except ValueError:
                            print(f"❌ Cannot parse '{user_input}' as {col_dtype}. Please try again.")
                            continue
                
                except KeyboardInterrupt:
                    print("\n❌ Filter value input cancelled")
                    return
            
            # Prompt user for INCLUSION or EXCLUSION mode
            print(f"\n🔄 FILTER MODE SELECTION")
            print("=" * 70)
            print(f"How should the filter be applied to column '{selected_column}'?")
            print("\nOptions:")
            print("   1. INCLUSION: Keep rows WHERE column is IN {filter_values}")
            print("      (Create a subset containing only the specified values)")
            print(f"      Example: Keep rows where '{selected_column}' = {filter_values[:2] if len(filter_values) > 1 else filter_values}")
            print()
            print("   2. EXCLUSION: Keep rows WHERE column is NOT IN {filter_values}")
            print("      (Create a subset excluding the specified values)")
            print(f"      Example: Remove rows where '{selected_column}' = {filter_values[:2] if len(filter_values) > 1 else filter_values}")
            print("=" * 70)
            
            filter_mode = None
            while filter_mode is None:
                try:
                    mode_choice = input("\nSelect filter mode (1 for INCLUSION, 2 for EXCLUSION, or 'cancel'): ").strip()
                    
                    if mode_choice.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Filter mode selection cancelled")
                        return
                    
                    if mode_choice == '1':
                        filter_mode = 'inclusion'
                        print("✅ Selected: INCLUSION mode (keep matching rows)")
                    elif mode_choice == '2':
                        filter_mode = 'exclusion'
                        print("✅ Selected: EXCLUSION mode (remove matching rows)")
                    else:
                        print("❌ Invalid selection. Please enter 1, 2, or 'cancel'")
                        
                except KeyboardInterrupt:
                    print("\n❌ Filter mode selection cancelled")
                    return
            
            # Apply filters to create subset based on selected mode
            print(f"\n📊 APPLYING FILTERS")
            print("=" * 70)
            
            if filter_mode == 'inclusion':
                print(f"Filtering table '{selected_table}' WHERE '{selected_column}' IN {filter_values}")
                subset_df = table_df[table_df[selected_column].isin(filter_values)]
            else:  # exclusion mode
                print(f"Filtering table '{selected_table}' WHERE '{selected_column}' NOT IN {filter_values}")
                subset_df = table_df[~table_df[selected_column].isin(filter_values)]
            
            print(f"✅ Filter applied successfully")
            print(f"   📊 Original rows: {len(table_df):,}")
            print(f"   🔍 Filtered rows: {len(subset_df):,}")
            print(f"   📈 Retention rate: {(len(subset_df)/len(table_df))*100:.2f}%")
            print(f"   🔄 Filter mode: {filter_mode.upper()}")
            
            if subset_df.empty:
                print("❌ No rows matched the filter criteria")
                return
            
            # Get output table name from user
            print(f"\n💾 OUTPUT TABLE CONFIGURATION")
            print("=" * 70)
            
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            default_table_name = f"{selected_table}_{filter_mode}_{timestamp}"
            
            user_table_name = input(f"Enter table name for subset (default: {default_table_name}): ").strip()
            output_table_name = user_table_name if user_table_name else default_table_name
            output_table_name = self._sanitize_table_name(output_table_name)
            
            # Check if output table already exists
            available_tables = self.get_all_tables()
            
            while output_table_name in available_tables:
                print(f"\n⚠️  Table '{output_table_name}' already exists!")
                
                # Ask user for a new name
                user_table_name = input(f"Enter a different table name (or 'cancel' to abort): ").strip()
                
                if user_table_name.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Table subsetting cancelled by user")
                    return
                
                output_table_name = self._sanitize_table_name(user_table_name)
                
                # Refresh available tables list
                available_tables = self.get_all_tables()
                
                if output_table_name in available_tables:
                    print(f"⚠️  Table '{output_table_name}' still exists. Please try a different name.")
                else:
                    print(f"✅ New table name accepted: '{output_table_name}'")
                
            print(f"✅ Output table name: '{output_table_name}'")
            
            # Save subset to database
            print(f"\n💾 SAVING SUBSET TO DATABASE")
            print("-" * 70)
            
            try:
                conn = sqlite3.connect(self.__chemspace_db)
                cursor = conn.cursor()
                
                # Create new table with same schema as original
                # Get original table schema
                cursor.execute(f"PRAGMA table_info({selected_table})")
                columns_info = cursor.fetchall()
                
                # Build CREATE TABLE statement
                column_defs = []
                for col_info in columns_info:
                    col_name = col_info[1]
                    col_type = col_info[2]
                    col_notnull = col_info[3]
                    col_default = col_info[4]
                    col_pk = col_info[5]
                    
                    col_def = f"{col_name} {col_type}"
                    
                    if col_notnull:
                        col_def += " NOT NULL"
                    
                    # Only add DEFAULT when non-empty; quote string defaults if needed
                    if col_default is not None and str(col_default).strip() != "":
                        default_str = str(col_default).strip()
                        if default_str.upper() == "NULL":
                            col_def += " DEFAULT NULL"
                        elif default_str.replace(".", "", 1).lstrip("-").isdigit():
                            col_def += f" DEFAULT {default_str}"
                        else:
                            if not (default_str.startswith("'") or default_str.startswith('"')):
                                default_str = f"'{default_str}'"
                            col_def += f" DEFAULT {default_str}"
                    
                    if col_pk:
                        col_def += " PRIMARY KEY AUTOINCREMENT"
                    
                    column_defs.append(col_def)
                
                # Add UNIQUE constraint if present
                unique_constraint = ""
                cursor.execute(f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{selected_table}'")
                original_sql = cursor.fetchone()
                if original_sql and 'UNIQUE' in original_sql[0]:
                    if 'UNIQUE(smiles, name)' in original_sql[0]:
                        unique_constraint = ", UNIQUE(smiles, name)"
                
                create_table_sql = f"CREATE TABLE {output_table_name} ({', '.join(column_defs)}{unique_constraint})"

                cursor.execute(create_table_sql)
                
                # Create indexes
                cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_smiles ON {output_table_name}(smiles)")
                cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{output_table_name}_name ON {output_table_name}(name)")
                
                # Insert subset data
                insert_columns = [col[1] for col in columns_info if col[1] != 'id']  # Exclude auto-increment id
                placeholders = ', '.join(['?'] * len(insert_columns))
                insert_sql = f"INSERT INTO {output_table_name} ({', '.join(insert_columns)}) VALUES ({placeholders})"
                
                inserted_count = 0
                for _, row in subset_df.iterrows():
                    try:
                        values = [row[col] for col in insert_columns]
                        cursor.execute(insert_sql, values)
                        inserted_count += 1
                    except sqlite3.IntegrityError:
                        continue  # Skip duplicates
                
                conn.commit()
                conn.close()
                
                print(f"✅ Successfully created subset table '{output_table_name}'")
                print(f"   📊 Rows inserted: {inserted_count:,}")
                print(f"   🔍 Filter criteria: {selected_column} {filter_mode.upper()} {filter_values}")
                print(f"   📈 Original table: {selected_table}")
                print(f"   🔄 Filter mode: {filter_mode.upper()}")
                
            except Exception as e:
                print(f"❌ Error creating subset table: {e}")
                return

        except Exception as e:
            print(f"❌ Error in subset_table method: {e}")


