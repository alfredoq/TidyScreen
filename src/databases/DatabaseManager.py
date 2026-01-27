import sqlite3
import os
import sys
import json

class DatabaseManager:
    
    def check_db(self, db):
        if os.path.exists(db):
            print("TidyScreen main database found! Continuing...")
            return True
        else:
            print("TidyScreen main database does NOT exist. Creating it...")
            sys.exit(1)
    
    def create_projects_database(self, db_path):
        """Create the projects database with the required table structure."""
        try:
            # Ensure the directory exists
            db_dir = os.path.dirname(db_path)
            os.makedirs(db_dir, exist_ok=True)
            
            # Connect to database (this will create it if it doesn't exist)
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Create projects table
            create_table_query = """
            CREATE TABLE IF NOT EXISTS projects (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT UNIQUE NOT NULL,
                path TEXT NOT NULL,
                description TEXT,
                created_date TEXT NOT NULL
            )
            """
            
            cursor.execute(create_table_query)
            
            # Create chem_filters table
            create_chem_filters_query = """
            CREATE TABLE IF NOT EXISTS chem_filters (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                filter_name TEXT UNIQUE NOT NULL,
                smarts TEXT NOT NULL
            )
            """
            
            cursor.execute(create_chem_filters_query)
            
            # Populate chem_filters table from JSON file
            config_dir = os.path.dirname(os.path.dirname(os.path.dirname(db_path)))
            json_file_path = os.path.join(config_dir, 'src', 'config', 'chem_filters.json')
            
            if os.path.exists(json_file_path):
                with open(json_file_path, 'r') as f:
                    filters_data = json.load(f)
                
                insert_filter_query = """
                INSERT OR IGNORE INTO chem_filters (filter_name, smarts) VALUES (?, ?)
                """
                
                for filter_item in filters_data:
                    cursor.execute(insert_filter_query, (filter_item['filter_name'], filter_item['smarts']))
                    
                print(f"Chemical filters populated from: {json_file_path}")
            else:
                print(f"Warning: Chemical filters JSON file not found at {json_file_path}")
            
            conn.commit()
            conn.close()
            
            print(f"Database created successfully at: {db_path}")
            return True
            
        except Exception as e:
            print(f"Error creating database: {e}")
            return False
            
    def refresh_chem_filters(self, db_path):
        """Refresh the chem_filters table by reading the JSON file again and adding new filters."""
        
        print("ENTRO")
        
        try:
            # Connect to database
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Create chem_filters table if it doesn't exist
            create_chem_filters_query = """
            CREATE TABLE IF NOT EXISTS chem_filters (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                filter_name TEXT UNIQUE NOT NULL,
                smarts TEXT NOT NULL
            )
            """
            cursor.execute(create_chem_filters_query)
            
            # Get the path to the JSON file
            config_dir = os.path.dirname(os.path.dirname(os.path.dirname(db_path)))
            json_file_path = os.path.join(config_dir, 'tidyscreen', 'config', 'chem_filters.json')
            
            if not os.path.exists(json_file_path):
                print(f"Warning: Chemical filters JSON file not found at {json_file_path}")
                conn.close()
                return False
            
            # Read the JSON file
            with open(json_file_path, 'r') as f:
                filters_data = json.load(f)
            
            # Get existing filter names to avoid duplicates
            cursor.execute("SELECT filter_name FROM chem_filters")
            existing_filters = set(row[0] for row in cursor.fetchall())
            
            # Insert new filters
            insert_filter_query = """
            INSERT OR IGNORE INTO chem_filters (filter_name, smarts) VALUES (?, ?)
            """
            
            new_filters_count = 0
            updated_filters_count = 0
            
            for filter_item in filters_data:
                filter_name = filter_item['filter_name']
                smarts = filter_item['smarts']
                
                if filter_name not in existing_filters:
                    # New filter
                    cursor.execute(insert_filter_query, (filter_name, smarts))
                    new_filters_count += 1
                    print(f"  ➕ Added new filter: {filter_name}")
                else:
                    # Check if SMARTS pattern has changed
                    cursor.execute("SELECT smarts FROM chem_filters WHERE filter_name = ?", (filter_name,))
                    current_smarts = cursor.fetchone()[0]
                    
                    if current_smarts != smarts:
                        # Update existing filter with new SMARTS pattern
                        cursor.execute(
                            "UPDATE chem_filters SET smarts = ? WHERE filter_name = ?", 
                            (smarts, filter_name)
                        )
                        updated_filters_count += 1
                        print(f"  🔄 Updated filter: {filter_name}")
            
            conn.commit()
            conn.close()
            
            # Summary
            print(f"\n✅ Chemical filters refresh completed:")
            print(f"   📊 Total filters in JSON: {len(filters_data)}")
            print(f"   ➕ New filters added: {new_filters_count}")
            print(f"   🔄 Filters updated: {updated_filters_count}")
            print(f"   📍 Source file: {json_file_path}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error refreshing chemical filters: {e}")
            return False

    def refresh_chem_reactions(self, db_path):
        """Refresh the chem_reactions table by reading the JSON file again and adding new reaction types."""
        try:
            # Connect to database
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Create chem_reactions table if it doesn't exist
            create_chem_reactions_query = """
            CREATE TABLE IF NOT EXISTS chem_reactions (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                reaction_name TEXT UNIQUE NOT NULL,
                smarts TEXT NOT NULL
            )
            """
            cursor.execute(create_chem_reactions_query)
            
            # Get the path to the JSON file
            config_dir = os.path.dirname(os.path.dirname(os.path.dirname(db_path)))
            json_file_path = os.path.join(config_dir, 'tidyscreen', 'config', 'chem_reactions.json')
            
            if not os.path.exists(json_file_path):
                print(f"Warning: Chemical reactions JSON file not found at {json_file_path}")
                conn.close()
                return False
            
            # Read the JSON file
            with open(json_file_path, 'r') as f:
                reactions_data = json.load(f)
            
            # Get existing reactions with their SMARTS patterns
            cursor.execute("SELECT reaction_name, smarts FROM chem_reactions")
            existing_reactions = {row[0]: row[1] for row in cursor.fetchall()}
            
            # Insert new reactions
            insert_reaction_query = """
            INSERT OR IGNORE INTO chem_reactions (reaction_name, smarts) VALUES (?, ?)
            """
            
            new_reactions_count = 0
            updated_reactions_count = 0
            skipped_reactions_count = 0
            
            for reaction_item in reactions_data:
                reaction_name = reaction_item['reaction_name']
                new_smarts = reaction_item['smarts']
                
                if reaction_name not in existing_reactions:
                    # New reaction - add it
                    cursor.execute(insert_reaction_query, (reaction_name, new_smarts))
                    new_reactions_count += 1
                    print(f"  ➕ Added new reaction: {reaction_name}")
                    
                else:
                    # Reaction exists - check if SMARTS pattern is different
                    current_smarts = existing_reactions[reaction_name]
                    
                    if current_smarts != new_smarts:
                        # SMARTS pattern has changed - ask user for confirmation
                        print(f"\n🔄 Reaction '{reaction_name}' already exists with different SMARTS:")
                        print(f"   📋 Current SMARTS: {current_smarts}")
                        print(f"   🆕 New SMARTS:     {new_smarts}")
                        print(f"   📍 Source: {json_file_path}")
                        
                        while True:
                            user_choice = input("\n❓ Do you want to replace the existing SMARTS pattern? (y/n/s=skip): ").strip().lower()
                            
                            if user_choice in ['y', 'yes']:
                                # U pdate existing reaction with new SMARTS pattern
                                cursor.execute(
                                    "UPDATE chem_reactions SET smarts = ? WHERE reaction_name = ?", 
                                    (new_smarts, reaction_name)
                                )
                                updated_reactions_count += 1
                                print(f"  ✅ Updated reaction: {reaction_name}")
                                break
                                
                            elif user_choice in ['n', 'no']:
                                # Keep existing SMARTS pattern
                                print(f"  ⏭️  Keeping existing SMARTS for: {reaction_name}")
                                skipped_reactions_count += 1
                                break
                                
                            elif user_choice in ['s', 'skip']:
                                # Skip this reaction entirely
                                print(f"  ⏭️  Skipped reaction: {reaction_name}")
                                skipped_reactions_count += 1
                                break
                                
                            else:
                                print("  ❌ Invalid choice. Please enter 'y' (yes), 'n' (no), or 's' (skip)")
                    else:
                        # SMARTS pattern is the same - no action needed
                        print(f"  ✅ Reaction '{reaction_name}' already up to date")
            
            conn.commit()
            conn.close()
            
            # Summary
            print(f"\n✅ Chemical reactions refresh completed:")
            print(f"   📊 Total reactions in JSON: {len(reactions_data)}")
            print(f"   ➕ New reactions added: {new_reactions_count}")
            print(f"   🔄 Reactions updated: {updated_reactions_count}")
            print(f"   ⏭️  Reactions skipped: {skipped_reactions_count}")
            print(f"   📍 Source file: {json_file_path}")
            
            return True
            
        except KeyboardInterrupt:
            print(f"\n\n⏹️  Operation cancelled by user")
            if 'conn' in locals():
                conn.close()
            return False
            
        except Exception as e:
            print(f"❌ Error refreshing chemical reactions: {e}")
            if 'conn' in locals():
                conn.close()
            return False

    def connect_db(self, db):
        try:
            conn = sqlite3.connect(db)
            print(f"Connected to database at {db} successfully.")
            return conn
        except sqlite3.Error as e:
            print(f"An error occurred while connecting to the database: {e}")
            sys.exit(1)
            
    def execute_query(self, conn, query, message):
        try:
            cursor = conn.cursor()
            cursor.execute(query)
            conn.commit()
            print(f"{message}")
        except sqlite3.Error as e:
            print(f"An error occurred while executing the query: {e}")
            sys.exit(1)

def connect_db(db):
    try:
        conn = sqlite3.connect(db)
        print(f"Connected to database at {db} successfully.")
        return conn
    except sqlite3.Error as e:
        print(f"An error occurred while connecting to the database: {e}")
        sys.exit(1)


def create_table_from_columns_dict(cursor: sqlite3.Cursor, table_name: str, columns_dict: dict, verbose=True) -> None:
    """
    Create a database table dynamically from a columns dictionary.
    
    Args:
        cursor (sqlite3.Cursor): Database cursor for executing SQL commands
        table_name (str): Name of the table to create
        columns_dict (dict): Dictionary mapping column names to their SQL types
                            Example: {'id': 'INTEGER PRIMARY KEY AUTOINCREMENT', 'name': 'TEXT NOT NULL'}
    
    Returns:
        None
    
    Raises:
        sqlite3.Error: If table creation fails
    """
    try:
        # Build the CREATE TABLE statement from columns_dict - This works if the table does not exist
        columns_sql = ",\n    ".join([f"{col_name} {col_type}" for col_name, col_type in columns_dict.items()])
        
        create_table_sql = f"""
            CREATE TABLE IF NOT EXISTS {table_name} (
                {columns_sql}
            )
        """
        
        cursor.execute(create_table_sql)
        
        if verbose:
            print(f"   ✅ Table '{table_name}' created/verified with {len(columns_dict)} columns")
        
    except sqlite3.Error as e:
        print(f"   ❌ Error creating table '{table_name}': {e}")
        raise


def update_legacy_table_columns(cursor: sqlite3.Cursor, table_name: str, columns_dict: dict, verbose=True) -> None:
    """
    Update legacy database tables by adding missing columns based on columns_dict reference.
    
    This method compares existing table columns against the expected columns_dict schema
    and adds any missing columns to maintain compatibility with newer database versions.
    
    Args:
        cursor (sqlite3.Cursor): Database cursor for executing SQL commands
        table_name (str): Name of the table to update
        columns_dict (dict): Dictionary mapping column names to their SQL types
                            Example: {'renumbering_dict': 'TEXT', 'status': "TEXT DEFAULT 'created'"}
    
    Returns:
        None
    
    Raises:
        sqlite3.Error: If column addition fails
    
    Example:
        columns_dict = {'renumbering_dict': 'TEXT', 'notes': 'TEXT'}
        self._update_legacy_table_columns(cursor, 'pdb_templates', columns_dict)
    """
    try:
        # Get existing columns from the table
        cursor.execute(f"PRAGMA table_info({table_name})")
        existing_columns = {row[1] for row in cursor.fetchall()}
        
        # Identify missing columns
        missing_columns = set(columns_dict.keys()) - existing_columns
        
        if not missing_columns:
            if verbose:
                print(f"   ✅ Table '{table_name}' is up to date - no missing columns")
                return
        
        if verbose:
            print(f"   🔄 Updating legacy table '{table_name}' - adding {len(missing_columns)} missing column(s)")
        
        # Add each missing column
        for col_name in missing_columns:
            col_type = columns_dict[col_name]
            try:
                alter_sql = f"ALTER TABLE {table_name} ADD COLUMN {col_name} {col_type}"
                cursor.execute(alter_sql)
                print(f"      ✓ Added column: {col_name} ({col_type})")
            except sqlite3.Error as e:
                print(f"      ⚠️  Could not add column '{col_name}': {e}")
        
        if verbose:
            print(f"   ✅ Legacy table '{table_name}' updated successfully")
        
    except sqlite3.Error as e:
        print(f"   ❌ Error updating legacy table '{table_name}': {e}")
        raise

def remove_legacy_table_columns(cursor: sqlite3.Cursor, table_name: str, columns_dict: dict, verbose=True) -> None:
    """
    Remove legacy columns from database tables that are no longer in the columns_dict reference.
    
    This method compares existing table columns against the expected columns_dict schema
    and removes any extra columns that are no longer needed. This is useful for cleaning up
    deprecated columns from legacy database versions.
    
    **WARNING**: This operation is destructive and cannot be easily reversed. Column data
    will be permanently deleted.
    
    Args:
        cursor (sqlite3.Cursor): Database cursor for executing SQL commands
        table_name (str): Name of the table to clean up
        columns_dict (dict): Dictionary mapping column names to their SQL types
                            Example: {'pdb_id': 'INTEGER PRIMARY KEY', 'pdb_name': 'TEXT'}
    
    Returns:
        None
    
    Raises:
        sqlite3.Error: If column removal fails
    
    Example:
        columns_dict = {'pdb_id': 'INTEGER PRIMARY KEY', 'pdb_name': 'TEXT'}
        self._remove_legacy_table_columns(cursor, 'pdb_templates', columns_dict)
    
    Note:
        SQLite does not support DROP COLUMN directly in older versions (< 3.35.0).
        This method uses the recommended approach of creating a new table and copying data.
    """
    try:
        # Get existing columns from the table
        cursor.execute(f"PRAGMA table_info({table_name})")
        existing_columns_info = cursor.fetchall()
        existing_columns = {row[1] for row in existing_columns_info}
        
        # Identify columns to remove
        columns_to_remove = existing_columns - set(columns_dict.keys())
        
        if not columns_to_remove:
            if verbose:
                print(f"   ✅ Table '{table_name}' has no extra columns to remove")
            return
        
        if verbose:
            print(f"   🗑️  Removing {len(columns_to_remove)} legacy column(s) from '{table_name}'")
            print(f"      Columns to remove: {', '.join(sorted(columns_to_remove))}")
        
        # Confirm destructive operation
        print(f"   ⚠️  WARNING: This will permanently delete column data!")
        confirm = input(f"   Continue with column removal? (yes/no, default: no): ").strip().lower()
        
        if confirm != 'yes':
            print(f"   ❌ Column removal cancelled by user")
            return
        
        # SQLite approach: Create new table with desired columns, copy data, rename
        # Get columns to keep (in original order)
        columns_to_keep = [
            (row[1], row[2])  # (column_name, column_type)
            for row in existing_columns_info 
            if row[1] in columns_dict
        ]
        
        # Build new table schema
        new_columns_sql = ",\n    ".join([
            f"{col_name} {columns_dict[col_name]}" 
            for col_name, _ in columns_to_keep
        ])
        
        # Column names for data transfer
        transfer_columns = ", ".join([col_name for col_name, _ in columns_to_keep])
        
        # Create temporary table with new schema
        temp_table_name = f"{table_name}_temp_migration"
        
        cursor.execute(f"""
            CREATE TABLE {temp_table_name} (
                {new_columns_sql}
            )
        """)
        
        # Copy data from old table to new table
        cursor.execute(f"""
            INSERT INTO {temp_table_name} ({transfer_columns})
            SELECT {transfer_columns}
            FROM {table_name}
        """)
        
        # Drop old table
        cursor.execute(f"DROP TABLE {table_name}")
        
        # Rename new table to original name
        cursor.execute(f"ALTER TABLE {temp_table_name} RENAME TO {table_name}")
        
        print(f"   ✅ Successfully removed legacy columns from '{table_name}'")
        for col_name in sorted(columns_to_remove):
            print(f"      ✓ Removed column: {col_name}")
        
    except sqlite3.Error as e:
        print(f"   ❌ Error removing legacy columns from '{table_name}': {e}")
        # Attempt rollback if temp table exists
        try:
            cursor.execute(f"DROP TABLE IF EXISTS {table_name}_temp_migration")
        except:
            pass
        raise
    
def insert_data_dinamically_into_table(cursor: sqlite3.Cursor, table_name: str, data_dict: dict) -> None:
    """
    Insert data dynamically into a database table.
    Args:
        cursor (sqlite3.Cursor): Database cursor for executing SQL commands
        table_name (str): Name of the table to insert data into 
        data (Dict[str, Any]): Dictionary mapping column names to their values
                                    Example: {'pdb_template_name': 'template1', 'pdb_model_name':
                                    'model1', 'project_name': 'MyProject'}  
    """

    # Build dynamic INSERT statement (excluding auto-increment and timestamp columns)
    insert_columns = [col for col in data_dict.keys()]
    placeholders = ', '.join(['?' for _ in insert_columns])
    columns_str = ', '.join(insert_columns)
    values_tuple = tuple(data_dict[col] for col in insert_columns)

    # Insert receptor entry using dynamic column list
    cursor.execute(f'''
    INSERT INTO {table_name} ({columns_str})
    VALUES ({placeholders})
    ''', values_tuple)

    # Retrieve the pdb_id of the newly inserted record
    pdb_template_name = data_dict.get('pdb_template_name', None)
    pdb_model_name = data_dict.get('pdb_model_name', None)

    template_id = cursor.lastrowid or cursor.execute(
        f"SELECT pdb_id FROM {table_name} WHERE pdb_template_name = ? AND pdb_model_name = ?", 
        (pdb_template_name, pdb_model_name)
    ).fetchone()[0]

    return template_id