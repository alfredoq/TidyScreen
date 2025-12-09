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
            
            # Create chem_filters table if it doesn't exist
            create_chem_filters_query = """
            CREATE TABLE IF NOT EXISTS chem_reactions (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                reaction_name TEXT UNIQUE NOT NULL,
                smarts TEXT NOT NULL
            )
            """
            cursor.execute(create_chem_filters_query)
            
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
            
            # Get existing filter names to avoid duplicates
            cursor.execute("SELECT reaction_name FROM chem_reactions")
            existing_reactions = set(row[0] for row in cursor.fetchall())
            
            # Insert new filters
            insert_reaction_query = """
            INSERT OR IGNORE INTO chem_reactions (reaction_name, smarts) VALUES (?, ?)
            """
            
            new_reactions_count = 0
            updated_reactions_count = 0
            
            for reaction_item in reactions_data:
                reaction_name = reaction_item['reaction_name']
                smarts = reaction_item['smarts']
                
                if reaction_name not in existing_reactions:
                    # New filter
                    cursor.execute(insert_reaction_query, (reaction_name, smarts))
                    new_reactions_count += 1
                    print(f"  ➕ Added new reaction: {reaction_name}")
                else:
                    # Check if SMARTS pattern has changed
                    cursor.execute("SELECT smarts FROM chem_reactions WHERE reaction_name = ?", (reaction_name,))
                    current_smarts = cursor.fetchone()[0]
                    
                    if current_smarts != smarts:
                        # Update existing reaction with new SMARTS pattern
                        cursor.execute(
                            "UPDATE chem_reactions SET smarts = ? WHERE reaction_name = ?", 
                            (smarts, reaction_name)
                        )
                        updated_filters_count += 1
                        print(f"  🔄 Updated reaction: {reaction_name}")
            
            conn.commit()
            conn.close()
            
            # Summary
            print(f"\n✅ Chemical reactions refresh completed:")
            print(f"   📊 Total reactions in JSON: {len(reactions_data)}")
            print(f"   ➕ New reactions added: {new_reactions_count}")
            print(f"   🔄 Reactions updated: {updated_reactions_count}")
            print(f"   📍 Source file: {json_file_path}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error refreshing chemical reactions: {e}")
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




