import sqlite3
import os
import sys

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
            conn.commit()
            conn.close()
            
            print(f"Database created successfully at: {db_path}")
            return True
            
        except Exception as e:
            print(f"Error creating database: {e}")
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




