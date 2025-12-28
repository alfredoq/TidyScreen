from tidyscreen.databases.DatabaseManager import DatabaseManager
import site
from datetime import datetime
import os
import shutil


class ProjectsManagement(DatabaseManager):
    
    def __init__(self):
        # Initialize parent class first
        super().__init__()
        
        self.__projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"

    def list_all_projects(self, print_output=True):
        """
        List all projects stored in the projects database.

        Args:
            print_output (bool): Whether to print formatted output. Default is True.
                               When True, only prints formatted output and returns None.
                               When False, returns list of dictionaries without printing.

        Returns:
            list or None: List of dictionaries when print_output=False, None when print_output=True
        """
        # Check if database exists, create it if it doesn't
        if not os.path.exists(self.__projects_db):
            if print_output:
                print(f"Database does not exist at: {self.__projects_db}")
                print("Creating new projects database...")

            if not self.create_projects_database(self.__projects_db):
                if print_output:
                    print("Failed to create database. Cannot proceed.")
                    return None
                return []

            if print_output:
                print("Database created successfully. No projects found yet.")
                return None
            return []

        try:
            # Validate database structure
            self.check_db(self.__projects_db)

            # Connect to database
            conn = self.connect_db(self.__projects_db)
            cursor = conn.cursor()

            # Query all projects
            query = "SELECT * FROM projects ORDER BY created_date DESC"
            cursor.execute(query)

            # Fetch all project data
            projects_data = cursor.fetchall()

            if not projects_data:
                if print_output:
                    print("No projects found in the database.")
                    conn.close()
                    return None
                conn.close()
                return []

            # Get column names for better display
            column_names = [description[0] for description in cursor.description]

            if print_output:
                # Table formatting
                col_widths = [max(len(str(col)), 12) for col in column_names]
                for row in projects_data:
                    for i, val in enumerate(row):
                        col_widths[i] = max(col_widths[i], len(str(val)) if val is not None else 4)

                # Print header
                header = " | ".join(f"{col:<{col_widths[i]}}" for i, col in enumerate(column_names))
                print("=" * (sum(col_widths) + 3 * (len(column_names) - 1)))
                print(header)
                print("-" * (sum(col_widths) + 3 * (len(column_names) - 1)))

            projects_list = []
            for project_row in projects_data:
                # Create dictionary for each project
                project_dict = {}
                for i, column_name in enumerate(column_names):
                    project_dict[column_name] = project_row[i]

                projects_list.append(project_dict)

                if print_output:
                    row_str = " | ".join(
                        f"{str(project_row[i]) if project_row[i] is not None else '':<{col_widths[i]}}"
                        for i in range(len(column_names))
                    )
                    print(row_str)

            if print_output:
                print("=" * (sum(col_widths) + 3 * (len(column_names) - 1)))

            # Close connection
            conn.close()

            # Return None when printing formatted output, return list when not printing
            if print_output:
                return None
            else:
                return projects_list

        except Exception as e:
            if print_output:
                print(f"Error reading projects database: {e}")
                return None
            return []


    def load_project_info(self, project_name):
        """
        Retrieve project information from the database and return as dictionary.
        """
        try:
            # Connect to the database using inherited method
            conn = self.connect_db(self.__projects_db)
            
            # Query for project information
            cursor = conn.cursor()
            query = "SELECT * FROM projects WHERE name = ?"
            cursor.execute(query, (project_name,))
            
            # Fetch the project data
            project_data = cursor.fetchone()
            
            if project_data:
                # Get column names to map data to attributes
                column_names = [description[0] for description in cursor.description]
                
                # Create dictionary with project information
                project_dict = {}
                for i, column_name in enumerate(column_names):
                    project_dict[column_name] = project_data[i]
                
                print(f"Project '{project_name}' loaded successfully from database.")
                project_dict['_project_exists'] = True
                
                # Close the connection
                conn.close()
                
                return project_dict
                
            else:
                print(f"Project '{project_name}' not found in database.")
                # Close the connection
                conn.close()
                
                # Return default attributes
                return {
                    'id': None,
                    'name': project_name,
                    'path': None,
                    'description': None,
                    'created_date': None,
                    '_project_exists': False
                }
                
        except Exception as e:
            print(f"Error loading project information: {e}")
            # Return default attributes in case of error
            return {
                'id': None,
                'name': project_name,
                'path': None,
                'description': None,
                'created_date': None,
                '_project_exists': False
            }
    
    def create_project_entry(self, name, path, description=None):
        """
        Create a new project entry in the database with name and path.
        
        Args:
            name (str): The name of the project
            path (str): The path to the project
            description (str, optional): The description of the project
            
        Returns:
            dict: Dictionary containing project creation results with keys:
                - success (bool): Whether the project was created successfully
                - project_id (int): The ID of the created project (if successful)
                - created_date (str): The creation date (if successful)
                - message (str): Status message
        """
        try:
            # Connect to the database using inherited method
            conn = self.connect_db(self.__projects_db)
            
            # Check if project already exists
            cursor = conn.cursor()
            check_query = "SELECT name FROM projects WHERE name = ?"
            cursor.execute(check_query, (name,))
            
            if cursor.fetchone():
                conn.close()
                return {
                    'success': False,
                    'project_id': None,
                    'created_date': None,
                    'message': f"Project '{name}' already exists in database."
                }
            
            # Insert new project entry
            created_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            
            insert_query = """
                INSERT INTO projects (name, path, description, created_date) 
                VALUES (?, ?, ?, ?)
            """
            
            cursor.execute(insert_query, (name, path, description, created_date))
            conn.commit()
            
            # Get the ID of the newly created project
            project_id = cursor.lastrowid
            
            # Close the connection
            conn.close()
            
            return {
                'success': True,
                'project_id': project_id,
                'created_date': created_date,
                'message': f"Project '{name}' created successfully with ID {project_id}"
            }
            
        except Exception as e:
            return {
                'success': False,
                'project_id': None,
                'created_date': None,
                'message': f"Error creating project entry: {e}"
            }

    def delete_project_entry(self, project_name, delete_directory=True):
        """
        Delete a project entry from the database and optionally remove its directory.
        
        Args:
            project_name (str): The name of the project to delete
            delete_directory (bool): Whether to delete the project directory as well
            
        Returns:
            dict: Dictionary containing deletion results with keys:
                - success (bool): Whether the project was deleted successfully
                - project_info (dict): Information about the deleted project (if successful)
                - directory_deleted (bool): Whether the directory was deleted
                - message (str): Status message
        """
        try:
            # First, get project information before deletion
            project_info = self.load_project_info(project_name)
            
            if not project_info.get('_project_exists', False):
                return {
                    'success': False,
                    'project_info': None,
                    'directory_deleted': False,
                    'message': f"Project '{project_name}' not found in database."
                }
            
            # Connect to the database using inherited method
            conn = self.connect_db(self.__projects_db)
            cursor = conn.cursor()
            
            # Delete project from database
            delete_query = "DELETE FROM projects WHERE name = ?"
            cursor.execute(delete_query, (project_name,))
            
            if cursor.rowcount == 0:
                conn.close()
                return {
                    'success': False,
                    'project_info': project_info,
                    'directory_deleted': False,
                    'message': f"Project '{project_name}' could not be deleted from database."
                }
            
            conn.commit()
            conn.close()
            
            # Handle directory deletion if requested and path exists
            directory_deleted = False
            directory_message = ""
            
            if delete_directory and project_info.get('path'):
                project_path = project_info['path']
                
                if os.path.exists(project_path):
                    try:
                        if os.path.isdir(project_path):
                            shutil.rmtree(project_path)
                            directory_deleted = True
                            directory_message = f" Directory '{project_path}' deleted successfully."
                        else:
                            directory_message = f" Path '{project_path}' is not a directory, skipped deletion."
                    except Exception as dir_error:
                        directory_message = f" Failed to delete directory '{project_path}': {dir_error}"
                else:
                    directory_message = f" Directory '{project_path}' does not exist, skipped deletion."
            elif delete_directory:
                directory_message = " No path specified in project, skipped directory deletion."
            
            return {
                'success': True,
                'project_info': project_info,
                'directory_deleted': directory_deleted,
                'message': f"Project '{project_name}' deleted successfully from database.{directory_message}"
            }
            
        except Exception as e:
            return {
                'success': False,
                'project_info': None,
                'directory_deleted': False,
                'message': f"Error deleting project entry: {e}"
            }
