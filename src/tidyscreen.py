import sys
import os
import json
from tidyscreen.databases.DatabaseManager import DatabaseManager
from tidyscreen.projects.projects_management import ProjectsManagement
import site 

def projects():
    """Will list all projects stored in the projects database."""
    projects_manager = ProjectsManagement()
    projects_manager.list_all_projects()  # This will print the formatted output

class ActivateProject:
    
    def __init__(self, name):
        # Initialize database manager as a component
        self.db_manager = DatabaseManager()
        
        self.__projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db" 
        self.name = name
        
        # Check if database exists, create it if it doesn't
        if not os.path.exists(self.__projects_db):
            print(f"Projects database does not exist at: {self.__projects_db}")
            print("Creating new projects database...")
            
            if not self.db_manager.create_projects_database(self.__projects_db):
                print("❌ Failed to create projects database. Cannot proceed.")
                # Set attributes to indicate project doesn't exist
                self._project_exists = False
                return
            
            print("✅ Projects database created successfully.")
        
        # Validate database structure
        self.db_manager.check_db(self.__projects_db)
        
        # Load project information from database using ProjectsManagement
        self.load_project_info()
        
    def load_project_info(self):
        """
        Load project information using the ProjectsManagement class.
        """
        projects_manager = ProjectsManagement()
        project_data = projects_manager.load_project_info(self.name)
        
        # Set attributes from the returned dictionary
        for key, value in project_data.items():
            setattr(self, key, value)
    
    def get_project_attributes(self):
        """
        Return a dictionary of all project attributes.
        """
        attributes = {}
        for attr_name in dir(self):
            if not attr_name.startswith('_') and not callable(getattr(self, attr_name)):
                attributes[attr_name] = getattr(self, attr_name)
        return attributes
    
    def project_exists(self):
        """
        Check if the project exists in the database.
        """
        return getattr(self, '_project_exists', False)
    
    # ChemSpace Integration Methods
    def _get_chemspace_db_path(self):
        """
        Get the path to the ChemSpace database for this project.
        """
        return os.path.join(self.path, 'chemspace', 'processed_data', 'chemspace.db')
    
    
    def get_chemspace_compounds(self, limit=None, flag_filter=None, name_pattern=None):
        """
        Retrieve compounds from the ChemSpace database with optional filtering.
        
        Args:
            limit (int, optional): Maximum number of compounds to return
            flag_filter (str, optional): Filter by specific flag value
            name_pattern (str, optional): Filter by name pattern (SQL LIKE)
            
        Returns:
            list: List of compound dictionaries
        """
        import sqlite3
        
        db_path = self._get_chemspace_db_path()
        
        if not os.path.exists(db_path):
            print("⚠️  ChemSpace database not found. Load some compounds first.")
            return []
        
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Build query with filters
            query = "SELECT * FROM compounds WHERE 1=1"
            params = []
            
            if flag_filter:
                query += " AND flag = ?"
                params.append(flag_filter)
            
            if name_pattern:
                query += " AND name LIKE ?"
                params.append(f"%{name_pattern}%")
            
            query += " ORDER BY id DESC"
            
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
            print(f"❌ Error retrieving compounds: {e}")
            return []
    
    def print_chemspace_summary(self):
        """
        Print a formatted summary of the ChemSpace database.
        """
        import sqlite3
        
        db_path = self._get_chemspace_db_path()
        
        if not os.path.exists(db_path):
            print("⚠️  ChemSpace database not found. Load some compounds first.")
            return
        
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Total compounds
            cursor.execute("SELECT COUNT(*) FROM compounds")
            total_compounds = cursor.fetchone()[0]
            
            # Compounds by flag
            cursor.execute("SELECT flag, COUNT(*) FROM compounds GROUP BY flag ORDER BY COUNT(*) DESC")
            flag_counts = dict(cursor.fetchall())
            
            # Unique source files
            cursor.execute("SELECT COUNT(DISTINCT source_file) FROM compounds")
            unique_files = cursor.fetchone()[0]
            
            conn.close()
            
            print("\n" + "="*60)
            print(f"CHEMSPACE SUMMARY - Project: {self.name}")
            print("="*60)
            print(f"📊 Total Compounds: {total_compounds}")
            print(f"📁 Source Files: {unique_files}")
            print(f"🗄️  Database Path: {db_path}")
            
            if flag_counts:
                print("\n🏷️  Compounds by Flag:")
                for flag, count in flag_counts.items():
                    flag_display = flag if flag else '(no flag)'
                    print(f"   {flag_display}: {count}")
            
            print("="*60)
            
        except Exception as e:
            print(f"❌ Error getting ChemSpace summary: {e}")

class CreateProject:
    
    def __init__(self, name, path, description=None):
        # Initialize database manager as a component
        self.db_manager = DatabaseManager()
        
        self.__projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db" 
        self.name = name
        self.description = description
        
        # Append project name as the last subdirectory
        self.path = os.path.join(path, name)
        
        # Check if database exists, create it if it doesn't
        if not os.path.exists(self.__projects_db):
            print(f"Projects database does not exist at: {self.__projects_db}")
            print("Creating new projects database...")
            
            if not self.db_manager.create_projects_database(self.__projects_db):
                print("❌ Failed to create projects database. Cannot proceed.")
                self._project_created = False
                return
            
            print("✅ Projects database created successfully.")
        
        # Validate database structure
        self.db_manager.check_db(self.__projects_db)
        
        # Create the project directory first
        self.create_project_directory()
        
        # Create the project entry in database
        self.create_project_entry()
    
    def load_project_template(self):
        """
        Load the project structure template from JSON configuration file.
        """
        try:
            # Get the path to the configuration file
            config_path = os.path.join(os.path.dirname(__file__), 'config', 'project_template.json')
            
            with open(config_path, 'r', encoding='utf-8') as file:
                template = json.load(file)
                return template.get('project_structure', {})
        except FileNotFoundError:
            print(f"⚠️  Warning: Project template not found at {config_path}")
            print("Creating basic project directory structure...")
            return self._get_default_structure()
        except json.JSONDecodeError as e:
            print(f"⚠️  Warning: Error parsing project template: {e}")
            print("Using default project directory structure...")
            return self._get_default_structure()
        except Exception as e:
            print(f"⚠️  Warning: Error loading project template: {e}")
            print("Using default project directory structure...")
            return self._get_default_structure()
    
    def _get_default_structure(self):
        """
        Return a default project structure if the JSON file is not available.
        """
        return {
            "folders": ["src", "tests", "docs"],
            "subfolders": {
                "src": ["models", "views", "controllers"],
                "tests": ["unit", "integration"]
            },
            "files": {
                "README.md": f"# {self.name}\n\n## Description\n{self.description or 'Project description'}\n",
                ".gitignore": "__pycache__/\n*.pyc\n*.pyo\n*.pyd\n.Python\nenv/\nvenv/\n",
                "requirements.txt": "# Add your project dependencies here\n"
            }
        }
    
    def create_project_directory(self):
        """
        Create the project directory structure on the filesystem based on JSON template.
        """
        try:
            # Create the main project directory
            os.makedirs(self.path, exist_ok=True)
            print(f"📁 Project directory created: {self.path}")
            
            # Load the project template
            template = self.load_project_template()
            
            # Create folders
            self._create_folders(template.get('folders', []))
            
            # Create subfolders
            self._create_subfolders(template.get('subfolders', {}))
            
            # Create template files
            self._create_template_files(template.get('files', {}))
            
            self._directory_created = True
            print(f"✅ Project structure created successfully!")
            
        except PermissionError:
            print(f"❌ Error: Permission denied when creating directory: {self.path}")
            self._directory_created = False
        except Exception as e:
            print(f"❌ Error creating project directory: {e}")
            self._directory_created = False
    
    def _create_folders(self, folders):
        """
        Create the main folders for the project.
        """
        for folder in folders:
            folder_path = os.path.join(self.path, folder)
            os.makedirs(folder_path, exist_ok=True)
            print(f"  📂 Created folder: {folder}")
    
    def _create_subfolders(self, subfolders):
        """
        Create subfolders within the main folders.
        """
        for parent_folder, child_folders in subfolders.items():
            for child_folder in child_folders:
                subfolder_path = os.path.join(self.path, parent_folder, child_folder)
                os.makedirs(subfolder_path, exist_ok=True)
                print(f"    📁 Created subfolder: {parent_folder}/{child_folder}")
    
    def _create_template_files(self, files):
        """
        Create template files with content.
        """
        for file_path, content in files.items():
            full_file_path = os.path.join(self.path, file_path)
            
            # Create directory if it doesn't exist (for files in subfolders)
            file_dir = os.path.dirname(full_file_path)
            if file_dir != self.path:
                os.makedirs(file_dir, exist_ok=True)
            
            # Replace template variables
            processed_content = content.replace('{{project_name}}', self.name)
            processed_content = processed_content.replace('{{project_description}}', self.description or 'Project description')
            
            # Write the file
            with open(full_file_path, 'w', encoding='utf-8') as file:
                file.write(processed_content)
            print(f"  📄 Created file: {file_path}")
    
    def create_project_entry(self):
        """
        Create a new project entry in the database with name and path using ProjectsManagement.
        """
        # Only proceed with database entry if directory was created successfully
        if not getattr(self, '_directory_created', False):
            print("❌ Project creation failed: Could not create project directory")
            self._project_created = False
            return
        
        # Use ProjectsManagement to handle project creation
        projects_manager = ProjectsManagement()
        result = projects_manager.create_project_entry(self.name, self.path, self.description)
        
        # Handle the result
        if result['success']:
            self.id = result['project_id']
            self.created_date = result['created_date']
            self._project_created = True
            
            print(result['message'])
            print(f"Project path: {self.path}")
        else:
            self._project_created = False
            print(result['message'])
    
    def directory_created_successfully(self):
        """
        Check if the project directory was created successfully.
        """
        return getattr(self, '_directory_created', False)

class DeleteProject:
    
    def __init__(self, name, delete_directory=True, confirm=True):
        # Initialize database manager as a component
        self.db_manager = DatabaseManager()
        
        self.__projects_db = f"{site.getsitepackages()[0]}/tidyscreen/projects_db/projects_database.db"
        self.name = name
        self.delete_directory = delete_directory
        self.confirm = confirm
        
        # Check if database exists, create it if it doesn't
        if not os.path.exists(self.__projects_db):
            print(f"Projects database does not exist at: {self.__projects_db}")
            print("Creating new projects database...")
            
            if not self.db_manager.create_projects_database(self.__projects_db):
                print("❌ Failed to create projects database. Cannot proceed.")
                self._project_deleted = False
                return
            
            print("✅ Projects database created successfully.")
        
        # Validate database structure
        self.db_manager.check_db(self.__projects_db)
        
        # Delete the project entry from database
        self.delete_project_entry()
    
    def delete_project_entry(self):
        """
        Delete a project entry from the database with confirmation using ProjectsManagement.
        """
        # Use ProjectsManagement to get project information first
        projects_manager = ProjectsManagement()
        
        # Load project info to show what will be deleted
        project_info = projects_manager.load_project_info(self.name)
        
        if not project_info.get('_project_exists', False):
            print(f"Project '{self.name}' not found in database.")
            self._project_deleted = False
            return
        
        # Show project information and ask for confirmation
        if self.confirm:
            print("\n" + "="*60)
            print("PROJECT DELETION CONFIRMATION")
            print("="*60)
            print(f"Project Name: {project_info.get('name', 'N/A')}")
            print(f"Project ID: {project_info.get('id', 'N/A')}")
            print(f"Project Path: {project_info.get('path', 'N/A')}")
            print(f"Description: {project_info.get('description', 'N/A')}")
            print(f"Created Date: {project_info.get('created_date', 'N/A')}")
            print("\nActions to be performed:")
            print(f"- Remove project record from database: YES")
            print(f"- Delete project directory: {'YES' if self.delete_directory else 'NO'}")
            
            if self.delete_directory and project_info.get('path'):
                print(f"- Directory to be deleted: {project_info['path']}")
                print("\n⚠️  WARNING: This will permanently delete all files in the project directory!")
            
            print("="*60)
            
            confirmation = input("\nAre you sure you want to delete this project? (yes/no): ").strip().lower()
            
            if confirmation not in ['yes', 'y']:
                print("Project deletion cancelled.")
                self._project_deleted = False
                return
        
        # Proceed with deletion
        result = projects_manager.delete_project_entry(self.name, self.delete_directory)
        
        # Handle the result
        if result['success']:
            self.id = result['project_info'].get('id')
            self.path = result['project_info'].get('path')
            self.description = result['project_info'].get('description')
            self.created_date = result['project_info'].get('created_date')
            self.directory_deleted = result['directory_deleted']
            self._project_deleted = True
            
            print(f"\n✅ {result['message']}")
            
        else:
            self._project_deleted = False
            print(f"\n❌ {result['message']}")
    
    def project_deleted_successfully(self):
        """
        Check if the project was deleted successfully.
        """
        return getattr(self, '_project_deleted', False)


