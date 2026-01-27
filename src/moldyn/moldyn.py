from tidyscreen import tidyscreen
import sys
import os

# Add the parent directory to path to import our local tidyscreen module
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir)

# Import from our local tidyscreen module
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject


class MolDyn:
    """
    MolDyn class for managing molecular dynamics assays within a project.
    Uses an ActivateProject object to access project information and database functionality.
    
    """
    def __init__(self, project_obj: ActivateProject):
        """
        Initialize MolDyn with an ActivateProject object.
        
        Args:
            project_obj (ActivateProject): An instantiated ActivateProject object
        """
        # Validate that we received a proper ActivateProject object
        if not isinstance(project_obj, ActivateProject):
            raise TypeError("MolDyn requires an ActivateProject object")
        
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
        
        # Set up receptors path within the project directory
        self.__receptor_path = os.path.join(self.path, 'docking/receptors')
        
        # Set up docking assays registers database path within the project directory
        self.__docking_registers_db = os.path.join(self.path, 'docking/docking_registers', 'registers.db')

        # Set up docking params registers database path within the project directory
        self.__docking_params_db = os.path.join(self.path, 'docking/params', 'params.db')

        # Set up molecular dynamics assays registers database path within the project directory
        self.__md_registers_db = os.path.join(self.path, 'dynamics/md_registers', 'md_registers.db')

        # Set up molecular dynamics assays folder path within the project directory
        self.__md_assays_folder = os.path.join(self.path, 'dynamics/md_assays')

        # Set up molecular dynamics params folder path within the project directory
        self.__md_params_folder = os.path.join(self.path, 'dynamics/md_params')

    
    def create_md_method():

        """
        Placeholder method for creating a molecular dynamics method.
        """
        pass