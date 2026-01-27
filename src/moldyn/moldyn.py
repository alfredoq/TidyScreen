from tidyscreen import tidyscreen
import sys
import os
from databases import DatabaseManager as dbm

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

    
    def create_md_method(self):
        """
        Uses a similar approach to docking methods creation in moldock.py
        to create molecular dynamics methods.
        
        Will prompt the user for the following method parameters:
        - MD engine: options: GROMACS, AMBER, NAMD, OpenMM
        
        The method will create a md_methods.db under project path in the dynamics/md_registers folder
        """
        try:
            import sqlite3
            import json
            from datetime import datetime
            
            print(f"🧬 CREATE MOLECULAR DYNAMICS METHOD")
            print("=" * 50)
            
            # Ensure MD directories exist
            md_registers_dir = os.path.dirname(self.__md_registers_db)
            os.makedirs(md_registers_dir, exist_ok=True)
            
            # Define available MD engines with details
            md_engines = {
                '1': {
                    'name': 'AMBER',
                    'description': 'Comprehensive MD suite for biomolecular simulations',
                    'requirements': ['AMBER installation', 'AmberTools', 'Force field parameters'],
                },
            }
            
            # Display available MD engines
            print(f"📋 Available MD Engines:")
            print("-" * 70)
            
            for key, engine in md_engines.items():
                print(f"{key}. {engine['name']}")
                print(f"   📝 {engine['description']}")
                print(f"   📋 Requirements: {', '.join(engine['requirements'])}")
                print("-" * 70)
            
            # Get user selection
            while True:
                try:
                    selection = input(f"\n🧬 Select MD engine or 'cancel': ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ MD method creation cancelled")
                        return None
                    
                    if selection in md_engines:
                        selected_engine = md_engines[selection]
                        print(f"\n✅ Selected: {selected_engine['name']}")
                        break
                    else:
                        print("❌ Invalid selection. Please enter a valid option.")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ MD method creation cancelled")
                    return None
            
            ## Set the method method name
            while True:
                try:
                    method_name = input(f"\n📝 Enter method name (or 'cancel'): ").strip()
                    
                    if method_name.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ MD method creation cancelled")
                        return None
                    
                    if not method_name:
                        print("❌ Method name cannot be empty")
                        continue
                    
                    # Validate method name (no special characters) in order to store in DB
                    if not method_name.replace('_', '').replace('-', '').replace(' ', '').isalnum():
                        print("❌ Method name can only contain letters, numbers, spaces, hyphens, and underscores")
                        continue
                    
                    break
                    
                except KeyboardInterrupt:
                    print("\n❌ MD method creation cancelled")
                    return None
            
            # Get optional description
            
            description = input(f"\n📄 Enter method description (optional): ").strip()
            if not description:
                description = f"Molecular dynamics method using {selected_engine['name']}"
            
            
            # Get system preparation options
            prep_params = self._get_system_preparation_parameters()
            


            # Save parameters to database 
            self._save_parameters_to_db(prep_params)

                
        except Exception as e:
            print(f"❌ Error in create_md_method: {e}")
            return None
    
    
    def _get_system_preparation_parameters(self):
        """Get system preparation parameters for MD simulations."""
        
        print(f"\n🧪 CONFIGURE SYSTEM PREPARATION PARAMETERS")
        print("-" * 70)
        
        params = {}
        
        # Query if implicit of explicit solvent simulation is required
        solvent_type = input("Solvent type (implicit/explicit) [default: explicit]: ").strip().lower() or 'explicit'
        
        if solvent_type not in ['implicit', 'explicit']:
            print("❌ Invalid solvent type. Please enter 'implicit' or 'explicit'.")
            return None
        
        params['solvent_type'] = solvent_type
        
        ## Set tleap parameters for system preparation
        tleap_params = self._set_tleap_params(params)
        params['tleap_params'] = tleap_params

        ## Set minimization parameters

        ## Set heating parameters

        ## Set equilibration parameters

        ## Set production parameters



        print(params)

        
            
        print(f"✅ System preparation parameters configured")
        
        return params
    
    def _set_tleap_params(self, params):
        """Set tleap parameters for system preparation."""
        print(f"\n🧪 TLEAP PARAMETERS CONFIGURATION")
        print("-" * 70)
        
        tleap_params = {}

        # Force field
        tleap_params['force_field'] = input("Force field (e.g., ff14SB, ff19SB) [default: ff14SB]: ").strip() or 'ff14SB'
        
        # Set small molecule parameters
        tleap_params['small_molecule_params'] = input("Small molecule parameters (e.g., gaff2, gaff) [default: gaff2]: ").strip() or 'gaff2'
        
        # Water model only for explicit solvent
        if params['solvent_type'] == 'explicit':
            tleap_params['water_model'] = input("Water model (e.g., TIP3P, OPC) [default: TIP3P]: ").strip() or 'TIP3P'
        
        # Ion model only for explicit solvent
        tleap_params['add_ions'] = input("Add ions? (yes/no) [default: yes]: ").strip().lower() or 'yes'
        if tleap_params['add_ions'] in ['yes', 'y']:
            tleap_params['ion_concentration'] = float(input("Ion concentration (M) [default: 0.15]: ").strip() or '0.15')
            tleap_params['positive_ion'] = input("Positive ion (e.g., Na+, K+) [default: Na+]: ").strip() or 'Na+'
            tleap_params['negative_ion'] = input("Negative ion (e.g., Cl-) [default: Cl-]: ").strip() or 'Cl-'
        
        ## Set SolvateBox parameters
        # Box shape
        tleap_params['box_type'] = input("Box type (cubic, dodecahedron, octahedron) [default: dodecahedron]: ").strip() or 'dodecahedron'
        tleap_params['box_distance_nm'] = float(input("Box minimum distance from protein (nm) [default: 10.0]: ").strip() or '10.0')
        tleap_params['wat_closeness'] = float(input("Water closeness parameter (Å) [default: 1.0]: ").strip() or '1.0')
        
        return tleap_params

    def _save_parameters_to_db(self, params):
        """Save system preparation parameters to the MD registers database."""
        try:
            import sqlite3
            import json
            import datetime

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()
            
            # Define dictionary to hold table structure info
            columns_dict = {
                'method_name': 'TEXT',
                'description': 'TEXT',
                'engine': 'TEXT',
                'parameters': 'TEXT',
                'created_date': 'TEXT',
            }

            # Import from moldock.py the _create_table_from_columns_dict helper function
            dbm.create_table_from_columns_dict(cursor, 'md_methods', columns_dict)

            # Import from moldock.py the _update_legacy_table_columns helper function
            dbm.update_legacy_table_columns(cursor, 'md_methods', columns_dict)

            # Import from moldock.py the _remove_legacy_table_columns helper function
            dbm.remove_legacy_table_columns(cursor, 'md_methods', columns_dict)

            # Serialize parameters to JSON string
            params_json = self._serialize_parameters(params)
            
            # Prepare data values dictionary (excluding auto-generated columns)
            data_dict = {
                'method_name': params.get('method_name', 'Unnamed_Method'),
                'description': params.get('description', ''),
                'engine': params.get('engine', ''),
                'parameters': params_json,
                'created_date': datetime.now().isoformat(),
            }

            # Import from moldock.py the _insert_data_into_table helper function
            # _insert_data_into_table is already imported at the top
            dbm.insert_data_dinamically_into_table(cursor, 'md_methods', data_dict)
            
            conn.commit()
            method_id = cursor.lastrowid
            conn.close()
            
            # Display method summary
            self._display_method_summary(method_id, params.get('method_name'), params.get('engine'), 
                                         params.get('description'), params, params_json, 
                                         self.__md_registers_db)
            
        except Exception as e:
            print(f"❌ Error saving parameters to database: {e}")


    def _serialize_parameters(self, parameters):
        """Serialize parameters dictionary to JSON string for database storage."""
        try:
            import json
            return json.dumps(parameters, indent=2)
        except Exception:
            return str(parameters)
    
    def _display_method_summary(self, method_id, method_name, engine_info, 
                                description, parameters, system_prep_params, db_path):
        """Display summary of created MD method."""
        try:
            print(f"\n🎯 MOLECULAR DYNAMICS METHOD CREATED SUCCESSFULLY")
            print("=" * 60)
            print(f"📋 Method ID: {method_id}")
            print(f"🏷️  Method Name: {method_name}")
            print(f"🧬 MD Engine: {engine_info['name']}")
            print(f"📝 Description: {description}")
            print(f"🗂️  Database: {db_path}")
            
            print(f"\n⚙️  ENGINE PARAMETERS:")
            for key, value in parameters.items():
                print(f"   • {key.replace('_', ' ').title()}: {value}")

            print(f"\n🧪 SYSTEM PREPARATION PARAMETERS:")
            if system_prep_params:
                for key, value in system_prep_params.items():
                    print(f"   • {key.replace('_', ' ').title()}: {value}")
            else:
                print("   • None")

            print(f"\n📋 REQUIREMENTS:")
            for req in engine_info['requirements']:
                print(f"   • {req}")

            print(f"\n💡 NEXT STEPS:")
            print("   1. Ensure all requirements are installed")
            print("   2. Use this method for MD simulations")
            print("   3. Prepare input structures for simulation")

            print("=" * 60)

        except Exception as e:
            print(f"⚠️  Error displaying method summary: {e}")