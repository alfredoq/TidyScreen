from tidyscreen import tidyscreen
import sys
import os
from databases import DatabaseManager as dbm
from molecule_management import ligand_management as lm


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
        self.__docking_registers_db = os.path.join(self.path, 'docking/docking_registers', 'docking_assays.db')

        # Set up docking params registers database path within the project directory
        self.__docking_params_db = os.path.join(self.path, 'docking/params', 'params.db')

        # Set up molecular dynamics assays registers database path within the project directory
        self.__md_registers_db = os.path.join(self.path, 'dynamics/md_registers', 'md_registers.db')
        
        # Set up molecular dynamics methods registers database path within the project directory
        self.__md_methods_db = os.path.join(self.path, 'dynamics/md_registers', 'md_methods.db')

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
            
            # Set the params directory for the MD method
            params = {}
            
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
                        params['engine'] = selected_engine['name']
                        
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
                    
                    params['method_name'] = method_name
                    
                    break
                    
                except KeyboardInterrupt:
                    print("\n❌ MD method creation cancelled")
                    return None
            
            # Get optional description
            
            description = input(f"\n📄 Enter method description (optional): ").strip()
            if not description:
                description = f"Molecular dynamics method using {selected_engine['name']}"
            
            params['description'] = description
            
            # Get system preparation options
            params = self._get_system_preparation_parameters(params)

            # Save parameters to database 
            params = self._save_parameters_to_db(params)
            
            # Display method summary
            self._display_method_summary(params)
                
        except Exception as e:
            print(f"❌ Error in create_md_method: {e}")
            return None
    
    def _get_system_preparation_parameters(self, params):
        """Get system preparation parameters for MD simulations."""
        
        print(f"\n🧪 CONFIGURE SYSTEM PREPARATION PARAMETERS")
        print("-" * 70)
        
        # Query if implicit of explicit solvent simulation is required
        solvent_type = input("Solvent type (implicit/explicit) [default: explicit]: ").strip().lower() or 'explicit'
        
        if solvent_type not in ['implicit', 'explicit']:
            print("❌ Invalid solvent type. Please enter 'implicit' or 'explicit'.")
            return None
        
        params['solvent_type'] = solvent_type
        
        ## Set tleap parameters for system preparation
        params = self._set_tleap_params(params)

        ## Set minimization parameters
        params = self._set_minimization_params(params)

        ## Set heating parameters
        params = self._set_heating_params(params)

        ## Set equilibration parameters
        params = self._set_equilibration_params(params)

        ## Set production parameters
        params = self._set_production_params(params)

            
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
            tleap_params['water_model'] = input("Water model (e.g., tip3p, opc, etc) [default: tip3p]: ").strip() or 'tip3p'
        
        if params['solvent_type'] == 'explicit':
            tleap_params['water_type'] = input("Water type (e.g., TIP3PBOX, TIP4PBOX, etc) [default: TIP3PBOX]: ").strip() or 'TIP3PBOX'

        # Ion model only for explicit solvent
        tleap_params['add_ions'] = input("Add ions? (yes/no) [default: yes]: ").strip().lower() or 'yes'
        if tleap_params['add_ions'] in ['yes', 'y']:
            tleap_params['ion_concentration'] = float(input("Ion concentration (M) [default: 0.15]: ").strip() or '0.15')
            tleap_params['positive_ion'] = input("Positive ion (e.g., Na+, K+) [default: Na+]: ").strip() or 'Na+'
            tleap_params['negative_ion'] = input("Negative ion (e.g., Cl-) [default: Cl-]: ").strip() or 'Cl-'
        
        ## Set SolvateBox parameters
        # Box shape
        tleap_params['box_type'] = input("Box type (Box, Oct) [default: Box]: ").strip() or 'Box'
        tleap_params['box_distance_nm'] = float(input("Box minimum distance from protein (nm) [default: 10.0]: ").strip() or '10.0')
        tleap_params['wat_closeness'] = float(input("Water closeness parameter (Å) [default: 1.0]: ").strip() or '1.0')
        
        params['tleap_params'] = tleap_params

        return params

    def _set_minimization_params(self, params):

        # Set parameters for energy minimization step
        print(f"\n🧪 ENERGY MINIMIZATION PARAMETERS CONFIGURATION")
        print("-" * 70)
        min1_params = {}
        min2_params = {}

        # Set parameters for first minimization step (with restraints on solute)
        print(f"\n🔹 First Minimization Step (with restraints on solute)")
        print("-" * 50)
        min1_params['min1_maxcyc'] = int(input("Maximum number of minimization cycles [default: 5000]: ").strip() or '5000')
        min1_params['min1_ncyc'] = int(input("Number of cycles for steepest descent [default: 2500]: ").strip() or '2500')
        if 'general_params' not in params:
            params['general_params'] = {}
        params['general_params']['cutoff'] = int(input("Cutoff distance for nonbonded interactions (Å) [default: 10]: ").strip() or '10')
        min1_params['min1_restraint_selector'] = input("Restraint selector for first minimization (e.g., '(@C,N,CA,O)' ; '(!:WAT & !@Na+)' ) [default: '(!:WAT & !@Na+ & !@Cl-)']: ").strip() or '(!:WAT & !@Na+ & !@Cl-)'
        min1_params['min1_restraint_wt'] = float(input("Restraint weight on solute (kcal/mol·Å²) [default: 10.0]: ").strip() or '10.0')
        
        # Assign first minimization parameters to min_params
        params['first_minimization'] = min1_params
        
        # Set parameters for second minimization step (without restraints)
        print(f"\n🔹 Second Minimization Step (without restraints)")
        print("-" * 50)
        min2_params['min2_maxcyc'] = int(input("Maximum number of minimization cycles [default: 5000]: ").strip() or '5000')
        min2_params['min2_ncyc'] = int(input("Number of cycles for steepest descent [default: 2500]: ").strip() or '2500')
        min2_params['min2_restraint_selector'] = input("Restraint selector for second minimization (e.g., '@C,N,CA,O' ; '!@WAT & !@Na+' ) [default: '@C,N,CA,O']: ").strip() or '@C,N,CA,O'
        min2_params['min2_restraint_wt'] = float(input("Restraint weight on solute (kcal/mol·Å²) [default: 1.0]: ").strip() or '1.0')
        
        # Assign second minimization parameters to min_params
        params['second_minimization'] = min2_params

        return params

    def _set_heating_params(self, params):

        # Set the parameters for the heating stage
        print(f"\n🧪 HEATING PARAMETERS CONFIGURATION")
        print("-" * 70)
        heating_params = {}
        heating_params['initial_temp'] = float(input("Initial temperature for heating (K) [default: 0]: ").strip() or '0')
        params['general_params']['target_temp'] = float(input("Target temperature for heating (K) [default: 300]: ").strip() or '300')
        params['general_params']['collision_freq'] = float(input("Collision frequency (ps^-1) [default: 1]: ").strip() or '1')
        params['general_params']['thermostat'] = int(input("Thermostat selection [default: 3]: ").strip() or '3')
        heating_params['heating_steps'] = int(input("Number of steps for heating [default: 50000]: ").strip() or '50000')
        heating_params['heating_restraint_selector'] = input("Restraint selector during heating (e.g., '(:* & !:WAT & !:Na+ & !:Cl-)', '@C,N,CA,O' ; '(!:WAT & !@Na+)' ) [default: '(:* & !:WAT & !:Na+ & !:Cl-)']: ").strip() or '(:* & !:WAT & !:Na+ & !:Cl-)'
        heating_params['heating_restraint_wt'] = float(input("Restraint weight during heating (kcal/mol·Å²) [default: 0.5]: ").strip() or '0.5')
        heating_params['heating_restart_write'] = int(input("Number of frame to write restart [default: 1000]: ").strip() or '1000')
        heating_params['heating_trajectory_write'] = int(input("Number of frame to write to trajectory [default: 1000]: ").strip() or '1000')
        heating_params['heating_output_write'] = int(input("Number of frame to write to output [default: 1000]: ").strip() or '1000')
        
        params['heating_params'] = heating_params
        
        return params

    def _set_equilibration_params(self, params):

        # Set the parameters for the equilibration stage
        print(f"\n🧪 EQUILIBRATION PARAMETERS CONFIGURATION")
        print("-" * 70)
        equilibration_params = {}
        equilibration_params['equilibration_steps'] = int(input("Number of steps for equilibration [default: 500000]: ").strip() or '500000')
        equilibration_params['equilibration_restraint_selector'] = input("Restraint selector during equilibration (e.g., '(:* & !:WAT & !:Na+ & !:Cl-)', '@C,N,CA,O' ; '(!:WAT & !@Na+)' ) [default: '@C,N,CA,O']: ").strip() or '@C,N,CA,O'
        equilibration_params['equilibration_restraint_wt'] = float(input("Restraint weight during equilibration (kcal/mol·Å²) [default: 0.5]: ").strip() or '0.5')
        equilibration_params['equilibration_restart_write'] = int(input("Number of frame to write restart [default: 1000]: ").strip() or '1000')
        equilibration_params['equilibration_trajectory_write'] = int(input("Number of frame to write to trajectory [default: 1000]: ").strip() or '1000')
        equilibration_params['equilibration_output_write'] = int(input("Number of frame to write to output [default: 1000]: ").strip() or '1000')
        
        params['equilibration_params'] = equilibration_params
        
        return params

    def _set_production_params(self, params):

        # Set the parameters for the production stage
        print(f"\n🧪 PRODUCTION PARAMETERS CONFIGURATION")
        print("-" * 70)
        production_params = {}
        production_params['production_steps'] = int(input("Number of steps for production [default: 5000000]: ").strip() or '5000000')
        production_params['production_restraint_selector'] = input("Restraint selector during production (e.g., '(:* & !:WAT & !:Na+ & !:Cl-)', '@C,N,CA,O' ; '(!:WAT & !@Na+)' ) [default: '@C,N,CA,O']: ").strip() or '@C,N,CA,O'
        production_params['production_restraint_wt'] = float(input("Restraint weight during production (kcal/mol·Å²) [default: 0.5]: ").strip() or '0.5')
        production_params['production_restart_write'] = int(input("Number of frame to write restart [default: 1000]: ").strip() or '1000')
        production_params['production_trajectory_write'] = int(input("Number of frame to write to trajectory [default: 1000]: ").strip() or '1000')
        production_params['production_output_write'] = int(input("Number of frame to write to output [default: 1000]: ").strip() or '1000')
        
        params['production_params'] = production_params
        
        return params

    def _save_parameters_to_db(self, params):
        """Save system preparation parameters to the MD registers database."""
        try:
            import sqlite3

            conn = sqlite3.connect(self.__md_methods_db)
            cursor = conn.cursor()
            
            # Define dictionary to hold table structure info
            columns_dict = {
                'method_id': 'INTEGER PRIMARY KEY AUTOINCREMENT',
                'method_name': 'TEXT',
                'description': 'TEXT',
                'engine': 'TEXT',
                'parameters': 'TEXT',
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
            }

            # Insert dinamically data into md_methods table
            dbm.insert_data_dinamically_into_table(cursor, 'md_methods', data_dict)
            
            conn.commit()
            method_id = cursor.lastrowid
            
            params['method_id'] = method_id
            
            conn.close()
            
            return params
           
        except Exception as e:
            print(f"❌ Error saving parameters to database: {e}")

    def _serialize_parameters(self, parameters):
        """Serialize parameters dictionary to JSON string for database storage."""
        try:
            import json
            return json.dumps(parameters, indent=2)
        except Exception:
            return str(parameters)
    
    def _display_method_summary(self, params):
        """Display summary of created MD method."""
        try:
            print(f"\n🎯 MOLECULAR DYNAMICS METHOD CREATED SUCCESSFULLY")
            print("=" * 60)
            print(f"\n🏷️  Method ID: {params.get('method_id', 'N/A')}")
            print(f"🏷️  Method Name: {params.get('method_name', 'Unnamed_Method')}")
            print(f"🧬 MD Engine: {params.get('engine', '')}")
            print(f"📝 Description: {params.get('description', '')}")
            
            print(f"\n⚙️  MD SPECIFIC PARAMETERS:")
            
            tleap_params = params.get('tleap_params', {})
            print(f"  🔬 TLEAP PARAMETERS:")
            for key, value in tleap_params.items():
                print(f"    - {key}: {value}")

            min1_params = params.get('minimization_params', {}).get('first_minimization', {})
            min2_params = params.get('minimization_params', {}).get('second_minimization', {})
            print(f"  🧪 MINIMIZATION PARAMETERS:")
            # Print first minimization step parameters
            print(f"    🔹 First Minimization Step:")
            for key, value in min1_params.items():
                print(f"      - {key}: {value}")
            # Print second minimization step parameters
            print(f"    🔹 Second Minimization Step:")
            for key, value in min2_params.items():
                print(f"      - {key}: {value}")

            # Print heating parameters
            heating_params = params.get('heating_params', {})
            print(f"  🔥 HEATING PARAMETERS:")
            for key, value in heating_params.items():
                print(f"    - {key}: {value}")

            # Print equilibration parameters
            equilibration_params = params.get('equilibration_params', {})
            print(f"  ❄️  EQUILIBRATION PARAMETERS:")
            for key, value in equilibration_params.items():
                print(f"    - {key}: {value}")

            # Print production parameters   
            production_params = params.get('production_params', {})
            print(f"  🎬 PRODUCTION PARAMETERS:")
            for key, value in production_params.items():
                print(f"    - {key}: {value}")

            print("=" * 60)

        except Exception as e:
            print(f"⚠️  Error displaying method summary: {e}")
    
    def list_md_methods(self):
        """List all molecular dynamics methods registered in the project and allow user to view details."""
        try:
            import sqlite3
            import json

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()
            
            # Check if md_methods table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_methods';")
            if not cursor.fetchone():
                print("❌ No MD methods found. The 'md_methods' table does not exist.")
                return
            
            # Fetch all methods
            cursor.execute("SELECT method_id, method_name, engine, description FROM md_methods;")
            methods = cursor.fetchall()
            
            if not methods:
                print("❌ No MD methods found in the database.")
                return
            
            print(f"\n🧬 MOLECULAR DYNAMICS METHODS IN PROJECT '{self.name}':")
            print("=" * 70)
            for method in methods:
                method_id, method_name, engine, description = method
                print(f"🏷️  Method ID: {method_id}")
                print(f"🏷️  Method Name: {method_name}")
                print(f"🧬 MD Engine: {engine}")
                print(f"📝 Description: {description}")
                print("-" * 70)

            # Prompt user to select a method to show details
            method_ids = [str(m[0]) for m in methods]
            while True:
                selection = input("\n🔎 Enter the Method ID to view details (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    conn.close()
                    return
                if selection in method_ids:
                    # Fetch full parameters for the selected method
                    cursor.execute("SELECT parameters FROM md_methods WHERE method_id = ?", (selection,))
                    row = cursor.fetchone()
                    if row:
                        try:
                            params = json.loads(row[0])
                        except Exception:
                            params = {}
                        params['method_id'] = int(selection)
                        self._display_method_summary(params)
                    else:
                        print("❌ Method not found.")
                    break
                else:
                    print("❌ Invalid Method ID. Please try again.")

            conn.close()
            
        except Exception as e:
            print(f"❌ Error listing MD methods: {e}")
    
    def delete_md_method(self):
        """Delete a molecular dynamics method from the project database."""
        try:
            import sqlite3

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()
            
            # Check if md_methods table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_methods';")
            if not cursor.fetchone():
                print("❌ No MD methods found. The 'md_methods' table does not exist.")
                return
            
            # Fetch all methods
            cursor.execute("SELECT method_id, method_name FROM md_methods;")
            methods = cursor.fetchall()
            
            if not methods:
                print("❌ No MD methods found in the database.")
                return
            
            print(f"\n🧬 MOLECULAR DYNAMICS METHODS IN PROJECT '{self.name}':")
            print("=" * 70)
            for method in methods:
                method_id, method_name = method
                print(f"🏷️  Method ID: {method_id}")
                print(f"🏷️  Method Name: {method_name}")
                print("-" * 70)

            # Prompt user to select a method to delete
            method_ids = [str(m[0]) for m in methods]
            while True:
                selection = input("\n🗑️  Enter the Method ID to delete (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    conn.close()
                    return
                if selection in method_ids:
                    # Confirm deletion
                    confirm = input(f"⚠️  Are you sure you want to delete Method ID {selection}? (yes/no): ").strip().lower()
                    if confirm in ['yes', 'y']:
                        cursor.execute("DELETE FROM md_methods WHERE method_id = ?", (selection,))
                        conn.commit()
                        print(f"✅ Method ID {selection} deleted successfully.")
                    else:
                        print("❌ Deletion cancelled.")
                    break
                else:
                    print("❌ Invalid Method ID. Please try again.")

            conn.close()
            
        except Exception as e:
            print(f"❌ Error deleting MD method: {e}")
            
    def perform_md_assay(self):
        """
        Will connect to the docking registers database and list available docking assays.
        """
        
        try:
            import sqlite3

            # Connect to docking registers database
            conn = sqlite3.connect(self.__docking_registers_db)
            cursor = conn.cursor()
            # Show docking assays available
            docking_assay_params_dict = self._select_docking_assay(cursor)
            conn.close()

            # Select unique ligand molecules in the selected docking assay
            docking_assay_params_dict = self._select_unique_ligands_in_docking_assay(docking_assay_params_dict)

            # Select pose id for the selected ligand to perform MD assay
            docking_assay_params_dict = self._select_ligand_pose_for_md_assay(docking_assay_params_dict)

            # Select a MD method to perform the MD assay
            md_parameters_dict = self._select_md_method()

            # Restore the selected docked pose using ligand_management module
            pdb_dict = lm.restore_single_docked_pose(docking_assay_params_dict.get('docking_results_db'), docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict.get('selected_pose_id'))

            # Query user for MD assay description
            assay_description = input("\n📝 Enter MD assay description (optional): ").strip()
            if not assay_description:
                assay_description = f"MD assay for ligand {docking_assay_params_dict.get('selected_ligand_name')} (pose {docking_assay_params_dict.get('selected_pose_id')})"

            # Create MD assay folder structure and prepare input files for the selected MD engine
            md_assay_id, md_assay_folder = self._create_md_assay_folder()

            # Save the docked pose PDB file in the MD assay folder
            selected_pose_pdb = lm.write_pdb_with_moldf(pdb_dict, docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict.get('selected_pose_id'), md_assay_folder)

            # Prepare the ligand tleap input file in the MD assay folder
            mol2_file, frcmod_file = lm.prepare_ligand_tleap_input_files(self.__chemspace_db, docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict, md_assay_folder)
            
            # Create complex .prmtop and .inpcrd files
            prmtop_file, inpcrd_file = self._prepare_complex_prmtop_inpcrd_for_md(mol2_file, frcmod_file, docking_assay_params_dict, md_assay_folder, md_assay_folder, md_parameters_dict, selected_pose_pdb)
            
            # Prepare md simulation input files
            self._prepare_md_simulation_input_files(md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file)
            
            # Prepare execution scripts for the MD assay
            self._prepare_md_execution_script(md_assay_folder)
            
            # Create MD register entry in the md_assays table
            self._create_md_assay_register_entry(md_assay_id, md_assay_folder, assay_description, docking_assay_params_dict.get('assay_id'), docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict.get('selected_pose_id'), md_parameters_dict)
            
            # Query if the user want to start the MD simulation now
            start_now = input("\n▶️  Do you want to start the MD simulation now? (yes/no) [default: no]: ").strip().lower() or 'no'
            if start_now in ['yes', 'y']:
                self._start_md_simulation(md_assay_folder)
            else:
                print("❌ MD simulation not started. You can start it later by running the execution script in the assay folder.")
                

        except Exception as e:
            print(f"❌ Error performing MD assay: {e}")
        
    def _select_docking_assay(self, cursor):
        
        # Check if docking_assays table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='docking_assays';")
            if not cursor.fetchone():
                print("❌ No docking assays found. The 'docking_assays' table does not exist.")
                return

            # Fetch all docking assays
            cursor.execute("SELECT assay_id, assay_name, table_name, receptor_info, notes, assay_folder_path FROM docking_assays;")
            assays = cursor.fetchall()

            if not assays:
                print("❌ No docking assays found in the database.")
                return

            print(f"\n🧬 DOCKING ASSAYS IN PROJECT '{self.name}':")
            print("=" * 70)
            for assay in assays:
                assay_id, assay_name, _, _, notes, _ = assay
                print(f"🏷️  Assay ID: {assay_id}")
                print(f"🏷️  Assay Name: {assay_name}")
                print(f"📝 Description: {notes}")
                print("-" * 70)

            # Prompt user to select an assay for further action (optional, extend as needed)
            assay_ids = [str(a[0]) for a in assays]
            while True:
                selection = input("\n🔎 Enter the Assay ID to proceed (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    return
                if selection in assay_ids:
                    print(f"✅ Selected Assay ID: {selection}")
                    
                    # Find the selected assay tuple
                    selected_assay = next((a for a in assays if str(a[0]) == selection), None)
                    if selected_assay:
                        # Get column names from cursor description
                        column_names = [desc[0] for desc in cursor.description]
                        # Create dictionary mapping column names to values
                        docking_assay_params_dict = dict(zip(column_names, selected_assay))
                        
                    else:
                        print("❌ Could not retrieve selected assay details.")
                    
                    
                    # Further actions can be implemented here
                    return docking_assay_params_dict
                else:
                    print("❌ Invalid Assay ID. Please try again.")
    
    def _select_unique_ligands_in_docking_assay(self, docking_assay_params_dict):
        try:
            import sqlite3

            # Get required parameters from docking_assay_params_dict
            assay_folder_path = docking_assay_params_dict.get('assay_folder_path')
            assay_name = docking_assay_params_dict.get('assay_name')
            # Connect to docking results database within the assay folder
            docking_results_db = os.path.join(assay_folder_path, 'results', f'{assay_name}.db')
            conn = sqlite3.connect(docking_results_db)
            cursor = conn.cursor()

            print(f"\n🧬 UNIQUE LIGANDS IN DOCKING ASSAY '{assay_name}':")
            print("=" * 70)
        # Fetch unique ligands in the selected docking assay
            cursor.execute("""
                SELECT DISTINCT LigName
                FROM Results;
            """)
            ligands = cursor.fetchall()

            if not ligands:
                print("❌ No ligands found for the selected docking assay.")
                return
            
            for idx, ligand in enumerate(ligands, start=1):
                lig_name = ligand[0]
                print(f"🏷️  Ligand ID: {idx}")
                print(f"🏷️  Ligand Name: {lig_name}")
                print("-" * 70)
            conn.close()
            
            # Prompt the used to select a ligand id and return the ligand name
            ligand_ids = [str(i) for i in range(1, len(ligands) + 1)]
            while True:
                selection = input("\n🔎 Enter the Ligand ID to proceed (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    return
                if selection in ligand_ids:
                    selected_ligand_name = ligands[int(selection) - 1][0]
                    print(f"✅ Selected Ligand: {selected_ligand_name}")
                    # Add relevan info to docking_assay_params_dict
                    docking_assay_params_dict['selected_ligand_name'] = selected_ligand_name
                    docking_assay_params_dict['docking_results_db'] = docking_results_db
                    
                    return docking_assay_params_dict
                else:
                    print("❌ Invalid Ligand ID. Please try again.")
            
        except Exception as e:
            print(f"❌ Error fetching unique ligands: {e}")
    
    def _select_ligand_pose_for_md_assay(self, docking_assay_params_dict):
        try:
            import sqlite3

            docking_results_db = docking_assay_params_dict.get('docking_results_db')
            assay_name = docking_assay_params_dict.get('assay_name')
            ligname = docking_assay_params_dict.get('selected_ligand_name')

            # Connect to docking results database within the assay folder
            conn = sqlite3.connect(docking_results_db)
            cursor = conn.cursor()

            print(f"\n🧬 POSES FOR LIGAND '{ligname}':")
            print("=" * 70)
        # Fetch poses for the selected ligand
            cursor.execute("""
                SELECT Pose_ID, LigName, pose_rank, docking_score, cluster_size
                FROM Results
                WHERE LigName = ?;
            """, (ligname,))
            poses = cursor.fetchall()

            if not poses:
                print("❌ No poses found for the selected ligand.")
                return
            
            for pose in poses:
                pose_id, ligname, pose_rank, docking_score, cluster_size = pose
                print(f"🏷️  Pose ID: {pose_id}")
                print(f"🏷️  Ligand Name: {ligname}")
                print(f"🏷️  Pose Rank: {pose_rank}")
                print(f"🏷️  Docking Score: {docking_score}")
                print(f"🏷️  Cluster Size: {cluster_size}")
                print("-" * 70)
            conn.close()
            
            # Prompt the user to select a pose id for MD assay
            pose_ids = [str(p[0]) for p in poses]
            while True:
                selection = input("\n🔎 Enter the Pose ID to proceed with MD assay (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    return
                if selection in pose_ids:
                    print(f"✅ Selected Pose ID: {selection} for MD assay")
                    # Add the selected pose id to docking_assay_params_dict
                    docking_assay_params_dict['selected_pose_id'] = int(selection)
                    
                    return docking_assay_params_dict
                else:
                    print("❌ Invalid Pose ID. Please try again.")
            
        except Exception as e:
            print(f"❌ Error fetching ligand poses: {e}")

    def _create_md_assay_folder(self):

        """Create molecular dynamics assay folder structure and prepare input files."""
        
        last_assay_id = self._get_last_md_assay_id()

        new_assay_id = last_assay_id + 1
        new_assay_folder = os.path.join(self.__md_assays_folder, f'md_assay_{new_assay_id}')
        try:
            os.makedirs(new_assay_folder, exist_ok=False)
            print(f"✅ Created MD assay folder: {new_assay_folder}")
            # Further input file preparation can be implemented here
            return new_assay_id, new_assay_folder
        except Exception as e:
            print(f"❌ Error creating MD assay folder: {e}")

    def _get_last_md_assay_id(self):
        """Retrieve the last MD assay ID from the MD registers database."""
        try:
            import sqlite3

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()

            # Check if md_assays table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_assays';")
            if not cursor.fetchone():
                print("❌ No MD assays found. The 'md_assays' table does not exist.")
                return 0

            # Fetch the last assay ID
            cursor.execute("SELECT MAX(assay_id) FROM md_assays;")
            row = cursor.fetchone()
            last_assay_id = row[0] if row and row[0] is not None else 0

            conn.close()
            return last_assay_id
        except Exception as e:
            print(f"❌ Error retrieving last MD assay ID: {e}")
            sys.exit(1)

    def _create_md_assay_register_entry(self, new_assay_id, new_assay_folder, assay_description, docking_assay_id=None, ligname=None, pose_id=None, md_parameters_dict=None):
        """Create a new entry in the md_assays table for the newly created MD assay."""
        try:
            import sqlite3

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()

            # Define dictionary to hold table structure info
            columns_dict = {
                'assay_id': 'INTEGER PRIMARY KEY AUTOINCREMENT',
                'md_assay': 'TEXT',
                'description': 'TEXT',
                'assay_folder_path': 'TEXT',
                'docking_assay_id': 'INTEGER',
                'ligand_name': 'TEXT',
                'pose_id': 'INTEGER',
                'md_parameters': 'TEXT',
            }
            # Create md_assays table if it doesn't exist
            dbm.create_table_from_columns_dict(cursor, 'md_assays', columns_dict, verbose=False)

            # Update legacy table columns 
            dbm.update_legacy_table_columns(cursor, 'md_assays', columns_dict, verbose=False)

            # Remove legacy table columns
            dbm.remove_legacy_table_columns(cursor, 'md_assays', columns_dict, verbose=False)

            # Prepare data values dictionary (excluding auto-generated columns)
            data_dict = {
                'md_assay': f'assay_{new_assay_id}',
                'description': assay_description,
                'assay_folder_path': new_assay_folder,
                'docking_assay_id': f'assay_{docking_assay_id}',
                'ligand_name': ligname,
                'pose_id': pose_id,
                'md_parameters': self._serialize_parameters(md_parameters_dict),
            }
            # Insert dinamically data into md_assays table
            dbm.insert_data_dinamically_into_table(cursor, 'md_assays', data_dict)
            conn.commit()
            conn.close()

        except Exception as e:
            print(f"❌ Error creating MD assay register entry: {e}")
            
    def _prepare_complex_prmtop_inpcrd_for_md(self, mol2_file, frcmod_file, assay_info, output_dir, output_file, md_parameters_dict, pdb_file):
        
        import subprocess
        import json

        # Parse receptor_info as a dictionary if it's a JSON string or already a dict
        receptor_info = assay_info.get('receptor_info', None)
        if isinstance(receptor_info, str):
            try:
                receptor_details = json.loads(receptor_info)
            except Exception:
                receptor_details = {}
        elif isinstance(receptor_info, dict):
            receptor_details = receptor_info
        else:
            receptor_details = {}

        # get the .pdbqt used for docking
        receptor_pdb = receptor_details.get('pdbqt_file', None)
        
        # Define the raw .pdb file of the receptor
        receptor_pdb_path = '/'.join(receptor_pdb.split('/')[:-1]) + '/receptor_checked.pdb'

        # # Create tleap input file to load receptor and ligand
        tleap_in_file = os.path.join(output_dir, "complex.in")

        # Defined output files
        prmtop_file = os.path.join(output_dir, 'complex.prmtop')
        inpcrd_file = os.path.join(output_dir, 'complex.inpcrd')

        tleap_params = md_parameters_dict.get('tleap_params', {})
        
        protein_ff = tleap_params.get('force_field', 'ff14SB')
        ligand_ff = tleap_params.get('small_molecule_params', 'gaff2')
        solvent_type = md_parameters_dict.get('solvent_type', 'explicit')
        if solvent_type == 'explicit':
            water_model = tleap_params.get('water_model', 'tip3p')
            water_type = tleap_params.get('water_type', 'TIP3PBOX')
            add_ions = tleap_params.get('add_ions', 'yes')
            ion_concentration = tleap_params.get('ion_concentration', 0.15)
            positive_ion = tleap_params.get('positive_ion', 'Na+')
            negative_ion = tleap_params.get('negative_ion', 'Cl-')
            box_type = tleap_params.get('box_type', 'Box')
            box_distance_nm = tleap_params.get('box_distance_nm', 1.0)
            wat_closeness = tleap_params.get('wat_closeness', 1.0)
        
        with open(tleap_in_file, 'w') as f:
            f.write(f"source leaprc.protein.{protein_ff}\n")
            f.write(f"source leaprc.{ligand_ff}\n")
            if solvent_type == 'explicit':
                f.write(f"source leaprc.water.{water_model}\n")
                f.write(f"HOH = WAT\n")
            f.write(f"rec = loadpdb {receptor_pdb_path}\n")
            f.write(f"UNL = loadmol2 {mol2_file}\n")
            f.write(f"loadamberparams {frcmod_file}\n")
            f.write(f"lig = loadpdb {pdb_file}\n")
            f.write(f"COM = combine {{rec lig}}\n")
            if solvent_type == 'explicit':
                f.write(f"solvate{box_type} COM {water_type} {box_distance_nm} {wat_closeness}\n")
                if add_ions in ['yes', 'y']:
                    f.write(f"addions COM {positive_ion} 0\n")
                    f.write(f"addions COM {negative_ion} 0\n")
            f.write(f"saveamberparm COM {prmtop_file} {inpcrd_file}\n")
            f.write("quit\n")

        # Run tleap to generate prmtop and inpcrd files
        tleap_command = f"tleap -f {tleap_in_file}"
        try:
            subprocess.run(tleap_command, shell=True, check=True, stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

        except Exception as e:
            print(f"[ERROR] Failed to run tleap for complex generation: {e}")
            sys.exit(1)

        return prmtop_file, inpcrd_file
        
    def _select_md_method(self):
        """Select a molecular dynamics method from the registered methods."""
        try:
            import sqlite3
            import json

            conn = sqlite3.connect(self.__md_methods_db)
            cursor = conn.cursor()
            
            # Check if md_methods table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_methods';")
            if not cursor.fetchone():
                print("❌ No MD methods found. The 'md_methods' table does not exist.")
                return
            
            # Fetch all methods
            cursor.execute("SELECT method_id, method_name, engine, description, parameters FROM md_methods;")
            methods = cursor.fetchall()
            
            if not methods:
                print("❌ No MD methods found in the database.")
                return
            
            print(f"\n🧬 MOLECULAR DYNAMICS METHODS IN PROJECT '{self.name}':")
            print("=" * 70)
            for method in methods:
                method_id, method_name, engine, description, _ = method
                print(f"🏷️  Method ID: {method_id}")
                print(f"🏷️  Method Name: {method_name}")
                print(f"🧬 MD Engine: {engine}")
                print(f"📝 Description: {description}")
                print("-" * 70)

            # Prompt user to select a method to show details
            method_ids = [str(m[0]) for m in methods]
            while True:
                selection = input("\n🔎 Enter the Method ID to select for MD assay (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    conn.close()
                    return
                if selection in method_ids:
                    print(f"✅ Selected Method ID: {selection} for MD assay")
                    
                    # Fetch the parameters column for the selected method
                    cursor.execute("SELECT parameters FROM md_methods WHERE method_id = ?", (selection,))
                    row = cursor.fetchone()
                    
                    if row:
                        try:
                            parameters_dict = json.loads(row[0])
                        except Exception:
                            parameters_dict = {}
                        print(f"\n📋 Parameters for selected MD method (as dictionary):")
                        return parameters_dict
                        
                    else:
                        print("❌ Could not retrieve parameters for the selected method.")
                    
                    conn.close()
                    return int(selection)
                else:
                    print("❌ Invalid Method ID. Please try again.")

        except Exception as e:
            print(f"❌ Error selecting MD method: {e}")

    def _prepare_md_simulation_input_files(self, md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file):

        """Prepare MD simulation input files based on the selected MD method and engine."""
        try:
            import os

            engine = md_parameters_dict.get('engine', 'AMBER')
            solvent_type = md_parameters_dict.get('solvent_type', 'explicit')

            if solvent_type == 'explicit':
                self._prepare_explicit_solvent_md_input_files(md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file)

            else:
                self._prepare_implicit_solvent_md_input_files(engine, md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file)

        except Exception as e:
            print(f"❌ Error preparing MD simulation input files: {e}")

    def _prepare_explicit_solvent_md_input_files(self, md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file):

        ## Prepare Min1 input file
        # Create min1 input file
        min1_in_file = os.path.join(md_assay_folder, "min1.in")

        with open(min1_in_file, 'w') as f:
            f.write(f"Minimization 1 stage\n")
            f.write(f"&cntrl\n")
            f.write(f"  imin=1,\n")
            f.write(f"  maxcyc={md_parameters_dict.get('first_minimization').get('min1_maxcyc')},\n")
            f.write(f"  ncyc={md_parameters_dict.get('first_minimization').get('min1_ncyc')},\n")
            f.write(f"  ntb=1,\n")
            f.write(f"  ntc=2,\n")
            f.write(f"  ntf=2,\n")
            f.write(f"  cut={md_parameters_dict.get('general_params').get('cutoff')},\n")
            f.write(f"  ntpr=100,\n")
            f.write(f"  iwrap=1,\n")
            f.write(f"  ntr=1,\n")
            f.write(f"  restraint_wt={md_parameters_dict.get('first_minimization').get('min1_restraint_wt')},\n")
            f.write(f"  restraintmask='{md_parameters_dict.get('first_minimization').get('min1_restraint_selector')}',\n")
            f.write(f"/\n")
            f.write(f"END\n")
            
        
        # ## Prepare Min2 input file
        min2_in_file = os.path.join(md_assay_folder, "min2.in")

        with open(min2_in_file, 'w') as f:
            f.write(f"Minimization 2 stage\n")
            f.write(f"&cntrl\n")
            f.write(f"  imin=1,\n")
            f.write(f"  maxcyc={md_parameters_dict.get('second_minimization').get('min2_maxcyc')},\n")
            f.write(f"  ncyc={md_parameters_dict.get('second_minimization').get('min2_ncyc')},\n")
            f.write(f"  ntb=1,\n")
            f.write(f"  ntc=2,\n")
            f.write(f"  ntf=2,\n")
            f.write(f"  cut={md_parameters_dict.get('general_params').get('cutoff')},\n")
            f.write(f"  ntpr=100,\n")
            f.write(f"  iwrap=1,\n")
            f.write(f"  ntr=1,\n")
            f.write(f"  restraint_wt={md_parameters_dict.get('second_minimization').get('min2_restraint_wt')},\n")
            f.write(f"  restraintmask='{md_parameters_dict.get('second_minimization').get('min2_restraint_selector')}',\n")
            f.write(f"/\n")
            f.write(f"END\n")

        # Prepare heating input file
        heating_in_file = os.path.join(md_assay_folder, "heating.in")
        
        with open(heating_in_file, 'w') as f:
            f.write(f"Heating stage\n")
            f.write(f"&cntrl\n")
            f.write(f"  imin=0,\n")
            f.write(f"  irest=0,\n")
            f.write(f"  ntx=1,\n")
            f.write(f"  ntb=2,\n")
            f.write(f"  ntp=1,\n")
            f.write(f"  cut={md_parameters_dict.get('general_params').get('cutoff')},\n")
            f.write(f"  ntr=1,\n")
            f.write(f"  ntc=2,\n")
            f.write(f"  ntf=2,\n")
            f.write(f"  temp0={md_parameters_dict.get('general_params').get('target_temp')},\n")
            f.write(f"  tempi={md_parameters_dict.get('heating_params').get('initial_temp')},\n")
            f.write(f"  ntt={md_parameters_dict.get('general_params').get('thermostat')},\n")
            f.write(f"  gamma_ln={md_parameters_dict.get('general_params').get('collision_freq')},\n")
            f.write(f"  nstlim={md_parameters_dict.get('heating_params').get('heating_steps')},\n")
            f.write(f"  dt=0.002,\n")
            f.write(f"  ntpr={md_parameters_dict.get('heating_params').get('heating_output_write')},\n")
            f.write(f"  ntwx={md_parameters_dict.get('heating_params').get('heating_trajectory_write')},\n")
            f.write(f"  ntwr={md_parameters_dict.get('heating_params').get('heating_restart_write')},\n")
            f.write(f"  restraint_wt={md_parameters_dict.get('heating_params').get('heating_restraint_wt')},\n")
            f.write(f"  restraintmask='{md_parameters_dict.get('heating_params').get('heating_restraint_selector')}',\n")
            f.write(f"  iwrap=1,\n")
            f.write(f"  nmropt=1,\n")
            f.write(f"&end\n")
            f.write(f"&wt\n")
            f.write(f"  TYPE = 'TEMP0',\n")
            f.write(f"  ISTEP1 = 0,\n")
            f.write(f"  ISTEP2 = {md_parameters_dict.get('heating_params').get('heating_steps') // 2},\n")
            f.write(f"  VALUE1 = 0,\n")
            f.write(f"  VALUE2 = {md_parameters_dict.get('general_params').get('target_temp')},\n")
            f.write(f"/\n")
            f.write(f"&wt TYPE = 'END'\n")
            f.write(f"/\n")
            f.write(f"END\n")

        # Prepare the equilibration input file
        equilibration_in_file = os.path.join(md_assay_folder, "equilibration.in")
        with open(equilibration_in_file, 'w') as f:
            f.write(f"Equilibration stage\n")
            f.write(f"&cntrl\n")
            f.write(f"  imin=0,\n")
            f.write(f"  irest=1,\n")
            f.write(f"  ntx=5,\n")
            f.write(f"  ntb=1,\n")
            f.write(f"  ntp=0,\n")
            f.write(f"  cut={md_parameters_dict.get('general_params').get('cutoff')},\n")
            f.write(f"  ntr=1,\n")
            f.write(f"  ntc=2,\n")
            f.write(f"  ntf=2,\n")
            f.write(f"  temp0={md_parameters_dict.get('general_params').get('target_temp')},\n")
            f.write(f"  tempi={md_parameters_dict.get('general_params').get('target_temp')},\n")
            f.write(f"  ntt={md_parameters_dict.get('general_params').get('thermostat')},\n")
            f.write(f"  gamma_ln={md_parameters_dict.get('general_params').get('collision_freq')},\n")
            f.write(f"  nstlim={md_parameters_dict.get('equilibration_params').get('equilibration_steps')},\n")
            f.write(f"  dt=0.002,\n")
            f.write(f"  ntpr={md_parameters_dict.get('equilibration_params').get('equilibration_output_write')},\n")
            f.write(f"  ntwx={md_parameters_dict.get('equilibration_params').get('equilibration_trajectory_write')},\n")
            f.write(f"  ntwr={md_parameters_dict.get('equilibration_params').get('equilibration_restart_write')},\n")
            f.write(f"  iwrap=1,\n")
            f.write(f"/\n")
            f.write(f"END\n")

        # Prepare the production MD input file
        production_in_file = os.path.join(md_assay_folder, "production.in")
        with open(production_in_file, 'w') as f:
            f.write(f"Production MD stage\n")
            f.write(f"&cntrl\n")
            f.write(f"  imin=0,\n")
            f.write(f"  irest=1,\n")
            f.write(f"  ntx=5,\n")
            f.write(f"  ntb=1,\n")
            f.write(f"  ntp=0,\n")
            f.write(f"  cut={md_parameters_dict.get('general_params').get('cutoff')},\n")
            f.write(f"  ntc=2,\n")
            f.write(f"  ntf=2,\n")
            f.write(f"  temp0={md_parameters_dict.get('general_params').get('target_temp')},\n")
            f.write(f"  ntt={md_parameters_dict.get('general_params').get('thermostat')},\n")
            f.write(f"  gamma_ln={md_parameters_dict.get('general_params').get('collision_freq')},\n")
            f.write(f"  nstlim={md_parameters_dict.get('production_params').get('production_steps')},\n")
            f.write(f"  dt=0.002,\n")
            f.write(f"  ntpr={md_parameters_dict.get('production_params').get('production_output_write')},\n")
            f.write(f"  ntwx={md_parameters_dict.get('production_params').get('production_trajectory_write')},\n")
            f.write(f"  ntwr={md_parameters_dict.get('production_params').get('production_restart_write')},\n")
            f.write(f"  iwrap=1,\n")
            f.write(f"/\n")
            f.write(f"END\n")

    def _prepare_md_execution_script(self, md_assay_folder):
        
        execution_script_path = os.path.join(md_assay_folder, "run_md.sh")
        with open(execution_script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("echo 'Running min1'\n")
            f.write("pmemd.cuda -O -i min1.in -o min1.out -p complex.prmtop -c complex.inpcrd -r min1.crd -ref complex.inpcrd\n")
            f.write("echo 'Running min2'\n")
            f.write("pmemd.cuda -O -i min2.in -o min2.out -p complex.prmtop -c complex.inpcrd -r min2.crd -ref min1.crd\n")
            f.write("echo 'Running heating'\n")
            f.write("pmemd.cuda -O -i heating.in -o heating.out -p complex.prmtop -c min2.crd -r heating.crd -ref min2.crd\n")
            f.write("echo 'Running equilibration'\n")
            f.write("pmemd.cuda -O -i equilibration.in -o equilibration.out -p complex.prmtop -c heating.crd -r equilibration.crd -ref heating.crd\n")
            f.write("echo 'Running production'\n")
            f.write("pmemd.cuda -O -i production.in -o production.out -p complex.prmtop -c equilibration.crd -r production.crd -ref equilibration.crd\n")

        os.chmod(execution_script_path, 0o755)
    
    def _start_md_simulation(self, md_assay_folder):
        
        # Query if the user wants to run in the foreground or background
        run_mode = input("\n⚙️  Do you want to run the MD simulation in the foreground or background? (fg/bg) [default: fg]: ").strip().lower() or 'fg'
        execution_script_path = os.path.join(md_assay_folder, "run_md.sh")
        try:
            import subprocess

            if run_mode in ['fg', 'foreground']:
                print("\n▶️  Starting MD simulation in the foreground...")
                subprocess.run([execution_script_path], check=True, cwd=md_assay_folder)

            elif run_mode in ['bg', 'background']:
                print("\n▶️  Starting MD simulation in the background...")
                subprocess.Popen([execution_script_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=md_assay_folder)

            else:
                print("❌ Invalid option. Please choose 'fg' or 'bg'.")
                return
        except Exception as e:
            print(f"❌ Error starting MD simulation: {e}")