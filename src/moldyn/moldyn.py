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
        tleap_params = self._set_tleap_params(params)
        params['tleap_params'] = tleap_params

        ## Set minimization parameters

        ## Set heating parameters

        ## Set equilibration parameters

        ## Set production parameters

            
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

            conn = sqlite3.connect(self.__md_registers_db)
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
            docking_assay_id, assay_name, assay_folder_path = self._select_docking_assay(cursor)
            conn.close()

            print(f"Docking Assay ID: {docking_assay_id}")
            
            # Select unique ligand molecules in the selected docking assay
            ligname, docking_results_db = self._select_unique_ligands_in_docking_assay(docking_assay_id, assay_name, assay_folder_path)

            # Select pose id for the selected ligand to perform MD assay
            pose_id = self._select_ligand_pose_for_md_assay(docking_results_db, assay_name, ligname)

            # Restore the selected docked pose using ligand_management module
            pdb_dict = lm.restore_single_docked_pose(docking_results_db, ligname, pose_id)

            # Query user for MD assay description
            assay_description = input("\n📝 Enter MD assay description (optional): ").strip()
            if not assay_description:
                assay_description = f"MD assay for ligand {ligname} (pose {pose_id})"

            # Create MD assay folder structure and prepare input files for the selected MD engine
            new_assay_id, new_assay_folder = self._create_md_assay_folder()

            # Create MD register entry in the md_assays table
            self._create_md_assay_register_entry(new_assay_id, new_assay_folder, assay_description, docking_assay_id, ligname, pose_id)

        except Exception as e:
            print(f"❌ Error performing MD assay: {e}")
        
    def _select_docking_assay(self, cursor):
        
        # Check if docking_assays table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='docking_assays';")
            if not cursor.fetchone():
                print("❌ No docking assays found. The 'docking_assays' table does not exist.")
                return

            # Fetch all docking assays
            cursor.execute("SELECT assay_id, assay_name, notes, assay_folder_path FROM docking_assays;")
            assays = cursor.fetchall()

            if not assays:
                print("❌ No docking assays found in the database.")
                return

            print(f"\n🧬 DOCKING ASSAYS IN PROJECT '{self.name}':")
            print("=" * 70)
            for assay in assays:
                assay_id, assay_name, notes, assay_folder_path = assay
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
                    # Further actions can be implemented here
                    return assay_id, assay_name, assay_folder_path
                else:
                    print("❌ Invalid Assay ID. Please try again.")
    
    def _select_unique_ligands_in_docking_assay(self, assay_id, assay_name, assay_folder_path):
        try:
            import sqlite3

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
                    return selected_ligand_name, docking_results_db
                else:
                    print("❌ Invalid Ligand ID. Please try again.")
            
        except Exception as e:
            print(f"❌ Error fetching unique ligands: {e}")
    
    def _select_ligand_pose_for_md_assay(self, docking_results_db, assay_name, ligname):
        try:
            import sqlite3

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
                    # Further actions to set up MD assay can be implemented here
                    return selection
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

    def _create_md_assay_register_entry(self, new_assay_id, new_assay_folder, assay_description, docking_assay_id=None, ligname=None, pose_id=None):
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
            }
            # Create md_assays table if it doesn't exist
            dbm.create_table_from_columns_dict(cursor, 'md_assays', columns_dict)

            # Update legacy table columns 
            dbm.update_legacy_table_columns(cursor, 'md_assays', columns_dict)

            # Remove legacy table columns
            dbm.remove_legacy_table_columns(cursor, 'md_assays', columns_dict)

            print(f"Docking Assay ID: {docking_assay_id}")

            # Prepare data values dictionary (excluding auto-generated columns)
            data_dict = {
                'md_assay': f'assay_{new_assay_id}',
                'description': assay_description,
                'assay_folder_path': new_assay_folder,
                'docking_assay_id': f'assay_{docking_assay_id}',
                'ligand_name': ligname,
                'pose_id': pose_id,
            }
            # Insert dinamically data into md_assays table
            dbm.insert_data_dinamically_into_table(cursor, 'md_assays', data_dict)
            conn.commit()
            conn.close()

        except Exception as e:
            print(f"❌ Error creating MD assay register entry: {e}")