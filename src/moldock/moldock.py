import ast
import os
import sqlite3
from turtle import pd
from tidyscreen import tidyscreen
import sys
from typing import Dict, List, Tuple, Optional, Any
from rdkit.Chem import AllChem
#from rdkit import Chem
import shutil

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
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject

class MolDock:
    """
    MolDock class for managing molecular docking assays within a project.
    Uses an ActivateProject object to access project information and database functionality.
    
    """
    
    def __init__(self, project_obj: ActivateProject):
        """
        Initialize MolDock with an ActivateProject object.
        
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
        
        # Set up receptors path within the project directory
        self.__receptor_path = os.path.join(self.path, 'docking/receptors')
        
        # Set up docking assays registers database path within the project directory
        self.__docking_registers_db = os.path.join(self.path, 'docking/docking_registers', 'registers.db')

        # Set up docking params registers database path within the project directory
        self.__docking_params_db = os.path.join(self.path, 'docking/params', 'params.db')

        # Ensure data directory exists
        data_dir = os.path.dirname(self.__chemspace_db)
        os.makedirs(data_dir, exist_ok=True)

    def create_docking_method(self):
        """
        Will prompt the user for the following method parameters:
        - Docking engine: options: AutoDockGPU, Autodock, Vina

        The method will create a docking_methods.db under project path in the docking/docking_registers folder
        """
        try:
            print(f"🧬 CREATE DOCKING METHOD")
            print("=" * 50)
            
            # Ensure docking directories exist
            docking_registers_dir = os.path.dirname(self.__docking_registers_db)
            os.makedirs(docking_registers_dir, exist_ok=True)
            
            # Define available docking engines with details
            docking_engines = {
                '1': {
                    'name': 'AutoDockGPU',
                    'description': 'GPU-accelerated AutoDock for high-throughput docking',
                    'requirements': ['GPU with CUDA support', 'AutoDockGPU installation'],
                },
                '2': {
                    'name': 'AutoDock',
                    'description': 'Classical AutoDock CPU-based docking',
                    'requirements': ['AutoDock installation', 'AutoDockTools'],
                },
                '3': {
                    'name': 'Vina',
                    'description': 'AutoDock Vina - fast and accurate docking',
                    'requirements': ['AutoDock Vina installation'],
                }
            }
            
            # Display available docking engines
            print(f"📋 Available Docking Engines:")
            print("-" * 70)
            
            for key, engine in docking_engines.items():
                print(f"{key}. {engine['name']}")
                print(f"   📝 {engine['description']}")
                print(f"   📋 Requirements: {', '.join(engine['requirements'])}")
                print("-" * 70)
            
            # Get user selection
            while True:
                try:
                    selection = input(f"\n🧬 Select docking engine (1-3) or 'cancel': ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Docking method creation cancelled")
                        return None
                    
                    if selection in docking_engines:
                        selected_engine = docking_engines[selection]
                        print(f"\n✅ Selected: {selected_engine['name']}")
                        break
                    else:
                        print("❌ Invalid selection. Please enter 1, 2, or 3")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ Docking method creation cancelled")
                    return None
            
            # Get method name
            while True:
                try:
                    method_name = input(f"\n📝 Enter method name (or 'cancel'): ").strip()
                    
                    if method_name.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Docking method creation cancelled")
                        return None
                    
                    if not method_name:
                        print("❌ Method name cannot be empty")
                        continue
                    
                    # Validate method name (no special characters)
                    if not method_name.replace('_', '').replace('-', '').replace(' ', '').isalnum():
                        print("❌ Method name can only contain letters, numbers, spaces, hyphens, and underscores")
                        continue
                    
                    break
                    
                except KeyboardInterrupt:
                    print("\n❌ Docking method creation cancelled")
                    return None
            
            # Get optional description
            description = input(f"\n📄 Enter method description (optional): ").strip()
            if not description:
                description = f"Docking method using {selected_engine['name']}"
            
            # Get additional parameters based on engine
            additional_params = self._get_engine_specific_parameters(selected_engine['name'])
                        
            # Get ligand preparation options
            ligand_prep_params = self._get_ligand_preparation_parameters()
            
            if additional_params is None:
                print("❌ Parameter configuration cancelled")
                return None
            
            # Create or connect to docking methods database
            methods_db_path = os.path.join(docking_registers_dir, 'docking_methods.db')
            
            try:
                import sqlite3
                conn = sqlite3.connect(methods_db_path)
                cursor = conn.cursor()
                
                # Create docking_methods table with ligand_prep_params column
                cursor.execute('''
                    CREATE TABLE IF NOT EXISTS docking_methods (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        method_name TEXT UNIQUE NOT NULL,
                        docking_engine TEXT NOT NULL,
                        description TEXT,
                        parameters TEXT,
                        ligand_prep_params TEXT,
                        created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                    )
                ''')
                
                # Check if method name already exists
                cursor.execute("SELECT COUNT(*) FROM docking_methods WHERE method_name = ?", (method_name,))
                if cursor.fetchone()[0] > 0:
                    print(f"⚠️  Method name '{method_name}' already exists")
                    while True:
                        overwrite = input("Overwrite existing method? (y/n): ").strip().lower()
                        if overwrite in ['y', 'yes']:
                            # Update existing method
                            cursor.execute('''
                                UPDATE docking_methods 
                                SET docking_engine = ?, description = ?, parameters = ?, ligand_prep_params = ?, 
                                    created_date = CURRENT_TIMESTAMP
                                WHERE method_name = ?
                            ''', (
                                selected_engine['name'],
                                description,
                                self._serialize_parameters(additional_params),
                                self._serialize_parameters(ligand_prep_params),
                                method_name
                            ))
                            print(f"✅ Updated existing docking method '{method_name}'")
                            break
                        elif overwrite in ['n', 'no']:
                            print("❌ Method creation cancelled - name already exists")
                            conn.close()
                            return None
                        else:
                            print("❌ Please answer 'y' or 'n'")
                            continue
                else:
                    # Insert new method
                    cursor.execute('''
                        INSERT INTO docking_methods 
                        (method_name, docking_engine, description, parameters, ligand_prep_params)
                        VALUES (?, ?, ?, ?, ?)
                    ''', (
                        method_name,
                        selected_engine['name'],
                        description,
                        self._serialize_parameters(additional_params),
                        self._serialize_parameters(ligand_prep_params)
                    ))
                    print(f"✅ Created new docking method '{method_name}'")
                
                method_id = cursor.lastrowid or cursor.execute("SELECT id FROM docking_methods WHERE method_name = ?", (method_name,)).fetchone()[0]
                
                conn.commit()
                conn.close()
                
                # Display method summary
                self._display_method_summary(method_id, method_name, selected_engine, description, additional_params, ligand_prep_params, methods_db_path)
                
                return {
                    'method_id': method_id,
                    'method_name': method_name,
                    'docking_engine': selected_engine['name'],
                    'description': description,
                    'parameters': additional_params,
                    'ligand_prep_params': ligand_prep_params,
                    'database_path': methods_db_path
                }
                
            except Exception as e:
                print(f"❌ Error creating docking methods database: {e}")
                return None
                
        except Exception as e:
            print(f"❌ Error in create_docking_method: {e}")
            return None

    def list_docking_methods(self):
        """
        Will read the docking_methods.db under project path in the docking/docking_registers folder and inform the user about available methods
        """
        
        import sqlite3
        import json
        
        try:
            # Path to docking methods database
            docking_registers_dir = os.path.dirname(self.__docking_registers_db)
            methods_db_path = os.path.join(docking_registers_dir, 'docking_methods.db')
            
            # Check if database exists
            if not os.path.exists(methods_db_path):
                print(f"❌ No docking methods database found")
                print(f"   Expected location: {methods_db_path}")
                print(f"\n💡 Create a docking method first using create_docking_method()")
                return None
            
            # Connect to database
            conn = sqlite3.connect(methods_db_path)
            cursor = conn.cursor()
            
            # Check if table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='docking_methods'")
            if not cursor.fetchone():
                print(f"❌ Docking methods table not found in database")
                conn.close()
                return None
            
            # Retrieve all docking methods
            cursor.execute('''
                SELECT id, method_name, docking_engine, description, parameters, 
                    ligand_prep_params, created_date
                FROM docking_methods
                ORDER BY created_date ASC
            ''')
            
            methods = cursor.fetchall()
            conn.close()
            
            # Check if any methods exist
            if not methods:
                print(f"📋 NO DOCKING METHODS FOUND")
                print(f"   Database exists but contains no methods")
                print(f"\n💡 Create a docking method using create_docking_method()")
                return None
            
            # Display methods
            print(f"\n🧬 AVAILABLE DOCKING METHODS")
            print("=" * 80)
            print(f"   Found {len(methods)} method(s) in database")
            print(f"   Database: {methods_db_path}")
            print("=" * 80)
            
            methods_list = []
            
            for idx, (method_id, method_name, engine, description, params_json, 
                    ligand_prep_json, created_date) in enumerate(methods, 1):
                
                # Parse JSON parameters
                try:
                    parameters = json.loads(params_json) if params_json else {}
                    ligand_prep_params = json.loads(ligand_prep_json) if ligand_prep_json else {}
                except json.JSONDecodeError:
                    parameters = {}
                    ligand_prep_params = {}
                
                # Store method info
                method_info = {
                    'id': method_id,
                    'method_name': method_name,
                    'docking_engine': engine,
                    'description': description,
                    'parameters': parameters,
                    'ligand_prep_params': ligand_prep_params,
                    'created_date': created_date
                }
                methods_list.append(method_info)
                
                # Display basic method details only
                print(f"\n{idx}. 🏷️  {method_name}")
                print(f"   {'─' * 76}")
                print(f"   📋 ID: {method_id}")
                print(f"   🧬 Engine: {engine}")
                print(f"   📝 Description: {description or 'No description'}")
                print(f"   📅 Created: {created_date}")
                print(f"   {'─' * 76}")
            
            print("=" * 80)
            print(f"\n💡 Type 'details <number>' to view engine and ligand preparation parameters")
            print(f"💡 Example: details 1")
            
            # Interactive details viewing
            while True:
                try:
                    user_input = input(f"\nEnter command (or press Enter to exit): ").strip()
                    
                    if not user_input:
                        break
                    
                    # Check if user wants details
                    if user_input.lower().startswith('details'):
                        parts = user_input.split()
                        if len(parts) < 2:
                            print(f"❌ Please specify a method number. Example: details 1")
                            continue
                        
                        try:
                            method_number = int(parts[1])
                            if method_number < 1 or method_number > len(methods_list):
                                print(f"❌ Invalid method number. Choose between 1 and {len(methods_list)}")
                                continue
                            
                            # Display detailed parameters for selected method
                            selected_method = methods_list[method_number - 1]
                            print(f"\n{'═' * 80}")
                            print(f"🔍 DETAILED PARAMETERS FOR: {selected_method['method_name']}")
                            print(f"{'═' * 80}")
                            
                            # Display engine parameters
                            if selected_method['parameters']:
                                print(f"\n⚙️  ENGINE PARAMETERS ({selected_method['docking_engine']}):")
                                print(f"   {'─' * 76}")
                                for key, value in selected_method['parameters'].items():
                                    param_display = key.replace('_', ' ').title()
                                    print(f"   • {param_display}: {value}")
                            else:
                                print(f"\n⚙️  Engine Parameters: None configured")
                            
                            # Display ligand preparation parameters
                            if selected_method['ligand_prep_params']:
                                print(f"\n🧪 LIGAND PREPARATION PARAMETERS:")
                                print(f"   {'─' * 76}")
                                for key, value in selected_method['ligand_prep_params'].items():
                                    param_display = key.replace('_', ' ').title()
                                    print(f"   • {param_display}: {value}")
                            else:
                                print(f"\n🧪 Ligand Preparation: Default settings")
                            
                            print(f"{'═' * 80}")
                            
                        except ValueError:
                            print(f"❌ Invalid number format. Please use: details <number>")
                            continue
                    else:
                        print(f"❌ Unknown command. Use 'details <number>' to view parameters")
                
                except KeyboardInterrupt:
                    print(f"\n\n👋 Exiting method list")
                    break
            
            return methods_list
            
        except sqlite3.Error as e:
            print(f"❌ Database error: {e}")
            return None
        except Exception as e:
            print(f"❌ Error listing docking methods: {e}")
            return None

    def _get_engine_specific_parameters(self, engine_name: str) -> Optional[Dict[str, Any]]:
        """
        Get engine-specific parameters from user input.
        
        Args:
            engine_name (str): Name of the selected docking engine
            
        Returns:
            Optional[Dict[str, Any]]: Engine parameters or None if cancelled
        """
        try:
            print(f"\n⚙️  CONFIGURE {engine_name.upper()} PARAMETERS")
            print("-" * 50)
            
            parameters = {}
            
            if engine_name == 'AutoDockGPU':
                # AutoDockGPU specific parameters
                
                parameters['heuristics'] = self._get_parameter_choice(
                    "Apply heuristics algorithm?",
                    [0, 1],
                    default=1
                )
                
                parameters['heurmax'] = self._input_parameter_choice(
                    "Asymptotic heuristics # evals limit (smooth limit: 10000000-20000000. Default: 12000000)",
                    default=12000000
                )
                
                
                parameters['autostop'] = self._get_parameter_choice(
                    "Automatic stopping criterion based on convergence",
                    [0, 1],
                    default=1
                )
                
                parameters['asfreq'] = self._input_parameter_choice(
                    "AutoStop testing frequency (in # of generations: 1-20. Default: 5)",
                    default=5
                )
                
                parameters['nrun'] = self._input_parameter_choice(
                    "LGA runs (5-500). Default: 20",
                    default=20
                )
                
                parameters['nev'] = self._input_parameter_choice(
                    "Score evaluations (max.) per LGA run",
                    default=2500000
                )
                
                parameters['ngen'] = self._input_parameter_choice(
                    "Generations (max.) per LGA run. Default: 42000",
                    default=42000
                )
                
                parameters['lsmet'] = self._get_parameter_choice(
                    "Local-search method",
                    ["ad"],
                    default="ad"
                )
                
                parameters['lsit'] = self._input_parameter_choice(
                    "Local-search iterations (max.)",
                    default=300
                )
                
                parameters['psize'] = self._input_parameter_choice(
                    "Population size. Default: 150",
                    default=150
                )
                
                parameters['mrat'] = self._input_parameter_choice(
                    "Mutation rate (%). Default: 2",
                    default=2
                )
                
                parameters['crat'] = self._input_parameter_choice(
                    "Crossover rate (%). Default: 80",
                    default=80
                )
                
                parameters['lsrat'] = self._input_parameter_choice(
                    "Local-search rate (%). Default: 100",
                    default=100
                )
                
                parameters['trat'] = self._input_parameter_choice(
                    "Tournament (selection) rate (%). Default: 60",
                    default=60
                )
                
                parameters['dmov'] = self._input_parameter_choice(
                    "Maximum LGA movement delta (Angstroms). Default: 6",
                    default=6
                )
                
                parameters['dang'] = self._input_parameter_choice(
                    "Maximum LGA angle delta (Degrees). Default: 90",
                    default=90
                )
                
                parameters['rholb'] = self._input_parameter_choice(
                    "Solis-Wets lower bound of rho parameter. Default=0.01",
                    default=0.01
                )
                
                parameters['lsmov'] = self._input_parameter_choice(
                    "Solis-Wets movement delta (A). Default=2",
                    default=2
                )
                
                parameters['lsang'] = self._input_parameter_choice(
                    "Solis-Wets angle delta (Degrees). Default=75",
                    default=75
                )
                
                parameters['cslim'] = self._input_parameter_choice(
                    "Solis-Wets cons. success/failure limit to adjust rho. Default=4",
                    default=4
                )
                
                parameters['stopstd'] = self._input_parameter_choice(
                    "AutoStop energy standard deviation tolerance (kcal/mol). Default=0.15",
                    default=0.15
                )
                
                parameters['contact_analysis'] = self._get_parameter_choice(
                    "Perform distance-based analysis",
                    [0, 1],
                    default=1
                )
                
                parameters['xmloutput'] = self._get_parameter_choice(
                    "Specify if xml output format is wanted",
                    [0, 1],
                    default=0
                )
                
                parameters['clustering'] = self._get_parameter_choice(
                    "Output clustering analysis in dlg and/or xml file ",
                    [0, 1],
                    default=1
                )
                
                parameters['output-cluster-poses'] = self._get_parameter_choice(
                    "Output up to a certain number of poses per cluster (0: all)",
                    [0, 1],
                    default=0
                )
                
                parameters['hsym'] = self._get_parameter_choice(
                    "Handle symmetry in RMSD calc.",
                    [0, 1],
                    default=1
                )
                
                parameters['rmstol'] = self._input_parameter_choice(
                    "RMSD clustering tolerance (Angstroms). Default=2",
                    default=2
                )
                
            elif engine_name == 'AutoDock':
                # AutoDock specific parameters
                parameters['search_algorithm'] = self._get_parameter_choice(
                    "Search algorithm",
                    ['Lamarckian GA', 'Simulated Annealing', 'Genetic Algorithm'],
                    default='Lamarckian GA'
                )
                
                parameters['num_runs'] = self._get_parameter_integer(
                    "Number of docking runs per ligand",
                    default=10, min_val=1, max_val=100
                )
                
                parameters['population_size'] = self._get_parameter_integer(
                    "GA population size",
                    default=150, min_val=50, max_val=500
                )
                
                parameters['max_generations'] = self._get_parameter_integer(
                    "Maximum GA generations",
                    default=27000, min_val=1000, max_val=100000
                )
                
                parameters['energy_evaluations'] = self._get_parameter_integer(
                    "Maximum energy evaluations",
                    default=2500000, min_val=100000, max_val=10000000
                )
                
            elif engine_name == 'Vina':
                # AutoDock Vina specific parameters
                
                parameters['sf_name'] = self._get_parameter_choice(
                    "Scoring function name to use",
                    ['vina', 'ad4'],
                    default='vina'
                )

                parameters['cpu'] = self._get_parameter_integer(
                    "Number of CPU to use (default: 0; use all of them)",
                    default=0, min_val=0, max_val=64
                )
                
                parameters['exhaustiveness'] = self._get_parameter_integer(
                    "Search exhaustiveness (1-32)",
                    default=8, min_val=1, max_val=32
                )
                
                parameters['n_poses'] = self._get_parameter_integer(
                    "Number of poses to generate in MC search (1-100)",
                    default=20, min_val=1, max_val=100
                )

                parameters['min_rmsd'] = self._get_parameter_float(
                    "minimum RMSD difference between poses in MC search (0.1-3)",
                    default=1, min_val=0.1, max_val=3
                )

                parameters['num_modes'] = self._get_parameter_integer(
                    "Number of binding modes to return",
                    default=9, min_val=1, max_val=20
                )

                parameters['max_evals'] = self._get_parameter_integer(
                    "Maximum number of evaluations in MC search (default: 0; use heuristics rules)",
                    default=0, min_val=0, max_val=10000000
                )
                
                parameters['energy_range'] = self._get_parameter_float(
                    "Energy range (kcal/mol)",
                    default=3.0, min_val=1.0, max_val=10.0
                )
                
                parameters['cpu_cores'] = self._get_parameter_integer(
                    "Number of CPU cores to use (0 = auto)",
                    default=0, min_val=0, max_val=64
                )
            
                parameters['n_poses'] = self._get_parameter_integer(
                    "Number of poses to write in output",
                    default=9, min_val=1, max_val=100
                )

                parameters['energy_range'] = self._get_parameter_float(
                    "Maximum energy difference from best pose (kcal/mol)",
                    default=3, min_val=0.1, max_val=10.0
                )
            
            return parameters if parameters else None
            
        except KeyboardInterrupt:
            print("\n❌ Parameter configuration cancelled")
            return None
        except Exception as e:
            print(f"❌ Error configuring parameters: {e}")
            return None

    def _get_parameter_choice(self, param_name: str, choices: List[str], default: str) -> str:
        """Get parameter choice from user with validation."""
        print(f"\n📋 {param_name}:")
        for i, choice in enumerate(choices, 1):
            marker = " (default)" if choice == default else ""
            print(f"  {i}. {choice}{marker}")
        
        while True:
            try:
                selection = input(f"Select {param_name.lower()} (1-{len(choices)} or press Enter for default): ").strip()
                
                if not selection:
                    return default
                
                try:
                    idx = int(selection) - 1
                    if 0 <= idx < len(choices):
                        return choices[idx]
                    else:
                        print(f"❌ Please enter 1-{len(choices)}")
                        continue
                except ValueError:
                    print("❌ Please enter a number")
                    continue
                    
            except KeyboardInterrupt:
                raise
    
    def _input_parameter_choice(self, param_name: str, default: str = None) -> str:
        """
        Prompt the user to input a parameter value directly (no predefined options).

        Args:
            param_name (str): Name or description of the parameter to prompt for.
            default (str, optional): Default value to use if the user presses Enter.

        Returns:
            str: The value entered by the user, or the default if provided and no input is given.
        """
        while True:
            try:
                print(f"\n📋 {param_name}:")
                prompt = f"Enter value for {param_name}"
                if default is not None:
                    prompt += f" (default: {default})"
                prompt += ": "
                value = input(prompt).strip()
                if not value:
                    if default is not None:
                        return default
                    else:
                        return None
                return value
            except KeyboardInterrupt:
                print("\n❌ Parameter input cancelled.")
                return None
    
    def _get_parameter_integer(self, param_name: str, default: int, min_val: int, max_val: int) -> int:
        """Get integer parameter from user with validation."""
        while True:
            try:
                value_str = input(f"📊 {param_name} ({min_val}-{max_val}, default: {default}): ").strip()
                
                if not value_str:
                    return default
                
                try:
                    value = int(value_str)
                    if min_val <= value <= max_val:
                        return value
                    else:
                        print(f"❌ Value must be between {min_val} and {max_val}")
                        continue
                except ValueError:
                    print("❌ Please enter a valid integer")
                    continue
                    
            except KeyboardInterrupt:
                raise

    def _get_parameter_float(self, param_name: str, default: float, min_val: float, max_val: float) -> float:
        """Get float parameter from user with validation."""
        while True:
            try:
                value_str = input(f"📊 {param_name} ({min_val}-{max_val}, default: {default}): ").strip()
                
                if not value_str:
                    return default
                
                try:
                    value = float(value_str)
                    if min_val <= value <= max_val:
                        return value
                    else:
                        print(f"❌ Value must be between {min_val} and {max_val}")
                        continue
                except ValueError:
                    print("❌ Please enter a valid number")
                    continue
                    
            except KeyboardInterrupt:
                raise

    def _serialize_parameters(self, parameters: Dict[str, Any]) -> str:
        """Serialize parameters dictionary to JSON string for database storage."""
        try:
            import json
            return json.dumps(parameters, indent=2)
        except Exception:
            return str(parameters)

    def _display_method_summary(self, method_id: int, method_name: str, engine_info: Dict,
                            description: str, parameters: Dict[str, Any], ligand_prep_params: Dict[str, Any], db_path: str) -> None:
        """Display summary of created docking method."""
        try:
            print(f"\n🎯 DOCKING METHOD CREATED SUCCESSFULLY")
            print("=" * 60)
            print(f"📋 Method ID: {method_id}")
            print(f"🏷️  Method Name: {method_name}")
            print(f"🧬 Docking Engine: {engine_info['name']}")
            print(f"📝 Description: {description}")
            print(f"🗂️  Database: {db_path}")
            
            print(f"\n⚙️  CONFIGURED PARAMETERS:")
            for key, value in parameters.items():
                print(f"   • {key.replace('_', ' ').title()}: {value}")

            print(f"\n🧪 LIGAND PREPARATION PARAMETERS:")
            if ligand_prep_params:
                for key, value in ligand_prep_params.items():
                    print(f"   • {key.replace('_', ' ').title()}: {value}")
            else:
                print("   • None")

            print(f"\n📋 REQUIREMENTS:")
            for req in engine_info['requirements']:
                print(f"   • {req}")

            print(f"\n💡 NEXT STEPS:")
            print("   1. Ensure all requirements are installed")
            print("   2. Use this method for docking experiments")
            print("   3. Configure receptor files for docking")

            print("=" * 60)

        except Exception as e:
            print(f"⚠️  Error displaying method summary: {e}")

    def dock_table(self, show_details: bool = True, 
                filter_by_type: Optional[str] = None,
                sort_by: str = 'name',
                show_all_tables: bool = False,
                clean_ligand_files: bool = True,
                notes: Optional[str] = None) -> Optional[Any]:
        """
        List tables available in the chemspace database, select a docking method, and allow selection for docking preparation.
        By default, only shows tables ready for docking (with SDF blobs).
        Creates a docking assay registry upon successful completion.
        
        Allows user to choose between:
        - Running in current console (foreground)
        - Running as background process (detached, survives SSH disconnect)
        
        Args:
            show_details (bool): Whether to show detailed information including compound counts
            filter_by_type (Optional[str]): Filter tables by type ('original', 'reaction_products', 'filtered_compounds', etc.)
            sort_by (str): Sort criteria ('name', 'compounds', 'type', 'date')
            show_all_tables (bool): If True, shows all tables regardless of docking readiness
            clean_ligand_files (bool): Whether to clean up ligand files after processing
            notes (Optional[str]): Optional notes to store in the assay registry; prompts if None
            
        Returns:
            Optional[Any]: Dictionary with selected table name, docking method, assay registry info and SMILES data if requested, or None if cancelled
        """
        try:
            print(f"🧬 MOLDOCK - Table and Method Selection for Docking")
            print("=" * 70)
            
            # Step 1: Select docking method first and return a dictionary of method parameters (selected_method)
            selected_method = self._select_docking_method()
            
            if not selected_method:
                print("❌ No docking method selected")
                return None
            
            print(f"\n✅ Selected docking method: '{selected_method['method_name']}' ({selected_method['docking_engine']})")
            print("-" * 70)
            
            # Step 1b: Prompt user for execution mode
            execution_mode = self._prompt_execution_mode()
            if execution_mode is None:
                print("❌ Docking preparation cancelled")
                return None
            
            # Step 2: Check if chemspace database exists
            if not os.path.exists(self.__chemspace_db):
                print(f"❌ ChemSpace database not found: {self.__chemspace_db}")
                print("   Please ensure ChemSpace has been initialized for this project.")
                return None
            
            # Step 3: Connect to chemspace database and get available tables
            try:
                import sqlite3
                conn = sqlite3.connect(self.__chemspace_db)
                cursor = conn.cursor()
                
                # Get all tables (reusing ChemSpace pattern)
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")
                available_tables = [row[0] for row in cursor.fetchall()]
                conn.close()
                
            except Exception as e:
                print(f"❌ Error accessing ChemSpace database: {e}")
                return None
            
            if not available_tables:
                print("📋 No tables found in ChemSpace database")
                print("   Load some compound data into ChemSpace first.")
                return None
            
            print(f"📊 Found {len(available_tables)} table(s) in ChemSpace database")
            
            # Step 4: Collect table information using existing patterns
            table_data = []
            ready_table_data = []
            not_ready_count = 0
            
            if show_details:
                print("📊 Analyzing tables for docking preparation...")
            
            for table_name in available_tables:
                try:
                    # Get compound count (reusing ChemSpace pattern)
                    table_info = self._get_table_info_for_docking(table_name)
                    
                    # Apply type filter if specified
                    if filter_by_type and table_info['type'] != filter_by_type:
                        continue
                    
                    # Separate ready vs not ready tables
                    if table_info['docking_status'] == 'ready for docking':
                        ready_table_data.append(table_info)
                    else:
                        not_ready_count += 1
                    
                    # Add to all tables list (for show_all_tables option)
                    table_data.append(table_info)
                    
                except Exception as e:
                    print(f"   ⚠️  Error analyzing table '{table_name}': {e}")
                    continue
            
            # Step 5: Show summary of table readiness
            total_analyzed = len(table_data)
            ready_count = len(ready_table_data)
            
            print(f"\n📋 TABLE READINESS SUMMARY:")
            print(f"   📊 Total tables analyzed: {total_analyzed}")
            print(f"   ✅ Ready for docking: {ready_count}")
            print(f"   ❌ Not ready for docking: {not_ready_count}")
            
            # Step 6: Determine which tables to display
            if not show_all_tables:
                # Default behavior: show only ready tables
                tables_to_display = ready_table_data
                
                if not ready_table_data:
                    print(f"\n⚠️  NO TABLES READY FOR DOCKING")
                    print(f"   All {total_analyzed} tables need 3D molecular preparation")
                    print(f"   💡 Suggestion: Use ChemSpace.generate_mols_in_table() to prepare molecules")
                    
                    # Offer to show all tables
                    if total_analyzed > 0:
                        print(f"\n🔄 OPTIONS:")
                        print(f"   1. Show all tables (including not ready)")
                        print(f"   2. Cancel and prepare molecules first")
                        
                        while True:
                            try:
                                choice = input(f"\n🧬 Choose option (1/2): ").strip()
                                
                                if choice == '1':
                                    print(f"\n📋 Showing all {total_analyzed} tables (including not ready)...")
                                    tables_to_display = table_data
                                    show_all_tables = True
                                    break
                                elif choice == '2' or choice.lower() in ['cancel', 'quit', 'exit']:
                                    print("❌ Docking preparation cancelled")
                                    return None
                                else:
                                    print("❌ Invalid choice. Please enter 1 or 2")
                                    continue
                                    
                            except KeyboardInterrupt:
                                print("\n❌ Docking preparation cancelled")
                                return None
                    else:
                        return None
                else:
                    print(f"\n✅ Showing {ready_count} table(s) ready for docking")
                    
            else:
                # Show all tables when explicitly requested
                tables_to_display = table_data
                print(f"\n📋 Showing all {total_analyzed} tables (ready + not ready)")
            
            if not tables_to_display:
                if filter_by_type:
                    print(f"❌ No tables found with type '{filter_by_type}'")
                else:
                    print("❌ No valid tables found for docking preparation")
                return None
            
            # Step 7: Sort tables (reusing existing logic)
            tables_to_display = self._sort_table_data_for_docking(tables_to_display, sort_by)
            
            # Step 8: Display table list for selection with readiness context
            selected_table = self._display_and_select_docking_table(
                tables_to_display, filter_by_type, show_all_tables, ready_count, not_ready_count
            )
            
            if selected_table:
                print(f"\n✅ Selected table for docking: '{selected_table}'")
                
                # Step 9: Show docking setup summary and get user choice
                user_choice = self._show_docking_setup_summary(selected_table, selected_method, tables_to_display)
                
                if user_choice == 1:
                    # Step 10: Create docking assay registry before returning results
                    print(f"\n📝 Creating docking assay registry...")
                    note_text = notes
                    if note_text is None:
                        try:
                            note_text = input(f"\n📝 Enter notes for this docking assay (press Enter to leave blank): ")
                        except KeyboardInterrupt:
                            print("\n❌ Docking preparation cancelled")
                            return None
                    if note_text is None:
                        note_text = ""
                    
                    # Get selected table info for registry
                    selected_table_info = next((t for t in tables_to_display if t['name'] == selected_table), None)
                    
                    # Create assay registry
                    assay_registry = self._create_docking_assay_registry(
                        selected_table, 
                        selected_method, 
                        selected_table_info,
                        user_choice,
                        note_text
                    )
                    
                    if not assay_registry:
                        print("⚠️  Warning: Could not create docking assay registry, but continuing...")
                    else:
                        print(f"✅ Docking assay registry created with ID: {assay_registry['assay_id']}")
                    
                    if user_choice == 1:
                        # Load SMILES and return complete setup
                        print(f"\n🧪 Loading SMILES from table '{selected_table}'...")
                        smiles_data = self._get_smiles_from_table(selected_table)
                        
                        if smiles_data:
                            print(f"✅ Loaded {len(smiles_data)} compounds with SMILES")
                            
                            # Update registry with loaded compound count
                            if assay_registry:
                                self._update_assay_registry_compound_count(assay_registry['assay_id'], len(smiles_data))
                            
                        else:
                            print("❌ No valid SMILES found in selected table")
                            
                            # Update registry status to indicate failure
                            if assay_registry:
                                self._update_assay_registry_status(assay_registry['assay_id'], 'failed', 'No valid SMILES found')
                            
                            return None
                        
                        # Select receptor before executing docking
                        print(f"\n🧬 Selecting receptor for docking...")
                        receptor_info = self._select_receptor_from_db()
                        
                        if not receptor_info:
                            print("❌ No receptor selected. Docking cancelled.")
                            return None
                        
                        print(f"✅ Selected receptor: '{receptor_info['pdb_name']}'")
                        
                        # Execute docking based on execution mode
                        if execution_mode == 'foreground':
                            # Run directly in current console
                            self._run_docking_foreground(selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info)
                        else:
                            # Run as background/detached process
                            self._run_docking_background(selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info)
                        
                    else:
                        print(f"❌ Docking engine '{selected_method['docking_engine']}' not yet implemented")
                    
                else:
                    return None
            else:
                print("❌ No table selected for docking preparation")
                return None
                
        except Exception as e:
            print(f"❌ Error in dock_table method: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _prompt_execution_mode(self) -> Optional[str]:
        """
        Prompt user to choose between foreground (current console) or background (detached) execution.
        
        Returns:
            Optional[str]: 'foreground', 'background', or None if cancelled
        """
        print(f"\n🚀 EXECUTION MODE SELECTION")
        print("=" * 70)
        print(f"How would you like to run the docking process?")
        print(f"\nOptions:")
        print(f"   1. FOREGROUND: Run in current console")
        print(f"      - Process will stop if SSH connection is closed")
        print(f"      - Real-time output and progress visible")
        print(f"      - Easier to debug if issues occur")
        print(f"")
        print(f"   2. BACKGROUND: Run as detached process")
        print(f"      - Process continues even if SSH connection is closed")
        print(f"      - Output redirected to log file")
        print(f"      - Better for long-running jobs over SSH")
        print(f"      - Free up console for other work")
        print("=" * 70)
        
        while True:
            try:
                choice = input(f"\n🧬 Choose execution mode (1/2, or 'cancel'): ").strip().lower()
                
                if choice in ['1', 'foreground', 'fg']:
                    print(f"✅ Selected: FOREGROUND execution")
                    return 'foreground'
                elif choice in ['2', 'background', 'bg']:
                    print(f"✅ Selected: BACKGROUND execution")
                    return 'background'
                elif choice in ['cancel', 'quit', 'exit', 'c']:
                    print(f"❌ Execution mode selection cancelled")
                    return None
                else:
                    print(f"❌ Invalid choice. Please enter 1, 2, or 'cancel'")
                    continue
                    
            except KeyboardInterrupt:
                print(f"\n❌ Execution mode selection cancelled")
                return None

    def _run_docking_foreground(self, 
                                selected_table: str, 
                                selected_method: Dict[str, Any], 
                                assay_registry: Dict[str, Any], 
                                clean_ligand_files: bool,
                                receptor_info: Dict[str, Any]
        ) -> None:
        """
        Run docking in foreground (current console).
        
        Args:
            selected_table (str): Name of the selected table
            selected_method (Dict): Docking method information
            assay_registry (Dict): Assay registry information
            clean_ligand_files (bool): Whether to clean ligand files
        """
        try:
            print(f"\n💻 Running docking in FOREGROUND mode...")
            print(f"   Process will run in current console")
            
            ## Dock using AutoDockGPU if selected    
            if selected_method['docking_engine'] == 'AutoDockGPU':
                self._dock_table_with_autodockgpu(selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info)
                docking_mode = 'dlg'

            ## Dock using Vina if selected    
            elif selected_method['docking_engine'] == 'Vina':
                self._dock_table_with_vina(selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info)
                docking_mode = 'vina'
            else:
                print(f"❌ Docking engine '{selected_method['docking_engine']}' not yet implemented")
                return

            ## After the docking run has been completed, execute the corresponding analysis
            results_db_file = self._process_docking_assay_results(assay_registry, docking_mode, clean_ligand_files, max_poses=10)
            print(f"✅ Docking results saved to: {results_db_file}")
            
        except Exception as e:
            print(f"❌ Error during foreground docking execution: {e}")
            import traceback
            traceback.print_exc()

    def _run_docking_background(self, 
                                selected_table: str, 
                                selected_method: Dict[str, Any], 
                                assay_registry: Dict[str, Any], 
                                clean_ligand_files: bool,
                                receptor_info: Dict[str, Any]
            ) -> None:
        """
        Run docking as background/detached process (survives SSH disconnect).
        Uses nohup to detach from parent process and output to log file.
        
        Args:
            selected_table (str): Name of the selected table
            selected_method (Dict): Docking method information
            assay_registry (Dict): Assay registry information
            clean_ligand_files (bool): Whether to clean ligand files
        """
        import subprocess
        import sys
        from datetime import datetime
        
        try:
            print(f"\n🔄 Launching docking as BACKGROUND process...")
            
            # Create log file path
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            log_file = os.path.join(
                self.path, 
                f'docking/logs/docking_assay_{assay_registry["assay_id"]}_{timestamp}.log'
            )
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            
            # Build Python command to run docking
            python_script = f"""
import sys
import os
sys.path.insert(0, '{self.path}')

from tidyscreen import tidyscreen
from tidyscreen.moldock.moldock import MolDock

try:
    # Recreate the project and MolDock objects
    project = tidyscreen.ActivateProject('{self.name}')
    moldock = MolDock(project)
    
    # Retrieve method info from assay registry
    assay_id = {assay_registry['assay_id']}
    
    # Prepare parameters
    selected_table = '{selected_table}'
    selected_method = {selected_method}
    assay_registry = {assay_registry}
    clean_ligand_files = {clean_ligand_files}
    receptor_info = {receptor_info}
    
    # Run docking
    moldock._run_docking_foreground(selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info)
    
    print(f'✅ Background docking process completed successfully')
    sys.exit(0)
    
except Exception as e:
    print(f'❌ Error in background docking process: {{e}}')
    import traceback
    traceback.print_exc()
    sys.exit(1)
    """
            
            # Write script to temporary file
            script_file = os.path.join(
                self.path,
                f'docking/scripts/docking_process_{assay_registry["assay_id"]}_{timestamp}.py'
            )
            os.makedirs(os.path.dirname(script_file), exist_ok=True)
            
            with open(script_file, 'w') as f:
                f.write(python_script)
            
            # Launch process using nohup (survives SSH disconnect)
            with open(log_file, 'w') as logf:
                process = subprocess.Popen(
                    [sys.executable, script_file],
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    #preexec_fn=os.setpgrp if hasattr(os, 'setpgrp') else None,  # Detach from parent process
                    start_new_session=True  # Windows-compatible detachment
                )
            
            process_id = process.pid
            
            print(f"✅ Background docking process started successfully!")
            print(f"   🆔 Process ID: {process_id}")
            print(f"   📝 Log file: {log_file}")
            print(f"   🔄 Assay ID: {assay_registry['assay_id']}")
            print(f"\n💡 You can safely close this SSH connection.")
            print(f"   The docking process will continue running in the background.")
            print(f"\n📊 To monitor progress, use:")
            print(f"   tail -f {log_file}")
            
        except Exception as e:
            print(f"❌ Error launching background docking process: {e}")
            import traceback
            traceback.print_exc()

    def _create_docking_assay_registry(self, table_name: str, docking_method: Dict[str, Any], 
                                    table_info: Optional[Dict[str, Any]], setup_type: int,
                                    notes: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        Create a docking assay registry entry in the docking_assays database and create associated folder.
        
        Args:
            table_name (str): Name of the selected table
            docking_method (Dict): Selected docking method information
            table_info (Optional[Dict]): Table information dictionary
            setup_type (int): Setup type (1=configuration_only, 2=ready_for_docking)
            notes (Optional[str]): Notes to store alongside the registry entry
            
        Returns:
            Optional[Dict[str, Any]]: Assay registry information or None if failed
        """
        try:
            note_text = notes if notes is not None else ""
            # Ensure docking registers directory exists
            docking_registers_dir = os.path.dirname(self.__docking_registers_db)
            os.makedirs(docking_registers_dir, exist_ok=True)
            
            # Create assays directory for storing docking results
            docking_assays_base_dir = os.path.join(self.path, 'docking', 'docking_assays')
            os.makedirs(docking_assays_base_dir, exist_ok=True)
            
            # Create assay registry database path
            assays_db_path = os.path.join(docking_registers_dir, 'docking_assays.db')
            
            import sqlite3
            import json
            from datetime import datetime
            
            conn = sqlite3.connect(assays_db_path)
            cursor = conn.cursor()
            
            # Create docking_assays table with assay_folder_path column
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS docking_assays (
                    assay_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    assay_name TEXT NOT NULL,
                    project_name TEXT,
                    table_name TEXT NOT NULL,
                    docking_method_id INTEGER,
                    docking_method_name TEXT NOT NULL,
                    docking_engine TEXT NOT NULL,
                    compound_count INTEGER DEFAULT 0,
                    prepared_molecules INTEGER DEFAULT 0,
                    setup_type TEXT,
                    status TEXT DEFAULT 'created',
                    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    last_modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    table_readiness_score INTEGER,
                    docking_status TEXT,
                    notes TEXT,
                    configuration TEXT,
                    assay_folder_path TEXT
                )
            ''')
            
            # First, insert a temporary record to get the assay_id
            cursor.execute('''
                INSERT INTO docking_assays (
                    assay_name, project_name, table_name, docking_method_id, docking_method_name,
                    docking_engine, compound_count, prepared_molecules, setup_type, status,
                    table_readiness_score, docking_status, notes, configuration, assay_folder_path
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'temp_assay',  # Temporary name - will be updated
                self.name,  # project name
                table_name,
                docking_method.get('method_id'),
                docking_method['method_name'],
                docking_method['docking_engine'],
                table_info['compound_count'] if table_info else 0,
                table_info['sdf_blob_count'] if table_info else 0,
                'configuration_only' if setup_type == 1 else 'ready_for_docking',
                'created',
                table_info['readiness_score'] if table_info else 0,
                table_info['docking_status'] if table_info else 'unknown',
                note_text,
                '{}',  # Temporary empty configuration - will be updated
                ''     # Temporary empty folder path - will be updated
            ))
            
            # Get the generated assay_id
            assay_id = cursor.lastrowid
            
            # Generate assay name using the assay_id
            assay_name = f"assay_{assay_id}"
            
            # Create unique assay folder using the new assay name
            assay_folder_path = self._create_assay_folder(docking_assays_base_dir, assay_name)
            
            if not assay_folder_path:
                print("⚠️  Warning: Could not create assay folder, but continuing...")
                assay_folder_path = ""
            
            # Prepare configuration data
            configuration_data = {
                'docking_method': docking_method,
                'table_info': {
                    'name': table_name,
                    'compound_count': table_info['compound_count'] if table_info else 0,
                    'has_sdf_blob': table_info['has_sdf_blob'] if table_info else False,
                    'sdf_blob_count': table_info['sdf_blob_count'] if table_info else 0,
                    'readiness_score': table_info['readiness_score'] if table_info else 0,
                    'docking_status': table_info['docking_status'] if table_info else 'unknown'
                },
                'setup_parameters': {
                    'setup_type': 'configuration_only' if setup_type == 1 else 'ready_for_docking',
                    'created_via': 'dock_table_method'
                },
                'assay_folder': {
                    'path': assay_folder_path,
                    'created': assay_folder_path != ""
                }
            }
            
            # Update the record with the correct assay_name, configuration, and folder path
            cursor.execute('''
                UPDATE docking_assays 
                SET assay_name = ?, configuration = ?, assay_folder_path = ?, 
                    notes = ?, last_modified = CURRENT_TIMESTAMP
                WHERE assay_id = ?
            ''', (
                assay_name,
                json.dumps(configuration_data, indent=2),
                assay_folder_path,
                note_text,
                assay_id
            ))
            
            conn.commit()
            conn.close()
            
            print(f"📝 Assay Registry Created:")
            print(f"   🆔 Assay ID: {assay_id}")
            print(f"   📋 Assay Name: {assay_name}")
            print(f"   📊 Table: {table_name}")
            print(f"   🧬 Method: {docking_method['method_name']} ({docking_method['docking_engine']})")
            print(f"   🔧 Setup Type: {'configuration_only' if setup_type == 1 else 'ready_for_docking'}")
            print(f"   📁 Assay Folder: {assay_folder_path}")
            print(f"   💾 Database: {assays_db_path}")
            
            return {
                'assay_id': assay_id,
                'assay_name': assay_name,
                'project_name': self.name,
                'table_name': table_name,
                'docking_method': docking_method,
                'setup_type': 'configuration_only' if setup_type == 1 else 'ready_for_docking',
                'status': 'created',
                'database_path': assays_db_path,
                'assay_folder_path': assay_folder_path,
                'created_date': datetime.now().isoformat(),
                'notes': note_text
            }
            
        except Exception as e:
            print(f"❌ Error creating docking assay registry: {e}")
            return None

    def _create_assay_folder(self, base_dir: str, assay_name: str) -> str:
        """
        Create a unique folder for the docking assay with proper structure.
        
        Args:
            base_dir (str): Base directory for docking assays
            assay_name (str): Name of the assay
            
        Returns:
            str: Path to the created assay folder, empty string if failed
        """
        try:
            # Clean assay name for folder creation (remove invalid characters)
            clean_name = "".join(c for c in assay_name if c.isalnum() or c in ('_', '-', '.')).strip()
            
            # Create unique folder path
            assay_folder_path = os.path.join(base_dir, clean_name)
            
            # Handle name conflicts by adding suffix
            original_path = assay_folder_path
            counter = 1
            while os.path.exists(assay_folder_path):
                assay_folder_path = f"{original_path}_{counter}"
                counter += 1
            
            # Create the assay folder and subdirectories
            os.makedirs(assay_folder_path, exist_ok=True)
            
            # Create standard subdirectories for docking workflow
            subdirs = [
                'ligands',          # For prepared ligand files (PDBQT, SDF)
                'receptors',        # For receptor files (PDB, PDBQT)
                'grid_files',       # For grid/map files (AutoDock grid files)
                'results',          # For docking output files
                'logs',             # For log files and error reports
                'analysis',         # For analysis results and plots
                'configurations'    # For configuration files (DPF, config files)
            ]
            
            for subdir in subdirs:
                subdir_path = os.path.join(assay_folder_path, subdir)
                os.makedirs(subdir_path, exist_ok=True)
            
            # Create a README file with assay information
            self._create_assay_readme(assay_folder_path, assay_name)
            
            print(f"   📁 Created assay folder: {assay_folder_path}")
            print(f"   📋 Subdirectories: {', '.join(subdirs)}")
            
            return assay_folder_path
            
        except Exception as e:
            print(f"⚠️  Error creating assay folder: {e}")
            return ""

    def _create_assay_readme(self, assay_folder_path: str, assay_name: str) -> None:
        """
        Create a README file in the assay folder with information about the structure.
        
        Args:
            assay_folder_path (str): Path to the assay folder
            assay_name (str): Name of the assay
        """
        try:
            from datetime import datetime
            
            readme_content = f"""# Docking Assay: {assay_name}

    Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    Project: {self.name}

    ## Folder Structure

    - **ligands/**: Prepared ligand files (PDBQT, SDF formats)
    - **receptors/**: Receptor structure files (PDB, PDBQT formats)  
    - **grid_files/**: Grid and map files for docking
    - **results/**: Docking output files and poses
    - **logs/**: Log files and error reports
    - **analysis/**: Analysis results, plots, and summaries
    - **configurations/**: Docking configuration files (DPF, config files)

    ## Usage Notes

    This folder structure is designed to organize all files related to this docking assay.
    Place input files in the appropriate directories and docking results will be saved
    to the results/ directory.

    ## File Naming Convention

    - Use consistent naming with assay ID prefix
    - Keep original compound identifiers where possible
    - Use standard file extensions (.pdbqt, .sdf, .pdb, .dlg, etc.)

    Generated by TidyScreen MolDock module.
    """
            
            readme_path = os.path.join(assay_folder_path, 'README.md')
            with open(readme_path, 'w') as f:
                f.write(readme_content)
                
        except Exception as e:
            print(f"⚠️  Could not create README file: {e}")

    def _update_assay_registry_compound_count(self, assay_id: int, compound_count: int) -> bool:
        """
        Update the compound count in the assay registry after SMILES are loaded.
        
        Args:
            assay_id (int): Assay registry ID
            compound_count (int): Actual number of loaded compounds
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            assays_db_path = os.path.join(os.path.dirname(self.__docking_registers_db), 'docking_assays.db')
            
            import sqlite3
            conn = sqlite3.connect(assays_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                UPDATE docking_assays 
                SET compound_count = ?, status = 'molecules_loaded', last_modified = CURRENT_TIMESTAMP
                WHERE assay_id = ?
            ''', (compound_count, assay_id))
            
            conn.commit()
            conn.close()
            
            print(f"   📊 Updated assay registry: {compound_count:,} compounds loaded")
            return True
            
        except Exception as e:
            print(f"⚠️  Error updating assay registry compound count: {e}")
            return False

    def _update_assay_registry_status(self, assay_id: int, status: str, notes: str = "") -> bool:
        """
        Update the status of an assay registry entry.
        
        Args:
            assay_id (int): Assay registry ID
            status (str): New status ('created', 'molecules_loaded', 'failed', 'completed')
            notes (str): Additional notes about the status change
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            assays_db_path = os.path.join(os.path.dirname(self.__docking_registers_db), 'docking_assays.db')
            
            import sqlite3
            conn = sqlite3.connect(assays_db_path)
            cursor = conn.cursor()
            
            # Get existing notes
            cursor.execute("SELECT notes FROM docking_assays WHERE assay_id = ?", (assay_id,))
            result = cursor.fetchone()
            existing_notes = result[0] if result else ""
            
            # Append new notes
            updated_notes = f"{existing_notes}\n{notes}" if existing_notes else notes
            
            cursor.execute('''
                UPDATE docking_assays 
                SET status = ?, notes = ?, last_modified = CURRENT_TIMESTAMP
                WHERE assay_id = ?
            ''', (status, updated_notes.strip(), assay_id))
            
            conn.commit()
            conn.close()
            
            print(f"   📝 Updated assay registry status: {status}")
            return True
            
        except Exception as e:
            print(f"⚠️  Error updating assay registry status: {e}")
            return False

    def get_assay_folder_path(self, assay_id: int) -> Optional[str]:
        """
        Retrieve the folder path for a specific assay ID.
        
        Args:
            assay_id (int): Assay registry ID
            
        Returns:
            Optional[str]: Path to assay folder or None if not found
        """
        try:
            assays_db_path = os.path.join(os.path.dirname(self.__docking_registers_db), 'docking_assays.db')
            
            if not os.path.exists(assays_db_path):
                return None
            
            import sqlite3
            conn = sqlite3.connect(assays_db_path)
            cursor = conn.cursor()
            
            cursor.execute("SELECT assay_folder_path FROM docking_assays WHERE assay_id = ?", (assay_id,))
            result = cursor.fetchone()
            
            conn.close()
            
            if result and result[0]:
                folder_path = result[0]
                # Verify folder still exists
                if os.path.exists(folder_path):
                    return folder_path
                else:
                    print(f"⚠️  Assay folder not found: {folder_path}")
                    return None
            
            return None
            
        except Exception as e:
            print(f"❌ Error retrieving assay folder path: {e}")
            return None

    def list_assay_folders(self) -> Optional[List[Dict[str, Any]]]:
        """
        List all docking assays and their folder paths for this project.
        
        Returns:
            Optional[List[Dict]]: List of assay information with folder paths
        """
        try:
            assays_db_path = os.path.join(os.path.dirname(self.__docking_registers_db), 'docking_assays.db')
            
            if not os.path.exists(assays_db_path):
                print("📋 No docking assays found")
                return None
            
            import sqlite3
            conn = sqlite3.connect(assays_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT assay_id, assay_name, table_name, docking_method_name, 
                    docking_engine, status, created_date, assay_folder_path, notes
                FROM docking_assays 
                WHERE project_name = ?
                ORDER BY created_date ASC
            ''', (self.name,))
            
            results = cursor.fetchall()
            conn.close()
            
            if not results:
                print("📋 No docking assays found for this project")
                return None
            
            assay_list = []
            print(f"\n📋 DOCKING ASSAYS FOR PROJECT: {self.name}")
            print("=" * 130)
            print(f"{'ID':<5} {'Assay Name':<25} {'Table':<20} {'Method':<15} {'Notes':<60}")
            print("=" * 130)
            
            for row in results:
                assay_id, assay_name, table_name, method_name, engine, status, created_date, folder_path, notes = row
                
                folder_exists = "✅ Yes" if folder_path and os.path.exists(folder_path) else "❌ No"
                note_preview = (notes or "").replace('\n', ' ')[:60]
                
                print(f"{assay_id:<5} {assay_name[:24]:<25} {table_name[:14]:<20} "
                    f"{method_name[:14]:<15} {note_preview:<60}")
                
                assay_list.append({
                    'assay_id': assay_id,
                    'assay_name': assay_name,
                    'table_name': table_name,
                    'docking_method_name': method_name,
                    'docking_engine': engine,
                    'status': status,
                    'created_date': created_date,
                    'assay_folder_path': folder_path,
                    'folder_exists': folder_path and os.path.exists(folder_path),
                    'notes': notes
                })
            
            print("=" * 130)
            print(f"Total assays: {len(assay_list)}")
            
            return assay_list
            
        except Exception as e:
            print(f"❌ Error listing assay folders: {e}")
            return None

    def _update_assay_registry_compound_count(self, assay_id: int, compound_count: int) -> bool:
        """
        Update the compound count in the assay registry after SMILES are loaded.
        
        Args:
            assay_id (int): Assay registry ID
            compound_count (int): Actual number of loaded compounds
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            assays_db_path = os.path.join(os.path.dirname(self.__docking_registers_db), 'docking_assays.db')
            
            import sqlite3
            conn = sqlite3.connect(assays_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                UPDATE docking_assays 
                SET compound_count = ?, status = 'molecules_loaded', last_modified = CURRENT_TIMESTAMP
                WHERE assay_id = ?
            ''', (compound_count, assay_id))
            
            conn.commit()
            conn.close()
            
            print(f"   📊 Updated assay registry: {compound_count:,} compounds loaded")
            return True
            
        except Exception as e:
            print(f"⚠️  Error updating assay registry compound count: {e}")
            return False

    def _select_docking_method(self) -> Optional[Dict[str, Any]]:
        """
        Select a docking method from available methods created by create_docking_method.
        
        Returns:
            Optional[Dict[str, Any]]: Selected docking method information or None if cancelled
        """
        try:
            print(f"🎯 SELECT DOCKING METHOD")
            print("=" * 50)
            
            # Check if docking methods database exists
            methods_db_path = os.path.join(os.path.dirname(self.__docking_registers_db), 'docking_methods.db')
            
            if not os.path.exists(methods_db_path):
                print(f"❌ No docking methods found")
                print(f"   Database not found: {methods_db_path}")
                print(f"   💡 Create a docking method first using create_docking_method()")
                return None
            
            # Connect to docking methods database
            try:
                import sqlite3
                import json
                conn = sqlite3.connect(methods_db_path)
                cursor = conn.cursor()
                
                # Get available docking methods, including ligand_prep_params
                cursor.execute('''
                    SELECT id, method_name, docking_engine, description, parameters, ligand_prep_params, created_date
                    FROM docking_methods
                    ORDER BY created_date ASC
                ''')
                methods = cursor.fetchall()
                conn.close()
                
            except Exception as e:
                print(f"❌ Error accessing docking methods database: {e}")
                return None
            
            if not methods:
                print("❌ No docking methods available")
                print("   💡 Create a docking method first using create_docking_method()")
                return None
            
            print(f"📋 Available Docking Methods ({len(methods)} found):")
            print("-" * 80)
            print(f"{'#':<3} {'Method ID':<10} {'Method Name':<20} {'Engine':<15} {'Description':<30} {'Created':<12}")
            print("-" * 80)
            
            # Display available methods
            method_list = []
            for i, method in enumerate(methods, 1):
                method_id, method_name, docking_engine, description, parameters_json, ligand_prep_params_json, created_date = method
                # Parse parameters for display
                try:
                    parameters = json.loads(parameters_json) if parameters_json else {}
                except:
                    parameters = {}
                try:
                    ligand_prep_parameters = json.loads(ligand_prep_params_json) if ligand_prep_params_json else {}
                except:
                    ligand_prep_parameters = {}
                # Format created date
                try:
                    from datetime import datetime
                    created_dt = datetime.strptime(created_date, '%Y-%m-%d %H:%M:%S')
                    created_str = created_dt.strftime('%Y-%m-%d')
                except:
                    created_str = created_date[:10] if created_date else 'Unknown'
                # Truncate description for display
                desc_short = (description[:27] + '...') if len(description) > 30 else description
                print(f"{i:<3}{method_id:>10} {method_name[:19]:>20} {docking_engine[:14]:>15} {desc_short:>30} {created_str:>12}")
                method_list.append({
                    'method_id': method_id,
                    'method_name': method_name,
                    'docking_engine': docking_engine,
                    'description': description,
                    'parameters': parameters,
                    'ligand_prep_parameters': ligand_prep_parameters,
                    'created_date': created_date
                })
            
            print("-" * 80)
            print("Commands: Enter method number, method name, 'details <num>' for info, or 'cancel'")
            
            # Method selection loop
            while True:
                try:
                    selection = input(f"\n🎯 Select docking method: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    
                    # Handle details command
                    if selection.lower().startswith('details '):
                        try:
                            detail_idx = int(selection.split()[1]) - 1
                            if 0 <= detail_idx < len(method_list):
                                self._show_docking_method_details(method_list[detail_idx])
                                continue
                            else:
                                print(f"❌ Invalid method number. Please enter 1-{len(method_list)}")
                                continue
                        except (IndexError, ValueError):
                            print("❌ Invalid details command. Use 'details <number>'")
                            continue
                    
                    # Try as number first
                    try:
                        method_idx = int(selection) - 1
                        if 0 <= method_idx < len(method_list):
                            selected_method = method_list[method_idx]
                            print(f"\n✅ Selected: '{selected_method['method_name']}' using {selected_method['docking_engine']}")
                            return selected_method
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(method_list)}")
                            continue
                            
                    except ValueError:
                        # Try as method name
                        matching_methods = [m for m in method_list if m['method_name'].lower() == selection.lower()]
                        if matching_methods:
                            selected_method = matching_methods[0]
                            print(f"\n✅ Selected: '{selected_method['method_name']}' using {selected_method['docking_engine']}")
                            return selected_method
                        else:
                            print(f"❌ Method '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Method selection cancelled")
                    return None
            
        except Exception as e:
            print(f"❌ Error in docking method selection: {e}")
            return None

    def _show_docking_method_details(self, method_info: Dict[str, Any]) -> None:
        """
        Show detailed information about a docking method.
        
        Args:
            method_info (Dict): Docking method information dictionary
        """
        try:
            print(f"\n🎯 DOCKING METHOD DETAILS: '{method_info['method_name']}'")
            print("=" * 60)
            print(f"🆔 Method ID: {method_info['method_id']}")
            print(f"🏷️  Method Name: {method_info['method_name']}")
            print(f"🧬 Docking Engine: {method_info['docking_engine']}")
            print(f"📝 Description: {method_info['description']}")
            print(f"📅 Created: {method_info['created_date']}")
            
            print(f"\n⚙️  CONFIGURED DOCKING PARAMETERS:")
            parameters = method_info.get('parameters', {})
            if parameters:
                for key, value in parameters.items():
                    param_name = key.replace('_', ' ').title()
                    print(f"   • {param_name}: {value}")
            else:
                print("   • No specific docking parameters configured")
                
            print(f"\n⚙️  CONFIGURED LIGANDS PREP PARAMETERS:")
            lig_parameters = method_info.get('ligand_prep_parameters', {})
            if lig_parameters:
                for key, value in lig_parameters.items():
                    param_name = key.replace('_', ' ').title()
                    print(f"   • {param_name}: {value}")
            else:
                print("   • No specific ligand preparation parameters configured")
            
            # Show engine-specific information
            engine = method_info['docking_engine']
            print(f"\n🔧 {engine.upper()} ENGINE INFO:")
            
            if engine == 'AutoDockGPU':
                print("   • Optimized for GPU acceleration")
                print("   • Best for high-throughput screening")
                print("   • Requires CUDA-compatible GPU")
            elif engine == 'AutoDock':
                print("   • Classical CPU-based docking")
                print("   • Highly configurable parameters")
                print("   • Suitable for detailed studies")
            elif engine == 'Vina':
                print("   • Fast and accurate docking")
                print("   • Easy to use and configure")
                print("   • Good balance of speed and accuracy")
            
            print("=" * 60)
            
        except Exception as e:
            print(f"❌ Error showing method details: {e}")

    def _show_docking_setup_summary(self, selected_table: str, selected_method: Dict[str, Any], 
                                table_data: List[Dict]) -> Optional[int]:
        """
        Show summary of selected table and docking method, and get user choice for next steps.
        
        Args:
            selected_table (str): Name of selected table
            selected_method (Dict): Selected docking method information
            table_data (List[Dict]): Table data list
            
        Returns:
            Optional[int]: User choice (1=configuration only, 2=load SMILES, None=cancel)
        """
        try:
            # Find selected table info
            table_info = next((t for t in table_data if t['name'] == selected_table), None)
            
            if not table_info:
                return None
            
            print(f"\n🎯 DOCKING SETUP SUMMARY")
            print("=" * 60)
            print(f"📋 Selected Table: '{selected_table}'")
            print(f"📊 Compounds: {table_info['compound_count']:,}")
            print(f"🎯 Table Status: {table_info['docking_readiness']['status']}")
            
            if table_info['has_sdf_blob'] and table_info['sdf_blob_count'] > 0:
                print(f"🚀 3D Molecules Ready: {table_info['sdf_blob_count']:,}")
                prep_coverage = (table_info['sdf_blob_count'] / table_info['compound_count']) * 100
                print(f"📈 Preparation Coverage: {prep_coverage:.1f}%")
            
            print(f"\n🔬 Selected Docking Method: '{selected_method['method_name']}'")
            print(f"🧬 Docking Engine: {selected_method['docking_engine']}")
            print(f"📝 Description: {selected_method['description']}")
            
            # Show key parameters
            parameters = selected_method.get('parameters', {})
            if parameters:
                print(f"⚙️  Key Parameters:")
                # Show most important parameters based on engine
                engine = selected_method['docking_engine']
                if engine == 'AutoDockGPU':
                    key_params = ['nruns', 'psize']
                elif engine == 'AutoDock':
                    key_params = ['num_runs', 'population_size', 'energy_evaluations']
                elif engine == 'Vina':
                    key_params = ['exhaustiveness', 'num_modes', 'cpu_cores']
                else:
                    key_params = list(parameters.keys())[:3]  # First 3 parameters
                
                for param in key_params:
                    if param in parameters:
                        param_name = param.replace('_', ' ').title()
                        print(f"   • {param_name}: {parameters[param]}")
            
            # Show readiness assessment
            readiness = table_info['docking_readiness']
            if readiness['score'] >= 70:
                print(f"\n✅ Setup Status: Ready for docking")
            else:
                print(f"\n⚠️  Setup Status: {readiness['status']}")
                if readiness['issues']:
                    print("   Issues to consider:")
                    for issue in readiness['issues'][:3]:  # Show first 3 issues
                        print(f"   • {issue}")
            
            print(f"\n🔄 NEXT STEPS:")
            print("   1. Load molecules and prepare for immediate docking")
            print("   2. Cancel setup")
            print("=" * 60)
            
            # Get user choice
            while True:
                try:
                    choice = input("🎯 Select next step (1/2): ").strip()
                    
                    if choice == '1':
                        print("🧪 Will load molecules for immediate docking preparation")
                        return 1
                    elif choice == '2' or choice.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Docking setup cancelled")
                        return None
                    else:
                        print("❌ Invalid choice. Please enter 1 or 2")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ Docking setup cancelled")
                    return None
            
        except Exception as e:
            print(f"⚠️  Error showing docking setup summary: {e}")
            return None

    def _show_docking_preparation_summary(self, selected_table: str, table_data: List[Dict]) -> Optional[int]:
        """
        Show summary of selected table for docking preparation and get user choice.
        
        Args:
            selected_table (str): Name of selected table
            table_data (List[Dict]): Table data list
            
        Returns:
            Optional[int]: User choice (1=table name, 2=SMILES list, None=cancel)
        """
        try:
            # Find selected table info
            table_info = next((t for t in table_data if t['name'] == selected_table), None)
            
            if not table_info:
                return None
            
            print(f"\n🎯 DOCKING PREPARATION SUMMARY")
            print("-" * 50)
            print(f"📋 Selected Table: '{selected_table}'")
            print(f"📊 Compounds: {table_info['compound_count']:,}")
            print(f"🎯 Readiness: {table_info['docking_readiness']['status']}")
            
            readiness = table_info['docking_readiness']
            if readiness['score'] >= 70:
                print("✅ Table is ready for docking preparation")
            elif readiness['score'] >= 50:
                print("🟡 Table has minor issues but can be used")
            else:
                print("⚠️  Table has significant issues for docking")
            
            print(f"\n🔄 PREPARATION OPTIONS:")
            print("   1. Return table name (for later processing)")
            print("   2. Load SMILES immediately (for molecule preparation)")
            print("   3. Cancel selection")
            print("-" * 50)
            
            # Get user choice
            while True:
                try:
                    choice = input("🧬 Select preparation option (1/2/3): ").strip()
                    
                    if choice == '1':
                        print("✅ Will return table name for later processing")
                        return 1
                    elif choice == '2':
                        print("🧪 Will load SMILES for immediate molecule preparation")
                        return 2
                    elif choice == '3' or choice.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Preparation cancelled")
                        return None
                    else:
                        print("❌ Invalid choice. Please enter 1, 2, or 3")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ Preparation cancelled")
                    return None
            
        except Exception as e:
            print(f"⚠️  Error showing preparation summary: {e}")
            return None

    def _get_smiles_from_table(self, table_name: str) -> Optional[List[Dict[str, Any]]]:
        """
        Extract SMILES and associated data from the selected table.
        
        Args:
            table_name (str): Name of the table to extract SMILES from
            
        Returns:
            Optional[List[Dict]]: List of compound dictionaries with SMILES and metadata
        """
        try:
            import sqlite3
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get table schema to identify available columns
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [row[1] for row in cursor.fetchall()]
            
            # Check for required and optional columns
            has_smiles = 'smiles' in [col.lower() for col in columns]
            if not has_smiles:
                print("❌ No SMILES column found in table")
                conn.close()
                return None
            
            # Build query to get available data
            base_columns = ['smiles']
            optional_columns = ['name', 'id', 'inchi_key', 'molecular_weight', 'logp']
            
            available_columns = []
            for col in columns:
                col_lower = col.lower()
                if col_lower == 'smiles':
                    available_columns.append('smiles')
                elif col_lower in ['name', 'compound_name', 'title']:
                    available_columns.append(f"{col} as name")
                elif col_lower in ['id', 'compound_id', 'mol_id']:
                    available_columns.append(f"{col} as id")
                elif col_lower in ['inchi_key', 'inchikey', 'inchi']:
                    available_columns.append(f"{col} as inchi_key")
                elif col_lower in ['molecular_weight', 'mw', 'mol_weight']:
                    available_columns.append(f"{col} as molecular_weight")
                elif col_lower in ['logp', 'clogp', 'alogp']:
                    available_columns.append(f"{col} as logp")
            
            # Execute query
            query = f"SELECT {', '.join(available_columns)} FROM {table_name} WHERE smiles IS NOT NULL AND smiles != ''"
            
            if TQDM_AVAILABLE:
                # First get count for progress bar
                cursor.execute(f"SELECT COUNT(*) FROM {table_name} WHERE smiles IS NOT NULL AND smiles != ''")
                total_count = cursor.fetchone()[0]
                
                cursor.execute(query)
                results = cursor.fetchall()
                
                progress_bar = tqdm(total=len(results), desc="Loading SMILES", unit="compounds")
            else:
                cursor.execute(query)
                results = cursor.fetchall()
                progress_bar = None
                print(f"📊 Processing {len(results)} compounds...")
            
            smiles_data = []
            processed_count = 0
            
            for row in results:
                try:
                    # Create compound dictionary
                    compound_dict = {}
                    
                    for i, col_info in enumerate(available_columns):
                        if ' as ' in col_info:
                            col_name = col_info.split(' as ')[1]
                        else:
                            col_name = col_info
                        
                        compound_dict[col_name] = row[i] if i < len(row) else None
                    
                    # Validate SMILES
                    smiles = compound_dict.get('smiles', '').strip()
                    if smiles and len(smiles) > 0:
                        # Add metadata for docking preparation
                        compound_dict['conf_rank'] = 0  # Default conformer rank
                        compound_dict['source_table'] = table_name
                        compound_dict['compound_index'] = processed_count
                        
                        smiles_data.append(compound_dict)
                        processed_count += 1
                    
                    if progress_bar:
                        progress_bar.update(1)
                    elif processed_count % 1000 == 0:
                        print(f"   📊 Processed: {processed_count:,} compounds")
                    
                except Exception as e:
                    print(f"   ⚠️  Error processing compound at index {len(smiles_data)}: {e}")
                    continue
            
            if progress_bar:
                progress_bar.close()
            
            conn.close()
            
            print(f"\n✅ Successfully loaded {len(smiles_data)} valid compounds")
            print(f"   📋 Table: {table_name}")
            print(f"   🧪 Valid SMILES: {len(smiles_data):,}")
            
            if len(smiles_data) == 0:
                print("❌ No valid SMILES found in table")
                return None
            
            return smiles_data
            
        except Exception as e:
            print(f"❌ Error loading SMILES from table: {e}")
            return None

    def _get_table_info_for_docking(self, table_name: str) -> Dict[str, Any]:
        """
        Get table information specifically for docking preparation.
        Reuses existing database access patterns.
        
        Args:
            table_name (str): Name of the table to analyze
            
        Returns:
            Dict[str, Any]: Table information dictionary
        """
        try:
            import sqlite3
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get basic table info
            cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            compound_count = cursor.fetchone()[0]
            
            # Check table schema for required columns
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [row[1] for row in cursor.fetchall()]
            has_smiles = 'smiles' in [col.lower() for col in columns]
            has_sdf_blob = 'sdf_blob' in [col.lower() for col in columns]
            
            # Count SDF blobs if column exists
            sdf_blob_count = 0
            if has_sdf_blob:
                cursor.execute(f"SELECT COUNT(*) FROM {table_name} WHERE sdf_blob IS NOT NULL")
                sdf_blob_count = cursor.fetchone()[0]
            
            # Get sample data to assess docking suitability
            if has_smiles and compound_count > 0:
                cursor.execute(f"SELECT smiles FROM {table_name} WHERE smiles IS NOT NULL LIMIT 5")
                sample_smiles = [row[0] for row in cursor.fetchall() if row[0] and row[0].strip()]
                valid_smiles_count = len(sample_smiles)
            else:
                valid_smiles_count = 0
                sample_smiles = []
            
            conn.close()
            
            # Classify table type (reusing existing logic from ChemSpace)
            table_type = self._classify_table_type_for_docking(table_name)
            
            # Assess docking readiness with SDF blob information
            docking_readiness = self._assess_docking_readiness(
                compound_count, has_smiles, valid_smiles_count, sample_smiles,
                has_sdf_blob, sdf_blob_count
            )
            
            return {
                'name': table_name,
                'type': table_type,
                'compound_count': compound_count,
                'has_smiles': has_smiles,
                'has_sdf_blob': has_sdf_blob,
                'sdf_blob_count': sdf_blob_count,
                'valid_smiles_count': valid_smiles_count,
                'docking_readiness': docking_readiness,
                'readiness_score': docking_readiness['score'],
                'docking_status': docking_readiness['docking_status'],
                'icon': self._get_table_icon_for_docking(table_type, docking_readiness['score'], has_sdf_blob, sdf_blob_count)
            }
            
        except Exception as e:
            return {
                'name': table_name,
                'type': 'error',
                'compound_count': 0,
                'has_smiles': False,
                'has_sdf_blob': False,
                'sdf_blob_count': 0,
                'valid_smiles_count': 0,
                'docking_readiness': {'score': 0, 'status': 'error', 'docking_status': 'NOT ready for docking', 'issues': [str(e)]},
                'readiness_score': 0,
                'docking_status': 'NOT ready for docking',
                'icon': '❌'
            }

    def _classify_table_type_for_docking(self, table_name: str) -> str:
        """
        Classify table type for docking purposes.
        Reuses existing classification logic from ChemSpace.
        
        Args:
            table_name (str): Name of the table
            
        Returns:
            str: Table type classification
        """
        name_lower = table_name.lower()
        
        # Reuse ChemSpace classification patterns
        if '_products_' in name_lower or 'reaction' in name_lower:
            return 'reaction_products'
        elif '_filtered_' in name_lower or 'filter' in name_lower:
            return 'filtered_compounds'
        elif 'stereoisomer' in name_lower or 'stereo' in name_lower:
            return 'stereoisomers'
        elif 'ersilia' in name_lower or 'prediction' in name_lower:
            return 'predicted_compounds'
        elif name_lower.startswith('temp_') or name_lower.startswith('tmp_'):
            return 'temporary'
        elif any(word in name_lower for word in ['workflow', 'batch', 'processed']):
            return 'processed'
        else:
            return 'original'

    def _assess_docking_readiness(self, compound_count: int, has_smiles: bool, 
                                valid_smiles_count: int, sample_smiles: List[str],
                                has_sdf_blob: bool = False, sdf_blob_count: int = 0) -> Dict[str, Any]:
        """
        Assess how ready a table is for docking preparation.
        
        Args:
            compound_count (int): Total number of compounds
            has_smiles (bool): Whether table has SMILES column
            valid_smiles_count (int): Number of valid SMILES in sample
            sample_smiles (List[str]): Sample SMILES strings
            has_sdf_blob (bool): Whether table has SDF blob column
            sdf_blob_count (int): Number of compounds with SDF blobs
            
        Returns:
            Dict[str, Any]: Readiness assessment
        """
        try:
            issues = []
            score = 100  # Start with perfect score and deduct points
            docking_status = "NOT ready for docking"  # Default status
            
            # Check basic requirements
            if compound_count == 0:
                issues.append("Table is empty")
                score = 0
            elif compound_count < 5:
                issues.append(f"Very few compounds ({compound_count})")
                score -= 20
            
            if not has_smiles:
                issues.append("No SMILES column found")
                score = 0
            elif valid_smiles_count == 0:
                issues.append("No valid SMILES found in sample")
                score -= 50
            elif valid_smiles_count < len(sample_smiles):
                issues.append("Some invalid SMILES detected")
                score -= 10
            
            # Check for SDF blob column (critical for docking readiness)
            if has_sdf_blob:
                if sdf_blob_count > 0:
                    docking_status = "ready for docking"
                    print(f"   ✅ Found {sdf_blob_count:,} prepared molecules with 3D coordinates")
                else:
                    issues.append("SDF blob column exists but no prepared molecules found")
                    docking_status = "NOT ready for docking"
                    score -= 30
            else:
                issues.append("No SDF blob column - molecules need 3D preparation")
                docking_status = "NOT ready for docking"
                score -= 40
            
            # Check SMILES complexity for docking suitability (only if we have valid SMILES)
            if sample_smiles and score > 0:
                complexity_assessment = self._assess_smiles_complexity(sample_smiles)
                if complexity_assessment['has_issues']:
                    issues.extend(complexity_assessment['issues'])
                    score -= complexity_assessment['score_penalty']
            
            # Determine overall status with docking-specific criteria
            if has_sdf_blob and sdf_blob_count > 0:
                if score >= 90:
                    status = 'excellent'
                elif score >= 70:
                    status = 'good'
                elif score >= 50:
                    status = 'fair'
                else:
                    status = 'poor'
            else:
                # Without SDF blobs, maximum status is 'needs_preparation'
                if score >= 70:
                    status = 'needs_preparation'
                elif score >= 50:
                    status = 'fair'
                elif score > 0:
                    status = 'poor'
                else:
                    status = 'unsuitable'
            
            return {
                'score': max(0, score),
                'status': status,
                'docking_status': docking_status,
                'issues': issues,
                'compound_count': compound_count,
                'has_smiles': has_smiles,
                'has_sdf_blob': has_sdf_blob,
                'sdf_blob_count': sdf_blob_count,
                'valid_smiles_sample': valid_smiles_count
            }
            
        except Exception as e:
            return {
                'score': 0,
                'status': 'error',
                'docking_status': 'NOT ready for docking',
                'issues': [f"Assessment error: {e}"],
                'compound_count': compound_count,
                'has_smiles': has_smiles,
                'has_sdf_blob': False,
                'sdf_blob_count': 0,
                'valid_smiles_sample': 0
            }

    def _assess_smiles_complexity(self, sample_smiles: List[str]) -> Dict[str, Any]:
        """
        Assess SMILES complexity for docking suitability.
        
        Args:
            sample_smiles (List[str]): Sample SMILES strings
            
        Returns:
            Dict[str, Any]: Complexity assessment
        """
        try:
            issues = []
            score_penalty = 0
            
            for smiles in sample_smiles[:3]:  # Check first 3 samples
                # Check for very simple molecules (might not be drug-like)
                if len(smiles) < 10:
                    issues.append("Very simple molecules detected")
                    score_penalty += 5
                
                # Check for very complex molecules (might be problematic for docking)
                elif len(smiles) > 200:
                    issues.append("Very complex molecules detected")
                    score_penalty += 10
                
                # Check for inorganic/problematic patterns
                problematic_patterns = ['[Na+]', '[K+]', '[Cl-]', '[SO4-2]', '[PO4-3]']
                if any(pattern in smiles for pattern in problematic_patterns):
                    issues.append("Inorganic/salt components detected")
                    score_penalty += 15
            
            return {
                'has_issues': len(issues) > 0,
                'issues': list(set(issues)),  # Remove duplicates
                'score_penalty': min(score_penalty, 30)  # Cap penalty
            }
            
        except Exception:
            return {
                'has_issues': True,
                'issues': ["Error assessing SMILES complexity"],
                'score_penalty': 10
            }

    def _get_table_icon_for_docking(self, table_type: str, readiness_score: int, 
                                has_sdf_blob: bool = False, sdf_blob_count: int = 0) -> str:
        """
        Get appropriate icon for table based on type, docking readiness, and SDF blob availability.
        
        Args:
            table_type (str): Table type
            readiness_score (int): Docking readiness score
            has_sdf_blob (bool): Whether table has SDF blob column
            sdf_blob_count (int): Number of compounds with SDF blobs
            
        Returns:
            str: Emoji icon
        """
        # Base icons by type (reusing ChemSpace patterns)
        type_icons = {
            'original': '📋',
            'reaction_products': '🧪',
            'filtered_compounds': '🔍',
            'stereoisomers': '🔄',
            'predicted_compounds': '🧠',
            'processed': '⚗️',
            'temporary': '📄',
            'error': '❌'
        }
        
        base_icon = type_icons.get(table_type, '📊')
        
        # Docking readiness indicators with SDF blob consideration
        if has_sdf_blob and sdf_blob_count > 0:
            # Has prepared 3D molecules - ready for docking
            if readiness_score >= 90:
                return f"{base_icon}🚀"  # Ready to dock - excellent
            elif readiness_score >= 70:
                return f"{base_icon}✅"  # Ready to dock - good
            else:
                return f"{base_icon}🟢"  # Ready to dock - fair quality
        else:
            # No 3D molecules prepared - needs preparation
            if readiness_score >= 70:
                return f"{base_icon}🔧"  # Needs 3D preparation
            elif readiness_score >= 50:
                return f"{base_icon}🟡"  # Needs preparation + has issues
            elif readiness_score > 0:
                return f"{base_icon}🟠"  # Poor quality, needs work
            else:
                return f"{base_icon}🔴"  # Unsuitable for docking

    def _sort_table_data_for_docking(self, table_data: List[Dict], sort_by: str) -> List[Dict]:
        """
        Sort table data for docking display.
        Reuses existing sorting logic with docking-specific priorities.
        
        Args:
            table_data (List[Dict]): List of table information dictionaries
            sort_by (str): Sort criteria
            
        Returns:
            List[Dict]: Sorted table data
        """
        try:
            if sort_by == 'readiness':
                # Sort by docking readiness score (highest first)
                return sorted(table_data, key=lambda x: x['readiness_score'], reverse=True)
            elif sort_by == 'compounds':
                return sorted(table_data, key=lambda x: x['compound_count'], reverse=True)
            elif sort_by == 'type':
                return sorted(table_data, key=lambda x: (x['type'], x['name'].lower()))
            else:  # sort by name (default)
                return sorted(table_data, key=lambda x: x['name'].lower())
        except Exception:
            return table_data

    def _display_and_select_docking_table(self, table_data: List[Dict], 
                                        filter_type: Optional[str],
                                        show_all_tables: bool = False,
                                        ready_count: int = 0,
                                        not_ready_count: int = 0) -> Optional[str]:
        """
        Display table list and handle user selection for docking.
        
        Args:
            table_data (List[Dict]): Table information data
            filter_type (Optional[str]): Applied filter type
            show_all_tables (bool): Whether showing all tables or only ready ones
            ready_count (int): Total number of ready tables
            not_ready_count (int): Total number of not ready tables
            
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        try:
            # Enhanced header showing filtering context
            if show_all_tables:
                print(f"📊 Displaying ALL tables: {len(table_data)} (Ready: {ready_count}, Not Ready: {not_ready_count})")
                if filter_type:
                    print(f"🔍 Additionally filtered by type: '{filter_type}'")
            else:
                print(f"✅ Displaying READY tables only: {len(table_data)} of {ready_count + not_ready_count} total")
                if filter_type:
                    print(f"🔍 Additionally filtered by type: '{filter_type}'")
            
            print("-" * 90)
            print(f"{'#':<3} {'Icon':<6} {'Table Name':<25} {'Type':<15} {'Compounds':<10} {'3D Mols':<8} {'Docking Status':<15}")
            print("-" * 90)
            
            for i, table_info in enumerate(table_data, 1):
                compounds_str = f"{table_info['compound_count']:,}"
                sdf_count_str = f"{table_info['sdf_blob_count']:,}" if table_info['has_sdf_blob'] else "0"
                docking_status = "✅ Ready" if table_info['docking_status'] == 'ready for docking' else "❌ Not Ready"
                
                print(f"{i:<3} {table_info['icon']:<6} {table_info['name'][:24]:<25} "
                    f"{table_info['type'][:14]:<15} {compounds_str:<10} {sdf_count_str:<8} {docking_status:<15}")
            
            print("-" * 90)
            
            # Show enhanced readiness legend
            print("📊 Icons: 🚀✅Ready to dock  🔧Needs 3D prep  🟡🟠Issues  🔴Unsuitable")
            print("💡 3D Mols: Number of molecules with prepared 3D coordinates (SDF blobs)")
            
            # Show filtering status and options
            if not show_all_tables and not_ready_count > 0:
                print(f"ℹ️  Note: {not_ready_count} not-ready tables are hidden. Use show_all_tables=True to see them.")
            
            print("\nCommands: Enter table number, table name, 'details <num>' for info, or 'cancel'")
            
            # Rest of the selection logic remains the same...
            while True:
                try:
                    selection = input(f"\n🧬 Select table for docking preparation: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    
                    # Handle details command
                    if selection.lower().startswith('details '):
                        try:
                            detail_idx = int(selection.split()[1]) - 1
                            if 0 <= detail_idx < len(table_data):
                                self._show_table_details_for_docking(table_data[detail_idx])
                                continue
                            else:
                                print(f"❌ Invalid table number. Please enter 1-{len(table_data)}")
                                continue
                        except (IndexError, ValueError):
                            print("❌ Invalid details command. Use 'details <number>'")
                            continue
                    
                    # Try as number first
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(table_data):
                            selected_info = table_data[table_idx]
                            
                            # Enhanced warning for tables not ready for docking
                            if selected_info['docking_status'] == 'NOT ready for docking':
                                print(f"\n⚠️  WARNING: Table '{selected_info['name']}' is NOT ready for docking")
                                if not selected_info['has_sdf_blob']:
                                    print("   🔧 Missing: 3D molecular coordinates (SDF blobs)")
                                    print("   💡 Suggestion: Use ChemSpace.generate_mols_in_table() to prepare 3D molecules")
                                elif selected_info['sdf_blob_count'] == 0:
                                    print("   📊 SDF blob column exists but no prepared molecules found")
                                
                                print("Issues:")
                                for issue in selected_info['docking_readiness']['issues']:
                                    print(f"   • {issue}")
                                
                                confirm = input("\nContinue with this table anyway? (y/n): ").strip().lower()
                                if confirm not in ['y', 'yes']:
                                    continue
                            else:
                                print(f"\n✅ Table '{selected_info['name']}' is ready for docking!")
                                print(f"   🚀 {selected_info['sdf_blob_count']:,} molecules with 3D coordinates available")
                            
                            return selected_info['name']
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(table_data)}")
                            continue
                            
                    except ValueError:
                        # Try as table name
                        matching_tables = [t for t in table_data if t['name'].lower() == selection.lower()]
                        if matching_tables:
                            return matching_tables[0]['name']
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error in table selection: {e}")
            return None

    def _show_table_details_for_docking(self, table_info: Dict[str, Any]) -> None:
        """
        Show detailed information about a table for docking preparation.
        
        Args:
            table_info (Dict): Table information dictionary
        """
        try:
            print(f"\n📋 DOCKING TABLE DETAILS: '{table_info['name']}'")
            print("=" * 60)
            print(f"🏷️  Type: {table_info['type']}")
            print(f"📊 Total Compounds: {table_info['compound_count']:,}")
            print(f"🧪 Has SMILES: {'Yes' if table_info['has_smiles'] else 'No'}")
            print(f"✅ Valid SMILES: {table_info['valid_smiles_count']} (sample)")
            
            # SDF blob information
            print(f"🔬 Has 3D Coordinates: {'Yes' if table_info['has_sdf_blob'] else 'No'}")
            if table_info['has_sdf_blob']:
                print(f"🚀 Prepared Molecules: {table_info['sdf_blob_count']:,}")
                if table_info['sdf_blob_count'] > 0:
                    prep_percentage = (table_info['sdf_blob_count'] / table_info['compound_count']) * 100
                    print(f"📈 Preparation Coverage: {prep_percentage:.1f}%")
            
            readiness = table_info['docking_readiness']
            print(f"\n🎯 DOCKING READINESS:")
            print(f"   Status: {readiness['docking_status']}")
            print(f"   Quality Score: {readiness['score']}/100")
            print(f"   Overall Status: {readiness['status']}")
            
            if readiness['issues']:
                print(f"   Issues:")
                for issue in readiness['issues']:
                    print(f"      • {issue}")
            else:
                print(f"   ✅ No issues detected")
            
            # Enhanced recommendations based on docking readiness
            print(f"\n💡 RECOMMENDATIONS:")
            if table_info['docking_status'] == 'ready for docking':
                print("   ✅ This table is ready for docking experiments")
                print("   🚀 You can proceed with receptor preparation and docking setup")
                print("   📊 All molecules have prepared 3D coordinates")
            else:
                print("   🔧 This table needs 3D molecular preparation before docking")
                
                if not table_info['has_sdf_blob']:
                    print("   📝 Steps needed:")
                    print("      1. Generate 3D conformers: chemspace.generate_mols_in_table()")
                    print("      2. This will create SDF blobs with 3D coordinates")
                    print("      3. Table will then be ready for docking")
                elif table_info['sdf_blob_count'] == 0:
                    print("   📝 SDF blob column exists but no molecules prepared")
                    print("      • Run: chemspace.generate_mols_in_table() to populate 3D coordinates")
                else:
                    partial_prep = (table_info['sdf_blob_count'] / table_info['compound_count']) * 100
                    print(f"   📊 Partial preparation: {partial_prep:.1f}% molecules ready")
                    print("      • Some molecules have 3D coordinates, others need preparation")
            
            print("=" * 60)
            
        except Exception as e:
            print(f"❌ Error showing table details: {e}")
    
    def save_pdb_model(self) -> Optional[Dict[str, Any]]:
        """
        The user is prompted to select a PDB file, which will be stored as a blob object in a table named pdb_files in the pdbs.db located in the 
        project_path/docking/receptors directory
        
        Returns:
            Optional[Dict[str, Any]]: Information about the saved PDB file or None if failed
        """
        import sqlite3
        import json
        from datetime import datetime
        
        try:
            # Prompt for PDB file path
            print(f"\n💾 SAVE PDB MODEL TO DATABASE")
            print("=" * 80)
            print(f"   This will store a PDB file as a blob in the database")
            print(f"   Location: {self.path}/docking/receptors/pdbs.db")
            print("=" * 80)
            
            # Get PDB file path
            while True:
                try:
                    pdb_file = input(f"\n📁 Enter path to PDB file: ").strip()
                    
                    if not pdb_file:
                        print(f"❌ No file path provided")
                        return None
                    
                    # Expand user path if needed
                    pdb_file = os.path.expanduser(pdb_file)
                    
                    # Check if file exists
                    if not os.path.exists(pdb_file):
                        print(f"❌ File not found: {pdb_file}")
                        retry = input(f"   Try again? (y/n, default: y): ").strip().lower()
                        if retry == 'n':
                            return None
                        continue
                    
                    # Check if it's a PDB file
                    if not pdb_file.lower().endswith('.pdb'):
                        print(f"⚠️  Warning: File doesn't have .pdb extension")
                        proceed = input(f"   Proceed anyway? (y/n, default: n): ").strip().lower()
                        if proceed != 'y':
                            continue
                    
                    break
                    
                except KeyboardInterrupt:
                    print(f"\n\n❌ Operation cancelled by user")
                    return None
            
            # Get PDB name/identifier
            default_name = os.path.splitext(os.path.basename(pdb_file))[0]
            while True:
                try:
                    pdb_name = input(f"\n🏷️  Enter name for this PDB (default: {default_name}): ").strip()
                    
                    if not pdb_name:
                        pdb_name = default_name
                    
                    # Check if name already exists
                    pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
                    
                    if os.path.exists(pdbs_db_path):
                        conn = sqlite3.connect(pdbs_db_path)
                        cursor = conn.cursor()
                        
                        # Check if table exists
                        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='pdb_files'")
                        if cursor.fetchone():
                            cursor.execute("SELECT COUNT(*) FROM pdb_files WHERE pdb_name = ?", (pdb_name,))
                            if cursor.fetchone()[0] > 0:
                                print(f"⚠️  PDB name '{pdb_name}' already exists in database")
                                overwrite = input(f"   Overwrite existing entry? (y/n, default: n): ").strip().lower()
                                conn.close()
                                if overwrite != 'y':
                                    continue
                        
                        conn.close()
                    
                    break
                    
                except KeyboardInterrupt:
                    print(f"\n\n❌ Operation cancelled by user")
                    return None
            
            # Get optional description
            try:
                description = input(f"\n📝 Enter description (optional, press Enter to skip): ").strip()
            except KeyboardInterrupt:
                print(f"\n   ⚠️  Description entry cancelled. Using empty description.")
                description = ""
            
            # Read PDB file content
            print(f"\n📖 Reading PDB file...")
            try:
                with open(pdb_file, 'rb') as f:
                    pdb_blob = f.read()
                
                file_size = len(pdb_blob)
                print(f"   ✓ File read successfully ({file_size:,} bytes)")
                
            except Exception as e:
                print(f"❌ Error reading PDB file: {e}")
                return None
            
            # Get basic file information
            file_info = {
                'original_path': pdb_file,
                'filename': os.path.basename(pdb_file),
                'file_size': file_size,
                'modified_date': datetime.fromtimestamp(os.path.getmtime(pdb_file)).isoformat()
            }
            
            # Create database and table
            os.makedirs(os.path.dirname(pdbs_db_path), exist_ok=True)
            
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            
            # Create pdb_files table if it doesn't exist
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS pdb_files (
                    file_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pdb_name TEXT UNIQUE NOT NULL,
                    project_name TEXT,
                    pdb_blob BLOB NOT NULL,
                    original_path TEXT,
                    filename TEXT,
                    file_size INTEGER,
                    description TEXT,
                    file_info TEXT,
                    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    last_modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    notes TEXT
                )
            ''')
            
            # Insert or replace PDB file
            print(f"\n💾 Saving to database...")
            
            cursor.execute('''
                INSERT OR REPLACE INTO pdb_files (
                    pdb_name, project_name, pdb_blob, original_path, filename,
                    file_size, description, file_info, notes
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                pdb_name,
                self.name,
                pdb_blob,
                pdb_file,
                os.path.basename(pdb_file),
                file_size,
                description,
                json.dumps(file_info, indent=2),
                f"PDB file saved from: {pdb_file}"
            ))
            
            file_id = cursor.lastrowid or cursor.execute(
                "SELECT file_id FROM pdb_files WHERE pdb_name = ?", 
                (pdb_name,)
            ).fetchone()[0]
            
            conn.commit()
            conn.close()
            
            # Display summary
            print(f"\n{'=' * 80}")
            print(f"✅ PDB MODEL SAVED SUCCESSFULLY")
            print(f"{'=' * 80}")
            print(f"   📋 File ID: {file_id}")
            print(f"   🏷️  Name: {pdb_name}")
            print(f"   📁 Original Path: {pdb_file}")
            print(f"   📊 File Size: {file_size:,} bytes")
            print(f"   📝 Description: {description or 'None'}")
            print(f"   💾 Database: {pdbs_db_path}")
            print(f"   📅 Saved: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"{'=' * 80}")
            
            print(f"\n💡 Use list_pdb_templates() to view all saved PDB files")
            
            return {
                'file_id': file_id,
                'pdb_name': pdb_name,
                'project_name': self.name,
                'original_path': pdb_file,
                'filename': os.path.basename(pdb_file),
                'file_size': file_size,
                'description': description,
                'database_path': pdbs_db_path,
                'saved_date': datetime.now().isoformat()
            }
            
        except sqlite3.Error as e:
            print(f"❌ Database error: {e}")
            return None
        except Exception as e:
            print(f"❌ Error saving PDB model: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def delete_pdb_model(self):
        """
        Will list the available pdb models in the database and prompt the user to select one to delete.
        """
        import sqlite3
        from datetime import datetime

        try:
            # Path to pdbs database
            pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
            
            # Check if database exists
            if not os.path.exists(pdbs_db_path):
                print(f"❌ No PDB database found")
                print(f"   Expected location: {pdbs_db_path}")
                return None
            
            # List available PDB models
            print(f"\n🗑️  DELETE PDB MODEL FROM DATABASE")
            print("=" * 80)
            
            models_list = self.list_pdb_models()
            
            if not models_list:
                print(f"❌ No PDB models available to delete")
                return None
            
            # Prompt user to select a model
            while True:
                try:
                    selection_input = input(f"\nEnter the model number to delete (or 'q' to quit): ").strip()
                    
                    if selection_input.lower() == 'q':
                        print(f"❌ Operation cancelled by user")
                        return None
                    
                    try:
                        model_number = int(selection_input)
                        if model_number < 1 or model_number > len(models_list):
                            print(f"❌ Invalid model number. Choose between 1 and {len(models_list)}")
                            continue
                        
                        selected_model = models_list[model_number - 1]
                        break
                        
                    except ValueError:
                        print(f"❌ Please enter a valid number or 'q' to quit")
                        continue
                        
                except KeyboardInterrupt:
                    print(f"\n\n❌ Operation cancelled by user")
                    return None
            
            # Confirm deletion
            print(f"\n⚠️  CONFIRM DELETION")
            print("=" * 80)
            print(f"   PDB Name: {selected_model['pdb_name']}")
            print(f"   Filename: {selected_model['filename']}")
            print(f"   File Size: {selected_model['file_size']:,} bytes" if selected_model['file_size'] else "   File Size: Unknown")
            print(f"   Created: {selected_model['created_date']}")
            print(f"   Description: {selected_model['description'] or 'None'}")
            print("=" * 80)
            
            confirm = input(f"\n⚠️  Are you sure you want to delete this PDB model? (yes/no, default: no): ").strip().lower()
            
            if confirm != 'yes':
                print(f"❌ Deletion cancelled by user")
                return None
            
            # Delete from database
            try:
                conn = sqlite3.connect(pdbs_db_path)
                cursor = conn.cursor()
                
                # Delete the record
                cursor.execute("DELETE FROM pdb_files WHERE file_id = ?", (selected_model['file_id'],))
                
                if cursor.rowcount == 0:
                    print(f"❌ Failed to delete PDB model from database")
                    conn.close()
                    return None
                
                conn.commit()
                conn.close()
                
                print(f"\n✓ PDB model '{selected_model['pdb_name']}' successfully deleted from database")
                print(f"   File ID: {selected_model['file_id']}")
                print(f"   Database: {pdbs_db_path}")
                
                return {
                    'status': 'deleted',
                    'pdb_name': selected_model['pdb_name'],
                    'file_id': selected_model['file_id'],
                    'deleted_date': datetime.now().isoformat()
                }
                
            except sqlite3.Error as e:
                print(f"❌ Database error while deleting PDB model: {e}")
                return None
                
        except Exception as e:
            print(f"❌ Error deleting PDB model: {e}")
            return None

    def write_pdb_model(self):
        """
        Will use the list_pdb_models to show the user available pdb files. The user is queried to select one model, after which the pdb file will be written to the project_path/docking/receptors directory using the original name
        
        Returns:
            Optional[Dict[str, Any]]: Information about the written PDB file or None if failed
        """
        import sqlite3
        from datetime import datetime
        
        try:
            # List available PDB models
            print(f"\n📁 RETRIEVE PDB MODEL FROM DATABASE")
            print("=" * 80)
            
            models_list = self.list_pdb_models()
            
            if not models_list:
                print(f"❌ No PDB models available to retrieve")
                return None
            
            # Prompt user to select a model
            while True:
                try:
                    selection_input = input(f"\nEnter the model number to retrieve (or 'q' to quit): ").strip()
                    
                    if selection_input.lower() == 'q':
                        print(f"❌ Operation cancelled by user")
                        return None
                    
                    try:
                        model_number = int(selection_input)
                        if model_number < 1 or model_number > len(models_list):
                            print(f"❌ Invalid model number. Choose between 1 and {len(models_list)}")
                            continue
                        
                        selected_model = models_list[model_number - 1]
                        break
                        
                    except ValueError:
                        print(f"❌ Invalid input. Please enter a number or 'q' to quit")
                        continue
                        
                except KeyboardInterrupt:
                    print(f"\n\n❌ Operation cancelled by user")
                    return None
            
            # Display selected model info
            print(f"\n{'=' * 80}")
            print(f"✓ Selected Model: {selected_model['pdb_name']}")
            print(f"  Original name: {selected_model['filename']}")
            print(f"  File size: {selected_model['file_size']:,} bytes" if selected_model['file_size'] else "  File size: Unknown")
            print(f"{'=' * 80}")
            
            # Retrieve the PDB blob from database
            print(f"\n📂 Retrieving PDB file from database...")
            
            pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT pdb_blob FROM pdb_files WHERE file_id = ?
            ''', (selected_model['file_id'],))
            
            result = cursor.fetchone()
            conn.close()
            
            if not result or not result[0]:
                print(f"❌ Failed to retrieve PDB blob from database")
                return None
            
            pdb_blob = result[0]
            
            # Determine output path and filename
            receptors_dir = os.path.join(self.path, 'docking', 'receptors')
            output_filename = selected_model['filename']  # Use original filename
            output_path = os.path.join(receptors_dir, output_filename)
            
            # Check if file already exists
            if os.path.exists(output_path):
                print(f"⚠️  File already exists: {output_path}")
                overwrite = input(f"   Overwrite? (y/n, default: n): ").strip().lower()
                if overwrite != 'y':
                    print(f"❌ Operation cancelled by user")
                    return None
            
            # Write PDB file to disk
            print(f"\n💾 Writing PDB file to disk...")
            try:
                os.makedirs(receptors_dir, exist_ok=True)
                
                with open(output_path, 'wb') as f:
                    f.write(pdb_blob)
                
                written_size = len(pdb_blob)
                print(f"   ✓ File written successfully ({written_size:,} bytes)")
                
            except Exception as e:
                print(f"❌ Error writing PDB file: {e}")
                return None
            
            # Display summary
            print(f"\n{'=' * 80}")
            print(f"✅ PDB FILE RETRIEVED AND WRITTEN SUCCESSFULLY")
            print(f"{'=' * 80}")
            print(f"   📋 Model ID: {selected_model['file_id']}")
            print(f"   🏷️  Model Name: {selected_model['pdb_name']}")
            print(f"   📁 Output Filename: {output_filename}")
            print(f"   📍 Output Path: {output_path}")
            print(f"   📊 File Size: {written_size:,} bytes")
            print(f"   📅 Retrieved: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"{'=' * 80}")
            
            print(f"\n💡 PDB file is ready for use in docking workflow")
            
            return {
                'file_id': selected_model['file_id'],
                'model_name': selected_model['pdb_name'],
                'filename': output_filename,
                'output_path': output_path,
                'file_size': written_size,
                'retrieved_date': datetime.now().isoformat()
            }
            
        except sqlite3.Error as e:
            print(f"❌ Database error: {e}")
            return None
        except Exception as e:
            print(f"❌ Error retrieving PDB file: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def list_pdb_models(self) -> Optional[List[Dict[str, Any]]]:
        """
        List the PDB models saved in the pdbs.db located in the 
        project_path/docking/receptors directory.
        
        Returns:
            Optional[List[Dict[str, Any]]]: List of PDB model information or None if failed
        """
        import sqlite3
        import json
        from datetime import datetime
        
        try:
            # Path to pdbs database
            pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
            
            # Check if database exists
            if not os.path.exists(pdbs_db_path):
                print(f"❌ No PDB database found")
                print(f"   Expected location: {pdbs_db_path}")
                print(f"\n💡 Save a PDB model first using save_pdb_model()")
                return None
            
            # Connect to database
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            
            # Check if pdb_files table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='pdb_files'")
            if not cursor.fetchone():
                print(f"❌ PDB files table not found in database")
                conn.close()
                return None
            
            # Retrieve all PDB models
            cursor.execute('''
                SELECT file_id, pdb_name, project_name, original_path, filename,
                    file_size, description, created_date, last_modified, notes
                FROM pdb_files
                ORDER BY created_date ASC
            ''')
            
            pdb_models = cursor.fetchall()
            conn.close()
            
            # Check if any models exist
            if not pdb_models:
                print(f"📋 NO PDB MODELS FOUND")
                print(f"   Database exists but contains no PDB models")
                print(f"\n💡 Save a PDB model using save_pdb_model()")
                return None
            
            # Display models
            print(f"\n📦 SAVED PDB MODELS")
            print("=" * 100)
            print(f"   Found {len(pdb_models)} PDB model(s) in database")
            print(f"   Database: {pdbs_db_path}")
            print("=" * 100)
            
            models_list = []
            
            for idx, (file_id, pdb_name, project_name, original_path, filename,
                    file_size, description, created_date, last_modified, notes) in enumerate(pdb_models, 1):
                
                # Store model info
                model_info = {
                    'file_id': file_id,
                    'pdb_name': pdb_name,
                    'project_name': project_name,
                    'original_path': original_path,
                    'filename': filename,
                    'file_size': file_size,
                    'description': description,
                    'created_date': created_date,
                    'last_modified': last_modified,
                    'notes': notes
                }
                models_list.append(model_info)
                
                # Display basic model details
                print(f"\n{idx}. 📦 {pdb_name}")
                print(f"   {'─' * 96}")
                print(f"   📋 ID: {file_id}")
                print(f"   📁 Filename: {filename}")
                print(f"   📍 Path: {original_path}")
                print(f"   📊 File Size: {file_size:,} bytes" if file_size else "   📊 File Size: Unknown")
                print(f"   📝 Description: {description or 'None'}")
                print(f"   📅 Created: {created_date}")
                
                # Display modification date if different from created date
                if last_modified and last_modified != created_date:
                    print(f"   🔄 Modified: {last_modified}")
                
                print(f"   {'─' * 96}")
            
            print("=" * 100)
            print(f"💡 Type 'details <number>' to view full information about a model")
            
            # Interactive details viewing
            while True:
                try:
                    user_input = input(f"\nEnter command (or press Enter to exit): ").strip()
                    
                    if not user_input:
                        break
                    
                    # Check if user wants details
                    if user_input.lower().startswith('details'):
                        parts = user_input.split()
                        if len(parts) < 2:
                            print(f"❌ Please specify a model number. Example: details 1")
                            continue
                        
                        try:
                            model_number = int(parts[1])
                            if model_number < 1 or model_number > len(models_list):
                                print(f"❌ Invalid model number. Choose between 1 and {len(models_list)}")
                                continue
                            
                            # Display detailed information for selected model
                            selected_model = models_list[model_number - 1]
                            print(f"\n{'═' * 100}")
                            print(f"🔍 DETAILED INFORMATION FOR: {selected_model['pdb_name']}")
                            print(f"{'═' * 100}")
                            print(f"   📋 File ID: {selected_model['file_id']}")
                            print(f"   🏷️  PDB Name: {selected_model['pdb_name']}")
                            print(f"   📁 Filename: {selected_model['filename']}")
                            print(f"   📍 Original Path: {selected_model['original_path']}")
                            print(f"   📦 Project: {selected_model['project_name']}")
                            
                            if selected_model['file_size']:
                                print(f"\n   📊 SIZE INFORMATION:")
                                print(f"      • File Size: {selected_model['file_size']:,} bytes")
                            
                            if selected_model['description']:
                                print(f"\n   📝 DESCRIPTION:")
                                print(f"      {selected_model['description']}")
                            
                            print(f"\n   📅 TIMESTAMPS:")
                            print(f"      • Created: {selected_model['created_date']}")
                            if selected_model['last_modified'] and selected_model['last_modified'] != selected_model['created_date']:
                                print(f"      • Last Modified: {selected_model['last_modified']}")
                            
                            if selected_model['notes']:
                                print(f"\n   📋 NOTES:")
                                print(f"      {selected_model['notes']}")
                            
                            print(f"{'═' * 100}")
                            
                        except ValueError:
                            print(f"❌ Invalid number format. Please use: details <number>")
                            continue
                    else:
                        print(f"❌ Unknown command. Use 'details <number>' to view full information")
                
                except KeyboardInterrupt:
                    print(f"\n\n👋 Exiting PDB models list")
                    break
            
            return models_list
            
        except sqlite3.Error as e:
            print(f"❌ Database error: {e}")
            return None
        except Exception as e:
            print(f"❌ Error listing PDB models: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def create_pdb_template(self, pdb_name: Optional[str] = None, verbose=True) -> Optional[Dict[str, Any]]:
        """
        Create a PDB file template from a PDB model stored in the database using save_pdb_model.
        The user is prompted to select a PDB model from the database, after which the PDB file 
        is analyzed and a register for the receptor is created in the pdb_templates table.
        
        Args:
            pdb_name (Optional[str]): Name for the receptor template. If None, user will be prompted.
            verbose (bool): Whether to display detailed output during processing.
            
        Returns:
            Optional[Dict[str, Any]]: PDB template registry information or None if failed
        """
        import sqlite3
        import tempfile
        from datetime import datetime
        
        try:
            print(f"🧬 CREATE PDB TEMPLATE FROM DATABASE")
            print("=" * 80)
            
            # Step 1: List and select PDB model from database
            print(f"\n📂 Retrieving PDB models from database...")
            models_list = self.list_pdb_models()
            
            if not models_list:
                print(f"❌ No PDB models available in database")
                return None
            
            # Prompt user to select a model
            while True:
                try:
                    selection_input = input(f"\nSelect PDB model number (or 'q' to quit): ").strip()
                    
                    if selection_input.lower() == 'q':
                        print(f"❌ PDB template creation cancelled")
                        return None
                    
                    try:
                        model_number = int(selection_input)
                        if model_number < 1 or model_number > len(models_list):
                            print(f"❌ Invalid model number. Choose between 1 and {len(models_list)}")
                            continue
                        
                        selected_model = models_list[model_number - 1]
                        break
                        
                    except ValueError:
                        print(f"❌ Invalid input. Please enter a number or 'q' to quit")
                        continue
                        
                except KeyboardInterrupt:
                    print(f"\n❌ PDB template creation cancelled")
                    return None
            
            # Step 2: Retrieve PDB blob from database
            print(f"\n📥 Retrieving PDB file from database...")
            pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT pdb_blob FROM pdb_files WHERE file_id = ?
            ''', (selected_model['file_id'],))
            
            result = cursor.fetchone()
            conn.close()
            
            if not result or not result[0]:
                print(f"❌ Failed to retrieve PDB blob from database")
                return None
            
            pdb_blob = result[0]
            
            # Step 3: Write blob to temporary file for analysis
            print(f"📝 Processing PDB file...")
            with tempfile.NamedTemporaryFile(mode='wb', suffix='.pdb', delete=False) as tmp_file:
                tmp_file.write(pdb_blob)
                temp_pdb_path = tmp_file.name
            
            try:
                # Step 4: Analyze PDB file
                print(f"🔬 Analyzing PDB structure...")
                pdb_analysis = self._analyze_pdb_file(temp_pdb_path)
                
                if not pdb_analysis:
                    print(f"❌ Failed to analyze PDB file")
                    os.unlink(temp_pdb_path)
                    return None
                
                if verbose:
                    self._display_pdb_analysis(pdb_analysis)
                
                # Step 5: Get PDB template name
                if pdb_name is None:
                    while True:
                        try:
                            default_name = f"{selected_model['pdb_name']}_template"
                            pdb_name = input(f"\n🏷️  Enter PDB template name (default: {default_name}): ").strip()
                            
                            if not pdb_name:
                                pdb_name = default_name
                            
                            if pdb_name.lower() in ['cancel', 'quit', 'exit']:
                                print("❌ PDB template creation cancelled")
                                os.unlink(temp_pdb_path)
                                return None
                            
                            # Validate template name
                            if not pdb_name.replace('_', '').replace('-', '').replace(' ', '').isalnum():
                                print("❌ Template name can only contain letters, numbers, spaces, hyphens, and underscores")
                                continue
                            
                            break
                            
                        except KeyboardInterrupt:
                            print("\n❌ PDB template creation cancelled")
                            os.unlink(temp_pdb_path)
                            return None
                
                # Step 6: Check if template name already exists
                existing_check = self._check_pdb_name_exists(pdb_name)
                if existing_check:
                    while True:
                        overwrite = input(f"\n⚠️  Template '{pdb_name}' already exists. Overwrite? (y/n): ").strip().lower()
                        if overwrite in ['y', 'yes']:
                            print("🔄 Will overwrite existing PDB template")
                            break
                        elif overwrite in ['n', 'no']:
                            print("❌ PDB template creation cancelled - name already exists")
                            os.unlink(temp_pdb_path)
                            return None
                        else:
                            print("❌ Please answer 'y' or 'n'")
                
                # Step 7: Create receptors directory structure
                receptors_base_dir = os.path.join(self.path, 'docking', 'receptors')
                os.makedirs(receptors_base_dir, exist_ok=True)
                
                # Create template-specific directory
                pdb_folder = self._create_receptor_folder(receptors_base_dir, pdb_name)
                
                if not pdb_folder:
                    print("❌ Failed to create PDB template folder")
                    os.unlink(temp_pdb_path)
                    return None
                
                # Step 8: Copy and process PDB file
                print(f"📂 Processing and saving PDB template...")
                processed_pdb_path, checked_pdb_path = self._process_pdb_file(temp_pdb_path, pdb_folder, pdb_analysis)
                
                if not processed_pdb_path:
                    print("❌ Failed to process PDB file")
                    os.unlink(temp_pdb_path)
                    return None
                
                # Step 9: Create template registry entry
                pdb_registry = self._create_pdb_registry_entry(
                    pdb_name, selected_model['original_path'], processed_pdb_path, pdb_folder, pdb_analysis, checked_pdb_path
                )
                
                if not pdb_registry:
                    print("❌ Failed to create PDB template registry")
                    os.unlink(temp_pdb_path)
                    return None
                
                if verbose:
                    # Step 10: Display success summary
                    self._display_pdb_creation_summary(pdb_registry)
                
                print(f"\n{'=' * 80}")
                print(f"✅ PDB TEMPLATE CREATED SUCCESSFULLY FROM DATABASE MODEL")
                print(f"{'=' * 80}")
                print(f"   Source Model: {selected_model['pdb_name']}")
                print(f"   Template Name: {pdb_name}")
                print(f"   Location: {pdb_folder}")
                print(f"{'=' * 80}\n")
                
                return pdb_registry
                
            finally:
                # Clean up temporary file
                if os.path.exists(temp_pdb_path):
                    os.unlink(temp_pdb_path)
                
        except Exception as e:
            print(f"❌ Error creating PDB template: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _analyze_pdb_file(self, pdb_file: str) -> Optional[Dict[str, Any]]:
        """
        Analyze PDB file structure and content, including detection and processing of alternate locations (altlocs)
        and gaps in the protein structure (missing residues).
        
        Args:
            pdb_file (str): Path to PDB file
            
        Returns:
            Optional[Dict[str, Any]]: PDB analysis results or None if failed
        """
        try:
            analysis = {
                'file_path': pdb_file,
                'file_size': os.path.getsize(pdb_file),
                'chains': set(),
                'residue_count': 0,
                'atom_count': 0,
                'hetero_count': 0,
                'has_waters': False,
                'has_ligands': False,
                'ligand_residues': set(),
                'resolution': None,
                'header_info': {},
                'altlocs': {},
                'has_altlocs': False,
                'gaps': {},  # New: gaps per chain
                'total_gaps': 0,  # New: total gaps
                'gap_residues': []  # New: residues involved in gaps
            }
            altlocs_per_residue = {}
            altloc_counts = {}
            residues_by_chain = {}

            with open(pdb_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    # Parse HEADER information
                    if line.startswith('HEADER'):
                        analysis['header_info']['title'] = line[10:50].strip()
                        analysis['header_info']['deposition_date'] = line[50:59].strip()
                        analysis['header_info']['pdb_id'] = line[62:66].strip()
                    
                    # Parse RESOLUTION
                    elif line.startswith('REMARK   2 RESOLUTION'):
                        try:
                            resolution_text = line[23:].strip()
                            if 'ANGSTROMS' in resolution_text:
                                resolution_value = float(resolution_text.split()[0])
                                analysis['resolution'] = resolution_value
                        except (ValueError, IndexError):
                            pass
                    
                    # Parse ATOM records
                    elif line.startswith('ATOM'):
                        analysis['atom_count'] += 1
                        chain_id = line[21].strip()
                        if chain_id:
                            analysis['chains'].add(chain_id)
                        # Track residues by chain for gap analysis
                        resnum = int(line[22:26].strip())
                        resname = line[17:20].strip()
                        residues_by_chain.setdefault(chain_id, set()).add(resnum)
                        analysis['residue_count'] += 1

                        # --- ALTLOC detection ---
                        altloc = line[16] if len(line) > 16 else " "
                        if altloc and altloc not in (" ", ""):
                            analysis['has_altlocs'] = True
                            res_uid = (chain_id, resname, str(resnum))
                            altlocs_per_residue.setdefault(res_uid, set()).add(altloc)
                            altloc_counts[altloc] = altloc_counts.get(altloc, 0) + 1

                    # Parse HETATM records
                    elif line.startswith('HETATM'):
                        analysis['hetero_count'] += 1
                        residue_name = line[17:20].strip()
                        if residue_name == 'HOH':
                            analysis['has_waters'] = True
                        else:
                            analysis['has_ligands'] = True
                            analysis['ligand_residues'].add(residue_name)
                        chain_id = line[21].strip()
                        if chain_id:
                            analysis['chains'].add(chain_id)
                        # --- ALTLOC detection for HETATM ---
                        altloc = line[16] if len(line) > 16 else " "
                        resnum = int(line[22:26].strip())
                        resname = line[17:20].strip()
                        if altloc and altloc not in (" ", ""):
                            analysis['has_altlocs'] = True
                            res_uid = (chain_id, resname, str(resnum))
                            altlocs_per_residue.setdefault(res_uid, set()).add(altloc)
                            altloc_counts[altloc] = altloc_counts.get(altloc, 0) + 1

            # Convert sets to lists for JSON serialization
            analysis['chains'] = sorted(list(analysis['chains']))
            analysis['ligand_residues'] = sorted(list(analysis['ligand_residues']))

            # Add altlocs summary
            if analysis['has_altlocs']:
                analysis['altlocs'] = {
                    'total_altlocs': sum(len(v) for v in altlocs_per_residue.values()),
                    'unique_altloc_ids': sorted(list(altloc_counts.keys())),
                    'altloc_counts': dict(sorted(altloc_counts.items())),
                    'residues_with_altlocs': len(altlocs_per_residue),
                }
            else:
                analysis['altlocs'] = {
                    'total_altlocs': 0,
                    'unique_altloc_ids': [],
                    'altloc_counts': {},
                    'residues_with_altlocs': 0,
                }

            # --- GAP ANALYSIS ---
            gap_count = 0
            gap_residues = []
            gaps_per_chain = {}
            for chain_id, resnums in residues_by_chain.items():
                sorted_resnums = sorted(resnums)
                chain_gaps = []
                for i in range(1, len(sorted_resnums)):
                    prev_res = sorted_resnums[i - 1]
                    curr_res = sorted_resnums[i]
                    if curr_res - prev_res > 1:
                        gap_count += 1
                        gap_range = (prev_res + 1, curr_res - 1)
                        chain_gaps.append({'start': prev_res, 'end': curr_res, 'gap_residues': list(range(gap_range[0], gap_range[1] + 1))})
                        gap_residues.extend([(chain_id, r) for r in range(gap_range[0], gap_range[1] + 1)])
                gaps_per_chain[chain_id] = chain_gaps
            analysis['gaps'] = gaps_per_chain
            analysis['total_gaps'] = gap_count
            analysis['gap_residues'] = gap_residues

            return analysis

        except Exception as e:
            print(f"❌ Error analyzing PDB file: {e}")
            return None

    def _process_pdb_file(self, pdb_file: str, pdb_folder: str, analysis: Dict[str, Any]) -> str:
        """
        Process and copy PDB file to receptor folder with chain, ligand, and water selection options.
        Now also manages alternate locations (altlocs) and water molecules if present.
        """
        try:
            import shutil

            # Copy original PDB to original folder
            original_dir = os.path.join(pdb_folder, 'original')
            original_pdb_path = os.path.join(original_dir, 'receptor_original.pdb')
            os.makedirs(original_dir, exist_ok=True)
            shutil.copy2(pdb_file, original_pdb_path)
            print(f"   📋 Copied original PDB file to: {os.path.relpath(original_pdb_path)}")

            # Create processed PDB path in processed folder
            processed_dir = os.path.join(pdb_folder, 'processed')
            os.makedirs(processed_dir, exist_ok=True)
            processed_pdb_path = os.path.join(processed_dir, 'receptor.pdb')
            shutil.copy2(original_pdb_path, processed_pdb_path)
            print(f"   🔄 Created working copy for processing: {os.path.relpath(processed_pdb_path)}")

            # --- ALTLOC MANAGEMENT ---
            altloc_handling = None
            altloc_map = {}
            if analysis.get('has_altlocs', False):
                altlocs = analysis.get('altlocs', {})
                print(f"\n⚠️  Alternate locations (altlocs) detected in structure.")
                print(f"   • Residues with altlocs: {altlocs.get('residues_with_altlocs', 0)}")
                print(f"   • Unique altloc IDs: {', '.join(altlocs.get('unique_altloc_ids', [])) or 'None'}")
                print(f"   • Total altloc atoms: {altlocs.get('total_altlocs', 0)}")
                print("   ⚠️  You must select which altlocs to keep in the processed file.")

                # Query user for altloc handling strategy
                altloc_handling, altloc_map = self._query_altloc_handling(original_pdb_path, analysis)
                if altloc_handling == "cancel":
                    print("   ❌ Altloc selection cancelled. Aborting processing.")
                    return ""

                # Apply altloc filtering to processed file
                self._filter_altlocs_in_pdb(processed_pdb_path, altloc_handling, altloc_map)
                print(f"   ✅ Applied altloc filtering: {altloc_handling}")

            # Process chains if multiple chains exist
            selected_chains = None
            if len(analysis['chains']) > 1:
                print(f"   🔗 Multi-chain structure detected: {', '.join(analysis['chains'])}")
                while True:
                    try:
                        choice = input(f"   🧬 Keep all chains or select specific chain? (all/select): ").strip().lower()
                        if choice in ['all', 'a']:
                            print("   ✅ Keeping all chains")
                            selected_chains = analysis['chains']
                            break
                        elif choice in ['select', 's']:
                            selected_chain = self._select_chain_interactive(analysis['chains'])
                            if selected_chain:
                                selected_chains = [selected_chain]
                                print(f"   ✅ Selected chain: {selected_chain}")
                                break
                            else:
                                print("   ⚠️  No chain selected, keeping all chains")
                                selected_chains = analysis['chains']
                                break
                        else:
                            print("   ❌ Please enter 'all' or 'select'")
                            continue
                    except KeyboardInterrupt:
                        print("\n   ⚠️  Using all chains by default")
                        selected_chains = analysis['chains']
                        break
            else:
                selected_chains = analysis['chains']

            # Process ligands if any exist
            selected_ligands = None
            if analysis['has_ligands'] and analysis['ligand_residues']:
                print(f"\n   💊 Ligand residues detected: {', '.join(analysis['ligand_residues'])}")
                ligand_details = self._analyze_ligands_in_pdb(original_pdb_path, analysis['ligand_residues'], selected_chains)
                if ligand_details:
                    self._display_ligand_analysis(ligand_details)
                    while True:
                        try:
                            print("\n   🧪 LIGAND OPTIONS:")
                            print("   1. Remove all ligands (clean receptor)")
                            print("   2. Keep all ligands")
                            print("   3. Select specific ligands to keep")
                            print("   4. Keep only co-crystallized ligands (exclude waters)")
                            ligand_choice = input("   💊 Select ligand handling (1-4): ").strip()
                            if ligand_choice == '1':
                                print("   🧹 Will remove all ligands from processed PDB")
                                selected_ligands = []
                                break
                            elif ligand_choice == '2':
                                print("   ✅ Will keep all ligands")
                                selected_ligands = list(analysis['ligand_residues'])
                                break
                            elif ligand_choice == '3':
                                selected_ligands = self._select_ligands_interactive(ligand_details)
                                if selected_ligands is not None:
                                    if selected_ligands:
                                        print(f"   ✅ Selected ligands: {', '.join(selected_ligands)}")
                                    else:
                                        print("   🧹 No ligands selected - will create clean receptor")
                                    break
                                else:
                                    print("   ⚠️  Ligand selection cancelled, keeping all ligands")
                                    selected_ligands = list(analysis['ligand_residues'])
                                    break
                            elif ligand_choice == '4':
                                non_water_ligands = [lig for lig in analysis['ligand_residues'] if lig != 'HOH']
                                selected_ligands = non_water_ligands
                                if non_water_ligands:
                                    print(f"   ✅ Will keep co-crystallized ligands: {', '.join(non_water_ligands)}")
                                else:
                                    print("   🧹 No co-crystallized ligands found - will create clean receptor")
                                break
                            else:
                                print("   ❌ Invalid choice. Please enter 1, 2, 3, or 4")
                                continue
                        except KeyboardInterrupt:
                            print("\n   ⚠️  Keeping all ligands by default")
                            selected_ligands = list(analysis['ligand_residues'])
                            break
                else:
                    selected_ligands = list(analysis['ligand_residues'])
            else:
                selected_ligands = []

            # --- WATER MANAGEMENT ---
            selected_waters = []
            water_details = self._analyze_waters_in_pdb(original_pdb_path, selected_chains)
            if water_details:
                print(f"\n   💧 Water molecules detected: {len(water_details)}")
                selected_waters = self._select_waters_interactive(water_details)
                if selected_waters is None:
                    print("   ⚠️  Water selection cancelled, will remove all waters.")
                    selected_waters = []
                elif selected_waters:
                    print(f"   ✅ Will keep {len(selected_waters)} water(s): {', '.join(selected_waters)}")
                else:
                    print("   🧹 Will remove all water molecules.")

            # Create processed PDB file with selected chains, ligands, and waters
            if (selected_chains != analysis['chains'] or
                selected_ligands != list(analysis['ligand_residues']) or
                (water_details and selected_waters != list(water_details.keys()))):
                print(f"   🔄 Applying filters to processed PDB (chains, ligands, waters)...")
                
                filtered_pdb_path = self._create_filtered_pdb(
                    processed_pdb_path, selected_chains, selected_ligands, selected_waters, analysis
                )
                
                if not filtered_pdb_path:
                    print("   ❌ Failed to create filtered PDB file")
                    return ""
                processed_pdb_path = filtered_pdb_path
            else:
                print(f"   ✅ No filtering needed - processed PDB ready")

            # Create reference ligand file if ligands were kept
            if selected_ligands:
                self._extract_reference_ligands(processed_pdb_path, pdb_folder, selected_ligands)

            # Create processing summary file
            self._create_processing_summary(pdb_folder, original_pdb_path, processed_pdb_path, selected_chains, selected_ligands, analysis)
            
            # --- Add missing atoms using tleap ---
            tleap_processed_file = self._add_missing_atoms_in_pdb_tleap(processed_pdb_path, analysis)
            
            # Process tleap_processed_file to add the element name column if missing
            receptor_moldf_dict = self._refine_receptor_model(tleap_processed_file, processed_pdb_path)

            # Write (overwrite) the tleap_processed_file
            from moldf import write_pdb
            
            write_pdb(receptor_moldf_dict, tleap_processed_file)
            
            print(f"✅ tleap processing complete. Output: {tleap_processed_file}")
            
            return processed_pdb_path, tleap_processed_file
            

        except Exception as e:
            print(f"❌ Error processing PDB file: {e}")
            return ""

    def _analyze_waters_in_pdb(self, pdb_file: str, selected_chains: Optional[List[str]] = None) -> Optional[Dict[str, Any]]:
        """
        Analyze water molecules (HOH) in the PDB file.
        Only include waters in chains defined in selected_chains if provided.
        Returns a dict: {water_id: {'chain':..., 'resnum':..., 'center':(x,y,z)}}
        """
        try:
            waters = {}
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('HETATM') and line[17:20].strip() == 'HOH':
                        chain_id = line[21].strip()
                        # Only include if chain is in selected_chains (if provided)
                        if selected_chains and chain_id not in selected_chains:
                            continue
                        resnum = line[22:26].strip()
                        atom_name = line[12:16].strip()
                        water_id = f"HOH_{chain_id}_{resnum}"
                        if water_id not in waters:
                            waters[water_id] = {
                                'chain': chain_id,
                                'resnum': resnum,
                                'atoms': [],
                                'coordinates': []
                            }
                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            waters[water_id]['coordinates'].append((x, y, z))
                        except ValueError:
                            pass
                        waters[water_id]['atoms'].append(atom_name)
            
            print(selected_chains)

            # Calculate center for each water
            for wid, info in waters.items():
                if info['coordinates']:
                    coords = info['coordinates']
                    center = tuple(sum(c[i] for c in coords) / len(coords) for i in range(3))
                    info['center'] = center
                else:
                    info['center'] = (0, 0, 0)
            return waters if waters else None
        except Exception as e:
            print(f"   ⚠️  Error analyzing waters: {e}")
            return None

    def _select_waters_interactive(self, water_details: Dict[str, Any]) -> Optional[List[str]]:
        """
        Interactive selection of water molecules to keep.
        Returns a list of water_ids to keep, or None if cancelled.
        Allows selection only by ResNum.
        """
        try:
            print("\n   💧 WATER SELECTION:")
            print("   " + "-" * 60)
            print(f"   {'#':<3} {'Water ID':<12} {'Chain':<6} {'ResNum':<7} {'Atoms':<7} {'Center (Å)':<20}")
            print("   " + "-" * 60)
            water_list = list(water_details.items())
            for i, (wid, info) in enumerate(water_list, 1):
                center = info['center']
                print(f"   {i:<3} {wid:<12} {info['chain']:<6} {info['resnum']:<7} {len(info['atoms']):<7} "
                      f"({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
            print("   " + "-" * 60)
            print("   Commands:")
            print("   • 'all' to keep all waters")
            print("   • 'none' to remove all waters")
            print("   • Comma-separated ResNums (e.g., 101,202) to select by residue number (all chains)")
            print("   • 'cancel' to abort water selection")
            while True:
                selection = input("   💧 Select waters to keep: ").strip().lower()
                if selection in ['all', 'a']:
                    return [wid for wid, _ in water_list]
                elif selection in ['none', 'n']:
                    return []
                elif selection in ['cancel', 'quit', 'exit']:
                    return None
                else:
                    try:
                        selected = set()
                        items = [x.strip() for x in selection.split(',') if x.strip()]
                        for item in items:
                            # Only allow selection by ResNum (must be alphanumeric)
                            for wid, info in water_list:
                                if info['resnum'] == item:
                                    selected.add(wid)
                        if selected:
                            return list(selected)
                        else:
                            print("   ❌ No valid selections made (use ResNum only)")
                            continue
                    except Exception:
                        print("   ❌ Invalid input. Use ResNums, 'all', or 'none'")
                        continue
        except KeyboardInterrupt:
            return None
        except Exception as e:
            print(f"   ⚠️  Error in water selection: {e}")
            return None

    def _create_filtered_pdb(self, pdb_file: str, selected_chains: List[str],
                                        selected_ligands: List[str], selected_waters: List[str],
                                        analysis: Dict[str, Any]) -> str:
        """
        Create a filtered PDB file with only selected chains, ligands, and waters.
        Adds a 'TER' card after each ligand/water for clarity, after the last ligand in selected_ligands,
        and inserts 'TER' cards between gaps as registered in the analysis dictionary.
        """
        try:
            temp_filtered_path = f"{pdb_file}.tmp"
            atoms_kept = 0
            ligand_atoms_kept = 0
            water_atoms_kept = 0
            last_ligand_res = None
            last_water_res = None

            # Prepare gap lookup: {chain_id: set of residue numbers that are gap starts}
            gap_starts = {}
            for chain_id, gaps in analysis.get('gaps', {}).items():
                gap_starts[chain_id] = set(gap['start'] for gap in gaps)

            last_chain = None
            last_resnum = None

            with open(pdb_file, 'r') as infile, open(temp_filtered_path, 'w') as outfile:
                for line in infile:
                    write_line = False
                    if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'REMARK')):
                        write_line = True
                    elif line.startswith('ATOM'):
                        chain_id = line[21].strip()
                        resnum = int(line[22:26].strip())
                        if not selected_chains or chain_id in selected_chains:
                            # Insert TER if a gap is detected before this residue (only once per gap)
                            if last_chain == chain_id and last_resnum is not None:
                                # Only insert TER when transitioning from the gap-start residue to the next residue
                                if chain_id in gap_starts and last_resnum in gap_starts[chain_id]:
                                    if resnum != last_resnum:
                                        outfile.write("TER\n")
                            write_line = True
                            atoms_kept += 1
                            last_chain = chain_id
                            last_resnum = resnum
                    elif line.startswith('TER'):
                        chain_id = line[21].strip()
                        if not selected_chains or chain_id in selected_chains:
                            write_line = True
                            atoms_kept += 1
                    elif line.startswith('HETATM'):
                        residue_name = line[17:20].strip()
                        chain_id = line[21].strip()
                        residue_number = line[22:26].strip()
                        res_uid = (chain_id, residue_name, residue_number)
                        water_id = f"HOH_{chain_id}_{residue_number}"
                        # Waters
                        if residue_name == 'HOH':
                            if selected_waters and water_id in selected_waters:
                                if last_water_res is None:
                                    outfile.write("TER\n")
                                write_line = True
                                water_atoms_kept += 1
                                if last_water_res and res_uid != last_water_res:
                                    outfile.write("TER\n")
                                last_water_res = res_uid
                        # Ligands
                        elif selected_ligands and residue_name in selected_ligands:
                            if not selected_chains or chain_id in selected_chains:
                                write_line = True
                                ligand_atoms_kept += 1
                                if last_ligand_res and res_uid != last_ligand_res:
                                    outfile.write("TER\n")
                                last_ligand_res = res_uid
                    if write_line:
                        outfile.write(line)
                # After all lines, if any water was written, add a TER card
                if last_water_res:
                    outfile.write("TER\n")
                outfile.write("END\n")
                outfile.truncate()
            import shutil
            shutil.move(temp_filtered_path, pdb_file)
            print(f"   ✅ Applied filtering to PDB file (chains, ligands, waters, gaps).")
            print(f"      • Chains kept: {', '.join(selected_chains) if selected_chains else 'all'}")
            print(f"      • Protein atoms: {atoms_kept:,}")
            if selected_ligands and ligand_atoms_kept > 0:
                ligand_names = [lig for lig in selected_ligands if lig != 'HOH']
                if ligand_names:
                    print(f"      • Ligand atoms: {ligand_atoms_kept:,} ({', '.join(ligand_names)})")
            if selected_waters and water_atoms_kept > 0:
                print(f"      • Water atoms: {water_atoms_kept:,}")
            elif selected_waters == []:
                print(f"      • All waters removed.")
            return pdb_file
        except Exception as e:
            print(f"   ❌ Error creating filtered PDB with waters: {e}")
            temp_filtered_path = f"{pdb_file}.tmp"
            if os.path.exists(temp_filtered_path):
                os.remove(temp_filtered_path)
            return ""

    def _query_altloc_handling(self, pdb_file: str, analysis: Dict[str, Any]) -> Tuple[str, dict]:
        """
        Query the user for how to handle altlocs: keep first, keep highest occupancy, or select manually.
        Returns the chosen strategy and a mapping if manual selection.
        """
        altlocs = analysis.get('altlocs', {})
        unique_altlocs = altlocs.get('unique_altloc_ids', [])
        print("\n   ALTLOC HANDLING OPTIONS:")
        print("   1. Keep only the first altloc (A/B/C...) for each residue")
        print("   2. Keep altloc with highest occupancy (if available)")
        print("   3. Select altloc manually for each residue with altlocs")
        print("   4. Cancel processing")
        while True:
            choice = input("   Select altloc handling (1-4): ").strip()
            if choice == '1':
                return "first", {}
            elif choice == '2':
                return "highest_occupancy", {}
            elif choice == '3':
                # Manual selection per residue
                altloc_map = self._manual_select_altlocs(pdb_file, analysis)
                return "manual", altloc_map
            elif choice == '4' or choice.lower() in ['cancel', 'quit', 'exit']:
                return "cancel", {}
            else:
                print("   ❌ Invalid choice. Please enter 1, 2, 3, or 4.")

    def _manual_select_altlocs(self, pdb_file: str, analysis: Dict[str, Any]) -> dict:
        """
        For each residue with altlocs, ask the user which altloc to keep.
        Returns a mapping: {res_uid: altloc_id}
        """
        altlocs = {}
        altlocs_per_residue = {}
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    chain_id = line[21].strip()
                    resname = line[17:20].strip()
                    resnum = line[22:26].strip()
                    altloc = line[16] if len(line) > 16 else " "
                    if altloc and altloc not in (" ", ""):
                        res_uid = (chain_id, resname, resnum)
                        altlocs_per_residue.setdefault(res_uid, set()).add(altloc)
        for res_uid, altloc_set in altlocs_per_residue.items():
            print(f"   Residue {res_uid} has altlocs: {', '.join(sorted(altloc_set))}")
            while True:
                sel = input(f"   Select altloc to keep for {res_uid} (or leave blank for first): ").strip().upper()
                if not sel:
                    altlocs[res_uid] = sorted(altloc_set)[0]
                    break
                elif sel in altloc_set:
                    altlocs[res_uid] = sel
                    break
                else:
                    print(f"   ❌ Invalid altloc. Choose from: {', '.join(sorted(altloc_set))}")
        return altlocs

    def _filter_altlocs_in_pdb(self, pdb_file: str, handling: str, altloc_map: dict) -> None:
        """
        Filter altlocs in the given pdb_file in-place according to the chosen strategy.
        Before writing the outfile, the altloc field is cleared (set to a space).
        """
        import tempfile
        tmp_path = pdb_file + ".altloc.tmp"
        with open(pdb_file, 'r') as infile, open(tmp_path, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    altloc = line[16] if len(line) > 16 else " "
                    chain_id = line[21].strip()
                    resname = line[17:20].strip()
                    resnum = line[22:26].strip()
                    res_uid = (chain_id, resname, resnum)
                    write_this = False
                    if altloc and altloc not in (" ", ""):
                        if handling == "first":
                            if not hasattr(self, "_altloc_first_seen"):
                                self._altloc_first_seen = {}
                            if res_uid not in self._altloc_first_seen:
                                self._altloc_first_seen[res_uid] = altloc
                            if altloc == self._altloc_first_seen[res_uid]:
                                write_this = True
                        elif handling == "highest_occupancy":
                            if not hasattr(self, "_altloc_occ"):
                                self._altloc_occ = {}
                            occ = float(line[54:60].strip() or 0)
                            if res_uid not in self._altloc_occ or occ > self._altloc_occ[res_uid][0]:
                                self._altloc_occ[res_uid] = (occ, line)
                        elif handling == "manual":
                            if res_uid in altloc_map and altloc == altloc_map[res_uid]:
                                write_this = True
                        else:
                            write_this = True
                    else:
                        # No altloc, always keep
                        write_this = True

                    if write_this:
                        # Clear altloc field (column 17, index 16)
                        new_line = line[:16] + " " + line[17:]
                        outfile.write(new_line)
                else:
                    outfile.write(line)
            # For highest_occupancy, after reading all lines, write the best for each residue
            if handling == "highest_occupancy" and hasattr(self, "_altloc_occ"):
                for occ, line in self._altloc_occ.values():
                    # Clear altloc field before writing
                    new_line = line[:16] + " " + line[17:]
                    outfile.write(new_line)
                del self._altloc_occ
            if hasattr(self, "_altloc_first_seen"):
                del self._altloc_first_seen
        # Replace original file
        import shutil
        shutil.move(tmp_path, pdb_file)

    def _create_processing_summary(self, receptor_folder: str, original_pdb_path: str, 
                                processed_pdb_path: str, selected_chains: List[str], 
                                selected_ligands: List[str], analysis: Dict[str, Any]) -> None:
        """
        Create a summary file documenting the processing steps applied, including altloc handling.
        
        Args:
            receptor_folder (str): Receptor folder path
            original_pdb_path (str): Path to original PDB file
            processed_pdb_path (str): Path to processed PDB file
            selected_chains (List[str]): Chains that were kept
            selected_ligands (List[str]): Ligands that were kept
            analysis (Dict): PDB analysis results
        """
        try:
            from datetime import datetime
            
            summary_path = os.path.join(receptor_folder, 'processing_summary.txt')
            
            with open(summary_path, 'w') as f:
                f.write(f"PDB PROCESSING SUMMARY\n")
                f.write(f"=" * 50 + "\n")
                f.write(f"Processing Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Project: {self.name}\n\n")
                
                f.write(f"ORIGINAL FILE:\n")
                f.write(f"  Path: {original_pdb_path}\n")
                f.write(f"  Size: {analysis['file_size']:,} bytes\n")
                f.write(f"  Chains: {', '.join(analysis['chains'])}\n")
                f.write(f"  Atoms: {analysis['atom_count']:,}\n")
                f.write(f"  Ligands: {', '.join(analysis['ligand_residues']) if analysis['ligand_residues'] else 'None'}\n")
                if analysis.get('resolution'):
                    f.write(f"  Resolution: {analysis['resolution']:.2f} Å\n")
                f.write("\n")
                
                f.write(f"PROCESSING APPLIED:\n")
                f.write(f"  Chains kept: {', '.join(selected_chains) if selected_chains else 'All'}\n")
                f.write(f"  Ligands kept: {', '.join(selected_ligands) if selected_ligands else 'None (clean receptor)'}\n")
                # Altloc summary
                if analysis.get('has_altlocs', False):
                    altlocs = analysis.get('altlocs', {})
                    f.write(f"  Altlocs present: Yes\n")
                    f.write(f"    • Residues with altlocs: {altlocs.get('residues_with_altlocs', 0)}\n")
                    f.write(f"    • Unique altloc IDs: {', '.join(altlocs.get('unique_altloc_ids', [])) or 'None'}\n")
                    f.write(f"    • Total altloc atoms: {altlocs.get('total_altlocs', 0)}\n")
                    # If altloc handling info is available, add it
                    altloc_handling_path = os.path.join(receptor_folder, 'processed', 'altloc_handling.txt')
                    if os.path.exists(altloc_handling_path):
                        f.write(f"    • Altloc handling: See processed/altloc_handling.txt\n")
                    else:
                        f.write(f"    • Altloc handling: Only one altloc per residue kept (see code for details)\n")
                else:
                    f.write(f"  Altlocs present: No\n")
                f.write("\n")
                
                f.write(f"PROCESSED FILE:\n")
                f.write(f"  Path: {processed_pdb_path}\n")
                
                # Get processed file stats
                if os.path.exists(processed_pdb_path):
                    processed_size = os.path.getsize(processed_pdb_path)
                    f.write(f"  Size: {processed_size:,} bytes\n")
                    
                    # Count atoms in processed file
                    with open(processed_pdb_path, 'r') as pdb:
                        processed_atoms = sum(1 for line in pdb if line.startswith('ATOM'))
                    f.write(f"  Protein atoms: {processed_atoms:,}\n")
                
                f.write("\n")
                f.write(f"FOLDER STRUCTURE:\n")
                f.write(f"  original/     - Original PDB file\n")
                f.write(f"  processed/    - Processed receptor file\n")
                f.write(f"  analysis/     - Reference ligands and analysis\n")
                f.write(f"  grid_files/   - Grid files (to be created)\n")
                f.write(f"  logs/         - Processing and docking logs\n")
            
            print(f"   📄 Created processing summary: {os.path.relpath(summary_path)}")
            
        except Exception as e:
            print(f"   ⚠️  Could not create processing summary: {e}")

    def _extract_reference_ligands(self, pdb_file: str, receptor_folder: str, selected_ligands: List[str]) -> None:
        """
        Extract reference ligands to separate files for docking reference.
        Adds a 'TER' card after the protein and after each ligand for clarity.
        
        Args:
            pdb_file (str): Processed PDB file path
            receptor_folder (str): Receptor folder path
            selected_ligands (List[str]): List of ligand residue names that were kept
        """
        try:
            if not selected_ligands or selected_ligands == ['HOH']:
                return  # No ligands to extract or only waters
            
            reference_dir = os.path.join(receptor_folder, 'analysis')
            os.makedirs(reference_dir, exist_ok=True)
            
            # Extract each ligand type to separate files
            ligands_extracted = []
            
            for ligand_name in selected_ligands:
                if ligand_name == 'HOH':
                    continue  # Skip water molecules
                
                ligand_file = os.path.join(reference_dir, f"reference_ligand_{ligand_name}.pdb")
                ligand_count = 0
                last_residue = None
                protein_section_done = False

                with open(pdb_file, 'r') as infile, open(ligand_file, 'w') as ligand_out:
                    # Write header
                    ligand_out.write(f"HEADER    REFERENCE LIGAND {ligand_name}\n")
                    ligand_out.write(f"REMARK    Extracted from processed receptor\n")
                    ligand_out.write(f"REMARK    Source: {os.path.basename(pdb_file)}\n")
                    ligand_out.write(f"REMARK    Project: {self.name}\n")
                    
                    for line in infile:
                        # Add 'TER' after the last ATOM line and before the first HETATM for this ligand
                        if not protein_section_done and line.startswith('HETATM'):
                            ligand_out.write("TER\n")
                            protein_section_done = True
                        if line.startswith('HETATM'):
                            residue_name = line[17:20].strip()
                            chain_id = line[21].strip()
                            residue_number = line[22:26].strip()
                            res_uid = (chain_id, residue_name, residue_number)
                            if residue_name == ligand_name:
                                ligand_out.write(line)
                                ligand_count += 1
                                # If residue changes, write TER
                                if last_residue and res_uid != last_residue:
                                    ligand_out.write("TER\n")
                                last_residue = res_uid
                    if ligand_count > 0:
                        ligand_out.write("TER\n")
                        ligand_out.write("END\n")
                
                if ligand_count > 0:
                    print(f"   📄 Extracted reference ligand: {ligand_name} ({ligand_count} atoms)")
                    ligands_extracted.append(ligand_name)
                else:
                    # Remove empty file
                    os.remove(ligand_file)
            
            if ligands_extracted:
                print(f"   📁 Reference ligands saved to: {os.path.relpath(reference_dir)}")
            
        except Exception as e:
            print(f"   ⚠️  Error extracting reference ligands: {e}")

    def _analyze_ligands_in_pdb(self, pdb_file: str, ligand_residues: List[str], selected_chains: Optional[List[str]] = None) -> Optional[Dict[str, Any]]:
        """
        Analyze ligands in PDB file to provide detailed information for selection.
        Only retain ligands in chains defined in selected_chains if provided.

        Args:
            pdb_file (str): Path to PDB file
            ligand_residues (List[str]): List of ligand residue names
            selected_chains (Optional[List[str]]): List of chain IDs to include

        Returns:
            Optional[Dict[str, Any]]: Detailed ligand analysis or None if failed
        """
        try:
            ligand_details = {}

            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('HETATM'):
                        residue_name = line[17:20].strip()
                        chain_id = line[21].strip()
                        residue_number = line[22:26].strip()
                        atom_name = line[12:16].strip()

                        # Only include ligands in selected_chains if provided
                        if selected_chains and chain_id not in selected_chains:
                            continue

                        if residue_name in ligand_residues:
                            # Create unique ligand identifier
                            ligand_id = f"{residue_name}_{chain_id}_{residue_number}"

                            if ligand_id not in ligand_details:
                                ligand_details[ligand_id] = {
                                    'residue_name': residue_name,
                                    'chain_id': chain_id,
                                    'residue_number': residue_number,
                                    'atoms': [],
                                    'atom_count': 0,
                                    'is_water': residue_name == 'HOH',
                                    'coordinates': []
                                }

                            # Extract coordinates
                            try:
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])
                                ligand_details[ligand_id]['coordinates'].append((x, y, z))
                            except ValueError:
                                pass

                            ligand_details[ligand_id]['atoms'].append(atom_name)
                            ligand_details[ligand_id]['atom_count'] += 1

            # Calculate ligand properties
            for ligand_id, details in ligand_details.items():
                if details['coordinates']:
                    # Calculate center of mass (approximate)
                    coords = details['coordinates']
                    center_x = sum(coord[0] for coord in coords) / len(coords)
                    center_y = sum(coord[1] for coord in coords) / len(coords)
                    center_z = sum(coord[2] for coord in coords) / len(coords)
                    details['center'] = (center_x, center_y, center_z)

                    # Calculate size (max distance from center)
                    max_dist = max(
                        ((coord[0] - center_x)**2 + (coord[1] - center_y)**2 + (coord[2] - center_z)**2)**0.5
                        for coord in coords
                    )
                    details['size'] = max_dist
                else:
                    details['center'] = (0, 0, 0)
                    details['size'] = 0

            return ligand_details if ligand_details else None

        except Exception as e:
            print(f"   ⚠️  Error analyzing ligands: {e}")
            return None

    def _display_ligand_analysis(self, ligand_details: Dict[str, Any]) -> None:
        """
        Display detailed ligand analysis for user review.
        
        Args:
            ligand_details (Dict): Detailed ligand information
        """
        try:
            print(f"\n   📊 LIGAND ANALYSIS:")
            print("   " + "-" * 70)
            print(f"   {'Ligand ID':<15} {'Type':<8} {'Chain':<6} {'ResNum':<7} {'Atoms':<7} {'Size':<8} {'Center (Å)':<15}")
            print("   " + "-" * 70)
            
            # Sort ligands by type (waters last) and then by size (largest first)
            sorted_ligands = sorted(
                ligand_details.items(),
                key=lambda x: (x[1]['is_water'], -x[1]['size'], x[1]['residue_name'])
            )
            
            for ligand_id, details in sorted_ligands:
                ligand_type = "Water" if details['is_water'] else "Ligand"
                center_str = f"({details['center'][0]:.1f}, {details['center'][1]:.1f}, {details['center'][2]:.1f})"
                
                print(f"   {ligand_id[:14]:<15} {ligand_type:<8} {details['chain_id']:<6} "
                    f"{details['residue_number']:<7} {details['atom_count']:<7} "
                    f"{details['size']:8.2f}Å {center_str:<15}")
            
            print("   " + "-" * 70)
            
            # Show summary
            total_ligands = len([l for l in ligand_details.values() if not l['is_water']])
            total_waters = len([l for l in ligand_details.values() if l['is_water']])
            
            print(f"   📊 Summary: {total_ligands} ligand(s), {total_waters} water molecule(s)")
            
            if total_ligands > 0:
                largest_ligand = max(
                    [l for l in ligand_details.values() if not l['is_water']],
                    key=lambda x: x['size']
                )
                print(f"   🎯 Largest ligand: {largest_ligand['residue_name']} ({largest_ligand['atom_count']} atoms, {largest_ligand['size']:.2f}Å)")
            
        except Exception as e:
            print(f"   ⚠️  Error displaying ligand analysis: {e}")

    def _select_ligands_interactive(self, ligand_details: Dict[str, Any]) -> Optional[List[str]]:
        """
        Interactive ligand selection interface.
        
        Args:
            ligand_details (Dict): Detailed ligand information
            
        Returns:
            Optional[List[str]]: List of selected ligand residue names or None if cancelled
        """
        try:
            # Separate waters from other ligands
            waters = {k: v for k, v in ligand_details.items() if v['is_water']}
            ligands = {k: v for k, v in ligand_details.items() if not v['is_water']}
            
            selected_residues = []
            
            # Handle non-water ligands first
            if ligands:
                print(f"\n   🧪 SELECT LIGANDS TO KEEP:")
                print("   " + "-" * 50)
                
                ligand_list = list(ligands.items())
                for i, (ligand_id, details) in enumerate(ligand_list, 1):
                    print(f"   {i}. {ligand_id} ({details['atom_count']} atoms, {details['size']:.2f}Å)")
                
                print(f"\n   Commands:")
                print(f"   • Enter numbers separated by commas (e.g., 1,3)")
                print(f"   • 'all' to keep all ligands")
                print(f"   • 'none' to remove all ligands")
                print(f"   • 'largest' to keep only the largest ligand")
                
                while True:
                    try:
                        selection = input(f"\n   💊 Select ligands to keep: ").strip().lower()
                        
                        if selection in ['cancel', 'quit', 'exit']:
                            return None
                        elif selection == 'all':
                            selected_residues.extend([details['residue_name'] for details in ligands.values()])
                            print(f"   ✅ Selected all ligands: {', '.join(set(selected_residues))}")
                            break
                        elif selection == 'none':
                            print("   🧹 No ligands selected")
                            break
                        elif selection == 'largest':
                            largest = max(ligands.values(), key=lambda x: x['size'])
                            selected_residues.append(largest['residue_name'])
                            print(f"   ✅ Selected largest ligand: {largest['residue_name']}")
                            break
                        else:
                            # Parse comma-separated numbers
                            try:
                                indices = [int(x.strip()) - 1 for x in selection.split(',')]
                                valid_selections = []
                                
                                for idx in indices:
                                    if 0 <= idx < len(ligand_list):
                                        ligand_id, details = ligand_list[idx]
                                        selected_residues.append(details['residue_name'])
                                        valid_selections.append(ligand_id)
                                    else:
                                        print(f"   ⚠️  Invalid selection: {idx + 1}")
                                
                                if valid_selections:
                                    print(f"   ✅ Selected: {', '.join(valid_selections)}")
                                    break
                                else:
                                    print("   ❌ No valid selections made")
                                    continue
                                    
                            except ValueError:
                                print("   ❌ Invalid input. Use numbers, 'all', 'none', or 'largest'")
                                continue
                                
                    except KeyboardInterrupt:
                        return None
            
            # Handle waters separately
            if waters:
                print(f"\n   💧 WATER MOLECULES FOUND: {len(waters)}")
                
                while True:
                    try:
                        water_choice = input(f"   Keep water molecules? (y/n/ask): ").strip().lower()
                        
                        if water_choice in ['y', 'yes']:
                            selected_residues.append('HOH')
                            print(f"   💧 Will keep {len(waters)} water molecules")
                            break
                        elif water_choice in ['n', 'no']:
                            print("   🧹 Will remove all water molecules")
                            break
                        elif water_choice == 'ask':
                            # Show water details and let user decide
                            print(f"   Water molecule details:")
                            water_list = list(waters.items())
                            for i, (water_id, details) in enumerate(water_list[:10], 1):  # Show first 10
                                center = details['center']
                                print(f"     {i}. {water_id} at ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
                            
                            if len(waters) > 10:
                                print(f"     ... and {len(waters) - 10} more water molecules")
                            
                            keep_waters = input(f"   Keep these water molecules? (y/n): ").strip().lower()
                            if keep_waters in ['y', 'yes']:
                                selected_residues.append('HOH')
                                print(f"   💧 Will keep {len(waters)} water molecules")
                            else:
                                print("   🧹 Will remove all water molecules")
                            break
                        else:
                            print("   ❌ Please answer 'y', 'n', or 'ask'")
                            continue
                            
                    except KeyboardInterrupt:
                        return None
            
            # Remove duplicates and return
            return list(set(selected_residues)) if selected_residues else []
            
        except Exception as e:
            print(f"   ❌ Error in ligand selection: {e}")
            return None

    def _display_pdb_analysis(self, analysis: Dict[str, Any]) -> None:
        """
        Display PDB file analysis results.
        
        Args:
            analysis (Dict): PDB analysis dictionary
        """
        try:
            print(f"\n📊 PDB FILE ANALYSIS")
            print("-" * 50)
            print(f"📁 File: {os.path.basename(analysis['file_path'])}")
            print(f"💾 Size: {analysis['file_size']:,} bytes")
            
            if analysis['header_info']:
                header = analysis['header_info']
                if header.get('title'):
                    print(f"🏷️  Title: {header['title']}")
                if header.get('pdb_id'):
                    print(f"🆔 PDB ID: {header['pdb_id']}")
                if header.get('deposition_date'):
                    print(f"📅 Date: {header['deposition_date']}")
            
            if analysis['resolution']:
                print(f"🔬 Resolution: {analysis['resolution']:.2f} Å")
            
            print(f"\n🧬 STRUCTURE CONTENT:")
            print(f"   ⚛️  Protein atoms: {analysis['atom_count']:,}")
            print(f"   🔗 Chains: {', '.join(analysis['chains'])} ({len(analysis['chains'])} total)")
            print(f"   📊 Residues: ~{analysis['residue_count']:,}")
            
            if analysis['hetero_count'] > 0:
                print(f"   🧪 Hetero atoms: {analysis['hetero_count']:,}")
                
                if analysis['has_waters']:
                    print(f"   💧 Contains water molecules")
                
                if analysis['has_ligands']:
                    print(f"   💊 Contains ligands: {', '.join(analysis['ligand_residues'])}")

            # Inform about alternate locations (altlocs)
            if analysis.get('has_altlocs', False):
                altlocs = analysis.get('altlocs', {})
                print(f"\n⚠️  Alternate locations (altlocs) detected:")
                print(f"   • Residues with altlocs: {altlocs.get('residues_with_altlocs', 0)}")
                print(f"   • Unique altloc IDs: {', '.join(altlocs.get('unique_altloc_ids', [])) or 'None'}")
                print(f"   • Total altloc atoms: {altlocs.get('total_altlocs', 0)}")
                print("   ⚠️  Consider handling altlocs (alternate conformations) during processing.")
            else:
                print(f"\n✅ No alternate location (altloc) atoms detected.")

            # Inform about gaps in the protein structure
            if analysis.get('total_gaps', 0) > 0:
                print(f"\n⚠️  Gaps detected in protein structure:")
                print(f"   • Number of gaps: {analysis['total_gaps']}")
                print(f"   • Residues involved in gaps: {len(analysis['gap_residues'])}")
                # Optionally, show a few gap details
                for chain_id, gaps in analysis.get('gaps', {}).items():
                    if gaps:
                        print(f"   • Chain {chain_id}: {len(gaps)} gap(s)")
                        for gap in gaps[:3]:  # Show up to 3 gaps per chain
                            print(f"      - Between residues {gap['start']} and {gap['end']}: missing {len(gap['gap_residues'])} residue(s)")
                        if len(gaps) > 3:
                            print(f"      ... {len(gaps) - 3} more gaps in chain {chain_id}")
            else:
                print(f"\n✅ No gaps detected in protein structure.")

            # Assessment
            print(f"\n✅ DOCKING SUITABILITY:")
            issues = []
            
            if analysis['atom_count'] == 0:
                issues.append("No protein atoms found")
            elif analysis['atom_count'] < 100:
                issues.append("Very small protein structure")
            
            if len(analysis['chains']) == 0:
                issues.append("No chain identifiers found")
            elif len(analysis['chains']) > 4:
                issues.append(f"Many chains ({len(analysis['chains'])}) - may need selection")
            
            if analysis['resolution'] and analysis['resolution'] > 3.0:
                issues.append(f"Low resolution ({analysis['resolution']:.1f} Å)")
            
            if issues:
                print("   ⚠️  Potential issues:")
                for issue in issues:
                    print(f"      • {issue}")
            else:
                print("   ✅ Structure appears suitable for docking")
            
            print("-" * 50)
            
        except Exception as e:
            print(f"⚠️  Error displaying PDB analysis: {e}")

    def _check_pdb_name_exists(self, receptor_name: str) -> bool:
        """
        Check if receptor name already exists in the database.
        
        Args:
            receptor_name (str): Receptor name to check
            
        Returns:
            bool: True if name exists, False otherwise
        """
        try:
            receptors_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
            
            if not os.path.exists(receptors_db_path):
                return False
            
            import sqlite3
            conn = sqlite3.connect(receptors_db_path)
            cursor = conn.cursor()
            
            cursor.execute("SELECT COUNT(*) FROM receptors WHERE receptor_name = ?", (receptor_name,))
            count = cursor.fetchone()[0]
            
            conn.close()
            return count > 0
            
        except Exception:
            return False

    def _create_receptor_folder(self, base_dir: str, receptor_name: str) -> str:
        """
        Create folder structure for receptor files.
        
        Args:
            base_dir (str): Base receptors directory
            receptor_name (str): Name of the receptor
            
        Returns:
            str: Path to created receptor folder
        """
        try:
            # Clean receptor name for folder creation
            clean_name = "".join(c for c in receptor_name if c.isalnum() or c in ('_', '-', '.')).strip()
            receptor_folder = os.path.join(base_dir, clean_name)
            
            # Create receptor folder and subdirectories
            os.makedirs(receptor_folder, exist_ok=True)
            
            subdirs = [
                'original',      # Original PDB file
                'processed',     # Processed receptor files (PDBQT, etc.)
                'grid_files',    # Grid/map files
                'analysis',      # Analysis results
                'logs'           # Processing logs
            ]
            
            for subdir in subdirs:
                os.makedirs(os.path.join(receptor_folder, subdir), exist_ok=True)
            
            print(f"   📁 Created receptor folder: {receptor_folder}")
            return receptor_folder
            
        except Exception as e:
            print(f"❌ Error creating receptor folder: {e}")
            return ""

    def _select_chain_interactive(self, chains: List[str]) -> Optional[str]:
        """
        Interactive chain selection for multi-chain structures.
        
        Args:
            chains (List[str]): Available chain identifiers
            
        Returns:
            Optional[str]: Selected chain ID or None if cancelled
        """
        try:
            print(f"   📋 Available chains: {', '.join(chains)}")
            
            while True:
                try:
                    selection = input(f"   🔗 Select chain ({'/'.join(chains)}): ").strip().upper()
                    
                    if selection in chains:
                        return selection
                    elif selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    else:
                        print(f"   ❌ Invalid chain. Please select from: {', '.join(chains)}")
                        continue
                        
                except KeyboardInterrupt:
                    return None
                    
        except Exception:
            return None

    def _extract_single_chain(self, pdb_file: str, chain_id: str) -> str:
        """
        Extract single chain from PDB file, including 'TER' cards associated with that chain.
        
        Args:
            pdb_file (str): Input PDB file path
            chain_id (str): Chain identifier to extract
            
        Returns:
            str: Path to single-chain PDB file
        """
        try:
            output_file = pdb_file.replace('.pdb', f'_chain_{chain_id}.pdb')
            with open(pdb_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    
                    if line.startswith(('ATOM', 'HETATM')):
                        if line[21].strip() == chain_id:
                            outfile.write(line)
                    elif line.startswith('TER'):
                        # For TER, the chain ID is at column 22 (index 21)
                        print(line.strip())   
                        if len(line) > 21 and line[21].strip() == chain_id:
                            outfile.write(line)
                    else:
                        # Keep header and other non-atom lines
                        outfile.write(line)
            
            return output_file
            
        except Exception as e:
            print(f"   ❌ Error extracting chain: {e}")
            return ""

    def _create_pdb_registry_entry(self, pdb_name: str, original_pdb_path: str, 
                                    processed_pdb_path: str, pdb_folder: str, 
                                    analysis: Dict[str, Any], checked_pdb_path) -> Optional[Dict[str, Any]]:
        """
        Create PDB templates registry entry in pdbs.db

        Args:
            pdb_name (str): Name of the PDB
            original_pdb_path (str): Path to original PDB file
            processed_pdb_path (str): Path to processed PDB file
            pdb_folder (str): Path to PDB folder
            analysis (Dict): PDB analysis results

        Returns:
            Optional[Dict[str, Any]]: Registry entry information or None if failed
        """
        try:
            # Create receptors database
            pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')

            import sqlite3
            import json
            from datetime import datetime

            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()

            # Create receptors table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS pdb_templates (
                    pdb_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pdb_name TEXT UNIQUE NOT NULL,
                    project_name TEXT,
                    original_pdb_path TEXT NOT NULL,
                    processed_pdb_path TEXT NOT NULL,
                    checked_pdb_path TEXT NOT NULL,
                    pdb_folder_path TEXT NOT NULL,
                    pdb_analysis TEXT,
                    chains TEXT,
                    resolution REAL,
                    atom_count INTEGER,
                    has_ligands BOOLEAN,
                    status TEXT DEFAULT 'created',
                    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    last_modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    notes TEXT
                )
            ''')

            # Prepare data for insertion
            pdb_analysis_json = json.dumps(analysis, indent=2)
            chains_str = ','.join(analysis['chains'])

            # Check for identical registry (excluding 'notes')
            cursor.execute('''
                SELECT pdb_id, pdb_name, project_name, original_pdb_path, processed_pdb_path, checked_pdb_path,
                       pdb_folder_path, pdb_analysis, chains, resolution, atom_count, has_ligands, status
                FROM pdb_templates
            ''')
            new_entry = (
                pdb_name,
                self.name,
                original_pdb_path,
                processed_pdb_path,
                checked_pdb_path,
                pdb_folder,
                pdb_analysis_json,
                chains_str,
                analysis.get('resolution'),
                analysis['atom_count'],
                analysis['has_ligands'],
                'created'
            )
            for row in cursor.fetchall():
                # Compare all fields except notes, created_date, last_modified, and pdb_id
                existing_entry = row[1:13]  # skip pdb_id
                if existing_entry == new_entry:
                    print(f"⚠️  An identical receptor registry already exists (PDB name: {row[1]}).")
                    print("   No new registry will be created.")
                    conn.close()
                    return None

            # Prompt for notes
            try:
                notes = input("Enter a note for this PDB file (Enter: empty note): ").strip()
                if not notes:
                    notes = ""
            except KeyboardInterrupt:
                print("\n   ⚠️  Note entry cancelled. Using empty note.")
                notes = ""

            # Insert or update receptor entry
            cursor.execute('''
                INSERT OR REPLACE INTO pdb_templates (
                    pdb_name, project_name, original_pdb_path, processed_pdb_path, checked_pdb_path,
                    pdb_folder_path, pdb_analysis, chains, resolution, atom_count,
                    has_ligands, status, notes
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                pdb_name,
                self.name,
                original_pdb_path,
                processed_pdb_path,
                checked_pdb_path,
                pdb_folder,
                pdb_analysis_json,
                chains_str,
                analysis.get('resolution'),
                analysis['atom_count'],
                analysis['has_ligands'],
                'created',
                notes if notes is not None else f"PDB model created from PDB file: {os.path.basename(original_pdb_path)}"
            ))

            pdb_id = cursor.lastrowid or cursor.execute(
                "SELECT pdb_id FROM pdb_templates WHERE pdb_name = ?", 
                (pdb_name,)
            ).fetchone()[0]

            conn.commit()
            conn.close()

            return {
                'pdb_id': pdb_id,
                'pdb_name': pdb_name,
                'project_name': self.name,
                'original_pdb_path': original_pdb_path,
                'processed_pdb_path': processed_pdb_path,
                'checked_pdb_path': checked_pdb_path,
                'pdb_folder_path': pdb_folder,
                'database_path': pdbs_db_path,
                'pdb_analysis': analysis,
                'created_date': datetime.now().isoformat(),
                'notes': notes
            }

        except Exception as e:
            print(f"❌ Error creating PDB registry: {e}")
            return None

    def _display_pdb_creation_summary(self, pdb_registry: Dict[str, Any]) -> None:
        """
        Display summary of created receptor.
        
        Args:
            receptor_info (Dict): Receptor registry information
        """
        try:
            print(f"\n✅ PDB FILE CREATED SUCCESSFULLY")
            print("=" * 60)
            print(f"🆔 PDB ID: {pdb_registry['pdb_id']}")
            print(f"🏷️  PDB Name: {pdb_registry['pdb_name']}")
            print(f"📁 Original PDB: {os.path.basename(pdb_registry['original_pdb_path'])}")
            print(f"📋 Processed PDB: {os.path.relpath(pdb_registry['processed_pdb_path'])}")
            print(f"📂 PDB Folder: {os.path.relpath(pdb_registry['pdb_folder_path'])}")
            print(f"💾 Database: {os.path.relpath(pdb_registry['database_path'])}")
            
            analysis = pdb_registry['pdb_analysis']
            print(f"\n📊 STRUCTURE INFO:")
            print(f"   ⚛️  Atoms: {analysis['atom_count']:,}")
            print(f"   🔗 Chains: {', '.join(analysis['chains'])}")
            if analysis.get('resolution'):
                print(f"   🔬 Resolution: {analysis['resolution']:.2f} Å")
            if analysis['has_ligands']:
                print(f"   💊 Ligands: {', '.join(analysis['ligand_residues'])}")
            
            print(f"\n💡 NEXT STEPS:")
            print("   1. Prepare receptor for docking (convert to PDBQT format)")
            print("   2. Generate grid files for your docking software")
            print("   3. Use this receptor in docking assays")
            print("   4. Configure binding site and docking parameters")
            
            print("=" * 60)
            
        except Exception as e:
            print(f"⚠️  Error displaying receptor summary: {e}")
        
    def _add_missing_atoms_in_receptor_tleap(self, receptor_file, receptor_info):
        """
        Will process using tleap a receptor_file by loading the .pdb and saving to a new file.
        Will use tleap as installed in the environment
        """
        import subprocess
        import tempfile
        import os
        import shutil

        try:
            # Check tleap is available
            tleap_path = shutil.which("tleap")
            if tleap_path is None:
                print("❌ tleap not found in PATH. Please install AmberTools and ensure tleap is available.")
                return None

            # Prepare output file path
            receptor_dir = os.path.dirname(receptor_file)
            output_pdb = os.path.join(receptor_dir, "receptor_checked.pdb")

            # Prepare tleap input script
            tleap_script = f"""
source leaprc.protein.ff14SB
mol = loadpdb "{receptor_file}"
savepdb mol "{output_pdb}"
quit
"""
            # Write tleap script to a temp file
            with tempfile.NamedTemporaryFile("w", delete=False, suffix=".in") as tf:
                tf.write(tleap_script)
                tleap_input_path = tf.name


            ## This first tleap run will ensure addition of missing heavy atoms
            print(f"🛠️  Running tleap to add missing atoms to receptor...")
            result = subprocess.run(
                ["tleap", "-f", tleap_input_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Clean up temp file
            os.remove(tleap_input_path)

            # Check for errors
            if result.returncode != 0:
                print("❌ tleap failed to process the receptor file.")
                print(result.stderr)
                return None

            # Check if output file was created
            if not os.path.exists(output_pdb):
                print("❌ tleap did not produce the expected output file.")
                return None

            # Hydrogen atoms will deleted using pdb4amber available in the environment
            intermediate_receptor_file = receptor_dir + "/receptor_checked.pdb"
            output_receptor_file_intermediate = receptor_dir + "/receptor_checked_woHs.pdb"

            result = subprocess.run(
                ["pdb4amber", f"{intermediate_receptor_file}", "-y", "-o", f"{output_receptor_file_intermediate}"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Move the intermediate file (without hydrogens) to the final output path
            shutil.move(output_receptor_file_intermediate, output_pdb)

#             # A second round of tleap processing will be performed on the receptor containing all heavy atoms

#             # Prepare tleap input script. The original output_pdb will be overwritten.
#             tleap_script = f"""
# source leaprc.protein.ff14SB
# mol = loadpdb "{output_receptor_file_intermediate}"
# savepdb mol "{output_pdb}"
# quit
# """

#             # Write tleap script to a temp file
#             with tempfile.NamedTemporaryFile("w", delete=False, suffix=".in") as tf:
#                 tf.write(tleap_script)
#                 tleap_input_path = tf.name

#             ## This second tleap run will ensure addition of hydrogens on all added heavy atoms
#             result = subprocess.run(
#                 ["tleap", "-f", tleap_input_path],
#                 stdout=subprocess.PIPE,
#                 stderr=subprocess.PIPE,
#                 text=True
#             )

#             # Clean up temp file
#             os.remove(tleap_input_path)

#             # Check for errors
#             if result.returncode != 0:
#                 print("❌ tleap failed to process the receptor file.")
#                 print(result.stderr)
#                 return None

#             # Check if output file was created
#             if not os.path.exists(output_pdb):
#                 print("❌ tleap did not produce the expected output file.")
#                 return None

            # Reassign chain IDs and residue numbers
            self._reassign_chain_and_resnums(receptor_file, output_pdb)

            print(f"✅ tleap processing complete. Output: {output_pdb}")
            return output_pdb

        except Exception as e:
            print(f"❌ Error running tleap for missing atom addition: {e}")
            return None

    def _add_missing_atoms_in_pdb_prody(self, receptor_file, receptor_info):

        import prody

        PDBname = receptor_file
        receptor_path = os.path.dirname(receptor_file)
        print(f"folder: {receptor_path}")
        final_file = receptor_file.replace(".pdb","_checked.pdb")
        
        ## Add missing atoms (including heavy and hydrogens) using prody
        print("🔨 Adding missing atoms...")
        structure = prody.addMissingAtoms(PDBname, method='pdbfixer')

        ## Rename the default output of prody 
        shutil.move(receptor_path + '/addH_receptor.pdb', final_file)

        return final_file

    def _add_missing_atoms_in_pdb_tleap(self, receptor_file, receptor_info):
        
        import tempfile
        import subprocess
        
        print(f"Processing pdb: \n {receptor_file} \n with tleap")
        
        ## Write a tleap input to load receptor_file
        try:
            tleap_path = shutil.which("tleap")

            # Fail is tleap is not available
            if tleap_path is None:
                print("❌ tleap not found in PATH. Please install AmberTools and ensure tleap is available.")
                return None
        
            # Process with tleap
            tleap_processed_file = receptor_file.replace(".pdb", "_checked.pdb")
            tleap_script = f"""
source leaprc.protein.ff14SB
source leaprc.water.tip3p
HOH = WAT
mol = loadpdb {receptor_file}
savepdb mol {tleap_processed_file}
quit
"""
            # Write tleap script to a temp file
            with tempfile.NamedTemporaryFile("w", delete=False, suffix=".in") as tf:
                tf.write(tleap_script)
                tleap_input_path = tf.name

            ## Run tleap to process the receptor
            print(f"🛠️  Running tleap to process receptor...")
            result = subprocess.run(
                ["tleap", "-f", tleap_input_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            ## Rename 'WAT' residues to 'HOH' in the tleap_processed_file
            with open(tleap_processed_file, 'r') as file:
                filedata = file.read()

            filedata = filedata.replace(' WAT ', ' HOH ')

            with open(tleap_processed_file, 'w') as file:
                file.write(filedata)

            # Clean up temp file
            os.remove(tleap_input_path)

            # Check for errors
            if result.returncode != 0:
                print("❌ tleap failed to process the receptor file.")
                print(result.stderr)
                return None

            # Check if output file was created
            if not os.path.exists(tleap_processed_file):
                print("❌ tleap did not produce the expected output file.")
                return None

            return tleap_processed_file
                    
        except Exception as e:
            print(f"Failed processing: {receptor_file} with tleap")
            print(e)
            return None
        
    def _reassign_chain_and_resnums(self, reference_pdb, target_pdb):
        """
        Create a dictionary in which the key is formed by residue name, the chain and the residue number in the reference pdb_file, and the value is formed by the residue name and residue number in the target pdb file.
        Then, process the target file so that each line is updated to use the chain id and residue number from the reference pdb (from the key), while keeping the residue name from the target pdb (from the value).
        """
        
        try:
            # Step 1: Build mapping dictionary by residue order
            ref_residues = []
            ref_coords = []
            with open(reference_pdb, 'r') as ref:
                for line in ref:
                    if line.startswith(('ATOM', 'HETATM')):
                        ref_resname = line[17:20].strip()
                        ref_chain = line[21].strip()
                        ref_resnum = line[22:26].strip()
                        ref_name = line[12:16].strip()
                        ref_occupancy = line[54:60].strip()
                        ref_tempFactor = line[60:66].strip().replace('\n', '')
                        ref_element = line[76:78].strip().replace('\n', '')
                        
                        # Extract atom coordinates
                        try:
                            ref_x = float(line[30:38].strip())
                            ref_y = float(line[38:46].strip())
                            ref_z = float(line[46:54].strip())
                        except ValueError:
                            ref_x, ref_y, ref_z = None, None, None
                        key = (
                            ref_resname, ref_chain, ref_resnum,
                            ref_name, ref_occupancy, ref_tempFactor, ref_element,
                            ref_x, ref_y, ref_z
                        )
                        if (ref_x, ref_y, ref_z) not in ref_coords:
                            ref_coords.append((ref_x, ref_y, ref_z))  # Add coordinates tuple
                        if key not in ref_residues:
                            ref_residues.append(key)

            tgt_residues = []
            tgt_coords = []  # New list to store (tgt_x, tgt_y, tgt_z) tuples
            with open(target_pdb, 'r') as tgt:
                for line in tgt:
                    if line.startswith(('ATOM', 'HETATM')):
                        tgt_resname = line[17:20].strip()
                        tgt_resnum = line[22:26].strip()
                        tgt_name = line[12:16].strip()
                        # Extract atom coordinates
                        try:
                            tgt_x = float(line[30:38].strip())
                            tgt_y = float(line[38:46].strip())
                            tgt_z = float(line[46:54].strip())
                        except ValueError:
                            tgt_x, tgt_y, tgt_z = None, None, None
                        key = (tgt_resname, tgt_resnum, tgt_name, tgt_x, tgt_y, tgt_z)
                        if (tgt_x, tgt_y, tgt_z) not in tgt_coords:
                            tgt_coords.append((tgt_x, tgt_y, tgt_z))  # Add coordinates tuple
                        if key not in tgt_residues and (tgt_x, tgt_y, tgt_z) in ref_coords:
                            tgt_residues.append(key)
                            
            mapping = {}
            
            # Build mapping by matching coordinates (x, y, z) between reference and target
            for ref_key in ref_residues:
                ref_x, ref_y, ref_z = ref_key[-3:]
                ref_coord_tuple = (ref_x, ref_y, ref_z)
                for tgt_key in tgt_residues:
                    tgt_x, tgt_y, tgt_z = tgt_key[-3:]
                    tgt_coord_tuple = (tgt_x, tgt_y, tgt_z)
                    if ref_coord_tuple == tgt_coord_tuple:
                        mapping[ref_key] = tgt_key
                        break

            # Step 2: Process target file and update chain id and residue number
            import shutil
            temp_out = target_pdb + ".chainfix"
            tgt_res_idx = 0
            ref_keys_list = list(mapping.keys())
            with open(target_pdb, 'r') as tgt, open(temp_out, 'w') as out:
                for line in tgt:
                    if line.startswith(('ATOM', 'HETATM')):
                        tgt_resname = line[17:20].strip()
                        tgt_resnum = line[22:26].strip()
                        tgt_name = line[12:16].strip()
                        # Extract atom coordinates
                        try:
                            tgt_x = float(line[30:38].strip())
                            tgt_y = float(line[38:46].strip())
                            tgt_z = float(line[46:54].strip())
                        except ValueError:
                            tgt_x, tgt_y, tgt_z = None, None, None
                        tgt_key = (tgt_resname, tgt_resnum, tgt_name, tgt_x, tgt_y, tgt_z)
                        
                        # Find the corresponding ref_key by matching tgt_key in mapping values
                        ref_key = None
                        
                        for k, v in mapping.items():
                            if v == tgt_key:
                                ref_key = k
                                break
                        
                        if ref_key is not None:
                            ref_resname, ref_chain, ref_resnum, ref_name, ref_occupancy, ref_tempFactor, ref_element, ref_x, ref_y, ref_z = ref_key
                            record = line[0:6]
                            atom_serial = line[6:11]
                            atom_name = line[12:16]
                            altLoc = " "
                            resName = line[17:20]
                            insertion = line[20:21]
                            chainID = ref_chain
                            resSeq = ref_resnum
                            x = line[30:38]
                            y = line[38:46]
                            z = line[46:54]
                            # Set occupancy and tempFactor to default values
                            occupancy = ref_occupancy
                            tempFactor = ref_tempFactor
                            # Use element and charge if present
                            element = ref_element
                            charge = line[78:80] if len(line) >= 80 else "  "
                            # Reconstruct the line in strict PDB format
                            # Format fields to strict PDB format (columns)
                            # Columns: 
                            #  1-6   Record name   (ATOM  / HETATM)
                            #  7-11  Atom serial   (right, 5)
                            # 13-16  Atom name     (left, 4)
                            # 17     altLoc        (1)
                            # 18-20  resName       (right, 3)
                            # 22     chainID       (1)
                            # 23-26  resSeq        (right, 4)
                            # 27     iCode         (1)
                            # 31-38  x             (8.3f)
                            # 39-46  y             (8.3f)
                            # 47-54  z             (8.3f)
                            # 55-60  occupancy     (6.2f)
                            # 61-66  tempFactor    (6.2f)
                            # 77-78  element       (right, 2)
                            # 79-80  charge        (right, 2)
                            new_line = (
                                f"{record:<6}{int(atom_serial):5d} "
                                f"{atom_name:<4}{altLoc:1}"
                                f"{resName:>3}{insertion:1}"
                                f"{chainID:1}"
                                f"{int(resSeq):4d}"
                                f"    "
                                f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
                                f"{float(occupancy):6.2f}{float(tempFactor):6.2f}"
                                f"          "
                                f"{element:>2}{charge:>2}\n"
                            )
                            
                            # Write the new_line to a temporary file for debugging
                            with open("debug_chainfix_output_heavy.txt", "a") as debug_out:
                                debug_out.write(new_line)
                            
                            out.write(new_line)
                            tgt_res_idx += 1
                        else:
                            # If no mapping, parse the line and reconstruct using PDB fields
                            # Parse fields from the original line
                            record = line[0:6]
                            atom_serial = line[6:11]
                            atom_name = line[12:16]
                            altLoc = line[16]
                            resName = line[17:20]
                            insertion = line[20:21]
                            #chainID = line[21:22]
                            chainID = chainID
                            #resSeq = line[22:26]
                            resSeq = resSeq
                            x = line[30:38]
                            y = line[38:46]
                            z = line[46:54]
                            # Set occupancy and tempFactor to default values
                            occupancy = f"{1.00:6.2f}"
                            tempFactor = f"{0.00:6.2f}"
                            # Use element and charge if present
                            element = line[76:78] if len(line) >= 78 else " H"
                            charge = line[78:80] if len(line) >= 80 else "  "
                            # Reconstruct the line in strict PDB format
                            new_line = (
                                f"{record:<6}{int(atom_serial):5d} "
                                f"{atom_name:<4}{altLoc:1}"
                                f"{resName:>3}{insertion:1}"
                                f"{chainID:1}"
                                f"{int(resSeq):4d}"
                                f"    "
                                f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
                                f"{float(occupancy):6.2f}{float(tempFactor):6.2f}"
                                f"          "
                                f"{element:>2}{charge:>2}\n"
                            )
                            
                            # Write the new_line to a temporary file for debugging
                            with open("debug_chainfix_output.txt", "a") as debug_out:
                                debug_out.write(new_line)
                            
                            out.write(new_line)
                    else:
                        out.write(line)
            shutil.move(temp_out, target_pdb)
            print("✅ Chain IDs and residue numbers reassigned in target PDB using reference mapping.")
            return mapping
        except Exception as e:
            print(f"   ⚠️  Error in residue mapping: {e}")
            return {}
    
    def create_receptor_for_docking(self, selection: int = None):
        """
        Will load the pdbs.db database from project_path/docking/receptors and prompt the user to select one of the PBD models available to create a pdbqt file for docking purposes using helper function to be created.
        If 'selection' is provided, it will be used directly; otherwise, the user will be prompted.
        """
        import sqlite3

        pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
        if not os.path.exists(pdbs_db_path):
            print(f"❌ No pdbs database found at {pdbs_db_path}")
            return None

        try:
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            cursor.execute('''
                SELECT pdb_id, pdb_name, processed_pdb_path, checked_pdb_path, pdb_folder_path, notes
                FROM pdb_templates
                ORDER BY created_date ASC
            ''')
            pdbs = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error loading PDB files database: {e}")
            return None

        if not pdbs:
            print("❌ No PDB files found in the database.")
            return None

        print("\n📋 Available PDB files:")
        print("=" * 60)
        for i, (rid, name, processed_pdb, checked_pdb, folder, notes) in enumerate(pdbs, 1):
            print(f"{i}. {name} (ID: {rid})")
            print(f"   Processed PDB: {os.path.basename(processed_pdb)}")
            print(f"   Checked PDB: {os.path.basename(checked_pdb)}")
            print(f"   Folder: {os.path.relpath(folder)}")
            print(f"   Notes: {notes}")
        print("=" * 60)

        # Prompt user to select a receptor only if selection is not provided
        if selection is None:
            while True:
                try:
                    selection_input = input("Select the pdb file by number (or 'cancel'): ").strip()
                    if selection_input.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ PDB file selection cancelled.")
                        return None
                    try:
                        idx = int(selection_input) - 1
                        if 0 <= idx < len(pdbs):
                            selected = pdbs[idx]
                            break
                        else:
                            print(f"❌ Invalid selection. Enter a number between 1 and {len(pdbs)}.")
                    except ValueError:
                        print("❌ Please enter a valid number.")
                except KeyboardInterrupt:
                    print("\n❌ PDB selection cancelled.")
                    return None
        else:
            idx = selection - 1
            if 0 <= idx < len(pdbs):
                selected = pdbs[idx]
            else:
                print(f"❌ Invalid selection. Enter a number between 1 and {len(pdbs)}.")
                return None

        pdb_id, pdb_name, processed_pdb_path, checked_pdb_path, rec_main_path, notes = selected

        # Select the pdb file to convert to .pdbqt
        pdb_to_convert = checked_pdb_path

        print(f"PDB to convert: {pdb_to_convert}")

        # Call helper to create PDBQT
        pdbqt_file, configs = self._create_pdbqt_from_pdb(pdb_to_convert)

        # Check if pdbqt receptor file contains a Zinc ion in order to apply AutoDock4Zn parameters
        has_ZN = self._check_ZN_presence(pdbqt_file, rec_main_path)

        if has_ZN:
            
            print("⚠️  Zn atom(s) detected in the receptor structure.")
            while True:
                zn_choice = input("Do you want to process with AutoDock4Zn parameters? (y/n): ").strip().lower()
                if zn_choice in ['y', 'yes']:
                    print("✅ AutoDock4Zn parameters will be applied.")
                    self._apply_ad4zn_params(pdbqt_file, rec_main_path)
                    break

                elif zn_choice in ['n', 'no']:
                    print("ℹ️  Skipping AutoDock4Zn processing.")
                    has_ZN = False
                    break
                else:
                    print("❌ Please answer 'y' or 'n'.")

            print("✅ AutoDock4Zn parameters to be applied.")

        # Call helper method to create receptor grids
        self._create_receptor_grids(pdbqt_file, configs, rec_main_path, has_ZN)

        if pdbqt_file:

            self._create_receptor_register(pdb_id, pdb_name, pdb_to_convert, pdbqt_file, configs)
            
            print(f"✅ PDBQT file created: {pdbqt_file}")
        
        else:
            print("❌ Failed to create PDBQT file.")
            return None
        
    def _create_pdbqt_from_pdb(self, pdb_to_convert):
        
        from meeko import MoleculePreparation, ResidueChemTemplates
        from meeko import Polymer

    
        destination_pdbqt = pdb_to_convert.replace(".pdb",".pdbqt")

        structure = self._create_prody_selection(pdb_to_convert)
        
        configs = self._prompt_recprep_opts()
        
        mk_prep = MoleculePreparation()
        chem_templates = ResidueChemTemplates.create_from_defaults()
        mypol = Polymer.from_prody(structure, chem_templates, mk_prep)

        self._write_pdbqt(mypol, destination_pdbqt)
        
        if os.path.exists(destination_pdbqt):
            return destination_pdbqt, configs
        else:
            print(f"❌ Failed to create PDBQT file at {destination_pdbqt}")
            return False

    def _create_prody_selection(self, pdb_to_convert):
        
        from prody import parsePDB
        
        structure = parsePDB(pdb_to_convert)
        
        return structure
    
    def _prompt_recprep_opts(self) -> Optional[dict]:
        """
        Prompt the user to enter docking box center coordinates and box size for receptor preparation.

        This method interactively asks the user for the X, Y, Z coordinates of the box center and the
        X, Y, Z dimensions (size) of the docking box. It validates the input and returns a dictionary
        with the box configuration, suitable for use in docking setup and registry.

        Returns:
            dict: Dictionary with keys 'center' and 'size', each mapping to a dict with 'x', 'y', 'z' values.
                Example: {'center': {'x': 10.0, 'y': 20.0, 'z': 30.0}, 'size': {'x': 20, 'y': 20, 'z': 20}}
            None: If the user cancels the input.

        Raises:
            None: All exceptions are handled internally with user feedback.
        """
        # Prompt for box center coordinates
        while True:
            try:
                x = input("Enter X coordinate for box center (e.g., 12.34): ").strip()
                y = input("Enter Y coordinate for box center (e.g., -56.78): ").strip()
                z = input("Enter Z coordinate for box center (e.g., 90.12): ").strip()
                try:
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    coords = {'x': round(x, 2), 'y': round(y, 2), 'z': round(z, 2)}
                    print(f"Coordinates entered: ({coords['x']:.2f}, {coords['y']:.2f}, {coords['z']:.2f})")
                except ValueError:
                    print("❌ Please enter valid numbers for all coordinates (e.g., 12.34)")
                    continue
                break
            except KeyboardInterrupt:
                print("\n❌ Coordinate entry cancelled.")
                return None

        # Prompt for box dimensions
        while True:
            try:
                x = input("Enter X box size (e.g., 20): ").strip()
                y = input("Enter Y box size (e.g., 20): ").strip()
                z = input("Enter Z box size (e.g., 20): ").strip()
                try:
                    x = int(x)
                    y = int(y)
                    z = int(z)
                    sizes = {'x': x, 'y': y, 'z': z}
                    print(f"Sizes entered: ({x}, {y}, {z})")
                except ValueError:
                    print("❌ Please enter valid integers for all box dimensions")
                    continue
                break
            except KeyboardInterrupt:
                print("\n❌ Box size entry cancelled.")
                return None

        # Prompt for dielectric constant
        while True:
            try:
                dielectric = input("Enter dielectric constant to use (default: -42): ").strip()
                if not dielectric:
                    dielectric = -42
                else:
                    try:
                        dielectric = int(dielectric)
                    except ValueError:
                        print("❌ Please enter a valid integer for the dielectric constant")
                        continue
                try:
                    dielectric = int(dielectric)
                    misc = {'dielectric': dielectric}
                    print(f"Dielectric constant: ({dielectric})")
                except ValueError:
                    print("❌ Please enter valid integers for dielectric constant")
                    continue
                break
            except KeyboardInterrupt:
                print("\n❌ Dielectric constant entry cancelled.")
                return None

        # Prompt for smooth constant
        while True:
            try:
                smooth = input("Enter smooth constant to use (default: 0.5): ").strip()
                if not smooth:
                    smooth = 0.5
                else:
                    try:
                        smooth = float(smooth)
                    except ValueError:
                        print("❌ Please enter a valid float for the smooth constant")
                        continue
                try:
                    smooth = float(smooth)
                    misc = {'smooth': smooth}
                    print(f"Smooth value: ({smooth})")
                except ValueError:
                    print("❌ Please enter valid float for smooth value")
                    continue
                break
            except KeyboardInterrupt:
                print("\n❌ Smooth value entry cancelled.")
                return None

        # Prompt for spacing constant
        while True:
            try:
                spacing = input("Enter spacing to use (default: 0.375): ").strip()
                if not spacing:
                    spacing = 0.375
                else:
                    try:
                        spacing = float(spacing)
                    except ValueError:
                        print("❌ Please enter a valid float for the spacing constant")
                        continue
                try:
                    spacing = float(spacing)
                    misc = {'spacing': spacing}
                    print(f"Spacing value: ({spacing})")
                except ValueError:
                    print("❌ Please enter valid float for spacing value")
                    continue
                break
            except KeyboardInterrupt:
                print("\n❌ Spacing value entry cancelled.")
                return None

        return {'center': coords, 'size': sizes, 'dielectric': dielectric, 'smooth': smooth, 'spacing': spacing}

    def _write_pdbqt(self, mypol, destination_pdbqt: str) -> None:
        """
        Write a PDBQT file from a Meeko Polymer object.

        This method uses Meeko's PDBQTWriterLegacy to generate the PDBQT string from the provided
        Polymer object and writes it to the specified destination file. Only the rigid PDBQT string
        is written (no flexible residues).

        Args:
            mypol: Meeko Polymer object representing the prepared receptor or ligand.
            destination_pdbqt (str): Path to the output .pdbqt file.

        Returns:
            None

        Raises:
            Exception: If writing the file fails.
        """
        try:
            from meeko import PDBQTWriterLegacy

            rigid_pdbqt_string, _ = PDBQTWriterLegacy.write_string_from_polymer(mypol)
            with open(destination_pdbqt, "w") as f:
                f.write(rigid_pdbqt_string)
            print(f"   ✅ PDBQT file written: {destination_pdbqt}")
        except Exception as e:
            print(f"   ❌ Error writing PDBQT file: {e}")
            raise

    def _create_receptor_grids(self, pdbqt_file, configs, rec_main_path, has_ZN):
        """
        Generate AutoGrid4 grid parameter files (.gpf) for the prepared receptor.

        This method writes the grid parameter file using Meeko's get_gpf_string, applies Zn parameters if needed,
        and runs AutoGrid4 to generate grid maps. Handles Zn-specific parameter file and grid settings.

        Args:
            pdbqt_file (str): Path to the receptor PDBQT file.
            configs (dict): Docking box configuration (center, size, dielectric, smooth, spacing).
            rec_main_path (str): Path to the receptor's main folder.
            has_ZN (bool): Whether the receptor contains a Zn atom.

        Returns:
            None
        """
        import os
        import shutil
        import subprocess
        import site
        from meeko.gridbox import get_gpf_string

        # Atom types
        ligand_types = ["HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S", "Cl", "Br", "I"]
        receptor_types = ligand_types[:-1] + (["ZN", "TZ", "Zn"] if has_ZN else ["Zn"]) # Iodine is removed from receptor types because grids fails to compute

        # Box center and dimensions
        grid_center = list(configs['center'].values())
        box_dims = list(configs['size'].values())

        # Prepare grid files directory
        grids_path = os.path.join(rec_main_path, 'grid_files')
        os.makedirs(grids_path, exist_ok=True)
        grids_file_path = os.path.join(grids_path, 'receptor.gpf')
        grid_log_file = grids_file_path.replace(".gpf", ".glg")

        # Write .gpf file
        with open(grids_file_path, "w") as f:
            gpf_string, _ = get_gpf_string(
                grid_center, box_dims, pdbqt_file, receptor_types, ligand_types,
                dielectric=configs['dielectric'], smooth=configs['smooth'], spacing=configs['spacing']
            )
            f.write(gpf_string)

            # If Zn present, copy AD4Zn.dat and add Zn grid parameters
            if has_ZN:
                try:
                    ad4zn_src = os.path.join(site.getsitepackages()[0], "tidyscreen/config/AD4Zn.dat")
                    shutil.copy(ad4zn_src, grids_path)
                    print(f"✅ Copied Zn parameters file: {ad4zn_src} → {grids_path}")
                except Exception as e:
                    print(f"⚠️  Error copying Zn parameters file: {e}")

                # Write Zn-specific grid parameters
                zn_params = [
                    "nbp_r_eps 0.25 23.2135 12 6 NA TZ",
                    "nbp_r_eps 2.1   3.8453 12 6 OA Zn",
                    "nbp_r_eps 2.25  7.5914 12 6 SA Zn",
                    "nbp_r_eps 1.0   0.0    12 6 HD Zn",
                    "nbp_r_eps 2.0   0.0060 12 6 NA Zn",
                    "nbp_r_eps 2.0   0.2966 12 6  N Zn"
                ]
                for line in zn_params:
                    f.write(line + "\n")

        # Prepend parameter file instruction if Zn present
        if has_ZN:
            try:
                with open(grids_file_path, 'r') as f:
                    original_content = f.read()
                with open(grids_file_path, 'w') as f:
                    f.write("parameter_file AD4Zn.dat\n" + original_content)
            except Exception as e:
                print(f"⚠️  Error prepending Zn parameter instruction: {e}")

        # Run AutoGrid4
        print(f"🛠️  Running AutoGrid4 with {os.path.basename(grids_file_path)} ...")
        try:
            result = subprocess.run(
                ["autogrid4", "-p", os.path.basename(grids_file_path), "-l", os.path.basename(grid_log_file)],
                cwd=grids_path,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            if result.returncode == 0:
                print(f"✅ AutoGrid4 completed successfully. Log: {grid_log_file}")
            else:
                print(f"❌ AutoGrid4 failed. See log for details: {grid_log_file}")
                print(result.stderr)
        except FileNotFoundError:
            print("❌ AutoGrid4 executable not found. Please ensure it is installed and in your PATH.")
        except Exception as e:
            print(f"❌ Error running AutoGrid4: {e}") 

    def _create_receptor_register(self, 
                                pdb_id: int,                 
                                pdb_name: str, 
                                pdb_to_convert: str,
                                pdbqt_file: str,
                                configs: dict, 
                                notes: str = None
                                ) -> bool:
        """
        Create a receptor register entry in the database for this project.

        This method stores metadata about a prepared receptor for docking, including:
        - The originating PDB file and its database ID
        - The processed PDBQT file path
        - The docking box configuration (center and size)
        - Optional user notes
        - Timestamp of creation

        The register is stored in the SQLite database:
            project_path/docking/receptors/receptors.db
        in a table called 'receptor_registers'.

        Args:
            pdb_id (int): ID for the PDB file originating the .pdbqt file (from pdbs.db)
            pdb_name (str): Name of the processed PDB file originating the .pdbqt file
            pdb_to_convert (str): Path to the PDB file used for conversion
            pdbqt_file (str): Path to the generated .pdbqt file
            configs (dict): Dictionary containing box information for docking (e.g., {'center': {...}, 'size': {...}})
            notes (str, optional): Optional notes about this receptor preparation

        Returns:
            bool: True if the register was created successfully, False otherwise

        Example:
            self._create_receptor_register(
                pdb_id=1,
                pdb_name="1abc",
                pdb_to_convert="/path/to/1abc_checked.pdb",
                pdbqt_file="/path/to/1abc_checked.pdbqt",
                configs={'center': {'x': 10, 'y': 20, 'z': 30}, 'size': {'x': 20, 'y': 20, 'z': 20}},
                notes="Prepared for Vina docking"
            )
        """
        try:
            import sqlite3
            import json
            from datetime import datetime

            receptors_db_path = os.path.join(self.path, 'docking', 'receptors', 'receptors.db')

            # Ensure the directory exists
            os.makedirs(os.path.dirname(receptors_db_path), exist_ok=True)

            conn = sqlite3.connect(receptors_db_path)
            cursor = conn.cursor()

            # Create table if it doesn't exist
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS receptor_registers (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    pdb_id INTEGER NOT NULL,
                    pdb_name TEXT NOT NULL,
                    pdb_to_convert TEXT NOT NULL,
                    pdbqt_file TEXT,
                    configs TEXT,
                    notes TEXT,
                    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')

            # Check for identical register (excluding 'notes')
            cursor.execute('''
                SELECT pdb_id, pdb_name, pdb_to_convert, pdbqt_file, configs
                FROM receptor_registers
            ''')
            new_entry = (
                pdb_id,
                pdb_name,
                pdb_to_convert,
                pdbqt_file,
                json.dumps(configs, indent=2) if configs is not None else None
            )
            for row in cursor.fetchall():
                # Compare all fields except notes and created_date
                if row == new_entry:
                    print(f"⚠️  An identical receptor register already exists for PDB name: {pdb_name}.")
                    print(f"   Existing register details:")
                    print(f"   - pdb_id: {row[0]}")
                    print(f"   - pdb_name: {row[1]}")
                    print(f"   - pdb_to_convert: {row[2]}")
                    print(f"   - pdbqt_file: {row[3]}")
                    print(f"   - configs: {row[4]}")
                    print("   No new register will be created.")
                    conn.close()
                    return False

            # Prompt for notes if not provided
            if notes is None:
                try:
                    notes = input("Enter a note for this .pdbqt file (Enter: empty note): ").strip()
                    if not notes:
                        notes = ""
                except KeyboardInterrupt:
                    print("\n   ⚠️  Note entry cancelled. Using empty note.")
                    notes = ""

            # Insert register entry
            cursor.execute('''
                INSERT INTO receptor_registers (
                    pdb_id, pdb_name, pdb_to_convert, pdbqt_file, configs, notes
                ) VALUES (?, ?, ?, ?, ?, ?)
            ''', (
                pdb_id,
                pdb_name,
                pdb_to_convert,
                pdbqt_file,
                json.dumps(configs, indent=2) if configs is not None else None,
                notes
            ))

            conn.commit()
            conn.close()
            print(f"   ✅ Receptor register created for PDB ID: '{pdb_id}' (PDBQT: {os.path.basename(pdbqt_file)})")
            return True

        except Exception as e:
            print(f"   ❌ Error creating receptor register: {e}")
            return False
        
    def _check_ZN_presence(self, pdbqt_file: str, rec_main_path: str) -> bool:
        """
        Check if a Zinc (ZN) atom is present in the given PDBQT file.

        This method scans the PDBQT file for lines containing a ZN atom (in columns 13-16)
        and returns True if at least one is found, otherwise False.

        Args:
            pdbqt_file (str): Path to the PDBQT file to check.
            rec_main_path (str): Path to the receptor's main folder (unused, kept for compatibility).

        Returns:
            bool: True if a ZN atom is present, False otherwise.

        Raises:
            None. All exceptions are handled internally with user feedback.
        """
        try:
            with open(pdbqt_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')) and ' ZN ' in line[12:16]:
                        return True
            return False
        except Exception as e:
            print(f"   ❌ Error checking Zn presence: {e}")
            return False

    def _apply_ad4zn_params(self, pdbqt_file: str, rec_main_path: str) -> None:
        """
        Apply AutoDock4Zn parameters to a receptor PDBQT file by running the zinc_pseudo.py script.

        This method uses the zinc_pseudo.py script (from the TidyScreen package) to add the TZ pseudo atom
        required for proper AutoDock4Zn scoring. It runs the script in the 'adt' conda environment and
        replaces the original PDBQT file with the processed one containing the TZ atom.

        Args:
            pdbqt_file (str): Path to the receptor PDBQT file to process.
            rec_main_path (str): Path to the receptor's main folder (not used, kept for compatibility).

        Returns:
            None

        Raises:
            None. All exceptions are handled internally with user feedback.
        """
        import shutil
        import subprocess
        import site
        import os

        try:
            # Locate the zinc_pseudo.py script in the TidyScreen package
            zn_script_path = os.path.join(site.getsitepackages()[0], "tidyscreen", "misc", "zinc_pseudo.py")
            output_pdbqt = pdbqt_file.replace(".pdbqt", "_tz.pdbqt")

            # Run the script in the 'adt' conda environment
            result = subprocess.run(
                ["conda", "run", "-n", "adt", "python", zn_script_path, "-r", pdbqt_file, "-o", output_pdbqt],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Check if output file was created
            if os.path.exists(output_pdbqt):
                shutil.move(output_pdbqt, pdbqt_file)
                if result.returncode == 0:
                    print(f"✅ Successfully applied AutoDock4Zn parameters to {os.path.basename(pdbqt_file)}")
                else:
                    print(f"⚠️  AutoDock4Zn script returned an error, but output file was created.")
                    print(result.stderr)
            else:
                print(f"❌ Failed to create TZ-modified PDBQT file. See error below:")
                print(result.stderr)

        except Exception as e:
            print(f"   ❌ Error applying AD4Zn parameters: {e}")
        
    def _dock_table_with_autodockgpu(self, selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info):
        """
        Perform molecular docking of compounds using AutoDockGPU.
        This method iterates through compounds stored in a database table and performs
        docking simulations against a specified receptor using AutoDockGPU. The workflow
        includes:
        1. Retrieving molecular structures from the sdf_blob column
        2. Converting molecules to PDBQT format using Meeko
        3. Executing AutoDockGPU with configured parameters
        4. Optionally cleaning up temporary files
        Parameters
        ----------
        selected_table : str
            The name of the database table containing compounds to dock.
            Expected to have columns: inchi_key and sdf_blob.
        selected_method : dict
            Configuration dictionary containing:
            - 'parameters' (dict): AutoDockGPU parameters (e.g., exhaustiveness, number of modes)
            - 'ligand_prep_params' (dict): Meeko ligand preparation parameters
        assay_registry : dict
            Registry dictionary containing:
            - 'assay_folder_path' (str): Base directory for assay results and intermediate files
        clean_ligand_files : bool
            If True, removes temporary SDF and PDBQT files after successful docking.
            If False, retains all intermediate files for inspection.
        Returns
        -------
        None
            Results are written to disk in the assay results directory.
        Notes
        -----
        - Requires receptor PDBQT and grid field files (.fld) to be pre-processed
        - Displays progress using tqdm if available, otherwise prints status every 100 compounds
        - Errors in individual compound docking are caught and logged; processing continues
        - Output PDBQT files are named using inchi_key for compound identification
        Raises
        ------
        No exceptions are raised; errors are logged to console and processing continues.
        Warnings
        --------
        - If no receptor PDBQT file is found, the method returns early without docking
        - If no receptor grid field file is found, the method returns early without docking
        - If no compounds with SDF blobs are found, the method returns early
        """
        
        import sqlite3
        import tempfile
        import os
        import subprocess
        from meeko import MoleculePreparation
        from rdkit import Chem

        # Try to import tqdm for progress bar
        try:
            from tqdm import tqdm
            TQDM_AVAILABLE = True
        except ImportError:
            TQDM_AVAILABLE = False

        print(f"\n🚀 Starting docking with AutoDockGPU for table: {selected_table}")

        # Retrieve docking method parameters
        method_params = selected_method.get('parameters', {})
        
        # Retrieve ligand preparation params
        ligand_prep_params = selected_method.get('ligand_prep_params', {})
        
        # Select the receptor to be used for docking
        receptor_main_path = self.__receptor_path + f"/{receptor_info.get('pdb_name', None)}"
        receptor_pdbqt_file = receptor_main_path + f"/processed/receptor_checked.pdbqt"
        fld_file = receptor_main_path + f"/grid_files/receptor_checked.maps.fld"
        assay_folder = assay_registry['assay_folder_path']
    
        if not receptor_pdbqt_file:
            print("❌ No receptor .pdbqt file found.")
            return
        
        if not fld_file:
            print("❌ No receptor .maps.fld file found.")
            return
            
        # Connect to chemspace database and get compounds
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            cursor.execute(f"SELECT inchi_key, sdf_blob FROM {selected_table} WHERE sdf_blob IS NOT NULL")
            compounds = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error retrieving compounds from table: {e}")
            return

        if not compounds:
            print("❌ No compounds with SDF blobs found in the selected table.")
            return

        print(f"🧪 {len(compounds)} compounds to dock.")

        # Prepare output directory for docking results
        docking_results_dir = f"{assay_folder}/results"

        # Progress bar setup
        progress_bar = None
        if TQDM_AVAILABLE:
            progress_bar = tqdm(total=len(compounds), desc="Docking compounds", unit="ligand")

        # Iterate and dock each compound
        for idx, (inchi_key, sdf_blob) in enumerate(compounds, 1):
            try:
                # Write SDF blob to temp file
                # Use inchi_key as the temporary file name for the SDF
                sdf_filepath = f"{assay_folder}/results/{inchi_key}.sdf"
                
                with open(sdf_filepath, "wb") as sdf_temp:
                    sdf_temp.write(sdf_blob)
                
                # Create an RDKit mol object from a sdf file
                mol = self._mol_from_sdf(sdf_filepath)
                
                ## Generate the ligand.pdbqt file
                ligand_pdbqt_filepath = f"{docking_results_dir}/{inchi_key}.pdbqt"
                self._pdbqt_from_mol(mol, ligand_pdbqt_filepath, inchi_key, ligand_prep_params)
                
                ## Execute autodockgpu
                self._execute_autodockgpu(ligand_pdbqt_filepath, fld_file, method_params)
                
                if clean_ligand_files:
                    # Remove the temporary SDF and .pdbqt files
                    os.remove(sdf_filepath)
                    os.remove(ligand_pdbqt_filepath)

            except Exception as e:
                print(f"   ❌ Error docking compound {inchi_key}: {e}")
                continue

            if progress_bar:
                progress_bar.update(1)
            elif idx % 100 == 0:
                print(f"   Docked {idx} compounds...")

        if progress_bar:
            progress_bar.close()

        print(f"\n🎉 Docking completed for table: {selected_table}")
        print(f"Results saved in: {docking_results_dir}")
        
    def _select_receptor_from_db(self) -> Optional[dict]:
        """
        Prompt the user to select a receptor from the receptors database.

        This method loads the receptors database from the project's docking/receptors folder,
        displays a list of available receptors, and prompts the user to select one by number.
        It returns the value stored at the pdbqt_file column.
            
        """
        
        import sqlite3
        import os

        receptors_db_path = os.path.join(self.path, 'docking', 'receptors', 'receptors.db')
        if not os.path.exists(receptors_db_path):
            print(f"❌ No receptors database found at {receptors_db_path}")
            return None

        try:
            conn = sqlite3.connect(receptors_db_path)
            cursor = conn.cursor()
            cursor.execute('''
                SELECT id, pdb_id, pdb_name, pdbqt_file, notes, created_date
                FROM receptor_registers
                ORDER BY created_date ASC
            ''')
            receptors = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error loading receptors database: {e}")
            return None

        if not receptors:
            print("❌ No receptors found in the database.")
            return None

        print("\n📋 Available Receptors:")
        print("=" * 70)
        for i, (rid, pdb_id, pdb_name, pdbqt_file, notes, created_date) in enumerate(receptors, 1):
            print(f"{i}. {pdb_name} (ID: {pdb_id})")
            print(f"   PDBQT: {os.path.basename(pdbqt_file) if pdbqt_file else 'N/A'}")
            print(f"   Notes: {notes}")
            print(f"   Created: {created_date}")
        print("=" * 70)

        while True:
            try:
                selection = input("Select a receptor by number (or 'cancel'): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Receptor selection cancelled.")
                    return None
                try:
                    idx = int(selection) - 1
                    if 0 <= idx < len(receptors):
                        selected = receptors[idx]
                        # Return as a dictionary for convenience
                        return {
                            'id': selected[0],
                            'pdb_id': selected[1],
                            'pdb_name': selected[2],
                            'pdbqt_file': selected[3],
                            'notes': selected[4],
                            'created_date': selected[5]
                        }
                    else:
                        print(f"❌ Invalid selection. Enter a number between 1 and {len(receptors)}.")
                except ValueError:
                    print("❌ Please enter a valid number.")
            except KeyboardInterrupt:
                print("\n❌ Receptor selection cancelled.")
                return None
    
    def _get_receptor_conditions(self, receptor_id):
        """
        Retrieve receptor configuration and metadata from the receptors database.
        This method queries the receptors.db SQLite database located in the project's
        docking/receptors directory to fetch detailed information about a specific receptor
        based on its unique identifier.
        Args:
            receptor_id: The unique identifier (int or str) of the receptor to retrieve.
        Returns:
            dict or None: A dictionary containing receptor information with the following keys:
                - id (int): Unique receptor identifier
                - pdb_id (str): PDB identifier code
                - pdb_name (str): Human-readable receptor name
                - pdb_to_convert (str): Path to the PDB file requiring conversion
                - pdbqt_file (str): Path to the converted PDBQT file
                - configs (dict or None): Parsed JSON configuration object, or None if invalid/absent
                - notes (str): Additional notes or comments about the receptor
                - created_date (str): Timestamp of receptor record creation
            Returns None if:
                - The receptors database file does not exist
                - No receptor with the given ID is found
                - An error occurs during database query execution
        Raises:
            Catches and handles sqlite3 and json.JSONDecodeError exceptions internally,
            printing error messages to stdout.
        Side Effects:
            Prints status messages to stdout:
            - ❌ Error message if database not found
            - ❌ Error message if receptor not found
            - ⚠️  Warning if JSON config parsing fails
            - ❌ Error message if database query fails
        Note:
            The configs field contains JSON-serialized configuration data. Malformed
            JSON is logged as a warning but does not raise an exception.
        """ 
        
    
        import sqlite3
        import json
        
        receptors_db_path = os.path.join(self.path, 'docking', 'receptors', 'receptors.db')
        
        if not os.path.exists(receptors_db_path):
            print(f"❌ No receptors database found at {receptors_db_path}")
            
            return None
        
        try:
            conn = sqlite3.connect(receptors_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT id, pdb_id, pdb_name, pdb_to_convert, pdbqt_file, configs, notes, created_date
                FROM receptor_registers
                WHERE id = ?
            ''', (receptor_id,))
            
            receptor = cursor.fetchone()
            conn.close()
            
            if not receptor:
                print(f"❌ No receptor found with ID: {receptor_id}")
                return None
            
            # Parse configs JSON if present
            configs = None
            if receptor[5]:
                try:
                    configs = json.loads(receptor[5])
                except json.JSONDecodeError:
                    print(f"⚠️  Warning: Could not parse configs for receptor ID {receptor_id}")
                    configs = None
            
            # Return as dictionary
            return {
                'id': receptor[0],
                'pdb_id': receptor[1],
                'pdb_name': receptor[2],
                'pdb_to_convert': receptor[3],
                'pdbqt_file': receptor[4],
                'configs': configs,
                'notes': receptor[6],
                'created_date': receptor[7]
            }
            
        except Exception as e:
            print(f"❌ Error retrieving receptor conditions: {e}")
            return None
        
    def _mol_from_sdf(self, sdf_file: str):
            """
            Load a molecule from an SDF file using RDKit, restoring atom names if present.

            This method reads the first molecule from the given SDF file using RDKit's SDMolSupplier.
            If the molecule contains an 'AtomNames' property, atom names are restored and set as
            RDKit atom properties for compatibility with downstream tools (e.g., Meeko).

            Args:
                sdf_file (str): Path to the SDF file.

            Returns:
                mol (rdkit.Chem.Mol): RDKit molecule object with atom names restored, or None if loading fails.

            Raises:
                None. All exceptions are handled internally with user feedback.
            """
            from rdkit.Chem import SDMolSupplier

            try:
                supplier = SDMolSupplier(sdf_file, removeHs=False)
                mol = next((m for m in supplier if m is not None), None)
                if mol is None:
                    print(f"   ❌ No valid molecule found in SDF file: {sdf_file}")
                    return None

                # Restore atom names if present
                if mol.HasProp("AtomNames"):
                    atom_names = mol.GetProp("AtomNames").split("|")
                    for i, atom_name in enumerate(atom_names):
                        if i < mol.GetNumAtoms():
                            atom = mol.GetAtomWithIdx(i)
                            atom.SetProp("_Name", atom_name)
                            atom.SetProp("atomLabel", atom_name)

                return mol

            except Exception as e:
                print(f"   ❌ Error reading SDF file '{sdf_file}': {e}")
                return None
            
    def _pdbqt_from_mol(self, mol, ligand_pdbqt_file, inchi_key, ligand_prep_params):
        
        from meeko import MoleculePreparation
        from meeko import MoleculeSetup
        from meeko import PDBQTWriterLegacy
        
        
        try:
            # Compose MoleculePreparation arguments from ligand_prep_params
            mk_prep = MoleculePreparation(**ligand_prep_params)
            mol_setup_list = mk_prep.prepare(mol, rename_atoms=True)
            mol_setup = mol_setup_list[0]

            pdbqt_string = PDBQTWriterLegacy.write_string(mol_setup)

            with open(ligand_pdbqt_file, "w") as f:
                f.write(pdbqt_string[0])

        except Exception as e:
            print(f"   ❌ Error creating PDBQT for molecule {inchi_key}: {e}")
            
    def _execute_autodockgpu(self, ligand_pdbqt_filepath, fld_file, method_params):
        import subprocess
        import os

        fld_file_path = os.path.dirname(fld_file)
        fld_filename = fld_file.split('/')[-1]

        # Build AutoDockGPU command
        cmd = [
            "autodock_gpu_128wi",
            "-ffile", fld_filename,
            "-lfile", ligand_pdbqt_filepath,
            "--heuristics", str(method_params.get('heuristics', 1)),
            "--heurmax", str(method_params.get('heurmax', 12000000)),
            "--autostop", str(method_params.get('autostop', 1)),
            "--asfreq", str(method_params.get('asfreq', 5)),
            "--nrun", str(method_params.get('nrun', 10)),
            "--nev", str(method_params.get('nev', 2500000)),
            "--ngen", str(method_params.get('gen', 42000)),
            "--lsmet", str(method_params.get('lsmet', "ad")),
            "--lsit", str(method_params.get('lsit', 300)),
            "--psize", str(method_params.get('psize', 150)),
            "--mrat", str(method_params.get('mrat', 2)),
            "--crat", str(method_params.get('crat', 80)),
            "--lsrat", str(method_params.get('lsrat', 100)),
            "--trat", str(method_params.get('trat', 60)),
            "--dmov", str(method_params.get('dmov', 6)),
            "--dang", str(method_params.get('dang', 90)),
            "--rholb", str(method_params.get('rholb', 0.01)),
            "--lsmov", str(method_params.get('lsmov', 2)),
            "--lsang", str(method_params.get('lsang', 75)),
            "--cslim", str(method_params.get('cslim', 4)),
            "--stopstd", str(method_params.get('stopstd', 0.15)),
            "--contact_analysis", str(method_params.get('contact_analysis', 1)),
            "--xmloutput", str(method_params.get('xmloutput', 0)),
            "--clustering", str(method_params.get('clustering', 1)),
            "--output-cluster-poses", str(method_params.get('output-cluster-poses', 0)),
            "--hsym", str(method_params.get('hsym', 1)),
            "--rmstol", str(method_params.get('rmstol', 2)),
        ]
        
        # Change working directory to docking_results_dir before running the command
        try:
            result = subprocess.run(cmd, cwd=fld_file_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                #print(f"   ✅ Docking completed for {ligand_pdbqt_filepath.split("/")[-1]}")
                pass
            else:
                print(f"   ❌ Docking failed for {ligand_pdbqt_filepath.split("/")[-1]}")
                print(result.stderr)
        except Exception as e:
            print(f"   ❌ Error running AutoDockGPU for {ligand_pdbqt_filepath.split("/")[-1]}: {e}")

    def _process_docking_assay_results(self, assay_registry, docking_mode, clean_ligand_files, max_poses):

        from ringtail import RingtailCore
        import os

        try:
            assay_folder = assay_registry['assay_folder_path']
            assay_id = assay_registry['assay_id']
            results_folder = f"{assay_folder}/results"
            results_db_file = f"{assay_folder}/results/assay_{assay_id}.db"
            rtc = RingtailCore(db_file = results_db_file, docking_mode=docking_mode )
            rtc.add_results_from_files(file_path = f"{assay_folder}/results", recursive = True, save_receptor = False, max_poses=max_poses)
            
            
            ## Add ligand input model if docking with docking output is .pdqbt type
            if docking_mode == 'vina':
                
                import glob
                
                # Get all .out.pdbqt files from the results directory
                out_pdbqt_files = glob.glob(os.path.join(assay_folder, 'results', '*_out.pdbqt'))

                # Add ligand input models to the ringtail results
                for out_file in out_pdbqt_files:
                    input_model = self._parse_pdbqt_output(out_file)
                    self._add_input_model_in_db(results_db_file, input_model, out_file)
                    
            ## Clean the .dlg files is requested
            
            if clean_ligand_files:
                if os.path.exists(results_db_file) and os.path.getsize(results_db_file) > 0:
                    for fname in os.listdir(results_folder):
                        if fname.endswith('.dlg'):
                            try:
                                os.remove(os.path.join(results_folder, fname))
                            except Exception as e:
                                print(f"Warning: could not delete {fname}: {e}")

            return results_db_file

        except Exception as error:
            print(error)
            print("Error processing with Ringtail")

    def list_pdb_templates(self) -> None:
        """
        Reads the 'pdbs.db' database located in the project_path/docking/receptors folder
        and prints out the following columns from the 'pdb_templates' table:
        - 'pdb_name'
        - 'pdb_analysis'
        - 'processed_pdb_path'
        - 'notes'
        Each field in the rows is printed in a new line.
        """
        import sqlite3
        pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
        if not os.path.exists(pdbs_db_path):
            print(f"❌ Database not found: {pdbs_db_path}")
            return

        try:
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            cursor.execute('''
                SELECT pdb_id, pdb_name, processed_pdb_path, notes
                FROM pdb_templates
                ORDER BY created_date ASC
            ''')
            rows = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error reading pdbs.db: {e}")
            return

        if not rows:
            print("📋 No PDB templates found in the database.")
            return

        print("\n📋 PDB Templates:")
        print("=" * 40)
        for pdb_id, pdb_name, processed_pdb_path, notes in rows:
            print(f"pdb_id: {pdb_id}")
            print(f"Name: {pdb_name}")
            print(f"Processed PDB Path: {processed_pdb_path}")
            print(f"Notes: {notes}")
            print("-" * 40)

    def remove_pdb_template(self):

        """
        Will list entried in the pdb_templates table of the pdbs.db database located in the project_path/docking/receptors folder, and query the user which one to delete. The registry is finally deleted from the database.
        """

        import sqlite3
        pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
        if not os.path.exists(pdbs_db_path):
            print(f"❌ Database not found: {pdbs_db_path}")
            return

        try:
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            cursor.execute('''
                SELECT pdb_id, pdb_name
                FROM pdb_templates
                ORDER BY created_date ASC
            ''')
            rows = cursor.fetchall()
        except Exception as e:
            print(f"❌ Error reading pdbs.db: {e}")
            return

        if not rows:
            print("📋 No PDB templates found in the database.")
            return

        print("\n📋 PDB Templates:")
        print("=" * 40)
        for i, (pdb_id, pdb_name) in enumerate(rows, 1):
            print(f"{i}. ID: {pdb_id}, Name: {pdb_name}")
        print("=" * 40)

        while True:
            try:
                selection = input("Select a PDB template to delete by number (or 'cancel'): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ PDB template deletion cancelled.")
                    return
                try:
                    idx = int(selection) - 1
                    if 0 <= idx < len(rows):
                        selected_id = rows[idx][0]
                        cursor.execute('DELETE FROM pdb_templates WHERE pdb_id = ?', (selected_id,))
                        conn.commit()
                        conn.close()
                        print(f"✅ Deleted PDB template with ID: {selected_id}")
                        return
                    else:
                        print(f"❌ Invalid selection. Enter a number between 1 and {len(rows)}.")
                except ValueError:
                    print("❌ Please enter a valid number.")
            except KeyboardInterrupt:
                print("\n❌ PDB template deletion cancelled.")
                return

    def list_receptor_models(self):
        """
        Reads the 'receptors.db' database located in the project_path/docking/receptors folder
        and prints out the following columns from the 'receptor_registers' table:
        - 'pdb_name'
        - 'pdb_to_convert'
        - 'pdbqt_file'
        - 'configs'
        - 'notes'
        Each field in the rows is printed in a new line.
        """
        import sqlite3
        receptors_db_path = os.path.join(self.path, 'docking', 'receptors', 'receptors.db')
        if not os.path.exists(receptors_db_path):
            print(f"❌ Database not found: {receptors_db_path}")
            return

        try:
            conn = sqlite3.connect(receptors_db_path)
            cursor = conn.cursor()
            cursor.execute('''
                SELECT pdb_name, pdb_to_convert, pdbqt_file, configs, notes
                FROM receptor_registers
                ORDER BY created_date ASC
            ''')
            rows = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error reading receptors.db: {e}")
            return

        if not rows:
            print("📋 No receptor models found in the database.")
            return

        print("\n📋 Receptor Models:")
        print("=" * 40)
        for pdb_name, pdb_to_convert, pdbqt_file, configs, notes in rows:
            print(f"PDB Name: {pdb_name}")
            print(f"PDB to Convert: {pdb_to_convert}")
            print(f"PDBQT File: {pdbqt_file}")
            print(f"Configs: {configs}")
            print(f"Notes: {notes}")
            print("-" * 40)
        
    def _get_ligand_preparation_parameters(self):
        """
        Get ligand preparation parameters to be used for Meeko procedures
        
        """    
        
        try:
            print(f"\n⚙️  CONFIGURE LIGAND PREPARATION PARAMETERS (MEEKO)")
            print("-" * 50)
            
            ligand_prep_parameters = {}
            
            ## Meeko specific parameters
            
            ligand_prep_parameters['merge_these_atom_types'] = self._input_parameter_choice(
                "Atom types to merge. Default: '(H,)'",
                default="(H,)"
            )
            
            ligand_prep_parameters['hydrate'] = self._get_parameter_choice(
                "Add water molecules to the structure for hydrated docking.",
                [False, True],
                default=False
            )
            
            ligand_prep_parameters['add_atom_types'] = self._input_parameter_choice(
                "Specify additional atom types to assign in JSON format, with SMARTS patterns and atom type names. (i.e.: '[{'smarts': '[#1][#6X3]:[#6X3]([#6])[#7]:[#7][#7][#6]', 'atype': 'HD'}]' for 1,2,3-triazole donor HD in 4 position. Default: None)",
                default=None
            )
            
            ligand_prep_parameters['charge'] = self._get_parameter_choice(
                "Choose the charge model.",
                ["gasteiger", "espaloma", "zero"],
                default="gasteiger"
            )
            
            return ligand_prep_parameters if ligand_prep_parameters else None
            
        except KeyboardInterrupt:
            print("\n❌ Ligands parameter configuration cancelled")
            return None
        except Exception as e:
            print(f"❌ Error configuring ligands parameters: {e}")
            return None
        pass
    
    def extract_docked_poses(self, assay=None):
        """
        Extract docked poses from a docking assay and save them as PDB files based on user-selected criteria.
        This method allows users to select a docking assay from the registry and extract poses based on three
        different selection criteria: lowest energy poses, most populated poses, or a combination of both.
        The selected poses are processed and saved as individual PDB files in an output directory.
        Parameters
        ----------
        assay : int, str, or None, optional
            The docking assay identifier to process. Can be:
            - None (default): User is prompted to select from available assays
            - int: Assay ID to match against the registry
            - str: Assay name to match against the registry
        Returns
        -------
        None
            The method does not return a value. Output is written to PDB files in the designated
            output directory based on the selection criteria.
        Raises
        ------
        None
            Errors are caught internally and reported via print statements. Returns None on failure.
        Notes
        -----
        - Requires a valid docking assay registry at 'docking/docking_registers/docking_assays.db'
        - Selection criteria options:
            1. Lowest energy poses: Selects poses with minimum binding energy
            2. Most populated poses: Selects poses with highest population/frequency
            3. Both criteria: Selects poses meeting both lowest energy and highest population thresholds
        - Output directories are created automatically in the results directory with descriptive names
        - Each extracted pose is saved as a separate PDB file named with ligand name and run number
        Examples
        --------
        >>> moldock_instance.extract_docked_poses()  # Interactive mode - prompts user selection
        >>> moldock_instance.extract_docked_poses(assay=1)  # Extract from assay with ID 1
        >>> moldock_instance.extract_docked_poses(assay="kinase_screen")  # Extract from named assay
        """
        
        import os
        import sqlite3

        # Path to docking assay registry
        registry_db = os.path.join(self.path, 'docking', 'docking_registers', 'docking_assays.db')
        if not os.path.exists(registry_db):
            print(f"❌ Docking assay registry not found at {registry_db}")
            return None

        # Connect to registry and fetch assays
        try:
            conn = sqlite3.connect(registry_db)
            cursor = conn.cursor()
            cursor.execute("SELECT assay_id, assay_name, docking_engine, assay_folder_path, status FROM docking_assays")
            assays = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error reading docking assay registry: {e}")
            return None

        if not assays:
            print("❌ No docking assays found in the registry.")
            return None

        # If assay is not provided, prompt user to select one
        selected_assay = None
        if assay is None:
            print("\n📋 Available Docking Assays:")
            print("=" * 60)
            for i, (aid, name, engine, folder, status) in enumerate(assays, 1):
                print(f"{i}. {name} (ID: {aid}, Engine: {engine}, Status: {status})\n   Folder: {folder}")
            print("=" * 60)
            while True:
                try:
                    choice = input("Select an assay by number (or press Enter to cancel): ").strip()
                    if not choice:
                        print("❌ Operation cancelled.")
                        return None
                    idx = int(choice) - 1
                    if 0 <= idx < len(assays):
                        selected_assay = assays[idx]
                        break
                    else:
                        print("Invalid selection. Try again.")
                except (ValueError, KeyboardInterrupt):
                    print("❌ Operation cancelled.")
                    return None
        else:
            # If assay is provided as ID or name, find it
            for a in assays:
                if (isinstance(assay, int) and a[0] == assay) or (isinstance(assay, str) and a[1] == assay):
                    selected_assay = a
                    break
            if not selected_assay:
                print(f"❌ Assay '{assay}' not found in registry.")
                return None

        assay_id, assay_name, docking_engine, assay_folder, status = selected_assay

        # Locate results directory
        results_dir = os.path.join(assay_folder, "results")
        if not os.path.exists(results_dir):
            print(f"❌ Results directory not found: {results_dir}")
            return None

        # Create the path were the database to be analyzed lives
        db_path = os.path.join(results_dir, f"assay_{assay_id}.db")

        ## Evaluate the docking engine in order to use the corresponding pose extraction method
        if docking_engine == "AutoDockGPU":
            self._extract_docked_poses_autodockgpu(db_path)
            
        elif docking_engine == "Vina":
            self._extract_docked_poses_autodockvina(db_path)

    def _extract_docked_poses_autodockgpu(self, db_path):
    
        ## Start the analysis
        print("Starting poses extraction for engine AutoDockGPU")
        print("")
        print("Select an option:")
        print("1 - Select lowest energy poses")
        print("2 - Select most populated poses")
        print("3 - Select most lowest energy AND populated poses (both criteria)")
        selection = input("Enter your choice (1, 2, or 3): ").strip()

        ## Prompt the user which selection to make
        while selection not in ("1", "2", "3"):
            print("Invalid selection. Only 1, 2, or 3 can be selected.")
            selection = input("Enter your choice (1, 2, or 3): ").strip()

        # Act according to user selection
        if selection == "1":
            print("You selected: lowest energy poses")
            # Create output directory for most stable poses
            output_dir = os.path.join(os.path.dirname(os.path.abspath(db_path)), "most_stable_poses")
            os.makedirs(output_dir, exist_ok=True)
            
            df = self._select_most_stable_poses(db_path)
        
        elif selection == "2":
            print("You selected: most populated poses")
            # Create output directory for most populated poses
            output_dir = os.path.join(os.path.dirname(os.path.abspath(db_path)), "most_populated_poses")
            os.makedirs(output_dir, exist_ok=True)
            df = self._select_most_populated_poses(db_path)

        elif selection == "3":
            print("You selected: most populated AND lowest energy poses")
            # Create output directory for most populated and stable poses
            output_dir = os.path.join(os.path.dirname(os.path.abspath(db_path)), "most_populated_and_stable_poses")
            os.makedirs(output_dir, exist_ok=True)
            df = self._select_most_populated_and_stable_poses(db_path)

        # Continue pprocessing the selected poses
        for row in df.iterrows():
            ligname = row[1]['LigName']
            pose_id = row[1]['Pose_ID']
            run_number = row[1]['run_number']
            
            input_model, pose_coords_json = self._retrieve_pose_info(db_path, pose_id)
            print(type(input_model))
            pdb_dict = self._process_input_model_for_moldf(input_model, pose_coords_json)
            self._write_pdb_with_moldf(pdb_dict, ligname, run_number, output_dir)

        if selection == "1":
            print(f"Most stable poses extracted and saved as PDB files to: \n \t {output_dir}")
        elif selection == "2":
            print(f"Most populated poses extracted and saved as PDB files to: \n \t {output_dir}")
        elif selection == "3":
            print(f"Most populated and stable poses extracted and saved as PDB files to: \n \t {output_dir}")    
    
    def _extract_docked_poses_autodockvina(self, db_path):
        
        ## Start the analysis
        print("Starting poses extraction for engine Vina")
        print("")
        print("Select an option:")
        print("1 - Extract lowest energy pose")
        print("2 - Extract all poses")
        selection = input("Enter your choice (1 or 2): ").strip()
        
        ## Prompt the user which selection to make
        while selection not in ("1", "2"):
            print("Invalid selection. Only 1 or 2 can be selected.")
            selection = input("Enter your choice (1 or 2): ").strip()
        
        # Act according to user selection
        if selection == "1":
            print("You selected: lowest energy poses")
            # Create output directory for most stable poses
            output_dir = os.path.join(os.path.dirname(os.path.abspath(db_path)), "most_stable_poses")
            os.makedirs(output_dir, exist_ok=True)
            
            df = self._select_most_stable_poses(db_path)

        elif selection == "2":
            print("You selected: all poses")
            # Create output directory for all poses
            output_dir = os.path.join(os.path.dirname(os.path.abspath(db_path)), "all_poses")
            os.makedirs(output_dir, exist_ok=True)
            
            df = self._select_all_poses(db_path)

        # Continue pprocessing the selected poses
        for row in df.iterrows():
            ligname = row[1]['LigName']
            pose_id = row[1]['Pose_ID']
            run_number = row[1]['run_number']
            
            input_model, pose_coords_json = self._retrieve_pose_info(db_path, pose_id)
            pdb_dict = self._process_input_model_for_moldf(input_model, pose_coords_json)
            self._write_pdb_with_moldf(pdb_dict, ligname, run_number, output_dir)

        print(f"Docked poses extracted to {output_dir}")

    def _select_most_stable_poses(self, db_path):
        """
        Select the most stable pose (lowest docking score) for each ligand
        from the Ringtail database located at db_path.
        Returns a DataFrame with LigName, Pose_ID, and docking_score.
        """
        import sqlite3
        import pandas as pd
        
        conn = sqlite3.connect(db_path)
        query = """
        SELECT LigName, Pose_ID, docking_score, run_number
        FROM Results AS R1
        WHERE docking_score = (
            SELECT MIN(docking_score)
            FROM Results AS R2
            WHERE R1.LigName = R2.LigName
        )
        """
        df = pd.read_sql(query, conn)
        conn.close()
        
        return df

    def _select_all_poses(self, db_path):
        """
        Select all the poses for each ligand
        from the Ringtail database located at db_path.
        Returns a DataFrame with LigName, Pose_ID, and docking_score.
        """
        import sqlite3
        import pandas as pd
        
        conn = sqlite3.connect(db_path)
        query = """
        SELECT LigName, Pose_ID, docking_score, run_number
        FROM Results 
        """
        df = pd.read_sql(query, conn)
        conn.close()
        
        return df

    def _select_most_populated_poses(self, db_path):
        """
        Select the most populated pose (highest value in the 'cluster_size' column) for each ligand
        from the Ringtail database located at db_path.
        Returns a DataFrame with LigName, Pose_ID, and docking_score.
        """

        import sqlite3
        import pandas as pd

        conn = sqlite3.connect(db_path)
        query = """
        SELECT LigName, Pose_ID, cluster_size, run_number
        FROM Results AS R1
        WHERE cluster_size = (
            SELECT MAX(cluster_size)
            FROM Results AS R2
            WHERE R1.LigName = R2.LigName
        )
        """
        df = pd.read_sql(query, conn)
        conn.close()
        
        return df

    def _select_most_populated_and_stable_poses(self, db_path):
        """
        Selects the most populated and stable docking poses for each ligand from a SQLite database.
        This method connects to the specified SQLite database, queries the 'Results' table,
        and retrieves the pose(s) for each ligand that have the largest cluster size (most populated)
        and the lowest docking score (most stable). The results are returned as a pandas DataFrame.
        Args:
            db_path (str): Path to the SQLite database file containing the docking results.
        Returns:
            pandas.DataFrame: A DataFrame containing the ligand name, pose ID, cluster size, docking score, and run number for the selected poses.
        """

        import sqlite3
        import pandas as pd

        conn = sqlite3.connect(db_path)
        query = """
        SELECT LigName, Pose_ID, cluster_size, docking_score, run_number
        FROM Results AS R1
        WHERE cluster_size = (
            SELECT MAX(cluster_size)
            FROM Results AS R2
            WHERE R1.LigName = R2.LigName
        )
        AND docking_score = (
            SELECT MIN(docking_score)
            FROM Results AS R3
            WHERE R1.LigName = R3.LigName
        )
        """
        df = pd.read_sql(query, conn)
        conn.close()    
        
        return df

    def _retrieve_pose_info(self, db_path, pose_id):
        """
        Retrieves ligand pose information from the database for a given pose ID.
        Connects to the specified SQLite database, queries the Ligands and Results tables
        to obtain the input model and ligand coordinates associated with the provided pose ID.
        Args:
            db_path (str): Path to the SQLite database file.
            pose_id (int or str): Identifier of the pose to retrieve.
        Returns:
            tuple: A tuple containing:
                - input_model (str): The input model data for the ligand.
                - pose_coords_json (str): JSON string of the ligand coordinates.
        Raises:
            sqlite3.DatabaseError: If there is an error connecting to or querying the database.
            TypeError: If the query does not return any results.
        """
       
        import sqlite3

        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Retrieve the info required to reconstruct the ligand pose
        cursor.execute("""
        SELECT L.input_model, R.ligand_coordinates FROM Ligands L
        JOIN Results R ON L.LigName = R.LigName
        WHERE R.Pose_ID = ?
        """, (pose_id,))
        row = cursor.fetchone()
        
        input_model, pose_coords_json = row
        
        return input_model, pose_coords_json

    def _process_input_model_for_moldf(self, input_model, coords_json):
        """
        Processes an input molecular model and coordinates to generate a DataFrame suitable for moldf.
        This method parses the input model (as a string representation of a list), extracts atom information,
        replaces atomic coordinates with those provided in `coords_json`, and formats the data into a DataFrame
        with standardized columns for molecular data. The resulting DataFrame is returned in a dictionary under
        the key '_atom_site', suitable for use with moldf.
        Args:
            input_model (str): String representation of a list containing atom records. Each record can be a string
                or a list/tuple, where the first element should be 'ATOM'.
            coords_json (str): JSON-encoded list of new coordinates (list of [x, y, z]) to replace the original
                coordinates in the input model.
        Returns:
            dict: A dictionary with a single key '_atom_site', whose value is a pandas DataFrame containing the
                processed atom data with updated coordinates and standardized columns.
        """

        import ast
        import pandas as pd
        import json


        input_model_list = ast.literal_eval(input_model)
        atom_rows = []
        for item in input_model_list:
            # If item is a string, split it and check if first element is 'ATOM'
            if isinstance(item, str):
                parts = item.split()
                if len(parts) > 0 and parts[0] == 'ATOM':
                    atom_rows.append(parts)
            # If item is a list/tuple, check if first element is 'ATOM'
            elif isinstance(item, (list, tuple)) and len(item) > 0 and item[0] == 'ATOM':
                atom_rows.append(item)
        
        
        df = pd.DataFrame(atom_rows)
        
        # Rename columns using a list of column names
        column_names = [
            'record_name', 'atom_number', 'atom_name', 'residue_name',
            'residue_number', 'x_coord', 'y_coord', 'z_coord',
            'occupancy', 'b_factor', 'charge', 'element_symbol'
        ]
        df.columns = column_names[:len(df.columns)]
        
        # Insert an empty column 'alt_loc' at the fourth position (index 3)
        df.insert(3, 'alt_loc', '')
        df.insert(5, 'chain_id', '')
        df.insert(7, 'insertion', '')
        df.insert(13, 'segment_id', '')

        # Move 'charge' column to the last position
        if 'charge' in df.columns:
            charge_col = df.pop('charge')
            df['charge'] = charge_col

        # Replace 'charge' column with empty items of type object
        if 'charge' in df.columns:
            df['charge'] = pd.Series([''] * len(df), dtype=object)

        ## Replace the coordinates with the ones from coords_json for the given pose
        coords_list = json.loads(coords_json)

        for i, coord in enumerate(coords_list):
            # Replace the original coordinates with the new ones corresponding to the pose
            df.at[i, 'x_coord'] = coord[0]
            df.at[i, 'y_coord'] = coord[1]
            df.at[i, 'z_coord'] = coord[2]  

            # Attempt to convert columns to appropriate types
        for col in df.columns:
            # Try to convert to int, then float, else keep as string
            try:
                df[col] = df[col].astype(int)
            except ValueError:
                try:
                    df[col] = df[col].astype(float)
                except ValueError:
                    pass

        # Return as dictionary for moldf
        return {'_atom_site': df}

    def _write_pdb_with_moldf(self,pdb_dict, ligname, run_number, output_dir):
        """
        Writes a PDB file using the provided molecular data dictionary and saves it to the specified output directory.

        Args:
            pdb_dict (dict): Dictionary containing the molecular structure data to be written to the PDB file.
            ligname (str): Name of the ligand, used as part of the output filename.
            run_number (int or str): Identifier for the run, used as part of the output filename.
            output_dir (str): Directory path where the output PDB file will be saved.

        Returns:
            None

        Side Effects:
            Creates a PDB file in the specified output directory with a filename formatted as "{ligname}_{run_number}.pdb".
        """
        from moldf import write_pdb
        import os

        output_file = os.path.join(output_dir, f"{ligname}_{run_number}.pdb")
        write_pdb(pdb_dict, output_file)

    def _dock_table_with_vina(self, selected_table, selected_method, assay_registry, clean_ligand_files, receptor_info):
        
        """
        Perform molecular docking of compounds using Vina.
        This method iterates through compounds stored in a database table and performs
        docking simulations against a specified receptor using Vina. The workflow
        includes:
        1. Retrieving molecular structures from the sdf_blob column
        2. Converting molecules to PDBQT format using Meeko

        """

        import sqlite3
        import tempfile
        import os
        import subprocess
        from meeko import MoleculePreparation
        from rdkit import Chem

        # Try to import tqdm for progress bar
        try:
            from tqdm import tqdm
            TQDM_AVAILABLE = True
        except ImportError:
            TQDM_AVAILABLE = False

        print(f"\n🚀 Starting docking with AutoDockGPU for table: {selected_table}")

        # Retrieve docking method parameters
        method_params = selected_method.get('parameters', {})
        
        # Retrieve ligand preparation params
        ligand_prep_params = selected_method.get('ligand_prep_params', {})
        
        # Select the receptor to be used for docking
        #selected_receptor = self._select_receptor_from_db()
        
        receptor_conditions = self._get_receptor_conditions(receptor_info.get('id', None))
        receptor_main_path = self.__receptor_path + f"/{receptor_info.get('pdb_name', None)}"
        receptor_pdbqt_file = receptor_main_path + f"/processed/receptor_checked.pdbqt"
        fld_file = receptor_main_path + f"/grid_files/receptor_checked.maps.fld"
        assay_folder = assay_registry['assay_folder_path']
    
        if not receptor_pdbqt_file:
            print("❌ No receptor .pdbqt file found.")
            return
        
        if not fld_file:
            print("❌ No receptor .maps.fld file found.")
            return
            
        # Connect to chemspace database and get compounds
        try:
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            cursor.execute(f"SELECT inchi_key, sdf_blob FROM {selected_table} WHERE sdf_blob IS NOT NULL")
            compounds = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error retrieving compounds from table: {e}")
            return

        if not compounds:
            print("❌ No compounds with SDF blobs found in the selected table.")
            return

        print(f"🧪 {len(compounds)} compounds to dock.")

        # Prepare output directory for docking results
        docking_results_dir = f"{assay_folder}/results"

        # Progress bar setup
        progress_bar = None
        if TQDM_AVAILABLE:
            progress_bar = tqdm(total=len(compounds), desc="Docking compounds", unit="ligand")

        ## If 'vina' scoring function is to be applied, compute the maps once for batch docking
        if method_params['sf_name'] == 'vina':          
            vina_rec_object = self._prepare_vina_receptor(method_params, receptor_pdbqt_file, receptor_conditions)

        # Iterate and dock each compound
        for idx, (inchi_key, sdf_blob) in enumerate(compounds, 1):
            try:
                # Write SDF blob to temp file
                # Use inchi_key as the temporary file name for the SDF
                sdf_filepath = f"{assay_folder}/results/{inchi_key}.sdf"
                
                with open(sdf_filepath, "wb") as sdf_temp:
                    sdf_temp.write(sdf_blob)
                
                # Create an RDKit mol object from a sdf file
                mol = self._mol_from_sdf(sdf_filepath)
                
                ## Generate the ligand.pdbqt file
                ligand_pdbqt_filepath = f"{docking_results_dir}/{inchi_key}.pdbqt"
                self._pdbqt_from_mol(mol, ligand_pdbqt_filepath, inchi_key, ligand_prep_params)

                ## Execute vina
                self._execute_vina(ligand_pdbqt_filepath, assay_registry, method_params, receptor_pdbqt_file, receptor_conditions, vina_rec_object=None)
                
                if clean_ligand_files:
                    # Remove the temporary SDF and .pdbqt files
                    os.remove(sdf_filepath)
                    os.remove(ligand_pdbqt_filepath)

            except Exception as e:
                print(f"   ❌ Error docking compound {inchi_key}: {e}")
                continue

            if progress_bar:
                progress_bar.update(1)
            elif idx % 100 == 0:
                print(f"   Docked {idx} compounds...")

        if progress_bar:
            progress_bar.close()

        print(f"\n🎉 Docking completed for table: {selected_table}")
        print(f"Results saved in: {docking_results_dir}")
    
    def _prepare_vina_receptor(self, method_params, receptor_pdbqt_file, receptor_conditions):
        
        from vina import Vina
        
        vina_rec_object = Vina(sf_name=method_params['sf_name'], verbosity=0)
        
        vina_rec_object.set_receptor(receptor_pdbqt_file)
        
        box_center = [receptor_conditions['configs']['center']['x'], receptor_conditions['configs']['center']['y'], receptor_conditions['configs']['center']['z']]
        
        box_size = [receptor_conditions['configs']['size']['x'], receptor_conditions['configs']['size']['y'], receptor_conditions['configs']['size']['z']]

        vina_rec_object.compute_vina_maps(center=box_center, box_size=box_size)        
        
        return vina_rec_object
    
    def _execute_vina(self, ligand_pdbqt_filepath, assay_registry, method_params, receptor_pdbqt_file, receptor_conditions, vina_rec_object):
        """
        Execute molecular docking using AutoDock Vina with specified scoring function.
        This method performs batch docking of a ligand against a receptor using either
        the 'vina' or 'ad4' (AutoDock 4) scoring functions. It handles ligand preparation,
        docking execution, and output file generation.
        Parameters
        ----------
        ligand_pdbqt_filepath : str
            Full file path to the ligand PDBQT file to be docked.
        assay_registry : dict or object
            Registry containing assay-related information and metadata.
        method_params : dict
            Dictionary containing docking method parameters:
            - 'sf_name' : str
                Scoring function name ('vina' or 'ad4').
            - 'exhaustiveness' : int
                Exhaustiveness parameter for docking search space exploration.
            - 'n_poses' : int
                Number of output poses to generate.
        receptor_pdbqt_file : str
            File path to the receptor PDBQT file.
        receptor_conditions : dict
            Dictionary containing receptor-related conditions:
            - 'pdbqt_file' : str
                Path to the processed receptor PDBQT file.
        vina_rec_object : Vina or None
            Pre-initialized Vina object for 'vina' scoring function, or None for 'ad4'.
        Returns
        -------
        None
            Results are written to output PDBQT files (ligand_*_out.pdbqt).
        Notes
        -----
        - For 'vina' scoring function: Uses pre-loaded receptor maps from vina_rec_object.
        - For 'ad4' scoring function: Initializes new Vina object and loads affinity maps
          from the grid_files directory.
        - Output files are generated with '_out.pdbqt' suffix, replacing '.pdbqt' in
          the original ligand filename.
        - Existing output files are overwritten if they exist.
        """
        

        from vina import Vina
        
        ## This will apply batch docking using the 'vina' scoring function
        if method_params['sf_name'] == 'vina' and vina_rec_object is not None:

            vina_rec_object.set_ligand_from_file(ligand_pdbqt_filepath)

            vina_rec_object.dock(exhaustiveness=method_params['exhaustiveness'], n_poses=method_params['n_poses'])

            results_file = ligand_pdbqt_filepath.replace('.pdbqt', '_out.pdbqt')

            vina_rec_object.write_poses(results_file, n_poses=method_params['n_poses'], overwrite=True)

    
        # This will apply docking using the 'ad4' scoring function
        elif method_params['sf_name'] == 'ad4':
            
            vina_rec_object = Vina(method_params['sf_name'])
            
            # Determine the grids files path using the receptor file as reference and load affinity maps for docking
            grid_maps_path = receptor_conditions['pdbqt_file'].replace('processed/receptor_checked.pdbqt', 'grid_files/receptor_checked')
            vina_rec_object.load_maps(grid_maps_path)
            
            vina_rec_object.set_ligand_from_file(ligand_pdbqt_filepath)
            
            vina_rec_object.dock(exhaustiveness=method_params['exhaustiveness'], n_poses=method_params['n_poses'])
            
            results_file = ligand_pdbqt_filepath.replace('.pdbqt', '_out.pdbqt')

            vina_rec_object.write_poses(results_file, n_poses=method_params['n_poses'], overwrite=True)
            
    def compute_roc_curve(self):
        
        """
        Prompt the user to select docking assays for binders (positives) and decoys (negatives)
        using list_assay_folders. Loads the best docking score per ligand from each assay’s
        Ringtail DB, computes an ROC curve/AUC, and writes ROC data to docking/analysis.
        """
        import os
        import sqlite3
        from datetime import datetime

        assays = self.list_assay_folders()
        if not assays:
            return None

        assay_by_id = {a['assay_id']: a for a in assays}

        print("\nSelect binder (positive) assay:")
        binder_id = self._prompt_assay("binders", assay_by_id)
        if binder_id is None:
            return None

        print("\nSelect decoy (negative) assay:")
        decoy_id = self._prompt_assay("decoys", assay_by_id)
        if decoy_id is None:
            return None

        binder_info = assay_by_id[binder_id]
        decoy_info = assay_by_id[decoy_id]
        
        binder_db = os.path.join(binder_info['assay_folder_path'], "results", f"assay_{binder_id}.db")
        decoy_db = os.path.join(decoy_info['assay_folder_path'], "results", f"assay_{decoy_id}.db")

        try:
            binder_scores = self._load_best_scores(binder_db)
            decoy_scores = self._load_best_scores(decoy_db)
        except Exception as e:
            print(f"❌ Error loading scores: {e}")
            return None

        if not binder_scores or not decoy_scores:
            print("❌ Missing scores: ensure both assays have populated results databases.")
            return None

        labels = [1] * len(binder_scores) + [0] * len(decoy_scores)
        scores = binder_scores + decoy_scores
        adjusted_scores = [-s for s in scores]  # higher is better binder

        paired = sorted(zip(adjusted_scores, labels), key=lambda x: x[0], reverse=True)
        pos_total = sum(labels)
        neg_total = len(labels) - pos_total
        tp = fp = 0
        roc_points = []

        for score, label in paired:
            if label == 1:
                tp += 1
            else:
                fp += 1
            tpr = tp / pos_total if pos_total else 0.0
            fpr = fp / neg_total if neg_total else 0.0
            roc_points.append((fpr, tpr, score, label))

        auc = 0.0
        prev_fpr = prev_tpr = 0.0
        for fpr, tpr, _, _ in roc_points:
            auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2
            prev_fpr, prev_tpr = fpr, tpr

        analysis_dir = os.path.join(self.path, "docking", "analysis")
        os.makedirs(analysis_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_path = os.path.join(analysis_dir, f"roc_assay_{binder_id}_vs_{decoy_id}_{timestamp}.csv")

        try:
            with open(csv_path, "w") as f:
                f.write("fpr,tpr,score,label\n")
                for fpr, tpr, score, label in roc_points:
                    f.write(f"{fpr},{tpr},{score},{label}\n")
        except Exception as e:
            print(f"⚠️ Could not write ROC CSV: {e}")
            csv_path = None

        print("\nROC analysis complete")
        print(f"Binders: {len(binder_scores)} | Decoys: {len(decoy_scores)} | AUC: {auc:.3f}")
        if csv_path:
            print(f"ROC curve data saved to: {csv_path}")


        # Create a plot of the ROC curve
        if csv_path:
            print(f"ROC curve data saved to: {csv_path}")
            
            # Create ROC plot using matplotlib
            try:
                import matplotlib.pyplot as plt
                import csv
                
                # Read the CSV to extract FPR and TPR for plotting
                fpr_list = []
                tpr_list = []
                with open(csv_path, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        fpr_list.append(float(row['fpr']))
                        tpr_list.append(float(row['tpr']))
                
                # Create the plot
                plt.figure(figsize=(8, 8))
                plt.plot(fpr_list, tpr_list, 'b-', linewidth=2, label=f'ROC curve (AUC = {auc:.3f})')
                plt.plot([0, 1], [0, 1], 'r--', linewidth=1, label='Random classifier')
                plt.xlim([0.0, 1.0])
                plt.ylim([0.0, 1.0])
                plt.xlabel('False Positive Rate', fontsize=12)
                plt.ylabel('True Positive Rate', fontsize=12)
                plt.title(f'ROC Curve - Assay {binder_id} vs {decoy_id}', fontsize=14)
                plt.legend(loc="lower right", fontsize=10)
                plt.grid(True, alpha=0.3)
                
                # Save the plot
                plot_path = csv_path.replace('.csv', '.png')
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"ROC plot saved to: {plot_path}")
                
            except ImportError:
                print("⚠️ matplotlib not available; ROC plot not created")
            except Exception as e:
                print(f"⚠️ Could not create ROC plot: {e}")

        return {
            "binder_assay_id": binder_id,
            "decoy_assay_id": decoy_id,
            "auc": auc,
            "roc_points": roc_points,
            "csv_path": csv_path,
        }
        
    def _prompt_assay(self, role: str, assay_by_id) -> Optional[int]:
        
        while True:
            choice = input(f"Enter assay ID to use as {role} (press Enter to cancel): ").strip()
            if not choice:
                print("❌ Operation cancelled.")
                return None
            try:
                aid = int(choice)
            except ValueError:
                print("Please enter a numeric assay ID.")
                continue
            if aid not in assay_by_id:
                print("Assay ID not found; choose one from the list above.")
                continue
            return aid
        
    def _load_best_scores(self, db_path: str) -> List[float]:
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Results DB not found: {db_path}")
        
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("""
            SELECT MIN(docking_score) AS best_score
            FROM Results
            GROUP BY LigName
        """)
        
        rows = cursor.fetchall()
        conn.close()
        return [float(r[0]) for r in rows if r[0] is not None]
    
    def _parse_pdbqt_output(self, filepath):
        """
        Will read a .pdbqt outfile generated by vina, and will parse the lines starting from 'MODEL 1' to 'ENDMDL'. Starting and ending lines are not included
        """    
        
        try:
            model = []
            with open(filepath, 'r') as f:
                
                in_model = False
                
                for line in f:
                    line = line.rstrip('\n')
                    
                    # Check if we're entering a MODEL block
                    if line.startswith('MODEL 1'):
                        in_model = True
                        continue
                    
                    # Check if we're exiting a MODEL block
                    if line.startswith('ENDMDL'):
                        in_model = False
                        continue
                    
                    # Collect lines within a MODEL block
                    if in_model:
                        model.append(line + "\n")

            return model[6:]  # Exclude the first 6 lines to simulate .dlg file registry
        
        except FileNotFoundError:
            print(f"❌ File not found: {filepath}")
            return None
        except Exception as e:
            print(f"❌ Error parsing PDBQT file: {e}")
            return None
        
    def _add_input_model_in_db(self, db, input_model, outfile):
        """
        Will connect to db and add the string in input model to the Ligands table, matching the row by outfile_prefix
        """
        
        import sqlite3
        import json

    
        outfile_prefix = os.path.basename(outfile).replace('.pdbqt', '')
        
        try:
            conn = sqlite3.connect(db)
            cursor = conn.cursor()
            
            # Update the Ligands table with the input_model for the matching LigName
            cursor.execute("""
                UPDATE Ligands 
                SET input_model = ? 
                WHERE LigName = ?
            """, (json.dumps(input_model), outfile_prefix))
            
            conn.commit()
            conn.close()
            
            return True
            
        except sqlite3.Error as e:
            print(f"❌ Database error while adding input model for {outfile_prefix}: {e}")
            return False
        except Exception as e:
            print(f"❌ Error adding input model for {outfile_prefix}: {e}")
            return False
    
    def _refine_receptor_model(self, tleap_processed_file, processed_pdb_path):
        
        # Minize the receptor model as processed with tleap
        minimized_pdb_file = self._minimize_receptor(tleap_processed_file)

        ## Create residue number matching dictionary
        renumbering_dict = self._construct_resnumbers_matching_dictionary(processed_pdb_path, minimized_pdb_file)

        ## Renumber the minimize pdb file using the dictioary
        renumbered_minimized_pdb_file = self._renumber_minimized_pdb_file(minimized_pdb_file, renumbering_dict)

        ## Add element column to the renumbered minimized pdb file 
        moldf_df = self._add_element_column_to_pdb_file(renumbered_minimized_pdb_file)
        
        ## Delete the renumbered file since it is not used any more
        os.remove(renumbered_minimized_pdb_file)

        # Create a moldf type dictionary with the modified df
        receptor_moldf_dict = {
            '_atom_site': moldf_df
        }
        
        return receptor_moldf_dict
        
    def _add_element_column_to_pdb_file(self, tleap_processed_file):
        
        from moldf import read_pdb
        from moldf import write_pdb
        import pandas as pd
        
        # Create a df from the tleap_receptor file
        df = pd.read_csv(tleap_processed_file, sep=r'\s+', engine='python', header=None)
        
        # Change the column names
        df.columns = ['record_name', 'atom_number', 'atom_name', 'residue_name', 'residue_number', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor']
        # Fill NaN values in the 'atom_number' column with forward fill and then add a cumulative sum to ensure unique atom numbers
        df['atom_number'] = df['atom_number'].ffill().add(df['atom_number'].isna().astype(int).cumsum()).astype(int)
        # Fill NaN values with a default value (e.g., 0) and change the residue_number column to integer type
        df['residue_number'] = df['residue_number'].ffill().fillna(0).astype(int)
        # Fill NaN atom names with an empty string
        df['atom_name'] = df['atom_name'].fillna('')
        # Fil NaN residue names with the same value as above
        df['residue_name'] = df['residue_name'].ffill().fillna('')

        ## Add missing columns according to moldf format
        df.insert(3, 'alt_loc', '')
        df.insert(5, 'chain_id', '')
        df.insert(7, 'insertion', '')
        df.insert(13, 'segment_id', '')
        df.insert(14, 'element_symbol', '')
        df.insert(15, 'charge', '')

        # Create a list of atom names in pdb file corresponding to carbon atoms
        carbon_atoms = ['C', 'CA', 'CB', 'CG', 'CD', 'CE', 'CZ', 'CH', 'CW', 'CV', 'CX', 'CY', 'CZ1', 'CZ2', 'CZ3', 'CE1', 'CE2', 'CE3', 'CD1', 'CD2', 'CG1', 'CG2', 'CH2', 'CH3']

        # Check the df for any atom names that are in the carbon_atoms list and change their element_symbol to 'C'
        for index, row in df.iterrows():
            if row['atom_name'].strip() in carbon_atoms:
                df.at[index, 'element_symbol'] = 'C'

        # Create a list of atom names in pdb file corresponding to nitrogen atoms
        nitrogen_atoms = ['N', 'NA', 'NB', 'NC', 'ND', 'NE', 'NZ', 'NH1', 'NH2', 'ND1', 'ND2', 'NE1', 'NE2']

        # Check the df for any atom names that are in the nitrogen_atoms list and change their element_symbol to 'N'
        for index, row in df.iterrows():
            if row['atom_name'].strip() in nitrogen_atoms:
                df.at[index, 'element_symbol'] = 'N'

        # Create a list of atom names in pdb file corresponding to hydrogen atoms
        hydrogen_atoms = ['H', 'HA', 'HB', 'HC', 'HD', 'HE', 'HZ', 'HH', 'HN', 'HW', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'HB1', 'HB2', 'HB3', 'HG1', 'HG2', 'HG3', 'HD1', 'HD2', 'HD3', 'HE1', 'HE2', 'HE3', 'HG', 'HH1', 'HH2', 'HD11', 'HD12', 'HD13', 'HE21', 'HE22', 'HE23', 'HD21', 'HD22', 'HD23', 'HZ1', 'HZ2', 'HZ3', 'HN1', 'HN2', 'HN3','HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HH11', 'HH12', 'HH13', 'HH21', 'HH22', 'HH23', 'HA2', 'HA3']

        # Check the df for any atom names that are in the hydrogen_atoms list and change their element_symbol to 'H'
        for index, row in df.iterrows():
            if row['atom_name'].strip() in hydrogen_atoms:
                df.at[index, 'element_symbol'] = 'H'

        # Create a list of atom names in pdb file corresponding to oxygen atoms
        oxygen_atoms = ['O', 'OA', 'OB', 'OC', 'OD', 'OE', 'OXT', 'OG', 'OG1', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OH']

        # Check the df for any atom names that are in the oxygen_atoms list and change their element_symbol to 'O'
        for index, row in df.iterrows():
            if row['atom_name'].strip() in oxygen_atoms:
                df.at[index, 'element_symbol'] = 'O'

        # Create a list of atom names in pdb file corresponding to sulfur atoms
        sulfur_atoms = ['S', 'SD', 'SG']

        # Check the df for any atom names that are in the sulfur_atoms list and change their element_symbol to 'S'
        for index, row in df.iterrows():
            if row['atom_name'].strip() in sulfur_atoms:
                df.at[index, 'element_symbol'] = 'S'

        # Create a list of atom names in pdb file corresponding to zinc atoms
        zinc_atoms = ['ZN', 'Z']

        # Check the df for any atom names that are in the zinc_atoms list and change their element_symbol to 'Zn'
        for index, row in df.iterrows():
            if row['atom_name'].strip() in zinc_atoms:
                df.at[index, 'element_symbol'] = 'ZN'
                
        # Return the moldf like dataframe
        return df
    
    def _minimize_receptor(self, tleap_processed_file):
        """
        Will minimize the receptor under processing using sander
        """
        
        import subprocess
        import os
        import tempfile
        import shutil

        print("Minimizing receptor")

        # Create minimization input file in the same directory as tleap_processed_file
        min_in_file = os.path.join(os.path.dirname(tleap_processed_file), "min.in")

        with open(min_in_file, 'w') as f:
            f.write("""Initial minimisation of the receptor
 &cntrl
  imin=1, maxcyc=50, ncyc=25,
  cut=16, ntb=0, igb=1,
 &end
        """)

        # Create a tleap input file in the same directory as tleap_processed_file
        tleap_in_file = os.path.join(os.path.dirname(tleap_processed_file), "minimize.in")

        with open(tleap_in_file, 'w') as f:
            f.write(f"""source leaprc.protein.ff14SB
source leaprc.water.tip3p
HOH = WAT
rec = loadpdb "{tleap_processed_file}"
saveamberparm rec "{os.path.join(os.path.dirname(tleap_processed_file), 'receptor.prmtop')}" "{os.path.join(os.path.dirname(tleap_processed_file), 'receptor.inpcrd')}"
quit
        """)

        # Run tleap to generate prmtop and inpcrd files
        tleap_command = f"tleap -f {tleap_in_file}"
        subprocess.run(tleap_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Set file to run sander to minimize the receptor
        prmtop_file = os.path.join(os.path.dirname(tleap_processed_file), "receptor.prmtop")
        inpcrd_file = os.path.join(os.path.dirname(tleap_processed_file), "receptor.inpcrd")
        min_out_file = os.path.join(os.path.dirname(tleap_processed_file), "min.out")
        min_rst_file = os.path.join(os.path.dirname(tleap_processed_file), "min.rst")
        
        # Run sander to minimize the receptor
        sander_command = f"sander -O -i {min_in_file} -o {min_out_file} -p {prmtop_file} -c {inpcrd_file} -r {min_rst_file}"
        subprocess.run(sander_command, shell=True, check=True)   

        # Convert minimized restart file back to PDB format
        minimized_pdb_file = tleap_processed_file.replace('.pdb', '_minimized.pdb')
        ambpdb_min_command = f"ambpdb -p {prmtop_file} -c {min_rst_file} > {minimized_pdb_file}"
        subprocess.run(ambpdb_min_command, shell=True, check=True)

        # Rename 'WAT' residues to 'HOH' in the minimized PDB file
        with open(minimized_pdb_file, 'r') as file:
            filedata = file.read()

        filedata = filedata.replace(' WAT ', ' HOH ')

        with open(minimized_pdb_file, 'w') as file:
            file.write(filedata)

        # Remove element column from minimized PDB file
        with open(minimized_pdb_file, 'r') as file:
            lines = file.readlines()

        with open(minimized_pdb_file, 'w') as file:
            for line in lines:
                if line.startswith(('ATOM', 'HETATM')):
                    # Keep only columns 0-76 (removes columns 77-80 which contain element symbol)
                    file.write(line[:76] + '\n')
                else:
                    file.write(line)

        return minimized_pdb_file
        
    def _pdb_to_moldf_df(self, pdb_file):
        
        from moldf import read_pdb
        from moldf import write_pdb
        import pandas as pd
        
        # Create a df from the tleap_receptor file
        df = pd.read_csv(pdb_file, sep=r'\s+', engine='python', header=None)
        # Change the column names
        df.columns = ['record_name', 'atom_number', 'atom_name', 'residue_name', 'residue_number', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor', 'element_symbol']
        # Fill NaN values in the 'atom_number' column with forward fill and then add a cumulative sum to ensure unique atom numbers
        df['atom_number'] = df['atom_number'].ffill().add(df['atom_number'].isna().astype(int).cumsum()).astype(int)
        # Fill NaN values with a default value (e.g., 0) and change the residue_number column to integer type
        df['residue_number'] = df['residue_number'].ffill().fillna(0).astype(int)
        # Fill NaN atom names with an empty string
        df['atom_name'] = df['atom_name'].fillna('')
        # Fil NaN residue names with the same value as above
        df['residue_name'] = df['residue_name'].ffill().fillna('')

        ## Add missing columns according to moldf format
        df.insert(3, 'alt_loc', '')
        df.insert(5, 'chain_id', '')
        df.insert(7, 'insertion', '')
        df.insert(13, 'segment_id', '')
        df.insert(15, 'charge', '')
        
        return df

    def _construct_resnumbers_matching_dictionary(self, processed_pdb_path, minimized_pdb_file):

        with open(processed_pdb_path, 'r') as f:
            crystallographic_lines = f.readlines()

        with open(minimized_pdb_file, 'r') as f:
            minimized_lines = f.readlines()

        # Construct a list of keys corresponding the crystallographic numbering
        crystal_key_list = []
        for line in crystallographic_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract residue name and number from crystallographic file
                cryst_res_name = line[17:20].strip()
                cryst_res_num = int(line[22:26].strip())
                key = f"{cryst_res_name}_{cryst_res_num}"
                if key not in crystal_key_list:
                    crystal_key_list.append(key)

        # Construct a list of keys corresponding the tleap renumbered sequence
        tleap_key_list = []
        for line in minimized_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract residue name and number from crystallographic file
                tleap_res_name = line[17:20].strip()
                tleap_res_num = int(line[22:26].strip())
                key = f"{tleap_res_name}_{tleap_res_num}"
                if key not in tleap_key_list:
                    tleap_key_list.append(key)

        # Create a dictionary matching the residue numbers between both lists
        renumbering_dict = dict(zip(tleap_key_list, crystal_key_list))

        return renumbering_dict

    def _renumber_minimized_pdb_file(self, minimized_pdb_file, renumbering_dict):

        """
        Will read the minimized_pdb_file which a pdb file, and the residue name and residue number fields will be evaluated, construction a key that will be searched in renumbering_dict. The value returned will be used to reassing residue number
        """
        
        try:
            import shutil
            
            # Create temporary output file
            temp_out = minimized_pdb_file + ".renumbered"
            
            with open(minimized_pdb_file, 'r') as infile, open(temp_out, 'w') as outfile:
                for line in infile:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Extract residue name and number from the minimized file
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26].strip())
                        
                        # Construct the key to search in renumbering_dict
                        key = f"{res_name}_{res_num}"
                        
                        # Look up the crystallographic numbering
                        if key in renumbering_dict:
                            # Parse the value to get the new residue number
                            crystal_value = renumbering_dict[key]
                            # The value format is "ResName_ResNum"
                            crystal_res_name, crystal_res_num = crystal_value.rsplit('_', 1)
                            new_res_num = int(crystal_res_num)
                            
                            # Reconstruct the line with the new residue number
                            # PDB format: columns 23-26 are residue number (right-justified, 4 characters)
                            new_line = line[:22] + f"{new_res_num:4d}" + line[26:]
                            outfile.write(new_line)
                        else:
                            # If key not found in dictionary, keep original line
                            outfile.write(line)
                    else:
                        # Non-ATOM/HETATM lines are written as-is
                        outfile.write(line)
            
            # Replace original file with renumbered version
            shutil.move(temp_out, minimized_pdb_file)
            print(f"✅ Residue numbers renumbered in {minimized_pdb_file}")
            
            return minimized_pdb_file
            
        except Exception as e:
            print(f"   ⚠️  Error renumbering minimized PDB file: {e}")
            return None