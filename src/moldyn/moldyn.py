from tidyscreen import tidyscreen
import sys
import os
from tidyscreen.databases import DatabaseManager as dbm
from tidyscreen.molecule_management import ligand_management as lm


# Add the parent directory to path to import our local tidyscreen module
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir)

# Import from our local tidyscreen module
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject
from tidyscreen.actionlog.action_logger import log_all_public_methods

@log_all_public_methods
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

    def _remap_project_path(self, stored_path):
        """
        Remap a stored absolute path to the current project location.

        When a project is imported (ActivateProject via import_existing_project())
        the paths stored inside SQLite blobs/columns still reference the directory
        the project lived in when those rows were written. This detects the first
        occurrence of a known top-level project subdirectory inside *stored_path*
        and replaces everything before it with self.path, leaving the relative
        suffix intact. Mirrors MolDock._remap_project_path().
        """
        if not stored_path or os.path.exists(stored_path):
            return stored_path
        for marker in ('/docking/', '/chemspace/', '/ml/', '/dynamics/'):
            pos = stored_path.find(marker)
            if pos != -1:
                relative = stored_path[pos + 1:]  # strip the leading '/'
                return os.path.join(self.path, relative)
        return stored_path

    @staticmethod
    def _prompt(message):
        """
        Clipboard-safe replacement for input() for prompts that contain emoji.

        Problems with plain input(emoji_prompt):
        1. Bracketed-paste mode: terminals wrap pasted text with ESC[200~/ESC[201~
           which appear as literal characters when readline does not strip them.
        2. Emoji in the prompt string: readline counts each emoji as 1 column but
           they render as 2, so its cursor model drifts — backspace and arrow keys
           land in the wrong position.
        3. Writing terminal escape sequences (e.g. \\x1b[?2004l) to stdout before
           input() can interfere with readline's own terminal initialisation and
           prevent it from recognising arrow-key sequences, causing ^[[D etc. to
           be inserted as literal text.

        Fix:
        - Print the emoji-containing message on its own line so the cursor is at
          column 0 before readline starts, giving it an accurate starting position.
        - Use readline.parse_and_bind() to suppress the bracketed-paste sentinels
          inside readline itself rather than fighting the terminal directly.
          This keeps readline's terminal setup intact and arrow keys work normally.
        - Strip any residual escape sequences from the result as a safety net.
        """
        import re as _re
        print(message)
        try:
            import readline as _rl
            # readline ≥ 8.1 exposes this variable directly; fall back to
            # binding the two sentinel sequences to a no-op on older versions.
            try:
                _rl.parse_and_bind('set enable-bracketed-paste off')
            except Exception:
                _rl.parse_and_bind(r'"\e[200~": ""')
                _rl.parse_and_bind(r'"\e[201~": ""')
        except ImportError:
            pass  # no readline — residual-strip below is the only defence
        raw = input('> ')
        raw = _re.sub(r'\x1b\[\d+~', '', raw)
        return raw.strip()

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
                    selection = self._prompt(f"\n🧬 Select MD engine or 'cancel': ")
                    
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
                    method_name = self._prompt(f"\n📝 Enter method name (or 'cancel'): ")

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
            description = self._prompt(f"\n📄 Enter method description (optional): ")
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

    def show_tleap_template_guide(self):
        """Print guidance on how to write a tleap template file for create_md_method_using_tleap_template()."""
        print("\n📂 TLEAP TEMPLATE FILE")
        print("=" * 70)
        print("Provide the path to a working tleap input file for this system.")
        print("The file will be stored verbatim and replayed at assay-setup time")
        print("with the following substitution tokens replaced automatically:")
        print()
        print("  REQUIRED tokens (must appear exactly once each):")
        print("  {{RECEPTOR_PDB}}  → absolute path to receptor_checked.pdb")
        print("                      do NOT hardcode the receptor filename")
        print("  {{LIGAND_PDB}}    → absolute path to the docked-pose PDB for this assay")
        print("  {{MOL2_FILE}}     → absolute path to the GAFF2 mol2 file for the ligand")
        print("  {{FRCMOD_FILE}}   → absolute path to the frcmod file for the ligand")
        print("  {{PRMTOP_OUT}}    → absolute path where tleap writes the output .prmtop")
        print("  {{INPCRD_OUT}}    → absolute path where tleap writes the output .inpcrd")
        print()
        print("  HOW TO WRITE THE TEMPLATE:")
        print("  • Load the protein force field first (e.g. source leaprc.protein.ff19SB)")
        print("  • Load the water/ion model matching the solvent type you will choose")
        print("    (e.g. source leaprc.water.tip3p  for explicit solvent)")
        print("  • Load any custom .lib files with:  loadoff <basename>.lib")
        print("    (attach those files when prompted; reference by basename only)")
        print("  • Load any custom .frcmod files with:  loadamberparams <basename>.frcmod")
        print("  • Load the ligand mol2 and frcmod with:")
        print("      lig = loadmol2 {{MOL2_FILE}}")
        print("      loadamberparams {{FRCMOD_FILE}}")
        print("  • Load the receptor and ligand PDBs and form the complex:")
        print("      rec = loadpdb {{RECEPTOR_PDB}}")
        print("      pose = loadpdb {{LIGAND_PDB}}")
        print("      complex = combine { rec lig }")
        print("    (or combine { rec pose } if the ligand comes from the PDB)")
        print("  • Add solvent box and ions as needed, e.g.:")
        print("      solvatebox complex TIP3PBOX 12.0")
        print("      addions complex Na+ 0")
        print("      addions complex Cl- 0")
        print("  • Write output using the tokens (NOT literal filenames):")
        print("      saveamberparm complex {{PRMTOP_OUT}} {{INPCRD_OUT}}")
        print("  • End the file with:  quit")
        print()
        print("  MINIMAL EXAMPLE (protein + small-molecule ligand, TIP3P explicit solvent):")
        print("  ─" * 35)
        print("  source leaprc.protein.ff19SB")
        print("  source leaprc.gaff2")
        print("  source leaprc.water.tip3p")
        print("  loadmol2 {{MOL2_FILE}}")
        print("  loadamberparams {{FRCMOD_FILE}}")
        print("  rec = loadpdb {{RECEPTOR_PDB}}")
        print("  pose = loadpdb {{LIGAND_PDB}}")
        print("  complex = combine { rec pose }")
        print("  solvatebox complex TIP3PBOX 12.0")
        print("  addions complex Na+ 0")
        print("  addions complex Cl- 0")
        print("  saveamberparm complex {{PRMTOP_OUT}} {{INPCRD_OUT}}")
        print("  quit")
        print("  ─" * 35)
        print()
        print("  CYCLODEXTRIN / CUSTOM FF EXAMPLE (add extra loadoff / loadamberparams):")
        print("  ─" * 35)
        print("  source leaprc.GLYCAM_06j-1")
        print("  source leaprc.gaff2")
        print("  source leaprc.water.tip3p")
        print("  loadoff HP2.lib")
        print("  loadamberparams HP2.frcmod")
        print("  UNL = loadmol2 {{MOL2_FILE}}")
        print("  loadamberparams {{FRCMOD_FILE}}")
        print("  rec = loadpdb {{RECEPTOR_PDB}}")
        print("  pose = loadpdb {{LIGAND_PDB}}")
        print("  complex = combine { rec pose }")
        print("  solvatebox complex TIP3PBOX 12.0")
        print("  addions complex Na+ 0")
        print("  addions complex Cl- 0")
        print("  saveamberparm complex {{PRMTOP_OUT}} {{INPCRD_OUT}}")
        print("  quit")
        print("  ─" * 35)
        print("=" * 70)

    def create_md_method_using_tleap_template(self):
        """
        Create an MD method where system preparation (prmtop/inpcrd) is driven
        by a user-supplied tleap input file rather than interactively-collected
        tleap parameters.  Intended for non-standard systems (e.g. cyclodextrins,
        glycans, cofactors) that require custom force-field loading or explicit
        bond declarations not captured by the generic wizard.

        The template file is stored verbatim in the method record.  At assay time
        the template is replayed with path substitutions for the ligand/receptor
        PDB and the output prmtop/inpcrd file names.

        AMBER simulation parameters (minimization, heating, equilibration,
        production) are still collected interactively — they are independent of
        the system-preparation template.
        """
        try:
            print("🧬 CREATE MD METHOD (TLEAP TEMPLATE)")
            print("=" * 50)

            md_registers_dir = os.path.dirname(self.__md_methods_db)
            os.makedirs(md_registers_dir, exist_ok=True)

            params = {}
            params['engine'] = 'AMBER'
            params['method_type'] = 'tleap_template'

            # --- Method name ---
            while True:
                try:
                    method_name = self._prompt("\n📝 Enter method name (or 'cancel'): ")
                    if method_name.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ MD method creation cancelled")
                        return None
                    if not method_name:
                        print("❌ Method name cannot be empty")
                        continue
                    if not method_name.replace('_', '').replace('-', '').replace(' ', '').isalnum():
                        print("❌ Method name can only contain letters, numbers, spaces, hyphens, and underscores")
                        continue
                    params['method_name'] = method_name
                    break
                except KeyboardInterrupt:
                    print("\n❌ MD method creation cancelled")
                    return None

            # --- Description ---
            description = self._prompt("\n📄 Enter method description (optional): ")
            if not description:
                description = "AMBER MD method using tleap template"
            params['description'] = description

            # --- Tleap template file ---
            self.show_tleap_template_guide()
            while True:
                try:
                    tleap_template_path = self._prompt("\n📂 Path to tleap template file (or 'cancel'): ")
                    if tleap_template_path.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ MD method creation cancelled")
                        return None
                    tleap_template_path = os.path.expanduser(tleap_template_path)
                    if not os.path.isfile(tleap_template_path):
                        print(f"❌ File not found: {tleap_template_path}")
                        continue
                    with open(tleap_template_path, 'r') as fh:
                        tleap_template_content = fh.read()
                    params['tleap_template_path'] = tleap_template_path
                    params['tleap_template_content'] = tleap_template_content
                    n_lines = len(tleap_template_content.splitlines())
                    print(f"✅ Template loaded: {os.path.basename(tleap_template_path)} ({n_lines} lines)")
                    break
                except KeyboardInterrupt:
                    print("\n❌ MD method creation cancelled")
                    return None

            # --- Custom parameter files (.lib / .frcmod) ---
            print("\n📦 CUSTOM PARAMETER FILES")
            print("-" * 70)
            print("Attach any .lib or .frcmod files required by the tleap template.")
            print("They will be stored in the method record and restored to the assay")
            print("folder before tleap runs, so the template can reference them by")
            print("basename only (e.g.  loadoff HP2.lib  or  loadamberparams my.frcmod).")
            print("-" * 70)
            custom_parameter_files = []
            try:
                attach = self._prompt("\n➕ Attach custom parameter files? (yes/no) [default: no]: ").lower() or 'no'
                if attach in ['yes', 'y']:
                    print("Enter file paths one at a time. Type 'done' when finished.")
                    while True:
                        try:
                            file_path_input = self._prompt("  📎 File path (or 'done'): ")
                            if file_path_input.lower() in ['done', '']:
                                break
                            file_path_input = os.path.expanduser(file_path_input)
                            if not os.path.isfile(file_path_input):
                                print(f"  ❌ File not found: {file_path_input}")
                                continue
                            with open(file_path_input, 'r') as fh:
                                file_content = fh.read()
                            entry = {
                                'filename': os.path.basename(file_path_input),
                                'content': file_content,
                            }
                            custom_parameter_files.append(entry)
                            print(f"  ✅ Attached: {entry['filename']}")
                        except KeyboardInterrupt:
                            print("\n❌ MD method creation cancelled")
                            return None
            except KeyboardInterrupt:
                print("\n❌ MD method creation cancelled")
                return None
            params['custom_parameter_files'] = custom_parameter_files
            if custom_parameter_files:
                print(f"✅ {len(custom_parameter_files)} custom file(s) attached.")
            else:
                print("ℹ️  No custom parameter files attached.")

            # --- Solvent type (needed to generate correct AMBER .in files) ---
            solvent_type = input("\nSolvent type (implicit/explicit) [default: explicit]: ").strip().lower() or 'explicit'
            if solvent_type not in ['implicit', 'explicit']:
                print("⚠️  Invalid solvent type, defaulting to 'explicit'")
                solvent_type = 'explicit'
            params['solvent_type'] = solvent_type

            # --- AMBER simulation parameters (independent of tleap template) ---
            params = self._set_minimization_params(params)
            params = self._set_heating_params(params)
            params = self._set_equilibration_params(params)
            params = self._set_production_params(params)

            # --- Persist to DB ---
            params = self._save_parameters_to_db(params)

            # --- Summary ---
            self._display_method_summary_tleap_template(params)

        except Exception as e:
            print(f"❌ Error in create_md_method_using_tleap_template: {e}")
            return None

    def _display_method_summary_tleap_template(self, params):
        """Display summary for a tleap-template-based MD method."""
        try:
            print(f"\n🎯 MOLECULAR DYNAMICS METHOD CREATED SUCCESSFULLY")
            print("=" * 60)
            print(f"\n🏷️  Method ID   : {params.get('method_id', 'N/A')}")
            print(f"🏷️  Method Name : {params.get('method_name', 'Unnamed_Method')}")
            print(f"🧬 MD Engine   : {params.get('engine', '')}")
            print(f"📝 Description : {params.get('description', '')}")
            print(f"🔖 Method type : tleap template")

            template_path = params.get('tleap_template_path', 'N/A')
            template_content = params.get('tleap_template_content', '')
            n_lines = len(template_content.splitlines()) if template_content else 0
            print(f"\n📂 TLEAP TEMPLATE:")
            print(f"  Source file : {template_path}")
            print(f"  Lines stored: {n_lines}")

            custom_files = params.get('custom_parameter_files', [])
            print(f"\n📦 ATTACHED PARAMETER FILES ({len(custom_files)}):")
            if custom_files:
                for cf in custom_files:
                    n = len(cf.get('content', '').splitlines())
                    print(f"  - {cf['filename']} ({n} lines)")
            else:
                print("  (none)")

            print(f"\n⚙️  SIMULATION PARAMETERS:")

            min1 = params.get('first_minimization', {})
            min2 = params.get('second_minimization', {})
            general = params.get('general_params', {})
            print(f"  🧪 MINIMIZATION:")
            print(f"    🔹 Step 1 — maxcyc={min1.get('min1_maxcyc')}, ncyc={min1.get('min1_ncyc')}, "
                  f"restraint_wt={min1.get('min1_restraint_wt')}")
            print(f"    🔹 Step 2 — maxcyc={min2.get('min2_maxcyc')}, ncyc={min2.get('min2_ncyc')}, "
                  f"restraint_wt={min2.get('min2_restraint_wt')}")

            heating = params.get('heating_params', {})
            print(f"  🔥 HEATING:")
            print(f"    {heating.get('heating_steps')} steps, "
                  f"{heating.get('initial_temp')} K → {general.get('target_temp')} K")

            equil = params.get('equilibration_params', {})
            print(f"  ❄️  EQUILIBRATION: {equil.get('equilibration_steps')} steps")

            prod = params.get('production_params', {})
            print(f"  🎬 PRODUCTION   : {prod.get('production_steps')} steps")

            print("=" * 60)

        except Exception as e:
            print(f"⚠️  Error displaying method summary: {e}")

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

        # For stub receptors the tleap topology setup (force field, water model,
        # box, ions) is already encoded in the template stored at import time.
        # Only the AMBER simulation parameters (minimisation, heating, etc.) are needed.
        stub_answer = input(
            "\nIs this method for a stub receptor imported via import_receptor_model()? (yes/no) [default: no]: "
        ).strip().lower() or 'no'
        uses_receptor_template = stub_answer in ('yes', 'y')

        if uses_receptor_template:
            print(
                "ℹ️  Tleap topology parameters will be skipped — "
                "they are provided by the receptor's stored template."
            )
            params['tleap_params'] = {'uses_receptor_template': True}
        else:
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

            conn = sqlite3.connect(self.__md_methods_db)
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
        import traceback

        md_assay_folder = None
        md_assay_id = None

        try:
            import sqlite3

            # Connect to docking registers database
            conn = sqlite3.connect(self.__docking_registers_db)
            cursor = conn.cursor()
            # Show docking assays available
            docking_assay_params_dict = self._select_docking_assay(cursor)
            conn.close()
            if docking_assay_params_dict is None:
                return

            # Select unique ligand molecules in the selected docking assay
            docking_assay_params_dict = self._select_unique_ligands_in_docking_assay(docking_assay_params_dict)
            if docking_assay_params_dict is None:
                return

            # Select pose id for the selected ligand to perform MD assay
            docking_assay_params_dict = self._select_ligand_pose_for_md_assay(docking_assay_params_dict)
            if docking_assay_params_dict is None:
                return

            # Select a MD method to perform the MD assay
            md_parameters_dict = self._select_md_method()
            if md_parameters_dict is None:
                return

            # Restore the selected docked pose using ligand_management module
            print(f"\n⚙️  Restoring docked pose...")
            pdb_dict = lm.restore_single_docked_pose(docking_assay_params_dict.get('docking_results_db'), docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict.get('selected_pose_id'))
            if pdb_dict is None:
                raise RuntimeError("Failed to restore docked pose — restore_single_docked_pose returned None")

            # Query user for MD assay description
            assay_description = self._prompt("\n📝 Enter MD assay description (optional): ")
            if not assay_description:
                assay_description = f"MD assay for ligand {docking_assay_params_dict.get('selected_ligand_name')} (pose {docking_assay_params_dict.get('selected_pose_id')})"

            # Create MD assay folder structure and prepare input files for the selected MD engine
            md_assay_id, md_assay_folder = self._create_md_assay_folder()
            if md_assay_id is None:
                return

            # Save the docked pose PDB file in the MD assay folder
            print(f"⚙️  Writing docked pose PDB...")
            selected_pose_pdb = lm.write_pdb_with_moldf(pdb_dict, docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict.get('selected_pose_id'), md_assay_folder)
            if selected_pose_pdb is None:
                raise RuntimeError("Failed to write docked pose PDB — write_pdb_with_moldf returned None")

            # Prepare the ligand tleap input file in the MD assay folder
            print(f"⚙️  Preparing ligand force-field files...")
            mol2_file, frcmod_file = lm.prepare_ligand_tleap_input_files(self.__chemspace_db, docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict, md_assay_folder)
            if mol2_file is None or frcmod_file is None:
                raise RuntimeError("Failed to prepare ligand tleap input files — prepare_ligand_tleap_input_files returned None")

            # Create complex .prmtop and .inpcrd files
            print(f"⚙️  Preparing topology and coordinate files with tleap...")

            # Resolve receptor_info (may be a raw JSON string in older records)
            import json as _json
            _receptor_info = docking_assay_params_dict.get('receptor_info', {})
            if isinstance(_receptor_info, str):
                try:
                    _receptor_info = _json.loads(_receptor_info)
                except Exception:
                    _receptor_info = {}

            # Derive receptor PDB path (needed by both template branches)
            _pdbqt = self._remap_project_path(_receptor_info.get('pdbqt_file', ''))
            _receptor_pdb = os.path.join(os.path.dirname(_pdbqt), 'receptor_checked.pdb') if _pdbqt else None

            if md_parameters_dict.get('method_type') == 'tleap_template':
                # MD method carries its own tleap template (explicit user choice)
                prmtop_file, inpcrd_file = self._prepare_complex_prmtop_inpcrd_from_template(
                    selected_pose_pdb, md_assay_folder, md_parameters_dict,
                    mol2_file=mol2_file, frcmod_file=frcmod_file,
                    receptor_pdb=_receptor_pdb,
                )
            else:
                # Check if the receptor is a stub (created by import_receptor_model).
                # If so, use the receptor's stored tleap_config with solvation intact.
                _stub_tleap_config = None
                _is_stub = False
                _template_name = _receptor_info.get('template_name', '')
                if _template_name:
                    import sqlite3 as _sqlite3
                    _pdbs_db = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
                    try:
                        _conn = _sqlite3.connect(_pdbs_db)
                        _row = _conn.execute(
                            "SELECT status FROM pdb_templates WHERE pdb_template_name = ?",
                            (_template_name,)
                        ).fetchone()
                        _conn.close()
                        if _row and _row[0] == 'imported':
                            _is_stub = True
                            _receptors_db = os.path.join(self.path, 'docking', 'receptors', 'receptors.db')
                            _model_name = _receptor_info.get('receptor_model_name', '')
                            _conn2 = _sqlite3.connect(_receptors_db)
                            _row2 = _conn2.execute(
                                "SELECT tleap_config, fps_tleap_config FROM receptor_models "
                                "WHERE receptor_model_name = ?",
                                (_model_name,)
                            ).fetchone()
                            _conn2.close()
                            if _row2:
                                _raw = _row2[0] or _row2[1]
                                if _raw:
                                    try:
                                        _stub_tleap_config = _json.loads(_raw)
                                    except Exception:
                                        pass
                    except Exception:
                        pass

                if _stub_tleap_config:
                    print(f"   ℹ️  Stub receptor detected — using stored tleap template "
                          f"(solvation commands kept for MD).")
                    _stub_md_params = {
                        'tleap_template_content': _stub_tleap_config.get('template_content', ''),
                        'custom_parameter_files': _stub_tleap_config.get('custom_parameter_files', []),
                    }
                    prmtop_file, inpcrd_file = self._prepare_complex_prmtop_inpcrd_from_template(
                        selected_pose_pdb, md_assay_folder, _stub_md_params,
                        mol2_file=mol2_file, frcmod_file=frcmod_file,
                        receptor_pdb=_receptor_pdb,
                    )
                else:
                    # Compatibility check: method created for stub receptors but receptor is standard.
                    _uses_receptor_template = md_parameters_dict.get('tleap_params', {}).get('uses_receptor_template', False)
                    if _uses_receptor_template and not _is_stub:
                        print("❌ Incompatible MD method and receptor:")
                        print("   The selected MD method was created for an imported stub receptor")
                        print("   (tleap topology parameters were not collected).")
                        print("   The selected docking assay uses a standard receptor that requires")
                        print("   explicit tleap parameters (force field, water model, box, ions).")
                        print("   Solution: create a new MD method with create_md_method() and")
                        print("   answer 'no' when asked about stub receptors.")
                        return
                    if _is_stub:
                        print("⚠️  Stub receptor detected but no tleap template stored.")
                        print("   Re-run import_receptor_model() to provide a template,")
                        print("   or use create_md_method_using_tleap_template() for this assay.")
                    prmtop_file, inpcrd_file = self._prepare_complex_prmtop_inpcrd_for_md(
                        mol2_file, frcmod_file, docking_assay_params_dict,
                        md_assay_folder, md_assay_folder, md_parameters_dict, selected_pose_pdb,
                    )

            # Prepare md simulation input files
            print(f"⚙️  Writing AMBER simulation input files...")
            self._prepare_md_simulation_input_files(md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file)

            # Prepare execution scripts for the MD assay
            print(f"⚙️  Writing execution script...")
            self._prepare_md_execution_script(md_assay_folder, md_parameters_dict)

            # Create MD register entry in the md_assays table
            self._create_md_assay_register_entry(md_assay_id, md_assay_folder, assay_description, docking_assay_params_dict.get('assay_id'), docking_assay_params_dict.get('selected_ligand_name'), docking_assay_params_dict.get('selected_pose_id'), md_parameters_dict)

            print(f"\n✅ MD assay prepared successfully in: {md_assay_folder}")

            # Query if the user want to start the MD simulation now
            start_now = input("\n▶️  Do you want to start the MD simulation now? (yes/no) [default: no]: ").strip().lower() or 'no'
            bg_info = None
            if start_now in ['yes', 'y']:
                bg_info = self._start_md_simulation(md_assay_folder)
            else:
                print("ℹ️  MD simulation not started. Run the execution script in the assay folder when ready.")

            return bg_info

        except Exception as e:
            print(f"\n❌ Error performing MD assay: {e}")
            traceback.print_exc()
            if md_assay_folder and os.path.exists(md_assay_folder):
                print(f"⚠️  Assay folder left on disk for inspection: {md_assay_folder}")
        
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
            assay_folder_path = self._remap_project_path(docking_assay_params_dict.get('assay_folder_path'))
            docking_assay_params_dict['assay_folder_path'] = assay_folder_path
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
        except FileExistsError:
            print(f"⚠️  MD assay folder already exists: {new_assay_folder}")
            confirm = input("🗑️  Delete existing folder and continue? (yes/no) [default: no]: ").strip().lower() or 'no'
            if confirm not in ['yes', 'y']:
                print("❌ MD assay creation cancelled.")
                return None, None
            import shutil
            shutil.rmtree(new_assay_folder)
            os.makedirs(new_assay_folder)
        except Exception as e:
            print(f"❌ Error creating MD assay folder: {e}")
            return None, None
        print(f"✅ Created MD assay folder: {new_assay_folder}")
        return new_assay_id, new_assay_folder

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

    def _create_md_assay_register_entry(self, new_assay_id, new_assay_folder, assay_description, docking_assay_id=None, ligname=None, pose_id=None, md_parameters_dict=None, receptor_template_name=None):
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
                'receptor_template_name': 'TEXT',
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
                'receptor_template_name': receptor_template_name,
                'md_parameters': self._serialize_parameters(md_parameters_dict),
            }
            # Insert dinamically data into md_assays table
            dbm.insert_data_dinamically_into_table(cursor, 'md_assays', data_dict)
            conn.commit()
            conn.close()

        except Exception as e:
            print(f"❌ Error creating MD assay register entry: {e}")

    def _check_disulfides_in_template(self, pdb_template_name):
        """
        Retrieve the confirmed disulfide_bonds list stored for pdb_template_name
        in pdbs.db (persisted inside the pdb_analysis JSON blob by
        MolDock.create_pdb_template()). Mirrors MolDock._check_disulfides_in_template().

        Returns:
            List[Dict]: disulfide bond entries (chain1, resnum1, chain2,
            resnum2, distance), or [] if none are recorded / on any error.
        """
        import sqlite3
        import json

        try:
            pdb_templates_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
            conn = sqlite3.connect(pdb_templates_db_path)
            cursor = conn.cursor()

            cursor.execute("SELECT pdb_analysis FROM pdb_templates WHERE pdb_template_name = ?", (pdb_template_name,))
            result = cursor.fetchone()
            conn.close()

            if result is None or result[0] is None:
                return []

            pdb_analysis = json.loads(result[0])
            return pdb_analysis.get('disulfide_bonds', []) or []

        except Exception as e:
            print(f"Error retrieving disulfide bonds for template '{pdb_template_name}': {e}")
            return []

    def _get_tleap_bond_lines(self, pdb_path, disulfide_bonds, unit_name, target_already_tleap_processed=False):
        """
        Compute each disulfide-bonded residue's tleap-internal residue index
        and return the corresponding tleap `bond` command lines for
        `unit_name`. Mirrors MolDock._get_tleap_bond_lines() — see that
        method for the full rationale.

        tleap's `unit.<idx>` addressing on a FRESH loadpdb is NOT a simple
        1-based ordinal count of distinct residues in file order — that was
        the original, incorrect assumption here and caused tleap parser
        errors on receptors with resSeq gaps (missing/unmodeled loop
        residues) or multiple chains. Empirically verified actual rule for a
        fresh load:

        - Within a single chain, tleap's internal residue index tracks the
          PDB resSeq value directly: when tleap detects a break in the chain
          (a resSeq gap), it jumps the internal counter by the same gap so
          the index keeps matching the literal resSeq for the rest of that
          chain. Net effect: for the first chain loaded, internal index ==
          literal resSeq.
        - At a chain boundary (new chain ID, e.g. separated by a TER record),
          the internal counter advances by exactly 1 from the previous
          chain's last index, regardless of the new chain's own resSeq
          numbering (it does NOT jump to align with it).

        `target_already_tleap_processed` selects which file the computed
        indices must address (see MolDock._get_tleap_bond_lines() for the
        full explanation of the True case, used when pdb_path is an
        earlier-stage file standing in for a different, already
        loadpdb+savepdb-processed file).

        Returns:
            List[str]: ready-to-write tleap script lines (newline-terminated)
        """
        if not disulfide_bonds:
            return []

        residue_keys = set()
        for bond in disulfide_bonds:
            residue_keys.add((bond['chain1'], bond['resnum1']))
            residue_keys.add((bond['chain2'], bond['resnum2']))

        resnum_counts = {}
        for _chain_id, resnum in residue_keys:
            resnum_counts[resnum] = resnum_counts.get(resnum, 0) + 1
        ambiguous_resnums = {r for r, n in resnum_counts.items() if n > 1}

        index_by_chain_resnum = {}
        index_by_resnum = {}
        seen_residue = None
        current_chain = None
        chain_base_index = 0
        chain_first_resnum = None
        tleap_index = 0
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    chain_id = line[21].strip()
                    try:
                        resnum = int(line[22:26])
                    except ValueError:
                        continue
                    res_uid = (chain_id, resnum)
                    if res_uid != seen_residue:
                        seen_residue = res_uid
                        if target_already_tleap_processed:
                            tleap_index += 1
                        elif chain_id != current_chain:
                            tleap_index = resnum if current_chain is None else tleap_index + 1
                            current_chain = chain_id
                            chain_base_index = tleap_index
                            chain_first_resnum = resnum
                        else:
                            tleap_index = chain_base_index + (resnum - chain_first_resnum)
                        index_by_chain_resnum.setdefault(res_uid, tleap_index)
                        index_by_resnum.setdefault(resnum, tleap_index)

        def _resolve(chain_id, resnum):
            idx = index_by_chain_resnum.get((chain_id, resnum))
            if idx is not None:
                return idx
            if resnum not in ambiguous_resnums:
                return index_by_resnum.get(resnum)
            return None

        bond_lines = []
        for bond in disulfide_bonds:
            idx1 = _resolve(bond['chain1'], bond['resnum1'])
            idx2 = _resolve(bond['chain2'], bond['resnum2'])
            if idx1 is None or idx2 is None:
                print(
                    f"   ⚠️  Could not locate disulfide partner residue(s) for "
                    f"CYX{bond['resnum1']}_{bond['chain1']}--CYX{bond['resnum2']}_{bond['chain2']} "
                    f"in {pdb_path}; skipping explicit tleap bond for this pair."
                )
                continue
            bond_lines.append(f"bond {unit_name}.{idx1}.SG {unit_name}.{idx2}.SG\n")

        return bond_lines

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

        # get the .pdbqt used for docking, remapped to the current project location
        # so that imported projects (whose original path differs) work correctly.
        receptor_pdb = self._remap_project_path(receptor_details.get('pdbqt_file', None))

        # Define the raw .pdb file of the receptor
        receptor_pdb_path = os.path.join(os.path.dirname(receptor_pdb), 'receptor_checked.pdb')

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
        
        tleap_log_file = os.path.join(output_dir, 'tleap.log')

        # Confirmed disulfide bridges (if any) for this receptor template
        pdb_template_name = receptor_details.get('template_name', None)
        disulfide_bonds = self._check_disulfides_in_template(pdb_template_name)
        disulfide_bond_lines = self._get_tleap_bond_lines(receptor_pdb_path, disulfide_bonds, 'rec')

        try:
            with open(tleap_in_file, 'w') as f:
                f.write(f"source leaprc.protein.{protein_ff}\n")
                f.write(f"loadamberparams frcmod.ions1lm_126_tip3p\n")
                f.write(f"source leaprc.{ligand_ff}\n")
                if solvent_type == 'explicit':
                    f.write(f"source leaprc.water.{water_model}\n")
                    f.write(f"HOH = WAT\n")
                f.write(f"rec = loadpdb {receptor_pdb_path}\n")
                for bond_line in disulfide_bond_lines:
                    f.write(bond_line)
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
        except OSError as e:
            raise RuntimeError(f"Failed to write tleap input file '{tleap_in_file}': {e}") from e

        # Run tleap to generate prmtop and inpcrd files
        try:
            result = subprocess.run(
                f"tleap -f {tleap_in_file}",
                shell=True,
                capture_output=True,
                text=True,
            )
        except subprocess.SubprocessError as e:
            raise RuntimeError(f"Failed to execute tleap: {e}") from e

        # Always save tleap output for later inspection
        try:
            with open(tleap_log_file, 'w') as lf:
                lf.write(result.stdout)
                if result.stderr:
                    lf.write("\n--- STDERR ---\n")
                    lf.write(result.stderr)
        except OSError:
            pass

        if result.returncode != 0:
            tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
            print(f"   tleap output (last lines):\n{tail}")
            raise RuntimeError(
                f"tleap exited with code {result.returncode}. "
                f"Full log: {tleap_log_file}"
            )

        # tleap exits 0 even when it fails to write output files — verify explicitly
        missing = [f for f in (prmtop_file, inpcrd_file) if not os.path.exists(f)]
        if missing:
            tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
            print(f"   tleap output (last lines):\n{tail}")
            raise RuntimeError(
                f"tleap completed but did not produce: {', '.join(os.path.basename(f) for f in missing)}. "
                f"Full log: {tleap_log_file}"
            )

        print(f"   ✓ tleap completed successfully (log: {tleap_log_file})")
        return prmtop_file, inpcrd_file

    def _prepare_complex_prmtop_inpcrd_from_template(self, selected_pose_pdb, output_dir, md_parameters_dict, mol2_file=None, frcmod_file=None, receptor_pdb=None):
        """
        Replay the tleap template stored in the MD method record, substituting
        path tokens, then run tleap to produce prmtop/inpcrd.

        Used by both perform_md_assay() (ligand+receptor) and
        perform_md_assay_on_receptor() (receptor only; selected_pose_pdb=None).

        Tokens replaced in the template:
          {{RECEPTOR_PDB}} → full path to receptor_checked.pdb inside output_dir (file is copied there; surrounding single quotes in the template are also stripped)
          {{LIGAND_PDB}}   → absolute path to the docked-pose PDB (skipped if None)
          {{MOL2_FILE}}    → absolute path to the GAFF2 mol2 file for the ligand
          {{FRCMOD_FILE}}  → absolute path to the frcmod file for the ligand
          {{PRMTOP_OUT}}   → complex.prmtop inside output_dir
          {{INPCRD_OUT}}   → complex.inpcrd inside output_dir

        Raises RuntimeError on any failure; tleap stdout+stderr are always
        saved to tleap.log inside output_dir for post-mortem inspection.
        """
        import subprocess

        template_content = md_parameters_dict.get('tleap_template_content', '')
        if not template_content:
            raise RuntimeError(
                "MD method has no 'tleap_template_content' stored. "
                "Re-create the method with create_md_method_using_tleap_template()."
            )

        prmtop_file = os.path.join(output_dir, 'complex.prmtop')
        inpcrd_file = os.path.join(output_dir, 'complex.inpcrd')
        tleap_in_file = os.path.join(output_dir, 'complex.in')
        tleap_log_file = os.path.join(output_dir, 'tleap.log')

        # Restore custom parameter files (.lib, .frcmod) into the assay folder so
        # tleap can find them by basename, as written in the template.
        custom_parameter_files = md_parameters_dict.get('custom_parameter_files', [])
        for cf in custom_parameter_files:
            dest = os.path.join(output_dir, cf['filename'])
            try:
                with open(dest, 'w') as fh:
                    fh.write(cf['content'])
                print(f"   ✓ Restored parameter file: {cf['filename']}")
            except OSError as e:
                raise RuntimeError(
                    f"Failed to restore parameter file '{cf['filename']}' to assay folder: {e}"
                ) from e

        # Copy the receptor PDB into the assay folder and resolve its full path
        # so that complex.in can reference it unambiguously without quotes.
        receptor_pdb_dest = None
        if receptor_pdb and os.path.exists(receptor_pdb):
            import shutil as _shutil
            receptor_pdb_dest = os.path.join(output_dir, os.path.basename(receptor_pdb))
            try:
                _shutil.copy2(receptor_pdb, receptor_pdb_dest)
                print(f"   ✓ Copied receptor PDB to assay folder: {receptor_pdb_dest}")
            except OSError as e:
                raise RuntimeError(
                    f"Failed to copy receptor PDB '{receptor_pdb}' to assay folder: {e}"
                ) from e
        elif receptor_pdb:
            raise RuntimeError(f"Receptor PDB not found: {receptor_pdb}")

        content = template_content
        if receptor_pdb_dest:
            # Replace quoted token first, then unquoted — both map to the full path without quotes.
            content = content.replace("'{{RECEPTOR_PDB}}'", receptor_pdb_dest)
            content = content.replace('{{RECEPTOR_PDB}}', receptor_pdb_dest)
        if selected_pose_pdb:
            content = content.replace('{{LIGAND_PDB}}', selected_pose_pdb)
        content = content.replace('{{PRMTOP_OUT}}', prmtop_file)
        content = content.replace('{{INPCRD_OUT}}', inpcrd_file)
        if mol2_file:
            content = content.replace('{{MOL2_FILE}}', mol2_file)
        if frcmod_file:
            content = content.replace('{{FRCMOD_FILE}}', frcmod_file)

        try:
            with open(tleap_in_file, 'w') as f:
                f.write(content)
        except OSError as e:
            raise RuntimeError(f"Failed to write tleap input file '{tleap_in_file}': {e}") from e

        try:
            result = subprocess.run(
                f"tleap -f {tleap_in_file}",
                shell=True,
                capture_output=True,
                text=True,
                cwd=output_dir,
            )
        except subprocess.SubprocessError as e:
            raise RuntimeError(f"Failed to execute tleap: {e}") from e

        try:
            with open(tleap_log_file, 'w') as lf:
                lf.write(result.stdout)
                if result.stderr:
                    lf.write("\n--- STDERR ---\n")
                    lf.write(result.stderr)
        except OSError:
            pass

        if result.returncode != 0:
            tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
            print(f"   tleap output (last lines):\n{tail}")
            raise RuntimeError(
                f"tleap exited with code {result.returncode}. "
                f"Full log: {tleap_log_file}"
            )

        missing = [f for f in (prmtop_file, inpcrd_file) if not os.path.exists(f)]
        if missing:
            tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
            print(f"   tleap output (last lines):\n{tail}")
            raise RuntimeError(
                f"tleap completed but did not produce: "
                f"{', '.join(os.path.basename(f) for f in missing)}. "
                f"Full log: {tleap_log_file}"
            )

        print(f"   ✓ tleap completed successfully via template (log: {tleap_log_file})")
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
            raise RuntimeError(f"Failed to prepare MD simulation input files: {e}") from e

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
        eq_restraint_selector = md_parameters_dict.get('equilibration_params').get('equilibration_restraint_selector', '')
        with open(equilibration_in_file, 'w') as f:
            f.write(f"Equilibration stage\n")
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
            f.write(f"  tempi={md_parameters_dict.get('general_params').get('target_temp')},\n")
            f.write(f"  ntt={md_parameters_dict.get('general_params').get('thermostat')},\n")
            f.write(f"  gamma_ln={md_parameters_dict.get('general_params').get('collision_freq')},\n")
            f.write(f"  nstlim={md_parameters_dict.get('equilibration_params').get('equilibration_steps')},\n")
            f.write(f"  dt=0.002,\n")
            f.write(f"  ntpr={md_parameters_dict.get('equilibration_params').get('equilibration_output_write')},\n")
            f.write(f"  ntwx={md_parameters_dict.get('equilibration_params').get('equilibration_trajectory_write')},\n")
            f.write(f"  ntwr={md_parameters_dict.get('equilibration_params').get('equilibration_restart_write')},\n")
            f.write(f"  iwrap=1,\n")
            if eq_restraint_selector:
                f.write(f"  ntr=1,\n")
                f.write(f"  restraint_wt={md_parameters_dict.get('equilibration_params').get('equilibration_restraint_wt')},\n")
                f.write(f"  restraintmask='{eq_restraint_selector}',\n")
            f.write(f"/\n")
            f.write(f"END\n")

        # Prepare the production MD input file
        production_in_file = os.path.join(md_assay_folder, "production.in")
        prod_restraint_selector = md_parameters_dict.get('production_params').get('production_restraint_selector', '')
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
            if prod_restraint_selector:
                f.write(f"  ntr=1,\n")
                f.write(f"  restraint_wt={md_parameters_dict.get('production_params').get('production_restraint_wt')},\n")
                f.write(f"  restraintmask='{prod_restraint_selector}',\n")
            f.write(f"/\n")
            f.write(f"END\n")

    def _prepare_md_execution_script(self, md_assay_folder, md_parameters_dict):
        import shutil as _shutil

        pmemd = _shutil.which("pmemd.cuda") or "pmemd.cuda"

        min1_ref = md_parameters_dict.get('first_minimization', {}).get('min1_restraint_selector', '')
        min2_ref = md_parameters_dict.get('second_minimization', {}).get('min2_restraint_selector', '')
        heat_ref = md_parameters_dict.get('heating_params', {}).get('heating_restraint_selector', '')
        eq_ref   = md_parameters_dict.get('equilibration_params', {}).get('equilibration_restraint_selector', '')
        prod_ref = md_parameters_dict.get('production_params', {}).get('production_restraint_selector', '')

        gp = md_parameters_dict.get('general_params', {})
        hp = md_parameters_dict.get('heating_params', {})

        # -ref <file> is a coordinate file path, not the mask.
        # Empty string means no -ref flag (no restraints for that stage).
        ref_file_min1 = '"complex.inpcrd"' if min1_ref else '""'
        ref_file_min2 = '"${min1_out}"'    if min2_ref else '""'
        ref_file_heat = '"${min2_out}"'    if heat_ref else '""'
        ref_file_eq   = '"${heating_out}"' if eq_ref   else '""'
        ref_file_prod = '"${eq_out}"'      if prod_ref else '""'

        execution_script_path = os.path.join(md_assay_folder, "run_md.sh")
        try:
            with open(execution_script_path, 'w') as f:

                # --- Header ---
                f.write("#!/bin/bash\n")
                f.write("set -uo pipefail\n")
                f.write('cd "$(dirname "$0")"\n')
                f.write(f"export PATH={os.path.dirname(pmemd)}:$PATH\n\n")
                f.write("MAX_RESTARTS=5\n")
                f.write("STAGE_OUT=\n\n")

                # --- Engine selection ---
                f.write("if command -v pmemd.cuda &>/dev/null; then\n")
                f.write("    MD_ENGINE=pmemd.cuda\n")
                f.write('    echo "Using pmemd.cuda"\n')
                f.write("else\n")
                f.write("    MD_ENGINE=sander\n")
                f.write('    echo "pmemd.cuda not found, falling back to sander"\n')
                f.write("fi\n\n")

                # --- run_stage function ---
                f.write("# run_stage STAGE BASE_IN RESTART_IN COORD_IN REF_COORD WRITE_TRAJ\n")
                f.write("#   Runs one MD stage.  On failure retries from the ntwr checkpoint\n")
                f.write("#   (written to the -r file) up to MAX_RESTARTS times.\n")
                f.write("#   REF_COORD=\"\" for stages without positional restraints.\n")
                f.write("#   WRITE_TRAJ=\"yes\" to emit -x <stage>.nc.\n")
                f.write("#   On success sets STAGE_OUT to the final .crd path.\n")
                f.write("run_stage() {\n")
                f.write('    local stage="$1" base_in="$2" restart_in="$3"\n')
                f.write('    local coord_in="$4" ref_coord="$5" write_traj="$6"\n')
                f.write("    local attempt=0 md_status suffix cur_in coord_out\n\n")
                f.write('    while [[ $attempt -le $MAX_RESTARTS ]]; do\n')
                f.write('        if [[ $attempt -eq 0 ]]; then\n')
                f.write('            suffix="${stage}"; cur_in="${base_in}"\n')
                f.write("        else\n")
                f.write('            suffix="${stage}_restart${attempt}"; cur_in="${restart_in}"\n')
                f.write("        fi\n")
                f.write('        coord_out="${suffix}.crd"\n\n')
                f.write("        md_status=0\n")
                f.write('        if [[ -n "${ref_coord}" && "${write_traj}" == "yes" ]]; then\n')
                f.write('            "${MD_ENGINE}" -O -i "${cur_in}" -o "${suffix}.out" \\\n')
                f.write('                -p complex.prmtop -c "${coord_in}" \\\n')
                f.write('                -ref "${ref_coord}" -r "${coord_out}" -x "${suffix}.nc" \\\n')
                f.write("                || md_status=$?\n")
                f.write('        elif [[ -n "${ref_coord}" ]]; then\n')
                f.write('            "${MD_ENGINE}" -O -i "${cur_in}" -o "${suffix}.out" \\\n')
                f.write('                -p complex.prmtop -c "${coord_in}" \\\n')
                f.write('                -ref "${ref_coord}" -r "${coord_out}" \\\n')
                f.write("                || md_status=$?\n")
                f.write('        elif [[ "${write_traj}" == "yes" ]]; then\n')
                f.write('            "${MD_ENGINE}" -O -i "${cur_in}" -o "${suffix}.out" \\\n')
                f.write('                -p complex.prmtop -c "${coord_in}" \\\n')
                f.write('                -r "${coord_out}" -x "${suffix}.nc" \\\n')
                f.write("                || md_status=$?\n")
                f.write("        else\n")
                f.write('            "${MD_ENGINE}" -O -i "${cur_in}" -o "${suffix}.out" \\\n')
                f.write('                -p complex.prmtop -c "${coord_in}" \\\n')
                f.write('                -r "${coord_out}" \\\n')
                f.write("                || md_status=$?\n")
                f.write("        fi\n\n")
                f.write('        if [[ $md_status -eq 0 ]]; then\n')
                f.write('            echo "  ${stage} done (attempt $((attempt+1))) → ${coord_out}"\n')
                f.write('            STAGE_OUT="${coord_out}"\n')
                f.write("            return 0\n")
                f.write("        fi\n\n")
                f.write('        if [[ -f "${coord_out}" && -s "${coord_out}" ]]; then\n')
                f.write('            echo "  ${stage} attempt $((attempt+1)) failed — checkpoint found, restarting"\n')
                f.write('            coord_in="${coord_out}"\n')
                f.write("            (( attempt++ )) || true\n")
                f.write("        else\n")
                f.write('            echo "  ${stage} attempt $((attempt+1)) failed — no checkpoint written, giving up"\n')
                f.write("            return 1\n")
                f.write("        fi\n")
                f.write("    done\n\n")
                f.write('    echo "  ${stage} exhausted ${MAX_RESTARTS} restart(s) — giving up"\n')
                f.write("    return 1\n")
                f.write("}\n\n")

                # --- heating_restart.in heredoc ---
                # heating.in uses irest=0/ntx=1 (fresh start from min2).
                # A restart resumes velocities from the checkpoint at target temp.
                f.write("# Restart variant of heating.in: resumes from checkpoint at target temp.\n")
                f.write("cat > heating_restart.in << 'HEATR'\n")
                f.write("Heating stage restart\n")
                f.write("&cntrl\n")
                f.write("  imin=0,\n")
                f.write("  irest=1,\n")
                f.write("  ntx=5,\n")
                f.write("  ntb=2,\n")
                f.write("  ntp=1,\n")
                f.write(f"  cut={gp.get('cutoff')},\n")
                f.write("  ntc=2,\n")
                f.write("  ntf=2,\n")
                f.write(f"  temp0={gp.get('target_temp')},\n")
                f.write(f"  tempi={gp.get('target_temp')},\n")
                f.write(f"  ntt={gp.get('thermostat')},\n")
                f.write(f"  gamma_ln={gp.get('collision_freq')},\n")
                f.write(f"  nstlim={hp.get('heating_steps')},\n")
                f.write("  dt=0.002,\n")
                f.write(f"  ntpr={hp.get('heating_output_write')},\n")
                f.write(f"  ntwx={hp.get('heating_trajectory_write')},\n")
                f.write(f"  ntwr={hp.get('heating_restart_write')},\n")
                if heat_ref:
                    f.write("  ntr=1,\n")
                    f.write(f"  restraint_wt={hp.get('heating_restraint_wt')},\n")
                    f.write(f"  restraintmask='{heat_ref}',\n")
                f.write("  iwrap=1,\n")
                f.write("/\n")
                f.write("END\n")
                f.write("HEATR\n\n")

                # --- Stage calls ---
                f.write('echo "--- min1 ---"\n')
                f.write(f'run_stage min1 min1.in min1.in complex.inpcrd {ref_file_min1} "" \\\n')
                f.write('    || { echo "FAILED: min1 — aborting"; exit 1; }\n')
                f.write('min1_out="${STAGE_OUT}"\n\n')

                f.write('echo "--- min2 ---"\n')
                f.write(f'run_stage min2 min2.in min2.in complex.inpcrd {ref_file_min2} "" \\\n')
                f.write('    || { echo "FAILED: min2 — aborting"; exit 1; }\n')
                f.write('min2_out="${STAGE_OUT}"\n\n')

                f.write('echo "--- heating ---"\n')
                f.write(f'run_stage heating heating.in heating_restart.in "${{min2_out}}" {ref_file_heat} "" \\\n')
                f.write('    || { echo "FAILED: heating — aborting"; exit 1; }\n')
                f.write('heating_out="${STAGE_OUT}"\n\n')

                f.write('echo "--- equilibration ---"\n')
                f.write(f'run_stage equilibration equilibration.in equilibration.in "${{heating_out}}" {ref_file_eq} yes \\\n')
                f.write('    || { echo "FAILED: equilibration — aborting"; exit 1; }\n')
                f.write('eq_out="${STAGE_OUT}"\n\n')

                f.write('echo "--- production ---"\n')
                f.write(f'run_stage production production.in production.in "${{eq_out}}" {ref_file_prod} yes \\\n')
                f.write('    || { echo "FAILED: production — aborting"; exit 1; }\n')
                f.write('echo "MD completed successfully → ${STAGE_OUT}"\n')

            os.chmod(execution_script_path, 0o755)
        except OSError as e:
            raise RuntimeError(f"Failed to write execution script '{execution_script_path}': {e}") from e
    
    def _start_md_simulation(self, md_assay_folder):

        # Query if the user wants to run in the foreground or background
        run_mode = input("\n⚙️  Do you want to run the MD simulation in the foreground or background? (fg/bg) [default: fg]: ").strip().lower() or 'fg'
        execution_script_path = os.path.join(md_assay_folder, "run_md.sh")
        try:
            import subprocess

            if run_mode in ['fg', 'foreground']:
                print("\n▶️  Starting MD simulation in the foreground...")
                subprocess.run([execution_script_path], check=True, cwd=md_assay_folder)
                return None

            elif run_mode in ['bg', 'background']:
                print("\n▶️  Starting MD simulation in the background...")
                process = subprocess.Popen([execution_script_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=md_assay_folder)
                print(f"   🆔 Process ID: {process.pid}")
                # No dedicated subprocess log file (output goes to per-stage
                # min/heating/equilibration/production .out files instead);
                # completion is judged from those via get_md_assay_status-style checks.
                return {'background': True, 'pid': process.pid, 'log_file': None}

            else:
                print("❌ Invalid option. Please choose 'fg' or 'bg'.")
                return None
        except Exception as e:
            print(f"❌ Error starting MD simulation: {e}")
            return None
            
    def perform_md_assay_on_receptor(self):
        """
        Configure one or more replica MD simulations starting from a registered
        receptor template (apo simulation, no ligand).  The user is asked how many
        replicas to create; each gets its own assay folder and registry entry.
        Auto-start is offered only for the first replica — remaining replicas
        require manual execution of their respective run scripts.
        """
        import traceback

        try:
            # Select receptor template from pdbs.db (shared across all replicas)
            template_dict = self._select_receptor_template()
            if template_dict is None:
                return

            # Select a MD method (shared across all replicas)
            md_parameters_dict = self._select_md_method()
            if md_parameters_dict is None:
                return

            # Query user for MD assay description (shared base description)
            assay_description = input("\n📝 Enter MD assay description (optional): ").strip()
            if not assay_description:
                assay_description = f"Receptor-only MD assay for template '{template_dict['pdb_template_name']}'"

            # Query number of replicas
            while True:
                try:
                    n_replicas_input = input("\n🔁 How many replicas do you want to create? [default: 1]: ").strip()
                    n_replicas = int(n_replicas_input) if n_replicas_input else 1
                    if n_replicas < 1:
                        print("❌ Number of replicas must be at least 1.")
                        continue
                    break
                except ValueError:
                    print("❌ Please enter a valid integer.")
                except KeyboardInterrupt:
                    print("\n❌ Cancelled.")
                    return

            # Optionally start the first replica automatically
            start_first = 'no'
            if n_replicas >= 1:
                start_first = input("\n▶️  Start replica 1 automatically after setup? (yes/no) [default: no]: ").strip().lower() or 'no'

            print(f"\n🔁 Creating {n_replicas} replica(s) for template '{template_dict['pdb_template_name']}'...")

            bg_info = None
            for replica_idx in range(1, n_replicas + 1):
                md_assay_folder = None
                md_assay_id = None
                replica_label = f"replica {replica_idx}/{n_replicas}"

                print(f"\n{'─' * 60}")
                print(f"  Setting up {replica_label}...")
                print(f"{'─' * 60}")

                try:
                    # Each replica gets its own assay folder and ID
                    md_assay_id, md_assay_folder = self._create_md_assay_folder()
                    if md_assay_id is None:
                        print(f"❌ Could not create assay folder for {replica_label}. Skipping.")
                        continue

                    # Append replica tag to description when more than one replica
                    replica_description = (
                        f"{assay_description} (replica {replica_idx})"
                        if n_replicas > 1
                        else assay_description
                    )

                    # Prepare receptor-only prmtop and inpcrd via tleap
                    print(f"⚙️  Preparing topology and coordinate files with tleap...")
                    if md_parameters_dict.get('method_type') == 'tleap_template':
                        prmtop_file, inpcrd_file = self._prepare_complex_prmtop_inpcrd_from_template(
                            selected_pose_pdb=None,
                            output_dir=md_assay_folder,
                            md_parameters_dict=md_parameters_dict,
                            receptor_pdb=template_dict['checked_pdb_path'],
                        )
                    else:
                        prmtop_file, inpcrd_file = self._prepare_receptor_only_prmtop_inpcrd_for_md(
                            template_dict['checked_pdb_path'], md_assay_folder, md_parameters_dict,
                            pdb_template_name=template_dict['pdb_template_name'],
                        )

                    # Prepare AMBER input files
                    print(f"⚙️  Writing AMBER simulation input files...")
                    self._prepare_md_simulation_input_files(md_parameters_dict, md_assay_folder, prmtop_file, inpcrd_file)

                    # Prepare execution shell script
                    print(f"⚙️  Writing execution script...")
                    self._prepare_md_execution_script(md_assay_folder, md_parameters_dict)

                    # Register the assay in md_assays table
                    self._create_md_assay_register_entry(
                        md_assay_id, md_assay_folder, replica_description,
                        receptor_template_name=template_dict['pdb_template_name'],
                        md_parameters_dict=md_parameters_dict,
                    )

                    print(f"✅ {replica_label.capitalize()} prepared in: {md_assay_folder}")

                    # Auto-start only for the first replica if the user requested it
                    if replica_idx == 1 and start_first in ['yes', 'y']:
                        bg_info = self._start_md_simulation(md_assay_folder)
                    else:
                        if replica_idx == 1:
                            print(f"ℹ️  Replica 1 not started. Run the execution script manually when ready.")
                        else:
                            print(f"ℹ️  {replica_label.capitalize()} requires manual start: {md_assay_folder}")

                except Exception as e:
                    print(f"\n❌ Error setting up {replica_label}: {e}")
                    traceback.print_exc()
                    if md_assay_folder and os.path.exists(md_assay_folder):
                        print(f"⚠️  Assay folder left on disk for inspection: {md_assay_folder}")

            print(f"\n{'═' * 60}")
            print(f"✅ Done — {n_replicas} replica(s) registered for template '{template_dict['pdb_template_name']}'.")
            print(f"{'═' * 60}")

            return bg_info

        except Exception as e:
            print(f"\n❌ Error performing receptor MD assay: {e}")
            traceback.print_exc()

    def _select_receptor_template(self):
        """
        Prompt the user to select a receptor template registered in pdb_templates
        (stored in docking/receptors/pdbs.db).

        Returns a dict with pdb_template_name, checked_pdb_path, template_folder_path,
        chains, and notes; or None if the user cancels or no templates exist.
        """
        import sqlite3

        pdbs_db_path = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')
        if not os.path.exists(pdbs_db_path):
            print(f"❌ No receptors database found at {pdbs_db_path}")
            return None

        try:
            conn = sqlite3.connect(pdbs_db_path)
            cursor = conn.cursor()
            cursor.execute("""
                SELECT pdb_id, pdb_template_name, checked_pdb_path,
                       template_folder_path, chains, notes
                FROM pdb_templates
                ORDER BY pdb_id ASC
            """)
            rows = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"❌ Error reading pdb_templates: {e}")
            return None

        if not rows:
            print("❌ No receptor templates found in pdb_templates.")
            return None

        print(f"\n🧬 RECEPTOR TEMPLATES IN PROJECT '{self.name}':")
        print("=" * 70)
        for i, (pdb_id, name, checked_pdb, folder, chains, notes) in enumerate(rows, 1):
            print(f"{i}. {name} (ID: {pdb_id})")
            print(f"   Chains: {chains}")
            print(f"   Checked PDB: {os.path.basename(checked_pdb) if checked_pdb else 'N/A'}")
            print(f"   Notes: {notes}")
            print("-" * 70)

        while True:
            try:
                selection = input("\n🔎 Select template by number (or 'cancel'): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Template selection cancelled.")
                    return None
                idx = int(selection) - 1
                if 0 <= idx < len(rows):
                    row = rows[idx]
                    print(f"✅ Selected template: '{row[1]}'")
                    return {
                        'pdb_id': row[0],
                        'pdb_template_name': row[1],
                        'checked_pdb_path': self._remap_project_path(row[2]),
                        'template_folder_path': self._remap_project_path(row[3]),
                        'chains': row[4],
                        'notes': row[5],
                    }
                else:
                    print(f"❌ Invalid selection. Enter a number between 1 and {len(rows)}.")
            except ValueError:
                print("❌ Please enter a valid number.")
            except KeyboardInterrupt:
                print("\n❌ Template selection cancelled.")
                return None

    def _prepare_receptor_only_prmtop_inpcrd_for_md(self, checked_pdb_path, output_dir, md_parameters_dict, pdb_template_name=None):
        """
        Build AMBER topology (prmtop) and coordinate (inpcrd) files for a receptor-only
        (apo) system using tleap.  Output files are named complex.prmtop / complex.inpcrd
        so that the shared execution script and input-file templates work unchanged.

        Raises RuntimeError on any failure so the caller can handle it with full context.
        tleap stdout+stderr are always saved to tleap.log inside the assay folder.
        """
        import subprocess

        tleap_params = md_parameters_dict.get('tleap_params', {})
        protein_ff = tleap_params.get('force_field', 'ff14SB')
        solvent_type = md_parameters_dict.get('solvent_type', 'explicit')

        prmtop_file = os.path.join(output_dir, 'complex.prmtop')
        inpcrd_file = os.path.join(output_dir, 'complex.inpcrd')
        tleap_in_file = os.path.join(output_dir, 'complex.in')
        tleap_log_file = os.path.join(output_dir, 'tleap.log')

        # Confirmed disulfide bridges (if any) for this receptor template
        disulfide_bonds = self._check_disulfides_in_template(pdb_template_name)
        disulfide_bond_lines = self._get_tleap_bond_lines(checked_pdb_path, disulfide_bonds, 'rec')

        try:
            with open(tleap_in_file, 'w') as f:
                f.write(f"source leaprc.protein.{protein_ff}\n")
                f.write(f"loadamberparams frcmod.ions1lm_126_tip3p\n")
                if solvent_type == 'explicit':
                    water_model = tleap_params.get('water_model', 'tip3p')
                    f.write(f"source leaprc.water.{water_model}\n")
                    f.write("HOH = WAT\n")
                f.write(f"rec = loadpdb {checked_pdb_path}\n")
                for bond_line in disulfide_bond_lines:
                    f.write(bond_line)
                if solvent_type == 'explicit':
                    water_type = tleap_params.get('water_type', 'TIP3PBOX')
                    box_type = tleap_params.get('box_type', 'Box')
                    box_distance_nm = tleap_params.get('box_distance_nm', 10.0)
                    wat_closeness = tleap_params.get('wat_closeness', 1.0)
                    f.write(f"solvate{box_type} rec {water_type} {box_distance_nm} {wat_closeness}\n")
                    add_ions = tleap_params.get('add_ions', 'yes')
                    if add_ions in ['yes', 'y']:
                        positive_ion = tleap_params.get('positive_ion', 'Na+')
                        negative_ion = tleap_params.get('negative_ion', 'Cl-')
                        f.write(f"addions rec {positive_ion} 0\n")
                        f.write(f"addions rec {negative_ion} 0\n")
                f.write(f"saveamberparm rec {prmtop_file} {inpcrd_file}\n")
                f.write("quit\n")
        except OSError as e:
            raise RuntimeError(f"Failed to write tleap input file '{tleap_in_file}': {e}") from e

        try:
            result = subprocess.run(
                f"tleap -f {tleap_in_file}",
                shell=True,
                capture_output=True,
                text=True,
            )
        except subprocess.SubprocessError as e:
            raise RuntimeError(f"Failed to execute tleap: {e}") from e

        # Always save tleap output for later inspection
        try:
            with open(tleap_log_file, 'w') as lf:
                lf.write(result.stdout)
                if result.stderr:
                    lf.write("\n--- STDERR ---\n")
                    lf.write(result.stderr)
        except OSError:
            pass

        if result.returncode != 0:
            tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
            print(f"   tleap output (last lines):\n{tail}")
            raise RuntimeError(
                f"tleap exited with code {result.returncode}. "
                f"Full log: {tleap_log_file}"
            )

        # tleap exits 0 even when it fails to write output files — verify explicitly
        missing = [f for f in (prmtop_file, inpcrd_file) if not os.path.exists(f)]
        if missing:
            tail = result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout
            print(f"   tleap output (last lines):\n{tail}")
            raise RuntimeError(
                f"tleap completed but did not produce: {', '.join(os.path.basename(f) for f in missing)}. "
                f"Full log: {tleap_log_file}"
            )

        print(f"   ✓ tleap completed successfully (log: {tleap_log_file})")
        return prmtop_file, inpcrd_file

    def compute_mmgbsa_on_trajectory(self):
        """
        Compute MM-GBSA binding free energy on completed production MD trajectories
        using AMBER's MMPBSA.py (single-trajectory approach).

        Workflow:
          1. Select a completed ligand-receptor MD assay
          2. Collect MM-GBSA parameters interactively
          3. Run ante-MMPBSA.py to generate gas-phase receptor/ligand/complex topologies
          4. Write the MMPBSA.py input namelist
          5. Write run_mmgbsa.sh execution script
          6. Query user: run now (fg/bg) or leave for manual execution
          7. If fg: parse, store, and display results immediately after completion
        """
        import traceback

        try:
            print("\n🔬 COMPUTE MM-GBSA ON MD TRAJECTORY")
            print("=" * 60)

            assay_info = self._select_completed_md_assay_for_mmgbsa()
            if assay_info is None:
                return

            assay_folder = assay_info['assay_folder_path']
            assay_id = assay_info['assay_id']
            ligand_name = assay_info.get('ligand_name', '')

            prmtop_file = os.path.join(assay_folder, 'complex.prmtop')
            trajectory_file = os.path.join(assay_folder, 'production.nc')
            for req_file in [prmtop_file, trajectory_file]:
                if not os.path.exists(req_file):
                    print(f"❌ Required file not found: {req_file}")
                    return

            mmgbsa_params = self._collect_mmgbsa_parameters(ligand_name)
            if mmgbsa_params is None:
                return

            mmgbsa_folder = os.path.join(assay_folder, 'mmgbsa')
            if os.path.exists(mmgbsa_folder):
                print(f"\n⚠️  MM-GBSA folder already exists: {mmgbsa_folder}")
                confirm = input("🗑️  Delete existing folder and continue? (yes/no) [default: no]: ").strip().lower() or 'no'
                if confirm not in ['yes', 'y']:
                    print("❌ MM-GBSA computation cancelled.")
                    return
                import shutil as _shutil
                _shutil.rmtree(mmgbsa_folder)
            os.makedirs(mmgbsa_folder)
            print(f"\n📂 MM-GBSA output folder: {mmgbsa_folder}")

            print(f"\n⚙️  Running ante-MMPBSA.py to generate gas-phase topologies...")
            com_prmtop, rec_prmtop, lig_prmtop = self._run_ante_mmpbsa(
                prmtop_file, mmgbsa_folder, mmgbsa_params
            )
            if com_prmtop is None:
                return

            print(f"\n⚙️  Writing MM-GBSA input file...")
            mmgbsa_in_file = self._write_mmgbsa_input_file(mmgbsa_folder, mmgbsa_params)

            print(f"\n⚙️  Writing MM-GBSA execution script...")
            script_path = self._prepare_mmgbsa_execution_script(
                mmgbsa_folder, com_prmtop, rec_prmtop, lig_prmtop,
                trajectory_file, mmgbsa_in_file,
                prmtop_file, mmgbsa_params
            )
            print(f"✅ MM-GBSA assay prepared successfully.")
            print(f"   Script: {script_path}")

            start_now = input(
                "\n▶️  Do you want to start the MM-GBSA computation now? (yes/no) [default: no]: "
            ).strip().lower() or 'no'
            if start_now in ['yes', 'y']:
                self._start_mmgbsa_computation(
                    mmgbsa_folder, script_path, assay_id, mmgbsa_params, assay_info
                )
            else:
                print(f"ℹ️  MM-GBSA computation not started.")
                print(f"   Run the script manually when ready: {script_path}")

        except Exception as e:
            print(f"\n❌ Error computing MM-GBSA: {e}")
            traceback.print_exc()

    def _select_completed_md_assay_for_mmgbsa(self):
        """
        List completed ligand-receptor MD assays and prompt the user to pick one.
        Completion is determined by the presence of 'Total wall time:' or
        'FINAL RESULTS' in production.out.
        Returns a dict with assay metadata or None if the user cancels.
        """
        import sqlite3

        try:
            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()

            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_assays';")
            if not cursor.fetchone():
                print("❌ No MD assays found (md_assays table does not exist).")
                conn.close()
                return None

            cursor.execute("""
                SELECT assay_id, md_assay, description, assay_folder_path,
                       ligand_name, pose_id
                FROM md_assays
                WHERE ligand_name IS NOT NULL
                ORDER BY assay_id ASC
            """)
            rows = cursor.fetchall()
            conn.close()

            if not rows:
                print("❌ No ligand-receptor MD assays found.")
                return None

            completed = []
            for assay_id, md_assay, description, folder, ligname, pose_id in rows:
                folder = self._remap_project_path(folder)
                if not folder or not os.path.exists(folder):
                    continue
                prod_out = os.path.join(folder, 'production.out')
                if not os.path.exists(prod_out):
                    continue
                try:
                    with open(prod_out, 'r') as fh:
                        content = fh.read()
                except OSError:
                    continue
                if 'Total wall time:' not in content and 'FINAL RESULTS' not in content:
                    continue
                completed.append({
                    'assay_id': assay_id,
                    'md_assay': md_assay,
                    'description': description,
                    'assay_folder_path': folder,
                    'ligand_name': ligname,
                    'pose_id': pose_id,
                })

            if not completed:
                print("❌ No completed ligand-receptor MD assays found.")
                print("   (production.out must contain 'Total wall time:' or 'FINAL RESULTS')")
                return None

            print(f"\n🧬 COMPLETED LIGAND-RECEPTOR MD ASSAYS:")
            print("=" * 70)
            for a in completed:
                print(f"  Assay ID    : {a['assay_id']}")
                print(f"  Name        : {a['md_assay']}")
                print(f"  Ligand      : {a['ligand_name']}  (pose {a['pose_id']})")
                print(f"  Description : {a['description']}")
                print(f"  Folder      : {a['assay_folder_path']}")
                print("-" * 70)

            assay_ids = [str(a['assay_id']) for a in completed]
            while True:
                selection = input("\n🔎 Enter the Assay ID to compute MM-GBSA (or 'cancel'): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    return None
                if selection in assay_ids:
                    selected = next(a for a in completed if str(a['assay_id']) == selection)
                    print(f"✅ Selected assay: {selected['md_assay']} (ID {selected['assay_id']})")
                    return selected
                print("❌ Invalid Assay ID. Please try again.")

        except Exception as e:
            print(f"❌ Error selecting MD assay: {e}")
            return None

    def _collect_mmgbsa_parameters(self, ligand_name):
        """
        Interactively collect MM-GBSA parameters from the user.
        Returns a dict of parameters, or None if the user cancels.
        """
        try:
            print(f"\n⚙️  MM-GBSA CONFIGURATION")
            print("-" * 60)

            params = {}

            # Ligand residue mask in the AMBER topology.
            # antechamber preserves the residue name from the input PDB (usually 'LIG'
            # for AutoDock-derived poses or whatever moldf wrote).  Ask the user to
            # confirm, defaulting to :LIG.
            default_lig_mask = ':UNL'
            print(f"\nThe ligand residue name in the AMBER topology is set by antechamber")
            print(f"from the input PDB residue name (commonly 'UNL' for docked poses).")
            lig_mask = input(
                "Ligand residue mask (e.g. :UNL, :LIG) [default: " + default_lig_mask + "]: "
            ).strip() or default_lig_mask
            params['ligand_mask'] = lig_mask

            # Residues to strip — must include water and ions but NOT the ligand.
            # Every residue name needs its own : prefix so ante-MMPBSA.py identifies
            # each as a residue selector (e.g. :WAT,:Na+,:Cl-).
            strip_mask = input(
                "Mask for residues to strip (solvent/ions) [default: :WAT,:Na+,:Cl-]: "
            ).strip() or ':WAT,:Na+,:Cl-'
            params['strip_mask'] = strip_mask

            # GB model
            print("\nGB model options:")
            print("  2 = Onufriev et al. (2000) model I")
            print("  5 = Onufriev et al. (2000) model II  (recommended for proteins)")
            print("  7 = Mongan et al. (2007)")
            print("  8 = Nguyen et al. (2013)")
            igb_str = input("GB model (igb) [default: 5]: ").strip() or '5'
            try:
                igb = int(igb_str)
                if igb not in [1, 2, 5, 7, 8]:
                    print("⚠️  Unsupported igb value, defaulting to 5")
                    igb = 5
            except ValueError:
                igb = 5
            params['igb'] = igb

            # Salt concentration
            saltcon_str = input("\nSalt concentration (M) [default: 0.15]: ").strip() or '0.15'
            try:
                params['saltcon'] = float(saltcon_str)
            except ValueError:
                params['saltcon'] = 0.15

            # Frame selection
            print("\nTrajectory frame selection:")
            try:
                params['startframe'] = int(input("  Start frame [default: 1]: ").strip() or '1')
                params['endframe'] = int(input("  End frame [default: 9999]: ").strip() or '9999')
                params['interval'] = int(input("  Interval (every Nth frame) [default: 1]: ").strip() or '1')
            except ValueError:
                params['startframe'] = 1
                params['endframe'] = 9999
                params['interval'] = 1

            # Whether to keep MMPBSA.py temporary files
            keep_str = input(
                "\nKeep MMPBSA.py intermediate files? (yes/no) [default: no]: "
            ).strip().lower() or 'no'
            params['keep_files'] = keep_str in ['yes', 'y']

            print(f"\n✅ MM-GBSA parameters configured:")
            print(f"   Ligand mask : {params['ligand_mask']}")
            print(f"   Strip mask  : {params['strip_mask']}")
            print(f"   GB model    : igb={params['igb']}")
            print(f"   Salt conc.  : {params['saltcon']} M")
            print(f"   Frames      : {params['startframe']} to {params['endframe']}, every {params['interval']}")
            return params

        except KeyboardInterrupt:
            print("\n❌ MM-GBSA parameter collection cancelled.")
            return None
        except Exception as e:
            print(f"❌ Error collecting MM-GBSA parameters: {e}")
            return None

    def _run_ante_mmpbsa(self, prmtop_file, mmgbsa_folder, mmgbsa_params):
        """
        Run ante-MMPBSA.py to produce gas-phase topology files for the complex,
        receptor, and ligand (solvent + ions stripped).
        Returns (com_prmtop, rec_prmtop, lig_prmtop) paths, or (None, None, None)
        on failure.
        """
        import subprocess
        import shutil

        amber_bin = os.path.expanduser('~/Programas/Amber26/ambertools26/bin')
        ante_mmpbsa = (
            shutil.which('ante-MMPBSA.py')
            or os.path.join(amber_bin, 'ante-MMPBSA.py')
            or 'ante-MMPBSA.py'
        )

        com_prmtop = os.path.join(mmgbsa_folder, 'com.prmtop')
        rec_prmtop = os.path.join(mmgbsa_folder, 'rec.prmtop')
        lig_prmtop = os.path.join(mmgbsa_folder, 'lig.prmtop')
        ante_log   = os.path.join(mmgbsa_folder, 'ante_mmpbsa.log')

        cmd = [
            ante_mmpbsa,
            '-p', prmtop_file,
            '-c', com_prmtop,
            '-r', rec_prmtop,
            '-l', lig_prmtop,
            '-s', mmgbsa_params['strip_mask'],
            '-n', mmgbsa_params['ligand_mask'],
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=mmgbsa_folder)
        except subprocess.SubprocessError as e:
            print(f"❌ Failed to execute ante-MMPBSA.py: {e}")
            return None, None, None

        try:
            with open(ante_log, 'w') as lf:
                lf.write(result.stdout)
                if result.stderr:
                    lf.write("\n--- STDERR ---\n")
                    lf.write(result.stderr)
        except OSError:
            pass

        if result.returncode != 0:
            print(f"❌ ante-MMPBSA.py failed (exit code {result.returncode})")
            print(f"   Log: {ante_log}")
            tail = (result.stdout + result.stderr)[-2000:]
            print(f"   Output:\n{tail}")
            return None, None, None

        missing = [f for f in [com_prmtop, rec_prmtop, lig_prmtop] if not os.path.exists(f)]
        if missing:
            print(f"❌ ante-MMPBSA.py did not produce: {', '.join(os.path.basename(f) for f in missing)}")
            print(f"   Log: {ante_log}")
            return None, None, None

        print(f"   ✓ ante-MMPBSA.py completed (log: {ante_log})")

        # Report atom counts and sanity-check the topology split.
        # com.prmtop must have more atoms than rec.prmtop and lig.prmtop individually,
        # and ideally com = rec + lig (single-trajectory approach).
        # Identical counts mean the ligand mask did not match any residue.
        n_com = self._get_prmtop_natom(com_prmtop)
        n_rec = self._get_prmtop_natom(rec_prmtop)
        n_lig = self._get_prmtop_natom(lig_prmtop)
        print(f"\n   Topology atom counts:")
        print(f"     com.prmtop (complex, gas phase) : {n_com if n_com is not None else 'unknown'} atoms")
        print(f"     rec.prmtop (receptor only)       : {n_rec if n_rec is not None else 'unknown'} atoms")
        print(f"     lig.prmtop (ligand only)         : {n_lig if n_lig is not None else 'unknown'} atoms")

        if None not in (n_com, n_rec, n_lig):
            if n_com == n_rec or n_com == n_lig or n_rec == n_lig:
                print(f"\n   ⚠️  WARNING: two or more topology files share the same atom count.")
                print(f"   This almost always means the ligand mask '{mmgbsa_params['ligand_mask']}'")
                print(f"   did not match any residue in the topology, so ante-MMPBSA.py could")
                print(f"   not separate receptor from ligand.")
                print(f"   Fix: inspect the residue names in the complex topology with")
                print(f"     cpptraj {prmtop_file}")
                print(f"     parminfo :*")
                print(f"   then re-run with the correct mask (e.g. ':LIG', ':UNL', ':MOL').")
                return None, None, None
            elif n_com == n_rec + n_lig:
                print(f"   ✓ Atom counts are consistent: {n_com} = {n_rec} (rec) + {n_lig} (lig)")
            else:
                print(f"\n   ⚠️  Atom counts do not satisfy com = rec + lig "
                      f"({n_com} ≠ {n_rec} + {n_lig}).")
                print(f"   This can happen when the strip mask also removes atoms that")
                print(f"   ante-MMPBSA.py needs to account for (e.g. crystal waters kept")
                print(f"   in the receptor).  Check {ante_log} for details.")

        return com_prmtop, rec_prmtop, lig_prmtop

    def _get_prmtop_natom(self, prmtop_file):
        """Read NATOM from the %FLAG NATOM section of an AMBER prmtop file."""
        try:
            with open(prmtop_file, 'r') as fh:
                lines = fh.readlines()
            for i, line in enumerate(lines):
                if '%FLAG NATOM' in line:
                    for j in range(i + 1, min(i + 4, len(lines))):
                        if lines[j].startswith('%FORMAT'):
                            continue
                        stripped = lines[j].strip()
                        if stripped:
                            return int(stripped.split()[0])
        except Exception:
            pass
        return None

    def _write_mmgbsa_input_file(self, mmgbsa_folder, mmgbsa_params):
        """Write the MMPBSA.py input namelist file and return its path."""
        mmgbsa_in = os.path.join(mmgbsa_folder, 'mmgbsa.in')
        with open(mmgbsa_in, 'w') as f:
            f.write("MM-GBSA calculation\n")
            f.write("&general\n")
            f.write(f"  startframe={mmgbsa_params['startframe']},\n")
            f.write(f"  endframe={mmgbsa_params['endframe']},\n")
            f.write(f"  interval={mmgbsa_params['interval']},\n")
            f.write("  verbose=2,\n")
            f.write("  netcdf=1,\n")
            f.write("/\n")
            f.write("&gb\n")
            f.write(f"  igb={mmgbsa_params['igb']},\n")
            f.write(f"  saltcon={mmgbsa_params['saltcon']},\n")
            f.write("/\n")
        print(f"   ✓ MM-GBSA input file written: {os.path.basename(mmgbsa_in)}")
        return mmgbsa_in

    def _prepare_mmgbsa_execution_script(self, mmgbsa_folder, com_prmtop, rec_prmtop,
                                           lig_prmtop, traj_file, mmgbsa_in,
                                           solvated_prmtop, mmgbsa_params):
        """
        Write run_mmgbsa.sh inside mmgbsa_folder and return its path.

        The script is fully self-contained and reproduces the complete pipeline:
          Step 1 — ante-MMPBSA.py: strips solvent/ions and splits the solvated
                   complex topology into gas-phase complex, receptor, and ligand
                   topology files.
          Step 2 — MMPBSA.py: computes MM-GBSA energies frame by frame and
                   writes the results summary and per-frame CSV.

        Including ante-MMPBSA.py in the script means the script can be re-run
        from scratch without needing any Python setup step.
        """
        import shutil as _shutil

        amber_bin = os.path.expanduser('~/Programas/Amber26/ambertools26/bin')
        amber_path_export = f"export PATH={amber_bin}:$PATH"

        ante_bin = _shutil.which('ante-MMPBSA.py') or os.path.join(amber_bin, 'ante-MMPBSA.py')
        mmpbsa_bin = _shutil.which('MMPBSA.py') or os.path.join(amber_bin, 'MMPBSA.py')

        script_path      = os.path.join(mmgbsa_folder, 'run_mmgbsa.sh')
        ante_log         = os.path.join(mmgbsa_folder, 'ante_mmpbsa.log')
        cpptraj_in       = os.path.join(mmgbsa_folder, 'strip_traj.in')
        cpptraj_log      = os.path.join(mmgbsa_folder, 'strip_traj.log')
        stripped_traj    = os.path.join(mmgbsa_folder, 'production_stripped.nc')
        results_dat      = os.path.join(mmgbsa_folder, 'mmgbsa_results.dat')
        energies_csv     = os.path.join(mmgbsa_folder, 'mmgbsa_energies.csv')
        mmpbsa_log       = os.path.join(mmgbsa_folder, 'mmpbsa.log')

        try:
            with open(script_path, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write("set -uo pipefail\n")
                f.write('cd "$(dirname "$0")"\n')
                f.write(f"{amber_path_export}\n\n")

                # --- Step 1: ante-MMPBSA.py ---
                strip_mask = mmgbsa_params['strip_mask']
                lig_mask   = mmgbsa_params['ligand_mask']
                f.write('echo "--- Step 1: ante-MMPBSA.py (gas-phase topology split) ---"\n')
                f.write("ante-MMPBSA.py \\\n")
                f.write(f"    -p {solvated_prmtop} \\\n")
                f.write(f"    -c {com_prmtop} \\\n")
                f.write(f"    -r {rec_prmtop} \\\n")
                f.write(f"    -l {lig_prmtop} \\\n")
                f.write(f"    -s {strip_mask} \\\n")
                f.write(f"    -n {lig_mask} \\\n")
                f.write(f"    > {ante_log} 2>&1\n\n")
                f.write('echo "   ante-MMPBSA.py done"\n\n')

                # --- Step 2: cpptraj — autoimage, center, and strip trajectory ---
                # The production trajectory contains the full solvated system; it must
                # be stripped with the same mask used for the topology so that frame
                # atom counts match com.prmtop when MMPBSA.py reads them.
                # autoimage fixes molecules broken across periodic boundaries.
                # center translates the solute (everything that is NOT stripped) to the
                # origin by center of mass before stripping the solvent.
                solute_mask = '!(' + strip_mask + ')'
                f.write('echo "--- Step 2: cpptraj (autoimage, center, strip trajectory) ---"\n')
                f.write(f"cat > {cpptraj_in} << 'CPPTRAJ_EOF'\n")
                f.write(f"parm {solvated_prmtop}\n")
                f.write(f"trajin {traj_file}\n")
                f.write("autoimage\n")
                f.write(f"center {solute_mask} mass origin\n")
                f.write("image origin center familiar\n")
                f.write(f"strip {strip_mask}\n")
                f.write(f"trajout {stripped_traj} netcdf\n")
                f.write("run\n")
                f.write("quit\n")
                f.write("CPPTRAJ_EOF\n\n")
                f.write(f"cpptraj -i {cpptraj_in} > {cpptraj_log} 2>&1\n")
                f.write('echo "   cpptraj trajectory processing done"\n\n')

                # --- Step 3: MMPBSA.py ---
                f.write('echo "--- Step 3: MMPBSA.py (MM-GBSA calculation) ---"\n')
                f.write(f"MMPBSA.py -O \\\n")
                f.write(f"    -i  {mmgbsa_in} \\\n")
                f.write(f"    -o  {results_dat} \\\n")
                f.write(f"    -eo {energies_csv} \\\n")
                f.write(f"    -cp {com_prmtop} \\\n")
                f.write(f"    -rp {rec_prmtop} \\\n")
                f.write(f"    -lp {lig_prmtop} \\\n")
                f.write(f"    -y  {stripped_traj} \\\n")
                f.write(f"    2>&1 | tee {mmpbsa_log}\n\n")
                f.write('echo "--- MM-GBSA computation finished ---"\n')

            os.chmod(script_path, 0o755)
        except OSError as e:
            raise RuntimeError(f"Failed to write MM-GBSA execution script: {e}") from e

        print(f"   ✓ Execution script written: {os.path.basename(script_path)}")
        return script_path

    def _start_mmgbsa_computation(self, mmgbsa_folder, script_path, assay_id,
                                   mmgbsa_params, assay_info):
        """
        Ask the user whether to run the MM-GBSA script in the foreground or background,
        then launch it accordingly.

        Foreground: blocks until completion, then parses mmgbsa_results.dat and stores
        results in the database.
        Background: launches the script detached; results must be collected manually.
        """
        import subprocess

        run_mode = input(
            "\n⚙️  Run in foreground or background? (fg/bg) [default: fg]: "
        ).strip().lower() or 'fg'

        results_dat  = os.path.join(mmgbsa_folder, 'mmgbsa_results.dat')
        energies_csv = os.path.join(mmgbsa_folder, 'mmgbsa_energies.csv')

        if run_mode in ['fg', 'foreground']:
            print("\n▶️  Starting MM-GBSA computation in the foreground...")
            try:
                proc = subprocess.run([script_path], cwd=mmgbsa_folder)
            except Exception as e:
                print(f"❌ Error running MM-GBSA script: {e}")
                return

            if proc.returncode != 0:
                print(f"❌ MM-GBSA script exited with code {proc.returncode}.")
                mmpbsa_log = os.path.join(mmgbsa_folder, 'mmpbsa.log')
                print(f"   Inspect the log for details: {mmpbsa_log}")
                return

            if not os.path.exists(results_dat):
                print(f"❌ MM-GBSA results file not found after run: {results_dat}")
                return

            # Parse, store, and display results only available after fg completion
            parsed = self._parse_mmpbsa_output(results_dat)
            if parsed is None:
                print(f"⚠️  Could not parse MMPBSA.py output — raw results in: {results_dat}")
                return

            parsed['results_file'] = results_dat
            parsed['energies_file'] = energies_csv if os.path.exists(energies_csv) else None
            self._store_mmgbsa_results(assay_id, parsed, mmgbsa_params, mmgbsa_folder)
            self._display_mmgbsa_results(parsed, assay_info)

        elif run_mode in ['bg', 'background']:
            print("\n▶️  Starting MM-GBSA computation in the background...")
            try:
                subprocess.Popen(
                    [script_path],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    cwd=mmgbsa_folder,
                )
            except Exception as e:
                print(f"❌ Error launching MM-GBSA script: {e}")
                return
            print(f"ℹ️  MM-GBSA running in background.")
            print(f"   Results will be written to: {results_dat}")
            print(f"   Log: {os.path.join(mmgbsa_folder, 'mmpbsa.log')}")
            print(f"   Re-run compute_mmgbsa_on_trajectory() after completion to parse and store results.")

        else:
            print("❌ Invalid option. Please choose 'fg' or 'bg'.")

    def _parse_mmpbsa_output(self, results_file):
        """
        Parse the MMPBSA.py summary output file.
        Returns a dict mapping energy component names to {average, std_dev} dicts,
        plus a 'DELTA_G_binding' key for the total binding free energy.
        Returns None if parsing fails entirely.
        """
        import re

        results = {}
        try:
            with open(results_file, 'r') as fh:
                content = fh.read()

            # Each line in the DELTA section looks like:
            #   COMPONENT    AVERAGE    STD_DEV    [STD_ERR_OF_MEAN]
            component_pattern = re.compile(
                r'^\s*(VDWAALS|EEL|EGB|ESURF|EPOL|ENPOLAR|EDISPER|EPB|ENPB|ECAVITY'
                r'|DELTA G gas|DELTA G solv|DELTA TOTAL)\s+'
                r'(-?\d+\.\d+)\s+(\d+\.\d+)',
                re.MULTILINE,
            )
            for m in component_pattern.finditer(content):
                results[m.group(1).strip()] = {
                    'average': float(m.group(2)),
                    'std_dev': float(m.group(3)),
                }

            # Some MMPBSA.py versions emit a dedicated "DELTA G binding" line.
            binding_match = re.search(
                r'DELTA G binding\s*=\s*(-?\d+\.\d+)\s+\+/-\s+(\d+\.\d+)',
                content,
            )
            if binding_match:
                results['DELTA_G_binding'] = {
                    'average': float(binding_match.group(1)),
                    'std_dev': float(binding_match.group(2)),
                }
            elif 'DELTA TOTAL' in results:
                # In single-trajectory MMGBSA, DELTA TOTAL IS the binding free energy.
                results['DELTA_G_binding'] = results['DELTA TOTAL']

            return results if results else None

        except Exception as e:
            print(f"⚠️  Error parsing MMPBSA.py output: {e}")
            return None

    def _store_mmgbsa_results(self, assay_id, results, mmgbsa_params, mmgbsa_folder):
        """
        Store MM-GBSA results as a JSON blob in the mmgbsa_results column of md_assays.
        The column is created if it does not yet exist.
        """
        import sqlite3
        import json

        try:
            payload = {
                'parameters': mmgbsa_params,
                'results': results,
                'mmgbsa_folder': mmgbsa_folder,
            }
            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()
            try:
                cursor.execute("ALTER TABLE md_assays ADD COLUMN mmgbsa_results TEXT")
            except Exception:
                pass  # column already exists
            cursor.execute(
                "UPDATE md_assays SET mmgbsa_results = ? WHERE assay_id = ?",
                (json.dumps(payload, indent=2), assay_id),
            )
            conn.commit()
            conn.close()
            print(f"   ✓ MM-GBSA results stored in database (assay_id={assay_id})")
        except Exception as e:
            print(f"⚠️  Error storing MM-GBSA results in database: {e}")

    def _display_mmgbsa_results(self, results, assay_info):
        """Print a formatted MM-GBSA results summary to the terminal."""
        print(f"\n{'=' * 60}")
        print(f"🔬 MM-GBSA RESULTS SUMMARY")
        print(f"{'=' * 60}")
        print(f"  Assay  : {assay_info.get('md_assay', 'N/A')} (ID {assay_info.get('assay_id', 'N/A')})")
        print(f"  Ligand : {assay_info.get('ligand_name', 'N/A')} (pose {assay_info.get('pose_id', 'N/A')})")
        print()

        for comp in ('VDWAALS', 'EEL', 'EGB', 'ESURF', 'DELTA G gas', 'DELTA G solv', 'DELTA TOTAL'):
            if comp in results:
                avg = results[comp]['average']
                std = results[comp]['std_dev']
                print(f"  {comp:<18} : {avg:>10.2f} ± {std:>7.2f} kcal/mol")

        binding = results.get('DELTA_G_binding')
        if binding:
            print()
            print(f"  {'─' * 48}")
            avg = binding['average']
            std = binding['std_dev']
            print(f"  {'ΔG binding':<18} : {avg:>10.2f} ± {std:>7.2f} kcal/mol")
            print(f"  {'─' * 48}")

        print()
        results_file = results.get('results_file', '')
        if results_file:
            print(f"  Full results   : {results_file}")
        energies_file = results.get('energies_file', '')
        if energies_file:
            print(f"  Per-frame data : {energies_file}")
        print(f"{'=' * 60}")

    def list_md_assays(self):
        """List all MD assays registered in the project."""
        try:
            import sqlite3

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()

            # Check if md_assays table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_assays';")
            if not cursor.fetchone():
                print("❌ No MD assays found. The 'md_assays' table does not exist.")
                return

            # Fetch all MD assays
            cursor.execute("SELECT * FROM md_assays;")
            assays = cursor.fetchall()

            if not assays:
                print("❌ No MD assays found.")
                return

            print("\n📋 MD Assays:")
            for assay in assays:
                print(f"  - {assay[1]} (ID: {assay[0]})")
                print(f"    Description: {assay[2]}")
                print(f"    Docking Assay ID: {assay[3]}")
                print(f"    Ligand Name: {assay[4]}")
                print(f"    Pose ID: {assay[5]}")
                print("-" * 50)

            # Select an assay to show details
            assay_ids = [str(a[0]) for a in assays]
            while True:
                selection = input("\n🔎 Enter the MD Assay ID to show details (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    return
                if selection in assay_ids:
                    selected_assay = next((a for a in assays if str(a[0]) == selection), None)
                    if selected_assay:
                        print(f"\n📋 Details for MD Assay '{selected_assay[1]}' (ID: {selected_assay[0]}):")
                        print(f"  Description: {selected_assay[2]}")    
                        print(f"  Assay folder path: {selected_assay[3]}")
                        print(f"  Docking Assay ID: {selected_assay[4]}")
                        print(f"  Ligand Name: {selected_assay[5]}")
                        print(f"  Pose ID: {selected_assay[6]}")
                        print(f"  MD Parameters: {selected_assay[7]}")

        except Exception as e:
            print(f"❌ Error listing MD assays: {e}")

    def delete_md_assay(self):
        """Delete an MD assay as created by perform_md_assay(), removing both its
        assay folder on disk and its registry entry in the md_assays table."""
        try:
            import sqlite3
            import shutil

            conn = sqlite3.connect(self.__md_registers_db)
            cursor = conn.cursor()

            # Check if md_assays table exists
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='md_assays';")
            if not cursor.fetchone():
                print("❌ No MD assays found. The 'md_assays' table does not exist.")
                conn.close()
                return

            # Fetch all MD assays
            cursor.execute("SELECT assay_id, md_assay, description, assay_folder_path FROM md_assays;")
            assays = [
                (assay_id, md_assay, description, self._remap_project_path(assay_folder_path))
                for assay_id, md_assay, description, assay_folder_path in cursor.fetchall()
            ]

            if not assays:
                print("❌ No MD assays found.")
                conn.close()
                return

            print(f"\n🧬 MOLECULAR DYNAMICS ASSAYS IN PROJECT '{self.name}':")
            print("=" * 70)
            for assay_id, md_assay, description, assay_folder_path in assays:
                print(f"🏷️  Assay ID: {assay_id}")
                print(f"🏷️  Assay Name: {md_assay}")
                print(f"📝 Description: {description}")
                print(f"📁 Folder: {assay_folder_path}")
                print("-" * 70)

            # Prompt user to select an assay to delete
            assay_ids = [str(a[0]) for a in assays]
            while True:
                selection = input("\n🗑️  Enter the Assay ID to delete (or 'cancel' to exit): ").strip()
                if selection.lower() in ['cancel', 'quit', 'exit']:
                    print("❌ Operation cancelled.")
                    conn.close()
                    return
                if selection in assay_ids:
                    selected_assay = next(a for a in assays if str(a[0]) == selection)
                    _, md_assay_name, _, assay_folder_path = selected_assay
                    confirm = input(
                        f"⚠️  Are you sure you want to delete Assay ID {selection} "
                        f"('{md_assay_name}') and its folder '{assay_folder_path}'? (yes/no): "
                    ).strip().lower()
                    if confirm in ['yes', 'y']:
                        if assay_folder_path and os.path.isdir(assay_folder_path):
                            try:
                                shutil.rmtree(assay_folder_path)
                                print(f"✅ Deleted assay folder: {assay_folder_path}")
                            except Exception as e:
                                print(f"⚠️  Could not delete assay folder '{assay_folder_path}': {e}")
                        else:
                            print(f"⚠️  Assay folder not found on disk, skipping: {assay_folder_path}")

                        cursor.execute("DELETE FROM md_assays WHERE assay_id = ?", (selection,))
                        conn.commit()
                        print(f"✅ Assay ID {selection} deleted from the registry.")
                    else:
                        print("❌ Deletion cancelled.")
                    break
                else:
                    print("❌ Invalid Assay ID. Please try again.")

            conn.close()

        except Exception as e:
            print(f"❌ Error deleting MD assay: {e}")
