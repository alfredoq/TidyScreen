###This section will import all the required modules
###DO NOT MODIFY THESE IMPORTS.
import sys
from pprint import pp as pp
sys.path.append('/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package')
##
## Imports related to the exploration of chemical space
from project_structure import generate_project_structure as gen_proj
import chemspace.actions_level_0 as acts_level_0 
import chemspace.actions_level_1 as acts_level_1 
import chemspace.general_procedures.operations_on_sql_database as db_ops
from  chemspace.general_procedures import smiles_processing as smiles_processing
from  chemspace.general_procedures import depiction_procedures as depic
from  chemspace.general_procedures import objects_processing as obj_proc
from  chemspace.general_procedures import blob_operations as blob_ops
from  chemspace.input_output_operations import actions_on_dbs as db_actions
from  chemspace.input_output_operations import i_o_checkings as i_o_chk
from  chemspace.input_output_operations import read_data_from_db as read_data
##
## Imports related to the preparation and execution of molecular docking assays
from docking.registers_management import receptor_models_registers as rec_regs
from docking.registers_management import docking_params_registers as reg_dock_parm
from docking.perform_docking import docking_assay as dock_ass
from docking.receptor_management import receptor_management as rec_mng
from docking.perform_docking import docking_configuration as dock_cfg
from docking.perform_docking import docking_assay as dock_assay
from docking.ligands_management import process_ligands_table as proc_ligs
from docking.docking_results_processing import process_docking_assay as dock_proc
##
##
### Step 0: Give the working project a name, create subfolder AND/OR export environment variables. Always leave this function activated since the global variables are exported by this function
project_base_dir = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform"
project_name = "sample"
raw_data_path, main_db_output_path, miscellaneous_files, receptor_models_registry_path,docking_assays_registry_path, docking_params_registry_path, docking_assays_storing_path, docking_raw_data_path = gen_proj.configure_project_structure(project_base_dir,project_name)
