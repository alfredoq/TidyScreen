import os
from termcolor import colored

current_dir = os.getcwd()

# Check if the execution has been done within the 'functions_execution_layer' folder:
if 'functions_execution_layer' in current_dir: # This will match if this script is executed within the folder
    previous_dir_list = current_dir.split('/')[:-1]
    package_name = previous_dir_list[-1]
    previous_dir = '/'.join(previous_dir_list)
    current_dir = previous_dir

import sys

sys.path.append(f'{current_dir}/drug_screening_package')
from project_structure import generate_project_structure as gen_struct


# This will rename all the system path of the modules
gen_struct.set_modules_path(f'{current_dir}/drug_screening_package')
project_name = input("Enter the name of the general project: ")
project_base_folder = input("Enter the project base folder (full path): ")
project_full_path = f'{project_base_folder}/{project_name}'

# Check if the folder exists. If TRUE, inform and exit to avoid overwritting. FALSE: create the project and proceed.
# Update
if os.path.exists(project_full_path):
    print(colored(f"The project {project_full_path} alredy exists. EXITING","red"))
    exit()
else:
    os.makedirs(f'{project_full_path}/functions_execution_layer')
    print(colored(f"Project folder {project_full_path} was created successfully.","green"))

with open(f"{project_full_path}/functions_execution_layer/{project_name}.py","w") as project_template:
    project_template.write('###This section will import all the required modules\n')
    project_template.write('###DO NOT MODIFY THESE IMPORTS.\n')
    project_template.write('import sys\n')
    project_template.write('from pprint import pp as pp\n')
    project_template.write(f"sys.path.append('{current_dir}/drug_screening_package')\n")
    project_template.write('##\n')
    project_template.write('## Imports related to the exploration of chemical space\n')
    project_template.write('from project_structure import generate_project_structure as gen_proj\n')
    project_template.write('import chemspace.actions_level_0 as acts_level_0 \n')
    project_template.write('import chemspace.actions_level_1 as acts_level_1 \n')
    project_template.write('import chemspace.general_procedures.operations_on_sql_database as db_ops\n')
    project_template.write('from  chemspace.general_procedures import smiles_processing as smiles_processing\n')
    project_template.write('from  chemspace.general_procedures import depiction_procedures as depic\n')
    project_template.write('from  chemspace.general_procedures import objects_processing as obj_proc\n')
    project_template.write('from  chemspace.general_procedures import blob_operations as blob_ops\n')
    project_template.write('from  chemspace.input_output_operations import actions_on_dbs as db_actions\n')
    project_template.write('from  chemspace.input_output_operations import i_o_checkings as i_o_chk\n')
    project_template.write('from  chemspace.input_output_operations import read_data_from_db as read_data\n')
    project_template.write('##\n')
    project_template.write('## Imports related to the preparation and execution of molecular docking assays\n')
    project_template.write('from docking.registers_management import receptor_models_registers as rec_regs\n')
    project_template.write('from docking.registers_management import docking_assays_registers as reg_dock_ass\n')
    project_template.write('from docking.registers_management import docking_params_registers as reg_dock_parm\n')
    project_template.write('from docking.perform_docking import docking_assay as dock_ass\n')
    project_template.write('from docking.receptor_management import receptor_management as rec_mng\n')
    project_template.write('from docking.perform_docking import docking_configuration as dock_cfg\n')
    project_template.write('from docking.perform_docking import docking_assay as dock_assay\n')
    project_template.write('from docking.ligands_management import process_ligands_table as proc_ligs\n')
    project_template.write('from docking.docking_results_processing import process_docking_assay as dock_proc\n')
    project_template.write('##\n')
    project_template.write('## Imports related to the machine learnign processing\n')
    project_template.write('from ml.fingerprints_colection import construct_fingerprints_database as cons_fps\n')
    project_template.write('## Imports related to the molecular dynamics simulations\n')
    project_template.write('from md.perform_md import md_assay as md_assay\n')
    project_template.write('from md.md_registers_management import md_assays_registers as md_assay_reg\n')
    project_template.write('from md.md_registers_management import md_conditions_registers as md_cond_reg\n')
    project_template.write('from md.md_analysis import md_analyses as md_analyses\n')
    project_template.write('from md.md_preparation_actions import ligand_preparation as lig_prep\n')
    
    project_template.write('##\n')
    project_template.write('##\n')
    project_template.write('### Step 0: Give the working project a name, create subfolder AND/OR export environment variables. Always leave this function activated since the global variables are exported by this function\n')
    project_template.write(f'project_base_dir = "{project_full_path}"\n')
    project_template.write(f'project_name = "{project_name}_project_results"\n')
    project_template.write('raw_data_path, main_db_output_path, miscellaneous_files, receptor_models_registry_path,docking_assays_registry_path, docking_params_registry_path, docking_assays_storing_path, docking_raw_data_path,docking_misc_data_path,ml_models_path,mds_path,qm_mm_mds_path = gen_proj.configure_project_structure(project_base_dir,project_name)\n')
    
