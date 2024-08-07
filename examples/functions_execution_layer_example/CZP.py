###This section will import all the required modules
###DO NOT MODIFY THESE IMPORTS.
import sys
from pprint import pp as pp
sys.path.append('/home/fredy/Desktop/test/TidyScreen/drug_screening_package')
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
from docking.registers_management import docking_assays_registers as reg_dock_ass
from docking.registers_management import docking_params_registers as reg_dock_parm
from docking.perform_docking import docking_assay as dock_ass
from docking.receptor_management import receptor_management as rec_mng
from docking.perform_docking import docking_configuration as dock_cfg
from docking.perform_docking import docking_assay as dock_assay
from docking.ligands_management import process_ligands_table as proc_ligs
from docking.docking_results_processing import process_docking_assay as dock_proc
##
## Imports related to the machine learnign processing
from ml.fingerprints_colection import construct_fingerprints_database as cons_fps
## Imports related to the molecular dynamics simulations
from md.perform_md import md_assay as md_assay
##
##
### Step 0: Give the working project a name, create subfolder AND/OR export environment variables. Always leave this function activated since the global variables are exported by this function
project_base_dir = "/home/fredy/Desktop/CZP"
project_name = "CZP_project_results"
raw_data_path, main_db_output_path, miscellaneous_files, receptor_models_registry_path,docking_assays_registry_path, docking_params_registry_path, docking_assays_storing_path, docking_raw_data_path,docking_misc_data_path,ml_models_path,mds_path,qm_mm_mds_path = gen_proj.configure_project_structure(project_base_dir,project_name)

## STEP 1: A custom csv file containing ligands in the SMILES format is read into a custom database and custom table.
#help(acts_level_0.process_raw_csv_file_pandarallel)
""" source_filename = f'{raw_data_path}/K777.csv'
destination_db = f'{main_db_output_path}/CZP_binders.db'
destination_table = source_filename.split('/')[-1].replace('.csv','')
retain_stereo = 0
acts_level_0.process_raw_csv_file_pandarallel(source_filename,destination_db,destination_table,retain_stereo) """


## STEP 2: This will enumerate all stereoisomers within a table.
#help(smiles_processing.enumerate_stereoisomers)
""" db_name = f'{main_db_output_path}/CZP_binders.db'
table_name = "K777"
smiles_processing.enumerate_stereoisomers(db_name,table_name) """


## STEP 3: Separate the ligands mantaining the S,S stereo
#help(db_ops.subset_molecules_by_column_value)
""" db_name = f'{main_db_output_path}/CZP_binders.db'
origin_table = "K777_stereo_enum"
dest_table = "K777_stereo_enum_SS"
column_name = "iso_config"
value = "SS"
action = "append"
db_ops.subset_molecules_by_column_value(db_name,origin_table,dest_table,column_name,value,action) """


## STEP 4: after a subseting of reactants is performed, it is useful to depict their structures.
#help(depic.generate_depiction_grid)
""" source_db = f'{main_db_output_path}/CZP_binders.db'
table_name = "K777_stereo_enum"
images_dest_dir =  f'{miscellaneous_files}/{table_name}'
max_mols_ppage = 25
depic.generate_depiction_grid_mol_id(source_db,table_name,images_dest_dir,max_mols_ppage) """


# STEP 5: Generate the .pdqbt files for the bioactive enantiomer of Nirmatrelvir.
#help(proc_ligs.process_ligands_table)
""" origin_db = f'{main_db_output_path}/CZP_binders.db'
table_name = "K777_stereo_enum"
dest_db = f'{main_db_output_path}/CZP_binders.db'
write_conformers = 0 # If activated ==1 - will write all conformers searched for ligand .pdbqt generation
conformer_rank = 0 # The rank of the conformer in the search to be used as ligand. 0-based index (0 is the lowest energy conformer)
proc_ligs.process_ligands_table(origin_db, table_name,dest_db,conformer_rank,write_conformers) """

# STEP 6: Creation of the receptor registry table to store all receptor models.
#help(rec_regs.create_receptor_registry_table)
""" rec_regs.create_receptor_registry_table(receptor_models_registry_path) """

# STEP 7: adding the receptor model files to the registry database. The models are read by folder name from 'docking_raw_data_path'
#help(rec_regs.append_receptor_registry)
""" rec_folder = f'{docking_raw_data_path}/2OZ2'
rec_regs.append_receptor_registry(receptor_models_registry_path,rec_folder) """

# STEP 8: Create of docking parameters registry table.
#help(reg_dock_parm.create_params_registry_table)
""" reg_dock_parm.create_params_registry_table(docking_params_registry_path) """

# STEP 9: Create a custom docking parameters condition to be applied
#help(reg_dock_parm.store_custom_param_dict)
""" reg_dock_parm.store_custom_param_dict(docking_params_registry_path) """

#STEP 10: Create a docking assays registry table
#help(reg_dock_ass.create_docking_registry_table)
""" reg_dock_ass.create_docking_registry_table(docking_assays_registry_path) """

# STEP 11: Configuration of a docking assay
#help(dock_ass.create_docking_assay)
ligands_db = "CZP_binders.db"
ligands_table_name = 'K777_stereo_enum_pdbqt_ligands'
rec_model_id = 1
docking_params_id = 1
hpc_assay_path = "/$PATH/TO/STORE/IN/HPC"
n_chunk = 10000
dock_ass.create_docking_assay(receptor_models_registry_path,docking_assays_registry_path,docking_params_registry_path,docking_assays_storing_path,f'{main_db_output_path}/{ligands_db}',ligands_table_name,rec_model_id,docking_params_id,hpc_assay_path,n_chunk)

