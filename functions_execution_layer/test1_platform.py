### This section will import all the required modules
### DO NOT MODIFY THESE IMPORTS.
import sys
from pprint import pp as pp
#sys.path.append('drug_screening_package')
sys.path.append('/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package')
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

## Imports related to the preparation and execution of molecular docking assays
from docking.registers_management import receptor_models_registers as rec_regs
from docking.registers_management import docking_params_registers as reg_dock_parm
from docking.perform_docking import docking_assay as dock_ass
from docking.receptor_management import receptor_management as rec_mng
from docking.perform_docking import docking_configuration as dock_cfg
from docking.perform_docking import docking_assay as dock_assay
from docking.ligands_management import process_ligands_table as proc_ligs
from docking.docking_results_processing import process_docking_assay as dock_proc

## Imports related to the machine learnign processing
from ml.fingerprints_colection import construct_fingerprints_database as cons_fps

### Step 0: Give the working project a name, create subfolder AND/OR export environment variables. Always leave this function activated since the global variables are exported by this function
project_base_dir = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform"
project_name = "TEST"
raw_data_path, main_db_output_path, miscellaneous_files, receptor_models_registry_path, docking_assays_registry_path, docking_params_registry_path, docking_assays_storing_path, docking_raw_data_path,ml_models_path = gen_proj.configure_project_structure(project_base_dir,project_name)

## STEP 1: the raw data included in the emolecules reactants files is processed.
#help(acts_level_0.process_emolecules_raw_file_pandarallel)
#source_filename = "emolecules_database.smi"
#destination_db = "emolecules_reagents.db"
#acts_level_0.process_emolecules_raw_file_pandarallel(f'{raw_data_path}/{source_filename}',f'{main_db_output_path}/{destination_db}')

## STEP 2: A custom csv file containing ligands in the SMILES format is read into a custom database and custom table.
#help(acts_level_0.process_raw_csv_file_pandarallel)
#source_filename = f'{raw_data_path}/MPRO_binders.csv'
#destination_db = f'{main_db_output_path}/MPRO_binders2.db'
#destination_table = "mpro_binders"
#acts_level_0.process_raw_csv_file_pandarallel(source_filename,destination_db,destination_table)

## STEP 3: after a subseting of reactants is performed, it is useful to depict their structures.
#help(depic.generate_depiction_grid)
#source_db = f'{main_db_output_path}/MPRO_binders2.db'
#table_name = "mpro_binders"
#images_dest_dir =  f'{miscellaneous_files}/{table_name}'
#max_mols_ppage = 25
#depic.generate_depiction_grid(source_db,table_name,images_dest_dir,max_mols_ppage)

# STEP 3: I will check if all MPRO_binders were successfully imported in the database
#help(db_actions.check_reactants_presence)
#source_db = "MPRO_binders2.db"
#table_name = "mpro_binders"
#source_filename = "MPRO_binders.csv"
#db_actions.check_reactants_presence(f"{main_db_output_path}/{source_db}",table_name,f'{raw_data_path}/{source_filename}')

# STEP 4: I will process a ligands table to generate .pdqbt files.
#help(proc_ligs.process_ligands_table)
#origin_db = "MPRO_binders.db"
#table_name = "mpro_binders"
#dest_db = "MPRO_binders.db"
#proc_ligs.process_ligands_table(f'{main_db_output_path}/{origin_db}', table_name, f'{main_db_output_path}/{dest_db}')

# STEP 5: Prior to initiation docking studies, a receptor registry table needs to be created to store all receptor models. If you want to create the table twice, the script will inform the action and exit
#help(rec_regs.create_receptor_registry_table)
#rec_regs.create_receptor_registry_table(receptor_models_registry_path)

# STEP 6: Now we need to add the receptor model files to the registry database. The models are read by folder name from 'docking_raw_data_path'
#help(rec_regs.append_receptor_registry)
#rec_folder = f'{docking_raw_data_path}/model_7vlp_pdb4amber'
#rec_regs.append_receptor_registry(receptor_models_registry_path,rec_folder)

# STEP 7: Now we need to create a docking parameters registry table.
#help(reg_dock_parm.create_params_registry_table)
#reg_dock_parm.create_params_registry_table(docking_params_registry_path)

# STEP 8: Now we need to create a custom docking parameters condition to be applied
#help(reg_dock_parm.store_custom_param_dict)
#reg_dock_parm.store_custom_param_dict(docking_params_registry_path)

# STEP 9: with all the above indicated steps completed, now we can configure a docking assay
#help(dock_ass.create_docking_assay)
#ligands_db = "MPRO_binders.db"
#ligands_table_name = "nirmatrelvir_selected_bioactive_stereo_pdbqt_ligands"
#rec_model_id = 1
#docking_params_id = 4
#dock_ass.create_docking_assay(receptor_models_registry_path,docking_assays_registry_path,docking_params_registry_path,docking_assays_storing_path,f'#{main_db_output_path}/{ligands_db}',ligands_table_name,rec_model_id,docking_params_id)

# STEP 10: Subset the Mpro binders based on the index in order to dock selected representative compounds.
#help(db_ops.subset_molecules_based_on_index)
#db_name = f'{main_db_output_path}/MPRO_binders.db'
#table_origin = "mpro_binders"
#table_dest = "nirmatrelvir_selected"
#db_ops.subset_molecules_based_on_index(db_name, table_origin, table_dest)

# STEP 11.1: Enumerate the stereoisomers on the subseted table
#help(smiles_processing.enumerate_stereoisomers)
#db_name = f'{main_db_output_path}/MPRO_binders.db'
#table_name = "nirmatrelvir_selected"
#smiles_processing.enumerate_stereoisomers(db_name,table_name)

# STEP 11.2: Enumerate the stereoisomers of a single selected molecule
#help(smiles_processing.enumerate_stereoisomers_selected_molecule)
#db_name_origin = f'{main_db_output_path}/MPRO_binders.db'
#db_name_dest = f'{main_db_output_path}/MPRO_binders.db'
#origin_table = "mpro_binders_selected"
#dest_table = "mpro_binders_nirmatrelvir_stereo_enum"
#row_id = 1
#smiles_processing.enumerate_stereoisomers_selected_molecule(db_name_origin, db_name_dest, origin_table, dest_table, row_id)

## STEP 12: after a subseting of reactants is performed, it is useful to depict their structures.
#help(depic.generate_depiction_grid)
#source_db = f'{main_db_output_path}/MPRO_binders.db'
#table_name = "nirmatrelvir_selected_stereo_enum"
#images_dest_dir =  f'{miscellaneous_files}/{table_name}'
#max_mols_ppage = 25
#depic.generate_depiction_grid(source_db,table_name,images_dest_dir,max_mols_ppage)

## STEP 13: filter specific stereoisomers from the enumerated table.
#help(db_ops.subset_molecules_by_column_value)
#db_name = f'{main_db_output_path}/MPRO_binders.db'
#origin_table = "nirmatrelvir_selected_stereo_enum"
#dest_table = "nirmatrelvir_selected_bioactive_stereo"
#column_name = "iso_config"
#value = "SSRSSS"
#action = "append"
#db_ops.subset_molecules_by_column_value(db_name,origin_table,dest_table,column_name,value,action)

## STEP 14: Isolate the bioactive Nirmatrelvir configuration
#help(db_ops.subset_molecules_by_column_value)
#db_name = f'{main_db_output_path}/MPRO_binders.db'
#origin_table = "mpro_binders_8ign_ligand_stereo_enum"
#dest_table = "mpro_binders_8ign_ligand_bioactive_isomer"
#column_name = "iso_config"
#value = "SSSSRS"
#action = "append"
#db_ops.subset_molecules_by_column_value(db_name,origin_table,dest_table,column_name,value,action)

# STEP 15: Generate the corresponding .pdqbt files for docking.
#help(proc_ligs.process_ligands_table)
#origin_db = f'{main_db_output_path}/MPRO_binders.db'
#table_name = "nirmatrelvir_selected_bioactive_stereo"
#dest_db = f'{main_db_output_path}/MPRO_binders.db'
#proc_ligs.process_ligands_table(origin_db, table_name, dest_db)

# STEP 16: configure a docking assay
#help(dock_ass.create_docking_assay)
#ligands_db = f'{main_db_output_path}/MPRO_binders.db'
#ligands_table_name = "mpro_binders_8ign_ligand_bioactive_isomer_pdbqt_ligands"
#rec_model_id = 1
#docking_params_id = 7
#dock_ass.create_docking_assay(receptor_models_registry_path,docking_assays_registry_path,docking_params_registry_path,docking_assays_storing_path,#ligands_db,ligands_table_name,rec_model_id,docking_params_id)

# STEP 17: perform and store the whole set of analysis for the docking assay
#help(dock_proc.perform_docking_assay_analysis)
#assay_nbr = 3
#dock_proc.perform_docking_assay_analysis(receptor_models_registry_path, docking_assays_registry_path, docking_assays_storing_path,#docking_raw_data_path, assay_nbr)

# STEP 18: Once the results have been processed, prepare the results for ml processing. The first step is to create a fingerprints database and table.
#help(cons_fps.create_fingerprints_table)
#cons_fps.create_fingerprints_table(ml_models_path)

# STEP 19: Store a positive fingerprint identified during the molecular docking assays;
#help(cons_fps.store_positive_fingerprint)
#assay_nbr = 1
#sub_pose = "LIENCHBZNNMNKG-OJFNHCPVSA-N_1"
#cons_fps.store_positive_fingerprint(ml_models_path, receptor_models_registry_path, docking_assays_registry_path,docking_assays_storing_path, assay_nbr, sub_pose)

# STEP 19: Store a negative fingerprint identified during the molecular docking assays;
#help(cons_fps.store_negative_fingerprint)
#assay_nbr = 1
#sub_pose = "LIENCHBZNNMNKG-OJFNHCPVSA-N_1"
#cons_fps.store_negative_fingerprint(ml_models_path, receptor_models_registry_path, docking_assays_registry_path,docking_assays_storing_path, assay_nbr, sub_pose)