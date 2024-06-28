### This section will import all the required modules
### DO NOT MODIFY THESE IMPORTS.
import sys
from pprint import pp as pp
#sys.path.append('../drug_screening_package/')
sys.path.append('/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package')
## Imports related to the exploration of chemical space
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
from docking.receptor_management import receptor_management as rec_mng
from docking.perform_docking import docking_configuration as dock_cfg
from docking.perform_docking import docking_assay as dock_assay

###
### Set some system variables respect to paths for your own computer.###

## Variables affecting the exploration of chemical space
raw_data_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/raw_data"
main_db_output_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data"
project_db_output_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects"
miscellaneous_files = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing"
blob_files = ""
## Variables affecting the execution of docking assays
receptor_models_registry_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_models_registers"
docking_assays_registry_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers"
docking_params_registry_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers"
docking_assays_storing_path = "/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays"

###################################
### SCREENING OF CHEMICAL SPACE ###
###################################

## STEP 1: the raw data included in the emolecules reactants files is processed.
#help(acts_level_0.process_emolecules_raw_file_pandarallel)
#acts_level_0.process_emolecules_raw_file_pandarallel(f'{raw_data_path}/emolecules_database.smi',f'{main_db_output_path}/Emolecules_reagents_5.db')

## STEP 2: an early step in the workflow is to add SMARTS notations to the database in order to filter/subset reactants
#help(acts_level_0.add_reactant_filter_to_db)
#acts_level_0.add_reactant_filter_to_db(f'{main_db_output_path}/Emolecules_reagents_4.db')

## STEP 3: a frequent action is to filter/subset reactants based on a filter alread stored within the database.
#help(db_ops.select_reactants_based_on_smarts_pandarallel)
#db_ops.select_reactants_based_on_smarts_pandarallel(f'{main_db_output_path}/Emolecules_reagents_4.db','emolecules_reagents','alpha_aminoacid')

## STEP 4: next step is to compute the molecular properties of a table of molecules, so it can be further filtered.
#help(smiles_processing.compute_molecule_properties_on_table)
#smiles_processing.compute_molecule_properties_on_table(f"{main_db_output_path}/Emolecules_reagents_4.db","alpha_aminoacid")

## STEP 4: after a subseting of reactants is performed, it is useful to depict their structures.
#help(depic.generate_depiction_grid)
#depic.generate_depiction_grid(f'{main_db_output_path}/Emolecules_reagents_4.db','alpha_aminoacid',f'{miscellaneous_files}/aminoacids',25)

## STEP 5: it is usefull to plot the chemical space of a table of compounds. In this case a two steps procedure need to be applied:
## STEP 5.1: compute the umap coordinates for a given table and store it
#help(smiles_processing.compute_chemspace_coord)
#smiles_processing.compute_chemspace_coord(f'{main_db_output_path}/Emolecules_reagents_4.db','alpha_aminoacid')

## STEP 5.2: plot the already computed coordinates for a given table
#help(depic.plot_table_chemspace)
#depic.plot_table_chemspace(f'{main_db_output_path}/Emolecules_reagents_4.db','alpha_aminoacid',f'{miscellaneous_files}/alpha_aminoacid_space.png')

## STEP 6: after inspecting the matched structures, it may be useful to exclude compounds compounds from a table using a SMARTS filter.
#help(db_ops.exclude_reactants_based_on_smarts_pandarallel)
#db_ops.exclude_reactants_based_on_smarts_pandarallel(f'{main_db_output_path}/Emolecules_reagents_4.db',"alpha_aminoacid","di_carboxylic_acid","aa_w_o_diacids")

## STEP ... Here multiple steps of post-filtering the table may be required depending on the expected reaction limitations.

## STEP 7: now a set of reactions can be included in the database as SMARTS based reactions schemes.
#help(acts_level_0.add_reaction_smarts_to_db)
#acts_level_0.add_reaction_smarts_to_db(f'{main_db_output_path}/Emolecules_reagents_3.db')

## STEP 8: Sometimes it is useful to manually add some reactants from a custom .csv file. In this case a couple of alchohols aimed to aminoacid acilation are included
#help(db_ops.insert_smiles_in_table_from_file)
#db_ops.insert_smiles_in_table_from_file("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/alcohols.csv")

## STEP 9: After having a set of reactions within the database, we can apply the corresponding schema: two types of reactions are available: 1 component transformation and 2 components reaction
# STEP 9.1: This is 1 component transformation, for example primary amine to azide transformation:
#help(smiles_processing.apply_reaction_1_component_pandarallel)
#smiles_processing.apply_reaction_1_component_pandarallel(f'{main_db_output_path}/Emolecules_reagents_3.db',"azide_from_primary_amine","alpha_aminoacid",f'{main_db_output_path}/Emolecules_reagents_3.db',"azide_from_aminoacids","replace")

# STEP 9.2: This is 2 component transformation, for example the acylation of aminoacids with alchols:
#help(smiles_processing.apply_reaction_2_components_pandarallel)
#smiles_processing.apply_reaction_2_components_pandarallel(f"{main_db_output_path}/Emolecules_reagents_3.db","CuAAC","azide_from_aminoacids","alkine_benzotriazole_blocks",f"{main_db_output_path}/Emolecules_reagents_3.db","triazoles")

## STEP 10: At some point in the screening campaign, a separate database needs to be constructed, to resemble the particularities of the research project (e.g. CZP, SARS-COV2 Mpro, etc.). In this case a function to retrieve all the necesary information into this new databased is used.
#help(acts_level_1.create_level_1_project)
#acts_level_1.create_level_1_project(f"{main_db_output_path}/Emolecules_reagents_3.db",f"{project_db_output_path}/new_sample_level1.db",["alpha_aminoacid","alcohols"],"reaction_types",["azide_from_primary_amine","acylation"])

## STEP 11: another useful action in the workflow is to store files within the database. In this case, the files are stored as BLOB objects within a dedicated table. So, if the database does not contain a table por storing purposes, the first step is to create it, and the second step is to store the object.
# STEP 11.1: Creation of the BLOB object table for storing miscellaneous files
#help(blob_ops.create_blob_objects_table)
#blob_ops.create_blob_objects_table(f'{main_db_output_path}/Emolecules_reagents_4.db',"stored_items")

# STEP 11.2: Once the table exists, we can store all kind of files into it:
#help(blob_ops.append_blob_object_to_table)
#blob_ops.append_blob_object_to_table(f'{main_db_output_path}/Emolecules_reagents_4.db',"stored_items","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/raw_data/core_scheme.mrv")

#STEP 11.3: Finally, the stored object can be retrieved from the database:
#help(blob_ops.retrieve_blob_object)
#blob_ops.retrieve_blob_object(f'{project_db_output_path}/new_sample_level1.db',"stored_items","reactions_scheme","/home/fredy/Desktop")

# STEP 12: A useful feature is to check the presence of molecules given in a .csv file within a certain table
#help(db_actions.check_reactants_presence)
#db_actions.check_reactants_presence(f"{main_db_output_path}/Emolecules_reagents_3.db","alpha_aminoacid",'/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/all_natural_aminoacids.csv')

# STEP 13: Miscellaneous actions: these are some useful action required during the workflow management.
#help(db_actions.copy_table_between_dbs)
#db_actions.copy_table_between_dbs(f"{main_db_output_path}/Emolecules_reagents_3.db",f"{main_db_output_path}/Emolecules_reagents_4.db","reactants_types")

# STEP 14: Subseting a table based on parameters limits
#help(db_ops.subset_molecules_by_property_range)
#db_ops.subset_molecules_by_property_range(f"{main_db_output_path}/Emolecules_reagents_4.db","alpha_aminoacid","alpha_aminoacid_LogP_subset","MolLogP",1,3)


###################################
### EXECUTION OF DOCKING ASSAYS ###
###################################

## Docking STEP 1: a first step to initiate a molecular docking campaign is to create a database registering and containing all receptor models. This is performed only once in a screening campaing. Afterwards, once the registry database is created, the user should use the function 'append_receptor_registry'.
#help(rec_regs.create_receptor_registry_table)
#rec_regs.create_receptor_registry_table(receptor_models_registry_path)

## Docking STEP 2: a receptor model registry should be added to the receptor registry table. This receptor model includes all the files within a folder that are required for a docking assay: grid maps, .fld, description in txt, etc etc.
#help(rec_regs.append_receptor_registry)
#rec_regs.append_receptor_registry(receptor_models_registry_path,"/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/receptors_models/model_1")

## Docking STEP 3: In some cases, a receptor model stored within the receptors registry database need to be retrieved into a folder for inspection or whatever reason. This can be done as follows:
#help(rec_mng.restore_receptor)
#rec_mng.restore_receptor(receptor_models_registry_path,"/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/tmp")

## Docking STEP 4: After a receptor/s model/s has/have been deposited in the corresponding database, it is time to construct a confined assay folder containing all the subfolders structure to contain: ligands, the receptor, dlg files.
#help(dock_cfg.create_docking_assay_folder)
#dock_cfg.create_docking_assay_folder(docking_assays_registry_path,docking_assays_storing_path,"/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects/sample.db","aminoacids_pdbqt_ligands",3)
#help(dock_assay.create_docking_assay)
#dock_assay.create_docking_assay(receptor_models_registry_path,docking_assays_registry_path,docking_params_registry_path,docking_assays_storing_path,"/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects/sample.db","aminoacids_pdbqt_ligands",3,1)