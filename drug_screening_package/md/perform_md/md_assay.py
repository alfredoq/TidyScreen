import sys
import sqlite3
import pandas as pd

MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import md.perform_md.md_configuration as md_conf
from docking.docking_results_processing import process_docking_assay as dock_proc
from md.md_preparation_actions import ligand_preparation as md_lig_prep
from md.md_preparation_actions import receptor_preparation as md_rec_prep
from md.md_preparation_actions import complex_preparation as md_com_prep

def retrieve_md_cond_set(md_params_registry_path,md_params_id):
    """
    This function will retrieve an md_cond_set stored as a dictionary in the corresponding conditions database
    
    ------
    Parameters:
    ------
    - md_params_registry_path: the full path to the database in which the md_conditions are stored
    - md_params_id: the id (integer) 
    """
    
    database_file = f'{md_params_registry_path}/md_params.db'
    conn = sqlite3.connect(database_file)
    sql = f"""SELECT *
             FROM md_params
             WHERE condition_id = {md_params_id}"""
    try: 
        df = pd.read_sql_query(sql,conn)
        
        for index, row in df.iterrows():
            md_cond_dict = eval(row[1])
            
        return md_cond_dict
    
    except Exception as error:
        print(error)
        
def create_md_assay(receptor_models_registry_path,docking_assays_storing_path,docking_assays_registry_path,mds_path,docking_assay_nbr,ligand_table_name,ligand_subpose_name,receptor_model_id,md_params_id,mendieta_assay_path):
    """
    This function will prepare the whole md assay folder given a specific binding pose of a ligand
    
    """
    # Prepare the MD assays path variables based on the general path.
    md_assays_registry_path = f'{mds_path}/mds_assays_registers'
    md_params_registry_path = f'{mds_path}/mds_params_registers'
    md_assays_storing_path =  f'{mds_path}/mds_assays'

    # Retrieve relevant data form the docking assays registries
    ligand_db,source_ligands_db,source_ligands_pdbqt_table,rec_model_id, receptor_pdb_file = md_conf.retrieve_assay_relevant_variables(docking_assays_registry_path,docking_assay_nbr,receptor_models_registry_path)
    
    ## This will create the assay folder structure and store the corresponding registry
    md_assay_folder = md_conf.create_md_assay_folder(md_assays_registry_path,md_assays_storing_path,ligand_db,ligand_table_name,ligand_subpose_name,receptor_model_id,md_params_id)

    # This will retrieve the md_parameters set from the corresponding registries database

    md_cond_dict = retrieve_md_cond_set(md_params_registry_path,md_params_id)

    ## These functions will retrieve and prepare the corresponding ligand
    ligand_docked_pdb_file = md_lig_prep.retrive_and_prepare_docked_pose(ligand_db,ligand_subpose_name,source_ligands_db,ligand_table_name,source_ligands_pdbqt_table,md_assay_folder,md_cond_dict)

    ## This functions will retrieve the corresponding receptor
    receptor_pdb_file = md_rec_prep.retrive_receptor_file(receptor_pdb_file,md_assay_folder)

    ## This function will prepare the complex between the ligand and receptor
    pdb_complex_file = md_com_prep.prepare_complex(md_assay_folder,receptor_pdb_file,ligand_docked_pdb_file)
    
    ## This functions will prepare the corresponding .prmtop and .inpcrd files for the md simulation.
    md_com_prep.write_tleap_input(md_assay_folder,pdb_complex_file,md_cond_dict)
    
    ## This function will prepare the md inputs
    md_conf.prepare_md_inputs(md_assay_folder,md_cond_dict)