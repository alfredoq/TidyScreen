
import sys

MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import md.perform_md.md_configuration as md_conf


def create_md_assay(receptor_models_registry_path,md_assays_registry_path,md_params_registry_path,md_assays_storing_path,ligand_db,ligand_table_name,ligand_subpose_name,rec_model_id,md_params_id,mendieta_assay_path):
    """
    This function will prepare the whole md assay folder given a specific binding pose of a ligand
    
    """

    ## This will create the assay folder structure and store the corresponding registry
    md_assay_folder = md_conf.create_md_assay_folder(md_assays_registry_path,md_assays_storing_path,ligand_db,ligand_table_name,ligand_subpose_name,rec_model_id,md_params_id)
