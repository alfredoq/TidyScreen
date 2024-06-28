import sys
import os
from termcolor import colored

MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import md.md_registers_management.md_assays_registers as md_regs


def create_md_assay_folder(md_assays_registry_path,md_assays_storing_path,ligand_db,ligand_table,ligand_subpose_name,rec_model_id,md_params_id):

    """
    This function will create a folder to store all the information required to run an md assay starting from a selected binding mode as obtained in a molecular docking assay.

    The function will query the last assay_id registered, and create the next consecutive folder.

    ------
    Parameters:
    ------

    ------
    Returns:
    ------
    A folder in which the information for a specific docking assay is stored.
    """

    last_value = md_regs.get_last_md_assay_id(md_assays_registry_path)

    # This will check is the database of registries exists, if not, will create it
    if last_value == None:
        md_regs.create_md_registry_table(md_assays_registry_path)
        print(colored("The docking assays registry tables did not existed. It was created","red"))
        last_value = 0

    md_assay_folder = f'{md_assays_storing_path}/md_assay_{last_value+1}'
    
    try: 
        ## This will append the md assay registry to the corresponding assays registries table
        description = input(colored("Provide a brief description of the MD assay: ","green"))
        md_regs.append_md_registry(md_assays_registry_path,ligand_db,ligand_table,ligand_subpose_name,rec_model_id,0,md_assay_folder,description,md_params_id)
        ## This will create the subfolders structure for the MD assay
        os.mkdir(md_assay_folder)
        print(colored(f"MD assay folder: md_assay_{last_value+1} CREATED","green"))
        
        return md_assay_folder
        
    except Exception as error:
        print(colored("Error creating the MD assay folder","red"))
        print(error)




