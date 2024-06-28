import sqlite3
import pandas as pd
import os.path
from termcolor import colored
import os


import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.general_procedures.operations_on_sql_database as sql_ops
import docking.registers_management.docking_assays_registers as dock_regs

def create_docking_assay_folder(docking_assays_registry_path,docking_assays_storing_path,ligands_db,ligands_table,rec_model_id,docking_params_id):
    """
    This function will create a folder to store all the information required to run a docking assay.

    The function will query the last assay_id registered, and create the next consecutive folder.

    ------
    Parameters:
    ------
    - docking_assays_registry_path: the full path to the folder in which docking assays registers are stored.
    - docking_assays_storing_path: the full path to the folder in which ALL docking assays folders are stored.
    - rec_model_id: the receptor model id that identifies the receptor to use. This value is consistent with the one stored in the receptors registries database.

    ------
    Returns:
    ------
    A folder in which the information for a specific docking assay is stored.
    """

    last_value = dock_regs.get_last_docking_assay_id(docking_assays_registry_path)
    
    # This will check is the database of registries exists, if not, will create it
    if last_value == None:
        dock_regs.create_docking_registry_table(docking_assays_registry_path)
        print(colored("The docking assays registry tables did not existed. It was created","red"))
        last_value = 0
    
    docking_assay_folder = f'{docking_assays_storing_path}/docking_assay_{last_value+1}'
    
    try: 
        ## This will append the docking assay registry to the corresponding assays registries table
        description = input(colored("Provide a brief description of the docking assay: ","magenta"))
        dock_regs.append_docking_registry(docking_assays_registry_path,ligands_db,ligands_table,rec_model_id,0,docking_assay_folder,description,docking_params_id)
        ## This will create the subfolders structure for the docking assay
        os.mkdir(docking_assay_folder)
        os.mkdir(f'{docking_assay_folder}/pdbqt_files')
        os.mkdir(f'{docking_assay_folder}/receptor')
        os.mkdir(f'{docking_assay_folder}/dlgs')
        print(colored(f"Docking assay folder: docking_assay_{last_value+1} CREATED","green"))
        
        return docking_assay_folder
        
    except Exception as error:
        print(colored("Error creating the docking assay folder","red"))
        print(error)



#################

if __name__ == '__main__':
    
    #create_docking_assay_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays/assay_1.db")
    
    #create_assay_database("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays",1,"/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/sample.db","alpha_aminoacid_pdbqt_ligands","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects/sample.db",'test',1)

    create_docking_assay_folder("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays")