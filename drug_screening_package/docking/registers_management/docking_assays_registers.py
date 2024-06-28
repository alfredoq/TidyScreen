import sqlite3
import pandas as pd
from termcolor import colored

import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.input_output_operations.actions_on_dbs as db_acts

def create_docking_registry_table(general_registries_folder):
    """
    This function will create a table designed to store information related to performed molecular docking campaigns.

    The table will contain the following columns:

        - assay_id: an integer values to store a specific docking assay;
        - db_name: the full path to the database containing the table of docked molecules ;
        - ligs_table_name: the name of the table localed in 'db_name' that contains the ligands in .pdbqt format.
        - rec_model_id: an integer indicating the receptor model used as stored within the receptors registers table;
        - stored_results: a flag indicating if the docking assay results are available in the assay table (0:no, 1:yes)
        - assay_folder: the path the folder were all the info to run the assay was created.
        - assay_description: a brief description of the objective of the specific docking assay.

    ATTENTION: if the table named 'docking_registries' is stored within the a database '{general_registries_folder/docking_registries.db} it will be replaced.

    The function should not be called again if the database/table has already been created, unless you want to clear the overall table of registries.
        
    ------
    Parameters:
    ------
    - general_registries_folder: the full path were the docking assays registry table is to be stored

    ------
    Returns:
    ------
    Creates an empty registry table. 
    """
    
    try: 
        conn = sqlite3.connect(f'{general_registries_folder}/docking_registries.db')

        ins = f'CREATE TABLE docking_registries ("assay_id" INTEGER, "db_name" TEXT, "ligs_table_name" TEXT, "rec_model_id" INTEGER,"docking_params_id" INTEGER, "stored_results" INTEGER, "assay_folder" TEXT, assay_description TEXT)'
        conn.execute(ins)
        print("Docking registers table created.")
        
    except Exception as error:
        print(error)
        print(colored("Error creating the database docking_registries.db","red"))

def erase_docking_assays_registry_table(general_registries_folder):
    """
    This function is used to completely delete all molecular docking assays registries stored in the general table.

    The table named 

    USE WITH CAUTION since a lot of information will be lost.
    
    Executing this function will recreate the corresponding table.
    ------
    Parameters:
    ------
    - docking_parameters_path: the full path were the docking registry table is to be stored
    """

    try: 
        db_acts.drop_table_from_db(f'{general_registries_folder}/docking_registries.db',"docking_registries")
        print(colored(f"Database 'docking_registries' within {general_registries_folder}/docking_registries.db has been DELETED ","green"))

    except Exception as error:
        print(error)
        print(colored("Problem deleting the registries table","red"))

def clear_docking_assays_registry_table(general_registries_folder):
    """
    This function will erase all the docking assays registries stored in the 'general_registries_folder' path within the table: 'docking_registries' stored in the 'docking_registries.db' database.

    ------
    Parameters:
    ------
    - general_registries_folder: the full path were the docking assays registry table is to be stored

    ------
    Returns:
    ------
    None
    """
    try: 
        conn = sqlite3.connect(f'{general_registries_folder}/docking_registries.db')
        cursor = conn.cursor()
        sql = "DELETE FROM docking_registries" # This will delete all the registries in the table.
        query = input(colored("Are you sure you want to delete all docking assay registries? (y/n): ","magenta"))
        if query == 'y':
            cursor.execute(sql)
            conn.commit()
            print(colored("All registries docking assays registries DELETED.","green"))
        else:
            exit()
    except Exception as error:
        print(error)
        print(colored("Deletion of docking assays registries FAILED","red"))

def get_last_docking_assay_id(docking_assays_registry_path):
    """
    This function will search for the last assay_id value in the 'docking_registries' table of the 'docking_registries.db' database under the path {general_registries_folder}.
    
    In case no registry value is found, '1' will be returned, assuming that the table is new or has been intentionally deleted.
    
    ------
    Parameters:
    ------
    - general_registries_folder: the full path to the folder in which the molecular docking registries are being stored.
    
    ------
    Returns:
    ------
    - last_value = The last incremental value stored as 'assay_id' for the docking runs.
    """
    
    try:
        conn = sqlite3.connect(docking_assays_registry_path+'/docking_registries.db')
        sql_assay = """SELECT assay_id 
                       FROM docking_registries;"""
    
        df = pd.read_sql_query(sql_assay,conn)
        assay_id_list = df["assay_id"].to_list()
        if len(assay_id_list) == 0: # This will hit in case no previous record in the table is present
            return 0
        else:
            last_value = assay_id_list[-1]
            return last_value
    
    except Exception as error:
        print(error)
        print(colored("Error retrieving the last assay_id value","red"))

def append_docking_registry(general_registries_folder,db_name, table_name, rec_model, stored_results, assay_folder,description,docking_params_id):
    """
    When preparing a docking assay, this function will append a registry with the required information into the registries table.
    
    ------
    Parameters:
    ------
    - general_registries_folder: the full path the folder in which the general docking registries are being stored.
    - db_name: the full path to the database containing the ligands to be docked.
    - table_name: the name of the table within 'db_name' that contains the ligands in the .pdbqt format.
    - rec_model: an integer identifying the receptor models that is going to be used in the docking assay.
    - stored_results: indicates wether the results of the docking has been processed and stored in the corresponding dabatase.
    - assay_folder: the full path to the folder in which the docking assay is to be stored.
    - description: a brief description of the objective of the specific docking assay
    
    ------
    Returns:
    ------
    A record in the register number is created with the provided information.
    """
    # Check if the assay_id already exists, if True, exit.
    
    try:
        conn = sqlite3.connect(general_registries_folder+'/docking_registries.db')
        sql_assay = """SELECT assay_id 
                       FROM docking_registries;"""
    
        # This will retrieve the last value of the assay_id in the registry table and increment by 1 for the next storage
        last_value = get_last_docking_assay_id(general_registries_folder)
        assay_id = last_value+1
    
        # Get all the values present in the table
        df = pd.read_sql_query(sql_assay,conn)
        assay_id_list = df["assay_id"].to_list()
        
        # Check if assay_id exists to prevent as a second factor repeating the id. Should never match because assay_is is automatically generated
        if assay_id in assay_id_list:
            print(colored(f"The docking assay_id {assay_id} already exists in the registry table. EXITING","red"))
            exit 
        else:
            raise # This will manually raise an exception to continue
    except:
        sql = f"""INSERT INTO docking_registries(assay_id,db_name,ligs_table_name,rec_model_id,docking_params_id,stored_results,assay_folder,assay_description) VALUES({assay_id},"{db_name}","{table_name}",{rec_model},{docking_params_id},{stored_results},"{assay_folder}","{description}");"""
        conn.execute(sql)
        conn.commit()
        print(colored(f"Docking assay_id {assay_id} is NEW, stored the registry","green"))

#################

if __name__ == '__main__':

    print("This is a test")

    #create_docking_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers")

    erase_docking_assays_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers")

    #append_docking_registry("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/sample.db","alpha_aminoacid_pdbqt_ligands",1,0)
    
    #last_value = get_last_docking_assay_id("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers")
    #print(last_value)

    #clear_docking_assays_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays_registers")