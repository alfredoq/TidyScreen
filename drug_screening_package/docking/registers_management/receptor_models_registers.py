import sqlite3
from termcolor import colored
import pandas as pd
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import docking.receptor_management.receptor_management as rec_mng

def create_receptor_registry_table(receptor_models_registry_path):
    """
    This function will create a database called 'receptors_models_registry.db' within the 'receptor_models_registry_path', with a table named 'receptor_models', which contains the following structure:

        - rec_model_id: an integer value that identifies the corresponding receptor model.
        - pdb_identifier: a string containing the pdb code of the receptor
        - receptor_object: a blob object containing the compresed file of the receptor structure and grid files.
        - description: a string to add text description to the receptor model.
        - pdb_file: the full path to the pdb file originating the corresponding .pdbqt receptor and grid files

    ------
    Parameters:
    ------
    - receptor_models_registry_path: the full path the folder in which the receptors models registries are stored.

    ------
    Returns:
    ------
    Creates a database and table to store receptor registry information.
    """
    conn = sqlite3.connect(f'{receptor_models_registry_path}/receptors_models_registry.db')

    try:
        ins = f'CREATE TABLE receptor_models ("rec_model_id" INTEGER, "pdb_identifier" TEXT, "receptor_object" BLOB, "description" TEXT, "pdb_file" TEXT)'
        conn.execute(ins)
        print(colored(f"Docking registers table created in: '{receptor_models_registry_path}'.","green"))
    except Exception as error:
        print(error)
        print(colored(f"Error creating receptors registries table in {receptor_models_registry_path}","red"))

def drop_receptor_registry_table(receptor_models_registry_path):
    """
    This function will drop the table 'receptor_models' from the database 'receptors_models_registry.db' stored in the 'receptor_models_registry_path'.

    ------
    Parameters:
    ------
    - receptor_models_registry_path: the full path the folder in which the receptors models registries are stored.

    ------
    Returns:
    ------
    Completely deletes the table 'receptor_models'.
    """
    try: 
        conn = sqlite3.connect(f'{receptor_models_registry_path}/receptors_models_registry.db')
        sql = "DROP TABLE receptor_models"
        query = input(colored("Are you sure to delete the table 'receptor_models'? (y/n): ","magenta"))
        if query == 'y':
            conn.execute(sql)
            conn.close()
            print(colored("SUCCESSFULLY deleted the table 'receptors_models'","green"))
        else: 
            print(colored("Canceling the deletion of table 'receptors_models'","red"))
            exit()
    except Exception as error:
        print(error)
        print(colored("Error deleting the table 'receptors_models'","red"))

def clear_receptor_registry_table(receptor_models_registry_path):
    """
    This function will erase all the records available in the table 'receptor_models' from the database 'receptors_models_registry.db' stored in the 'receptor_models_registry_path'.

    ------
    Parameters:
    ------
    - receptor_models_registry_path: the full path the folder in which the receptors models registries are stored.

    ------
    Returns:
    ------
    Completely deletes the records stored in the 'receptor_models'.
    """

    try: 
        conn = sqlite3.connect(f'{receptor_models_registry_path}/receptors_models_registry.db')
        cursor = conn.cursor()
        sql = "DELETE FROM receptor_models" # This will delete all the registries in the table.
        query = input(colored("Are you sure you want to clear all recetor model registries? (y/n): ","magenta"))
        if query == 'y':
            cursor.execute(sql)
            conn.commit()
            print(colored("All receptor models registries DELETED.","green"))
        else:
            exit()
    except Exception as error:
        print(error)
        print(colored("Deletion of receptor models registries FAILED","red"))

def get_last_receptor_id(receptor_models_registry_path):

    conn = sqlite3.connect(f'{receptor_models_registry_path}/receptors_models_registry.db')

    sql_model = """SELECT rec_model_id 
                       FROM receptor_models;"""

    try:     
        df = pd.read_sql_query(sql_model,conn)
        model_id_list = df["rec_model_id"].to_list()
        if len(model_id_list) == 0: # This will hit in case no previous record in the table is present
            return 0
        else:
            last_value = model_id_list[-1]
            return last_value

    except Exception as error:
        print(error)
        print(colored("Error retrieveing the last receptor model ID value","red"))

def append_receptor_registry(receptor_models_registry_path,receptor_folder):
    """
    This function will add a receptor registry to the 'receptor_models' table located in a database called 'receptors_models_registry.db' stored in the 'receptor_models_registry_path'.

    The register will be automatically added following a incremental 'rec_model_id' number.

    ------
    Parameters:
    ------
    - receptor_models_registry_path: this is the full path to the folder in whichn the 'receptors_models_registry.db' database is found.
    - receptor_folder: the name of the folder within the docking raw data path in which all files related to the receptor model are stored.

    ------
    Returns:
    ------
    A new registry is added to the '${receptor_models_registry_path}/receptors_models_registry.db' database.
    """
    last_value = get_last_receptor_id(receptor_models_registry_path)
    rec_mng.store_receptor(f'{receptor_models_registry_path}/receptors_models_registry.db',receptor_folder,last_value)

#################

if __name__ == '__main__':

    print("This is a test")

    #create_receptor_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_models_registers")

    #drop_receptor_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_models_registers")

    append_receptor_registry("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_models_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/receptors_models/model_1")

    #clear_receptor_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_models_registers")
