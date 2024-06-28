import sqlite3
from termcolor import colored
import xtarfile as tarfile
import pandas as pd
import tarfile 


import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import docking.receptor_management.blob_operations as blob_ops
import docking.input_output_operations.actions_on_dbs as acts_dbs

def create_receptors_table(db_name):
    """
    This function will create a table specifically designed to store receptors models.
    
    This function will create a fixed name table called 'receptors_models'.
    
    If the table already exists within the database, it will stop de execution,
    
    The receptors table is composed by the following rows:

        - rec_nbr: this field will be used as an integer value identifying the receptor model.
        - identifier: a value identifying the receptor stored in the row.
        - receptor_file: a BLOB object field that contains a .zip file with all the files associated to the receptor that are used in the docking procedure.
        - description: a brief comment of the receptor stored in the file.
    
    ------
    Parameters
    ------
    - db_name: the full path to the database file in which the table is going to be created.
    - table_name
    """
    
    conn = sqlite3.connect(db_name)
    
    try: # This will end the execution if the table is found
        sql_1 = f"""SELECT * 
                    FROM receptors_models'"""
        conn.execute(sql_1)
        print('Table already exists, please ', colored("append the object.","red"))
    
    except: # This will proceed with the table creation in case the table is not found
        print(colored("Creating ","green") + "new BLOB object table for storing receptors.")
        
        sql_2 = f"""CREATE TABLE receptors_models 
                    (rec_nbr INTEGER, identifier TEXT NOT NULL, blob_item BLOB NOT NULL, filename TEXT NOT NULL, description TEXT NOT NULL);"""      
        conn.execute(sql_2)

def store_receptor(db_name,from_path,last_value):
    """
    This function will store a receptor folder as a gzip file within a 'receptor_models' table.
    
    The function will ask for the following information:
    
        - rec_nbr: this is an integer that will be used to identify the receptor model.
        - identifier: a one_word description identifying the receptor model object within the database.
        - description: a brief description of the receptor model to be stored.
    
    ------
    Parameters:
    ------
    - db_name: the full path to the database file were the receptor folder is to be stored.
    - from_path: the full path to the folder containing all the receptor information that is intended to be stored within the datbase.
    
    ------
    Returns
    ------
    The receptor folder is saved in compressed format to the specified database file.
    """ 

    conn = sqlite3.connect(db_name)
    
    #First determine if 'receptors_models' exists within the database
    try:
        sql = f"""SELECT * 
              FROM receptor_models;"""
        df = pd.read_sql_query(sql,conn)
        print("Table 'receptor_models' found in database.")
        
        ## Continue with the model storage procedure
        # Create the tar file within the /tmp di for storage
        model_name = from_path.split('/')[-1]
        tar = tarfile.open(f"/tmp/{model_name}.tar.gz", "w:gz")
        #print(from_path)
        tar.add(from_path, arcname=model_name)
        #print(model_name)
        tar.close()
        
        rec_nbr = last_value+1
        pdb_file = input("Enter the name of the pdb file to store: ")
        pdb_file_path = f'{from_path}/{pdb_file}'
        identifier = pdb_file.replace(".pdb","")
        description = input("Enter a brief description of the features of the receptor model: ")    
    
        # Create the blob object for storing
        blob_obj = blob_ops.create_blob_object(f"/tmp/{model_name}.tar.gz")
        # Use the corresponding function to store the receptor object
        
        print(db_name)

        blob_ops.append_receptor_blob_object_to_table(db_name,blob_obj,rec_nbr,identifier,model_name+'.tar.gz',description,pdb_file_path)
        
        print(colored("Receptor object stored SUCCESSFULLY.","green"))
            
    except:
        print(colored("The table 'receptor_models' does not exists within the database. Create it first.","red"))
        exit()

def restore_receptor_manual(receptor_models_registry_path,to_path):
    """
    This function will restore a saved receptor model to the specified path.
    
    The function will ask for the 'rec_nbr' as the receptor model ID, which will be looked up in the 'receptor_models' table.
    
    ------
    Parameters:
    ------
    - db_name: the full path to the database file in which a table named 'receptors_models' exists.
    - rec_nbr: the integer value pointing to the receptor models and that is present in the 'rec_nbr' column of the 'receptors_models' table.
    """
    db_name = f'{receptor_models_registry_path}/receptors_models_registry.db'
    try: 
        
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
    except Exception as error:
        print(error)
        print(colored("Error connecting to receptors registry database","red"))
    
    # Indicate the receptor number
    while True:
        try:
            rec_nbr = int(input("Enter the receptor number to restore: "))
            # Check if 'rec_nbr' already exists in the 'receptors_models' table
            tag = acts_dbs.check_receptor_blob_item_index_presence(db_name,rec_nbr)
            if tag == 1:
                print(colored(f"Successfully retrieved rec_nbr {rec_nbr}.","green"))
                break
            else:
                print(colored("The receptor number entered does not exists. Check the input.","red"))
                continue
        except:
            print(colored("You need to enter an integer value","red"))

    #blob_ops.retrieve_blob_object(db_name,'receptor_models','rec_model_id',2,3,rec_nbr,to_path)
    blob_ops.retrieve_blob_object(db_name,'receptor_models','rec_model_id',2,rec_nbr,to_path)
    
def restore_receptor_acc(receptor_models_registry_path,to_path,rec_nbr):
    """
    This function will restore a saved receptor model to the specified path.
    
    The function will ask for the 'rec_nbr' as the receptor model ID, which will be looked up in the 'receptor_models' table.
    
    ------
    Parameters:
    ------
    - db_name: the full path to the database file in which a table named 'receptors_models' exists.
    - rec_nbr: the integer value pointing to the receptor models and that is present in the 'rec_nbr' column of the 'receptors_models' table.
    """
    db_name = f'{receptor_models_registry_path}/receptors_models_registry.db'
    try: 
        
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
    except Exception as error:
        print(error)
        print(colored("Error connecting to receptors registry database","red"))

    file_full_path = blob_ops.retrieve_blob_object(db_name,'receptor_models','rec_model_id',2,rec_nbr,to_path)
    
    return file_full_path

###############
    

if __name__ == "__main__":

    print("This is a local test")
        
    #create_receptors_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/receptors_models/receptors_models.db")
    
    #store_receptor("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/receptors_models/receptors_models.db","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_management/model_1")

    #restore_receptor("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects/sample.db",'/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/receptor_management/restore_rec')