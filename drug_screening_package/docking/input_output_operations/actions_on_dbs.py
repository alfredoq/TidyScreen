import sqlite3
import pandas as pd
from termcolor import colored
import tarfile
from glob import glob
import shutil
import os
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

def check_receptor_blob_item_index_presence(db_name,rec_nbr):
    """
    This function will check if the 'rec_nbr' integer value entered to store an item in the 'receptors_models' table exist.
    
    If the 'rec_nbr' exists, it will return '1', otherwise it will return '0'.
    
    ------
    Parameters
    ------
    - db_name: the full path to the database files in which a 'receptors_models' table exists, and contains a column named 'rec_nbr'.
    - rec_nbr: the integer to be searched in the corresponding field.
    
    ------
    Returns
    ------
    Will return '1' if 'rec_nbr' already exists, otherwise it will return '0'.
    """
    
    conn = sqlite3.connect(db_name)
    
    sql = """SELECT rec_model_id
              FROM receptor_models;"""
              
    try:
        df = pd.read_sql_query(sql,conn)
        list_of_rec_nbr = df['rec_model_id'].tolist()
        if rec_nbr in list_of_rec_nbr:
            return 1
        else:
            return 0
    except Exception as error:
        print(error)
        print(colored("ERROR searching the 'rec_model_id' value"))
        exit()
        
        
def uncompress_receptor_files(file,target_folder):
    """
    This function is used when a receptor BLOB object is retrieved from a database, and the intention is to uncompress the corresponding file and also extract all the corresponding file directly into the 'receptors' folder.
    
    ------
    Parameters:
    ------
    - file: the full path to the file to be extracted.
    - target_folder: the place in which the files will be extracted.
    
    ------
    Returns:
    ------
    The whole set of compressed file within the path of input file
    """
    
    try:
        ## This will uncompress the receptor in the given location
        file = tarfile.open(file)
        file.extractall(f'{target_folder}')
        # This list comprehension will return the corresponding directory
        dir = [x for x in glob(f"{target_folder}/*") if os.path.isdir(x)]
        all_files = os.listdir(dir[0])
        for file in all_files:
            shutil.move(f'{dir[0]}/{file}',f'{target_folder}/{file}')
        # Delete the old tree file    
        shutil.rmtree(dir[0])
    except Exception as error:
        print(error)
        print(colored("Problem extraction the receptor object at funct 'uncompress_receptor_files()'","red"))
        

###############
    

if __name__ == "__main__":

    print("This is a local test")
    
    uncompress_receptor_files("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays/docking_assay_10/receptor/rec_model_id_3.tar.gz","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays/docking_assay_10/receptor")
    
    