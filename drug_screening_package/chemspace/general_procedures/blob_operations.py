import sqlite3 
from termcolor import colored
import pandas as pd
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

def create_blob_object(db_name,table_name):
    '''
    This function will create from zero a table specifically used to store BLOB (Binary Large OBject) items, such as files, images, etc. 
    
    The new Table will contain 3 fields, as follows:
    
        - col_1: object_identifier - Type: text
        - col_2: blob_item - Type: blob
        - col_3: description - Type: text
        
    The function will initially check if the table to be created exists within the given database. If TRUE, it will exit with a message, otherwise, will create the corresponding table.
    
    ------
    Parameters:
    ------
    - db_name: the full path to the database in which the table storing the BLOB object is to be created
    - table_name: the name of the table that will be inserted in db_name for storing the BLOB objects.    
    
    ------
    Returns
    ------
    A table destinated to store files as BLOB object is created in the specified database.
    '''
    conn = sqlite3.connect(db_name)
    
    try: # This will end the execution if the table is found
        sql_1 = f"""SELECT * 
                    FROM {table_name}"""
        conn.execute(sql_1)
        print('Table already exists, please ', colored("append the object.","red"))
    
    except: # This will proceed with the table creation in case the table is not found
        print(colored("Creating ","green") + "new BLOB object table")
        sql_2 = f"""CREATE TABLE {table_name} 
                    (obj_id TEXT NOT NULL, blob_item BLOB NOT NULL, file TEXT NOT NULL, description TEXT NOT NULL);"""      
        conn.execute(sql_2)

def read_blob_object(filename):
    """
    This function will create a binary object from a filename provided that is aimed to be stored as a BLOB object within the database.
    ------
    Parameters:
    ------
    - filename: the filename in the corresponding format to be stored in the table containing BLOB objects.
    """
    
    with open(filename, 'rb') as file:
        blobData = file.read()
    
    file.close()
    
    return blobData

def write_blob_object(data,dest_file):
    """
    This function will create a binary object from a filename provided that is aimed to be stored as a BLOB object within the database.
    ------
    Parameters:
    ------
    - data: the information usually associated to a BLOB object
    - dest_file: the full path to the file of the destination file to be retrieved.
    """
    
    try: 
    
        with open(f'{dest_file}', 'wb') as file:
            file.write(data)
        file.close()
        print(f"File {dest_file} stored SUCCESSFULLY")
    
    except:
        print(f"Problem storing {dest_file}. Check for errors.")

def check_id_in_table(db_name,table,field_to_be_searched,value):
    """Will check if the provided obj_id tag already exists in the 'field_to_be_searched' and returns a key value to proceed with a corresponding action.
    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the table to be searched.
    - table: the name of the table containing the field to be compared 
    - field_to_be_searched: the name of the filed which will be searched for the presence of the provided value.
    - value: the string to be searched in the corresponding field, and for which a matching will be evaluated.
    """
    conn = sqlite3.connect(db_name)
    sql = f"""SELECT {field_to_be_searched}
              FROM {table};"""

    query_list = list(pd.read_sql_query(sql,conn)[field_to_be_searched])
    
    if value in query_list:
        print('The entered ID is', colored('already present in the database','red'))
        return 0 
    else:
        print('The entered ID is', colored('is new, adding to DB','green'))
        return 1

def append_blob_object_to_table(db_name,table_name,filename_to_blob):
    """
    This function will store a BLOB object to an already existing given table, which contains the following structure:
    
        - col_1: object_identifier - Type: text
        - col_2: blob_item - Type: blob
        - col_3: description - Type: text
    ------
    Parameters:
    ------
    - db_name: the database in which the BLOB object is to be stored.
    - table_name: the name of the table that will hold the BLOB object, and which complies with the above mentioned fields structure.
    - blob_object: the full path to the blob_object to be stored.

    ------
    Returns
    ------
    A file is written as a BLOB object into the database.
    """
    object_identifier = input("Provide an object identifier: ")
    description = input("Provide brief object description: ")

    filename = filename_to_blob.split('/')[-1]

    # check if the 'object_identifier' already exists in the target table
    check_value = check_id_in_table(db_name,table_name,"obj_id",object_identifier)
    
    if check_value == 0:
        print('The entered ID is', colored('already present in the database','red'))
        exit()
    
    blob_object = read_blob_object(filename_to_blob)
    
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    sql = f"""INSERT INTO {table_name} (obj_id, blob_item, file, description)
              VALUES (?,?,?,?);"""             
    data_tuple = (object_identifier, blob_object, filename, description)
    
    try: 
        cursor.execute(sql,data_tuple)
        conn.commit()
        print("BLOB data added successfully")
        cursor.close()
    except sqlite3.Error as error:
        print("Failed to insert blob data into sqlite table", colored(f"{error}","red"))

def retrieve_blob_object(db_name,table_name,obj_id,dest_dir):
    """
    This function will retrieve a blob object that is stored within the DB, an outputs its again to the destination directory
    ------
    Parameters:
    ------
    - db_name: the full path the database containing the corresponding DB storing the BLOB object to be retrieved.
    - table_name: the name of the table containing the BLOB object to be retrieved.
    - obj_id: the value used as identifier of the BLOB object and that is associated to a field name: 'obj_id'
    - dest_dir: the full path to the location were the store file is to be written.

    ------
    Returns
    ------
    A file is read from the database and written to the disk for further use.
    """

    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    sql = f"""SELECT * 
             FROM {table_name}
             WHERE obj_id = '{obj_id}'
          """
    cursor.execute(sql)
    query_output = cursor.fetchall()
    
    for row in query_output:
        data = row[1]
        file_full_path = f'{dest_dir}/{row[2]}'

    write_blob_object(data,file_full_path)
    print("BLOB object retrieved SUCCESSFULLY")


###############

if __name__ == "__main__":
    print("This is a local test")
    
    #create_blob_objects_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","blob_objs")
    
    #check_id_in_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","blob_objs","obj_id","a")
    
    #append_blob_object_to_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","blob_objs","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/platform_diagram.drawio.jpg")
    
    #retrieve_blob_object("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","blob_objs","diagram","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/blob")
    
    read_blob_object("/tmp/RYSXZQIOOACQRG-UHFFFAOYSA-N.pdbqt")