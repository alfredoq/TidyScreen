import sqlite3
import xtarfile as tarfile
from termcolor import colored
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)


def create_blob_object(filename):
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

def append_receptor_blob_object_to_table(db_name,blob_obj,rec_nbr,identifier,model_name,description,pdb_file_path):
    """
    This function will store a BLOB object to a receptor table named 'receptors_models'.
    
    The following information will be stored:
    
        - col_1: rec_nbr - Type: int
        - col_2: identifier - Type: text
        - col_3: tar file - Type: blob
        - col_4: description - Type: text
        - col 5: pdb_file_path - Type: text
    ------
    Parameters:
    ------
    - db_name: the database in which the BLOB object is to be stored.
    - blob_obj: the BLOB object to be stored within the 'receptors_models' table.
    - rec_nbr: a unique integer value identifying the receptor model.
    - identifier: a string identification of the receptor model.
    - description: a brief description of the features of the receptor model.
    - pdb_file_path: the full path to the pdb file that originates the receptor and grid files.
    
    ------
    Returns
    ------
    All the info including the blob object is appended to the 'receptors_models' table.
   
    """
        
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    sql = f"""INSERT INTO receptor_models (rec_model_id, pdb_identifier, receptor_object, description, pdb_file)
              VALUES (?,?,?,?,?);"""             
    data_tuple = (rec_nbr, identifier, blob_obj, description, pdb_file_path)
    
    try: 
        cursor.execute(sql,data_tuple)
        conn.commit()
        cursor.close()
        print(colored("BLOB data added successfully","green"))
    except sqlite3.Error as error:
        print(error)
        print("Failed to insert blob data into sqlite table", colored(f"{error}","red"))
        
def append_ligand_blob_object_to_table(db_name,table_name,filename):
    """
    This function will store a BLOB object of a ligand .pdbqt file to a given table.
    
    The following information will be stored:
    
        - col_1: smiles - Type: TEXT
        - col_2: inchi_key - Type: TEXT
        - col_3: .pdbqt file - Type: blob
        - col_4: filename - Type: TEXT
    ------
    Parameters:
    ------
    - db_name: the full path to the database in which the BLOB object is to be stored.
    - table_name: the name of the table within the database in which the .pdbqt BLOB object is to be stored.
    - blob_obj: the BLOB .pdbqt file object to be stored within the table.
    - file: The filename of the object to be stored within the database.
    
    ------
    Returns
    ------
    All the info including the blob object is appended to the corresponding ligands database within the table.
   
    """
        
    file_blob_object = create_blob_object(filename)
    
    print(file_blob_object)
    
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    sql = f"""INSERT INTO {table_name} (pdbqt_file)
              VALUES (?);"""             
    data_tuple = (file_blob_object)
    
    try: 
        cursor.execute(sql,data_tuple)
        conn.commit()
        cursor.close()
        print(colored("BLOB data added successfully","green"))
    except sqlite3.Error as error:
        print(error)
        print("Failed to insert blob data into sqlite table", colored(f"{error}","red"))
        
    cursor.close()
        
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

def retrieve_blob_object(db_name,table_name,id_col_name,blob_item_position,id_value,dest_dir):
    """
    This function will retrieve a blob object that is stored within the DB, an outputs its again to the destination directory
    ------
    Parameters:
    ------
    - db_name: the full path the database containing the corresponding DB storing the BLOB object to be retrieved.
    - table_name: the name of the table containing the BLOB object to be retrieved.
    - id_col_name: the name of the column in which the 'id_value' is going to be searched.
    - blob_item_position: the index were the BLOB item is located within the cursor search.
    - filename_position: the index in which the stored filename is located within the cursos search.
    - filename_col: the name of the column in which the name of the stored file is indicated.
    - id_value: the value used as identifier of the BLOB object and that is associated to a field name used as identifier.
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
             WHERE {id_col_name} = '{id_value}'
          """
    
    cursor.execute(sql)
    query_output = cursor.fetchall()
    
    for row in query_output:
        data = row[blob_item_position]
        file_full_path = f'{dest_dir}/rec_model_id_{id_value}.tar.gz'
        #print(file_full_path)
        try:
            write_blob_object(data,file_full_path)
            #tar = tarfile.open(data,file_full_path)
            #tar.extract(dest_dir)
            print(colored(f"SUCCESSFULLY retrieved the receptor BLOB object to: \n {dest_dir}","green"))
            return file_full_path
        
        except Exception as error:
            print(error)
            print(colored("ERROR retrieving the receptor BLOB object.","red"))


    #write_blob_object(data,file_full_path)
    #print("BLOB object retrieved SUCCESSFULLY")


##################3
    
    
    
if __name__ == "__main__":

    print("This is a local test")
    
    #append_ligand_blob_object_to_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/sample.db", "alpha_aminoacid_pdbqt_ligands","/tmp/RYSXZQIOOACQRG-UHFFFAOYSA-N.pdbqt")
    
    #create_blob_object("/tmp/RYSXZQIOOACQRG-UHFFFAOYSA-N.pdbqt")
    