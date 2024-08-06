import sqlite3
import pandas as pd
from termcolor import colored

def create_md_registry_table(general_md_registries_folder):
    """
    This function will create a table designed to store information related to performed molecular dynamics campaigns.

    The table will contain the following columns:

        - assay_id: an integer values to store a specific docking assay;
        - db_name: the full path to the database containing the table in which the docked molecule to be simulated is located ;
        - lig_table_name: the name of the table located in 'db_name' that contains the ligands pose to be simulated.
        - ligand_subpose_name: the ligand subpose name as obtrained from molecular docking that is going to be subjected to MD studies
        - rec_model_id: an integer indicating the receptor model used in the docking studies and that is the template for the MD simulation;
        - md_params_id: the set of parameters affecting the simulation, and that is referenced to the corresponding registry database.
        - stored_results: a flag indicating if the md assay results are available in the assay table (0:no, 1:yes)
        - assay_folder: the path the folder were all the info to run the assay was created.
        - assay_description: a brief description of the objective of the specific docking assay.

    ATTENTION: if the table named 'md_registries' is stored within the a database '{general_md_registries_folder/md_registries.db} it will be replaced.

    The function should not be called again if the database/table has already been created, unless you want to clear the overall table of registries.
        
    ------
    Parameters:
    ------
    - general_md_registries_folder: the full path were the molecular dynamics assays registry table is to be stored.

    ------
    Returns:
    ------
    Creates an empty registry table. 
    """
    
    try: 
        conn = sqlite3.connect(f'{general_md_registries_folder}/md_registries.db')

        ins = f'CREATE TABLE md_registries ("assay_id" INTEGER, "db_name" TEXT, "lig_table_name" TEXT, "lig_subpose_name" TEXT, "rec_model_id" INTEGER,"md_params_id" INTEGER, "stored_results" INTEGER, "assay_folder" TEXT, assay_description TEXT)'
        conn.execute(ins)
        print(colored("Molecular dynamics registers table created.","green"))
        
    except Exception as error:
        print(error)
        print(colored("Error creating the database docking_registries.db","red"))

def get_last_md_assay_id(md_assays_registry_path):
    """
    This function will search for the last assay_id value in the 'md_registries' table of the 'md_registries.db' database under the path {general_md_registries_folder}.
    
    In case no registry value is found, '1' will be returned, assuming that the table is new or has been intentionally deleted.
    
    ------
    Parameters:
    ------
    - general_md_registries_folder: the full path to the folder in which the molecular dynamics registries are being stored.
    
    ------
    Returns:
    ------
    - last_value = The last incremental value stored as 'assay_id' for the docking runs.
    """
    
    try:
        conn = sqlite3.connect(md_assays_registry_path+'/md_registries.db')
        sql_assay = """SELECT assay_id 
                       FROM md_registries;"""
    
        df = pd.read_sql_query(sql_assay,conn)
        assay_id_list = df["assay_id"].to_list()
        if len(assay_id_list) == 0: # This will hit in case no previous record in the table is present
            return 0
        else:
            last_value = assay_id_list[-1]
            return last_value
    
    except Exception as error:
        print(error)
        print(colored("Error retrieving the last MD assay_id value. Check if registries table already exists","red"))

def append_md_registry(general_registries_folder,db_name, table_name, ligand_subpose_name, rec_model, stored_results, assay_folder,description,md_params_id):
    """
    When preparing an md assay, this function will append a registry with the required information into the registries table.
    
    ------
    Parameters:
    ------
    - general_registries_folder: the full path the folder in which the general md registries are being stored.
    - db_name: the full path to the database containing the ligand to be docked.
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
        conn = sqlite3.connect(general_registries_folder+'/md_registries.db')
        sql_assay = """SELECT assay_id 
                       FROM md_registries;"""
    
        # This will retrieve the last value of the assay_id in the registry table and increment by 1 for the next storage
        last_value = get_last_md_assay_id(general_registries_folder)
        assay_id = last_value+1
    
        # Get all the values present in the table
        df = pd.read_sql_query(sql_assay,conn)
        assay_id_list = df["assay_id"].to_list()
        
        # Check if assay_id exists to prevent as a second factor repeating the id. Should never match because assay_is is automatically generated
        if assay_id in assay_id_list:
            print(colored(f"The md assay_id {assay_id} already exists in the registry table. EXITING","red"))

            exit 
        else:
            raise # This will manually raise an exception to continue
    except:
        sql = f"""INSERT INTO md_registries (assay_id,db_name,lig_table_name,lig_subpose_name,rec_model_id,md_params_id,stored_results,assay_folder,assay_description) VALUES({assay_id},"{db_name}","{table_name}","{ligand_subpose_name}",{rec_model},{md_params_id},{stored_results},"{assay_folder}","{description}");"""
        conn.execute(sql)
        conn.commit()
        print(colored(f"MD assay_id {assay_id} is NEW, stored the registry","green"))

