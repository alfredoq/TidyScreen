import sqlite3
from tqdm import tqdm # import tqdm
from pandarallel import pandarallel # import pandarallel
from rdkit.Chem import rdmolops
from rdkit import Chem
from rdkit import RDLogger
from termcolor import colored
RDLogger.DisableLog('rdApp.*')
import pandas as pd
import time
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.general_procedures.smiles_processing as general_proc
import chemspace.general_procedures.operations_on_sql_database as sql_ops
import chemspace.input_output_operations.i_o_checkings as i_o_check

"""
This module contains all the functions to process a raw Emolecule Database File, and the
objective is to create an initial SQL database containing 2 columns:

- SMILES: the smiles notation of each reactant, with all stereo a cofactors removed;
- InchiKey: the value as computed by the RDKit package.

It is expected that upon comparison of InchiKey values, each rectant in the database
is unique.
"""

def process_emolecules_raw_file(file,db_name,passed_reactants_table,repeated_reactants_table,error_reactants_table,action):    
    clean_smiles_list = []
    clean_inchikey_list = []
    repeated_smiles_list = []
    repeated_inchikey_list = []
    smiles_error_list = []
    
    counter = 0
    with open(file,'r') as input_file:
        num_lines = sum(1 for _ in input_file) - 2 # Get the number of lines in file
        input_file.seek(0) # Go to the start of the line
        next(input_file) # Will skip the first line (header)       
        actual_time = time.time()
        for line in input_file:               
            counter+=1
            new_line = line.split()
            smiles = new_line[0]
            # check first pair of lists
            clean_smiles_list,clean_inchikey_list,actual_time=i_o_check.check_lists_size(clean_smiles_list,clean_inchikey_list,50000,passed_reactants_table,num_lines,counter,actual_time,db_name)
            # check second pair of lists
            repeated_smiles_list,repeated_inchikey_list,actual_time=i_o_check.check_lists_size(repeated_smiles_list,repeated_inchikey_list,50000,repeated_reactants_table,num_lines,counter,actual_time,db_name)
            try:     
                clean_smiles, inchi_key = general_proc.generate_clean_smiles(smiles)

                if inchi_key not in clean_inchikey_list:
                    clean_smiles_list.append(clean_smiles)
                    clean_inchikey_list.append(inchi_key)
                else:
                    repeated_smiles_list.append(clean_smiles)
                    repeated_inchikey_list.append(inchi_key)
            except:
                smiles_error_list.append(smiles)
                continue
    input_file.close()
    
    # Join the data into a dataframe for further storing into database
    df = i_o_check.merge_2_lists_into_df(clean_smiles_list,clean_inchikey_list,"SMILES","Inchi_Key")
    df_repeated = i_o_check.merge_2_lists_into_df(repeated_smiles_list,repeated_inchikey_list,"SMILES","Inchi_Key")
    df_errors = i_o_check.convert_1_list_into_df(smiles_error_list,"SMILES")
    # This will store the remaining of data not affected by the block backup
    sql_ops.store_df_to_sql(db_name,df,passed_reactants_table,action)
    sql_ops.store_df_to_sql(db_name,df_repeated,repeated_reactants_table,action)
    sql_ops.store_df_to_sql(db_name,df_errors,error_reactants_table,action)

def process_emolecules_raw_file_tqdm(file,db_name,passed_reactants_table,action):
    
    raw_df = pd.read_csv(file,sep=' ')
    
    tqdm.pandas() # initialize tqdm for pandas
        
    raw_df["SMILES"] = raw_df["isosmiles"].progress_apply(lambda x: general_proc.generate_clean_smiles_2(x))

    raw_df["inchi_key"] = raw_df["SMILES"].progress_apply(lambda x: general_proc.compute_inchi_key(x))

    #raw_df.drop_duplicates('inchi_key',keep='last')
    
    clean_df = raw_df[["SMILES","inchi_key"]]

    sql_ops.store_df_to_sql(db_name,clean_df,passed_reactants_table,action)

def process_emolecules_raw_file_pandarallel(file,db_name):
    """
    This function will read a raw 'emolecules' files containing available reactants stored in the SMILES format. 
    
    The file needs to contain header values indicating the name of each files, which are separated by a 'space'. 

    It is important to note that the original 'emolecules' file contains the SMILES strings in a field originally containing the header: 'isosmiles', and this keyword will be retrieved to perform the file parsing.

    This function will make use of a dataframe feature implemented in the 'pandarallel' library in order to split into multiple processors the parsing of the information. Consquently, when using this function all the available CPU cores will be assigned to process the information.

    If the database file already exists, the table 'emolecules_reagents' will be replaced.

    --- Parametes ---
    - file: the full path to the raw file containing the 'emolecules' raw information
    - db_name: this full path and filename of the database containing the parsed information.

    --- Returns ---
    A database containing the corresponding information.
    """
    # Read the raw datafile into a dataframe
    raw_df = pd.read_csv(file,sep=' ')
    
    # initialize pandarallel for pandas and perform some SMILES processing
    pandarallel.initialize(progress_bar=True)
    print(colored("Parsing SMILES string.","green"))
    raw_df["SMILES"] = raw_df["isosmiles"].parallel_apply(lambda x: general_proc.generate_clean_smiles_2(x))
    print(colored("Computing inchi_key values.","green"))
    raw_df["inchi_key"] = raw_df["SMILES"].parallel_apply(lambda x: general_proc.compute_inchi_key(x))
    print(colored("Droping duplicated molecules.","green"))
    raw_df.drop_duplicates('inchi_key',keep='last')
    
    clean_df = raw_df[["SMILES","inchi_key"]]

    # Create a mol_id index and make it the first column
    clean_df["mol_id"] = clean_df.index
    first_column = clean_df.pop('mol_id')
    clean_df.insert(0, 'mol_id', first_column)

    # Store the final info to the corresponding database
    sql_ops.store_df_to_sql(db_name,clean_df,"emolecules_reagents","replace")

def process_raw_csv_file_pandarallel(file,db_name,table_name,retain_stereo):
    """
    This function will read a raw csv file containing SMILES molecules available in the first field. 

    The file needs to contain header values indicating the name of each files, which are separated by a ',' value. 

    It is important to note that the original csv file contains the SMILES strings in a field originally containing the header: 'SMILES', since this keyword will be retrieved to perform the file parsing.

    This function will make use of a dataframe feature implemented in the 'pandarallel' library in order to split into multiple processors the parsing of the information. Consquently, when using this function all the available CPU cores will be assigned to process the information.

    If the database file already exists, the table 'emolecules_reagents' will be replaced.

    --- Parametes ---
    - file: the full path to the raw csv file containing the raw information
    - db_name: this full path and filename of the database containing the parsed information.
    - table_name: the name of the table in which the parsed molecules will be stored. If the table already exists, molecules will be appended.
    - retain_stereo: flag (0/1: False:True) to retaing the stereochemistry explicitly indicated in the SMILES input.

    --- Returns ---
    A database containing the corresponding information.
    """
    # Read the raw datafile into a dataframe
    print(file)
    
    try: 
        raw_df = pd.read_csv(file)
        # initialize pandarallel for pandas and perform some SMILES processing
        pandarallel.initialize(progress_bar=True) 
        print(colored("Parsing SMILES string.","green"))
        if retain_stereo == 0: # Flag to delete existing stereochemistry
            print("COMIENZO")
            raw_df["SMILES"] = raw_df["SMILES"].parallel_apply(lambda x: general_proc.generate_clean_smiles_2(x))
            print(colored(" Deleting existing stereochemistry. Computing inchi_key values.","green"))
        if retain_stereo == 1: # Flag to retain existing stereochemistry
            raw_df["SMILES"] = raw_df["SMILES"].parallel_apply(lambda x: general_proc.generate_clean_smiles_2_retain_stereo(x))
            print(colored(" Mantaining existing stereochemistry. Computing inchi_key values.","green"))
        elif retain_stereo != 0 and retain_stereo != 1:
            print(retain_stereo)
            print(colored("The flag to manage existing stereochistetry is incorrect. EXITING","red"))
            exit()
        
        raw_df["inchi_key"] = raw_df["SMILES"].parallel_apply(lambda x: general_proc.compute_inchi_key(x))
        
        print(colored("Droping duplicated molecules.","green"))
        
        
        raw_df.drop_duplicates('inchi_key',keep='last')
        
        clean_df = raw_df[["SMILES","inchi_key"]]

        # Create a mol_id index and make it the first column
        clean_df["mol_id"] = clean_df.index
        first_column = clean_df.pop('mol_id')
        clean_df.insert(0, 'mol_id', first_column)

        # Store the final info to the corresponding database
        print(colored(f"Storing molecules to database {db_name}.","green"))
        sql_ops.store_df_to_sql(db_name,clean_df,table_name,"append")
        print(colored(f"Removing duplicated molecules in {db_name}.","green"))
        sql_ops.purge_duplicates_in_table(db_name,table_name,"inchi_key")
        
        print(colored(f"SUCCESSFULLY FINISHED PROCESING THE FILE {file}.","green"))

    except Exception as error:
        print(error)
        print(colored(f"ERROR importing the file {file}"))

def add_reactant_filter_to_db(db_name):
    """
    This function will include a molecule filter based on the SMARTS nomenclature into the database file, under the table: 'reactants_types'.

    Each filter will be appended to preexisting filters, and upon execution, the following information will be asked to the user:

        - Reactant: this is the name/identification of the reactant that matches the provided SMARTS.
        - SMARTS: the corresponding pattern that matches the desired molecules.
        - User: the name of the person including the corresponding filter.

    --- Parameters ---
    - db_name: the full path to the already created databsase in which the SMART filter for molecules is to be stored.
    
    --- Returns ---
    A table named 'reactants_types' containing the provided information.
    """
    reactant_name = input("Name the family of reactants: ")
    reactant_smarts = input("Reactant SMARTS: ")
    user = input("User who added the SMARTS?: ")
    dictionary = {"Reactant" : [reactant_name],
                  "SMARTS" : [reactant_smarts],
                  "User" : [user],
                }
    try: 
        df = pd.DataFrame.from_dict(dictionary)
        sql_ops.store_df_to_sql(db_name,df,"reactants_types","append")
        print("Added reactant succesfully")
    except:
        print("Error while adding reactant")
        
def add_reaction_smarts_to_db(db_name):
    """"
    This function will add a reaction modeled in SMARTS notation to the provided database. 
    
    The function will query for: 
    
        - Name of the reaction:
        - Reaction SMARTS: 
        - User who provided the reaction: 
        
    All this information will be stored in a table called 'reaction_types'.
    
    ------
    Parameters:
    ------
    - db_name: the full path and name of the database in which the reaction SMARTS will be inserted.

    ------
    Returns
    ------
    The provided information will be stored within the database in a table called 'reaction_types'.
    """

    reaction_name = input("Name the the reaction: ")
    reaction_smarts = input("Reaction SMARTS: ")
    user = input("User who added the SMARTS?: ")
    dictionary = {"Reaction" : [reaction_name],
                  "SMARTS" : [reaction_smarts],
                  "User" : [user],
                }

    df = pd.DataFrame.from_dict(dictionary)
    sql_ops.store_df_to_sql(db_name,df,"reaction_types","append")
    print(f"Reaction {reaction_name} added successfully")

def create_image_from_reactants(db_name,table):
    """
    This function will create a set of images of a reactants tables, with the number of plotted images being configured by the passed arguments.
    """
    conn = sqlite3.connect(db_name)
    query_df = pd.read_sql_query(f'SELECT * FROM {table}',conn,index_col=None)
    
    for index, inchi_key in query_df.iterrows():
        print(inchi_key["inchi_key"])
        
# Use to test the functions
if __name__ == '__main__':
    
    print("Launching local tests")
    #process_emolecules_raw_file(file,db_name,"reagents_pass","reagents_dup","reagents_err","append")

    #add_reactant_filter_to_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db")

    #process_emolecules_raw_file_tqdm("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/raw_data/sample_1.smi","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_2.db","reagents_pass","append")

    #process_emolecules_raw_file_pandarallel("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/raw_data/emolecules_database.smi","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_2.db","reagents_pass","append")

    #add_reactant_filter_to_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db")

    #add_reactant_filter_to_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db")

    #add_reaction_smarts_to_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db")

    #add_reaction_smarts_to_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db")
    
