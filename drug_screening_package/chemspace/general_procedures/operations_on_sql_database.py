
import sqlite3
import re
from tqdm import tqdm
from pandarallel import pandarallel
import readline # to improve the performance of input()
import os
from termcolor import colored
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.input_output_operations.i_o_checkings as i_o_check
import chemspace.general_procedures.smiles_processing as smiles_proc
import chemspace.general_procedures.operations_on_sql_database as db_ops
import chemspace.input_output_operations.actions_on_dbs as db_actions
import chemspace.general_procedures.smiles_processing as general_proc

def store_df_to_sql(db_name,df,table_name,action):
    """
    This function will store a dataframe that is provided within an sqlite3 database.
    -----
    Input Params:
    -----
    - db_name: the name of the database in which to store the table. It is does not exists, will create it, otherwise the table will be included among the rest of tables in the schema.
    - df: the dataframe to be stored as a table within the database.
    - table: the name of the table in which to drop the dataframe.
    - action: In case the 'table' already exists, how to manage the data. Options: "append", "replace".
    - value: the target value that will match the record to be droped.
    """
    
    connection = sqlite3.connect(f'{db_name}')
    df.to_sql(con=connection, name=table_name,if_exists=action,index=None)

def delete_record(db_name,table,column,value):
    """
    This function will drop a record included in a specific database. First, the function will check if the searched value (target to drop) exists, and will inform a message in case it doesnt.
    -----
    Input Params:
    -----
    - db_name: the full path and name of the sqlite database to be searched and modified.
    - table: the name of the table in which a value included in a record is to be droped.
    - column: the field in which to search the target to drop.
    - value: the target value that will match the record to be droped.
    """
    
    try: 
        connection = sqlite3.connect(db_name)
        sql = f"DELETE FROM {table} WHERE {column} == '{value}';"
        #sql = f"SELECT * FROM {table} WHERE {column} == '{value}';"
        print(sql)
        connection.execute(sql)
        #print(query)
        connection.commit()
        print("The record has been deleted successfully!")
    
    except:
        print("For some reason the record could not be deleted!")
    
    connection.close() 

def select_reactants_based_on_smarts(db_name,reactants_types_table,reactant,table_to_filter):
    """"
    This function will generate and store a table to the given database, that is populated with the 'SMILES' and 'inchi_key' corresponding to the reactants included in the 'table_to_filter' filter that matches the 'SMARTS' notation stored in the 'reactants_types_table' for the corresponding 'reactant'
    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the information to be processed.
    - reactants_types_table: the name of the table in the 'db_name' containig the type of reactants and the association with the corresponding SMARTS notation. A field named 'smarts' must exist in this table.
    - reactant: this is the name of the type of molecule (reactant) that is to be searched based on the SMARTS notation.
    - table_to_filter: the name of the table containing all the SMILES notation of the molecules to be compared with the SMARTS
    """
    
    database_path = '/'.join(db_name.split('/')[:-1])
    conn =  sqlite3.connect(f'{db_name}')
    sql_search_op = f'SELECT SMARTS FROM {reactants_types_table} WHERE Reactant == "{reactant}";'
    query_result = conn.execute(sql_search_op)
    target_smarts = query_result.fetchone()[0]
    
    cursor = conn.cursor()
    inchi_key_hit_list = []
    smiles_key_hit_list = []
    
    counter = 0
    for row in cursor.execute(f'SELECT SMILES, Inchi_Key FROM {table_to_filter};'):
        smiles = row[0]
        inchi_key = row[1]
        
        hit_key = check_smiles_vs_smarts(smiles, target_smarts)
        
        if hit_key == 1:
            if inchi_key not in inchi_key_hit_list:
                counter+=1
                print(f"Search HIT: {counter}")
                inchi_key_hit_list.append(inchi_key)
                smiles_key_hit_list.append(smiles)
            
                # This will trigger the write to the DB in case the list exceeds the allowed limit size        
                if len(inchi_key_hit_list) > 10000:    
                    df =i_o_check.merge_2_lists_into_df(smiles_key_hit_list,inchi_key_hit_list,"SMILES","Inchi_Key")
                    smiles_key_hit_list = [] # Will flush the list
                    inchi_key_hit_list = [] # Will flush the lis
                    # Will dump the dataframe to a file for further storing in the database
                    if not os.path.isfile(f'{database_path}/temp.csv'):
                        df.to_csv(f'{database_path}/temp.csv', header=True, index=None)
                    else:
                        df.to_csv(f'{database_path}/temp.csv', mode='a',header=False, index=None)
                
    cursor.close() # Close the cursor for the final writting operation
    df = i_o_check.merge_2_lists_into_df(smiles_key_hit_list,inchi_key_hit_list,"SMILES","inchi_key")
    df_previous = pd.read_csv(f'{database_path}/temp.csv')
    store_df_to_sql(db_name,df,reactant,"replace")
    store_df_to_sql(db_name,df_previous,reactant,"append")
    os.remove(f'{database_path}/temp.csv') # Will remove the temp file containing the data
    conn.close()

def select_reactants_based_on_smarts_tqdm(db_name,reactants_types_table,reactant,table_to_filter):
    """
    This function will generate and store a table to the given database, that is populated with the 'SMILES' and 'inchi_key' corresponding to the reactants included in the 'table_to_filter' filter that matches the 'SMARTS' notation stored in the 'reactants_types_table' for the corresponding 'reactant'. The procedure of filtering is based on the library 'tqdm'
    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the information to be processed.
    - reactants_types_table: the name of the table in the 'db_name' containig the type of reactants and the association with the corresponding SMARTS notation. A field named 'smarts' must exist in this table.
    - reactant: this is the name of the type of molecule (reactant) that is to be searched based on the SMARTS notation.
    - table_to_filter: the name of the table containing all the SMILES notation of the molecules to be compared with the SMARTS
    """
    conn =  sqlite3.connect(f'{db_name}')
    sql_search_op = f'SELECT SMARTS FROM {reactants_types_table} WHERE Reactant == "{reactant}";'
    query_result = conn.execute(sql_search_op)
    target_smarts = query_result.fetchone()[0]

    sql_search_op = f'SELECT * FROM {table_to_filter};'
    df_to_filter = pd.read_sql_query(sql_search_op,conn)

    tqdm.pandas() # Instatiate the tqdm package for posterior running.
    df_to_filter["hit_key"] = df_to_filter['SMILES'].progress_apply(lambda x: check_smiles_vs_smarts(x,target_smarts))

    # Subset the SMILES matching the SMARTS string
    df_hits = df_to_filter[df_to_filter["hit_key"] == 1].copy()
    # Drot the columns for hits tags
    df_hits.drop(columns=['hit_key'],axis=1,inplace=True)
    # Save the dataframe to the database
    store_df_to_sql(db_name,df_hits,reactant,"replace")

def select_reactants_based_on_smarts_pandarallel(db_name,table_to_filter,filter_name):
    """
    This function will search and find molecules stored a SMILES in a certain table within the provided database, with a matching criteria based on SMARTS notation.

    By default reactants SMARTS will be read form a table named 'reactants_types'
    
    All the matched molecules will be stored to table whose name is equal to the reactant keyword provided.

    In addition to the SMILES notation, the computed inchi_key will also be stored within the resulting table.

    This function is based on the 'pandarallel' library, so the whole search will be run as a multiprocessor event, ocuppying all available CPU cores. 

    If a previous filtering has been done, and consequently a table named after the searched molecule already exists within the database file, it will be replaced by the new results. 

    If duplicated reactants matching the smarts structure are found based on 'inchi_key' comparison, only one instance is retained.

    A column name 'available' is created. In the initial subseting the tag is set to 0. The availability is managed by a separated function.

    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the source table containing the SMILES notation which is to be filtered based on a SMARTS pattern.
    - table_to_filter: the name of the table containing all the SMILES notation of the molecules to be compared with the SMARTS
    - filter_name: this is the identifier of the type of molecules that is to be searched based on the SMARTS notation.
    
    ------
    Returns 
    ------
    A new table is created within the database and named after the filter_name name provided.
    """

    conn =  sqlite3.connect(f'{db_name}')
    sql_search_op = f'SELECT SMARTS FROM reactants_types WHERE Reactant == "{filter_name}";'
    query_result = conn.execute(sql_search_op)
    target_smarts = query_result.fetchone()[0]

    sql_search_op = f'SELECT * FROM {table_to_filter};'
    df_to_filter = pd.read_sql_query(sql_search_op,conn)

    pandarallel.initialize(progress_bar=True) # Instatiate the pandarallel package for posterior running.
    df_to_filter["hit_key"] = df_to_filter['SMILES'].parallel_apply(lambda x: check_smiles_vs_smarts(x,target_smarts))
    # Will tag initial compound availability as 0
    df_to_filter["available"] = df_to_filter['SMILES'].parallel_apply(lambda x: 0)

    # Subset the SMILES matching the SMARTS string
    df_hits = df_to_filter[df_to_filter["hit_key"] == 1].copy()
    # Drop the columns for hits tags
    df_hits.drop(columns=['hit_key'],axis=1,inplace=True)
    # Delete duplicates by inchi_key
    df_hits.drop_duplicates(subset='inchi_key', keep="last", inplace=True)
    # Save the dataframe to the database
    store_df_to_sql(db_name,df_hits,filter_name,"replace")

def exclude_reactants_based_on_smarts_pandarallel(db_name,table_to_filter,filter_name,output_table_name):
    """
    This function will search and EXCLUDE molecules stored as SMILES in a certain table within the provided database, with a matching criteria based on SMARTS notation.

    By default, the SMARTS of the reactant to exclude will be read from the table 'reactants_types'.
    
    All the matched molecules will be excluded from the input table and afterwards saved to a new table whose name is equal to the name of the SMARTS pattern provided.

    In addition to the SMILES notation, the computed inchi_key will also be stored within the resulting table.

    This function is based on the 'pandarallel' library, so the whole search will be run as a multiprocessor event, ocuppying all available CPU cores. 

    If a previous filtering has been done, and consequently a table named after the searched molecule already exists within the database file, it will be replaced by the new results. 
       
    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the source table containing the SMILES notation which is to be filtered based on a SMARTS pattern.
    - table_to_filter: the name of the table containing all the SMILES notation of the molecules to be compared with the SMARTS
    - filter_name: this is the identifier of the type of molecules that is to be searched based on the SMARTS notation.
    - output_table_name: the name of the table to which the compunds without excluded ones will be stored.
    
    ------
    Returns 
    ------
    A new table is created within the database and named after the 'output_table_name' parameter.
    """

    conn =  sqlite3.connect(f'{db_name}')
    sql_search_op = f'SELECT SMARTS FROM reactants_types WHERE Reactant == "{filter_name}";'
    query_result = conn.execute(sql_search_op)
    target_smarts = query_result.fetchone()[0]

    sql_search_op = f'SELECT * FROM {table_to_filter};'
    df_to_filter = pd.read_sql_query(sql_search_op,conn)

    pandarallel.initialize(progress_bar=True) # Instatiate the pandarallel package for posterior running.
    df_to_filter["hit_key"] = df_to_filter['SMILES'].parallel_apply(lambda x: check_smiles_vs_smarts(x,target_smarts))

    # Subset the SMILES matching the SMARTS string - Here only 'hit_key' == 0 is wanted because we need to EXCLUDE the compounds matching the target SMARTS notation.
    df_hits = df_to_filter[df_to_filter["hit_key"] == 0].copy()
    # Drop the columns for hits tags
    df_hits.drop(columns=['hit_key'],axis=1,inplace=True)
    # Save the dataframe to the database
    store_df_to_sql(db_name,df_hits,output_table_name,"replace")

def exclude_reactants_based_on_number_of_functional_groups(db_name,table_to_filter,filter_name,output_table_name,max_count,write_excluded):
    """
    This function will exclude molecules from a table based on the number of times a functional group indicated byt the provided 'filter_name' is present in the molecule.

    The number of functional group ocurrences can be explicitly indicated.

    The compounds passing the filter will be written to 'output_table_name', while if the options to store the excluded compounds (e.g. for inspection purposes) is selected, a table '{output_table_name_excluded}' will be created.

    The output table will contain a column named '{filter_name}_count' indicating the number of matching occurrences.

    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the source table containing the SMILES notation which is to be filtered based on a SMARTS pattern.
    - table_to_filter: the name of the table containing all the SMILES notation of the molecules to be compared with the SMARTS
    - filter_name: this is the identifier of the type of molecules that is to be searched based on the SMARTS notation.
    - output_table_name: the name of the table to which the compounds matching the limit of functional groups will be written.
    - max_count: the maximum (inclusive) number of functional groups tolerated within the molecule.
    - write_excluded: indicate (1:yes, 0:no) if a table named '{output_table_name_excluded}' is to be created to store excluded molecules.
    
    ------
    Returns:
    ------
    A table containing the molecules containing a number of functional groups below the number indicated.

    """
    conn =  sqlite3.connect(f'{db_name}')
    sql_search_op = f'SELECT SMARTS FROM reactants_types WHERE Reactant == "{filter_name}";'
    query_result = conn.execute(sql_search_op)
    target_smarts = query_result.fetchone()[0]

    sql_search_op = f'SELECT * FROM {table_to_filter};'
    df_to_filter = pd.read_sql_query(sql_search_op,conn)

    pandarallel.initialize(progress_bar=True) # Instatiate the pandarallel package for posterior running.
    df_to_filter[filter_name+'_count'] = df_to_filter['SMILES'].parallel_apply(lambda x: get_number_of_matches(x,target_smarts))

    df_clean = df_to_filter[df_to_filter[filter_name+'_count'] <= max_count]
    df_excluded = df_to_filter[df_to_filter[filter_name+'_count'] > max_count]

    store_df_to_sql(db_name,df_clean,output_table_name,"replace")

    if write_excluded == 1:
        store_df_to_sql(db_name,df_excluded,output_table_name+'_excluded',"replace")


def exclude_reactants_based_on_number_of_functional_groups_dict(db_name,table_to_filter,filter_dict,output_table_name,write_excluded):
    """
    This function will exclude molecules from a table based on the number of times a functional group indicated by the provided 'filter_name' is present in the molecule. In this case, a serie of filters are provided as a dictionary

    The number of functional group ocurrences can be explicitly indicated.

    The compounds passing the filter will be written to 'output_table_name', while if the options to store the excluded compounds (e.g. for inspection purposes) is selected, a table '{output_table_name_excluded}' will be created.

    The output table will contain a column named '{filter_name}_count' indicating the number of matching occurrences.

    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the source table containing the SMILES notation which is to be filtered based on a SMARTS pattern.
    - table_to_filter: the name of the table containing all the SMILES notation of the molecules to be compared with the SMARTS
    - filter_dict: A dictionary containing a key value pair corresponding to the 'filter_name:max_count' values
    - output_table_name: the name of the table to which the compounds matching the limit of functional groups will be written.
    - write_excluded: indicate (1:yes, 0:no) if a table named '{output_table_name_excluded}' is to be created to store excluded molecules.
    
    ------
    Returns:
    ------
    A table containing the molecules containing a number of functional groups below the number indicated.

    """
    
    counter = 1
    
    for item in filter_dict.items():
    
        # Extract the parameters from the dictionary
        filter_name = item[0]
        max_count = item[1]

        # Extract the SMARTS for the corresponding filter name
        conn =  sqlite3.connect(f'{db_name}')
        sql_search_op = f'SELECT SMARTS FROM reactants_types WHERE Reactant == "{filter_name}";'
        query_result = conn.execute(sql_search_op)
        target_smarts = query_result.fetchone()[0]

        # This will search the cpds in the first iteration
        if counter == 1:
            sql_search_op = f'SELECT * FROM {table_to_filter};'
            df_to_filter = pd.read_sql_query(sql_search_op,conn)

        else:
            df_to_filter = df_clean    

        pandarallel.initialize(progress_bar=True) # Instatiate the pandarallel package for posterior running.
        df_to_filter[filter_name+'_count'] = df_to_filter['SMILES'].parallel_apply(lambda x: get_number_of_matches(x,target_smarts))

        df_clean = df_to_filter[df_to_filter[filter_name+'_count'] <= max_count]
        
        df_excluded = df_to_filter[df_to_filter[filter_name+'_count'] > max_count]

        # Add the counter after the first iteration
        counter+=1

    
    # After finishing the last iteration output the final df to a table
    store_df_to_sql(db_name,df_clean,output_table_name,"replace")

    if write_excluded == 1:
        store_df_to_sql(db_name,df_excluded,output_table_name+'_excluded',"replace")

        

def check_smiles_vs_smarts(smiles, target_smarts):
    """
    This function will compare a SMILES input with a given SMARTS pattern, and will return a 0/1 tag in case the match is FALSE/TRUE.

    --- Parameters ---
    - smiles: a SMILES string
    - target_samarts: a SMARTS pattern used to evaluate the structure matching

    --- Returns ---
    A 0/1 tag in case the match is FALSE/TRUE.
    
    """
    try: 
        mol = Chem.MolFromSmiles(smiles)
        smarts_mol_filter = Chem.MolFromSmarts(target_smarts)
                
        if mol.HasSubstructMatch(smarts_mol_filter):
            return 1
        else:
            return 0
    except:
        pass

def get_number_of_matches(smiles, target_smarts):
    """
    This function will search the number of instances a certain substructure provided as 'target_smarts' exists in a molecule given its SMILES string.

    Te function will return the number of instances at which the group is present.

    ------
    Parameters:
    ------
    - smiles: the SMILES string corresponding to the molecule to be searched.
    - target_smarts: the functional group given as the SMARTS notation that is going to be searched.

    ------
    Returns:
    ------
    An integer indicating the number of times the corresponding functional group is present.
    """

    try: 
        mol = Chem.MolFromSmiles(smiles)
        smarts_mol_filter = Chem.MolFromSmarts(target_smarts)
        group_numbers = len(mol.GetSubstructMatches(smarts_mol_filter))
        
        return group_numbers
    
    except:
        pass

def read_database_record(db_name,table_name,field_header_to_retrieve,field_header_as_key,key):
    """
    This function will connect to a database, and retrieve the cell value as provided by paramenters
    ------
    Parameters
    ------
    - db_name: the full path to the database to be searched.
    - table_name: the name of the table containing the record to be retrieved.
    - field_header_to_retrieve: the field name were the record to retrieve exists.
    - field_header_as_key: the field name containing the key used to identify the record to extract.
    - key: the key value identifying the record to be extracted.
    
    """

    conn = sqlite3.connect(db_name)
    sql = f"""SELECT {field_header_to_retrieve}
              FROM {table_name}
              WHERE {field_header_as_key} = '{key}';"""

    query = conn.execute(sql).fetchone()[0]

    return query

def read_database_field(db_name,table_name,field_name):
    """
    This function will return as a pandas dataframe all the values of the field stored in the database as indicated by the search paramenters.
    ------
    Parameters
    ------
    - db_name: the name of the database to be searched.
    - table_name: the name of the table containing the field to be extracted.
    - field_name: the name of the field whose content will be returned as a pandas dataframe.
    """

    conn =sqlite3.connect(db_name)
    sql = f"""SELECT {field_name}
              FROM {table_name};"""

    query_df = pd.read_sql_query(sql,conn)
    
    return query_df

def read_database_table(db_name,table_name):
    """
    This function will read the whole set of fields from a table and will return it as a dataframe with column names corresponding to the original values.
    ------
    Parameters
    ------
    - db_name: the full path and name of the database to be searched
    - table_name : the name of the field whose content will be returned as a pandas dataframe.

    ------
    Returns:
    ------
    - query_df: a dataframe containing of the fields of the selected database table, included as columns in the dataframe.
    """

    conn =sqlite3.connect(db_name)
    sql = f"""SELECT *
              FROM {table_name};"""

    query_df = pd.read_sql_query(sql,conn)
    
    return query_df

def insert_smiles_in_table_from_file(db_name,table_name,smiles_file):
    """
    This function will insert one or more SMILES strings provided in an input file within a table in the database.
    
    The file provided may contain any number of columns, but MUST contain a column named 'SMILES'.

    In addition, the table stored will contain the inchi_key values.

    ------
    Parameters:
    ------
    - dbname: the full path to the database in which the record(s) is(are) to be inserted.
    - table_name: the table in which the SMILES strings and correspongins 'inchi_key' values are to be stored. 
    - smiles_file: the file containing a header named "SMILES" followed by one or more lines containing the corresponding smiles string
    
    ------
    Returns
    ------
    A table containing the SMILES strings and inchi_key value is stored within the indicated database.
    """
    
    df = pd.read_csv(smiles_file)
    
    # Instatiate the pandarallel package for posterior running
    pandarallel.initialize(progress_bar=True) 

    df["SMILES_clean"] = df["SMILES"].parallel_apply(lambda x: smiles_proc.generate_clean_smiles_2(x))
    df["inchi_key"] = df["SMILES_clean"].parallel_apply(lambda x: smiles_proc.compute_inchi_key(x))

    df_store = df[["SMILES_clean","inchi_key"]]
    df_store = df_store.rename(columns={"SMILES_clean":"SMILES"})
    db_ops.store_df_to_sql(db_name,df_store,table_name,"replace")

def add_smiles_in_table_from_file(db_name,table_name,smiles_file):
    """
    This function will add one or more SMILES strings provided in an input file within a preexisting table in the database.
    
    The file provided may contain any number of columns, but MUST contain a column named 'SMILES'.

    In addition, the table stored will contain the inchi_key values.

    ------
    Parameters:
    ------
    - dbname: the full path to the database in which the record(s) is(are) to be inserted.
    - table_name: the table in which the SMILES strings and correspongins 'inchi_key' values are to be stored. 
    - smiles_file: the file containing a header named "SMILES" followed by one or more lines containing the corresponding smiles string
    
    ------
    Returns
    ------
    A table containing the SMILES strings and inchi_key value is stored within the indicated database.
    """
    
    df = pd.read_csv(smiles_file)
    
    # Instatiate the pandarallel package for posterior running
    pandarallel.initialize(progress_bar=True) 

    df["SMILES_clean"] = df["SMILES"].parallel_apply(lambda x: smiles_proc.generate_clean_smiles_2(x))
    df["inchi_key"] = df["SMILES_clean"].parallel_apply(lambda x: smiles_proc.compute_inchi_key(x))

    df_store = df[["SMILES_clean","inchi_key"]]
    df_store = df_store.rename(columns={"SMILES_clean":"SMILES"})
    db_ops.store_df_to_sql(db_name,df_store,table_name,"append")
    
def subset_molecules_based_on_index(db_name,table_origin,table_dest):
    """
    This function will extract and write to a subset table containing a set of molecules in SMILES format that are extracted from another table by indication the corresponding index as obtained from the depiction grids. The function will systematically ask for indexes names and construct the corresponding query list until the stopping keyword is entered.
    ------
    Params
    ------
    - db_name: the full path to the database were the retrieving and writing of the corresponding compounds.
    - table_origin: the name of the table were the compounds are to be extracted from.
    - table_dest: the name of the destination table were the compounds are to be written.
    """
    index_of_cps_list = []
    while True:
        try: 
            input_index = int(input("Enter the row id to extract (Enter -1 to end addition): "))
            if input_index == -1:
                print("Finishing list contruction")
                print(f"List constructed includes rows: \n {index_of_cps_list}")
                break
            elif input_index not in index_of_cps_list:
                index_of_cps_list.append(input_index)
            elif input_index in index_of_cps_list:
                print(colored("The cpd row id already exist, enter a new value","red"))    
        except:    
            print(colored("Enter an integer index","red"))
            
    conn = sqlite3.connect(db_name)
    
    for item in index_of_cps_list:
        extracted_df = db_actions.extract_by_row_id(db_name,table_origin,item)
        store_df_to_sql(db_name,extracted_df,table_dest,"append")
        print(colored(f"Compound written to: {table_dest} in database"))

def subset_available_molecules(db_name,table_name):
    """
    This function will extract all the molecules within the given table using as pattern a column 'available' == 1.
    
    The function will search for a field/column named 'available', which should exist to terminate the execution.
    
    All the molecules matching that criteria will be stored within the same database under the name 'table_name'_avail. 
    
    If the destination table already exists, it will be replaced.
    
    ------
    Parameters
    ------
    - db_name: the full path to the database containing the table to be searched for available molecules.
    - table_name: the name of the table containing a field named 'available' to be searched.
    
    ------
    Returns
    ------
    A table named 'table_name'_avail will be written to the same database containing molecules tagged with 'available' == 1.
    """

    conn = sqlite3.connect(db_name)
    
    sql =  f"""SELECT *
               FROM {table_name}
               WHERE available == 1;"""
               
    try:
        df_tagged = pd.read_sql_query(sql,conn)
        store_df_to_sql(db_name,df_tagged,table_name+'_avail',"replace")
        print("Extraction of available compounds SUCCESSFUL")
        
    except:
        print("No tagged molecules were found. Check the action.")

def subset_molecules_by_property_range(db_name,table_origin,table_dest,parameter,lower_limit,upper_limit):
    """
    Given a table in which molecular properties are present for the corresponding molecules, this function will subset molecules based on a provided lower:upper range.

    ------
    Parameters:
    ------
    - db_name: the database name in which the table containing the precalculated molecular descriptors exists.
    - table_origin: the original table to be processed.
    - table_dest: the table within the same dabatase to which the subseted molecules will be written.
    - parameter: the target descriptor to be used for subseting by properties.
    - lower_limit: the lower (inclusive) value of the parameter to be used.
    - upper_limit: the higher (inclusive) value of the parameter to be used.
    """

    conn = sqlite3.connect(db_name)
    sql = f"""SELECT *
              FROM {table_origin}
              WHERE {parameter} >= {lower_limit}
              AND {parameter} <= {upper_limit};"""

    try:

        df = pd.read_sql_query(sql,conn)
        store_df_to_sql(db_name,df,table_dest,"replace")
        print(f"Table {table_origin} SUCCESFULLY subseted to {table_dest}")

    except Exception as error:
        print(error)
        print("Problem subseting the provided table")

def subset_molecules_by_set_of_properties_range(db_name,table_origin,table_dest,values_list):
    """
    This function will filter a table containing molecules with calculated properties, using multiple conditions as indicated in a properties list of lists. An example dictionary is:

        values_list = [["A",[">",30]],
               ["A",["<",50]],
               ["B",["<",20]],
               ["C",[">",70]],
              ]
 
    In order to efficiently define the values_list it is advisable to be assisted by the histograms plotted for every property

    ------
    Parameters
    ------
    db_name: the full path to the database containing the folder to be filtered.
    table_origin: the name of the table to be filtered.
    table_dest: tne name of the table in which the filtered compounds are to be stored.
    values_list: the values_list containing the multiple filtering criteria. Pay special attention to the sample dictionary format.

    ------
    Returns
    ------
    table_dest is written to the same database. If exists, it will be replaced by the new sunseting.
    """

    conn = sqlite3.connect(db_name)
    sql = f"SELECT * FROM {table_origin}"

    try:
        df = pd.read_sql_query(sql,conn)
        #query = ' and '.join([f'({k} {v[0]} {v[1]})' for k, v in values_dict.items()])
        query = ' and '.join([f'{v[0]} {v[1][0]} {v[1][1]}' for v in values_list])
        subset_df = df.query(query)
        store_df_to_sql(db_name,subset_df,table_dest,"replace")
        print(f"Table {table_origin} SUCCESFULLY subseted to {table_dest}")

    except Exception as error:
        print(error)
        print("Problem subseting the provided table")

def subset_molecules_by_column_value(db_name,table_origin,table_dest,column_name,value,action):
    """
    Given a table in which molecular properties are present for the corresponding molecules, this function will subset molecules based on a given column and provided value.

    ------
    Parameters:
    ------
    - db_name: the database name in which the table containing the precalculated molecular descriptors exists.
    - table_origin: the original table to be processed.
    - table_dest: the table within the same dabatase to which the subseted molecules will be written.
    - column_name: the target descriptor to be used for subseting by properties.
    - value: the string or integer value to use for filtering the corresponding compounds.
    - action: if the dest_table exists, indicate the action to apply: "replace", "append"
    """

    conn = sqlite3.connect(db_name)
    sql = f"""SELECT *
              FROM {table_origin}
              WHERE {column_name} = '{value}';"""

    try:
        df = pd.read_sql_query(sql,conn)
        store_df_to_sql(db_name,df,table_dest,action)
        print(f"Table {table_origin} SUCCESFULLY subseted to {table_dest}")

    except Exception as error:
        print(error)
        print("Problem subseting the provided table")

def purge_duplicates_in_table(db_name,table_name,field_key):
    """
    This function will evaluate the content of a table within a database, and provided a target column, will evaluate whether there are duplicates in it, an in that case will delete the duplicated record.
    ------
    Parameters:
    ------
    - db_name: the full path the database file where the duplicates checkings will be evaluated.
    - table_name: the table in which the search for duplicates will occur.
    - field_key: This is the name of the field used to search duplicates in order to delete duplicated rows.

    ------
    Returns
    ------
    The processed table without duplicates will replace the original table in the database.
    """

    conn = sqlite3.connect(db_name)
    sql = f"""SELECT *
              FROM '{table_name}';"""
    
    df = pd.read_sql_query(sql,conn)
    df.drop_duplicates(subset=[field_key],inplace=True)

    db_ops.store_df_to_sql(db_name,df,table_name,"replace")
        
    print("Duplicates dropped SUCCESFULLY")

def mark_molecule_as_available(db_name,table_name,cpd_availability_file):
    """
    This function will add a tag = 1 in the column names 'available' of a table containing molecules.

    The availability is read from a .csv file containing a SMILES string.

    Prior to tagging compounds, it is advisable to check the presence of the compounds to be tagged in the target table using the corresponding function.

    ------
    Parameters
    ------
    - db_name: the full path to the database in which the tagging is supossed to happen
    - table_name: the name of the table in which the compounds to be tagged exists.
    - cpd_availability_file: a .csv file with a header SMILES including the compounds to be tagged.

    ------
    Returns:
    ------

    The 1 tag is added to the compound in the target table
    """

    conn = sqlite3.connect(db_name)

    sql = f"""SELECT *
             FROM {table_name};"""
    
    df = pd.read_sql_query(sql,conn)

    df_tags = pd.read_csv(cpd_availability_file)
       
    pandarallel.initialize() # Instatiate the pandarallel package for posterior running.
    df_tags["SMILES_clean"] = df_tags["SMILES"].parallel_apply(lambda x: general_proc.generate_clean_smiles_2(x))
    df_tags["inchi_key"] = df_tags["SMILES"].parallel_apply(lambda x: general_proc.compute_inchi_key(x))
    df_tags["available"] = df_tags["SMILES"].parallel_apply(lambda x: 1)
    df_tags.drop(columns=["SMILES"],inplace=True)
    df_tags.rename(columns={"SMILES_clean":"SMILES"},inplace=True)
    df_append = df_tags[["SMILES","inchi_key","available"]]
    df_merged = pd.concat([df,df_append],axis=0)
    df_merged.drop_duplicates(subset='inchi_key',keep="last",inplace=True)

    store_df_to_sql(db_name,df_merged,table_name,"replace")

###############
if __name__ == "__main__":
    print("This is a test")
    
    #select_reactants_based_on_smarts("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reactants_types","alpha_aminoacid","reagents_pass")
   
    #select_reactants_based_on_smarts("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","reactants_types","primary_aliphatic_amine","reagents_pass")
    
    #select_reactants_based_on_smarts("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reactants_types","alpha_aminoacid","reagents_pass")

    #select_reactants_based_on_smarts_tqdm("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reactants_types","alpha_aminoacid","reagents_pass")

    #select_reactants_based_on_smarts_tqdm("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reactants_types","primary_aliphatic_amine","reagents_pass")

    #select_reactants_based_on_smarts_pandarallel("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","reactants_types","primary_aliphatic_amine","reagents_pass")

    #select_reactants_based_on_smarts_pandarallel("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reagents_pass","reactants_types","alpha_aminoacid")
    
    #exclude_reactants_based_on_number_of_functional_groups("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/Emolecules_reagents_4.db","alpha_aminoacid","carboxylic_acid","filter_test",1,1)

    #delete_record("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","reaction_types",'Reaction','primary_alyphatic_amines')
    
    #read_database_record("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reaction_types","SMARTS","Reaction","azide_from_amine")
    
    #read_database_field("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alpha_aminoacid","SMILES")

    #df = read_database_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alpha_aminoacid")

    #insert_smiles_in_table_from_file("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/alcohols.csv")

    #insert_smiles_in_table_from_file("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/alcohols.csv")

    #add_smiles_in_table_from_file("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/alcohols.csv")

    #insert_smiles_in_table_from_file("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/alcohols.csv")

    #subset_molecules_based_on_index("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","alpha_aminoacid","alpha_aminoacid_subset")
    
    #subset_molecules_based_on_index("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alpha_aminoacid","alpha_aminoacid_subset")
    
    #purge_duplicates_in_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alpha_aminoacid_subset",1)

    #purge_duplicates_in_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","alcohols","inchi_key")

    #mark_molecule_as_available("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alpha_aminoacid","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/files_for_testing/cpd_availability.csv")
    
    #subset_available_molecules("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alpha_aminoacid")

    #subset_molecules_by_column_value("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/Emolecules_reagents_4.db","alpha_aminoacid","aa_subset","NumHAcceptors",4)