import sqlite3
import pandas as pd
from termcolor import colored
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.input_output_operations.i_o_checkings as i_o_check
#import general_procedures.operations_on_sql_database as sql_ops
import chemspace.general_procedures.smiles_processing as smiles_proc
import chemspace.general_procedures.operations_on_sql_database as ops


def drop_table_from_db(db_name,table_name):
    """
    This function will drop an entire table from a database
    """
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    try:
        cursor.execute(f'DROP TABLE {table_name};')
        conn.commit
        cursor.close()
        print(f"The table {table_name} was dropped succesfully")
    except Exception as error:
        print(colored(f"Drop of the table {table_name} failed","red"))
        print(error)
        exit()
    
def drop_record_from_db(db_name,table_name,field_name,value):
    """
    This function will drop an entire table from a database
    """
    conn = sqlite3.connect(db_name)
    try:
        sql = f'''DELETE FROM {table_name} 
                  WHERE {field_name} = "{value}"'''
        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()
        print("Record dropped successfully")
    except:
        print("Failed to drop record")

def modify_cell_in_db(db_name,table_name,field_name,old_value,new_value):
    """
    This function will modify a cell in a table
    """
    conn = sqlite3.connect(db_name)
    try:
        
        sql = f'''UPDATE {table_name} 
                  SET {field_name} = "{new_value}" 
                  WHERE {field_name} = "{old_value}"'''
        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()
        print("Record modified successfully")
        cur.close()
    except:
        print("Failed to modify record")

def rename_table_in_db(db_name,origin,destination):
    """
    This function will rename with 'destination' a cell in a table that originally contained "origin"
    """
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    # This uses a two step renaming through a 'temp_table' because SQL is case insensitive, and there are cases where just the case is wanted to be modified.
    sql_1 = f'''ALTER TABLE {origin} 
               RENAME TO temp_table;'''
    cursor.execute(sql_1)

    sql_2 = f'''ALTER TABLE temp_table 
                RENAME TO {destination};'''
    cursor.execute(sql_2)
    conn.commit
    cursor.close()
    
def rename_colname_in_db(db_name,table,origin,destination):
    """
    This function will rename with 'destination' a cell in a table that originally contained "origin"
    """
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    # This uses a two step renaming through a 'temp_column' because SQL is case insensitive, and there are cases where just the case is wanted to be modified.
    try: 
        sql_1 = f'''ALTER TABLE {table} 
                    RENAME COLUMN {origin} 
                    TO temp_column;'''
        cursor.execute(sql_1)

        sql_2 = f'''ALTER TABLE {table} 
                    RENAME COLUMN temp_column 
                    TO {destination};'''
        cursor.execute(sql_2)
        conn.commit
        cursor.close()
        print("Column updated succesfully")
    except:
        print("Error in updating column")

def extract_joined_data_by_key(db_name,table_1,table_2,key,field_to_extract):
    """
    This function will return data based on a join of key of two tables
    ------
    Parameters
    ------
    - db_name: is the database were the information to be extracted is stored
    - table_1: is the table were the field_to_extract exists
    - table_2: is the table were the keys for joining exists
    - key: this is the key present in both tables in order to perform the inner join.
    - field_to_extract: this is the information that is to be extracted based on the inner join procedure
    """
    conn = sqlite3.connect(db_name)
    sql = f'''SELECT DISTINCT {table_1}.{field_to_extract}, {table_1}.{key} 
              FROM {table_1} 
              INNER JOIN {table_2} ON {table_2}.{key} = {table_1}.{key};'''
    query_df = pd.read_sql_query(sql,conn)
    #print(sql)
    #print(len(query_df))
    
    smiles_list = []
    inchi_key_list = []
    
    for index, row in query_df.iterrows():        
        if row.iloc[1] not in inchi_key_list:
    
            smiles_list.append(row.iloc[0])
            inchi_key_list.append(row.iloc[1])
    
    df = i_o_check.merge_2_lists_into_df(smiles_list,inchi_key_list,"SMILES","inchi_key")
    
    return df

def extract_row_by_value(db_name,table_name,target_column,value):
    """
    This function will extract a record matching a specific value, and will return it as a dataframe for further storing if required.
    ------
    Parameters
    ------
    - db_name: the full path of the database were the information to be extracted is stored.
    - table_name: the table were the record to extract exists.
    - target_column: The name of the column were the provided 'value' matching the record is stored.
    - value: the string corresponding to the value for the record to be extracted.
    """
    conn = sqlite3.connect(db_name)
    sql = f'''SELECT * 
            FROM {table_name} 
            WHERE {target_column} = "{value}";'''
    
    query_df = pd.read_sql_query(sql,conn)

    return query_df

def extract_by_row_id(db_name,table_name,row_id):
    """
    This function will extract a record matching a specific index, and will return it as a dataframe for further storing if required.
    ------
    Parameters
    ------
    - db_name: the full path of the database were the information to be extracted is stored.
    - table_name: the table were the record to extract exists.
    - index: the value corresponding to the index to be extracted.
    """
    conn = sqlite3.connect(db_name)
    sql = f'''SELECT * 
            FROM {table_name} 
            WHERE rowid = {row_id};'''
    
    query_df = pd.read_sql_query(sql,conn)

    return query_df

def check_reactants_presence(db_name, target_inchi_table,csv_file):
    """
    This function will check the presence of reactants in a database based on a .csv file provided containing a 'SMILES' column. The corresponding inchi_key are computer and evaluated for presence.
    ------
    Parameters
    ------
    - db_name: the full path to the database to be evaluated.
    - target_inchi_table: the table containing 'SMILES:inchi_key' pair to be evaluated. It should be noted that a column named 'inchi_key' must be available in this target table.
    - csv_file: the full path to the .csv containing a 'SMILES:name' pair to be evaluated for presence in the database. The separator should be a ',' and the header is 'SMILES,name'.
    """

    conn = sqlite3.connect(db_name)
    sql = f""" SELECT inchi_key
              FROM {target_inchi_table};
            """           
    query_result = list(conn.execute(sql).fetchall())
    query_result_list = [i[0] for i in query_result] # will generate a list with all the resulting inchi_key
    
    print(len(query_result_list))
    to_search_df = pd.read_csv(csv_file)
    
    for index, row in to_search_df.iterrows():
        try: 
            target_inchi_key = smiles_proc.compute_inchi_key(row["SMILES"])
            print(target_inchi_key)
            if target_inchi_key in query_result_list:
                print(f'The compound ' + colored(f'{row["SMILES"]}','green') + ' is present')
            else:
                print(f'The compound ' + colored(f'{row["SMILES"]}','red') + ' is NOT present')
        except:
            print(f"Error generation inchi_key for {row['SMILES']}")

def create_2_fields_table(db_name,table_name,field_1_name,field_2_name,field_1_type,field_2_type):
    """
    This function will create a table containing two field with the given fields name and field atributes.
    ------
    Parameters:
    ------
    - db_name: the full path to the database in which the table containing two fields is to be created.
    - table_name: the name of the correspoding table to be created
    - field_1_name: the name of col_1
    - field_2_name: the name of col_2
    - field_1_type: the data type of col_1
    - field_2_type: the data type of col_2
    """

    conn = sqlite3.connect(db_name)
    sql = f"""CREATE TABLE {table_name}
              ({field_1_name} {field_1_type},
               {field_2_name} {field_2_type});"""

    
    try:
        conn.execute(sql)
        print("Table SUCCESFULLY created")
    except:
        print("ERROR creating the table")

def copy_table_between_dbs(origin_db_name,destination_db_name,table_name):
    """
    This function will copy an entire table from Database_1 (origin) to Database_2 (destination), will all the corresponding contents. 
    
    It is important to note that if the table to be copied from DB_origin to DB_destination already exists, the table in the destination will be replaced.

    --- Parameters ---
    - origin_db_name: the full path of the database FROM which the table is to be copied.
    - destination_db_name: the full path of the database TO which the table is to be copied.
    - table_name: the name of the table that is to be copied.

    --- Returns ---
    A copy of the requested table into the destination database.
    """

    conn = sqlite3.connect(origin_db_name)
    sql = f"""SELECT *
             FROM {table_name};"""

    try: 
        df = pd.read_sql_query(sql,conn)
        ops.store_df_to_sql(destination_db_name,df,table_name,"replace")
        
        print(f"Table '{table_name}' copied SUCESSFULY")

    except:
        print(f"ERROR copying the table: '{table_name}'")

def merge_two_tables(db_name,table_1,table_2,new_table_name):
    """
    This function will receive two tables, merge them and save the resulting new table to the same database under the provided 'new_table_name'.

    There is no need that both tables contain the same number of columns or equal headers.

    In case column number of name are different, de 'None' value will be inserted.
    
    ------
    Parameters
    ------
    - db_name: the full path to the database file in which the two tables to merge exists.
    - table_1: name of table_1 destinated to the merging.
    - table_2: name of table_2 destinated to the merging.
    -new_table_name: the name of the merged table that will be saved to the database file.

    ------
    Returns
    ------
    A new table containing the merged information is saved in the same database.
    """

    conn = sqlite3.connect(db_name)

    sql_1 = f"""SELECT *
             FROM {table_1};"""
    
    sql_2 = f"""SELECT *
             FROM {table_2};"""
    
    try: 
        df_1 = pd.read_sql_query(sql_1,conn)
        df_2 = pd.read_sql_query(sql_2,conn)
        df_merged = pd.concat([df_1,df_2],axis=0)
        ops.store_df_to_sql(db_name,df_merged,new_table_name,"replace")
        print("Successfully merged the tables.")

    except:
        print("Problem reading the tables. Check the parameters.")
    
def extract_table_as_csv(db_name,table_name,dest_file):
    """
    This function will export a table and it whole contents as a .csv file.
    
    The function is intended for sharing tables as single files between users of the package.
    
    This function is combined with the 'create_table_from_csv()' function to restore the table from a .csv source file.
    
    ------
    Parameters
    ------
    - db_name: the full path to the database file from which the table to be exported is located.
    - table_name: the name of the table to be exported as a .csv file.
    - dest_file: the full path to the file were the .csv file is to be written.
    
    ------
    Returns
    ------
    A .csv file with al the table contents is created in the specified path.
    """

    conn = sqlite3.connect(db_name)
    
    sql = f"""SELECT * 
             FROM {table_name};"""
             
    try: 
        df = pd.read_sql_query(sql,conn)
        df.to_csv(f'{dest_file}',sep=',',header=True,index=None)
        print(f"Saving of table: {table_name} SUCCESSFUL")
    
    except Exception as error:
        print(error)
        print(f'Extraction of table: {table_name} FAILED.')

def create_table_from_csv(db_name,origin_file,dest_table_name):
    """
    This function will create a table within the specified database taking as input a .csv file.
    
    For compatibility purposes, it is recommended that the .csv file has been crated by the 'extract_table_as_csv()' function of this toolkit.
    
    If the table name already exists in the database, the creation will FAIL to avoid unwanted overwriting of information.
    
    ------
    Parameters
    ------
    - db_name: the full path to the database file in which the .csv file will be inserted a new table.
    - origin_file: the full path to the .csv file containing the info to create the table in the database.
    - dest_table_name: the name of the table to be created with the corresponding information.
    """

    try: 
        df = pd.read_csv(origin_file)
        conn = sqlite3.connect(db_name)
        ops.store_df_to_sql(db_name,df,dest_table_name,"fail")
        print("Table created SUCCESSFULLY")
    except:
        print("Table NOT created. Check for errors.")

def replicate_table_empty_state(table,src_conn,dest_conn):
    """
    This function will prepare an empty table by emulating the structure (columns and columns types) of a given table.

    This is useful to copy a table between different databases in a SQL row wise manner.

    ------
    Parameters:
    ------
    - table: the name of the table to be copied form the source connection to the destination connection.
    - src_conn: an SQL connection to the source database as provided by the main function.
    - dest_conn: an SQL connection to the destination database as provided by the main function.
    
    ------
    Returns
    ------
    An empty table with the required structure for a row wise copy is returned in the dest_conn
    """

    sc = src_conn.execute(f'SELECT * FROM {table}')
    ins = None
    for row in sc.fetchall():
        if not ins:
            cols = tuple([k for k in row.keys() if k != 'id'])
            ins = "EXIT"
    
    ## The PRAGMA instruction will get all the column types
    ins_2 = f'PRAGMA table_info({table})'
    sc.execute(ins_2)
    typeInfos = sc.fetchall()
    columnTypes = [item[2] for item in typeInfos] 
    
    ins = None
    
    for index, item in enumerate(cols):
        type = columnTypes[index]

        if index == 0:
            ins = f'CREATE TABLE {table} ("{item}" {type})'
            dest_conn.execute(ins)
            print("Destination table created")
        else:
            ins = f'ALTER TABLE {table} ADD {item} {type}'
            dest_conn.execute(ins)
            continue

def establish_row_factory_connection(db_name):
    """
    This function will establish an return a connection object to the given database, and return what is named a 'dictionary cursor', which returns dictionaries instead of tuples of a normal connection

    This type of connection is used by the function that copies tables in a SQL row based way.

    ------
    Parameters:
    ------
    - db_name: the full path to the database were the row_factory connection is to be established.

    ------
    Returns:
    ------
    A 'row factory' connection type.
    """
    conn = sqlite3.connect(db_name)
    # Let rows returned be of dict/tuple type
    conn.row_factory = sqlite3.Row
    
    return conn

def copy_table_between_dbs(src_db,dest_db,table):
    """
    This function will copy one table between two different databases: src_db, dest_db.

    The name of the table to be copied will be preserved, and the structure and contents of the table will be replicated.

    Is the table already exists in the 'dest_db' it will be replaced.

    ------
    Parameters:
    ------
    - src_db: the source database from which the table will be retrieved.
    - dest_db: the destination database in which the table will be copied.
    - table: the name of the table in 'src_db' that will be copied to 'dest_db'

    ------
    Returns:
    ------
    An identical copy of the table into the 'dest_db'
    
    """
    # Create the 'factory type' connection
    src_conn = establish_row_factory_connection(src_db)
    dest_conn = establish_row_factory_connection(dest_db)

    # Replicate an empty table in dest_db for copying in a row-based fashion
    replicate_table_empty_state(table,src_conn,dest_conn)

    # Copy iteratively the source rows into destination

    sc = src_conn.execute(f'SELECT * FROM {table}')
    ins = None
    dc = dest_conn.cursor()

    for row in sc.fetchall():
        if not ins:
            cols = tuple([k for k in row.keys() if k != 'id'])
            ins = 'INSERT OR REPLACE INTO %s %s VALUES (%s)' % (table, cols,','.join(['?'] * len(cols)))
        c = [row[c] for c in cols]
        dc.execute(ins, c)

    dest_conn.commit()

    print(colored(f"Destination table {table} SUCESSFULLY copied.","green"))

### Section to test functions 
    
if __name__ == "__main__":
    print("Write test function")
    
    #drop_table_from_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","aa_w/o_diacids")

    #rename_table_in_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","alchohols","alcohols")
    
    #rename_colname_in_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reagents_pass","Inchi_key","inchi_key")
    
    #drop_record_from_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","reactants_types","Reactant","di_carboxylic_acid")
    
    #modify_cell_in_db("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","reactants_types","Reactant","para_di_substituted_benzotriazole","para_di_subs_benzotriazole")
    
    #extract_joined_data_by_key("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reagents_pass","aminoacids","inchi_key","SMILES")
    
    #extract_row_by_value("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reaction_types","Reaction", "azide_from_amine")
    
    #df = extract_by_row_id("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","azides",1)
    
    #check_reactants_presence("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/Emolecules_reagents.db","alpha_aminoacid","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/all_natural_aminoacids.csv")

    #create_2_fields_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alcohols","SMILES","inchi_key","TEXT","TEXT")

    #copy_table_between_dbs("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents.db","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_2.db","reaction_types")

    #merge_two_tables("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","alpha_aminoacid","alcohols","aa_alcohols")
    
    #extract_table_as_csv("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/Emolecules_reagents_2.db","reaction_types","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/reactions_1.csv")
    
    #create_table_from_csv("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/alpha_aminoacid_copy.csv","table_created")

    copy_table_between_dbs("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/sample.db","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/Emolecules_reagents_4.db","alpha_aminoacid_pdbqt_ligands")