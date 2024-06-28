import sqlite3
import pandas as pd
import time
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.general_procedures.operations_on_sql_database as sql_ops

def check_lists_size(list1,list2,lines_limit,table_name,num_lines,counter,actual_time,db_name):
    actual_lines = len(list1)
    if actual_lines > lines_limit:
        time_elapsed = time.time()
        time_taken = round((time_elapsed-actual_time),2)
        print(f"Hice un backup: {counter}/{num_lines} - {time_taken} segs")
        df =merge_2_lists_into_df(list1,list2,"SMILES","Inchi_Key")
        sql_ops.store_df_to_sql(db_name,df,table_name,"append")
        list1 = []
        list2 = []
        return list1, list2,time_elapsed
    else:
        return list1, list2, actual_time

def check_1_list_size(list,lines_limit,table_name,db_name,action):
    actual_lines = len(list)
    if actual_lines > lines_limit:
        df =convert_1_list_into_df(list,"Inchi_Key")
        sql_ops.store_df_to_sql(db_name,df,table_name,action)
        list = []
        return list
    else:
        return list

def merge_2_lists_into_df(list1,list2,colname1,colname2):
    # Join the data into a dataframe for further storing into database
    df = pd.DataFrame(zip(list1,list2),columns=[colname1,colname2])
    
    return df

def merge_3_lists_into_df(list1,list2,list3,colname1,colname2,colname3):
    # Join the data into a dataframe for further storing into database
    df = pd.DataFrame(zip(list1,list2,list3),columns=[colname1,colname2,colname3])
    
    return df

def merge_4_lists_into_df(list1,list2,list3,list4,colname1,colname2,colname3,colname4):
    # Join the data into a dataframe for further storing into database
    df = pd.DataFrame(zip(list1,list2,list3,list4),columns=[colname1,colname2,colname3,colname4])
    
    return df

def convert_1_list_into_df(list,colname):
    # Join the data into a dataframe for further storing into database
    df = pd.DataFrame(zip(list),columns=[colname])
    
    return df

