import sqlite3
from termcolor import colored
import pandas as pd
import json
from pprint import pprint
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

def create_params_registry_table(docking_params_registry_path):
    """
    This function will create a database named 'docking_parameters.db', which contains a table named 'docking_params'.
    
    This database is intended to store all the parameters required by AutoDock-GPU using a dictionary based format.
    
    The 'docking_params' contains the following structure:
    
        - condition_id: an integer describing the set of conditions to be used
        - condition_dict: a dictionary containing all the available docking parameters
        - condition_description: a brief description of the parameters set
    
    ------
    Parameters:
    ------
    - docking_params_registry_path: the full path to the folder in which the 'docking_parameters.db' database is stored.
    
    ------
    Returns:
    ------
    Creates a database and table to store the parameters to be applied to docking assays.
    """
    
    try: 
        conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')

        ins = f'CREATE TABLE docking_params ("condition_id" INTEGER, "condition_dict" TEXT, condition_description TEXT)'
        conn.execute(ins)
        print(colored("Docking conditions registers table SUCCESSFULLY created.","green"))
        
    except Exception as error:
        print(error)
        print(colored("Error creating the database 'docking_parameters.db'","red"))

def drop_params_registry_table(docking_params_registry_path):
    """
    This function will delete the table named: 'docking_params' stored in the database 'docking_parameters.db' under the 'docking_params_registry_path' path.
    """
    try: 
        conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')
        sql = "DROP TABLE docking_params"
        conn.execute(sql)
        conn.close()
        print(colored("SUCCESSFULLY delete the table 'docking_params'","green"))
    except Exception as error:
        print(error)
        print(colored("Error deleting the table 'docking_params'","red"))

def clear_params_registry_table(docking_params_registry_path):
    """
    This function will erase all the molecular docking paramters stored within the table named 'docking_params' in the 'docking_parameters.db' under the 'docking_params_registry_path' path.

    ------
    Parameters:
    ------
    - docking_params_registry_path: the full path to the folder in which the 'docking_parameters.db' database is stored.

    ------
    Returns:
    ------
    None
    """
    try: 
        conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')
        cursor = conn.cursor()
        sql = "DELETE FROM docking_params" # This will delete all the registries in the table.
        query = input(colored("Are you sure you want to delete all docking params registries? (y/n): ","magenta"))
        if query == 'y':
            cursor.execute(sql)
            conn.commit()
            print(colored("All registries corresponding to docking params were DELETED.","green"))
        else:
            exit()
    except Exception as error:
        print(error)
        print(colored("Deletion of docking params registries FAILED","red"))

def create_default_dict():
    """
    This function will create and return dictionary containing all available docking paramentes for AutoDock-GPU.
    
    The corresponding dictionary will be returns for use at set conditions as default or for further modifications.
    
    ------
    Parameters:
    ------
    None
    
    ------
    Returns:
    ------
    A dictionary containing all available parameters with default values.
    """      

    general_opts_dict = {'files_ops':{},
                 'conversion_ops':{},
                 'output_ops':{},
                 'setup_ops':{},
                 'search_ops':{},
                 'scoring_ops':{},
                 }

    files_ops_dict = {'--lfile':'LFILE',
                       '--ffile':'FLDFILE',
                       '--flexres':'FLEXFILE',
                       '--filelist':'BATCHFILE',
                       '--import_dpf':"DPFFILE",
                       '--xraylfile': "XRAYLIGFILE",
                       }

    conversion_ops_dict = {'--xml2dlg':'None',
                        }
   
    output_ops_dict = {'--resname':'ligand basename',
                       '--contact_analysis':0,
                       '--xmloutput':1,
                       '--dlgoutput':1,
                       '--dlg2stdout':0,
                       '--rlige':0,
                       '--gfpop':0,
                       '--npdb':0,
                       '--gbest':0,
                       '--clustering':1,
                       '--output-cluster-poses':0,
                       '--hsym':1,
                       '--rmstol':2,
                       }

    setup_ops_dict = {'--devnum':1,
                      '--loadxml':None,
                      '--seed':'time',
                     }
    
    search_ops_dict = {'--heuristics':1,
                       '--heurmax':12000000,
                       '--autostop':1,
                       '--asfreq':5,
                       '--nrun':20,
                       '--nev':2500000,
                       '--ngen':42000,
                       '--lsmet':'ad',
                       '--lsit':300,
                       '--psize':150,
                       '--mrat':2,
                       '--crat':80,
                       '--lsrat':100,
                       '--trat':60,
                       '--dmov':6,
                       '--dang':90,
                       '--rholb':0.01,
                       '--lsmov':2,
                       '--lsang':75,
                       '--cslim':4,
                       '--stopstd':0.15,
                       '--initswgens':0
                       }

    scoring_ops_dict = {'--derivtype':None,
                        '--modpair':None,
                        '--ubmod':0,
                        '--smooth':0.5,
                        '--elecmindist':0.01,
                        '--modqp':0,
                        }
    
    # Configure the general options dictionary with the custom values
    general_opts_dict['files_ops'] = files_ops_dict
    general_opts_dict['conversion_ops'] = conversion_ops_dict
    general_opts_dict['output_ops'] = output_ops_dict
    general_opts_dict['setup_ops'] = setup_ops_dict
    general_opts_dict['search_ops'] = search_ops_dict
    general_opts_dict['scoring_ops'] = scoring_ops_dict

    return general_opts_dict

def check_param_exists(parameter):
    """
    This function will check if a given parameter exists within the set of parameters available in the AutoDock-GPU parameters dictionary.
    
    After checking, a flag '0' will be returned if the parameter does not exist, and '1' if the parameter exists.
    
    If the parameter is found, the parameter category in which it is found will be also returned.
    ------
    Parameters:
    ------
    - parameter: the name of the parameters to check.
    
    ------
    Returns:
    ------
    A flag 0/1 indicating if the parameter does not exist/exists.
    """
    default_dict = create_default_dict()

    for key in default_dict.keys():
        print(key)        
        for sub_key in default_dict[key].keys():
            # Check if the parameter matches
            if parameter == sub_key:
                print(colored(f"Parameter '{parameter}' found in category '{key}'","green"))
                return key, sub_key, 1
            else:
                continue
        
    return 0,0,0 # if the parameter does not exist, all zero values will be returned

def create_custom_params_dict():
    """
    This function will store a new registry in the params_registry table.
    
    A key to be modified is requested, after which the corresponding parameter is requested.
    
    The base params_dict is created using the 'create_default_dict', after which the requester parameters are modified.
    
    ------
    Parameters:
    ------
    - docking_params_registry_path: the full path to the folder in which the parameters database is stored.
    
    ------
    Returns
    ------
    A new registry is written to the registry table only in case the created parameters set is new.
    """ 
    
    # This will create a base params dictionary to be modified
    custom_dict = create_default_dict()
    
    while True:
        parameter = input("Enter the parameter to modify (enter 'DONE' to finish): ")
        
        if parameter == 'DONE':
            print(colored("Finished entering parameters.","green"))
            break
        # Check if the flag exists
        #check_flag = check_param_exists(parameter)
        key,sub_key,check_flag = check_param_exists(parameter)
        
        if check_flag == 0:
            print(colored(f"The parameter '{parameter}' DOES NOT EXIST. Enter the param again.","red"))
            continue
        elif check_flag == 1:
            value = input(f"Provide '{parameter}' parameter value: ")
            custom_dict[key][sub_key]=value
            continue
        
    return custom_dict

def get_last_param_record_id(conn):
    """
    This function will return the last value stored under the column 'condition_id' under the 'docking_parameters.db' database stored within 'docking_params_registry_path'.
    
    ------
    Parameters:
    ------
    - docking_params_registry_path: the full path to the place were 'docking_parameters.db' is located.
    - conn: an sql connection to the 'docking_parameters.db' database.
    ------
    Returns:
    ------
    An integer corresponding to the last value stored under the column 'condition_id'.
    """

    try:
        #conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')
        sql_assay = """SELECT condition_id 
                       FROM docking_params;"""
    
        df = pd.read_sql_query(sql_assay,conn)
        assay_id_list = df["condition_id"].to_list()
        if len(assay_id_list) == 0: # This will hit in case no previous record in the table is present
            return 0
        else:
            last_value = assay_id_list[-1]
            return last_value
    
    except Exception as error:
        print(error)
        print(colored("Error retrieving the 'condition_id' value","red"))
        exit()

def check_params_dict_existence(custom_dict,conn):
    """
    This function will check if a certain dictionary containing docking parameters already exists within the 'docking_parameters.db' located in the 'docking_params_registry_path'.
    
    If the same 'params_dict' already exists, the a flag 1 will be returned.
    
    ------
    Parameters: 
    ------
    - docking_params_registry_path: the full path to the place were 'docking_parameters.db' is located.
    - custom_dict: the dictionary of parameters to be searched in 'docking_parameters.db'.
    - conn: an sql connection to the 'docking_parameters.db' database.
    """

    sql = """SELECT *
             FROM docking_params;   
            """
    df = pd.read_sql_query(sql,conn)
    
    custom_dict_str = str(custom_dict) # Convert the custo
       
    for index, row in df.iterrows():
        current_dict = eval(row[1]) # Here the eval() function evaluates the input string (stored dict) and returns a dictionary
                
        if custom_dict == current_dict: # This will check if the custom_dict already exists
            print("hallado")
            print(colored(f"The custom_dictionary already exists under the index: {row[0]}. Exiting.","red"))
            exit()
    
def store_custom_param_dict(docking_params_registry_path):
    """ 
    This function will store a custom parameters directory into the 'docking_parameters.db' stored within 'docking_params_registry_path'
    
    ------
    Parameters:
    ------
    - docking_params_registry_path: the full path to the place were 'docking_parameters.db' is located.
    
    ------
    Returns:
    ------
    A record within the 'docking_params_registry_path'/'docking_parameters.db' is stored.
    """
    # The following storing of the registry will only occur if the 'check_params_dict_existence()' does not stops the execution.
    try: 
        conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')
        # Get the last value to append the new record
        last_value = get_last_param_record_id(conn)
        # Create the dictionary of values to store starting from the default set
        custom_dict = create_custom_params_dict()
        # Check if the 'custom_dict' is already stored within the database
        check_params_dict_existence(custom_dict,conn)
        description = input(colored("Enter a brief description of the conditions set: ","magenta"))
        
        # Store the corresponding registry to the database
        sql = f"""INSERT INTO docking_params(condition_id,condition_dict,condition_description) 
                  VALUES({last_value+1},"{custom_dict}","{description}");"""
        conn.execute(sql)
        conn.commit()
        print(colored(f"Docking params condition_id {last_value+1} is NEW, stored the registry","green"))
    
    except Exception as error:
        print(error)
        print(colored("Addition of the parameters conditions set FAILED","red"))
        exit()
    
def retrieve_docking_conditions_set(docking_params_registry_path,condition_id):
    """
    This function will return a dictionary object containing
    
    """
    
    conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')

    sql = f"""SELECT *
          FROM docking_params;          
          """
    try:
        df = pd.read_sql_query(sql,conn)
        row = df.loc[df["condition_id"] == condition_id]
        stored_dict = row["condition_dict"].iloc[0]
        conditions_dict = eval(stored_dict)
        return conditions_dict
        
    except Exception as error:
        print(error)
        print(colored(f"The condition_id: {condition_id} does not exists","red"))

def compare_params_dicts(dict_1,dict_2):
    """
    This function will compare to dictionaries, and return the modified key_value pairs ('differences string')in a unified string so as to serve as docking input.

    Typically, the first dictionary is the default dictionary of conditions, while the second one is a custom conditions dictionary.

    ------
    Parameters:
    ------
    - dict_1: the first dictionary that is taken as reference for comparison.
    - dict_2: the second dictionary that is takes as reference for comparison.

    ------
    Returns
    ------
    A 'differences_strins' containing the parameteres that differs from the default in the custom_dict
    """
            
    # This will detect parameters differences between 
    differences_string = ''
    
    for key in dict_1.keys():
        for sub_key in dict_1[key].keys():
            if dict_1[key][sub_key] != dict_2[key][sub_key]:
                differences_string = differences_string + f'{sub_key} {dict_2[key][sub_key]} '

    return differences_string


#################

if __name__ == '__main__':

    print("This is a test")
    
    #create_params_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers")
    
    #create_default_dict()
    
    #create_custom_params_dict("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers")
    
    #last_value = get_last_param_record_id("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers")
    
    #store_custom_param_dir("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers")

    #retrieve_docking_conditions_set("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers", 1)
    
    #compare_dictionaries("dymmy1","dymmy2")

    #clear_params_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers")

    #drop_params_registry_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_params_registers")