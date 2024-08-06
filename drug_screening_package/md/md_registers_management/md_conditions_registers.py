from termcolor import colored
import sqlite3
import pandas as pd

def create_md_cond_registry_table(md_params_registry_path):
    """
    This function will create a database named 'md_params.db', which contains a table named 'md_params'.
    
    This database is intended to store all the parameters required by Amber in order to perform the intended molecular dynamics simulations.
    
    The 'md_params' contains the following structure:
    
        - condition_id: an integer describing the set of conditions to be used
        - condition_dict: a dictionary containing all the set of parameters to be used during the MD simulations.
        - condition_description: a brief description of the parameters set
    
    ------
    Parameters:
    ------
    - md_params_registry_path: the full path to the folder in which the 'md_conditions.db' database is stored.
    
    ------
    Returns:
    ------
    Creates a database and table to store the simulation conditions to be applied to the molecular dynamics assays.
       
    """
    
    try: 
        conn = sqlite3.connect(f'{md_params_registry_path}/md_params.db')

        ins = f'CREATE TABLE md_params ("condition_id" INTEGER, "condition_dict" TEXT, condition_description TEXT)'
        conn.execute(ins)
        print(colored("Docking conditions registers table SUCCESSFULLY created.","green"))
        
    except Exception as error:
        print(error)
        print(colored("Error creating the database 'md_params.db'","red"))
        
def get_last_md_cond_id(conn):
    """
    This function will return the last value stored under the column 'condition_id' as stored in the 'md_params.db' located in  'md_params_registry_path'.
    
    ------
    Parameters:
    ------
    - md_params_registry_path: the full path to the place were 'md_params.db' is located.
    - conn: an sql connection to the 'md_params.db' database.
    ------
    Returns:
    ------
    An integer corresponding to the last value stored under the column 'condition_id'.
    """

    try:
        #conn = sqlite3.connect(f'{docking_params_registry_path}/docking_parameters.db')
        sql_assay = """SELECT condition_id 
                       FROM md_params;"""
    
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

def create_default_md_cond_dict():
    """
    This function will create and return dictionary containing the standard simulation conditions as required for a typical explicit solvent MD simulation excuted using pmemd.cuda included in the Amber software package.
    
    The corresponding dictionary will be returned for use as default set conditions or for further modifications by the corresponding parameters modifications function.
    
    ------
    Parameters:
    ------
    None
    
    ------
    Returns:
    ------
    A dictionary containing all available parameters with default values.
    """      

    general_opts_dict = {'lig_antechamber_parms':{},
                         'tleap_params':{},
                         'min1_params' : {},
                         'min2_params' : {},
                         'heat1_params' : {},
                         'heat2_params' : {},
                         'equi_params' : {},
                         'prod_params' : {},
                        }
    
    lig_antechamber_parms = {'charge_system' : 'gas',
                             'clean_files' : 'yes',
                             'residue_name' : 'UNL',
                        }
    tleap_params = {'protein_forcefield': 'leaprc.protein.ff14SB',
                 'ligand_forcefield':'leaprc.gaff2',
                 'water_forcefield': 'leaprc.water.tip3p',
                 'solvate_geometry': 'SolvateOct',
                 'solvent_model': 'TIP3PBOX',
                 'solvent_box_min_dist': 10,
                 'solvent_min_dist': 1.0,
                 }
    
    min1_params = {'imin' : 1 ,
                  'maxcyc' : 5000,
                  'ncyc' : 2500,
                  'ntb' : 1,
                  'ntr' : 1,
                  'cut' : 10.0,
                  'restraintmask' : f"(:* & !:WAT)",
                  'restraint_wt' : 500
                  }
    
    min2_params = {'imin' : 1 ,
                  'maxcyc' : 5000,
                  'ncyc' : 2500,
                  'ntb' : 1,
                  'ntr' : 0,
                  'cut' : 10.0,
                  }
    
    heat1_params = {'mds' : {'imin' : 0 ,
                             'irest' : 0, 
                             'ntx' : 1, 
                             'ntb' : 2, 
                             'ntp' : 1, 
                             'cut' : 10.0, 
                             'ntr' : 1, 
                             'ntc' : 2, 
                             'ntf' : 2, 
                             'tempi' : 0.0, 
                             'temp0' : 310.0, 
                             'ntt' : 3, 
                             'gamma_ln' : 1.0, 
                             'nstlim' : 200000, 
                             'dt' : 0.002, 
                             'ntpr' : 100, 
                             'ntwx' : 100, 
                             'ntwr' : 1000, 
                             'restraintmask' : f"(:* & !:WAT & !:Na+ & !:Cl-)", 
                             'restraint_wt' : 50.0, 
                             'nmropt' : 1,
                             },
                    'temp_change' : {'TYPE' : 'TEMP0',
                                     'ISTEP1' : 0,
                                     'ISTEP2' : 100000,
                                     'VALUE1' : 0.0,
                                     'VALUE2' : 310.0,
                                    }
                  }
    
    heat2_params = {'mds' : {'imin' : 0 ,
                             'irest' : 1, 
                             'ntx' : 5, 
                             'ntb' : 2, 
                             'ntp' : 1, 
                             'cut' : 10.0, 
                             'ntr' : 1, 
                             'ntc' : 2, 
                             'ntf' : 2, 
                             'tempi' : 310.0, 
                             'temp0' : 310.0, 
                             'ntt' : 3, 
                             'gamma_ln' : 1.0, 
                             'nstlim' : 100000, 
                             'dt' : 0.002, 
                             'ntpr' : 100, 
                             'ntwx' : 100, 
                             'ntwr' : 1000, 
                             'restraintmask' : f"(@CA,C,N,O & !:WAT | :{lig_antechamber_parms['residue_name']})", 
                             'restraint_wt' : 50.0, 
                             },
                  }
    
    equi_params = {'mds' : {'imin' : 0 ,
                             'irest' : 1, 
                             'ntx' : 5, 
                             'ntb' : 1, 
                             'ntp' : 0, 
                             'cut' : 10.0, 
                             'ntr' : 1, 
                             'ntc' : 2, 
                             'ntf' : 2, 
                             'tempi' : 310.0, 
                             'temp0' : 310.0, 
                             'ntt' : 3, 
                             'gamma_ln' : 1.0, 
                             'nstlim' : 100000, 
                             'dt' : 0.002, 
                             'ntpr' : 100, 
                             'ntwx' : 100, 
                             'ntwr' : 1000, 
                             'restraintmask' : f"(@CA,C,N,O & !:WAT)", 
                             'restraint_wt' : 50, 
                             },
                  }
    
    prod_params = {'mds' : {'imin' : 0 ,
                             'irest' : 1, 
                             'ntx' : 5, 
                             'ntb' : 1, 
                             'ntp' : 0, 
                             'cut' : 10.0, 
                             'ntr' : 1, 
                             'ntc' : 2, 
                             'ntf' : 2, 
                             'tempi' : 310.0, 
                             'temp0' : 310.0, 
                             'ntt' : 3, 
                             'gamma_ln' : 1.0, 
                             'nstlim' : 5000000, 
                             'dt' : 0.002, 
                             'ntpr' : 100, 
                             'ntwx' : 100, 
                             'ntwr' : 1000, 
                             'restraintmask' : f"(@CA,C,N,O & !:WAT)", 
                             'restraint_wt' : 50, 
                             },
                  }
    
    # Configure the general options dictionary with the custom values
    general_opts_dict['lig_antechamber_parms'] = lig_antechamber_parms
    general_opts_dict['tleap_params'] = tleap_params
    general_opts_dict['min1_params'] = min1_params
    general_opts_dict['min2_params'] = min2_params
    general_opts_dict['heat1_params'] = heat1_params
    general_opts_dict['heat2_params'] = heat2_params
    general_opts_dict['equi_params'] = equi_params
    general_opts_dict['prod_params'] = prod_params
    

    return general_opts_dict

def check_param_exists(parameter):
    """
    This function will check if a given parameter exists within the set of parameters available in parameters dictionary available to perform MD simulations.
    
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
    default_dict = create_default_md_cond_dict()

    for key in default_dict.keys():
        for sub_key in default_dict[key].keys():
            # Check if the parameter matches
            if parameter == sub_key:
                print(colored(f"Parameter '{parameter}' found in category '{key}'","green"))
                return key, sub_key, 1
            else:
                continue
        
    return 0,0,0 # if the parameter does not exist, all zero values will be returned

def create_custom_md_cond_dict(md_params_registry_path):
    """
    This function will store a new registry in the params_registry table.
    
    A key to be modified is requested, after which the corresponding parameter is requested.
    
    The base params_dict is created using the 'create_default_dict', after which the requested parameters are modified.
    
    ------
    Parameters:
    ------
    - md_params_registry_path: the full path to the folder in which the parameters database is stored.
    
    ------
    Returns
    ------
    A new registry is written to the registry table only in case the created parameters set is new.
    """ 
    
    # This will create a base params dictionary to be modified
    custom_dict = create_default_md_cond_dict()
    
    while True:
        parameter = input("Enter the parameter to modify (enter 'DONE' to finish): ")
        
        if parameter == 'DONE':
            print(colored("Finished entering parameters.","green"))
            break
        # Check if the flag exists
        key,sub_key,check_flag = check_param_exists(parameter)
        
        if check_flag == 0:
            print(colored(f"The parameter '{parameter}' DOES NOT EXIST. Enter the param again.","red"))
            continue
        elif check_flag == 1:
            value = input(f"Provide '{parameter}' parameter value: ")
            custom_dict[key][sub_key]=value
            continue
        
    return custom_dict

def check_params_dict_existence(custom_dict,conn):
    """
    This function will check if a certain dictionary containing MD parameters already exists within the 'md_params.db' located in the 'md_params_registry_path'.
    
    If the same 'params_dict' already exists, the a flag 1 will be returned.
    
    ------
    Parameters: 
    ------
    - custom_dict: the dictionary of parameters to be searched in 'md_params.db'.
    - conn: an sql connection to the 'md_params.db' database.
    """

    sql = """SELECT *
             FROM md_params;   
            """
    df = pd.read_sql_query(sql,conn)
    
    custom_dict_str = str(custom_dict) # Convert the custo
       
    for index, row in df.iterrows():
        current_dict = eval(row[1]) # Here the eval() function evaluates the input string (stored dict) and returns a dictionary
                
        if custom_dict == current_dict: # This will check if the custom_dict already exists
            print(colored(f"The custom_dictionary already exists under the index: {row[0]}. Exiting.","red"))
            exit()

def append_md_cond_registry(md_params_registry_path):
    """
    This function will store a custom md parameters conditions set into the 'md_params.db' stored within 'md_params' table
    
    ------
    Parameters:
    ------
    - md_params_registry_path: the full path to the place were 'md_params.db' is located.
    
    ------
    Returns:
    ------
    A record within the 'md_params.db' is stored.
    
    """
    
    # The following storing of the registry will only occur if the 'check_params_dict_existence()' does not stops the execution.
    try: 
        conn = sqlite3.connect(f'{md_params_registry_path}/md_params.db')
        # Get the last value to append the new record
        last_value = get_last_md_cond_id(conn)
    
        # Create the dictionary of values to store starting from the default set
        custom_dict = create_custom_md_cond_dict(md_params_registry_path)
        
        # Check if the 'custom_dict' is already stored within the database
        check_params_dict_existence(custom_dict,conn)
    
        # Provide a brief description of the parameters set
        description = input(colored("Enter a brief description of the conditions set: ","magenta"))
        
        # Store the corresponding registry to the database
        sql = f"""INSERT INTO md_params(condition_id,condition_dict,condition_description) 
                  VALUES({last_value+1},"{custom_dict}","{description}");"""
        conn.execute(sql)
        conn.commit()
        print(colored(f"Docking params condition_id {last_value+1} is NEW, stored the registry","green"))
        
    except Exception as error:
        print(error)
        print(colored("Addition of the parameters conditions set FAILED","red"))
        exit()
        

