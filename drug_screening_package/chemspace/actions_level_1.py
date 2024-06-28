import sqlite3
import pandas as pd
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.general_procedures.operations_on_sql_database as sql_ops
import chemspace.input_output_operations.actions_on_dbs as actions_dbs


def create_level_1_project(level_0_db_file,level_1_db_file,reactants_table_list,reactions_table,reactions_list):
    """
    This function will retrieve reactants and reactions in  order to create a new Level_1 project database that will be created in the 'project_db_output_path' path. 
    
    All the required information for the project will be retrieved. 
    
    The function will ask for a corresponding 'proyect description'.
    
    ------
    Parameters:
    ------
    - level_0_db_file: the level_0 database in which reactants and reactions are stored.
    - level_1_db_file: the name of the level_1_db that will be created
    - reactants_table_list: provide a list containing items corresponding to the reactants to retrieve form the level_0_db. E.g. ['reactant_1', 'reactant_2',...]
    - reactions_table: the name of the table in the 'level_0_db_file' in which the reactions to be used in the project are stored.
    - reactions_list: provide a list containing items corresponding to the reaction SMARTS to retrieve form the level_0_db. E.g. ['smarts_1', 'smarts_2',...]
    ------
    Returns
    ------
    A new database for a specific reasearch project is created, contaning molecules and reaction tables.
    """
    
    for reactant_table in reactants_table_list:
        actions_dbs.copy_table_between_dbs(level_0_db_file,level_1_db_file,reactant_table)

    counter = 1
    for reaction in reactions_list:
        df = actions_dbs.extract_row_by_value(level_0_db_file,reactions_table,"Reaction",reaction)
        if counter == 1:
            sql_ops.store_df_to_sql(level_1_db_file,df,"reaction_types","replace")
        else:
            sql_ops.store_df_to_sql(level_1_db_file,df,"reaction_types","append")
        counter+=1
    print("Reactions copied SUCCESFULLY")

###############
if __name__ == "__main__":
    print("This is a test") 
    
    #create_level_1_project("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/level_1_projects/new_sample_level1.db",["alpha_aminoacid","alcohols"],"reaction_types",["azide_from_primary_amine","acylation"])
