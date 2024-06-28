import sqlite3
import pandas as pd
import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

def read_reactants(db_name,reactant_table):
    """
    This function  will read a table containing reactants in a database, and return a df for further use in the workflow.
    """
    
    conn = sqlite3.connect(db_name)
    query = pd.read_sql_query(f"SELECT * FROM {reactant_table}")