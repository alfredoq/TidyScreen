
import sqlite3
import os
import pandas as pd
import rdkit.Chem as Chem
from termcolor import colored
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D
from rdkit.Chem import rdCIPLabeler
import sys
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import chemspace.general_procedures.objects_processing as obj_proc


def generate_depiction_grid(db_name, table_name, output_path, max_mols_ppage):
    """
    This function will create one/multiple .png documents depicting grids of molecules corresponding to the table and column of SMILES read as input.

    The number of generated .png files will corredpond to the ratio: SMILES rows / max_mols_ppage

    In order to run this function, the table provide for plotting has to contain at least two columns named: "SMILES" and "inchi_key".

    The SMILES notation will be used for plotting purposes, while 'inchi_key' will serve as compound label.
    ------
    Parameters
    ------
    - db_name: the full path of the database cotaining the SMILES to be depicted in the grid.
    - table_name: the table included in the 'db-name' that contains the SMILES strings to depict.
    - output_path: the full path were the generated .png files of the depicted images will be saved.    
    - max_mols_ppage: the maximum number of molecules that will be depipected per page.

    ------
    Returns 
    ------
    A serie of .png files that are stored in the 'output_path'
    """
    
    try: 

        # First check if the destination directory alredy exists. If True, exit informing.
        if os.path.isdir('new_folder'):
            print(colored(f"The destination folder '{output_path}' already exists. EXITING","red"))
            exit()
        else:

            # Create the destination directory
            os.makedirs(output_path)            
            print(colored(f"SUCESSFULLY created the output folder: {output_path}","green"))
            # Continue with the processing of the data
            conn = sqlite3.connect(db_name)
            sql = f'''SELECT SMILES, inchi_key
                      FROM {table_name};'''

            df_query = pd.read_sql_query(sql,conn)
            smiles_list = df_query["SMILES"].to_list()
            labels_list = df_query["inchi_key"].to_list()
            indexes_list = [f'cpd_row_id: {x+1}' for x in range(len(smiles_list))]
        
            smiles_list_chunks = obj_proc.split_list_in_chunks(smiles_list, max_mols_ppage)
            labels_list_chunks = obj_proc.split_list_in_chunks(labels_list, max_mols_ppage)
            indexes_list_chunks = obj_proc.split_list_in_chunks(indexes_list, max_mols_ppage)
    
            counter = 0 # used to number figures with chunks
            for index, list_item in enumerate(smiles_list_chunks):
                df_chunk  = pd.DataFrame({'SMILES':list_item})
                df_chunk.to_csv(f"{output_path}/{table_name}_{counter}.csv",index=None,header=None)
                mols_item = [Chem.MolFromSmiles(smile) for smile in list_item]
                img = MolsToGridImage(mols=mols_item, legends=indexes_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
                img.save(f"{output_path}/{table_name}_{counter}.png")
                counter+=1
    
            print(colored(f"SUCESSFULLY saved images in folder: {output_path}","green"))

    except Exception as error:
        print(error)
        print(colored(f"ERROR depicting molecules in table: '{table_name}'.","red"))

def generate_depiction_grid_mol_id(db_name, table_name, miscellaneous_files, max_mols_ppage):
    """
    This function will create one/multiple .png documents depicting grids of molecules corresponding to the table and column of SMILES read as input. The label of the compound will be assigned from the 'mol_id' column.

    The number of generated .png files will corredpond to the ratio: SMILES rows / max_mols_ppage

    In order to run this function, the table provide for plotting has to contain at least two columns named: "SMILES" and "inchi_key".

    The SMILES notation will be used for plotting purposes, while 'inchi_key' will serve as compound label.
    ------
    Parameters
    ------
    - db_name: the full path of the database cotaining the SMILES to be depicted in the grid.
    - table_name: the table included in the 'db-name' that contains the SMILES strings to depict.
    - miscellaneous_files: the full path to the misc files (provided by the general script)
    - max_mols_ppage: the maximum number of molecules that will be depipected per page.

    ------
    Returns 
    ------
    A serie of .png files that are stored in the 'output_path'
    """
    output_path = f'{miscellaneous_files}/{table_name}'

    try: 

        # First check if the destination directory alredy exists. If True, exit informing.
        if os.path.isdir('new_folder'):
            print(colored(f"The destination folder '{output_path}' already exists. EXITING","red"))
            exit()
        else:
            # Create the destination directory
            os.makedirs(output_path)            
            print(colored(f"SUCESSFULLY created the output folder: {output_path}","green"))
            # Continue with the processing of the data
            conn = sqlite3.connect(db_name)
            sql = f'''SELECT mol_id, SMILES, inchi_key
                      FROM '{table_name}';'''

            df_query = pd.read_sql_query(sql,conn)
            smiles_list = df_query["SMILES"].to_list()
            labels_list = df_query["mol_id"].to_list()
            inchi_labels_list = df_query["inchi_key"].to_list()
            # This will generate a list of labels only with mol_id to be depicted
            labels_list_str = [f'mol_id: {x}' for x in (labels_list)]
            # This will generate a list of labels combining mol_id and the corresponding inchi_key to be depicted
            inchi_labels_list_str= [f'mol_id: {x} - inchi_key: {y}' for x,y in zip(labels_list,inchi_labels_list)]
            
            indexes_list = [f'cpd_row_id: {x+1}' for x in range(len(smiles_list))]
        
            smiles_list_chunks = obj_proc.split_list_in_chunks(smiles_list, max_mols_ppage)
            # Return chunks of labels containing only the mol_id
            labels_list_chunks = obj_proc.split_list_in_chunks(labels_list_str, max_mols_ppage)
            # Return chunks of labels containing both the mol_id and inchi_key value
            inchi_labels_list_chunks = obj_proc.split_list_in_chunks(inchi_labels_list_str, max_mols_ppage)
               
            counter = 0 # used to number figures with chunks
            for index, list_item in enumerate(smiles_list_chunks):
                df_chunk  = pd.DataFrame({'SMILES':list_item})
                df_chunk.to_csv(f"{output_path}/{table_name}_{counter}.csv",index=None,header=None)
                mols_item = [Chem.MolFromSmiles(smile) for smile in list_item]
                ## This will output the plot containing molecules only with mol_id label 
                #img = MolsToGridImage(mols=mols_item, legends=labels_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
                ## This will output the plot containing molecules both with mol_id and inchi_key labels
                img = MolsToGridImage(mols=mols_item, legends=inchi_labels_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
                img.save(f"{output_path}/{table_name}_{counter}.png")
                counter+=1
    
            print(colored(f"SUCESSFULLY saved images in folder: {output_path}","green"))

    except Exception as error:
        print(error)
        print(colored(f"ERROR depicting molecules in table: '{table_name}'.","red"))

def plot_table_chemspace(db_name,table_name,output_file):
    """
    This function will generate a plot file (.png) in which the x-axis corresponds to the the UMAP-1 value of the chemspace, while the y-axis corresponds to UMAP-2. The corresponding file will be saved in the folder were the db_name lives
    --- Parameter ---
    - db_name: the full path the the database in which the table to be analyzed lives. 
    - table_name: the name of the table in which the two axes (UMAP-1 and UMAP-2) exists for plotting purposes.
    - output_file: the full path to the .png file showing the chemical space.

    --- Returns ---
    - png file: will create a file named 'chemspace.png' in the analyzed database folder. 
    """

    conn = sqlite3.connect(db_name)
    sql = f'''SELECT *
            FROM {table_name};'''

    try:
        df = pd.read_sql_query(sql,conn)
        figure(figsize=(15, 10), dpi=300)
        plt.scatter(df["UMAP-1"],df["UMAP-2"])
        plt.grid()
        plt.title("Chemical space plot")
        plt.xlabel("UMAP-1")
        plt.ylabel("UMAP-2")
        plt.savefig(output_file)
        
    except:
        print("ERROR in retrieving the chemspace info. Does is exists in the db table?")

##########################################
## Generate the tests for the functions ##
##########################################

if __name__ == "__main__":

    print("Running tests")

    #db_name = "/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/level_1_projects/sample_level1.db"
    #db_name = "/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db"
    #db_name = "/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db"

    #generate_depiction_grid(db_name,'para_di_subs_benzotriazole',"SMILES","inchi_key",25,"/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/compounds_depictions")

    #plot_table_chemspace(db_name,"alpha_aminoacid")
