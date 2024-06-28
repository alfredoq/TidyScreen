import os
from termcolor import colored
import sys
import rdkit.Chem as Chem
from rdkit.Chem.Draw import MolsToGridImage, rdMolDraw2D
from rdkit.Chem import rdCIPLabeler
import pandas as pd

MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.general_procedures.objects_processing as obj_proc

def generate_depiction_grid_with_results_labels(output_path,df,max_mols_ppage):
    """
    Given a dataframe contining docking results, this function will will create grid depictions in the 'docking_misc_data_path' path under the subfolded 'docking_assay_{#}'
    
    """
    
    # First check if the destination directory already exists. If True, exit informing.
        
    if os.path.isdir(output_path):
        print(colored(f"The destination folder '{output_path}' already exists. EXITING","red"))
        exit()
    else:
        # Create the destination directory
        os.makedirs(output_path)            
        print(colored(f"SUCESSFULLY created the output folder: {output_path}","green"))
    
        #ligname_list = df["LigName"].to_list() 
        ligname_list = df["sub_pose"].to_list()
        smiles_list = df["ligand_smile"].to_list()
        clusters_list = df["cluster_size"].to_list()
        score_list = df["docking_score"].to_list()
        
        labels_list = []
        counter=1
        for id, cluster in enumerate(clusters_list):
            label = f'rank_{counter}_{ligname_list[id]} \n size= {cluster} - score: {score_list[id]}' 
            labels_list.append(label)
            counter+=1
            
        
        smiles_list_chunks = obj_proc.split_list_in_chunks(smiles_list, max_mols_ppage)
        labels_list_chunks = obj_proc.split_list_in_chunks(labels_list, max_mols_ppage)
        
        counter = 0 # used to number figures with chunks
        for index, list_item in enumerate(smiles_list_chunks):
            df_chunk  = pd.DataFrame({'SMILES':list_item})
            df_chunk.to_csv(f"{output_path}/molecules_chunk_{counter}.csv",index=None,header=None)
            mols_item = [Chem.MolFromSmiles(smile) for smile in list_item]
            img = MolsToGridImage(mols=mols_item, legends=labels_list_chunks[index], molsPerRow=5,subImgSize=(400,250))
            img.save(f"{output_path}/molecules_chunk_{counter}.png")
            counter+=1
    
        print(colored(f"SUCESSFULLY saved images in folder: {output_path}","green"))
        
        
        
    