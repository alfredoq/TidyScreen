from ringtail import RingtailCore
from pprint import pprint
from termcolor import colored
import sqlite3
import pandas as pd
import os
import MDAnalysis as mda
import prolif as plf
from glob import glob
from pandarallel import pandarallel # import pandarallel
from openbabel import openbabel
import shutil
import numpy as np
import pickle
from sklearn.ensemble import RandomForestClassifier
import xtarfile as tarfile
from pathlib import Path
import matplotlib.pyplot as plt 

import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
from chemspace.general_procedures import operations_on_sql_database as sql_ops
from docking.docking_results_processing import results_depiction_procedures as depict_res
from docking.ligands_management import process_ligands_table as lig_proc
from md.perform_md import md_configuration as md_conf

def create_db_from_dlgs(docking_assays_storing_path, assay_nbr, max_poses):
    """
    This function will create an sql database from a folder containing all the docking .dlg files.

    ------
    Parameters:
    ------
    - docking_assays_storing_path: the full path to the folder in which all docking assays are run and stored in folders of type: 'docking_assay_#'.
    - assay_nbr: the number of the docking assays to be processed.
    - max_poses : the maximum number of poses clusters to retain from the docking file.

    ------
    Returns:
    ------
    A database named 'docking_assay_#' as generated after processing with ringtail.
    """

    # Check if the results database alredy exists. If True, stop execution
    if os.path.isfile(f'{docking_assays_storing_path}/docking_assay_{assay_nbr}'):
        print(colored(f"The database file: docking_assay_{assay_nbr} already exists. Exiting.","red"))
        exit()

    try:
        opts = RingtailCore.get_defaults()
        #Modify the corresponding values as needed
        opts['rman_opts']['values']['file_sources']['file_path']['path'] = [[f'{docking_assays_storing_path}/docking_assay_{assay_nbr}']]
        opts['rman_opts']['values']['file_sources']['file_path']['recursive'] = True
        opts['rman_opts']['values']['max_poses'] = max_poses
        opts['storage_opts']['values']['db_file'] = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
        opts['storage_opts']['values']['output_all_poses'] = True
        with RingtailCore(opts_dict=opts) as rt_core:
            rt_core.add_results()

    except Exception as error:
        print(error)
        print(colored("Error creating docking results database","red"))

    # Rename the corresponding docking subposes
    try:
        print(colored("Naming indidual subposes within table","green"))
        add_docking_subposes_nbr(docking_assays_storing_path, assay_nbr)
        print(colored(f"\n SUCCESFULLY created database: docking_assay_{assay_nbr}.db","green"))
    except Exception as error:
        print(error)
        print(colored("Problem renaming the docking subposes in database: docking_assay_{assay_nbr}.db","red"))

def create_db_from_dlgs_custom_path(dest_folder):
    """
    This function will create an sql database from a folder containing all the docking .dlg files in a provided path

    """

    try:
        opts = RingtailCore.get_defaults()
        #Modify the corresponding values as needed
        opts['rman_opts']['values']['file_sources']['file_path']['path'] = [[f'{dest_folder}']]
        opts['rman_opts']['values']['file_sources']['file_path']['recursive'] = True
        opts['rman_opts']['values']['max_poses'] = 10
        opts['storage_opts']['values']['db_file'] = f'{dest_folder}/results.db'
        opts['storage_opts']['values']['output_all_poses'] = True
        with RingtailCore(opts_dict=opts) as rt_core:
            rt_core.add_results()

    except Exception as error:
        print(error)
        print(colored("Error creating docking results database","red"))

    # Rename the corresponding docking subposes
    try:
        print(colored("Naming indidual subposes within table","green"))
        add_docking_subposes_nbr_custom_path(dest_folder)
        print(colored(f"\n SUCCESFULLY created database: results.db","green"))
    except Exception as error:
        print(error)
        print(colored("Problem renaming the docking subposes in database: docking_assay_{assay_nbr}.db","red"))

def add_docking_subposes_nbr(docking_assays_storing_path, assay_nbr):
    """
    This function will process an already prepared database containing docking results and that was processed using ringtail.

    The objective of this function is, given the 'Results' table in the resulting database, which contains subposes of a same compound name stored in the 'LigName' field, this will add a column named 'sub_pose_nbr' to the name.

    The values corresponding to the 'sub_pose_nbr' column corresponds to the unique identification of a molecule pose within the whole molecular docking assay.

    ------
    Parameters:
    ------
    - docking_assays_storing_path: the full path to the folder in which all docking assays are run and stored in folders of type: 'docking_assay_#'.
    - assay_nbr: the number of the docking assays to be processed.

    ------
    Returns:
    ------
    A modified 'LigName' field as stored in the databased corresponding to results processed with ringtail.
    """    
    # Connect to the target database
    try:
        database_file = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
        conn = sqlite3.connect(database_file)
    except Exception as error:
        print(error)
        print(colored(f"Error conecting to database: docking_assay{assay_nbr}.db","red"))

    # Retrieve the 'Results' table as a dataframe for internal processing
        
    try:
        sql = """SELECT *
               FROM Results;"""
        
        df = pd.read_sql(sql,conn)

        # Add the '_#' value to the corresponding name

        unique_names_list = []
        for index, row in df.iterrows():
            ligname = row['LigName']
            if ligname not in unique_names_list:
                counter=1
                df.at[index,'sub_pose'] =  ligname+'_'+str(counter)
                unique_names_list.append(ligname)
            else:
                counter+=1 
                df.at[index,'sub_pose'] =  ligname+'_'+str(counter)

        # Reorder df columns and store to the database replacing the original one
        df = df[['Pose_ID', 'LigName', 'sub_pose','receptor', 'pose_rank', 'run_number','docking_score', 'leff', 'deltas', 'cluster_rmsd', 'cluster_size','reference_rmsd', 'energies_inter', 'energies_vdw', 'energies_electro','energies_flexLig', 'energies_flexLR', 'energies_intra','energies_torsional', 'unbound_energy', 'nr_interactions', 'num_hb','about_x', 'about_y', 'about_z', 'trans_x', 'trans_y', 'trans_z',
        'axisangle_x', 'axisangle_y', 'axisangle_z', 'axisangle_w', 'dihedrals','ligand_coordinates', 'flexible_res_coordinates']]
        
        sql_ops.store_df_to_sql(database_file,df,"Results","replace")

    except Exception as error:
        print(error)

def add_docking_subposes_nbr_custom_path(dest_folder):
    """
    This function will process an already prepared database containing docking results and that was processed using ringtail.

    """    
    # Connect to the target database
    try:
        database_file = f'{dest_folder}/results.db'
        conn = sqlite3.connect(database_file)
    except Exception as error:
        print(error)
        print(colored(f"Error conecting to database: results.db in {dest_folder}","red"))

    # Retrieve the 'Results' table as a dataframe for internal processing
        
    try:
        sql = """SELECT *
               FROM Results;"""
        
        df = pd.read_sql(sql,conn)

        # Add the '_#' value to the corresponding name

        unique_names_list = []
        for index, row in df.iterrows():
            ligname = row['LigName']
            if ligname not in unique_names_list:
                counter=1
                df.at[index,'sub_pose'] =  ligname+'_'+str(counter)
                unique_names_list.append(ligname)
            else:
                counter+=1 
                df.at[index,'sub_pose'] =  ligname+'_'+str(counter)

        # Reorder df columns and store to the database replacing the original one
        df = df[['Pose_ID', 'LigName', 'sub_pose','receptor', 'pose_rank', 'run_number','docking_score', 'leff', 'deltas', 'cluster_rmsd', 'cluster_size','reference_rmsd', 'energies_inter', 'energies_vdw', 'energies_electro','energies_flexLig', 'energies_flexLR', 'energies_intra','energies_torsional', 'unbound_energy', 'nr_interactions', 'num_hb','about_x', 'about_y', 'about_z', 'trans_x', 'trans_y', 'trans_z',
        'axisangle_x', 'axisangle_y', 'axisangle_z', 'axisangle_w', 'dihedrals','ligand_coordinates', 'flexible_res_coordinates']]
        
        sql_ops.store_df_to_sql(database_file,df,"Results","replace")

    except Exception as error:
        print(error)

def create_table_with_lowest_energy_poses(docking_assays_storing_path, assay_nbr):
    """
    This function will process the 'Results' table generated by ringtail analysis on the corresponding .dlg files, an return a table named 'Results_lowest_energy' containing only the docking positions corresponding to the lowest energy cluster.
    
    ------
    Parameters:
    ------
    - docking_assays_storing_path: the full path to the folder in which all docking assays are run and stored in folders of type: 'docking_assay_#'.
    - assay_nbr: the number of the docking assays to be processed.

    ------
    Returns:
    ------
    A table named 'Results_lowest_energy' within the dockig assay database.
    """
    assay_db = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    
    try: 
        conn = sqlite3.connect(assay_db)

        sql = """SELECT * 
                 FROM Results;"""
    
        df = pd.read_sql_query(sql,conn)
        df.drop_duplicates('LigName', inplace=True, keep="first")

        sql_ops.store_df_to_sql(assay_db,df,"Results_lowest_energy","replace")

    except Exception as error:
        print(error)
        print(colored("Error preparing the table with lowest energy clusters","red"))

def print_ringtail_default_ops():
    """"
    This simple function will only print the available option in the Ringtail opts_dict
    """
    opts = RingtailCore.get_defaults()
    pprint(opts)

def split_sdf_file(output_path,filename):
    """
    This function will split a multimolecule .sdf file into submolecules

    ------
    Parameters:
    ------
    - filename: the name of the .sdf file to be splitted into the corresponding molecules.

    ------
    Returns:
    ------
    Multiple .sdf subfiles with the corresponding splitted molecules.
    """
    filename_prefix = filename.split('/')[-1].replace('.sdf','')
    string = ''
    counter = 1
    for line in open(filename):
        string+=line
        if '$$$$' in line:
            counter+=1
            string = ''
        output = filename_prefix+"_"+str(counter)+".sdf"
        output_file = f'{output_path}/{output}'
        with open(output_file, "w") as file:
            file.write(string)
            file.close()
    
    # This code will delete the last file generated by the loop, since the multimolecule .sdf file always ends with a '$$$$' string activating the writting of an empty file that needs to be deleted.
    os.remove(output_file)
    # This will delete the original multimolecule filename
    os.remove(filename)
    
def extract_selected_index_from_sdf_file(ligname,database,output_path,index,rank):
    """
    This function will split a multimolecule .sdf file into submolecules

    ------
    Parameters:
    ------
    - filename: the name of the .sdf file to be splitted into the corresponding molecules.

    ------
    Returns:
    ------
    Multiple .sdf subfiles with the corresponding splitted molecules.
    """
    
    # Extract the corresponding .sdf file from the database
    
    try: 
        opts = RingtailCore.get_defaults()
        #Modify the corresponding values as needed
        opts['filters']['values']['ligand_name'] = [ligname]
        opts['storage_opts']['values']['db_file'] = f'{database}'
        opts['storage_opts']['values']['output_all_poses'] = True
        opts['out_opts']['values']['export_sdf_path'] = f'{output_path}/'
        
        with RingtailCore(opts_dict=opts) as rt_core:
            rt_core.filter()
            rt_core.write_molecule_sdfs()

        # This will split the multimolecules .sdf files
        filename = f'{output_path}/{ligname}.sdf'

    except Exception as error:
        print(error)
       
    filename_prefix = filename.split('/')[-1].replace('.sdf','')
    string = ''
    counter = 1

    for line in open(filename):
        string+=line
        if '$$$$' in line:
            counter+=1
            string = ''
        
        # Disable writing the rank part of the file in case it is not needed. Flag = 0
        
        if rank == 0:
            output = f'{filename_prefix}_{str(counter)}.sdf'
                
        else:
            output = f'rank_{rank}_{filename_prefix}_{str(counter)}.sdf'

        output_file = f'{output_path}/{output}'
        
        if index == counter:
            with open(output_file, "w") as file:
                file.write(string)
                file.close()
    
            # Will save the corresponding filename in to return it in case it is needed for further processing
            extracted_sdf_filename = output_file.split('/')[-1]

    # This will delete the original multimolecule filename
    os.remove(filename)
    # Return the written sdf filename for further use
    
    return extracted_sdf_filename

def extract_single_sdf_molecule(ligname,database,output_path):
    """
    This function will extract a single sdf file based on the corresponding 'molname'
    """
    try: 
        opts = RingtailCore.get_defaults()
        #Modify the corresponding values as needed
        opts['filters']['values']['ligand_name'] = [ligname]
        opts['storage_opts']['values']['db_file'] = f'{database}'
        opts['storage_opts']['values']['output_all_poses'] = True
        opts['out_opts']['values']['export_sdf_path'] = f'{output_path}/'
        
        with RingtailCore(opts_dict=opts) as rt_core:
            rt_core.filter()
            rt_core.write_molecule_sdfs()

        # This will split the multimolecules .sdf files
        full_path_filename = f'{output_path}/{ligname}.sdf'
        
        split_sdf_file(output_path,full_path_filename)
        
    except Exception as error:
        print(error)

def extract_all_sdf_molecules(docking_assays_storing_path,assay_nbr):
    """
    This function will recursively extract .sdf files from a database created from a docking assay.

    Aditionally, this function will create the corresponding .pdb files from the .sdf docked poses.

    ------
    Parameters:
    ------
    - docking_assays_storing_path: the full path to the folder in which all docking assays are run and stored in folders of type: 'docking_assay_#'.
    - assay_nbr: the number of the docking assays to be processed.
    """
    
    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    
    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands'

    Path(output_path).mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(database)

    sql = f"""SELECT LigName
              FROM Ligands;"""

    ligands_name_df = pd.read_sql_query(sql,conn)

    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        extract_single_sdf_molecule(ligname,database,output_path)
    
    # This function will generate the corresponding pdb files from the generated sdf docked poses
    generate_pdb_from_sdf(output_path)

def extract_all_pdb_poses(docking_assays_storing_path,assay_nbr,max_poses):
    """
    This function will extract and stored the docked poses
    
    ------
    Parameters
    ------
    - docking_assays_storing_path: the full path the folder in which all docking assays are stored
    - assay_number: the number of the docking assay to be processed
    - max_poses: the maximum number of docked poses to be processed
    """

    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(database)
    
    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands'
    Path(output_path).mkdir(parents=True, exist_ok=True)
    
    # Retrieve the names of ligand in order to get the .dlg to process and that contains the ligands
    sql = f"""SELECT LigName
              FROM Ligands;"""

    ligands_name_df = pd.read_sql_query(sql,conn)

    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        dlg_file_to_parse = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/dlgs/{ligname}.dlg'
        parse_dlg_file(database,dlg_file_to_parse,max_poses)

def parse_dlg_file(database,dlg_file,max_poses):
    
    """"
    This function will parse a .dlg file in order to generate one .pdb file for each detected cluster
    
    ------
    Parameters:
    ------
    - database: the full path to the database in which the pdb files are to be stored
    - dlg_file: the full path to the .dlg file to be parsed for .pdb docked poses
    - max_poses: the maximun number of poses to be extracted from the .dlg file
    """

    ligand_prefix = dlg_file.split('/')[-1].replace('.dlg','')
    trigger = 0 
    activation_keywords =  ['CLUSTERING', 'HISTOGRAM'] # This is the opening line of the Clustering Histogram
    shutdown_keywords = ['RMSD', 'TABLE']
    number_of_poses_to_extract = 10

    # This loop will construct the list of docked poses to store based on the clustering histogram
    counter = 1
    pose_rank_list = []
    pose_id_list = []
    for line in open(dlg_file):
        new_line = line.rsplit()
    
        if new_line == activation_keywords:
            trigger = 1
        elif new_line == shutdown_keywords:
            trigger = 0
    
        if trigger == 1 and len(new_line) == 10 and counter <= number_of_poses_to_extract:
            pose_rank_list.append(int(new_line[0]))
            pose_id_list.append(new_line[4])
            counter+=1

    # This loop will extract each .pdbqt file of the stored docking pose and convert it into the corresponding pdb file.

    for pose_rank in pose_rank_list:
        output_file = f'/tmp/{ligand_prefix}_{pose_rank}.pdb'
        pose_id = pose_id_list[pose_rank-1]
        trigger = 0
        activation_keywords =  ['DOCKED:', 'USER', 'Run', '=', pose_id] # This is the opening line of the .pdbqt run to extract
        shutdown_keywords = ['DOCKED:', 'ENDMDL']
        with open(output_file,'w') as pdb_output:
            for line in open(dlg_file):
                new_line = line.rsplit()
                if new_line == activation_keywords:
                    trigger = 1
                elif new_line == shutdown_keywords:
                    trigger = 0
        
                if trigger == 1 and new_line[1] == 'ATOM': # This extract the lines constituting the .pdb file
                    # This will output the .pdb file with the corresponding format
                    #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n")
                    atom_field = "HETATM"
                    pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")
        
        pdb_output.close()
        
        # Compress the .pdb file
        # Create the tar file within the /tmp di for storage
        model_name = output_file.split('/')[-1]
        tar = tarfile.open(f"/tmp/{model_name}.tar.gz", "w:gz")
        tar.add(output_file, arcname=model_name)
        tar.close()
        
        pdb_blob_obj = lig_proc.create_blob_object(f"/tmp/{model_name}.tar.gz")
        store_pdb_file_in_database(database,output_file,pdb_blob_obj)
        
        # Clean the generated files in the /tmp dir
        os.remove(output_file)
        os.remove(f"/tmp/{model_name}.tar.gz")

def parse_dlg_file_by_cluster(docking_assays_registry_path,assay_id,receptor_models_registry_path,ligand_main_name,cluster_number):
    """
    This function will parse a .dlg file and extract all the poses in the given cluster in order to generate one .pdb file for each run.
        
    """
    ligand_db,source_ligands_db,source_ligands_pdbqt_table, rec_model_id, receptor_pdb_file = md_conf.retrieve_assay_relevant_variables(docking_assays_registry_path,assay_id,receptor_models_registry_path)

    dlg_path_list = ligand_db.split('/')[0:-1]
    assay_path = '/'.join(dlg_path_list)
    dlg_file = f'{assay_path}/dlgs/{ligand_main_name}.dlg'
    
    # This will get all the runs corresponding to the indicated cluster number
    run_list = []
    for line in open(dlg_file):
        new_line = line.rsplit()
        if "RANKING" in new_line and int(new_line[0]) == cluster_number:
            run_list.append(new_line[2])
    
    ## This loop will extract each .pdbqt file of the stored docking pose and convert it into the corresponding pdb file.

    # Create the directory to store the cluster .pdb files in the assay_dir
    output_dir = f'{assay_path}/pdbs_{ligand_main_name}_cluster_{cluster_number}' 
    isExist = os.path.exists(output_dir)
    if not isExist:
       # Create a new directory because it does not exist
        os.makedirs(output_dir)
    else:
        print(colored("The directory already exists","red"))
        exit()

    counter = 1
    for run in run_list:
        output_file = f'/{output_dir}/{ligand_main_name}_{run}.pdb'
        trigger = 0
        activation_keywords =  ['DOCKED:', 'USER', 'Run', '=', run] # This is the opening line of the .pdbqt run to extract
        shutdown_keywords = ['DOCKED:', 'ENDMDL']
        with open(output_file,'w') as pdb_output:
            for line in open(dlg_file):
                new_line = line.rsplit()
                if new_line == activation_keywords:
                    trigger = 1
                elif new_line == shutdown_keywords:
                    trigger = 0
        
                if trigger == 1 and new_line[1] == 'ATOM': # This extract the lines constituting the .pdb file
                    # This will output the .pdb file with the corresponding format
                    ## pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n") ## Deprecated due to parsing error
                    #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:>1} {new_line[6]:>11}{new_line[7]:>8}{new_line[8]:>8}\n")
                    atom_field = "HETATM"
                    pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")

        pdb_output.close()
        
        # Compress the .pdb file
        # Create the tar file within the /tmp for storage
        file_name = output_file.split('/')[-1].replace('.pdb','')
        tar = tarfile.open(f"/{output_dir}/{file_name}.tar.gz", "w:gz")
        tar.add(output_file, arcname=f'{file_name}.pdb')
        tar.close()
        
        pdb_blob_obj = lig_proc.create_blob_object(f"/{output_dir}/{file_name}.tar.gz")
        
        store_pdb_file_in_database_by_cluster(ligand_db,file_name,pdb_blob_obj,ligand_main_name,cluster_number,counter)
        
        counter+=1

def extract_1_per_cluster_from_dlg(docking_assays_storing_path,assay_nbr):
    """
    This function will parse a .dlg file and will extract 1 .pdb per identified cluster into a folder named 'cluster_pdb_files' located within the docking assay folder
    
    """
    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(database)

    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands_from_dlg_1_per_cluster'
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Retrieve the names of ligand in order to get the .dlg to process and that contains the ligands
    sql = f"""SELECT LigName
              FROM Ligands;"""

    ligands_name_df = pd.read_sql_query(sql,conn)

    activation_keywords =  ['CLUSTERING', 'HISTOGRAM'] # This is the opening line of the Clustering Histogram
    shutdown_keywords = ['RMSD', 'TABLE']
    number_of_poses_to_extract = 999

    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        dlg_file_to_parse = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/dlgs/{ligname}.dlg'
        
        trigger = 0 
        
        # This loop will construct the list of docked poses to store based on the clustering histogram
        counter = 1
        pose_rank_list = []
        pose_id_list = []

        for line in open(dlg_file_to_parse):
            new_line = line.rsplit()
            if new_line == activation_keywords:
                trigger = 1
            elif new_line == shutdown_keywords:
                trigger = 0
    
            if trigger == 1 and len(new_line) == 10 and counter <= number_of_poses_to_extract and "#" in new_line[-1]:
                pose_rank_list.append(int(new_line[0]))
                pose_id_list.append(new_line[4])
                counter+=1        

        # This will extract all the pre-identified pdb files and store them into the indicated folder
        for pose_rank in pose_rank_list:
            output_file = f'{output_path}/{ligname}_{pose_rank}.pdb'
            pose_id = pose_id_list[pose_rank-1]
            trigger_2 = 0
            activation_keywords_2 =  ['DOCKED:', 'USER', 'Run', '=', pose_id] # This is the opening line of the .pdbqt run to extract
            shutdown_keywords_2 = ['DOCKED:', 'ENDMDL']
            with open(output_file,'w') as pdb_output:
                for line in open(dlg_file_to_parse):
                    new_line = line.rsplit()
                    if new_line == activation_keywords_2:
                        trigger_2 = 1
                    elif new_line == shutdown_keywords_2:
                        trigger_2 = 0
        
                    if trigger_2 == 1 and new_line[1] == 'ATOM': # This extract the lines constituting the .pdb file
                        # This will output the .pdb file with the corresponding format
                        #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n")
                        atom_field = "HETATM"
                        pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")
        
            pdb_output.close()

def extract_1_per_cluster_from_dlg_custom_path(dest_folder):
    """
    This function will parse a .dlg file and will extract 1 .pdb per identified cluster into a folder named 'cluster_pdb_files' located within the docking assay folder
    
    """
    # Create a connection to the already created database
    database = f'{dest_folder}/results.db'
    conn = sqlite3.connect(database)

    # Configure the output path to store corresponding .sdf files    
    output_path = f'{dest_folder}/docked_ligands_from_dlg_1_per_cluster'
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Retrieve the names of ligand in order to get the .dlg to process and that contains the ligands
    sql = f"""SELECT LigName
              FROM Ligands;"""

    ligands_name_df = pd.read_sql_query(sql,conn)

    activation_keywords =  ['CLUSTERING', 'HISTOGRAM'] # This is the opening line of the Clustering Histogram
    shutdown_keywords = ['RMSD', 'TABLE']
    number_of_poses_to_extract = 999

    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        dlg_file_to_parse = f'{dest_folder}/dlgs/{ligname}.dlg'
        
        trigger = 0 
        
        # This loop will construct the list of docked poses to store based on the clustering histogram
        counter = 1
        pose_rank_list = []
        pose_id_list = []

        for line in open(dlg_file_to_parse):
            new_line = line.rsplit()
            if new_line == activation_keywords:
                trigger = 1
            elif new_line == shutdown_keywords:
                trigger = 0
    
            if trigger == 1 and len(new_line) == 10 and counter <= number_of_poses_to_extract and "#" in new_line[-1]:
                pose_rank_list.append(int(new_line[0]))
                pose_id_list.append(new_line[4])
                counter+=1        
                
        # This will extract all the pre-identified pdb files and store them into the indicated folder
        for pose_rank in pose_rank_list:
            output_file = f'{output_path}/{ligname}_{pose_rank}.pdb'
            pose_id = pose_id_list[pose_rank-1]
            trigger_2 = 0
            activation_keywords_2 =  ['DOCKED:', 'USER', 'Run', '=', pose_id] # This is the opening line of the .pdbqt run to extract
            shutdown_keywords_2 = ['DOCKED:', 'ENDMDL']
            with open(output_file,'w') as pdb_output:
                for line in open(dlg_file_to_parse):
                    new_line = line.rsplit()
                    if new_line == activation_keywords_2:
                        trigger_2 = 1
                    elif new_line == shutdown_keywords_2:
                        trigger_2 = 0
        
                    if trigger_2 == 1:
                        if new_line[1] == 'ATOM' or new_line[1] == 'HETATM': # This extract the lines constituting the .pdb file
                            # This will output the .pdb file with the corresponding format
                            #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n")
                            atom_field = "HETATM"
                            pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")
        
            pdb_output.close()
            
def extract_lowest_energy_conformation_from_dlg(docking_assays_storing_path,assay_nbr):
    """
    This function will parse all .dlg files in a docking experiment and will extract the lowest energy conformation as a .pdb file for each .dlg file.
    
    """
    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(database)

    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands_from_dlg_lowest_energy'
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Retrieve the names of ligand in order to get the .dlg to process and that contains the ligands
    sql = f"""SELECT LigName
              FROM Ligands;"""

    ligands_name_df = pd.read_sql_query(sql,conn)

    activation_keywords =  ['CLUSTERING', 'HISTOGRAM'] # This is the opening line of the Clustering Histogram
    shutdown_keywords = ['RMSD', 'TABLE']
    number_of_poses_to_extract = 1

    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        dlg_file_to_parse = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/dlgs/{ligname}.dlg'
        
        trigger = 0 
        
        # This loop will construct the list of docked poses to store based on the clustering histogram
        counter = 1
        pose_rank_list = []
        pose_id_list = []

        for line in open(dlg_file_to_parse):
            new_line = line.rsplit()
            if new_line == activation_keywords:
                trigger = 1
            elif new_line == shutdown_keywords:
                trigger = 0
    
            if trigger == 1 and len(new_line) == 10 and counter <= number_of_poses_to_extract and "#" in new_line[-1]:
                pose_rank_list.append(int(new_line[0]))
                pose_id_list.append(new_line[4])
                counter+=1        

        # This will extract all the pre-identified pdb files and store them into the indicated folder
        for pose_rank in pose_rank_list:
            output_file = f'{output_path}/{ligname}_{pose_rank}.pdb'
            pose_id = pose_id_list[pose_rank-1]
            trigger_2 = 0
            activation_keywords_2 =  ['DOCKED:', 'USER', 'Run', '=', pose_id] # This is the opening line of the .pdbqt run to extract
            shutdown_keywords_2 = ['DOCKED:', 'ENDMDL']
            with open(output_file,'w') as pdb_output:
                for line in open(dlg_file_to_parse):
                    new_line = line.rsplit()
                    if new_line == activation_keywords_2:
                        trigger_2 = 1
                    elif new_line == shutdown_keywords_2:
                        trigger_2 = 0
        
                    if trigger_2 == 1 and new_line[1] == 'ATOM': # This extract the lines constituting the .pdb file
                        # This will output the .pdb file with the corresponding format
                        #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n")
                        atom_field = "HETATM"
                        pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")
        
            pdb_output.close()            

def extract_most_populated_cluster_conformation_from_dlg(docking_assays_storing_path,assay_nbr):
    """
    This function will process a docking assay and will return the most populated cluster of conformations by the parsing the .dlg files stored in the assay.
    
    """

    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(database)

    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands_from_dlg_most_populated_clusters'
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Retrieve the names of ligand in order to get the .dlg to process and that contains the ligands
    sql = f"""SELECT LigName
              FROM Ligands;"""

    ligands_name_df = pd.read_sql_query(sql,conn)

    activation_keywords =  ['CLUSTERING', 'HISTOGRAM'] # This is the opening line of the Clustering Histogram
    shutdown_keywords = ['RMSD', 'TABLE']
    
    for index, row in ligands_name_df.iterrows():
        ligname = row['LigName']
        dlg_file_to_parse = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/dlgs/{ligname}.dlg'
        
        trigger = 0 
        
        # This loop will construct the list of docked poses to store based on the clustering histogram
        pose_rank_list = []
        pose_id_list = []
        cluster_size_list = []
        
        for line in open(dlg_file_to_parse):
            new_line = line.rsplit()
    
            if new_line == activation_keywords:
                trigger = 1
            elif new_line == shutdown_keywords:
                trigger = 0
    
            if trigger == 1 and len(new_line) == 10:
                cluster_size_list.append(len(new_line[9]))
                pose_rank_list.append(int(new_line[0]))
                pose_id_list.append(new_line[4])


        df = pd.DataFrame({"Rank":pose_rank_list,
                           "Pose_ID":pose_id_list,
                           "Cluster_size":cluster_size_list})

        sorted_df = df.sort_values("Cluster_size",ascending=False)

        most_populated_pose_id = sorted_df["Pose_ID"].iloc[0]
    
        output_file = f'/{output_path}/{ligname}_pose_{most_populated_pose_id}.pdb'
        trigger = 0
        activation_keywords_2 =  ['DOCKED:', 'USER', 'Run', '=', most_populated_pose_id] # This is the opening line of the .pdbqt run to extract
        shutdown_keywords_2 = ['DOCKED:', 'ENDMDL']
        
        with open(output_file,'w') as pdb_output:
            for line in open(dlg_file_to_parse):
                new_line = line.rsplit()
                if new_line == activation_keywords_2:
                    trigger = 1
                elif new_line == shutdown_keywords_2:
                    trigger = 0
        
                if trigger == 1 and new_line[1] == 'ATOM': # This extract the lines constituting the .pdb file
                    # This will output the .pdb file with the corresponding format
                    #pdb_output.write(f"{new_line[1]}{new_line[2]:>7}{'':<2}{new_line[3]:<4}{new_line[4]:<8}{new_line[5]:<7}{new_line[6]:<8}{new_line[7]:<8}{new_line[8]:<8}\n")
                    atom_field = "HETATM"
                    pdb_output.write(f"{atom_field:<6}{new_line[2]:>5}{'':<1}{new_line[3]:>4}{'':<1}{new_line[4]:>3}{'':<2}{new_line[5]:>4}{'':<4}{new_line[6]:>8}{new_line[7]:>8}{new_line[8]:>8}\n")
        pdb_output.close()

def store_pdb_file_in_database(database,pdb_file,pdb_blob_obj):
    """
    This function will store a pdb_file as a blob object within a database and under the table Results_pdb_files
    
    ------
    Parameters
    ------
    - database: the full path to the database in which the .pdb files corresponding to the docked poses are to be stored.
    - pdb_file: the name of the .pdb_file to be stored.
    - pdb_blob_obj: the blob object corresponding to the docked .pdb sub_pose
    """
    # Connect to the database
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    
    # Check if the table 'Results_pdb_files' exists, otherwhise, create it.
    try: 
        cur.execute('SELECT * FROM Results_pdb_files')
                
    except:
        # The table does not exists, so create it
        ins = f'CREATE TABLE Results_pdb_files ("sub_pose" TEXT, "blob_pdb_file" BLOB)'
        conn.execute(ins)
        print(colored("Docked poses table has been created","green"))
   
    # Store the sub_pose name and pdb BLOB object
    sub_pose_name = pdb_file.split('/')[-1].replace('.pdb','')
    sql = f"""INSERT INTO Results_pdb_files('sub_pose','blob_pdb_file') VALUES(?,?);"""
    data_tuple = (sub_pose_name,pdb_blob_obj)
    cur.execute(sql,data_tuple)
    conn.commit()
    cur.close()

def store_pdb_file_in_database_by_cluster(database,pdb_file,pdb_blob_obj,ligand_name,cluster_number,counter):
    """
    This function will store a pdb_file as a blob object within a database and under a table named after the cluster to which the results were extracted (e.g. pdb_files_{ligand_name}_cluster_{cluster_number})
    
    ------
    Parameters
    ------
    - database: the full path to the database in which the .pdb files corresponding to the docked poses are to be stored.
    - pdb_file: the name of the .pdb_file to be stored.
    - pdb_blob_obj: the blob object corresponding to the docked .pdb sub_pose
    - ligand_name: the name of the ligand to store
    - cluster_number: the number of the cluster originating the .pdb poses to be stored.
    - counter: if equal to 1, will create the table for storing, if higher, will pass.
    """
    
    # Connect to the database
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    
    # Check if the table 'Results_pdb_files_cluster_#' exists. If so exit to avoid rewriting.
    ligand_name_fixed = ligand_name.replace('-','_') # To avois sql errors with '-'
    
    if counter == 1:
        try: 
            cur.execute(f'SELECT * FROM pdb_files_{ligand_name_fixed}_cluster_{cluster_number}')
            #cur.execute(f'SELECT * FROM pdb_files_AAAAA_cluster_{cluster_number}')
            # Exit if the previous command is successful
            print(colored(f"The table pdb_files_{ligand_name_fixed}_cluster_{cluster_number} already exists","red"))
            exit()
                
        except:
            # The table does not exists, so create it
            ins = f'CREATE TABLE pdb_files_{ligand_name_fixed}_cluster_{cluster_number} ("sub_pose" TEXT, "blob_pdb_file" BLOB)'
            conn.execute(ins)
            print(colored("Docked poses table has been created","green"))
            table_exists = 1
        
    ## Store the sub_pose name and pdb BLOB object
    sub_pose_name = pdb_file.split('/')[-1].replace('.pdb','')
    sql = f"""INSERT INTO pdb_files_{ligand_name_fixed}_cluster_{cluster_number}('sub_pose','blob_pdb_file') VALUES(?,?);"""
    data_tuple = (sub_pose_name,pdb_blob_obj)
    cur.execute(sql,data_tuple)
    conn.commit()
    cur.close()
        
def generate_pdb_from_sdf(output_path):
    """This script will create pdb files starting from sdf files
    
    ------
    Parameters:
    ------
    - output_path: this is the full path to the folder in which one or more .sdf files are stored.
    
    ------
    Returns:
    ------
    The corresponding pdb files are created.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "pdb")
    
    sdf_files_list = glob(f'{output_path}/*.sdf')
    
    for file in sdf_files_list:
        pdb_file = file.replace(".sdf",".pdb")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol,file)
        obConversion.WriteFile(mol,pdb_file)

def generate_pdb_from_sdf_single_file(output_path,sdf_file):
    
    print(sdf_file)

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "pdb")
    
    pdb_file = sdf_file.replace(".sdf",".pdb")
   
   
    mol = openbabel.OBMol()
    input_file= f'{output_path}/{sdf_file}'
    output_file = f'{output_path}/{pdb_file}'
    obConversion.ReadFile(mol,input_file)
    obConversion.WriteFile(mol,output_file)

    # This will return the full path to the generate .pdb file in case it is required for futher processing
    return output_file

def extract_lowest_energy_cluster_pdb(docking_assays_storing_path,assay_nbr,extract_pdb_flag):
    """
    This function will extract the .pdb files corresponding to the lowest energy cluster and store them into a the 'docked_ligands_lowest_energy' folder

    ------
    Parameters:
    ------
    - docking_assays_storing_path: the full path to the folder in which all docking assays are run and stored in folders of type: 'docking_assay_#'.
    - assay_nbr: the number of the docking assays to be processed.
    - extract_pdb_flag: activate wether the extraction of the pdb files should be written to 'docked_ligands_lowest_energy' (0/1)
    """

    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    
    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands_lowest_energy'
   
    
    conn = sqlite3.connect(database)

    sql = f"""SELECT *
              FROM Results
              WHERE pose_rank = 1;"""

    lowest_energy_ligands_df = pd.read_sql_query(sql,conn)
    sorted_df = lowest_energy_ligands_df.sort_values(by=['cluster_size','docking_score'], ascending=[False,True])

    # Store the lowest energy df to a table in the database
    sql_ops.store_df_to_sql(f'{database}',sorted_df,"Results_lowest_energy","replace")

    # If selected, extract the lowest energy .pdb file to the corresponding directory

    if extract_pdb_flag == 1:

        Path(output_path).mkdir(parents=True, exist_ok=True)

        for index, row in lowest_energy_ligands_df.iterrows():
            ligname = row['LigName']
            ## This will extract the lowest energy binding pose from the sdf file, which is the first one: i.e. index 1
            extract_selected_index_from_sdf_file(ligname,database,output_path,1,1)
    
        ## This will generate the .pdb file for all the extracted files
        generate_pdb_from_sdf(output_path)
    
def extract_most_populated_cluster_pdb(docking_assays_storing_path,assay_nbr,extract_most_pop_pdb_flag):
    """"
    This function will extract the .pdb file corresponding to the most populated cluster into a folder named 'docked_ligands_most_populated'
    
    ------
    Parameters:
    ------
    - docking_assays_storing_path: the full path to the folder in which all docking assays are run and stored in folders of type: 'docking_assay_#'.
    - assay_nbr: the number of the docking assays to be processed.
    """

    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands_most_populated'
        
    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(database)

    sql = """SELECT *
             FROM Results;"""

    df = pd.read_sql_query(sql,conn)
    sorted_df = df.sort_values(by=['cluster_size','docking_score'], ascending=[False,True])
    sorted_df_unique = sorted_df.drop_duplicates(subset=['LigName'],keep='first')
    sorted_df_unique_final = sorted_df_unique.sort_values(by=['cluster_size','docking_score'], ascending=[False,True])
    
    # Store the most popuplated clusters df to a table in the database
    sql_ops.store_df_to_sql(f'{database}',sorted_df_unique_final,"Results_most_populated","replace")

    if extract_most_pop_pdb_flag == 1:
        
        Path(output_path).mkdir(parents=True, exist_ok=True)
        
        rank = 1
        for index, row in sorted_df_unique.iterrows():
            ligname = row['LigName']
            index_to_extract = int(row['sub_pose'].split('_')[-1])
            # This will extract the lowest energy binding pose from the sdf file, which is the first one: i.e. index 1
            extract_selected_index_from_sdf_file(ligname,database,output_path,index_to_extract,rank)
            rank+=1
    
        # This will generate the .pdb file for all the extracted files
        generate_pdb_from_sdf(output_path)

def create_receptor_mda_object(receptor_models_registry_path,docking_assays_registry_path,assay_nbr):
    """
    This function will receive the full path to a protein receptor and prepare an MDAnalysis (mda) object for further processing results.

    ------
    Parameters:
    ------
    - receptor_models_registry_path: The full path to the folder in which receptor models registry database is stored.
    - docking_assays_registry_path: The full path to the folder in which docking assays registry database is stored.
    - assay_nbr: the corresponding docking assay_nbr that is being processed.

    ------
    Returns:
    ------
    An MDAnalysis (mda) object corresponding to the receptor.

    """
    # Retrieve the rec_model_id as stored in the docking assays database
    try: 
        conn = sqlite3.connect(f'{docking_assays_registry_path}/docking_registries.db')
        cursor = conn.cursor()
        sql = f"""SELECT rec_model_id
                  FROM docking_registries
                  WHERE assay_id = {assay_nbr};"""
        cursor.execute(sql)
        rec_model_id = cursor.fetchall()[0][0]
        cursor.close()
        conn.close()

    except Exception as error:
        print(error)
        print(colored(f"Error retrieving the rec_model_id corrsponding to assay: {assay_nbr}.","red"))

    # Retrieve the filepath to the .pdb file associated to the docking assay and stored in the receptors registry database
    try: 
        conn = sqlite3.connect(f'{receptor_models_registry_path}/receptors_models_registry.db')
        cursor = conn.cursor()
        sql = f"""SELECT pdb_file
                  FROM receptor_models
                WHERE rec_model_id = {rec_model_id};"""
    
        cursor.execute(sql)
        pdb_file = cursor.fetchall()[0][0]
        cursor.close()
        conn.close()

    except Exception as error:
        print(error)
        print(colored(f"Error retrieving the pdb filepath corresponding to {rec_model_id} from the registry database.","red"))


    # Create the corresponding mda receptor object from the corresponding pdb file
    try:
        u = mda.Universe(pdb_file)
        receptor_plf = plf.Molecule.from_mda(u)
        return receptor_plf
    
    except Exception as error:
        print(error)
        print(colored(f"Error creacting the MDA universe object for receptor model: {rec_model_id}.","red"))

def compute_single_prolif_interaction_pattern(ligand_name,receptor_plf):
    """"
    This function will compute the corresponding prolif interaction pattern for a given ligand:receptor pair.
    
    ------
    Parameters:
    ------
    - ligand_name: the name of the ligand to be used in the interaction fingerpring analysis.
    - receptor_plf: an MDAnalysis receptor object already computed, corresponding to the receptor used in the docking campaign.
    
    ------
    Returns:
    ------
    A dataframe containing a single ligand interaction fingerprint as calculated by ProLIF.
    """
    ligand_plf = plf.sdf_supplier(ligand_name)
    # use default interactions
    list_of_interactions = ['Anionic','CationPi','Cationic','EdgeToFace','FaceToFace','HBAcceptor','HBDonor','Hydrophobic','MetalAcceptor','MetalDonor','PiCation','PiStacking','VdWContact','XBAcceptor','XBDonor']    
    # Compute the list of interactions
    fp = plf.Fingerprint(list_of_interactions)
    fp.run_from_iterable(ligand_plf, receptor_plf)
    df = fp.to_dataframe(drop_empty=False)
    df_proc = process_prolif_raw_dataframe(df)
    df_proc['sub_pose'] = ligand_name.split('/')[-1].replace('.sdf','')
    return df_proc

def compute_all_prolif_interaction_patterns(receptor_models_registry_path,docking_assays_registry_path,assay_nbr):

    receptor_plf = create_receptor_mda_object(receptor_models_registry_path,docking_assays_registry_path,assay_nbr)

    # Connect to the docking assays registries an retrieve to full path to the assay folder
    conn = sqlite3.connect(f'{docking_assays_registry_path}/docking_registries.db')
    cursor = conn.cursor()
    sql = f"""SELECT assay_folder
              FROM docking_registries
              WHERE assay_id = {assay_nbr};"""
    cursor.execute(sql)
    assay_folder = cursor.fetchall()[0][0]
    cursor.close()
    conn.close()

    print(assay_folder)

    # Connect to the processed result dataframe in order to retrieve the ligands to be processed
    conn = sqlite3.connect(f'{assay_folder}/docking_assay_{assay_nbr}.db')
    sql = f"""SELECT sub_pose
              FROM Results;"""

    df = pd.read_sql_query(sql,conn)

    sdf_files_path = f'{assay_folder}/docked_ligands'
    list_of_ligands = glob(f"{sdf_files_path}/*.sdf")
    
    counter = 1
    ligands_prefix_list = [] # This list will be used to append the ligand names to the final Prolif dataframe
    for ligand in list_of_ligands:
        ligand_prefix = ligand.split('/')[-1].replace('.sdf','')
        ligands_prefix_list.append(ligand_prefix)
        try: 
            df = compute_single_prolif_interaction_pattern(ligand,receptor_plf)
            if counter == 1:
                df_complete = df
                counter += 1
            else:
                df_complete = pd.concat([df_complete,df])
    
        except Exception as erro:
            continue

    ## This will make the column 'ligand' the first one and store it into the corresponding results database
    first_column = df_complete.pop('sub_pose')
    df_complete = df_complete.fillna(0) # This will fill missing values with a 0
    df_complete = df_complete.astype(int)
    df_complete_ordered = pd.concat([first_column,df_complete],axis=1)
    sql_ops.store_df_to_sql(f'{assay_folder}/docking_assay_{assay_nbr}.db',df_complete_ordered,"prolif_interactions","replace")

def compute_selected_prolif_interaction_patterns(receptor_models_registry_path,docking_assays_registry_path,docking_raw_data_path,assay_nbr,table_name,sub_pose):
    """
    Given a table name and a subpose name, this function will compute a single prolif interaction pattern, and will store it to a table named 'selected_prolif_interactions_fingerprints'
    """
    receptor_plf = create_receptor_mda_object(receptor_models_registry_path,docking_assays_registry_path,assay_nbr)

    # Connect to the docking assays registries an retrieve to full path to the assay folder
    
    conn = sqlite3.connect(f'{docking_assays_registry_path}/docking_registries.db')
    
    cursor = conn.cursor()
    sql = f"""SELECT assay_folder
              FROM docking_registries
              WHERE assay_id = {assay_nbr};"""
    cursor.execute(sql)
    assay_folder = cursor.fetchall()[0][0]
    cursor.close()
    conn.close()

    print(assay_folder)

    # Connect to the processed result dataframe in order to retrieve the ligands to be processed
    results_database = f'{assay_folder}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(results_database)
    sql = f"""SELECT LigName, sub_pose
              FROM {table_name};"""

    df = pd.read_sql_query(sql,conn)

    # If not exists, create the folder to store extracted sdf file for computing the fingerprint. Otherwise continue
    sdf_output_path = f'{assay_folder}/selected_docked_ligands'
        
    if not os.path.exists(sdf_output_path):
        os.makedirs(sdf_output_path)
    
    # This will extract the selected subpose to the corresponding folder
    index = int(sub_pose.split('_')[-1])
    ligname = sub_pose.split('_')[0]
    extract_selected_index_from_sdf_file(ligname,results_database,sdf_output_path,index,0) # extract the sdf subpose and do not append a rank index

    ligand_name = f'{sdf_output_path}/{sub_pose}.sdf'
    
    # Generate the corresponding .pdb file.
    generate_pdb_from_sdf_single_file(sdf_output_path,f'{sub_pose}.sdf')

    # Try to retrieve the 'prolif_interactions' table within the database to append the fingerprint
    try:
        sql = f"""SELECT *
              FROM prolif_interactions;"""

        df_fingerprints_existing = pd.read_sql_query(sql,conn)
        
        # Check if the sub pose has already been computed
        sub_pose_names = df_fingerprints_existing['sub_pose'].to_list()
        
        if sub_pose in sub_pose_names:
            print(colored(f"The sub pose: {sub_pose} has already been processes and is present. EXITING","red"))
            exit()
        
        fingerprint_df = compute_single_prolif_interaction_pattern(ligand_name,receptor_plf)
        first_column = fingerprint_df.pop('sub_pose')
        df_fingerprint_ordered = pd.concat([first_column,fingerprint_df],axis=1)
        df_complete = pd.concat([df_fingerprints_existing,df_fingerprint_ordered])
        df_complete = df_complete.fillna(0) # This will fill missing values with a 0
        sql_ops.store_df_to_sql(f'{results_database}',df_complete,"prolif_interactions","replace")

    except Exception as error:
        print("Starting the fingerprints table")
        fingerprint_df = compute_single_prolif_interaction_pattern(ligand_name,receptor_plf)
        first_column = fingerprint_df.pop('sub_pose')
        df_fingerprint_ordered = pd.concat([first_column,fingerprint_df],axis=1)
        df_complete = df_fingerprint_ordered.fillna(0) # This will fill missing values with a 0
        sql_ops.store_df_to_sql(f'{results_database}',df_fingerprint_ordered,"prolif_interactions","replace")

    print(colored(f"FINISHED extraction of fingerprint for sub pose: {sub_pose}.","green"))        

def process_prolif_raw_dataframe(df):
    """
    This function will process a dataframe generated by Prolif (i.e. multilevel headers) and will return a processed dataframe in a format compatible with its storing withine de cadd platform.
    
    ------
    Parameters:
    ------
    - df: the input dataframe as generated by the Prolif package.
    
    ------
    Returns:
    df_proc: the dataframe in a format compatible with its storing.
    ------
    """
    df.droplevel(level=0,axis=1).head()
    df_final_grouped = df.copy()
    df_final_grouped.columns = df_final_grouped.columns.map('{0[1]}-{0[2]}'.format)
    
    return df_final_grouped

def perform_ligands_histograms_analysis(docking_assays_storing_path, assay_nbr,cluster_number_threshold):
    """
    This function will save an histogram of the docked_score vs. cluster_size for each ligand processed in a database
    """
    
    database_file = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    output_dir = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/ligs_histograms'
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(database_file)

    sql_inst = "SELECT LigName, docking_score, cluster_size FROM Results;"
    df = pd.read_sql(sql_inst,conn)

    unique_name_list = []
    dock_score_list = []
    clust_size_list = []
    
    number_of_rows = df.shape[0] 
    counter = 1
    
    df.to_csv("/home/fredy/Desktop/ver.csv",index=False)

    for index, row in df.iterrows():
        lig_name = row['LigName']
        
        if lig_name not in unique_name_list: 
            # This part of the code will append the values for the first iteration in the first ligand
            if len(dock_score_list) == 0:
                try:
                    dock_score_list.append(row['docking_score'])
                    clust_size_list.append(row['cluster_size'])        
                    unique_name_list.append(lig_name)
                    #continue
                except Exception as error:
                    print("ERROR")
            
            # Once a new ligand in found, save the plot and reset the lists
            else:
                # This will create the previous ligname variable to store the corresponding data
                previous_lig_name = unique_name_list[-1]
                # This will save the corresponding plot prior to delete del lists
                plt.bar(dock_score_list, clust_size_list, color ='blue', width = 0.1)
                plt.xlabel("docking score")
                plt.ylabel("cluster_size")
                plt.title(f"Clustering analysis for {previous_lig_name}")
                plt.savefig(f"/{output_dir}/{previous_lig_name}.png")
                # Clear the plot for the next iteration
                plt.clf()
                ## Save a csv with the plotted values
                df = pd.DataFrame({'docking_score':dock_score_list,
                                   'cluster_size':clust_size_list})
                df.to_csv(f"/{output_dir}/{previous_lig_name}.csv",index=False,header=False)
                
                # Reset the lists
                dock_score_list = []
                clust_size_list = []
                
                # Store the values and append the LigName
                dock_score_list.append(row['docking_score'])
                clust_size_list.append(row['cluster_size'])        
                unique_name_list.append(lig_name)

        else: 
            # This part of the code will continue to append the values for the same ligand
            try:
                dock_score_list.append(row['docking_score'])
                clust_size_list.append(row['cluster_size'])        
            except Exception as error:
                print(error)

    
        # This will save the plot for the last ligand
        if counter == number_of_rows:
            # This will save the corresponding plot prior to delete del lists
            plt.bar(dock_score_list, clust_size_list, color ='blue', width = 0.1)
            plt.xlabel("docking score")
            plt.ylabel("cluster_size")
            plt.title(f"Clustering analysis for {lig_name}")
            plt.savefig(f"/{output_dir}/{lig_name}.png")
            plt.clf() # Clear the last plot
            previous_lig_name = unique_name_list[-1]
            
            #Save a csv with the plotted values
            df = pd.DataFrame({'docking_score':dock_score_list,
                               'cluster_size':clust_size_list})
            df.to_csv(f"{output_dir}/{previous_lig_name}.csv",index=False,header=False)
        
        # Update the counter 
        counter+=1

    ## Extract the preferred binder from the whole analysis
    extract_preferred_binders(conn,cluster_number_threshold,output_dir)

def extract_preferred_binders(connection,cluster_number_threshold,output_path):
    """
    This function will create a .csv file containing the SMILES and inchi_key of the compounds resulting from a docking campaign in which the number of conformations in the lowest energy cluster is above a certain threshold limit.

    ------
    Parameters
    ------
    - connection: a connection already established to a database corresponding to a docking result processing.
    - cluster_number_threshold: the threshold corresponding to the minumum number of cluster conformatios to be detected in the lowest energy cluster.
    - output_path: The full path were the output .csv file is to be written.

    ------
    Returns
    ------
    A .csv file
    """

    sql = "SELECT * FROM Results WHERE pose_rank = 1"
    df = pd.read_sql(sql,connection)
    df_preferred_binders = df[df["cluster_size"] > cluster_number_threshold]    
    preferred_binders_ligname_list = df_preferred_binders["LigName"].to_list()

    sql_ligands = "SELECT * FROM Ligands"
    df_ligands = pd.read_sql(sql_ligands,connection)

    selected_smiles_list = []
    selected_ligname_list = []
    for index, row in df_ligands.iterrows():
        ligname = row["LigName"]
        if ligname in preferred_binders_ligname_list:
            smiles = str(row["ligand_smile"]).strip()
            selected_ligname = str(row["LigName"]).strip()
            selected_smiles_list.append(smiles)
            selected_ligname_list.append(selected_ligname)

    df_preferred_binders_smiles = pd.DataFrame({"SMILES":selected_smiles_list,
                                                "name":selected_ligname_list})
    df_preferred_binders_smiles.to_csv(f"{output_path}/selected_binders.csv",header=True,index=False)

def perform_docking_assay_analysis(receptor_models_registry_path,docking_assays_registry_path,docking_assays_storing_path,docking_raw_data_path,assay_nbr,max_poses,compute_all_prolif_flag,cluster_number_threshold):
    """
    This function will perform the whole set of analysis for molecular docking assay.
    ------
    Parameters:
    ------
    - receptor_models_registry_path: the full path to the folder in which the project receptor models registries are stored.
    - docking_assays_registry_path: the full path to the folder in which the molecular docking assays registries are stored.
    - docking_assays_storing_path: the full path to the folder in which the molecular docking assays are stored.
    - docking_raw_data_path: the full path to the folder in which the project raw data intended for docking is stored.
    - assay_nbr: the number of the docking assay as stored in the corresponding registries that is to be analyzed.
    - max_poses : the maximum number of poses clusters to retain from the docking file.
    - compute_all_prolif_flag: Flag (0/1) to compute the prolif interactions on all the database ligands poses. Warning: SPACE and TIME.
    - extract_pdb_poses_from_cluster: this flag controls if a certain cluster number needs to be processed and store all the .pdb file included in it as a separate .pdb files in a table named: 'Results_pdb_files_cluster_#'
    
    ------
    Returns:
    ------
    The following information:
    - A folder named 'docked_ligands' stored within the general docking assay folder. Contains .sdf and .pdb files of all the molecules docked within the analyzed assay, including all the stored docking poses.
    - A database named 'docking_assay_#.db' containing all the ringtail structure and information associated to the docking assay analyzed.
    - A table named 'prolif_interactions' stored within the 'docking_assay_#.db' containing all the interaction fingerprints as calculated by the package ProLIF.
    """
    
    try: 
        ## Step_1: Create the corresponding ringtail database:
        create_db_from_dlgs(docking_assays_storing_path, assay_nbr, max_poses)
        print(colored("Finished creating database from .dlg files","green"))
            
        ## Perform the ploting and storing of the corresponding histograms
        perform_ligands_histograms_analysis(docking_assays_storing_path, assay_nbr, cluster_number_threshold)
            
        ## Extract one .pdb file for each cluster
        extract_1_per_cluster_from_dlg(docking_assays_storing_path,assay_nbr)
        
        ## Extract the lowest energy conformation for each ligand
        extract_lowest_energy_conformation_from_dlg(docking_assays_storing_path,assay_nbr)
        
        ## Extract the most populated cluster of conformations
        extract_most_populated_cluster_conformation_from_dlg(docking_assays_storing_path,assay_nbr)
        
        if compute_all_prolif_flag == 1:
            # Step_5: Compute and store the corresponding ProLIF interactions for the whole assay:
            compute_all_prolif_interaction_patterns(receptor_models_registry_path,docking_assays_registry_path,docking_raw_data_path,assay_nbr)
    
    
        print(colored(f"Docking analysis of assay: {assay_nbr} finished SUCCESFULLY.","green"))
    
    except Exception as error:
        print(error)
        print(colored(f"ERROR in docking analysis of assay: {assay_nbr}.","red"))

def process_docking_isolated_folder(dest_folder):
    """
    This function will apply all the docking processes to a custom folder in which the docking results are stored in the directories:

        $PATH/dlgs
        $PATH/pdbqt_files
        $PATH/receptor

    """

    create_db_from_dlgs_custom_path(dest_folder)
    extract_1_per_cluster_from_dlg_custom_path(dest_folder)

def extract_lowest_energy_pdb_file_range(docking_assays_storing_path,assay_nbr,start,stop):
    """
    This function will extract a pdb files range from the table 'Results_lowest_energy'
    """

    # Create a connection to the already created database
    database = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    
    # Configure the output path to store corresponding .sdf files    
    output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/results_lowest_energy_{start}_{stop}'
    
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    else: 
        print(colored("The directory containing the docked ligand files already exists. Stoping execution","red"))
        exit()

    conn = sqlite3.connect(database)

    sql = f"""SELECT LigName, sub_pose
              FROM Results_lowest_energy;"""

    ligands_name_df = pd.read_sql_query(sql,conn)
    ligands_name_df_range = ligands_name_df[start-1:stop+1] # the subtraction and addition is to compatibilize with the 0-based indexing

    for index, row in ligands_name_df_range.iterrows():
        ligname = row['LigName']
        index = int(row['sub_pose'].split('_')[1])
        print(index)
        extract_selected_index_from_sdf_file(ligname,database,output_path,index,0)
    
    # This function will generate the corresponding pdb files from the generated sdf docked poses
    generate_pdb_from_sdf(output_path)

def plot_grid_favored_binders(docking_assays_storing_path,docking_misc_data_path,assay_nbr,cluster_thr,docking_score_thr):
    """
    Given a table of processes docking results, this function will plot the favored binders (i.e. cpds with a higher number of cluster members and low docked energy) to a folder in the '/misc' subfolders.
    """
    db_file = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(db_file)
    
    sql = """SELECT ligand_smile, Results_lowest_energy.LigName, sub_pose, docking_score, cluster_size
             FROM Results_lowest_energy
             INNER JOIN Ligands ON Results_lowest_energy.LigName = Ligands.LigName;"""

    df = pd.read_sql_query(sql,conn)
    df_favored = df[df['cluster_size'] > cluster_thr]
    df_favored2 = df_favored[df_favored['docking_score'] < docking_score_thr]
    
    # Check if some compounds passed the threshold
    
    if len(df_favored2) == 0:
        print(colored("No compounds passed the indicated cluster and docking score values threshold value. EXITING","red"))
        print(colored(f"Passing cluster threshold: {len(df_favored)}. EXITING","red"))
        print(colored(f"Passing docking score threshold: {len(df_favored2)}. EXITING","red"))
        exit()
    
    else:
        print(colored(f"The number of cpds passing clustersize threshold is: {len(df_favored)}","green"))
        print(colored(f"The number of cpds passing docking_score threshold is: {len(df_favored2)}","green"))
        action = input(colored("Continue with this threshold value? (y/n): ", "green"))
        
        if action == 'y':
            sorted_df = df_favored2.sort_values(['cluster_size','docking_score'], ascending=[False,True])
            # generate the depiction grids
            output_path = f'{docking_misc_data_path}/docking_assay_{assay_nbr}'
            depict_res.generate_depiction_grid_with_results_labels(output_path,sorted_df,25) 
            # store the favored ligands dataframe to the database
            sql_ops.store_df_to_sql(f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docking_assay_{assay_nbr}.db',sorted_df,"favored_ligands","replace")
            
            # Create the .pdb files for the ligands passing the above mentiones threshold values
            pdb_output_path = f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/preferred_ligands'
            
            if os.path.isdir(pdb_output_path):
                print(colored(f"The destination folder '{pdb_output_path}' for preferred ligands already exists. EXITING","red"))
                exit()
            else:
                ## Store the favored ligands .pdb file within a specific folder
            
                if os.path.isdir(pdb_output_path):
                    print(colored(f"The destination folder '{pdb_output_path}' for preferred ligands already exists. EXITING","red"))
                    exit()
                else:
                    # Create the destination directory
                    os.makedirs(pdb_output_path)                   
                    print(colored(f"SUCESSFULLY created the output folder: {pdb_output_path}","green"))
                    # Exctract the .sdf and .pdb files for the selected df
                    create_pdb_from_df(db_file,sorted_df,pdb_output_path)
            
            
            #    ligname_list = sorted_df["sub_pose"].to_list()
            #    counter=1
            #    for ligname in ligname_list:
            #        src=f'{docking_assays_storing_path}/docking_assay_{assay_nbr}/docked_ligands/{ligname}.pdb'
            #        dest= f'{preferred_ligands_path}/rank_{counter}_{ligname}.pdb'
            #        shutil.copyfile(src,dest)
            #        counter+=1

        else:
            print(colored("EXIT for new threshold value input.","red"))
            exit()

def create_pdb_from_df(db_file,df,output_path):
    """
    This function will receive a df of selected binders and will output the corresponding .pdb files to a folder named 'selected_binders'
    """
        
    rank = 1
    for index, row in df.iterrows():
        ligname = row['LigName']
        index = int(row['sub_pose'].split('_')[-1])
        extract_selected_index_from_sdf_file(ligname,db_file,output_path,index,rank)
        rank+=1
    
    # This function will generate the corresponding pdb files from the generated sdf docked poses
    generate_pdb_from_sdf(output_path)

def predict_binder_class(ml_models_path,receptor_models_registry_path,docking_assays_registry_path,docking_raw_data_path,assay_nbr,table_name,classifier_model,pickle_file,extract_pdb,pdb_files_chunk_size):
    """
    This function use an already saved classification model (in the pickle format) and will classiffy a docked pose as a binder (1) or non-binder (0) class
    """

    receptor_plf = create_receptor_mda_object(receptor_models_registry_path,docking_assays_registry_path,assay_nbr)

    conn = sqlite3.connect(f'{docking_assays_registry_path}/docking_registries.db')
    
    cursor = conn.cursor()
    sql = f"""SELECT assay_folder
              FROM docking_registries
              WHERE assay_id = {assay_nbr};"""
    cursor.execute(sql)
    assay_folder = cursor.fetchall()[0][0]
    cursor.close()
    conn.close()

    # Load the DUMMY fingreprint record to normalize columns to whole receptor
    ml_database = f'{ml_models_path}/ml_model_devel/interaction_fingerprints_data.db'
    conn = sqlite3.connect(ml_database)
    cursor = conn.cursor()
    sql = """SELECT fingerprint_dict
             FROM positive_binders
             WHERE sub_pose == 'DUMMY';"""
    cursor.execute(sql)
    dummy_dict_str = cursor.fetchall()[0][0]
    dummy_dict = eval(dummy_dict_str)
    dummy_df = pd.DataFrame(dummy_dict, index=['i',])
    cursor.close()
    conn.close()

    # If not exists, create the folder to store extracted sdf file for computing the fingerprint. Otherwise continue
    sdf_output_path = f'{assay_folder}/tmp'
        
    if not os.path.exists(sdf_output_path):
        os.makedirs(sdf_output_path)
    
    # Connect to the processed result dataframe in order to retrieve the ligands to be processed
    results_database = f'{assay_folder}/docking_assay_{assay_nbr}.db'
    conn = sqlite3.connect(results_database)
    sql = f"""SELECT *
              FROM {table_name};"""

    df = pd.read_sql_query(sql,conn)
    conn.close()


    # load model
    with open(f"{ml_models_path}/misc/pickle_model_files/{pickle_file}",'rb') as f:
        loaded_classifier = pickle.load(f)

    # Create a list to store the class_classification code
    class_classification_list = []

    for index, row in df.iterrows():
        sub_pose = row['sub_pose']
        # This will extract the selected subpose to the corresponding folder
        index = int(sub_pose.split('_')[-1])
        ligname = sub_pose.split('_')[0]
        extract_selected_index_from_sdf_file(ligname,results_database,sdf_output_path,index,0) # extract the sdf subpose and do not append a rank index
        ligand_name = f'{sdf_output_path}/{sub_pose}.sdf'
        # Compute the fingerprint
        fingerprint_df = compute_single_prolif_interaction_pattern(ligand_name,receptor_plf)
        fingerprint_df.replace({False: 0.0, True: 1.0}, inplace=True)
        ## Concatenate the dummy dict with the corresponding calculated df
        df_concat = pd.concat(([dummy_df,fingerprint_df]))
        ## Delete the first row of the dataframe that corresponds to the the DUMMY record
        df_concat.drop(index=df_concat.index[0], axis=0, inplace=True)
        ## Fill missing values in the row
        df_concat.fillna(0,inplace=True)
        df_concat = df_concat.drop('sub_pose', axis=1)
        df_concat.to_csv("/home/fredy/Desktop/fingerprint.csv",header=None, index=False)     
        X_concat = df_concat.to_numpy()[0].reshape(1, -1)
        
        binder_class = loaded_classifier.predict(X_concat)[0]
        # Append a to a classification list in order to store in the database
        class_classification_list.append(binder_class)
        
        if binder_class == 1:
            print(colored(f"Binder class: {binder_class}","green"))
        else:
            print(colored(f"Binder class: {binder_class}","red"))

    # Add the class_classification to the dataframe for further storage
    df["binder_class"] = class_classification_list

    df_positives = df[df["binder_class"] == 1]
    df_negatives = df[df["binder_class"] == 0]

    # Replace the original table with the one containg the class_classification

    sql_ops.store_df_to_sql(results_database,df_positives,f'{table_name}_predicted_positive',"replace")
    sql_ops.store_df_to_sql(results_database,df_negatives,f'{table_name}_predicted_negative',"replace")

    # Flag to indicate wether the pdb files needs to be extracted
    if extract_pdb == 1:
        # Extract the corresponding pdb files for the possitive and negative tables
        extract_pdb_for_whole_table(results_database,f'{table_name}_predicted_positive',f'{assay_folder}/{table_name}_predicted_positive_pdb_files',pdb_files_chunk_size)
        extract_pdb_for_whole_table(results_database,f'{table_name}_predicted_negative',f'{assay_folder}/{table_name}_predicted_negative_pdb_files',pdb_files_chunk_size)
    
        # Delete the folder containing the corresponding sdf files
        shutil.rmtree(sdf_output_path)

    else:
        # Delete the folder containing the corresponding sdf files
        shutil.rmtree(sdf_output_path)

def extract_pdb_for_whole_table(db_name,table_name,main_output_path,pdb_files_chunk_size):
    """
    This function will receive a database and a table name, and will extract the whole set of pdb into a subfolder named 'table_name'_pdb_files
    """

    conn = sqlite3.connect(db_name)
    sql = f"""SELECT LigName, sub_pose
             FROM '{table_name}';"""

    df = pd.read_sql_query(sql,conn)

    # Create the folder to output the corresponding files
    output_path = f'{main_output_path}/part_1'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    rank = 0 # Do not append rank to the corresponding files.
    counter = 1
    part = 2
    for index, row in df.iterrows():
        # If the counter == the chunksize save to a new folder
        if counter % pdb_files_chunk_size == 0:
            # Generate the .pdb files for the folder prior to updating the part number
            generate_pdb_from_sdf(output_path)
            # Delete the .sdf files in the folder
            delete_sdf_files_in_folder(output_path)
            # Update the output path to the next part
            output_path = f'{main_output_path}/part_{part}'
            part+=1
            if not os.path.exists(output_path):
                os.makedirs(output_path)
        # Continue with the processing
        ligname = row['LigName']
        index = int(row['sub_pose'].split('_')[-1])
        extract_selected_index_from_sdf_file(ligname,db_name,output_path,index,rank)
        counter+=1
    
    # This function will generate the corresponding pdb files from the generated sdf docked poses
    generate_pdb_from_sdf(output_path)
    # Delete the .sdf files in the folder
    delete_sdf_files_in_folder(output_path)

def delete_sdf_files_in_folder(folder):
    list_of_sdf_files = glob(f'{folder}/*.sdf')
    for file in list_of_sdf_files:
        os.remove(file)


#################

if __name__ == '__main__':

    print("This is a test")

    #create_db_from_dlgs("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/docking_assays", 1)

    #print_ringtail_default_ops()

    #split_sdf_file("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/docking_assays/docking_assay_1/sdf_parts","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/docking_assays/docking_assay_1/docked_ligands/AXGVWFZIMFSNCN-UHFFFAOYSA-N.sdf")

    #extract_all_sdf_molecules("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/docking_assays", 1)

    #create_receptor_mda_object("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/receptor_models_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/docking_assays_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/raw_data,",1)

    #compute_all_prolif_interaction_patterns("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/receptor_models_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/docking_assays_registers","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/docking/raw_data",1)