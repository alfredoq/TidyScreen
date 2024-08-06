import sys
import os
from termcolor import colored
from ringtail import RingtailCore
import sqlite3
import pandas as pd

MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import md.md_registers_management.md_assays_registers as md_regs
from docking.docking_results_processing import process_docking_assay as dock_proc

def create_md_assay_folder(md_assays_registry_path,md_assays_storing_path,ligand_db,ligand_table,ligand_subpose_name,rec_model_id,md_params_id):

    """
    This function will create a folder to store all the information required to run an md assay starting from a selected binding mode as obtained in a molecular docking assay.

    The function will query the last assay_id registered, and create the next consecutive folder.

    ------
    Parameters:
    ------

    ------
    Returns:
    ------
    A folder in which the information for a specific docking assay is stored.
    """

    last_value = md_regs.get_last_md_assay_id(md_assays_registry_path)

    # This will check is the database of registries exists, if not, will create it
    if last_value == None:
        md_regs.create_md_registry_table(md_assays_registry_path)
        print(colored("The docking assays registry tables did not existed. It was created","red"))
        last_value = 0

    md_assay_folder = f'{md_assays_storing_path}/md_assay_{last_value+1}'
    
    try: 
        ## This will append the md assay registry to the corresponding assays registries table
        description = input(colored("Provide a brief description of the MD assay: ","green"))
        md_regs.append_md_registry(md_assays_registry_path,ligand_db,ligand_table,ligand_subpose_name,rec_model_id,0,md_assay_folder,description,md_params_id)
        ## This will create the subfolders structure for the MD assay
        os.mkdir(md_assay_folder)
        print(colored(f"MD assay folder: md_assay_{last_value+1} CREATED","green"))
        
        return md_assay_folder
        
    except Exception as error:
        print(colored("Error creating the MD assay folder","red"))
        print(error)

def retrieve_assay_relevant_variables(docking_assays_registry_path,assay_id,receptor_models_registry_path):

    """
    Given a docking assay ID, this function will return relevant variables specific to the docking assay performed.
    """
    docking_assays_database = f'{docking_assays_registry_path}/docking_registries.db'
    receptor_models_database = f'{receptor_models_registry_path}/receptors_models_registry.db'

    try: 
        # Open the connection to the database of docking assays
        conn = sqlite3.connect(docking_assays_database)
        sql = f"""SELECT * FROM docking_registries
                WHERE assay_id = {assay_id}"""
        df = pd.read_sql_query(sql,conn)

        ligand_db = df['assay_folder'][0]+'/'+df['assay_folder'][0].split('/')[-1]+'.db'
        rec_model_id = int(df['rec_model_id'][0])
        source_ligands_db = df['db_name'][0]
        source_ligands_pdbqt_table = df['ligs_table_name'][0]

        # Close the connection to the database of docking assays
        conn.close()

        # Open the connection to the database of receptor models
        conn = sqlite3.connect(receptor_models_database)
        sql = f"""SELECT * FROM receptor_models
                WHERE rec_model_id = {rec_model_id}"""
        
        df = pd.read_sql_query(sql,conn)
        receptor_pdb_file = df['pdb_file'][0]

        return ligand_db,source_ligands_db,source_ligands_pdbqt_table, rec_model_id, receptor_pdb_file

    except Exception as error:
        print(error)

def prepare_md_inputs(md_assay_folder,md_cond_dict):
    """
    This function will prepare the MD inputs (i.e. min1.in, ...) based on the assay parameters stored in corresponding database, and that are provided as a dictionary.
    
    ------
    Parameters:
    ------
    - md_assay_folder: the full path to the folder in which the md_assay if being stored.
    - md_cond_dict: the set of md parameters retrieved from the corresponding database.
    """

    # Prepare min1.in
    min1_dict = md_cond_dict["min1_params"]
    
    with open(f'{md_assay_folder}/min1.in','w') as min1:
        min1.write("Initial minimization of the system with restraints \n")
        min1.write("&cntrl \n")
        for param, value in min1_dict.items():
            # Check if param is 'restraintmask' to add a ' character not stored in the dict
            if param == 'restraintmask':
                min1.write(f" {param} = '{value}', \n")    
            else:
                min1.write(f' {param} = {value}, \n')
        min1.write("&end \n")
        min1.write("END \n")
    min1.close()
        
    # Prepare min2.in
    min2_dict = md_cond_dict["min2_params"]
    
    with open(f'{md_assay_folder}/min2.in','w') as min2:
        min2.write("Second minimization of the system with no restraints \n")
        min2.write("&cntrl \n")
        for param, value in min2_dict.items():
            min2.write(f' {param} = {value}, \n')
        min2.write("&end \n")
        min2.write("END \n")
    min2.close()
        
    # Prepare heat1.in
    heat1_dict = md_cond_dict["heat1_params"]
    
    with open(f'{md_assay_folder}/heat1.in','w') as heat1:
        heat1.write("First heating stage of the solvent, with restraints in the solute \n")
        heat1.write("&cntrl \n")
        # Write the params corresponding to the md section
        for param, value in heat1_dict["mds"].items():
            if param == 'restraintmask':
                heat1.write(f" {param} = '{value}', \n")    
            else:
                heat1.write(f' {param} = {value}, \n')
        heat1.write("&end \n")
        heat1.write("&wt \n")
        # Write the params corresponding to the temp change section
        for param, value in heat1_dict["temp_change"].items():
            if param == 'TYPE':
                heat1.write(f" {param} = '{value}', \n")
            else:
                heat1.write(f" {param} = {value}, \n")
        # End the input file
        heat1.write("/ \n")
        heat1.write("&wt TYPE = 'END' \n")
        heat1.write("/ \n")
        heat1.write("END \n")
    heat1.close()
    
    # Prepare heat2.in
    heat2_dict = md_cond_dict["heat2_params"]
    
    with open(f'{md_assay_folder}/heat2.in','w') as heat2:
        heat2.write("Second heating stage of the solvent and protein, with restraints in the ligand \n")
        heat2.write("&cntrl \n")
        # Write the params corresponding to the md section
        for param, value in heat2_dict["mds"].items():
            if param == 'restraintmask':
                heat2.write(f" {param} = '{value}', \n")    
            else:
                heat2.write(f' {param} = {value}, \n')
        heat2.write("&end \n")
        heat2.write("END \n")
    heat2.close()
    
    # Prepare equi.in
    equi_dict = md_cond_dict["equi_params"]
    
    with open(f'{md_assay_folder}/equi.in','w') as equi:
        equi.write("Equilibration stage in which only the protein backbone is restrained \n")
        equi.write("&cntrl \n")
        # Write the params corresponding to the md section
        for param, value in equi_dict["mds"].items():
            if param == 'restraintmask':
                equi.write(f" {param} = '{value}', \n")    
            else:
                equi.write(f' {param} = {value}, \n')
        equi.write("&end \n")
        equi.write("END \n")
    equi.close()
    
    # Prepare prod.in
    prod_dict = md_cond_dict["prod_params"]
    
    with open(f'{md_assay_folder}/prod.in','w') as prod:
        prod.write("Production stage in which only the protein backbone is restrained \n")
        prod.write("&cntrl \n")
        # Write the params corresponding to the md section
        for param, value in prod_dict["mds"].items():
            if param == 'restraintmask':
                prod.write(f" {param} = '{value}', \n")    
            else:
                prod.write(f' {param} = {value}, \n')
        prod.write("&end \n")
        prod.write("END \n")
    prod.close()
    
    # Write the local MD execution file
    with open(f'{md_assay_folder}/md_local_exec.sh','w') as execution:
        execution.write('pmemd.cuda -O -i min1.in -o min1.out -p complex.prmtop -c complex.inpcrd -r min1.crd -ref complex.inpcrd \n')
        execution.write('pmemd.cuda -O -i min2.in -o min2.out -p complex.prmtop -c min1.crd -r min2.crd  \n')
        execution.write('pmemd.cuda -O -i heat1.in -o heat1.out -p complex.prmtop -c min2.crd -r heat1.crd -ref min2.crd -x heat1.nc \n')
        execution.write('pmemd.cuda -O -i heat2.in -o heat2.out -p complex.prmtop -c heat1.crd -r heat2.crd -ref heat1.crd -x heat2.nc \n')
        execution.write('pmemd.cuda -O -i equi.in -o equi.out -p complex.prmtop -c heat2.crd -r equi.crd -ref heat2.crd -x equi.nc\n')
        execution.write('pmemd.cuda -O -i prod.in -o prod.out -p complex.prmtop -c equi.crd -r prod.crd -ref equi.crd -x prod.nc \n')
        
