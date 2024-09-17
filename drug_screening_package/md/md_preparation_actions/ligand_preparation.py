import sys
import subprocess
import sqlite3
import tarfile 
import os

MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import md.perform_md.md_configuration as md_conf
from docking.docking_results_processing import process_docking_assay as dock_proc

def retrieve_docked_subpose_pdb_file(ligand_table_name,ligand_db,ligand_subpose_name,output_folder):
    """
    This function will retrieve the corresponding sub_pose from the database, and will write it to a destination path
    
    ------
    Parameters:
    ------
    - ligand_table_name: the name of the table were the .pdb of the ligand to be retrieved is stored
    - ligand_db: the database as processed with ringtail in which the docked pose to be extracted is present
    - ligand_subpose: the name of the subpose, which is searched in the 'Results' table unde the column 'sub_pose'.
    - output_path: the full path to the place in which the retrieved .pdb file is to be stored.
     
    """
    conn = sqlite3.connect(ligand_db)
    
    cursor = conn.cursor()
    
    sql = f"""SELECT * 
             FROM {ligand_table_name}
             WHERE sub_pose = '{ligand_subpose_name}'
          """
    cursor.execute(sql)
    query_output = cursor.fetchall()
    
    for row in query_output:
        data = row[1]
    
    dest_file = f'{output_folder}/{ligand_subpose_name}.tar.gz'
    with open(dest_file, 'wb') as file:
        file.write(data)
    file.close()
    print(f"File {dest_file} stored SUCCESSFULLY")
    
    # This will extract the retrieved tar.gz file and delete it afterwards
    file = tarfile.open(f'{output_folder}/{ligand_subpose_name}.tar.gz')
    file.extract(f'{ligand_subpose_name}.pdb',f'{output_folder}')
    file.close()
    # Delete the compressed file restored from the database
    os.remove(f'{output_folder}/{ligand_subpose_name}.tar.gz')

    # return the ligand .pdb filename for further processing
    return f'{output_folder}/{ligand_subpose_name}.pdb'

def retrieve_reference_pdb_file(md_assay_folder,source_ligands_db,source_ligands_pdbqt_table,ligand_inchi_key):
    
    # This will retrieve the corresponding database containing ligands
    conn = sqlite3.connect(source_ligands_db)
    
    cursor = conn.cursor()
    
    sql = f"""SELECT * 
             FROM {source_ligands_pdbqt_table}
             WHERE inchi_key = '{ligand_inchi_key}'
          """
    cursor.execute(sql)
    query_output = cursor.fetchall()
    
    for row in query_output:
        stored_pdb_filename = row[4]
        stored_pdb_blob = row[5]
        
    dest_file = f'{md_assay_folder}/{stored_pdb_filename}'
    with open(dest_file, 'wb') as file:
        file.write(stored_pdb_blob)
    file.close()
    print(f"File {dest_file} stored SUCCESSFULLY")
    
     # This will extract the retrieved tar.gz file and delete it afterwards
    reference_pdb_file = stored_pdb_filename.replace('.tar.gz','.pdb')
    file = tarfile.open(f'{md_assay_folder}/{stored_pdb_filename}')
    file.extract(f'{reference_pdb_file}',f'{md_assay_folder}')
    file.close()
    os.remove(f'{md_assay_folder}/{stored_pdb_filename}')

    return f'{md_assay_folder}/{reference_pdb_file}'

def restore_all_reference_pdb_files_to_path(source_db,ligands_pdbqt_table,output_path):
    """
    This function will retrieve all the .pdb files containing all Hs atoms and correctly numbered by RDKit as stored in a {source_ligands_pdbqt_table} table, and will output the .pdb files the {output_path}. These .pdb files are consistant with the .pdb poses extracted by the docking analysis tools, and can be further used to prepare MD parameters files.
    
    """

    # This will retrieve the corresponding database containing ligands
    conn = sqlite3.connect(source_db)
    
    cursor = conn.cursor()
    
    sql = f"SELECT * FROM {ligands_pdbqt_table}"
    
    cursor.execute(sql)
    query_output = cursor.fetchall()
    
    
    for row in query_output:
        stored_pdb_filename = row[4]
        print(stored_pdb_filename)
        stored_pdb_blob = row[5]
    
        dest_file = f'{output_path}/{stored_pdb_filename}'
        with open(dest_file, 'wb') as file:
            file.write(stored_pdb_blob)
        file.close()
        print(f"File {dest_file} stored SUCCESSFULLY")
    
        # This will extract the retrieved tar.gz file and delete it afterwards
        reference_pdb_file = stored_pdb_filename.replace('.tar.gz','.pdb')
        file = tarfile.open(f'{output_path}/{stored_pdb_filename}')
        file.extract(f'{reference_pdb_file}',f'{output_path}')
        file.close()
        os.remove(f'{output_path}/{stored_pdb_filename}')

def process_pdb_with_antechamber(file,output_path,md_cond_dict):
    """
    This function will receive a pdb file and process it with antechamber in order to generate the corresponding .prepin and .frcmod files.

    ------
    Parameters
    ------
    -file: the full path to the .pdb file to be processed.
    - output_path: the full path to the folder in which the generated files are to be written.
    """

    file_prefix = file.split('/')[-1].replace('.pdb','')

    # Remove the unwanted keywords in the starting .pdb fileto avoid antechamber errors
    unwanted_keywords = ['CONECT','MASTER','END','COMPND','AUTHOR']
    
    with open(f'{output_path}/temp.pdb','w') as temp_pdb_file:
        for line in open(file):
            new_line = line.rstrip()
            if new_line.split()[0] in unwanted_keywords:
                continue
            else:
                temp_pdb_file.write(f'{new_line}\n')
    
    # Set some conditions to run antechamber using the cond_dictionary provided

    charge = md_cond_dict['lig_antechamber_parms']['charge_system']
    clean_files = md_cond_dict['lig_antechamber_parms']['clean_files']
    residue_name = md_cond_dict['lig_antechamber_parms']['residue_name']
            
    # Prepare the .prepin file
    bash_command = f'antechamber -i {output_path}/temp.pdb -fi pdb -o {output_path}/{file_prefix}.prepin -fo prepi -c {charge} -pf {clean_files} -rn {residue_name}'
    process = subprocess.Popen(bash_command.split(),stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Prepare the .frcmod file
    bash_command = f'parmchk2 -i {output_path}/{file_prefix}.prepin -f prepi -o {output_path}/{file_prefix}.frcmod'
    process = subprocess.Popen(bash_command.split(),stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Return the full path to the .prepin file for further use
    return f'{output_path}/{file_prefix}.prepin'

def prepare_pdb_from_prepin(file,md_assay_folder):
    """
    This function will use tleap in order to generate a .pdb file consistant with the .prepin and .frcmod files previously generated.

    ------
    Parameters:
    ------
    - file: the full path to the .prepin file to be processed

    """

    file_prefix = file.split('/')[-1].replace('.prepin','')

    # This will write a tleap input file to process the corresponding file
    with open(f'{md_assay_folder}/lig_tleap_{file_prefix}.in','w') as ligand_tleap_file:
        ligand_tleap_file.write("source leaprc.gaff2 \n") 
        ligand_tleap_file.write(f"loadamberprep {file} \n") 
        ligand_tleap_file.write(f"loadamberparams {md_assay_folder}/{file_prefix}.frcmod \n")
        ligand_tleap_file.write(f"savepdb MOL {md_assay_folder}/{file_prefix}_tleap.pdb \n")
        ligand_tleap_file.write(f"quit \n")

    # execute the tleap file prepared above
    bash_command = f'tleap -f {md_assay_folder}/lig_tleap_{file_prefix}.in'
    process = subprocess.Popen(bash_command.split(),stdout=subprocess.PIPE)
    output, error = process.communicate()

def retrive_and_prepare_docked_pose(ligand_db,ligand_subpose_name,source_ligands_db,ligand_table_name,source_ligands_pdbqt_table,md_assay_folder,md_cond_dict):
    """
    Given a docked sub_pose that is stored within a database as a docking result, this function will extract the corresponding .pdb file and prepare the ligand for molecular dynamics simulations.

    ------
    Parameters:
    ------
    - ligand_db: the database as processed with ringtail in which the docked pose to be extracted is present.
    - ligand_subpose: the name of the subpose, which is searched in the 'Results' table unde the column 'sub_pose'. 
    - source_ligands_db: the full path of the database containing the .pdbqt and associated .pdb file consistant with docked poses.
    - ligand_table_name: the name of the table in which the .pdb file to extract is stored
    - source_ligands_pdbqt_table: the name of the table within 'source_ligands_db' that contains the original .pdbqt and .pdb files.
    - md_assay_folder: the full path to the folder in which the MD assay information and files is to be stored.
    - md_cond_dict: the set of conditions (as a dictionary) to be applied to the MD simulations

    ------
    Returns:
    ------

    """
 
    # Retrieve the .pdb file of the docking sub_pose
    ligand_docked_pdb_file = retrieve_docked_subpose_pdb_file(ligand_table_name,ligand_db,ligand_subpose_name,md_assay_folder)
 
    # Retrieve the reference .pdb file used to prepare the corresponding .pdbqt prior to docking assays
    ligand_inchi_key = ligand_subpose_name.split('_')[0]
    reference_pdb_file = retrieve_reference_pdb_file(md_assay_folder,source_ligands_db,source_ligands_pdbqt_table,ligand_inchi_key)
    
    ## Process the corresponding .pdb with antechamber
    prepin_file = process_pdb_with_antechamber(reference_pdb_file,md_assay_folder,md_cond_dict)

    return ligand_docked_pdb_file