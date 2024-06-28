import sqlite3
from termcolor import colored
import xtarfile as tarfile
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from pandarallel import pandarallel # import pandarallel
from tqdm import tqdm # import tqdm
import tarfile
import os

import sys
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import chemspace.general_procedures.smiles_processing as smi_proc
import docking.receptor_management.blob_operations as blob_ops

def create_ligands_blob_table(db_name,table_name):
    """
    This function will create a table specifically designed to store ligands .pdbqt files for docking campaigns.

    The resulting table name derives from the name of the original table containing the SMILES strings, with '_pdbqt_ligands'.

    The existence of the table is checked when starting the execution. If the table to be created exists, the script will stop.

    The table contains the following string:
        - SMILES: (TEXT) the string corresponding to the processed compound.
        - inchi_key: (TEXT) a the corresponding identifier associated to the ligand file.
        - stored_file: (TEXT) the name of the file to be stored.
        - blob_item: (BLOB) the ligand file stored as a BLOB item.

    ------
    Parameters:
    ------
    - db_name: the fullpath to the database file in which the table to store the ligands .pdbqt files is to be created.
    - table_name: the table containing the SMILES string to be processed to generate the ligands .pdbqt files.

    ------
    Returns:
    ------
    A table containing all the info including the stored .pdbqt files.
    """    
    conn = sqlite3.connect(db_name)

    try:
        sql = f"""CREATE TABLE '{table_name}' 
                  (SMILES TEXT NOT NULL, 
                   inchi_key TEXT NOT NULL, 
                   stored_file TEXT NOT NULL, 
                   blob_item BLOB NOT NULL);"""      
        
        conn.execute(sql)
        print(colored("SUCCESFULLY created ligands BLOB objects table","green"))
    
    except:
        print(colored("The table containing ligands already exists. Revise.","red"))    
        sys.exit()

def create_blob_object(filename):
    """
    This function will create a binary object from a filename provided that is aimed to be stored as a BLOB object within the database.
    ------
    Parameters:
    ------
    - filename: the filename in the corresponding format to be stored in the table containing BLOB objects.

    ------
    Returns:
    ------
    A blob object containing the ligand .pdbqt file in a compressed (.tar.gz) format.
    """
    with open(filename, 'rb') as file:
        blobData = file.read()
    file.close()
    
    return blobData

def compute_inchi_key(smiles):
    """
    This function will be compute an inchi_key from a SMILES string input.

    ------
    Parameters:
    ------
    - smiles: the string input to compute the inchi_key.

    ------
    Returns
    ------
    An inchi_key value for further use.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        return inchi_key
    except:
        pass

def compute_pdqbt_file(smiles):
    dictio = {
    "ATOM_PARAMS": {
    "defaults+4H_triazole": [
    {"smarts": "[#1]", "atype": "H",},
    {"smarts": "[#1][#7,#8,#9,#15,#16]","atype": "HD"},
    {"smarts": "[#5]", "atype": "B"},
    {"smarts": "[C]", "atype": "C"},
    {"smarts": "[c]", "atype": "A"},
    {"smarts": "[#7]", "atype": "NA"},
    {"smarts": "[#8]", "atype": "OA"},
    {"smarts": "[#9]", "atype": "F"},
    {"smarts": "[#12]", "atype": "Mg"},
    {"smarts": "[#14]", "atype": "Si"},
    {"smarts": "[#15]", "atype": "P"},
    {"smarts": "[#16]", "atype": "S"},
    {"smarts": "[#17]", "atype": "Cl"},
    {"smarts": "[#20]", "atype": "Ca"},
    {"smarts": "[#25]", "atype": "Mn"},
    {"smarts": "[#26]", "atype": "Fe"},
    {"smarts": "[#30]", "atype": "Zn"},
    {"smarts": "[#35]", "atype": "Br"},
    {"smarts": "[#53]", "atype": "I"},
    {"smarts": "[#7X3v3][a]", "atype": "N", "comment": "pyrrole, aniline"},
    {"smarts": "[#7X3v3][#6X3v4]", "atype": "N", "comment": "amide"},
    {"smarts": "[#7+1]", "atype": "N", "comment": "ammonium, pyridinium"},
    {"smarts": "[SX2]", "atype": "SA", "comment": "sulfur acceptor"},
    {"smarts": "[#1][#6X3]:[#6X3]([#6X4])[#7]:[#7][#7][#6X4]", "atype": "HD", "comment": "4,5-H in 1,2,3-triazole"},
    ]
    }
    }

    try: 
        mol = Chem.MolFromSmiles(smiles)
        mol_hs = Chem.AddHs(mol)
        props = AllChem.ETKDG()
        props.pruneRmsThresh = 0.25
        props.useRandomCoords = True
        props.numThreads = 1
        # Create the object containing the conformers
        confs = AllChem.EmbedMultipleConfs(mol_hs,50,props)
        ps = AllChem.MMFFGetMoleculeProperties(mol_hs)

        # Prepare the confromers
        AllChem.MMFFOptimizeMoleculeConfs(mol_hs, mmffVariant='MMFF94', maxIters=10)
        conformers_energies_dict={} # Empty dictionary to store conformers 

        for conf in confs:
            try:
                filename='conformer-' + str(conf + 1) + '.pdb'
                ff = AllChem.MMFFGetMoleculeForceField(mol_hs,ps,confId=conf)
                ff.Minimize()
                energy_value = ff.CalcEnergy()
                conformers_energies_dict[conf] = energy_value
                # Generates a sorted dictionary containing the conformers
            except:
                continue

        conformers_energies_dict_sorted=sorted(conformers_energies_dict.items(), key=lambda x: x[1])
        lowest_energy_conformer_nbr = conformers_energies_dict_sorted[0][0]
        selected_conformer = Chem.MolToMolBlock(mol_hs,confId=lowest_energy_conformer_nbr)
        selected_mol = Chem.MolFromMolBlock(selected_conformer, removeHs=False)
        
        preparator = MoleculePreparation(atom_type_smarts=dictio,merge_these_atom_types=['H'])
        mol_setups = preparator.prepare(selected_mol)
    
        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
   
        pdbqt_string_ready = pdbqt_string.replace('UNL','MOL')

        inchi_key = compute_inchi_key(smiles)

        with open(f'/tmp/{inchi_key}.pdbqt','w') as pdbqt_file:
            pdbqt_file.write(pdbqt_string_ready)

        destination_dir = '/tmp'
        pdbqt_file = f'{destination_dir}/{inchi_key}.pdbqt'
        blob_pdbqt_file = f'{pdbqt_file}.tar.gz'
    
        tar = tarfile.open(blob_pdbqt_file, "w:gz")
        tar.add(pdbqt_file,arcname=f"{inchi_key}")
        tar.close()
    
        return blob_pdbqt_file

    except:
        return "ERROR_COMPUTING_PDBQT"
        pass

def append_ligand_blob_object_to_table(db_name,table_name,smiles):
    
    inchi_key = compute_inchi_key(smiles)
    blob_pdbqt_file = compute_pdqbt_file(smiles).split('/')[-1]

    if blob_pdbqt_file == "ERROR_COMPUTING_PDBQT":
        file_blob_object = "ERROR"
    else:
        file_blob_object = create_blob_object('/tmp/'+blob_pdbqt_file)

    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    sql = f"""INSERT INTO '{table_name}' (SMILES, inchi_key, stored_file, blob_item)
              VALUES (?,?,?,?);"""             
    
    data_tuple = (smiles,inchi_key,blob_pdbqt_file,file_blob_object)
    
    try: 
        cursor.execute(sql,data_tuple)
        conn.commit()
        cursor.close()
        #print(colored("BLOB data added successfully","green"))
    except sqlite3.Error as error:
        print(error)
        print("Failed to insert blob data into sqlite table", colored(f"{error}","red"))
        cursor.close()

def process_ligands_table_single_proc(origin_db,table_name,dest_database):
    conn = sqlite3.connect(origin_db)
    sql = f"""SELECT SMILES, inchi_key
              FROM {table_name};"""
    
    try: 
        df = pd.read_sql_query(sql,conn)
        
        # Create the table to store the corresponding ligands
        dest_table_name = table_name+'_pdbqt_ligands'
        create_ligands_blob_table(dest_database, dest_table_name)
        
        # Initiatialize tqdm for running
        tqdm.pandas() # initialize tqdm for pandas
        try: 
            df['SMILES'].progress_apply(lambda smiles: append_ligand_blob_object_to_table(dest_database,dest_table_name,smiles))
        except Exception as error:
            print(error)
            pass
    
    except Exception as error:
        print(error)
        print(colored("Error","red"))

def process_ligands_table(origin_db,table_name,dest_database):
    
    conn = sqlite3.connect(origin_db)
    sql = f"""SELECT SMILES, inchi_key
              FROM '{table_name}';"""
    
    try: 
        df = pd.read_sql_query(sql,conn)
        
        # Create the table to store the corresponding ligands
        dest_table_name = table_name+'_pdbqt_ligands'
        create_ligands_blob_table(dest_database, dest_table_name)
        
        # Initiatialize tqdm for running
        #tqdm.pandas() # initialize tqdm for pandas
        pandarallel.initialize(progress_bar=True)
        try: 
            df['SMILES'].parallel_apply(lambda smiles: append_ligand_blob_object_to_table(dest_database,dest_table_name,smiles))
        except Exception as error:
            print(error)
            pass
        print(colored(" PDQBT files processes and storeed as BLOB data successfully","green"))
    except Exception as error:
        print(error)
        print(colored("Error","red"))

def write_blob_object(data,dest_file):
    try: 
        with open(f'{dest_file}', 'wb') as file:
            file.write(data)
        file.close()
    
    except:
        print(f"Problem storing {dest_file}. Check for errors.")

def restore_ligands_from_table(db_name,table_name,dest_dir):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    sql = f"""SELECT * 
             FROM '{table_name}';"""

    cursor.execute(sql)
    query_output = cursor.fetchall()

    for row in query_output:
        smiles = row[0]
        inchi_key = row[1]
        data = row[3]
        dest_file = dest_dir+'/'+inchi_key+'.pdbqt.tar.gz'
        write_blob_object(data,dest_file)
        ## This will extract the file in the specified folder
        file = tarfile.open(dest_file)
        file.extractall(dest_dir)
        file.close()
        # This will rename the stored file to .pdbqt
        os.rename(f'{dest_dir}/{inchi_key}',f'{dest_dir}/{inchi_key}.pdbqt')
        # This will delete de stored compress version of the file
        os.remove(dest_file)


#################

if __name__ == '__main__':

    print("This is a test")

    process_ligands_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects/sample.db","aminoacids","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/level_1_projects/sample.db")