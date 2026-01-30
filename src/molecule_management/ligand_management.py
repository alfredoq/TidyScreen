

def restore_single_docked_pose(results_db, ligname, pose_id):

    # Retrieve input_model and pose_coords_json
    input_model, pose_coords_json = _retrieve_pose_info(results_db, pose_id)

    # Process input_model and pose_coords_json to create pdb_dict
    pdb_dict = _process_input_model_for_moldf(input_model, pose_coords_json)

    return pdb_dict


def _retrieve_pose_info(db_path, pose_id):
        """
        Retrieves ligand pose information from the database for a given pose ID.
        Connects to the specified SQLite database, queries the Ligands and Results tables
        to obtain the input model and ligand coordinates associated with the provided pose ID.
        Args:
            db_path (str): Path to the SQLite database file.
            pose_id (int or str): Identifier of the pose to retrieve.
        Returns:
            tuple: A tuple containing:
                - input_model (str): The input model data for the ligand.
                - pose_coords_json (str): JSON string of the ligand coordinates.
        Raises:
            sqlite3.DatabaseError: If there is an error connecting to or querying the database.
            TypeError: If the query does not return any results.
        """
       
        import sqlite3

        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Retrieve the info required to reconstruct the ligand pose
        cursor.execute("""
        SELECT L.input_model, R.ligand_coordinates FROM Ligands L
        JOIN Results R ON L.LigName = R.LigName
        WHERE R.Pose_ID = ?
        """, (pose_id,))
        row = cursor.fetchone()

        input_model, pose_coords_json = row
        
        return input_model, pose_coords_json

def _process_input_model_for_moldf(input_model, coords_json):
        """
        Processes an input molecular model and coordinates to generate a DataFrame suitable for moldf.
        This method parses the input model (as a string representation of a list), extracts atom information,
        replaces atomic coordinates with those provided in `coords_json`, and formats the data into a DataFrame
        with standardized columns for molecular data. The resulting DataFrame is returned in a dictionary under
        the key '_atom_site', suitable for use with moldf.
        Args:
            input_model (str): String representation of a list containing atom records. Each record can be a string
                or a list/tuple, where the first element should be 'ATOM'.
            coords_json (str): JSON-encoded list of new coordinates (list of [x, y, z]) to replace the original
                coordinates in the input model.
        Returns:
            dict: A dictionary with a single key '_atom_site', whose value is a pandas DataFrame containing the
                processed atom data with updated coordinates and standardized columns.
        """

        import ast
        import pandas as pd
        import json


        input_model_list = ast.literal_eval(input_model)
        atom_rows = []
        for item in input_model_list:
            # If item is a string, split it and check if first element is 'ATOM'
            if isinstance(item, str):
                parts = item.split()
                if len(parts) > 0 and parts[0] == 'ATOM':
                    atom_rows.append(parts)
            # If item is a list/tuple, check if first element is 'ATOM'
            elif isinstance(item, (list, tuple)) and len(item) > 0 and item[0] == 'ATOM':
                atom_rows.append(item)
        
        
        df = pd.DataFrame(atom_rows)
        
        # Rename columns using a list of column names
        column_names = [
            'record_name', 'atom_number', 'atom_name', 'residue_name',
            'residue_number', 'x_coord', 'y_coord', 'z_coord',
            'occupancy', 'b_factor', 'charge', 'element_symbol'
        ]
        df.columns = column_names[:len(df.columns)]
        
        # Insert an empty column 'alt_loc' at the fourth position (index 3)
        df.insert(3, 'alt_loc', '')
        df.insert(5, 'chain_id', '')
        df.insert(7, 'insertion', '')
        df.insert(13, 'segment_id', '')

        # Move 'charge' column to the last position
        if 'charge' in df.columns:
            charge_col = df.pop('charge')
            df['charge'] = charge_col

        # Replace 'charge' column with empty items of type object
        if 'charge' in df.columns:
            df['charge'] = pd.Series([''] * len(df), dtype=object)

        ## Replace the coordinates with the ones from coords_json for the given pose
        coords_list = json.loads(coords_json)

        for i, coord in enumerate(coords_list):
            # Replace the original coordinates with the new ones corresponding to the pose
            df.at[i, 'x_coord'] = coord[0]
            df.at[i, 'y_coord'] = coord[1]
            df.at[i, 'z_coord'] = coord[2]  

            # Attempt to convert columns to appropriate types
        for col in df.columns:
            # Try to convert to int, then float, else keep as string
            try:
                df[col] = df[col].astype(int)
            except ValueError:
                try:
                    df[col] = df[col].astype(float)
                except ValueError:
                    pass

        # Return as dictionary for moldf
        return {'_atom_site': df}
    
    

def write_pdb_with_moldf(pdb_dict, ligname, run_number, output_dir):
        """
        Writes a PDB file using the provided molecular data dictionary and saves it to the specified output directory.

        Args:
            pdb_dict (dict): Dictionary containing the molecular structure data to be written to the PDB file.
            ligname (str): Name of the ligand, used as part of the output filename.
            run_number (int or str): Identifier for the run, used as part of the output filename.
            output_dir (str): Directory path where the output PDB file will be saved.

        Returns:
            None

        Side Effects:
            Creates a PDB file in the specified output directory with a filename formatted as "{ligname}_{run_number}.pdb".
        """
        from moldf import write_pdb
        import os

        output_file = os.path.join(output_dir, f"{ligname}_{run_number}.pdb")
        write_pdb(pdb_dict, output_file)

        return output_file
    
def prepare_ligand_tleap_input_files(chemspace_db, ligname, assay_info, output_dir):
        
        import subprocess
        import sys
        import os
        
        # Retrieve original ligand in .sdf format and prepare a .pdb file
        pdb_file = _convert_ligand_sdf_into_pdb(chemspace_db, ligname, assay_info, output_dir)
        
        if pdb_file is None:
            print(f"[ERROR] Could not convert ligand {ligname} SDF to PDB")
            sys.exit(1)
            
        # Prepare the .prepin and .frcmod files using antechamber and parmchk2
        mol2_file = os.path.join(output_dir, f"{ligname}.mol2")
        frcmod_file = os.path.join(output_dir, f"{ligname}.frcmod")
      
        ## Compute ligand espaloma charges
        espaloma_output_file = _compute_ligand_espaloma_charges(pdb_file)
        
        # Run antechamber
        antechamber_command = f"antechamber -i {pdb_file} -fi pdb -o {mol2_file} -fo mol2 -c rc -cf {espaloma_output_file} "
        
        try:
            subprocess.run(antechamber_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        except Exception as e:
            print(f"[ERROR] Failed to run antechamber: {e}")
            return None

        parmchk2_command = f"parmchk2 -i {mol2_file} -f mol2 -o {frcmod_file}"
        try:
            subprocess.run(parmchk2_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        except Exception as e:
            print(f"[ERROR] Failed to run parmchk2: {e}")
            return None

        return mol2_file, frcmod_file
    
def _convert_ligand_sdf_into_pdb(chemspace_db, ligname, assay_info, output_dir):
    """
    Convert a ligand SDF file into a PDB file using RDKit.

    Args:
        ligname (str): Name of the ligand (used for file naming).
        assay_info (dict): Assay information dictionary containing paths and metadata.
        output_dir (str): Directory where the output PDB file should be saved.

    Returns:
        str: Path to the generated PDB file, or None if conversion failed.
    """
    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import sqlite3
    from moldf import read_pdb
    from moldf import write_pdb

    # Construct input SDF and output PDB file paths
    sdf_file = os.path.join(output_dir, f"{ligname}.sdf")
    pdb_file = os.path.join(output_dir, f"{ligname}.pdb")

    # Retrieve the ligand in .sdf format from the chemspace.db corresponding table
    chemspace_table = assay_info.get('table_name', None)
    if not chemspace_table:
        print("[ERROR] ChemSpace table name not provided in assay_info.")
        return None

    try:
        conn = sqlite3.connect(chemspace_db)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT sdf_blob FROM {} WHERE inchi_key = ?".format(chemspace_table),
            (ligname,)
        )
        row = cursor.fetchone()
        conn.close()
        if not row or not row[0]:
            print(f"[ERROR] Ligand {ligname} not found in table {chemspace_table}.")
            return None
        with open(sdf_file, "wb") as f:
            f.write(row[0])
    except Exception as e:
        print(f"[ERROR] Failed to retrieve ligand SDF from ChemSpace: {e}")
        return None

    if not os.path.isfile(sdf_file):
        print(f"[ERROR] SDF file not found: {sdf_file}")
        return None

    # Read molecule from SDF
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        print(f"[ERROR] Failed to read molecule from SDF: {sdf_file}")
        return None

    # Restore atom names
    if mol.HasProp("AtomNames"):
        atom_names = mol.GetProp("AtomNames").split("|")
        
        for i, atom_name in enumerate(atom_names):
            if i < mol.GetNumAtoms():
                atom = mol.GetAtomWithIdx(i)
                atom.SetProp("_Name", atom_name)
                atom.SetProp("atomLabel", atom_name)

    # Write to PDB
    try:
        Chem.MolToPDBFile(mol, pdb_file)
    
    except Exception as e:
        print(f"[ERROR] Failed to write PDB file: {e}")
        return None
    
    # Load the PDB file into a moldf dataframe
    try:
        # Insert a blank line at at the  pdb_file for moldf correct parsing
        with open(pdb_file, 'r') as original:
            data = original.read()
        with open(pdb_file, 'w') as modified:
            modified.write('\n' + data) 
        
        # Parse the file using moldf
        moldf_dict = read_pdb(pdb_file)
        moldf_df = moldf_dict['_atom_site']
        
        # Check if the length of 'atom_names' list is equal to the number of rows in moldf_df (number of atoms should be the same)
        if len(atom_names) != len(moldf_df):
            print(f"[ERROR] Number of atom names ({len(atom_names)}) does not match number of atoms in PDB ({len(moldf_df)}).")
            sys.exit(1)
        else:
            # Assign atom names to the moldf dataframe
            moldf_df['atom_name'] = atom_names
            # Assign empty value to 'insertion' column
            moldf_df['insertion'] = ''
        
        new_moldf_dict = {
            '_atom_site': moldf_df
        }
        
        # Write (overwrite) the tleap_processed_file
        from moldf import write_pdb
        
        write_pdb(new_moldf_dict, pdb_file)
        
        # Remove one space at column 20 of pdb_file
        with open(pdb_file, 'r') as original:
            data = original.readlines()
        with open(pdb_file, 'w') as modified:
            for line in data:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Remove one space at column 23 (index 22)
                    new_line = line[:22] + line[23:]
                    modified.write(new_line)
                else:
                    modified.write(line)
        
        return pdb_file
        
        
    except Exception as e:
        print(f"[ERROR] Failed to load PDB file into moldf dataframe: {e}")
        return None
    
    
def _compute_ligand_espaloma_charges(pdb_file):
    from rdkit import Chem
    from espaloma_charge import charge
    import numpy as np
    
    # Read the ligand PDB file into an RDKit molecule
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        print(f"[ERROR] Could not read molecule from PDB: {pdb_file}")
        return None
    
    espaloma_charges = charge(mol)
    
    # Write espaloma charges to a file with 8 columns, values separated by a single space
    espaloma_output_file = pdb_file.replace('.pdb', '_espaloma_charges.txt')
    with open(espaloma_output_file, 'w') as f:
        for i, val in enumerate(espaloma_charges):
            f.write(f"{val:.6f}")
            # Write a space after each value except the last in a line
            if (i + 1) % 8 == 0:
                f.write("\n")
            else:
                f.write(" ")
        # If the last line is not complete, add a newline
        if len(espaloma_charges) % 8 != 0:
            f.write("\n")

    return espaloma_output_file