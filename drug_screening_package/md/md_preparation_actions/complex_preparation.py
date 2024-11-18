import subprocess

def prepare_complex(md_assay_folder,receptor_file,ligand_file):
    """
    This function will combine a recetor and ligand .pdb files into the corresponding complex for further processing.
    
    ------
    Parameters
    ------
    - receptor_file: the full path to the receptor used in docking and that is consistant with docked sub_pose to be used.
    - ligand_file: the full path to the docked sub_pose to be studied.
    
    """
    ligand_prefix = ligand_file.split('/')[-1].replace('.pdb','')
    complex_pdb_file = f'complex_{ligand_prefix}.pdb'
    with open(f'{md_assay_folder}/{complex_pdb_file}','w') as complex_file:
        # append receptor lines to the complex file
        for line in open(receptor_file):
            new_line = line.rstrip().split()
            
            if new_line[0] == 'ATOM':
                complex_file.write(line)
            
        # append a TER card to the complex file
        complex_file.write("TER \n")
        
        # append ligand lines to the complex file
        
        for line in open(ligand_file):
            new_line = line.rstrip().split()
            
            if new_line[0] == 'ATOM' or new_line[0] == 'HETATM':
                complex_file.write(line)
        
        # append a TER and END card to the complex file
        complex_file.write("TER \n")
        complex_file.write("END")
        
    complex_file.close()
    
    return f'{md_assay_folder}/{complex_pdb_file}'

def write_tleap_input(md_assay_folder,pdb_complex_file,md_cond_dict):
    """"
    This function will write a tleap input file used to prepare to .prmtop and .inpcrd files corresponding to the input complex.
    
    ------
    Parameters
    ------
    - md_assay_folder: The full path to the assay folder in which the md studies are being stored.
    - pdb_complex_file: the full path to the pdb file corresponding to the ligand protein complex to be simulated.
    - md_cond_dict: the set of parameters to be used for the MD simulations. Provided as a dictionary of parameters.
    """
        
    ligand_base_prefix = pdb_complex_file.split('/')[-1].split('_')[1]
    
    tleap_file = f'{md_assay_folder}/tleap.in'
    
    protein_forcefield = md_cond_dict['tleap_params']['protein_forcefield']
    ligand_forcefield = md_cond_dict['tleap_params']['ligand_forcefield']
    water_forcefield = md_cond_dict['tleap_params']['water_forcefield']
    solvate_geometry = md_cond_dict['tleap_params']['solvate_geometry']
    solvent_model = md_cond_dict['tleap_params']['solvent_model']
    solvent_box_min_dist = md_cond_dict['tleap_params']['solvent_box_min_dist']
    solvent_min_dist = md_cond_dict['tleap_params']['solvent_min_dist']
    
    
    with open(tleap_file,'w') as tleap_input:
        tleap_input.write(f"source {protein_forcefield} \n")
        tleap_input.write(f"source {ligand_forcefield} \n")
        tleap_input.write(f"source {water_forcefield} \n")
        tleap_input.write(f"loadamberprep {md_assay_folder}/{ligand_base_prefix}.prepin \n")
        tleap_input.write(f"loadamberparams {md_assay_folder}/{ligand_base_prefix}.frcmod \n")
        tleap_input.write(f"COM = loadpdb {pdb_complex_file} \n")
        tleap_input.write(f"addions COM Na+ 0 \n")
        tleap_input.write(f"addions COM Cl- 0 \n")
        tleap_input.write(f"{solvate_geometry} COM {solvent_model} {solvent_box_min_dist} {solvent_min_dist} \n")
        tleap_input.write(f"saveamberparm COM  {md_assay_folder}/complex.prmtop {md_assay_folder}/complex.inpcrd \n")
        tleap_input.write(f"quit \n")
    
    tleap_input.close()
    
    # Prepare the .frcmod file
    bash_command = f'tleap -f {tleap_file}'
    process = subprocess.Popen(bash_command.split(),stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    