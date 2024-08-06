import shutil

def retrive_receptor_file(file,output_folder):
    """
    This function will retieve the receptor file based on the docking registries
    """

    receptor_file = file.split('/')[-1]

    shutil.copyfile(file, f'{output_folder}/{receptor_file}')
    
    receptor_ready_file = delete_receptor_hidrogens(file,output_folder)
    
    return receptor_ready_file


def delete_receptor_hidrogens(file,output_folder):
    """"
    This function will delete all Hidrogens present in a pdb file. The deletion of the hidrogens is required for further processing the receptor template by tleap, since otherwise error in Hs names occurs
    """
    
    hydrogen_types = ["HE","HN","HZ","HG","HA","HA1","HA2","HB","HN1","HN2","HN3","HA","HB1","HB2","HB3","HE1","HE2","HE3","HZ1","HZ2","HZ3","HG1","HG2","HD1","HD2","HH1","HH2","HH3","1HG1","2HG1","3HG1","1HH1","2HH1","1HE2","2HE2","1HG2","2HG2","3HG2","1HH2","2HH2","1HD2","2HD2","3HD2","1HD1","2HD1","3HD1","3HD2","HH","1HZ3","2HZ3","3HZ3","HO","HC"]
    
    with open(f'{output_folder}/receptor_ready.pdb','w') as receptor_output:
        for line in open(file):
            new_line = line.rstrip().split()
        
            # Check if a TER card is present and mantain it
            if new_line[0] == "TER":
                receptor_output.write(f'{line}')
                continue
            
            try: 
                # Check if the atom field (#2) contains a Hs and if so do not write the line to the output file
                if new_line[2] not in hydrogen_types and len(new_line) > 3:
                    receptor_output.write(f'{line}')
            except:
                continue
 
        receptor_output.close()
        
    return f'{output_folder}/receptor_ready.pdb'
            