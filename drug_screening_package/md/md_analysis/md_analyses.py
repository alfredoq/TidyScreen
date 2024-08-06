import MDAnalysis
import MDAnalysis.analysis.rms
from MDAnalysis.analysis import distances
import os
from termcolor import colored
import matplotlib.pyplot as plt
import subprocess

def create_analysis_folder(md_assay_folder):
    """
    This function will create a folder named 'analysis_results' in which all the files associated to the analyses will be stored.
    
    """

    results_folder = f'{md_assay_folder}/analysis_results'

    isExist = os.path.exists(results_folder)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(results_folder)
        print(colored("Analysis folder created","green"))
    
    return results_folder

def check_if_solvated_exists(md_assay_folder):
    """
    This function will check if the solvated trajectories corresponding to the heating, equilibration and production stages exists. In that case the function will return a 1 flag that is used to process the striping of the corresponding trajectories.
    """

    heat_1_exists = os.path.isfile(f"{md_assay_folder}/heat1.nc")
    heat_2_exists = os.path.isfile(f"{md_assay_folder}/heat2.nc")
    equi_exists = os.path.isfile(f"{md_assay_folder}/equi.nc")
    prod_exists = os.path.isfile(f"{md_assay_folder}/prod.nc")

    # If all the files exist, return a flag = 1
    if heat_1_exists and heat_2_exists and equi_exists and prod_exists:
        return 1
    else:
        return 0

def strip_trajectories(md_assay_folder):
    """"
    This function will be activated when the presence of solvated trajectories is detected. In that case, the striping process included in the function body will be activated.
    """

    print(colored("START stripping trajectories","green"))

    # Generate the striped .prmtop file

    with open(f'{md_assay_folder}/strip_prmtop.in','w') as strip_prmtop_file:
        strip_prmtop_file.write(f"parm {md_assay_folder}/complex.prmtop \n")
        strip_prmtop_file.write("parmstrip :WAT,Na+,Cl- \n")
        strip_prmtop_file.write(f"parmwrite out {md_assay_folder}/complex_strip.prmtop \n")
        strip_prmtop_file.write("go \n")
        strip_prmtop_file.write("quit")
    strip_prmtop_file.close()

    # Execute the stripping of the prmtop file
    bash_command = f'cpptraj -i {md_assay_folder}/strip_prmtop.in'
    process = subprocess.Popen(bash_command.split(),stdout=subprocess.PIPE)
    output, error = process.communicate()

    # Strip the trajectory set
    trajectories_list = ["heat1.nc","heat2.nc","equi.nc","prod.nc"]
    for trajectory in trajectories_list:
        traj_prefix = trajectory.replace('.nc','')
        with open(f'{md_assay_folder}/strip_{traj_prefix}.in','w') as strip_file:
            strip_file.write(f'parm {md_assay_folder}/complex.prmtop \n')
            strip_file.write(f'trajin {md_assay_folder}/{trajectory} \n')
            strip_file.write(f'autoimage \n')
            strip_file.write(f'strip :WAT,Na+,Cl- \n')
            strip_file.write(f'trajout {md_assay_folder}/{traj_prefix}_strip.nc \n')
            strip_file.write(f'go \n')
            strip_file.write(f'quit')
        strip_file.close()

        # Execute the stripping of the trajectories
        bash_command = f'cpptraj -i {md_assay_folder}/strip_{traj_prefix}.in'
        process = subprocess.Popen(bash_command.split(),stdout=subprocess.PIPE)
        output, error = process.communicate()

        # Delete the solvated trajectory
        os.remove(f'{md_assay_folder}/{trajectory}')


    print(colored("END stripping trajectories","green"))

def load_trajectory_in_mda(md_assay_folder,trajectory,strip_trajectories_flag):
    
    """"
    This function will load the trajectory indicated using MDAnalysis, and further return an Universe object.
    
    ------
    Params:
    ------
    - md_assay_folder: the full path to the folder in which the trajectory to be process is stored.
    - traj_prefix: the prefix of the file corresponding to the trajectory to be processed: i.e. equi, prod, etc.
    - strip_trajectories_flag: if == 1 a striped trajectory will be loaded, otherwise the solvated one is loaded
    
    """
    
    if strip_trajectories_flag == 1:
        prmtop_file = f'{md_assay_folder}/complex_strip.prmtop'
    else:
        prmtop_file = f'{md_assay_folder}/complex.prmtop'
    
    u = MDAnalysis.Universe(prmtop_file,trajectory)
    return u

def compute_backbone_rmsd(md_assay_folder,results_folder,strip_trajectories_flag):
    """
    This function will compute the backbone rmsd values for all the segments of simulation, using as reference structure the first frame of the heating stage.
    
    """
    # The trajectories to be processed
    if strip_trajectories_flag == 1:
        trajectories = ['heat1_strip.nc','heat2_strip.nc','equi_strip.nc','prod_strip.nc']
        # This will load the reference trajectory 
        ref_traj = 'heat1_strip.nc'
    else:
        trajectories = ['heat1.nc','heat2.nc','equi.nc','prod.nc']
        # This will load the reference trajectory 
        ref_traj = 'heat1.nc'

    ref = load_trajectory_in_mda(md_assay_folder,f'{md_assay_folder}/{ref_traj}',strip_trajectories_flag)
    
    for trajectory in trajectories:
        try:
            
            traj = f'{md_assay_folder}/{trajectory}'
            traj_prefix = trajectory.replace('.nc','')
            # Load the trajectory
            u = load_trajectory_in_mda(md_assay_folder,traj,strip_trajectories_flag)
            # Compute RMSD for the backbone
            R = MDAnalysis.analysis.rms.RMSD(u, ref,select="backbone")
            R.run()
            rmsd = R.results.rmsd.T
            rmsd_values = rmsd[2]
            time = rmsd[1]
        
            # Create a plot and store it
            plt.plot(time,rmsd_values)
            plt.title(f"Backbone RMSD corresponding to {trajectory}")
            plt.savefig(f'{results_folder}/RMSD-{traj_prefix}.png')
            plt.clf() # This will clear the plot for the next iteration
    
        except:
            continue
    
    print(colored("Receptor Backbone RMSD plots saved","green"))

def compute_pair_distances(md_assay_folder,results_folder,residue_atoms_pairs,strip_trajectories_flag):
    """
    This function will compute, plot and store pairs of distances
    """

    # The trajectories to be processed
    if strip_trajectories_flag == 1:
        trajectories = ['heat1_strip.nc','heat2_strip.nc','equi_strip.nc','prod_strip.nc']
    else:
        trajectories = ['heat1.nc','heat2.nc','equi.nc','prod.nc']
    
    for trajectory in trajectories:
        try:
            traj = f'{md_assay_folder}/{trajectory}'
            traj_prefix = trajectory.replace('.nc','')
            # Load the trajectory
            u = load_trajectory_in_mda(md_assay_folder,traj,strip_trajectories_flag)
        
            counter_plot = 1        
            
            # Clear these variables
            previous_res1 = ''
            previous_res2 = ''
            
            for int_num, pair_dict in residue_atoms_pairs.items():
                counter = 0
                for resnum, resname in pair_dict.items():
                    if counter == 0:
                        value_1 = resnum
                        value_2 = resname
                        counter+=1
                    if counter == 1:
                        value_3 = resnum
                        value_4 = resname
            
                # Compute the distance for a pair of atoms
                time_list, dist_list = calculate_distance_between_atoms(u,value_1,value_2,value_3,value_4)
                # Create a plot and store them for pairs of atoms
                plt.plot(time_list,dist_list)
                plt.title(f"Distances for pair: {pair_dict}")
                plt.savefig(f'{results_folder}/Atom-Atom-distance-{traj_prefix}_{counter_plot}.png')
                plt.clf() # This will clear the plot for the next iteration
                
                # Check if variables have changed, and in that case compute the corresponding differences
                if value_1 != previous_res1 and value_3 != previous_res2:
                    # Compute the distance for a pair of residues
                    time_list, dist_list = calculate_distance_between_residues(u,value_1,value_3)
                    # Create a plot and store them for pairs of atoms
                    plt.plot(time_list,dist_list)
                    plt.title(f"Distances for pair: {value_1}:{value_3}")
                    plt.savefig(f'{results_folder}/Residues-{value_1}-{value_3}-distance-{traj_prefix}.png')
                    plt.clf() # This will clear the plot for the next iteration
                
                ## Assign resname to variable to check if they have changed. If not, the res-res distances will not be computed again since they are the same
                previous_res1 = value_1
                previous_res2 = value_3
                
                # Increase the plot counter
                counter_plot+=1
        
        except Exception as error:
            continue

def calculate_distance_between_atoms(u,res1_nbr,atom_name_1,res2_nbr,atom_name_2):
    """
    This will compute and return the distance as a list
    """

    time_list = []
    distances_list = []
    
    for ts in u.trajectory:
        time = ts.time
        atom1 = u.select_atoms(f"resid {res1_nbr} and name {atom_name_1}")
        atom2 = u.select_atoms(f"resid {res2_nbr} and name {atom_name_2}")
        com1 = atom1.center_of_mass()
        com2 = atom2.center_of_mass()
        dist = distances.distance_array(com1,com2)

        time_list.append(time)
        distances_list.append(dist[0])

    return time_list, distances_list

def calculate_distance_between_residues(u,res1_nbr,res2_nbr):
    """
    This will compute and return the distance as a list for the COM of two residues
    """

    time_list = []
    distances_list = []
    
    for ts in u.trajectory:
        time = ts.time
        res1 = u.select_atoms(f"resid {res1_nbr}")
        res2 = u.select_atoms(f"resid {res2_nbr}")
        com1 = res1.center_of_mass()
        com2 = res2.center_of_mass()
        dist = distances.distance_array(com1,com2)

        time_list.append(time)
        distances_list.append(dist[0])

    return time_list, distances_list

def check_temperature(md_assay_folder,results_folder):
    """
    This function will plot the temperatures of the different simulation stages
    
    """    
    # Outputs to be processed
    stages = ['heat1.out','heat2.out','equi.out','prod.out']
    
    for stage in stages:
        time_list = []
        temp_list = []
        # Add a try condition to avoid errors in case a later stage has not been completed
        try: 
            for line in open(f'{md_assay_folder}/{stage}'):
                new_line = line.rsplit()
                # The following will end data acquisition after finding the word 'AVERAGES'
                if 'A' in new_line and 'V' in new_line and 'E' in new_line:
                    break # End this internal for loop
                
                if "NSTEP" in new_line and "TEMP(K)" in new_line:
                    time_list.append(float(new_line[5]))
                    temp_list.append(float(new_line[8]))
            
            # Create a plot and store it
            plt.plot(time_list,temp_list)
            plt.title(f"Temperature profile vs time")
            plt.savefig(f'{results_folder}/Temp-vs-Time-{stage}.png')        
            plt.clf() # Clear the figure
            
        except:
            continue
   
def perform_overall_md_analyses(mds_path,md_assay_number,residue_atoms_pairs,strip_trajectories_flag):
    """
    This function will perform a full set of analyses on the trajectories obtained for a given MD assay
    
    ------
    Parameters:
    ------
    - mds_path: the general path in wich the md assays are to be stored
    - md_assay_number: the number of the MD assay to be processed
    - residue_atoms_pairs: a dictionary containing the distances to be computed, with each distance consisting of res_nbr:atom_name pairs
    - strip_trajectories_flag: if 1, the solvated trajectories will be stripped and deleted afterwards
    """
    
    md_assay_folder = f'{mds_path}/mds_assays/md_assay_{md_assay_number}'
    
    ## Create the folder to store the results
    results_folder = create_analysis_folder(md_assay_folder)    
    ## Check simulation temperature
    check_temperature(md_assay_folder,results_folder)
    ## Check the presence of solvated trajectories
    solvated_exists = check_if_solvated_exists(md_assay_folder)
    ## Action to take if solvated exist
    
    if solvated_exists == 1 and strip_trajectories_flag == 1:
        strip_trajectories(md_assay_folder)
    
    ## First analysis that will occurr on the stripped trajectory
    compute_backbone_rmsd(md_assay_folder,results_folder,strip_trajectories_flag)
    ## Compute paired distances
    compute_pair_distances(md_assay_folder,results_folder,residue_atoms_pairs,strip_trajectories_flag)
