import sqlite3
import tarfile
import sys
import math
import subprocess
from glob import glob
from termcolor import colored
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)

import docking.ligands_management.process_ligands_table as lig_proc
import docking.registers_management.docking_assays_registers as dock_regs
import docking.registers_management.docking_params_registers as dock_parm_regs
import docking.perform_docking.docking_configuration as dock_cfgs
from docking.receptor_management import receptor_management as rec_mng
from docking.input_output_operations import actions_on_dbs as db_acts

def prepare_docking_exec_script_OBSOLETE(docking_assay_folder,custom_string):
    """
    This script will prepare the execution script to be run on bash for a whole docking assay.
    """
    list_of_ligands = glob(f'{docking_assay_folder}/pdbqt_files/*.pdbqt')
    fld_file = glob(f'{docking_assay_folder}/receptor/*.fld')[0].split('/')[-1]
    
    with open(f'{docking_assay_folder}/docking_execution.sh','a') as exec_file:
        counter = 1
        for ligand in list_of_ligands:
            # Append instruction at start of file to register initial time
            if counter == 1:
                exec_file.write("start=$(date +%s)\n")
                exec_file.write("date > timing.txt\n")
                counter+=1
                                    
            ligand_file = ligand.split('/')[-1]
            ligand_prefix = ligand_file.replace('.pdbqt','')
            exec_file.write(f'autodock_gpu_128wi --lfile ./pdbqt_files/{ligand_file} --ffile ./receptor/{fld_file} {custom_string}\n')
            #exec_file.write(f'autodock_gpu_128wi --xml2dlg ./pdbqt_files/{ligand_prefix}.xml --contact_analysis 1\n')
            exec_file.write(f'mv ./pdbqt_files/{ligand_prefix}.dlg ./dlgs\n')
        print(colored("SUCESSFULLY written the docking execution script","green"))

        # Append the final timings instruction
        exec_file.write("end=$(date +%s)\n")
        exec_file.write("date >> timing.txt\n")
        exec_file.write('echo "It took $(($end - $start)) seconds to finish." >> timing.txt\n')

    subprocess.call(['chmod','777',f'{docking_assay_folder}/docking_execution.sh'])

def prepare_docking_exec_script(docking_assay_folder,custom_string,n_chunk):
    """
    This script will prepare the execution script to be run on bash for a whole docking assay.
    """
    list_of_ligands = glob(f'{docking_assay_folder}/pdbqt_files/*.pdbqt')
    fld_file = glob(f'{docking_assay_folder}/receptor/*.fld')[0].split('/')[-1]
    # Divide the originial lists in sublists 
    divided_list = [list_of_ligands[i * n_chunk:(i + 1) * n_chunk] for i in range((len(list_of_ligands) + n_chunk - 1) // n_chunk )]  
    chunk = 1
    for list_of_ligands in divided_list:
        counter = 1
        with open(f'{docking_assay_folder}/docking_execution_{chunk}.sh','a') as exec_file:
            for ligand in list_of_ligands:
        
         # Append instruction at start of file to register initial time
                if counter == 1:
                    exec_file.write("start=$(date +%s)\n")
                    exec_file.write(f"date > timing_{chunk}.txt\n")
                    counter += 1 # This will exit for the loop for this iteration

                ligand_file = ligand.split('/')[-1]
                ligand_prefix = ligand_file.replace('.pdbqt','')
                exec_file.write(f'autodock_gpu_128wi --lfile ./pdbqt_files/{ligand_file} --ffile ./receptor/{fld_file} {custom_string}\n')
                exec_file.write(f'mv ./pdbqt_files/{ligand_prefix}.dlg ./dlgs\n')

            # Append the final timings instruction
            exec_file.write("end=$(date +%s)\n")
            exec_file.write(f"date >> timing_{chunk}.txt\n")
            exec_file.write(f'echo "It took $(($end - $start)) seconds to finish." >> timing_{chunk}.txt\n')
            # Change execution permissions for the generated file
            subprocess.call(['chmod','777',f'{docking_assay_folder}/docking_execution_{chunk}.sh'])
        
        # Update the chunk value
        chunk +=1

    print(colored("SUCESSFULLY written the docking execution script","green"))

def prepare_docking_exec_script_mendieta(docking_assay_folder,custom_string,mendieta_assay_path):
    """
    This script will prepare the execution script to be run on bash for a whole docking assay for execution in Mendieta cluster.
    """

    list_of_ligands = glob(f'{docking_assay_folder}/pdbqt_files/*.pdbqt')
    fld_file = glob(f'{docking_assay_folder}/receptor/*.fld')[0].split('/')[-1]
    assay_folder = docking_assay_folder.split('/')[-1]

    with open(f'{docking_assay_folder}/docking_execution_SLURM.sh','a') as mendieta_exec_file:
        mendieta_exec_file.write("#!/bin/bash\n")
        mendieta_exec_file.write("### Las líneas #SBATCH configuran los recursos de la tarea\n")
        mendieta_exec_file.write("### (aunque parezcan estar comentadas)\n")
        mendieta_exec_file.write("### Cola\n")
        mendieta_exec_file.write("#SBATCH --partition=multi\n")
        mendieta_exec_file.write("### Tiempo de ejecucion. Formato dias-horas:minutos. Maximo una semana.\n")
        mendieta_exec_file.write("#SBATCH --time 2-0:00\n")
        mendieta_exec_file.write("### Nombre de la tarea\n")
        mendieta_exec_file.write("#SBATCH --job-name=test\n")
        mendieta_exec_file.write("### Cantidad de nodos a usar\n")
        mendieta_exec_file.write("#SBATCH --nodes=1\n")
        mendieta_exec_file.write("### GPUs por nodo (<= 2)\n")
        mendieta_exec_file.write("### OJO: Todos los procesos en un nodo ven ambas GPU con gpu:2!\n")
        mendieta_exec_file.write("###      Llamar al programa con /opt/mendieta/bin/split_cpu.sh para dividirlas.\n")
        mendieta_exec_file.write("#SBATCH --gres=gpu:1\n")
        mendieta_exec_file.write("### Procesos por nodo\n")
        mendieta_exec_file.write("#SBATCH --ntasks-per-node=1\n")
        mendieta_exec_file.write("### Cores por proceso (OpenMP/Pthreads/etc)\n")
        mendieta_exec_file.write("### Recordar exportar OMP_NUM_THREADS/MKL_NUM_THREADS/etc con el mismo valor\n")
        mendieta_exec_file.write("#SBATCH --cpus-per-task=10\n")
        mendieta_exec_file.write("### Environment setup\n")
        mendieta_exec_file.write(". /etc/profile\n")
        mendieta_exec_file.write("### Environment modules\n")
        mendieta_exec_file.write("module load gcc/11.2.0\n")
        mendieta_exec_file.write("module load cuda/11.6.2\n")
        mendieta_exec_file.write(f"export DOCK_RUN_PATH='{mendieta_assay_path}/{assay_folder}'\n")
        mendieta_exec_file.write("### Ejecutar la tarea\n")
        mendieta_exec_file.write("### NOTA: srun configura MVAPICH2 y MPICH con lo puesto arriba,\n")
        mendieta_exec_file.write("### no hay que llamar a mpirun.\n")
        mendieta_exec_file.write("srun docking_execution_mendieta.sh\n")


    with open(f'{docking_assay_folder}/docking_execution_mendieta.sh','a') as exec_file:
        counter = 1
        for ligand in list_of_ligands:
            if counter == 1:
                exec_file.write("#!/bin/bash\n")
                # Append instruction at start of file to register initial time
                exec_file.write("start=$(date +%s)\n")
                exec_file.write("date > timing_mendieta.txt\n")
                counter+=1
            ligand_file = ligand.split('/')[-1]
            ligand_prefix = ligand_file.replace('.pdbqt','')
            exec_file.write(f'autodock_gpu_128wi --lfile $DOCK_RUN_PATH/pdbqt_files/{ligand_file} --ffile $DOCK_RUN_PATH/receptor/{fld_file} {custom_string}\n')
            exec_file.write(f'mv $DOCK_RUN_PATH/pdbqt_files/{ligand_prefix}.dlg $DOCK_RUN_PATH/dlgs\n')
        
        # Append the final timings instruction
        exec_file.write("end=$(date +%s)\n")
        exec_file.write("date >> timing_mendieta.txt\n")
        exec_file.write('echo "It took $(($end - $start)) seconds to finish." >> timing_mendieta.txt\n')
        
        subprocess.call(['chmod','777',f'{docking_assay_folder}/docking_execution_mendieta.sh'])

        print(colored("SUCESSFULLY written the docking execution script Mendieta","green"))

def prepare_docking_exec_script_hpc(docking_assay_folder,custom_string,hpc_assay_path,n_chunk):
    """
    This script will prepare the execution script to be run on bash for a whole docking assay for execution in Mendieta cluster.
    """

    list_of_ligands = glob(f'{docking_assay_folder}/pdbqt_files/*.pdbqt')
    fld_file = glob(f'{docking_assay_folder}/receptor/*.fld')[0].split('/')[-1]
    assay_folder = docking_assay_folder.split('/')[-1]

    # Divide the originial lists in sublists 
    divided_list = [list_of_ligands[i * n_chunk:(i + 1) * n_chunk] for i in range((len(list_of_ligands) + n_chunk - 1) // n_chunk )]  
    number_of_chunks = len(divided_list)

    for chunk in range(number_of_chunks):
        part = chunk+1
        with open(f'{docking_assay_folder}/docking_execution_SLURM_{part}.sh','a') as mendieta_exec_file:
            mendieta_exec_file.write("#!/bin/bash\n")
            mendieta_exec_file.write("### Las líneas #SBATCH configuran los recursos de la tarea\n")
            mendieta_exec_file.write("### (aunque parezcan estar comentadas)\n")
            mendieta_exec_file.write("### Cola\n")
            mendieta_exec_file.write("#SBATCH --partition=multi\n")
            mendieta_exec_file.write("### Tiempo de ejecucion. Formato dias-horas:minutos. Maximo una semana.\n")
            mendieta_exec_file.write("#SBATCH --time 2-0:00\n")
            mendieta_exec_file.write("### Nombre de la tarea\n")
            mendieta_exec_file.write(f"#SBATCH --job-name=part_{part}\n")
            mendieta_exec_file.write("### Cantidad de nodos a usar\n")
            mendieta_exec_file.write("#SBATCH --nodes=1\n")
            mendieta_exec_file.write("### GPUs por nodo (<= 2)\n")
            mendieta_exec_file.write("### OJO: Todos los procesos en un nodo ven ambas GPU con gpu:2!\n")
            mendieta_exec_file.write("###      Llamar al programa con /opt/mendieta/bin/split_cpu.sh para dividirlas.\n")
            mendieta_exec_file.write("#SBATCH --gres=gpu:1\n")
            mendieta_exec_file.write("### Procesos por nodo\n")
            mendieta_exec_file.write("#SBATCH --ntasks-per-node=1\n")
            mendieta_exec_file.write("### Cores por proceso (OpenMP/Pthreads/etc)\n")
            mendieta_exec_file.write("### Recordar exportar OMP_NUM_THREADS/MKL_NUM_THREADS/etc con el mismo valor\n")
            mendieta_exec_file.write("#SBATCH --cpus-per-task=10\n")
            mendieta_exec_file.write("### Environment setup\n")
            mendieta_exec_file.write(". /etc/profile\n")
            mendieta_exec_file.write("### Environment modules\n")
            mendieta_exec_file.write("module load gcc/11.2.0\n")
            mendieta_exec_file.write("module load cuda/11.6.2\n")
            mendieta_exec_file.write(f"export DOCK_RUN_PATH='{hpc_assay_path}/{assay_folder}'\n")
            mendieta_exec_file.write("### Ejecutar la tarea\n")
            mendieta_exec_file.write("### NOTA: srun configura MVAPICH2 y MPICH con lo puesto arriba,\n")
            mendieta_exec_file.write("### no hay que llamar a mpirun.\n")
            mendieta_exec_file.write(f"srun docking_execution_mendieta_{part}.sh\n")

    chunk = 1
    for list_of_ligands in divided_list:
        with open(f'{docking_assay_folder}/docking_execution_hpc_{chunk}.sh','a') as exec_file:
            counter = 1
            for ligand in list_of_ligands:
                if counter == 1:
                    exec_file.write("#!/bin/bash\n")
                    # Append instruction at start of file to register initial time
                    exec_file.write("start=$(date +%s)\n")
                    exec_file.write(f"date > timing_mendieta_{chunk}.txt\n")
                    counter+=1
                ligand_file = ligand.split('/')[-1]
                ligand_prefix = ligand_file.replace('.pdbqt','')
                exec_file.write(f'autodock_gpu_128wi --lfile $DOCK_RUN_PATH/pdbqt_files/{ligand_file} --ffile $DOCK_RUN_PATH/receptor/{fld_file} {custom_string}\n')
                exec_file.write(f'mv $DOCK_RUN_PATH/pdbqt_files/{ligand_prefix}.dlg $DOCK_RUN_PATH/dlgs\n')
        
            # Append the final timings instruction
            exec_file.write("end=$(date +%s)\n")
            exec_file.write(f"date >> timing_hpc_{chunk}.txt\n")
            exec_file.write('echo "It took $(($end - $start)) seconds to finish." >> timing_hpc.txt\n')
            subprocess.call(['chmod','777',f'{docking_assay_folder}/docking_execution_hpc_{chunk}.sh'])

        # Update the chunk number for the next operation
        chunk+=1

    print(colored("SUCESSFULLY written the docking execution script in HPC","green"))

def create_docking_assay(receptor_models_registry_path,docking_assays_registry_path,docking_params_registry_path,docking_assays_storing_path,ligands_db,ligands_table_name,rec_model_id,docking_params_id,hpc_assay_path,n_chunk):
    """
    This function will create a tar.gz file containing all the information (ligands, receptor, parameters, execution script) required to perform a docking assay.

    ------
    Parameters:
    ------
    - receptor_models_registry_path: the full path to the database containing the receptor models.
    - docking_assays_registry_path: the full path to the folder in which the docking assays registries is stored.
    - docking_assays_storing_path: the path in which all the docking assays are stored
    - ligands_db: the full path to the database in which the ligands in the .pdbqt format are stored.
    - ligands_table_name: the table name containing the ligands in the .pdbqt format.
    - receptor_model_id: the rec_id stored in the corresponding receptors database.
    - docking_params_id: the id of the parameters set as stored in the docking_parameters database.
    - mendieta_assay_path: the folder in Mendieta were the docking assays are stored.
    - n_chunk: the number of ligands that will be processed in each script submission
    
    """
    ## This will create the assay folder structure and store the corresponding registry
    docking_assay_folder = dock_cfgs.create_docking_assay_folder(docking_assays_registry_path,docking_assays_storing_path,ligands_db,ligands_table_name,rec_model_id,docking_params_id)
    ## This will restore the ligands to the corresponding folder
    lig_proc.restore_ligands_from_table(ligands_db,ligands_table_name,f'{docking_assay_folder}/pdbqt_files')
    ## This will restore the receptor model to the corresponding folder
    file_full_path = rec_mng.restore_receptor_acc(receptor_models_registry_path,f'{docking_assay_folder}/receptor',rec_model_id)
    db_acts.uncompress_receptor_files(file_full_path,f'{docking_assay_folder}/receptor')
    ## Create a default params dict as reference of docking conditions
    default_dict = dock_parm_regs.create_default_dict()
    ## Retrieve the dictionary of docking conditions
    params_dict = dock_parm_regs.retrieve_docking_conditions_set(docking_params_registry_path,docking_params_id)
    ## Compare the dicts and return the required custom parameters
    custom_params_string = dock_parm_regs.compare_params_dicts(default_dict,params_dict)
    ## Prepare the execution script to run the whole docking assay
    #prepare_docking_exec_script(docking_assay_folder,custom_params_string)
    ## Prepare the execution script to run the whole docking assay
    prepare_docking_exec_script(docking_assay_folder,custom_params_string,n_chunk)
    ## Prepare the execution script to run the whole docking assay in MENDIETA
    prepare_docking_exec_script_hpc(docking_assay_folder,custom_params_string,hpc_assay_path,n_chunk)

#################

if __name__ == '__main__':

    print("This is a test")

    #create_docking_assay()
    
    #prepare_docking_exec_script("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/files_for_testing/docking_module/docking_assays/docking_assay_15","nada")