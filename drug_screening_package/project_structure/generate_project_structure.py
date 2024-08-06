import os
from termcolor import colored
from pathlib2 import Path
import pathlib
from glob import glob

def configure_project_structure(project_base_dir,project_name):
    """
    This function will generate all the folder structure required to run a whole VS capaign.

    The function will also export all environment variables as required by the rest of functions of the CADD package

    The folling folders structure are created:

        - '/project_name
            |- chemspace
                |- raw_data
                |- processed_data
                |- misc
            |- docking
                |- receptor_models_registers
                |- docking_assays_registers
                |- docking_params_registers
                |- docking_assays
                |- raw_data : contains receptor preparation folders aimed to be stored in the registry database
            |- ml
                |- ml_model_devel

    ------
    Parameters:
    ------
    - project_base_dir: this is the full path in which the project folder will be created.
    - project_name: the name of the folder containing the project and corresponding substructures.
    ------
    Returns:
    ------
    The folder/subfoldes structure is created if not existing.
    """

    if os.path.isdir(f'{project_base_dir}/{project_name}'):
        print(colored(f"The project named {project_name} already exists. Exporting environment variables","green"))
    else:
        # Directories corresponding to ChemSpace
        os.makedirs(f'{project_base_dir}/{project_name}/chemspace/raw_data')
        os.makedirs(f'{project_base_dir}/{project_name}/chemspace/processed_data')
        os.makedirs(f'{project_base_dir}/{project_name}/chemspace/misc')
        # Directories corresponding to Docking
        os.makedirs(f'{project_base_dir}/{project_name}/docking/receptor_models_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/docking/docking_assays_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/docking/docking_params_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/docking/docking_assays')
        os.makedirs(f'{project_base_dir}/{project_name}/docking/raw_data')
        os.makedirs(f'{project_base_dir}/{project_name}/docking/misc')
        # Directories corresponding to Machine Learning
        os.makedirs(f'{project_base_dir}/{project_name}/ml/ml_model_devel')
        os.makedirs(f'{project_base_dir}/{project_name}/ml/misc/positive_binders')
        os.makedirs(f'{project_base_dir}/{project_name}/ml/misc/negative_binders')
        os.makedirs(f'{project_base_dir}/{project_name}/ml/misc/receptors')
        os.makedirs(f'{project_base_dir}/{project_name}/ml/misc/fingerprint_csv_files')
        os.makedirs(f'{project_base_dir}/{project_name}/ml/misc/pickle_model_files')
        # Directories corresponding to Molecular dynamics
        os.makedirs(f'{project_base_dir}/{project_name}/mds/mds_assays')
        os.makedirs(f'{project_base_dir}/{project_name}/mds/mds_assays_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/mds/mds_params_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/mds/misc')
        # Directories corresponding to QM/MM Molecular dynamics
        os.makedirs(f'{project_base_dir}/{project_name}/qm-mm-mds/qm-mm-mds_assays')
        os.makedirs(f'{project_base_dir}/{project_name}/qm-mm-mds/qm-mm-mds_assays_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/qm-mm-mds/qm-mm-mds_params_registers')
        os.makedirs(f'{project_base_dir}/{project_name}/qm-mm-mds/misc')

    ### This will export the required global environment variables
    ## Variables affecting ChemSpace
    raw_data_path = f"{project_base_dir}/{project_name}/chemspace/raw_data"
    main_db_output_path = f'{project_base_dir}/{project_name}/chemspace/processed_data'
    miscellaneous_files = f'{project_base_dir}/{project_name}/chemspace/misc'
    ## Variables affecting Docking
    receptor_models_registry_path = f'{project_base_dir}/{project_name}/docking/receptor_models_registers'
    docking_assays_registry_path = f'{project_base_dir}/{project_name}/docking/docking_assays_registers'
    docking_params_registry_path = f'{project_base_dir}/{project_name}/docking/docking_params_registers'
    docking_assays_storing_path = f'{project_base_dir}/{project_name}/docking/docking_assays'
    docking_raw_data_path = f'{project_base_dir}/{project_name}/docking/raw_data'
    docking_misc_data_path = f'{project_base_dir}/{project_name}/docking/misc'
    ## Variables affecting Machine Learning
    ml_models_path = f'{project_base_dir}/{project_name}/ml'
    ## Variables affecting Molecular Dynamics
    mds_path = f'{project_base_dir}/{project_name}/mds'
    ## Variables affecting QM/MM-Molecular Dynamics
    qm_mm_mds_path = f'{project_base_dir}/{project_name}/qm-mm-mds'

    return raw_data_path, main_db_output_path, miscellaneous_files, receptor_models_registry_path, docking_assays_registry_path, docking_params_registry_path, docking_assays_storing_path, docking_raw_data_path,docking_misc_data_path,ml_models_path,mds_path,qm_mm_mds_path

def set_modules_path(path):
    """
    This function will update the path of the package to the corresponding paths.

    ------
    Parameters:
    ------
    - path: the full path to the place in which the module exists
    """

    list_of_files = glob(f'{path}/**/*.py', recursive=True)

    for file in list_of_files:

        # The following if statement will prevent the writing of the patching script
        if "generate_project_structure.py" in file:
            continue    

        else: 

            try: 
            
                # This will find the status of the line to replace using as marker 'MODULE_PATH = ' to point to the variable assignment
                with open(file,'r') as file_to_check:
                    for line in file_to_check:
                        if 'MODULE_PATH = ' in line:
                            line_to_replace = line.rstrip()
            
                # Now execute the replacement of the corresponding line
                file = Path(rf"{file}")
                data = file.read_text()
                data = data.replace(line_to_replace,f"MODULE_PATH = '{path}'")
                file.write_text(data)

            except:
                continue

