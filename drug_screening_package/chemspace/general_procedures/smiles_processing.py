import sqlite3
import pandas as pd
from chemplot import Plotter
from termcolor import colored
from pandarallel import pandarallel # import pandarallel
from rdkit.Chem import rdmolops
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit import RDLogger
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import matplotlib.pyplot as plt
from pathlib import Path
RDLogger.DisableLog('rdApp.*')
import sys
import pandas as pd
from tqdm import tqdm
MODULE_PATH = '/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/drug_screening_package'
sys.path.append(MODULE_PATH)
import chemspace.general_procedures.operations_on_sql_database as db_ops
import chemspace.input_output_operations.i_o_checkings as db_i_o

def generate_clean_smiles(smiles):
    
    try:
        sanitized_smiles = sanitize_smiles_errors(smiles)
        smiles_free_stereo = remove_smiles_stereo(sanitized_smiles)
        largest_mol_smiles = get_largest_fragment(smiles_free_stereo)
        # Add here potential processing steps for cleaning/sanitizing the SMILES notation 
        #####
        #
        clean_smiles = largest_mol_smiles # The last cleaning step should be assigned to 'clean_smiles'
        inchi_key = compute_inchi_key(largest_mol_smiles)
        return clean_smiles, inchi_key
    except:
        pass

def generate_clean_smiles_2(smiles):
    try: 
        sanitized_smiles = sanitize_smiles_errors(smiles)
        smiles_free_stereo = remove_smiles_stereo(sanitized_smiles)
        largest_mol_smiles = get_largest_fragment(smiles_free_stereo)
        # Add here potential processing steps for cleaning/sanitizing the SMILES notation 
        #####
        #
        clean_smiles = largest_mol_smiles # The last cleaning step should be assigned to 'clean_smiles'
        return clean_smiles
    
    except Exception as error:
        print(error)

def generate_clean_smiles_2_retain_stereo(smiles):
    sanitized_smiles = sanitize_smiles_errors(smiles)
    largest_mol_smiles = get_largest_fragment(sanitized_smiles)
    # Add here potential processing steps for cleaning/sanitizing the SMILES notation 
    #####
    #
    clean_smiles = largest_mol_smiles # The last cleaning step should be assigned to 'clean_smiles'
    return clean_smiles

def remove_smiles_stereo(smiles):
    #smiles_free_stereo = smiles.replace("@","").replace("[","").replace("]","").replace("H","")
    smiles_free_stereo = smiles.replace("@","")
    if "[CH]" in smiles_free_stereo:
        smiles_free_stereo = smiles_free_stereo.replace("[CH]","C")
    return smiles_free_stereo

def get_largest_fragment(smiles):   
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_frags = rdmolops.GetMolFrags(mol, asMols = True)
        largest_mol = max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())
        largest_mol_smiles = Chem.MolToSmiles(largest_mol)    
        return largest_mol_smiles
    except:
        pass

def compute_inchi_key(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        inchi_key = Chem.MolToInchiKey(mol)
        return inchi_key
    except:
        pass

def depict_smiles(smiles):
    try: 
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol)
        return img
    except:
        print(f"Error procesando: {smiles}")

def sanitize_smiles_errors(smiles):
    
    if "[nH]" in smiles:
        sanitized_smiles = smiles.replace("[nH]","[NH]")
        return sanitized_smiles
    else:
        return smiles
    
def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol    

def smiles_from_mol(mol):
    smiles = Chem.MolToSmiles(mol)
    return smiles

def reaction_from_smarts(smarts):
    reaction_mol = AllChem.ReactionFromSmarts(smarts)
    return reaction_mol

def react_1_component(reaction_mol,smiles):
    """
    This function will receive a reaction_mol (i.e. the smarts reaction already converted to the RDKit format), will apply the chemical transformation to the corresponding SMILES input and will afterwards convert the product to the SMILES format, which will be returned by the fucntion.
    --- Parameters ---
    - reaction_mol: the 1 component reaction to be applied by the function, which has already been transformed to the RDKit format.
    - smiles: the smiles notation of the reactant to be transformed.

    --- Returns ---
    - product_smiles: the SMILES notation of the transformed reactant.
    """
    mol = mol_from_smiles(smiles)
    try:
        product_mol = reaction_mol.RunReactants((mol,))[0][0]
        product_smiles = smiles_from_mol(product_mol)
        return product_smiles
    
    except:
        pass

def react_2_components(reaction_mol,smiles_1,smiles_2):
    """
    This function will receive a reaction_mol (i.e. the smarts reaction already converted to the RDKit format), will apply the chemical transformation to the corresponding SMILES input and will afterwards prepare the corresponding product to the SMILES format, which will be returned by the function.
    --- Parameters ---
    - reaction_mol: the 1 component reaction to be applied by the function, which has already been transformed to the RDKit format.
    - smiles_1: the smiles notation of the first reactant to be transformed.
    - smiles_2: the smiles notation of the first reactant to be transformed.

    --- Returns ---
    - product_smiles: the SMILES notation of the transformed reactant.
    """
    
    try:
        mol_1 = mol_from_smiles(smiles_1)
        mol_2 = mol_from_smiles(smiles_2)
        product_mol = reaction_mol.RunReactants((mol_1,mol_2))[0][0]
        product_smiles = smiles_from_mol(product_mol)
        return product_smiles
    
    except:
        pass

def apply_reaction_1_component(origin_db_name,reaction_table,reaction_name,reactants_table,destination_db_name,product_table_name,sql_action):
    """
    This function will receive a reaction instruction in the SMARTS format in which only a 1 component transformation is performed, and will a apply it to a corresponding SMILES set and will return a reaction products
    ------
    Parameters
    ------
    """    
    reaction_smarts = db_ops.read_database_record(origin_db_name,reaction_table,"SMARTS","Reaction",reaction_name)    
    reaction_mol = reaction_from_smarts(reaction_smarts)
    reactants_smiles_df = db_ops.read_database_field(origin_db_name,reactants_table,"SMILES")
    reactants_inchi_key_df = db_ops.read_database_field(origin_db_name,reactants_table,"inchi_key")
       
    reactant_inchi_key_list = []
    product_smiles_list = []
    product_inchi_key_list = []

    counter = 0
    for index, row in reactants_smiles_df.iterrows():
        smiles = row.iloc[0]
        mol = mol_from_smiles(smiles)
        product_mol = reaction_mol.RunReactants((mol,))[0][0]
        product_smiles = smiles_from_mol(product_mol)
        product_inchi_key = compute_inchi_key(product_smiles)
        product_smiles_list.append(product_smiles)
        product_inchi_key_list.append(product_inchi_key)
        reactant_inchi_key_list.append(reactants_inchi_key_df.iloc[counter].values[0]) # Extract the corresponding inchi_key from the dataframe
        counter+=1
        
    product_df = db_i_o.merge_3_lists_into_df(product_smiles_list,product_inchi_key_list,reactant_inchi_key_list,"SMILES","inchi_key","origin_inchi_key")
    
    db_ops.store_df_to_sql(destination_db_name,product_df,product_table_name,sql_action)

def apply_reaction_1_component_tqdm(origin_db_name,reaction_table,reaction_name,reactants_table,destination_db_name,product_table_name,sql_action):
    """
    This function will receive a reaction instruction in the SMARTS format in which only a 1 component transformation is performed, and will a apply it to a corresponding SMILES set and will return a reaction products. 
    
    In this particular case, the function will process the corresponding dataframe using the 'tqdm' library. 
    ------
    Parameters
    ------
    - origin_db_name: this is the full path to the db file containing the one component reaction to be applied.
    - reaction_table: the name of the table containing the SMARTS based one component reaction to be applied.
    - reaction_name: the identifier of the reaction to be applied, and that is present in the corresponding table.
    - reactants_table: the name of the table containing the molecules in SMILES format that will be subjected to a one-component transformation.
    - destination_db_name: this is the full path to the db file were the products are to be stored.
    - product_table_name: the name of the table in the destination database that will hold the corresponding products.
    - sql_action: the instruction ('append','replace') that will be applied to the destination database in case the product table already exists.

    ------
    Returns
    ------
    A table containing a set of one-component reaction table.
    """    
    # Read and prepare the reaction to apply
    reaction_smarts = db_ops.read_database_record(origin_db_name,reaction_table,"SMARTS","Reaction",reaction_name)    
    reaction_mol = reaction_from_smarts(reaction_smarts)
    chemicals_df = db_ops.read_database_table(origin_db_name,reactants_table)
    
    tqdm.pandas() # Instatiate the tqdm package for posterior running.
    chemicals_df['product_smiles'] = chemicals_df['SMILES'].progress_apply(lambda x: react_1_component(reaction_mol,x))
    chemicals_df['product_inchi_key'] = chemicals_df['product_smiles'].progress_apply(lambda x: compute_inchi_key(x))
    chemicals_df.drop(['SMILES'],axis=1,inplace=True) # Drop the reactant SMILES column
    chemicals_df.rename(columns={'product_smiles':"SMILES"},inplace=True)
    chemicals_df.rename(columns={'inchi_key':'inchi_key_reactant_1'},inplace=True)
    chemicals_df.rename(columns={'product_inchi_key':"inchi_key"},inplace=True)
    chemicals_df = chemicals_df[['SMILES','inchi_key','inchi_key_reactant_1']] # Reorder the dataframe columns
    
    # Store products to the database
    db_ops.store_df_to_sql(destination_db_name,chemicals_df,product_table_name,sql_action)

def apply_reaction_1_component_pandarallel(origin_db_name,reaction_name,reactants_table,destination_db_name,product_table_name,sql_action):
    """
    This function will receive a reaction instruction in the SMARTS format in which only a 1 component transformation is performed, and will a apply it to a corresponding SMILES set and will return a reaction products. 
    
    By default, the reaction SMARTS will be read from a preexisting table called 'reaction_types'

    The function will also compute the inchi_key of the product.

    In this particular case, the function will process the corresponding dataframe using the 'pandarallel' library. 
    ------
    Parameters
    ------
    - origin_db_name: this is the full path to the db file containing the one component reaction to be applied.
    - reaction_name: the identifier of the reaction to be applied, and that is present in the corresponding table.
    - reactants_table: the name of the table containing the molecules in SMILES format that will be subjected to a one-component transformation.
    - destination_db_name: this is the full path to the db file were the products are to be stored.
    - product_table_name: the name of the table in the destination database that will hold the corresponding products.
    - sql_action: the instruction ('append','replace') that will be applied to the destination database in case the product table already exists.

    ------
    Returns
    ------
    A table containing a set of one-component reaction table.
    """    
    # Read and prepare the reaction to apply
    reaction_smarts = db_ops.read_database_record(origin_db_name,'reaction_types',"SMARTS","Reaction",reaction_name)    
    reaction_mol = reaction_from_smarts(reaction_smarts)
    chemicals_df = db_ops.read_database_table(origin_db_name,reactants_table)
    
    # Initialize the pandarallel library for pandas
    pandarallel.initialize(progress_bar=True)
    
    chemicals_df['product_smiles'] = chemicals_df['SMILES'].parallel_apply(lambda x: react_1_component(reaction_mol,x))
    chemicals_df.drop(['SMILES'],axis=1,inplace=True) # Drop the reactant SMILES column
    chemicals_df.drop(['inchi_key'],axis=1,inplace=True) # Drop the reactant SMILES column
    chemicals_df['inchi_key'] = chemicals_df['product_smiles'].parallel_apply(lambda x: compute_inchi_key(x))
    chemicals_df.rename(columns={'product_smiles':"SMILES"},inplace=True)
    chemicals_df = chemicals_df[['SMILES','inchi_key']] # Reorder the dataframe columns
    
    # Create a mol_id index and make it the first column
    chemicals_df["mol_id"] = chemicals_df.index
    first_column = chemicals_df.pop('mol_id')
    chemicals_df.insert(0, 'mol_id', first_column)
    
    # Store products to the database
    db_ops.store_df_to_sql(destination_db_name,chemicals_df,product_table_name,sql_action)

def apply_reaction_2_components(origin_db_name,reaction_table,reaction_name,reactants1_table,reactants2_table,destination_db_name,product_table_name,sql_action,store_errors):
    """
    This function will receive a reaction instruction in the SMARTS format in which a 2 components transformation is performed, and will a apply it to a corresponding SMILES set and will return a reaction products. 
    
    
    The combinatorial reaction between reactant_1_table and reactant_2_table will be performed. In this case, the reactant_1 and reactant_2 need to match the order of reactants in which the reaction SMARTS is defined
    ------
    Parameters
    ------
    - origin_db_name: the full path of the database in which the reactants will be read from
    - reaction_table: the name of the table were the reaction to be applied will be read by the corresponding name.
    - reaction_name: the name of the 2 component reaction to be applied in order to obtain the corresponding product. This must exist in the 'reaction_table'.
    - reactants1_table: the name of the table containing the first reactants of the corresponding reaction. Pay attention to the order of reactants as defined in the reaction.
    - reactants2_table: the name of the table containing the second reactants of the corresponding reaction. Pay attention to the order of reactants as defined in the reaction.
    - destination_db_name: the full path of the database to which the corresponding products are to be written.
    - product_table_name: the name of the table were the products are to be written.
    - sql_action: action to apply if the products table already exists: options: 'REPLACE'; 'APPEND'"
    - store_errors: in case a reaction error occurs, store the involved SMILES in a separate table for further analysis and checkings. Options: '0','1'. 
    """    
    reaction_smarts = db_ops.read_database_record(origin_db_name,reaction_table,"SMARTS","Reaction",reaction_name)    
    reaction_mol = reaction_from_smarts(reaction_smarts)
    reactants_1_smiles_df = db_ops.read_database_field(origin_db_name,reactants1_table,"SMILES")
    reactants_2_smiles_df = db_ops.read_database_field(origin_db_name,reactants2_table,"SMILES")
       
    reactants_1_inchi_key_list = []
    reactants_2_inchi_key_list = []
    product_smiles_list = []
    product_inchi_key_list = []

    smiles_errors_reactant_1_list = []
    smiles_errors_reactant_2_list = []

    for index_1, row_1 in reactants_1_smiles_df.iterrows():
        smiles_1 = row_1.iloc[0]
        inchi_key_1 = compute_inchi_key(smiles_1)
    
        mol_1 = mol_from_smiles(smiles_1)
        for index, row_2 in reactants_2_smiles_df.iterrows():
            smiles_2 = row_2.iloc[0]
            inchi_key_2 = compute_inchi_key(smiles_2)
            try:
                mol_2 = mol_from_smiles(smiles_2)
                
                product_mol = reaction_mol.RunReactants((mol_1,mol_2))[0][0]
                product_smiles = smiles_from_mol(product_mol)
                product_inchi_key = compute_inchi_key(f"{product_smiles}")
                product_smiles_list.append(product_smiles)
                product_inchi_key_list.append(product_inchi_key)
                reactants_1_inchi_key_list.append(inchi_key_1) # Extract the corresponding inchi_key from the dataframe
                reactants_2_inchi_key_list.append(inchi_key_2) # Extract the corresponding inchi_key from the dataframe
            except Exception as error:
                print("ERROR")
                smiles_errors_reactant_1_list.append(smiles_1)
                smiles_errors_reactant_2_list.append(smiles_2)
                continue
      
    product_df = db_i_o.merge_4_lists_into_df(product_smiles_list,product_inchi_key_list,reactants_1_inchi_key_list,reactants_2_inchi_key_list,"SMILES","inchi_key","reactant_1_inchi_key","reactant_2_inchi_key")
    
    if store_errors == 1:
        if len(smiles_errors_reactant_1_list) > 0:
            errors_table = f'ERRORS_{product_table_name}'
            errors_df = db_i_o.merge_2_lists_into_df(smiles_errors_reactant_1_list,smiles_errors_reactant_2_list,"SMILES_1","SMILES_2")
            db_ops.store_df_to_sql(destination_db_name,errors_df,errors_table,sql_action)

    db_ops.store_df_to_sql(destination_db_name,product_df,product_table_name,sql_action)

def apply_reaction_2_components_pandarallel(origin_db_name,reaction_name,reactants1_table,reactants2_table,destination_db_name,product_table_name):
    """
    This function will receive a reaction instruction in the SMARTS format in which a 2 components transformation is performed.
    
    By default the reaction SMARTS will be read from a table called 'reaction_types'
    
    The reaction will be applied to a corresponding SMILES set and will return the reaction products in a column named 'SMARTS'. 
    
    The combinatorial reaction between reactant_1_table and reactant_2_table will be performed. 
    
    In this case, the reactant_1 and reactant_2 need to match the order of reactants in which the reaction SMARTS is defined

    Also, the inchi_key values for each specie (product, reactant_1 and reactant_2) will be computed, and stored together with a description of the reaction applied.
    ------
    Parameters
    ------
    - origin_db_name: the full path of the database in which the reactants will be read from
    - reaction_table: the name of the table were the reaction to be applied will be read by the corresponding name.
    - reaction_name: the name of the 2 component reaction to be applied in order to obtain the corresponding product. This must exist in the 'reaction_table'.
    - reactants1_table: the name of the table containing the first reactants of the corresponding reaction. Pay attention to the order of reactants as defined in the reaction.
    - reactants2_table: the name of the table containing the second reactants of the corresponding reaction. Pay attention to the order of reactants as defined in the reaction.
    - destination_db_name: the full path of the database to which the corresponding products are to be written.
    - product_table_name: the name of the table were the products are to be written.

    ------
    Returns
    ------ 
    A table containing all the generated information is saved into the destination database.
    """    
    reaction_smarts = db_ops.read_database_record(origin_db_name,'reaction_types',"SMARTS","Reaction",reaction_name)    
    reaction_mol = reaction_from_smarts(reaction_smarts)
    print
    reactants_1_smiles_df = db_ops.read_database_field(origin_db_name,reactants1_table,"SMILES")
    reactants_2_smiles_df = db_ops.read_database_field(origin_db_name,reactants2_table,"SMILES")
    df_products = pd.DataFrame() # Will store the all combinations of products
    df_products_iteration = pd.DataFrame() # Will temporary products through iteration
    # Determine which dataframe is larger to optimize pandarallel running
    reactants_1_length = len(reactants_1_smiles_df)
    reactants_2_length = len(reactants_2_smiles_df)
    
    # Initialize the pandarallel library for pandas
    pandarallel.initialize(progress_bar=True)

    try:
        if reactants_1_length > reactants_2_length:
            counter = 0
            for index, row in reactants_2_smiles_df.iterrows(): # the for loop is applied to the smaller dataframe (reactant_2_smiles)
                reactant_2_smiles = row.iloc[0]
                # Apply the 2 components reaction to the whole dataframe
                df_products_iteration["SMILES"] = reactants_1_smiles_df["SMILES"].parallel_apply(lambda reactant_1_smiles: react_2_components(reaction_mol,reactant_1_smiles,reactant_2_smiles))
                # Compute the inchi_key for the product
                df_products_iteration["inchi_key"] = df_products_iteration["SMILES"].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
                ## Create a joined registry for the 2 reactants
                #df_products_iteration["SMILES_reactants"] = reactants_1_smiles_df["SMILES"].parallel_apply(lambda reactant_1_smiles: reactant_1_smiles + ">>" + reactant_2_smiles + f">>Obtained through reaction: {reaction_name}")
                ## Split the registry of the reactants
                #df_products_iteration[["SMILES_reactant_1","SMILES_reactant_2","Description"]] = df_products_iteration["SMILES_reactants"].str.split(">>",expand=True)
                ## Drop the reactants merged column
                #df_products_iteration.drop(columns=["SMILES_reactants"],inplace=True)
                ## Compute reactants inchi_key_values
                #df_products_iteration["reactant_1_inchi_key"] = df_products_iteration["SMILES_reactant_1"].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
                #df_products_iteration["reactant_2_inchi_key"] = df_products_iteration["SMILES_reactant_2"].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
                ## Drop reactants SMILES columns
                #df_products_iteration.drop(columns=["SMILES_reactant_1","SMILES_reactant_2"],inplace=True)
                ## Sort df column order
                #df_products_iteration = df_products_iteration[["SMILES","inchi_key","reactant_1_inchi_key","reactant_2_inchi_key","Description"]]
                # Concatenate the overall dataframe for each iteration
                df_products = pd.concat([df_products,df_products_iteration],axis=0)
                print(f"Iteration {counter}")
                counter+=1

            # Create a mol_id index and make it the first column
            df_products["mol_id"] = df_products.index
            first_column = df_products.pop('mol_id')
            df_products.insert(0, 'mol_id', first_column)

            # Store the final dataframe to the corresponding database
            db_ops.store_df_to_sql(destination_db_name,df_products,product_table_name,"replace")
            print("Reaction finished SUCCESSFULLY")
    
        else: # the for loop is applied to the smaller dataframe (reactant_1_smiles)
            counter = 0
            for index, row in reactants_1_smiles_df.iterrows():
                reactant_1_smiles = row.iloc[0]
                # Apply the 2 components reaction to the whole dataframe
                df_products_iteration["SMILES"] = reactants_2_smiles_df["SMILES"].parallel_apply(lambda reactant_2_smiles: react_2_components(reaction_mol,reactant_1_smiles,reactant_2_smiles))
                # Compute the inchi_key for the product
                df_products_iteration["inchi_key"] = df_products_iteration["SMILES"].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
                ## Create a joined registry for the 2 reactants
                #df_products_iteration["SMILES_reactants"] = reactants_2_smiles_df["SMILES"].parallel_apply(lambda reactant_2_smiles: reactant_1_smiles + ">>" + reactant_2_smiles + f">>Obtained through reaction: {reaction_name}")
                ## Split the registry of the reactants
                #df_products_iteration[["SMILES_reactant_1","SMILES_reactant_2","Description"]] = df_products_iteration["SMILES_reactants"].str.split(">>",expand=True)
                ## Drop the reactants merged column
                #df_products_iteration.drop(columns=["SMILES_reactants"],inplace=True)
                ## Compute reactants inchi_key_values
                #df_products_iteration["reactant_1_inchi_key"] = df_products_iteration["SMILES_reactant_1"].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
                #df_products_iteration["reactant_2_inchi_key"] = df_products_iteration["SMILES_reactant_2"].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
                ## Drop reactants SMILES columns
                #df_products_iteration.drop(columns=["SMILES_reactant_1","SMILES_reactant_2"],inplace=True)
                ## Sort df column order
                #df_products_iteration = df_products_iteration[["SMILES","inchi_key","reactant_1_inchi_key","reactant_2_inchi_key","Description"]]
                # Concatenate the overall dataframe for each iteration
                df_products = pd.concat([df_products,df_products_iteration],axis=0)
                print(f"Iteration {counter}")
                counter+=1

            # Create a mol_id index and make it the first column
            df_products["mol_id"] = df_products.index
            first_column = df_products.pop('mol_id')
            df_products.insert(0, 'mol_id', first_column)

            # Store the final dataframe to the corresponding database
            db_ops.store_df_to_sql(destination_db_name,df_products,product_table_name,"replace")
            print("Reaction finished SUCCESSFULLY")

    except Exception as error:
        print(error)

def compute_chemspace_coord(db_name,table_name,calc_type):
    """
    Given a table containing a column named 'SMILES', this function will compute the x, y coordinates to plot the corresponding chemical space.

    The coordinate are computed using the umap algorithm as implemented in the library 'Chemplot'

    --- Parameters ---
    - db_name: the full path of the database were the table to be used for chemical space computacion is located
    - table_name: the table in which a SMILES colum will be retrieved in order to compute the chemical space
    - calc_type: the type of calculation used to compute the chemical space: i.e. 'structural' or 'tailored'

    --- Returns ---
    The corresponding x,y coordinates are writen to the given table and stored within the dataframe.
    """

    conn = sqlite3.connect(db_name)
    sql_action = f"""SELECT *
                    FROM {table_name};"""

    df = pd.read_sql_query(sql_action,conn)
    print("START computing the chemical space")
    
    # Check if the table already contains chemical space calculation data
    columns_list = list(df.columns)
    if "UMAP-1" in columns_list and "UMAP-2" in columns_list:
        print(colored(f"The chemspace on table {table_name} has already been computed. Exiting.","red"))
        exit()
    
    # Compute the chemical space and perform some processing
    cp = Plotter.from_smiles(df['SMILES'], target=df.index, target_type="C", sim_type=calc_type)
    chemical_space_data = cp.umap(random_state=500).drop(columns=['target'])
    df_with_chem_space = pd.concat([df,chemical_space_data],axis=1)
    df_with_chem_space["calc_type"] = calc_type

    # Store the calculation to the database
    db_ops.store_df_to_sql(db_name,df_with_chem_space,table_name,"replace")
    print("Succesfully computed and stored de chemical space info.")

def compute_molecule_properties(smiles):
    """
    This function will receive a SMILES string and will compute a set of hardcoded molecular properties.

    The function will return a dictionary containing key:value pairs for the corresponding computed properties.

    ------
    Parameters
    ------
    - smiles: the SMILES notation of the molecule whose properties are to be computed.

    ------
    Returns
    ------
    A dictionary containing the key:value pairs for the corresponding computed properties
    """
    molecular_properties_to_compute = ['MolWt', 'MolLogP', 'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds']
    integer_values = ['NumHAcceptors', 'NumHDonors', 'NumRotatableBonds']
    
    properties_dict = {}
    
    try:

        mol = Chem.MolFromSmiles(smiles)

        for nm,fn in Descriptors._descList:
            if nm in molecular_properties_to_compute:
                value = fn(mol)
                properties_dict[nm] = round(value,2)

        chiral_ctrs_nbr = count_chiral_centers(smiles)
        properties_dict["chiral_ctrs_nbr"] = int(chiral_ctrs_nbr)

        return properties_dict
    
    except:
        # Return a dictionary with values 999
        for property in molecular_properties_to_compute:
            properties_dict[property] = round(999,2)
        
        properties_dict["chiral_ctrs_nbr"] = int(999)

        # Return the fake dictionary if the SMILES fails
        return properties_dict
        
        #pass

def compute_molecule_properties_on_table(db_name,table_name):
    """
    This function will receive a table, and compute the molecular properties for each SMILES included in it.

    ------
    Parameters
    ------
    - db_name: the full path to the database in which the table to be processed exists.
    - table_name: the name of the table to be processed.
    - molecular_properties_to_compute: the list of molecular properties to compute.
    - integer_values: the propoerties that are going to be converted into integer values.

    ------
    Returns
    ------
    Columns containing the calculated molecular descriptors are added to the original table.
    """
    
    integer_values = ['NumHAcceptors', 'NumHDonors', 'NumRotatableBonds','chiral_ctrs_nbr']

    conn = sqlite3.connect(db_name)
    sql = f"""SELECT * 
               FROM {table_name};"""

    try:
        df = pd.read_sql_query(sql,conn)

        ## This will check if the table already contains the molecular descriptors. In that case, it will stop execution.
        for column in df.columns:
            if column in integer_values:
                print("The computation on the table has already been done in this table. Exiting")
                exit()

        # Initialize the pandarallel library for pandas
        pandarallel.initialize(progress_bar=True)
        df['properties_dict'] = df['SMILES'].parallel_apply(lambda x: compute_molecule_properties(x))
        df_props = pd.concat([df.drop(['properties_dict'],axis=1),df['properties_dict'].apply(pd.Series)],axis=1)

        # This will process the integer values columns
        for header in integer_values:
            df_props = df_props.astype({header:int})    

        ## Drop all rows containing a 999 values originated in a SMILES problem
        #df_props = df_props.drop(df_props[df_props["chiral_ctrs_nbr"] == 999].index)

        db_ops.store_df_to_sql(db_name,df_props,table_name,"replace")
        print(" Succesfully computed and stored de molecular properties info.")

    except Exception as error:
        print(error)
        print(f"Error retrieving the table {table_name}. Check for errors.")

def plot_histograms_of_properties(db_name,table_name,n_bins,miscellaneous_files):
    """
    Given a table in which molecular properties have been computed for the stored molecules, this function will plot a serie of histograms corresponding to each property present in the table.
    
    ------
    Parameters:
    ------
    - db_name: the full path to the database in which the table containing the properties is stored.
    - table_name: the name of the table containing the computed properties
    - n_bins: the number of bins included in the histograms
    - miscellaneous_files: the full path to the folder in which the Chemspace miscellaneous files are stored.
    
    ------
    Returns
    ------
    A serie of histograms stored in de /misc/{table_name} folder
    """

    conn = sqlite3.connect(db_name)
    sql = f"SELECT * FROM {table_name}"
    df = pd.read_sql_query(sql,conn)
    
    number_of_cpds = df.shape[0]
    
    exclude_columns_list = ["Pose_ID","mol_id","SMILES","inchi_key","available","LigName","sub_pose","receptor"]

    hist_output_dir = f'{miscellaneous_files}/{table_name}_props_hist'
    Path(hist_output_dir).mkdir(parents=True, exist_ok=True)

    for column in df.columns:
        #print(column)
        if column not in exclude_columns_list:
            try:
                counts, edges, bars = plt.hist(df[column],bins=n_bins,edgecolor='k')
                plt.bar_label(bars, rotation=90)
                plt.xlabel(column)
                plt.ylabel("Cpds number")
                plt.title(f"Number of cpds: {number_of_cpds}")
                plt.savefig(f"{hist_output_dir}/hist_{column}.png")
                plt.clf()
            except:
                continue

def return_enumerated_stereoisomers(smiles):
    """
    This function will compute the possible stereoisomers combinations given a SMILES string.

    The function will calculate and return a dictionary, containing ordered lists with the corresponding stereoisomer SMILES, overal stereoconfig and number of total stereoisomers for the specific compound.

    ------
    Parameters:
    ------
    - smiles: a SMILES string that is the input to the computation.
    """
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_hs = Chem.AddHs(mol)
        isomers = list(EnumerateStereoisomers(mol_hs))
        isomers_list = [] # The stereoisomer SMILES will be appended here.
        stereo_config_list_1 = [] # The stereo config will be appended here.

        # Compute each stereoisomer property
        for isomer in isomers:
            isomer_no_hs = Chem.RemoveHs(isomer)
            isomer_smiles = Chem.MolToSmiles(isomer_no_hs)
            stereo_config = Chem.FindMolChiralCenters(isomer,force=True,
            includeUnassigned=True,useLegacyImplementation=True)
            isomers_list.append(isomer_smiles)
            stereo_config_list_1.append(stereo_config)
   
        number_of_isomers = len(isomers_list)
        number_of_optical_isomers = count_chiral_centers(smiles)
        number_of_isomers_list = []

        # This will create a list with n-times the number of isomers for further exploding the list in a dataframe.
        for nbr in range(number_of_isomers):
            number_of_isomers_list.append(number_of_optical_isomers)

        # This will clean the stereo_config_list_1 to retain only the config in 'SS' (or whatever config) formate.
        stereo_config_list_2 = []
    
        for isomer_nbr in range(number_of_isomers):
            stereo_config = ''        
            chiral_centers_nbr = len(stereo_config_list_1[isomer_nbr])
            for chiral_center_nbr in range(chiral_centers_nbr):
                center_config = (stereo_config_list_1[isomer_nbr][chiral_center_nbr][1])    
                stereo_config = stereo_config + center_config

            # Construct the final config list.
            stereo_config_list_2.append(stereo_config)

        results_dict = {'isomers_smiles':isomers_list,'isomers_config':stereo_config_list_2,'chiral_ctrs_nbrs':number_of_isomers_list}

        return results_dict
    
    except: # Will proceed to the next function call 
        pass

def enumerate_stereoisomers(db_name,table_name):
    """
    This function will process a whole table in order to exhaustively enumerate the SMILES strings in a column.

    A new table named 'table_name_stereo_enum' will be created to contain the results.
    ------
    Parameters:
    ------
    - db_name: the full path to the database containing the table in which the SMILES strings are to be enumerated respect to stereoisomers.
    - table_name: the name of the table to be enumerated.
    
    ------
    Returns:
    ------
    A new table named 'table_name_stereo_enum'
    """

    conn = sqlite3.connect(db_name)

    sql = f"""SELECT SMILES
          FROM {table_name};"""

    try: 
        df = pd.read_sql_query(sql,conn)
        pandarallel.initialize(progress_bar=True)
        df["stereo_info"] = df["SMILES"].parallel_apply(lambda SMILES: return_enumerated_stereoisomers(SMILES))
        
        df_stereo = df['stereo_info'].apply(pd.Series)
        
        df_stereo_2 = df_stereo.explode(['isomers_smiles', 'isomers_config', 'chiral_ctrs_nbrs'],).reset_index().rename(columns={"isomers_smiles":"SMILES","isomers_config":"iso_config","chiral_ctrs_nbrs":"chrl_ctrs_nbrs"})
        
        df_stereo_2['inchi_key'] = df_stereo_2['SMILES'].parallel_apply(lambda SMILES: compute_inchi_key(SMILES))
        df_stereo_2 = df_stereo_2[["SMILES","inchi_key","iso_config","chrl_ctrs_nbrs"]]

        # Drop duplicated columns in case some kind of enumeration bug happens on the SMILES column
        df_stereo_2.drop_duplicates(subset=['inchi_key','iso_config','chrl_ctrs_nbrs'],inplace=True)
        # Create a mol_id index and make it the first column
        df_stereo_2["mol_id"] = df_stereo_2.index
        first_column = df_stereo_2.pop('mol_id')
        df_stereo_2.insert(0, 'mol_id', first_column) 

        # Store the final dataframe with the stereo info to the database.
        db_ops.store_df_to_sql(db_name,df_stereo_2,table_name+'_stereo_enum',"replace")
        print(colored(" The enumerated database was written SUCCESSFULLY.", "green"))
        
    except Exception as error:
        print(error)
        print("Error")

def count_chiral_centers(smiles):
    """
    This function will count and return the number of chiral centers for a given SMILES string:

    ------
    Parameters:
    ------
    - smiles: the smiles strings of the molecule to be processed
    
    ------
    Returns:
    ------
    An integer value corresponding to the number of stereocenters.
    """

    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_hs = Chem.AddHs(mol)
        # This will compute all stereoisomers: optical plus spatial
        isomers = list(EnumerateStereoisomers(mol_hs))

        # This will isolate the count of optical isomers:
        for isomer in isomers:
            stereo_config = Chem.FindMolChiralCenters(isomer,force=True,includeUnassigned=True,useLegacyImplementation=True)
            chiral_centers_nbr = len(stereo_config) 

        return chiral_centers_nbr
    
    except Exception as error:
        print(error)
        print(colored(f"Error computing the number of stereocenters for {smiles}.","red"))

def enumerate_stereoisomers_selected_molecule(db_name_origin,db_name_dest,origin_table,dest_table,row_id):
    """
    This function will enumerate all stereoisomers of a single selected molecule by using the corresponding 'row_id'.

    ------
    Parameters:
    ------
    - db_name_origin: the full path to the database were the table containing the molecule to be enumerated is present.
    - db_name_dest: the full path to the database were the the enumerated molecule is to be written.
    - origin_table: the table in which the 'row_id' value is to be searched.
    - dest_table: the table in which the enumerated enantiomers are to be stored.
    - row_id: the index value of the compound to enumerate.
    """

    try:
        conn = sqlite3.connect(db_name_origin)
        sql = f"""SELECT rowid, SMILES
                 FROM {origin_table}
                WHERE rowid = {row_id};"""
    
        df = pd.read_sql_query(sql,conn)
        df["stereo_info"] = df["SMILES"].apply(lambda SMILES: return_enumerated_stereoisomers(SMILES))
        df_stereo = df['stereo_info'].apply(pd.Series)
        df_stereo_2 = df_stereo.explode(['isomers_smiles', 'isomers_config', 'isomers_nbrs'],).reset_index().rename(columns={"isomers_smiles":"SMILES","isomers_config":"iso_config","isomers_nbrs":"iso_nbr"})
        df_stereo_2['inchi_key'] = df_stereo_2['SMILES'].apply(lambda SMILES: compute_inchi_key(SMILES))
        df_stereo_2 = df_stereo_2[["SMILES","inchi_key","iso_config","iso_nbr"]]

        # Store the final dataframe with the stereo info to the database.
        db_ops.store_df_to_sql(db_name_dest,df_stereo_2,dest_table,"replace")
        print(colored(" The enumerated database was written SUCCESSFULLY.", "green"))

    except Exception as error:
        print(error)
        print(colored(f"Error enumerating the stereoisomers of ligand in row: {row_id}","red"))

##########################################
## Generate the tests for the functions ##
##########################################

if __name__ == "__main__":
    
    print("This is a test")

    #apply_reaction("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reaction_types","azide_from_amine","alpha_aminoacid","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","azides","replace")
    
    #apply_reaction_1_component("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reaction_types","azide_from_amine","alpha_aminoacid","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","A","append")
    
    #apply_reaction_1_component_tqdm("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reaction_types","azide_from_amine","alpha_aminoacid","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","azides","replace")

    #apply_reaction_1_component_pandarallel("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/sample.db","azide_from_amine","alpha_aminoacid","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/sample.db","azides","replace")

    #apply_reaction_2_components("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","reaction_types","acylation","alpha_aminoacid","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/sample.db","alkyl_alpha_aminoacid","replace",1)

    #apply_reaction_2_components_pandarallel("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","reaction_types","acylation","alpha_aminoacid","alcohols","/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","alkyl_alpha_aminoacid")

    #compute_chemspace_coord("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_3.db","alpha_aminoacid")    

    #compute_molecule_properties_on_table("/home/fredy/MisDocumentos/Diseno-de-Scripts/Drug_Screening_Platform/processed_data/Emolecules_reagents_4.db","alpha_aminoacid")

    #return_enumerated_stereoisomers("CC(O)CC(Br)CC(O)(F)C")

    #enumerate_stereoisomers("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/processed_data/Emolecules_reagents_4.db","alpha_aminoacid")

    #enumerate_stereoisomers_selected_molecule("/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/chemspace/processed_data/MPRO_binders.db","/home/fredy/MisDocumentos/Diseno-de-Scripts/CADD_platform/TEST/chemspace/processed_data/MPRO_binders.db","mpro_binders_selected","nirmatrelvir_enumerated_stereo",2)    

