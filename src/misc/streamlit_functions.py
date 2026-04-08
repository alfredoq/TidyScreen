from io import StringIO
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

def read_database_as_dataframe(db_path, table_name):
    """
    Read a table from a SQLite database and return it as a pandas DataFrame.

    Args:
        db_path (str): Path to the SQLite database file.
        table_name (str): Name of the table to read.

    Returns:
        pd.DataFrame: DataFrame containing the table data.
    """
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
        conn.close()
        return df
    except Exception as e:
        print(f"❌ Error reading database: {e}")
        return pd.DataFrame()
    
def get_tables_info(db_path):
    """
    Return a DataFrame with all table names and number of rows for a given SQLite database.

    Args:
        db_path (str): Path to the SQLite database.

    Returns:
        pd.DataFrame: DataFrame with columns ['table', 'rows'].
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        # Get all table names
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]
        info = []
        for table in tables:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            nrows = cursor.fetchone()[0]
            info.append({'table': table, 'rows': nrows})
        conn.close()
        return pd.DataFrame(info)
    except Exception as e:
        return pd.DataFrame([{'table': 'Error', 'rows': str(e)}])
    
def get_docking_assay_registers(db_path):
    """
    Retrieve columns 'assay_name', 'table_name', 'docking_method_name', 'compound_count', and 'notes'
    from the 'docking_assays' table in the given SQLite database.

    Args:
        db_path (str): Path to the SQLite database.

    Returns:
        pd.DataFrame: DataFrame with the specified columns, or empty DataFrame if not found.
    """
    columns = ['assay_name', 'table_name', 'docking_method_name', 'compound_count', 'notes', 'receptor_info']
    try:
        conn = sqlite3.connect(db_path)
        query = f"SELECT {', '.join(columns)} FROM docking_assays"
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df
    except Exception as e:
        return None
    
def get_docking_results(db_path):
    """
    Retrieve columns 'LigName', 'pose_rank', 'run_number', 'docking_score', and 'cluster_size'
    from the 'Results' table in the given SQLite database.

    Args:
        db_path (str): Path to the SQLite database.

    Returns:
        pd.DataFrame: DataFrame with the specified columns, or empty DataFrame if not found.
    """
    columns = ['LigName', 'pose_rank', 'run_number', 'docking_score', 'cluster_size']
    try:
        conn = sqlite3.connect(db_path)
        query = f"SELECT {', '.join(columns)} FROM Results"
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df
    except Exception as e:
        
        return None

def get_mmpbsa_results(db_path):
    
    # check if the MMPBSA_results table exists
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='processed_mmgbsa_decomposition_json';")
        if cursor.fetchone() is None:
            return None  # Table does not exist
        
        query = "SELECT pose_id FROM processed_mmgbsa_decomposition_json"
                
        df = pd.read_sql_query(query, conn)
        
        conn.close()
        return df
    
    except Exception as e:
        return None

def get_mmpbsa_data_for_pose(db_path, pose_id):
    
    # check if the MMPBSA_results table exists
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='processed_mmgbsa_decomposition_json';")
        if cursor.fetchone() is None:
            return None  # Table does not exist
        
        # Retrieve all processed MMGBSA JSON entries
        cursor.execute('''
            SELECT p.pose_id, r.LigName, p.data
            FROM processed_mmgbsa_decomposition_json p
            JOIN Results r ON p.pose_id = r.Pose_ID WHERE p.pose_id = ?
        ''', (pose_id,))
        
        records = cursor.fetchall()
        
        if not records:
            return None  # No data found for the given pose_id
        # Convert the retrieved data into a DataFrame
        # Reconstruct DataFrame from JSON
        
        for pose_id, ligname, mmgbsa_json in records:
        
            mmgbsa_df = pd.read_json(StringIO(mmgbsa_json), orient='split')
            
        conn.close()
        return mmgbsa_df
        
        
    
    except Exception as e:
        return None
    
    

def construct_hist_for_ligand(df, ligand_name):
    """
    Filter the dataframe for the given ligand name in the 'LigName' column and plot a barplot of 'pose_rank' vs 'cluster_size'.

    Args:
        df (pd.DataFrame): The input DataFrame containing at least 'LigName', 'pose_rank', and 'cluster_size' columns.
        ligand_name (str): The ligand name to filter by.

    Returns:
        matplotlib.figure.Figure: The matplotlib Figure object containing the barplot.
    """
    if df is None or df.empty or 'LigName' not in df.columns or 'cluster_size' not in df.columns or 'pose_rank' not in df.columns:
        raise ValueError("DataFrame must contain 'LigName', 'pose_rank', and 'cluster_size' columns and not be empty.")
    filtered = df[df['LigName'] == ligand_name]
    if filtered.empty:
        raise ValueError(f"No rows found for ligand: {ligand_name}")
    fig, ax = plt.subplots()
    ax.bar(filtered['pose_rank'], filtered['cluster_size'], color='skyblue', edgecolor='black')
    ax.set_title(f"{ligand_name}: Pose Rank vs Cluster Size", size=8)
    ax.set_xlabel('Pose Rank', size=8)
    ax.set_ylabel('Cluster Size', size=8)
    
    fig.set_size_inches(3, 3)  # width, height in inches
    
    return fig
    
    
def get_pdbs_info(pdbs_db_path):
    """
    Retrieve all entries from the 'pdbs' table in the given SQLite database.

    Args:
        pdbs_db_path (str): Path to the SQLite database containing the 'pdbs' table.    
    Returns:
        pd.DataFrame: DataFrame with all columns from the 'pdbs' table,
        
    """
    
    try:
        conn = sqlite3.connect(pdbs_db_path)
        query_pdb_models = "SELECT file_id, pdb_model_name, project_name, original_path, filename, description, notes FROM pdb_models"
        df_pdb_models = pd.read_sql_query(query_pdb_models, conn)
        conn.close()
        return df_pdb_models
    
    except Exception as e:
        return None
    
    
def get_pdb_templates_info(pdbs_db_path):
    
    try:
        conn = sqlite3.connect(pdbs_db_path)
        query_templates_models = "SELECT pdb_id, pdb_template_name, pdb_model_name, project_name, template_folder_path, pdb_analysis, his_names, has_ligands, renumbering_dict, notes FROM pdb_templates"
        df_templates_models = pd.read_sql_query(query_templates_models, conn)
        conn.close()
        return df_templates_models
    
    except Exception as e:
        return None
    
    
def get_docking_models_info(receptors_db_path):
    
    try:
        conn = sqlite3.connect(receptors_db_path)
        query_docking_models = "SELECT id, pdb_id, receptor_model_name, template_name, pdb_model_name, pdb_to_convert, configs, notes FROM receptor_models"
        df_receptors_models = pd.read_sql_query(query_docking_models, conn)
        conn.close()
        return df_receptors_models
    
    except Exception as e:
        return None