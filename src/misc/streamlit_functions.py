from io import StringIO, BytesIO
import sqlite3
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

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
    
def get_docking_methods(db_path):
    """
    Retrieve all entries from the 'docking_methods' table.

    Args:
        db_path (str): Path to the docking_methods.db SQLite database.

    Returns:
        pd.DataFrame: DataFrame with all docking method records, or None on error.
    """
    columns = ['id', 'method_name', 'docking_engine', 'description', 'parameters', 'ligand_prep_params', 'created_date']
    try:
        conn = sqlite3.connect(db_path)
        query = f"SELECT {', '.join(columns)} FROM docking_methods"
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df
    except Exception as e:
        return None

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
    Also includes 'mmgbsa_total_energy' and 'mmgbsa_gas_energy' when present.

    Args:
        db_path (str): Path to the SQLite database.

    Returns:
        pd.DataFrame: DataFrame with the specified columns, or empty DataFrame if not found.
    """
    base_columns = ['Pose_ID', 'LigName', 'pose_rank', 'run_number', 'docking_score', 'cluster_size']
    optional_columns = ['mmgbsa_total_energy', 'mmgbsa_gas_energy']
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(Results)")
        existing = {row[1] for row in cursor.fetchall()}
        columns = base_columns + [c for c in optional_columns if c in existing]
        df = pd.read_sql_query(f"SELECT {', '.join(columns)} FROM Results", conn)
        conn.close()
        return df
    except Exception as e:
        return None

def get_extracted_poses_info(results_db_path):
    """
    Check for extracted docked poses in the four known output subdirectories adjacent to the results DB.

    Always returns all four entries. Each entry has:
      - 'directory': folder name
      - 'count': number of PDB files (0 if folder absent or empty)
      - 'path': full path to the folder
      - 'active': True if the folder exists and contains at least one PDB file
    """
    import os
    import glob

    results_dir = os.path.dirname(os.path.abspath(results_db_path))
    pose_dirs = ["most_stable_poses", "most_populated_poses", "most_populated_and_stable_poses", "all_poses"]
    entries = []
    for name in pose_dirs:
        full_path = os.path.join(results_dir, name)
        if os.path.isdir(full_path):
            pdb_files = glob.glob(os.path.join(full_path, "*.pdb"))
            count = len(pdb_files)
            active = count > 0
        else:
            count = 0
            active = False
        entries.append({"directory": name, "count": count, "path": full_path, "active": active})
    return entries

def get_prolif_all_column_names(db_path, table_name):
    """
    Return the sorted union of all column names across every pose stored in table_name.
    This ensures filter checkboxes cover interaction types that may be absent in a single pose.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(f"SELECT data FROM {table_name}")
        rows = cursor.fetchall()
        conn.close()
        all_cols = set()
        for (data_json,) in rows:
            df = pd.read_json(data_json, orient='split')
            all_cols.update(df.columns.tolist())
        return sorted(all_cols)
    except Exception:
        return []

def get_all_prolif_fingerprints(db_path, table_name):
    """
    Retrieve and concatenate ProLIF fingerprint DataFrames for all poses in the given table.
    A 'pose_id' column is prepended to identify each row's origin.
    Returns a combined DataFrame or None if the table is empty / on error.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(f"SELECT pose_id, data FROM {table_name} ORDER BY pose_id")
        rows = cursor.fetchall()
        conn.close()
        if not rows:
            return None
        frames = []
        for pose_id, data_json in rows:
            df = pd.read_json(data_json, orient='split')
            df.insert(0, 'pose_id', pose_id)
            frames.append(df)
        return pd.concat(frames, ignore_index=True)
    except Exception:
        return None

def get_prolif_tables(db_path):
    """
    Return a list of ProLIF fingerprint table names present in the results database.
    Tables follow the pattern 'processed_prolif_fps_json_condition_<N>'.
    Returns an empty list if none are found.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'processed_prolif_fps_json_condition_%';")
        tables = [row[0] for row in cursor.fetchall()]
        conn.close()
        return tables
    except Exception:
        return []

def get_prolif_poses(db_path, table_name):
    """
    Return a DataFrame with pose_id and computation_mode for all entries in the given ProLIF table.
    """
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(f"SELECT pose_id, computation_mode FROM {table_name} ORDER BY pose_id", conn)
        conn.close()
        return df
    except Exception:
        return None

def get_prolif_fingerprint_for_pose(db_path, table_name, pose_id):
    """
    Retrieve and reconstruct the ProLIF fingerprint DataFrame for a given pose_id.
    Returns a DataFrame or None if not found.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(f"SELECT data FROM {table_name} WHERE pose_id = ?", (pose_id,))
        row = cursor.fetchone()
        conn.close()
        if row is None:
            return None
        return pd.read_json(row[0], orient='split')
    except Exception:
        return None

def get_pose_ids_for_directory(db_path, directory_name):
    """
    Return the list of Pose_IDs from the Results table that correspond to the given
    extracted-poses directory (mirrors the logic used when extracting PDB files).
    """
    _queries = {
        "most_stable_poses": """
            SELECT Pose_ID FROM Results AS R1
            WHERE docking_score = (
                SELECT MIN(docking_score) FROM Results AS R2
                WHERE R1.LigName = R2.LigName
            )
        """,
        "most_populated_poses": """
            SELECT Pose_ID FROM Results AS R1
            WHERE cluster_size = (
                SELECT MAX(cluster_size) FROM Results AS R2
                WHERE R1.LigName = R2.LigName
            )
        """,
        "most_populated_and_stable_poses": """
            SELECT Pose_ID FROM Results AS R1
            WHERE cluster_size = (
                SELECT MAX(cluster_size) FROM Results AS R2
                WHERE R1.LigName = R2.LigName
            ) AND docking_score = (
                SELECT MIN(docking_score) FROM Results AS R3
                WHERE R1.LigName = R3.LigName
            )
        """,
        "all_poses": "SELECT Pose_ID FROM Results",
    }
    query = _queries.get(directory_name)
    if query is None:
        return []
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(query)
        ids = [row[0] for row in cursor.fetchall()]
        conn.close()
        return ids
    except Exception:
        return []

def get_filename_to_pose_id_map(db_path):
    """
    Return a dict mapping PDB filename -> pose info dict for every row in Results.
    Each value has keys 'pose_id', 'ligname', and 'pose_rank'.
    Filenames follow the convention '{LigName}_{run_number}.pdb' used during pose extraction.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT Pose_ID, LigName, run_number, pose_rank FROM Results")
        rows = cursor.fetchall()
        conn.close()
        return {
            f"{ligname}_{run_number}.pdb": {"pose_id": pose_id, "ligname": ligname, "pose_rank": pose_rank}
            for pose_id, ligname, run_number, pose_rank in rows
        }
    except Exception:
        return {}

def get_pose_labels_for_pose_ids(db_path, pose_ids):
    """
    Return a dict mapping pose_id -> 'LigName_poserank' for the given pose_ids,
    looked up from the Results table.
    Falls back to str(pose_id) for any pose_id not found.
    """
    if not pose_ids:
        return {}
    try:
        placeholders = ",".join("?" * len(pose_ids))
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            f"SELECT Pose_ID, LigName, pose_rank FROM Results WHERE Pose_ID IN ({placeholders})",
            pose_ids
        )
        rows = cursor.fetchall()
        conn.close()
        label_map = {pid: str(pid) for pid in pose_ids}
        for pose_id, ligname, pose_rank in rows:
            label_map[pose_id] = f"{ligname}_{pose_rank}"
        return label_map
    except Exception:
        return {pid: str(pid) for pid in pose_ids}

def get_prolif_tables_by_pose_ids(db_path, pose_ids):
    """
    Return ProLIF table names that contain at least one entry whose pose_id is in pose_ids.
    """
    if not pose_ids:
        return []
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'processed_prolif_fps_json_condition_%';")
        all_tables = [row[0] for row in cursor.fetchall()]
        placeholders = ",".join("?" * len(pose_ids))
        matching = []
        for table in all_tables:
            cursor.execute(f"SELECT COUNT(*) FROM {table} WHERE pose_id IN ({placeholders})", pose_ids)
            if cursor.fetchone()[0] > 0:
                matching.append(table)
        conn.close()
        return matching
    except Exception:
        return []

def get_prolif_poses_by_pose_ids(db_path, table_name, pose_ids):
    """
    Return a DataFrame with pose_id for entries in table_name whose pose_id is in pose_ids.
    """
    if not pose_ids:
        return None
    try:
        placeholders = ",".join("?" * len(pose_ids))
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(
            f"SELECT pose_id FROM {table_name} WHERE pose_id IN ({placeholders}) ORDER BY pose_id",
            conn, params=pose_ids
        )
        conn.close()
        return df
    except Exception:
        return None

def get_prolif_all_column_names_by_pose_ids(db_path, table_name, pose_ids):
    """
    Return the sorted union of all column names for entries matching the given pose_ids.
    """
    if not pose_ids:
        return []
    try:
        placeholders = ",".join("?" * len(pose_ids))
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(f"SELECT data FROM {table_name} WHERE pose_id IN ({placeholders})", pose_ids)
        rows = cursor.fetchall()
        conn.close()
        all_cols = set()
        for (data_json,) in rows:
            df = pd.read_json(data_json, orient='split')
            all_cols.update(df.columns.tolist())
        return sorted(all_cols)
    except Exception:
        return []

def get_all_prolif_fingerprints_by_pose_ids(db_path, table_name, pose_ids):
    """
    Retrieve and concatenate ProLIF fingerprint DataFrames for the given pose_ids.
    A 'pose_id' column is prepended to identify each row's origin.
    Returns a combined DataFrame or None if no matching entries / on error.
    """
    if not pose_ids:
        return None
    try:
        placeholders = ",".join("?" * len(pose_ids))
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            f"SELECT pose_id, data FROM {table_name} WHERE pose_id IN ({placeholders}) ORDER BY pose_id",
            pose_ids
        )
        rows = cursor.fetchall()
        conn.close()
        if not rows:
            return None
        frames = []
        for pose_id, data_json in rows:
            df = pd.read_json(data_json, orient='split')
            df.insert(0, 'pose_id', pose_id)
            frames.append(df)
        return pd.concat(frames, ignore_index=True)
    except Exception:
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
    
    
def get_electrostatic_profile_for_pose(mmpbsa_data, threshold=-1.0):
    
    df_electrostatic = mmpbsa_data[["residue", "ele"]].copy()
    df_electrostatic = df_electrostatic[df_electrostatic["ele"] < threshold]
    
    return df_electrostatic
    

def get_vdw_profile_for_pose(mmpbsa_data, threshold=-1.0):
    
    df_vdw = mmpbsa_data[["residue", "vdw"]].copy()
    df_vdw = df_vdw[df_vdw["vdw"] < threshold]
    
    return df_vdw

def get_polsolv_profile_for_pose(mmpbsa_data, threshold=-1.0):
    
    df_polsolv = mmpbsa_data[["residue", "polar_solvation"]].copy()
    df_polsolv = df_polsolv[df_polsolv["polar_solvation"] < threshold]
    
    return df_polsolv

def get_nonpolsolv_profile_for_pose(mmpbsa_data, threshold=-1.0):
    
    df_nonpolsolv = mmpbsa_data[["residue", "nonpolar_solvation"]].copy()
    df_nonpolsolv = df_nonpolsolv[df_nonpolsolv["nonpolar_solvation"] < threshold]
    
    return df_nonpolsolv

def get_gas_profile_for_pose(mmpbsa_data, threshold=-1.0):
    
    df_gas = mmpbsa_data[["residue", "gas"]].copy()
    df_gas = df_gas[df_gas["gas"] < threshold]
    
    return df_gas

def get_total_profile_for_pose(mmpbsa_data, threshold=-1.0):
    
    df_total = mmpbsa_data[["residue", "total"]].copy()
    df_total = df_total[df_total["total"] < threshold]
    
    return df_total

def create_mmpbsa_component_plot(df, x_col, y_col, title, xlabel, ylabel):
    """
    Create a lollipop plot from the given DataFrame.
    """
    fig, ax = plt.subplots()
    
    # Create vertical lines (stems)
    ax.vlines(df[x_col], 0, df[y_col], colors='skyblue', linewidth=1)
    
    # Create dots (lollipop heads)
    ax.scatter(df[x_col], df[y_col], color='darkblue', s=15, zorder=3)
    
    ax.set_title(title, size=8)
    ax.set_xlabel(xlabel, size=8)
    ax.set_ylabel(ylabel, size=8)
    # set the size of the x-axis labels to 8
    ax.tick_params(axis='x', labelsize=4, rotation=90)
    ax.tick_params(axis='y', labelsize=4, rotation=0)
    ax.grid(axis='y', alpha=0.3)

    fig.set_size_inches(6, 3)
    return fig


def save_positive_binder(project_path: str, assay_name: str, pose_file: str, directory: str, pose_full_path: str) -> str:
    """
    Save a pose flagged as a positive binder to the positive_binders.db database.

    The database is created at <project_path>/ml/training_sets/positive_binders.db.
    Returns a status string: 'saved', 'duplicate', 'conflict', or 'error:<message>'.
    'conflict' means the pose is already registered as a negative binder.

    Args:
        project_path (str): Root path of the active project.
        assay_name (str): Name of the docking assay.
        pose_file (str): PDB filename of the pose.
        directory (str): Pose folder name (e.g. 'most_stable_poses').
        pose_full_path (str): Absolute path to the PDB file.
    """
    try:
        db_dir = os.path.join(project_path, "ml", "training_sets")
        os.makedirs(db_dir, exist_ok=True)

        neg_db_path = os.path.join(db_dir, "negative_binders.db")
        if os.path.exists(neg_db_path):
            neg_conn = sqlite3.connect(neg_db_path)
            neg_conn.execute("""
                CREATE TABLE IF NOT EXISTS negative_binders (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    assay_name TEXT NOT NULL,
                    pose_file TEXT NOT NULL,
                    directory TEXT NOT NULL,
                    pose_full_path TEXT NOT NULL,
                    flagged_at TEXT NOT NULL,
                    UNIQUE(assay_name, pose_file, directory)
                )
            """)
            row = neg_conn.execute(
                "SELECT 1 FROM negative_binders WHERE assay_name=? AND pose_file=?",
                (assay_name, pose_file)
            ).fetchone()
            neg_conn.close()
            if row:
                return "conflict"

        db_path = os.path.join(db_dir, "positive_binders.db")
        conn = sqlite3.connect(db_path)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS positive_binders (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                assay_name TEXT NOT NULL,
                pose_file TEXT NOT NULL,
                directory TEXT NOT NULL,
                pose_full_path TEXT NOT NULL,
                flagged_at TEXT NOT NULL,
                UNIQUE(assay_name, pose_file, directory)
            )
        """)
        conn.commit()
        existing = conn.execute(
            "SELECT 1 FROM positive_binders WHERE assay_name=? AND pose_file=?",
            (assay_name, pose_file)
        ).fetchone()
        if existing:
            conn.close()
            return "duplicate"
        try:
            conn.execute(
                "INSERT INTO positive_binders (assay_name, pose_file, directory, pose_full_path, flagged_at) VALUES (?, ?, ?, ?, ?)",
                (assay_name, pose_file, directory, pose_full_path, datetime.now().isoformat())
            )
            conn.commit()
            status = "saved"
        except sqlite3.IntegrityError:
            status = "duplicate"
        conn.close()
        return status
    except Exception as e:
        return f"error:{e}"


def save_negative_binder(project_path: str, assay_name: str, pose_file: str, directory: str, pose_full_path: str) -> str:
    """
    Save a pose flagged as a negative binder to the negative_binders.db database.

    The database is created at <project_path>/ml/training_sets/negative_binders.db.
    Returns a status string: 'saved', 'duplicate', 'conflict', or 'error:<message>'.
    'conflict' means the pose is already registered as a positive binder.

    Args:
        project_path (str): Root path of the active project.
        assay_name (str): Name of the docking assay.
        pose_file (str): PDB filename of the pose.
        directory (str): Pose folder name (e.g. 'most_stable_poses').
        pose_full_path (str): Absolute path to the PDB file.
    """
    try:
        db_dir = os.path.join(project_path, "ml", "training_sets")
        os.makedirs(db_dir, exist_ok=True)

        pos_db_path = os.path.join(db_dir, "positive_binders.db")
        if os.path.exists(pos_db_path):
            pos_conn = sqlite3.connect(pos_db_path)
            pos_conn.execute("""
                CREATE TABLE IF NOT EXISTS positive_binders (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    assay_name TEXT NOT NULL,
                    pose_file TEXT NOT NULL,
                    directory TEXT NOT NULL,
                    pose_full_path TEXT NOT NULL,
                    flagged_at TEXT NOT NULL,
                    UNIQUE(assay_name, pose_file, directory)
                )
            """)
            row = pos_conn.execute(
                "SELECT 1 FROM positive_binders WHERE assay_name=? AND pose_file=?",
                (assay_name, pose_file)
            ).fetchone()
            pos_conn.close()
            if row:
                return "conflict"

        db_path = os.path.join(db_dir, "negative_binders.db")
        conn = sqlite3.connect(db_path)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS negative_binders (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                assay_name TEXT NOT NULL,
                pose_file TEXT NOT NULL,
                directory TEXT NOT NULL,
                pose_full_path TEXT NOT NULL,
                flagged_at TEXT NOT NULL,
                UNIQUE(assay_name, pose_file, directory)
            )
        """)
        conn.commit()
        existing = conn.execute(
            "SELECT 1 FROM negative_binders WHERE assay_name=? AND pose_file=?",
            (assay_name, pose_file)
        ).fetchone()
        if existing:
            conn.close()
            return "duplicate"
        try:
            conn.execute(
                "INSERT INTO negative_binders (assay_name, pose_file, directory, pose_full_path, flagged_at) VALUES (?, ?, ?, ?, ?)",
                (assay_name, pose_file, directory, pose_full_path, datetime.now().isoformat())
            )
            conn.commit()
            status = "saved"
        except sqlite3.IntegrityError:
            status = "duplicate"
        conn.close()
        return status
    except Exception as e:
        return f"error:{e}"


def get_binders_registry(project_path: str, binder_type: str) -> "pd.DataFrame":
    """
    Load the positive or negative binders registry from the ml/training_sets database.

    Args:
        project_path (str): Root path of the active project.
        binder_type (str): Either 'positive' or 'negative'.

    Returns:
        pd.DataFrame or None if the database/table does not exist.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", f"{binder_type}_binders.db")
    table = f"{binder_type}_binders"
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(f"SELECT * FROM {table} ORDER BY flagged_at DESC", conn)
        conn.close()
        return df
    except Exception:
        return None


def remove_binder(project_path: str, binder_type: str, assay_name: str, pose_file: str, directory: str) -> str:
    """
    Remove a single entry from the positive or negative binders registry.

    Deletes the row matching (assay_name, pose_file, directory) from the corresponding table.
    Returns 'removed', 'not_found', or 'error:<message>'.

    Args:
        project_path (str): Root path of the active project.
        binder_type (str): Either 'positive' or 'negative'.
        assay_name (str): Name of the docking assay.
        pose_file (str): PDB filename of the pose.
        directory (str): Pose folder name.
    """
    try:
        db_path = os.path.join(project_path, "ml", "training_sets", f"{binder_type}_binders.db")
        table = f"{binder_type}_binders"
        if not os.path.exists(db_path):
            return "not_found"
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            f"DELETE FROM {table} WHERE assay_name = ? AND pose_file = ? AND directory = ?",
            (assay_name, pose_file, directory)
        )
        conn.commit()
        removed = cursor.rowcount > 0
        conn.close()
        return "removed" if removed else "not_found"
    except Exception as e:
        return f"error:{e}"


def clear_binders(project_path: str, binder_type: str) -> str:
    """
    Delete all entries from the positive or negative binders registry.

    Returns 'cleared', 'not_found', or 'error:<message>'.

    Args:
        project_path (str): Root path of the active project.
        binder_type (str): Either 'positive' or 'negative'.
    """
    try:
        db_path = os.path.join(project_path, "ml", "training_sets", f"{binder_type}_binders.db")
        table = f"{binder_type}_binders"
        if not os.path.exists(db_path):
            return "not_found"
        conn = sqlite3.connect(db_path)
        conn.execute(f"DELETE FROM {table}")
        conn.commit()
        conn.close()
        return "cleared"
    except Exception as e:
        return f"error:{e}"


def consolidate_training_set(project_path: str, training_set_id: str, notes: str = "") -> str:
    """
    Snapshot the current positive and negative binders registries into a named training set.

    The snapshot is stored in <project_path>/ml/training_sets/training_sets_snapshots.db
    using two tables:
      - training_set_snapshots: one row per snapshot (id, counts, notes, timestamp)
      - training_set_entries:   one row per binder entry, tagged with training_set_id and binder_type

    Returns a status string: 'saved', 'duplicate', or 'error:<message>'.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): Unique identifier for this training set snapshot.
        notes (str): Optional free-text notes to attach to the snapshot.
    """
    try:
        df_pos = get_binders_registry(project_path, "positive")
        df_neg = get_binders_registry(project_path, "negative")
        n_pos = len(df_pos) if df_pos is not None and not df_pos.empty else 0
        n_neg = len(df_neg) if df_neg is not None and not df_neg.empty else 0

        if n_pos == 0 and n_neg == 0:
            return "error:No binders registered in either registry"

        db_dir = os.path.join(project_path, "ml", "training_sets")
        os.makedirs(db_dir, exist_ok=True)
        db_path = os.path.join(db_dir, "training_sets_snapshots.db")
        conn = sqlite3.connect(db_path)

        conn.execute("""
            CREATE TABLE IF NOT EXISTS training_set_snapshots (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                training_set_id TEXT NOT NULL UNIQUE,
                created_at TEXT NOT NULL,
                n_positives INTEGER NOT NULL,
                n_negatives INTEGER NOT NULL,
                notes TEXT
            )
        """)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS training_set_entries (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                training_set_id TEXT NOT NULL,
                binder_type TEXT NOT NULL,
                assay_name TEXT NOT NULL,
                pose_file TEXT NOT NULL,
                directory TEXT NOT NULL,
                pose_full_path TEXT NOT NULL,
                flagged_at TEXT NOT NULL
            )
        """)
        conn.commit()

        try:
            conn.execute(
                "INSERT INTO training_set_snapshots (training_set_id, created_at, n_positives, n_negatives, notes) VALUES (?, ?, ?, ?, ?)",
                (training_set_id, datetime.now().isoformat(), n_pos, n_neg, notes)
            )
            conn.commit()
        except sqlite3.IntegrityError:
            conn.close()
            return "duplicate"

        for binder_type, df in [("positive", df_pos), ("negative", df_neg)]:
            if df is not None and not df.empty:
                for _, row in df.iterrows():
                    conn.execute(
                        "INSERT INTO training_set_entries (training_set_id, binder_type, assay_name, pose_file, directory, pose_full_path, flagged_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
                        (training_set_id, binder_type, row["assay_name"], row["pose_file"],
                         row["directory"], row["pose_full_path"], row["flagged_at"])
                    )
        conn.commit()
        conn.close()
        return "saved"
    except Exception as e:
        return f"error:{e}"


def get_training_set_entries(project_path: str, training_set_id: str) -> "pd.DataFrame":
    """
    Load all binder entries for a given consolidated training set.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): ID of the training set snapshot to load.

    Returns:
        pd.DataFrame with columns (binder_type, assay_name, pose_file, directory,
        pose_full_path, flagged_at), ordered positive-first then negative,
        or None on error.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(
            """SELECT binder_type, assay_name, pose_file, directory, pose_full_path, flagged_at
               FROM training_set_entries
               WHERE training_set_id = ?
               ORDER BY binder_type DESC, assay_name, pose_file""",
            conn,
            params=(training_set_id,)
        )
        conn.close()
        return df if not df.empty else None
    except Exception:
        return None


def get_training_set_snapshots(project_path: str) -> "pd.DataFrame":
    """
    Load the list of all consolidated training set snapshots.

    Args:
        project_path (str): Root path of the active project.

    Returns:
        pd.DataFrame with snapshot metadata, or None if none exist yet.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(
            "SELECT training_set_id, created_at, n_positives, n_negatives, notes FROM training_set_snapshots ORDER BY created_at ASC",
            conn
        )
        conn.close()
        return df
    except Exception:
        return None


def get_training_set_fingerprint_status(project_path: str) -> "dict":
    """
    Return a dict mapping training_set_id → fingerprint status string.

    Status values:
      '✅ Complete'      — every entry in the snapshot has at least one fingerprint
      '⚠️ Partial (X/Y)' — some but not all entries have fingerprints
      '❌ None'           — no fingerprints computed for this snapshot

    Entries are counted as distinct (assay_name, pose_file, directory) tuples so
    that multiple prolif_conditions runs for the same pose do not inflate the count.

    Args:
        project_path (str): Root path of the active project.

    Returns:
        dict[str, str]: Mapping of training_set_id to status string.
                        Empty dict if the database or tables do not exist.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return {}
    try:
        conn = sqlite3.connect(db_path)
        totals = pd.read_sql_query(
            "SELECT training_set_id, COUNT(*) AS total FROM training_set_entries GROUP BY training_set_id",
            conn,
        )
        computed = pd.read_sql_query(
            """SELECT training_set_id, COUNT(DISTINCT assay_name || '|' || pose_file || '|' || directory) AS done
               FROM training_set_fingerprints GROUP BY training_set_id""",
            conn,
        )
        conn.close()
    except Exception:
        return {}

    merged = totals.merge(computed, on="training_set_id", how="left")
    merged["done"] = merged["done"].fillna(0).astype(int)

    status = {}
    for _, row in merged.iterrows():
        ts_id = row["training_set_id"]
        total = int(row["total"])
        done = int(row["done"])
        if done == 0:
            status[ts_id] = "❌ None"
        elif done >= total:
            status[ts_id] = "✅ Complete"
        else:
            status[ts_id] = f"⚠️ Partial ({done}/{total})"
    return status


def get_training_set_prolif_conditions_used(project_path: str) -> "dict":
    """
    Return a dict mapping training_set_id → human-readable ProLIF conditions string.

    For each training set that has fingerprints, collects the distinct
    prolif_conditions_id values and resolves their descriptions from the
    ProLIF_Conditions table in the docking params DB.  Falls back to showing
    the raw IDs when the descriptions cannot be retrieved.

    Args:
        project_path (str): Root path of the active project.

    Returns:
        dict[str, str]: {training_set_id: "ID 1: desc1; ID 2: desc2"} or
                        "❌ None" when no fingerprints exist.
                        Empty dict if the snapshots database does not exist.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return {}

    try:
        conn = sqlite3.connect(db_path)
        rows = pd.read_sql_query(
            "SELECT DISTINCT training_set_id, prolif_conditions_id FROM training_set_fingerprints ORDER BY training_set_id, prolif_conditions_id",
            conn,
        )
        conn.close()
    except Exception:
        return {}

    if rows.empty:
        return {}

    # Try to resolve descriptions from docking params DB
    params_db = os.path.join(project_path, "docking", "params", "params.db")
    desc_map = {}
    try:
        if os.path.exists(params_db):
            conn = sqlite3.connect(params_db)
            conds = pd.read_sql_query("SELECT id, description FROM ProLIF_Conditions", conn)
            conn.close()
            desc_map = dict(zip(conds["id"], conds["description"]))
    except Exception:
        pass

    result = {}
    for ts_id, group in rows.groupby("training_set_id"):
        parts = []
        for cid in group["prolif_conditions_id"]:
            desc = desc_map.get(cid)
            parts.append(f"ID {cid}: {desc}" if desc else f"ID {cid}")
        result[ts_id] = "; ".join(parts)
    return result


def get_training_set_fingerprints_for_pose(project_path: str, training_set_id: str, assay_name: str, pose_file: str, directory: str) -> "list[dict] | None":
    """
    Retrieve all computed ProLIF fingerprints for a specific pose in a training set snapshot.

    A pose may have multiple fingerprints if computed under different ProLIF conditions.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): ID of the training set snapshot.
        assay_name (str): Name of the docking assay.
        pose_file (str): PDB filename of the pose.
        directory (str): Pose folder name (e.g. 'most_stable_poses').

    Returns:
        List of dicts with keys 'prolif_conditions_id' and 'fingerprint' (dict of interaction→bool),
        ordered by prolif_conditions_id. Returns None if no fingerprints exist for this pose.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        rows = pd.read_sql_query(
            """SELECT prolif_conditions_id, fingerprint_json
               FROM training_set_fingerprints
               WHERE training_set_id=? AND assay_name=? AND pose_file=? AND directory=?
               ORDER BY prolif_conditions_id""",
            conn,
            params=(training_set_id, assay_name, pose_file, directory),
        )
        conn.close()
        if rows.empty:
            return None
        return [
            {"prolif_conditions_id": row["prolif_conditions_id"], "fingerprint": json.loads(row["fingerprint_json"])}
            for _, row in rows.iterrows()
        ]
    except Exception:
        return None


def get_all_training_set_interactions(project_path: str, training_set_id: str, prolif_conditions_id: int) -> "list[str]":
    """
    Return a sorted list of all unique interaction keys present across every pose
    in a training set snapshot for a given ProLIF conditions ID.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): ID of the training set snapshot.
        prolif_conditions_id (int): ProLIF conditions identifier.

    Returns:
        Sorted list of interaction name strings. Empty list if none found.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return []
    try:
        conn = sqlite3.connect(db_path)
        rows = pd.read_sql_query(
            "SELECT fingerprint_json FROM training_set_fingerprints WHERE training_set_id=? AND prolif_conditions_id=?",
            conn,
            params=(training_set_id, prolif_conditions_id),
        )
        conn.close()
        all_keys: set = set()
        for _, row in rows.iterrows():
            fp = json.loads(row["fingerprint_json"])
            all_keys.update(fp.keys())

        def _seq_sort_key(name: str):
            first = name.split("_")[0]
            try:
                return (int(first), name)
            except ValueError:
                return (0, name)

        return sorted(all_keys, key=_seq_sort_key)
    except Exception:
        return []


def get_training_set_interaction_frequencies(project_path: str, training_set_id: str, prolif_conditions_id: int) -> "pd.DataFrame | None":
    """
    Compute the frequency (count and percentage) of each interaction across all poses
    in a training set snapshot for a given ProLIF conditions ID.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): ID of the training set snapshot.
        prolif_conditions_id (int): ProLIF conditions identifier.

    Returns:
        DataFrame with columns ['Interaction', 'Count', 'Frequency (%)'], sorted descending by Count.
        Returns None if no data is found.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        rows = pd.read_sql_query(
            "SELECT fingerprint_json FROM training_set_fingerprints WHERE training_set_id=? AND prolif_conditions_id=?",
            conn,
            params=(training_set_id, prolif_conditions_id),
        )
        conn.close()
        if rows.empty:
            return None
        n_poses = len(rows)
        counts: dict = {}
        for _, row in rows.iterrows():
            fp = json.loads(row["fingerprint_json"])
            for k, v in fp.items():
                if v:
                    counts[k] = counts.get(k, 0) + 1
        if not counts:
            return None
        def _seq_sort_key(name: str):
            first = name.split("_")[0]
            try:
                return (int(first), name)
            except ValueError:
                return (0, name)

        df = pd.DataFrame([
            {"Interaction": k, "Count": v, "Frequency (%)": round(100.0 * v / n_poses, 1)}
            for k, v in counts.items()
        ])
        df = df.iloc[sorted(range(len(df)), key=lambda i: _seq_sort_key(df["Interaction"].iloc[i]))].reset_index(drop=True)
        return df
    except Exception:
        return None


def restore_binders_from_snapshot(project_path: str, training_set_id: str) -> "dict | str":
    """
    Repopulate the Positive and Negative Binders tables from a training set snapshot.

    Entries already present in the destination table (matched on assay_name, pose_file,
    directory) are skipped; only genuinely new records are inserted.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): ID of the training set snapshot to restore from.

    Returns:
        dict with keys pos_inserted, pos_skipped, neg_inserted, neg_skipped,
        or a string 'error:<message>' on failure.
    """
    try:
        df = get_training_set_entries(project_path, training_set_id)
        if df is None or df.empty:
            return "error:No entries found for this snapshot"

        db_dir = os.path.join(project_path, "ml", "training_sets")
        os.makedirs(db_dir, exist_ok=True)

        counts = {"pos_inserted": 0, "pos_skipped": 0, "neg_inserted": 0, "neg_skipped": 0}

        for binder_type in ("positive", "negative"):
            db_path = os.path.join(db_dir, f"{binder_type}_binders.db")
            table = f"{binder_type}_binders"
            subset = df[df["binder_type"] == binder_type]
            if subset.empty:
                continue

            ins_key = "pos_inserted" if binder_type == "positive" else "neg_inserted"
            skip_key = "pos_skipped" if binder_type == "positive" else "neg_skipped"

            conn = sqlite3.connect(db_path)
            conn.execute(f"""
                CREATE TABLE IF NOT EXISTS {table} (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    assay_name TEXT NOT NULL,
                    pose_file TEXT NOT NULL,
                    directory TEXT NOT NULL,
                    pose_full_path TEXT NOT NULL,
                    flagged_at TEXT NOT NULL,
                    UNIQUE(assay_name, pose_file, directory)
                )
            """)
            conn.commit()

            for _, row in subset.iterrows():
                cursor = conn.execute(
                    f"INSERT OR IGNORE INTO {table} (assay_name, pose_file, directory, pose_full_path, flagged_at) VALUES (?, ?, ?, ?, ?)",
                    (row["assay_name"], row["pose_file"], row["directory"], row["pose_full_path"], row["flagged_at"])
                )
                if cursor.rowcount > 0:
                    counts[ins_key] += 1
                else:
                    counts[skip_key] += 1

            conn.commit()
            conn.close()

        return counts
    except Exception as e:
        return f"error:{e}"


def delete_training_set_snapshots(project_path: str, training_set_ids: list) -> str:
    """
    Delete one or more training set snapshots and all their associated entries.

    Removes matching rows from both training_set_snapshots and training_set_entries
    for each id in training_set_ids.
    Returns 'deleted', 'not_found', or 'error:<message>'.

    Args:
        project_path (str): Root path of the active project.
        training_set_ids (list[str]): List of training_set_id values to delete.
    """
    if not training_set_ids:
        return "not_found"
    try:
        db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
        if not os.path.exists(db_path):
            return "not_found"
        placeholders = ",".join("?" * len(training_set_ids))
        conn = sqlite3.connect(db_path)
        conn.execute(f"DELETE FROM training_set_snapshots WHERE training_set_id IN ({placeholders})", training_set_ids)
        conn.execute(f"DELETE FROM training_set_entries WHERE training_set_id IN ({placeholders})", training_set_ids)
        conn.execute(f"DELETE FROM training_set_fingerprints WHERE training_set_id IN ({placeholders})", training_set_ids)
        conn.commit()
        conn.close()
        return "deleted"
    except Exception as e:
        return f"error:{e}"


def get_table_columns(db_path: str, table_name: str) -> list:
    """
    Return the list of column names for a given table in a SQLite database.

    Args:
        db_path (str): Path to the SQLite database.
        table_name (str): Name of the table.

    Returns:
        list[str]: Column names, or empty list on error.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(f"SELECT * FROM [{table_name}] LIMIT 0")
        columns = [desc[0] for desc in cursor.description]
        conn.close()
        return columns
    except Exception:
        return []


def read_table_columns_as_dataframe(db_path: str, table_name: str, columns: list) -> "pd.DataFrame":
    """
    Read selected columns from a SQLite table and return as a pandas DataFrame.

    Args:
        db_path (str): Path to the SQLite database.
        table_name (str): Name of the table.
        columns (list): List of column names to include.

    Returns:
        pd.DataFrame: DataFrame with the selected columns, or empty DataFrame on error.
    """
    try:
        col_expr = ", ".join(f"[{c}]" for c in columns)
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(f"SELECT {col_expr} FROM [{table_name}]", conn)
        conn.close()
        return df
    except Exception as e:
        return pd.DataFrame()


def create_subset_table(db_path: str, source_table: str, new_table: str, row_indices: list, columns: list = None) -> str:
    """
    Create a new table in a SQLite database as a subset of rows from an existing table.

    Copies the rows identified by their 0-based positional indices (as returned by
    pandas iloc) from source_table into a new table called new_table. Only the
    columns listed in `columns` are written; if `columns` is None all columns are
    preserved.  The new table must not already exist.

    Returns 'created', 'duplicate', or 'error:<message>'.

    Args:
        db_path (str): Path to the SQLite database.
        source_table (str): Name of the source table.
        new_table (str): Name for the new subset table.
        row_indices (list[int]): 0-based row positions to include in the subset.
        columns (list[str], optional): Column names to include. Defaults to all columns.
    """
    if not row_indices:
        return "error:No rows selected"
    try:
        conn = sqlite3.connect(db_path)
        existing = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?", (new_table,)
        ).fetchone()
        if existing:
            conn.close()
            return "duplicate"
        if columns:
            col_expr = ", ".join(f"[{c}]" for c in columns)
            df_full = pd.read_sql_query(f"SELECT {col_expr} FROM [{source_table}]", conn)
        else:
            df_full = pd.read_sql_query(f"SELECT * FROM [{source_table}]", conn)
        df_subset = df_full.iloc[row_indices]
        df_subset.to_sql(new_table, conn, if_exists="fail", index=False)
        conn.close()
        return "created"
    except Exception as e:
        return f"error:{e}"



def depict_table_to_images(db_path: str, table_name: str,
                           max_molecules: int = 25,
                           molecules_per_image: int = 25,
                           mol_image_size: tuple = (300, 300),
                           legend_col: str = None) -> list:
    """
    Generate in-memory molecule grid images from a chemspace table for Streamlit display.

    Mirrors the logic of ChemSpace.depict_table but returns PIL Images instead of
    saving files to disk, so they can be rendered inline with st.image().

    Args:
        db_path (str): Path to the chemspace SQLite database.
        table_name (str): Name of the table to depict.
        max_molecules (int): Maximum number of molecules to render (-1 for all).
        molecules_per_image (int): Number of molecules per grid image.
        mol_image_size (tuple): (width, height) of each molecule cell in the grid.
        legend_col (str): Column to use as molecule label. Defaults to 'id' if present,
                          otherwise the first column.

    Returns:
        list[PIL.Image.Image]: List of grid images ready for st.image(). Empty on error.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError:
        raise ImportError("RDKit is required for depiction. Install with: conda install -c conda-forge rdkit")

    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(f"SELECT * FROM [{table_name}]", conn)
        conn.close()
    except Exception as e:
        raise RuntimeError(f"Could not read table '{table_name}': {e}")

    if df.empty:
        return []

    if 'smiles' not in df.columns:
        raise ValueError(f"Table '{table_name}' has no 'smiles' column.")

    total = len(df)
    if max_molecules == -1 or max_molecules >= total:
        subset = df
    else:
        subset = df.head(max_molecules)

    # Determine legend column: use caller's choice, else 'id', else first column
    if legend_col and legend_col in subset.columns:
        effective_legend_col = legend_col
    elif 'id' in subset.columns:
        effective_legend_col = 'id'
    else:
        effective_legend_col = subset.columns[0]

    # Calculate grid column count from molecules_per_image
    grid_cols = max(1, int(molecules_per_image ** 0.5))

    images = []
    batch_size = molecules_per_image
    rows_list = list(subset.itertuples(index=False))

    for batch_start in range(0, len(rows_list), batch_size):
        batch = rows_list[batch_start:batch_start + batch_size]
        mols = []
        legends = []

        for row in batch:
            smiles = getattr(row, 'smiles', None)
            if not smiles:
                continue
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                continue
            mols.append(mol)
            legend_val = str(getattr(row, effective_legend_col, ''))[:40]
            legends.append(legend_val)

        if not mols:
            continue

        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=grid_cols,
            subImgSize=mol_image_size,
            legends=legends,
            returnPNG=False
        )
        images.append(img)

    return images


def get_prolif_conditions_records(project_path: str) -> "list[dict] | None":
    """
    Load all ProLIF conditions stored in params.db.

    Returns a list of dicts with keys ``id``, ``description``, and
    ``conditions`` (the parsed JSON dict), ordered by id ascending.
    Returns None when the database or table does not exist.
    """
    db_path = os.path.join(project_path, 'docking', 'params', 'params.db')
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='ProLIF_Conditions'"
        )
        if not cursor.fetchone():
            conn.close()
            return None
        cursor.execute(
            "SELECT id, description, conditions FROM ProLIF_Conditions ORDER BY id ASC"
        )
        rows = cursor.fetchall()
        conn.close()
        if not rows:
            return None
        return [
            {"id": r[0], "description": r[1], "conditions": json.loads(r[2])}
            for r in rows
        ]
    except Exception:
        return None


def export_prolif_conditions_to_file(project_path: str, record_id: int, output_path: str) -> str:
    """
    Export a ProLIF conditions record from params.db to a portable JSON file.

    Mirrors the behaviour of MolDock.export_prolif_conditions().
    Returns 'exported', or 'error:<message>' on failure.

    Args:
        project_path (str): Root path of the active project.
        record_id (int): ID of the ProLIF_Conditions record to export.
        output_path (str): Destination file path for the JSON export.
    """
    db_path = os.path.join(project_path, 'docking', 'params', 'params.db')
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT description, conditions FROM ProLIF_Conditions WHERE id = ?",
            (record_id,)
        )
        row = cursor.fetchone()
        conn.close()
        if row is None:
            return f"error:No record found for ID {record_id}"
        description, conditions_json = row
        export_payload = {
            "description": description,
            "conditions": json.loads(conditions_json),
        }
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(export_payload, f, indent=2)
        return "exported"
    except Exception as e:
        return f"error:{e}"


def generate_vmd_script(pose_path: str, ref_pdb_path: str = None) -> str:
    """
    Generate a VMD Tcl script that loads a docked pose and, optionally, a
    reference structure.

    Args:
        pose_path:    Absolute path to the docked pose PDB file.
        ref_pdb_path: Absolute path to the reference PDB file (optional).

    Returns:
        VMD Tcl script as a plain string.
    """
    lines = [
        "# VMD script generated by TidyScreen",
        f"# Docked pose: {os.path.basename(pose_path)}",
        "",
        "# --- Docked pose (Licorice, coloured by element) ---",
        f'mol new "{pose_path}" type pdb waitfor all',
        "mol delrep 0 top",
        "mol representation Licorice 0.300000 12.000000 12.000000",
        "mol color Element",
        "mol selection {all}",
        "mol material Opaque",
        "mol addrep top",
    ]
    if ref_pdb_path:
        lines += [
            "",
            f"# --- Reference structure: {os.path.basename(ref_pdb_path)} (CPK, coloured by element) ---",
            f'mol new "{ref_pdb_path}" type pdb waitfor all',
            "mol delrep 0 top",
            "mol representation CPK 1.000000 0.300000 12.000000 12.000000",
            "mol color Element",
            "mol selection {all}",
            "mol material Opaque",
            "mol addrep top",
        ]
    lines += [
        "",
        "display resetview",
    ]
    return "\n".join(lines)


def export_training_set_fingerprints_as_csv_bytes(project_path: str, training_set_id: str) -> "bytes | None":
    """
    Build a CSV of all precomputed fingerprints for the given training set snapshot
    and return it as UTF-8 bytes suitable for st.download_button.

    Returns None if no fingerprints are found or the database is missing.
    """
    db_path = os.path.join(project_path, 'ml', 'training_sets', 'training_sets_snapshots.db')
    if not os.path.exists(db_path):
        return None
    try:
        conn = sqlite3.connect(db_path)
        rows = pd.read_sql_query(
            """SELECT assay_name, pose_file, label, fingerprint_json
               FROM training_set_fingerprints
               WHERE training_set_id = ?
               ORDER BY binder_type DESC, assay_name, pose_file""",
            conn, params=(training_set_id,)
        )
        conn.close()
    except Exception:
        return None

    if rows.empty:
        return None

    ligpose_col = (
        rows['pose_file'].str.replace('.pdb', '', regex=False) + '_' + rows['assay_name']
    ).reset_index(drop=True)
    label_col = rows['label'].reset_index(drop=True)

    fp_dicts = [json.loads(fp) for fp in rows['fingerprint_json']]
    fp_df = pd.DataFrame(fp_dicts).fillna(0).astype(int).reset_index(drop=True)

    out_df = pd.concat(
        [ligpose_col.rename('ligpose'), fp_df, label_col.rename('label')],
        axis=1
    )

    buffer = StringIO()
    out_df.to_csv(buffer, index=False)
    return buffer.getvalue().encode('utf-8')


def export_pdb_model(pdbs_db_path, file_id, output_dir):
    """
    Export a PDB model blob from the database to a file on disk.

    Args:
        pdbs_db_path (str): Path to the pdbs.db SQLite database.
        file_id (str): The file_id of the PDB model to export.
        output_dir (str): Directory where the PDB file will be written.

    Returns:
        tuple[bool, str]: (success, message) where message is the output path on success
                          or an error description on failure.
    """
    try:
        conn = sqlite3.connect(pdbs_db_path)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT filename, pdb_blob FROM pdb_models WHERE file_id = ?",
            (int(file_id),)
        )
        row = cursor.fetchone()
        conn.close()

        if row is None or row[1] is None:
            return False, f"No PDB blob found for file_id '{file_id}'"

        filename, pdb_blob = row
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, filename)

        with open(output_path, "wb") as f:
            f.write(pdb_blob)

        return True, output_path

    except Exception as e:
        return False, str(e)


def populate_results_with_mmgbsa_energies(db_path):
    """
    Read all per-pose MMGBSA decomposition data from processed_mmgbsa_decomposition_json,
    sum the 'total' and 'gas' columns for each pose, and write the results into
    mmgbsa_total_energy and mmgbsa_gas_energy columns of the Results table.

    Columns are added to Results via ALTER TABLE if they do not already exist.

    Returns:
        tuple[int, list[str]]: (number of rows updated, list of error messages)
    """
    errors = []
    updated = 0
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Ensure destination columns exist in Results
        cursor.execute("PRAGMA table_info(Results)")
        existing_columns = {row[1] for row in cursor.fetchall()}
        for col_name, col_type in [('mmgbsa_total_energy', 'REAL'), ('mmgbsa_gas_energy', 'REAL')]:
            if col_name not in existing_columns:
                cursor.execute(f"ALTER TABLE Results ADD COLUMN {col_name} {col_type}")

        # Read all decomposition entries
        cursor.execute("SELECT pose_id, data FROM processed_mmgbsa_decomposition_json")
        rows = cursor.fetchall()

        for pose_id, mmgbsa_json in rows:
            try:
                mmgbsa_df = pd.read_json(StringIO(mmgbsa_json), orient='split')
                total_energy = round(float(mmgbsa_df['total'].sum()), 3)
                gas_energy = round(float(mmgbsa_df['gas'].sum()), 3)
            except Exception as e:
                errors.append(f"Pose {pose_id}: could not parse MMGBSA data — {e}")
                total_energy, gas_energy = -1, -1

            cursor.execute(
                "UPDATE Results SET mmgbsa_total_energy = ?, mmgbsa_gas_energy = ? WHERE Pose_ID = ?",
                (total_energy, gas_energy, pose_id)
            )
            updated += 1

        conn.commit()
        conn.close()

    except Exception as e:
        errors.append(f"Database error: {e}")

    return updated, errors
