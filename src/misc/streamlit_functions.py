from io import StringIO, BytesIO
import sqlite3
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from datetime import datetime


def _db_mtime(db_path):
    """Modification time of db_path, used as a cache-busting key so cached
    reads below auto-invalidate whenever the underlying SQLite file changes."""
    try:
        return os.path.getmtime(db_path)
    except OSError:
        return None


_TS_META_TABLE = "_ts_table_meta"


def _sync_table_triggers(conn, table_names):
    """
    Ensure every table in table_names has a row in _ts_table_meta plus AFTER
    INSERT/UPDATE/DELETE triggers that bump that row's version on every write.

    SQLite has no built-in per-table "last modified" primitive, and a
    whole-file mtime can't tell which table changed. TidyScreen tables are
    routinely modified in place (e.g. inchi_key/sdf_blob backfills via
    UPDATE ... WHERE id=?) without changing row counts, so row-count-based
    signatures aren't reliable either. Triggers are the only cheap way to
    capture every INSERT/UPDATE/DELETE regardless of which code path made it.

    version is a monotonically incrementing counter rather than a
    timestamp: wall-clock time (e.g. strftime('%s','now')) only has
    1-second resolution, so two writes to the same table within the same
    second would produce identical "last modified" values and silently
    fail to invalidate the cache.
    """
    cursor = conn.cursor()
    cursor.execute(f"CREATE TABLE IF NOT EXISTS [{_TS_META_TABLE}] (table_name TEXT PRIMARY KEY, version INTEGER)")
    cursor.execute("SELECT name FROM sqlite_master WHERE type='trigger'")
    existing_triggers = {row[0] for row in cursor.fetchall()}
    for table in table_names:
        trigger_names = [f"_ts_trg_{table}_{op.lower()}" for op in ("INSERT", "UPDATE", "DELETE")]
        # DROP TABLE also drops that table's triggers, but not its
        # _ts_table_meta row. If the table was dropped and recreated under
        # the same name, its triggers are missing here even though a stale
        # version row still exists — trusting that row would return a
        # cached row count from before the drop. Treat "triggers missing"
        # (true for brand-new tables too) as "must invalidate" and bump the
        # version unconditionally instead of leaving a stale/absent row.
        if not existing_triggers.issuperset(trigger_names):
            cursor.execute(
                f"""
                INSERT INTO [{_TS_META_TABLE}] (table_name, version) VALUES (?, 1)
                ON CONFLICT(table_name) DO UPDATE SET version = version + 1
                """,
                (table,)
            )
            for trigger_name, op in zip(trigger_names, ("INSERT", "UPDATE", "DELETE")):
                if trigger_name in existing_triggers:
                    continue
                cursor.execute(f"""
                    CREATE TRIGGER [{trigger_name}]
                    AFTER {op} ON [{table}]
                    BEGIN
                        INSERT INTO [{_TS_META_TABLE}] (table_name, version)
                        VALUES ('{table}', 1)
                        ON CONFLICT(table_name) DO UPDATE SET version = version + 1;
                    END;
                """)
        else:
            cursor.execute(
                f"INSERT OR IGNORE INTO [{_TS_META_TABLE}] (table_name, version) VALUES (?, 0)",
                (table,)
            )


def get_table_signatures(db_path):
    """
    Return an ordered {table_name: version} dict for every table in db_path,
    where version is incremented by a trigger (see _sync_table_triggers)
    whenever that specific table is written to.

    Use this instead of _db_mtime as the cache-busting key for per-table
    reads below, so modifying one table doesn't invalidate every other
    table's cache. Deliberately not @st.cache_data'd: it must run fresh on
    every rerun to notice changes, but it only touches small metadata
    tables/queries (no scanning of user table contents), so it's cheap.

    Connects with isolation_level=None (autocommit) so every statement
    (including the CREATE TRIGGER/INSERT calls in _sync_table_triggers)
    commits immediately instead of sitting in one long-lived transaction —
    and the connection is always closed via `finally`, even if a Streamlit
    rerun interrupts this function mid-call (Streamlit's rerun signal isn't
    a plain Exception, so a bare try/except here wouldn't catch it and the
    connection/lock would otherwise leak). Without both of these, a
    long-held write lock on chemspace.db can make later reads/writes fail
    with "database is locked", since this function runs on every page rerun.

    Returns {} on error (e.g. db_path doesn't exist yet).
    """
    conn = None
    try:
        conn = sqlite3.connect(db_path, isolation_level=None, timeout=10)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")
        tables = [row[0] for row in cursor.fetchall() if row[0] != _TS_META_TABLE]
        _sync_table_triggers(conn, tables)
        cursor.execute(f"SELECT table_name, version FROM [{_TS_META_TABLE}]")
        version_by_table = dict(cursor.fetchall())
        return {table: version_by_table.get(table) for table in tables}
    except Exception:
        return {}
    finally:
        if conn is not None:
            conn.close()


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
    
@st.cache_data(show_spinner=False)
def _get_table_row_count(db_path, table_name, mtime=None):
    """Cached per (db_path, table_name, mtime) — see get_tables_info."""
    conn = None
    try:
        conn = sqlite3.connect(db_path, timeout=10)
        cursor = conn.cursor()
        cursor.execute(f"SELECT COUNT(*) FROM [{table_name}]")
        return cursor.fetchone()[0]
    except Exception:
        return None
    finally:
        if conn is not None:
            conn.close()


def get_tables_info(db_path, signatures=None):
    """
    Return a DataFrame with all table names and number of rows for a given SQLite database.

    Each table's row count is cached individually, keyed by that table's own
    last-modified signature (see get_table_signatures) rather than a single
    whole-file mtime — so writing to one table no longer forces a COUNT(*)
    rescan of every other table on the next rerun.

    Args:
        db_path (str): Path to the SQLite database.
        signatures (dict, optional): {table_name: last_modified} from
            get_table_signatures(db_path). Computed internally if omitted.

    Returns:
        pd.DataFrame: DataFrame with columns ['table', 'rows'].
    """
    if signatures is None:
        signatures = get_table_signatures(db_path)
    try:
        info = []
        for table, mtime in signatures.items():
            nrows = _get_table_row_count(db_path, table, mtime=mtime)
            if nrows is not None:
                info.append({'table': table, 'rows': nrows})
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

def get_md_method_registers(db_path):
    """
    Retrieve all entries from the 'md_methods' table in md_methods.db.

    Args:
        db_path (str): Path to the md_methods.db SQLite database.

    Returns:
        pd.DataFrame: DataFrame with all MD method records, or None on error.
    """
    columns = ['method_id', 'method_name', 'description', 'engine', 'parameters']
    try:
        conn = sqlite3.connect(db_path)
        query = f"SELECT {', '.join(columns)} FROM md_methods"
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df
    except Exception:
        return None

def get_md_assay_status(assay_folder_path):
    """
    Determine the run status of an MD assay by inspecting output files.

    Returns one of:
      "⬜ Not started"  — no AMBER output files exist yet
      "🔄 Running"      — at least one output file exists but production has not finished
      "✅ Completed"    — production.out exists and contains 'FINAL RESULTS'
    """
    if not assay_folder_path or not os.path.isdir(assay_folder_path):
        return "⬜ Not started"

    out_files = ['min1.out', 'min2.out', 'heating.out', 'equilibration.out', 'production.out']
    any_out = any(os.path.exists(os.path.join(assay_folder_path, f)) for f in out_files)

    if not any_out:
        return "⬜ Not started"

    production_out = os.path.join(assay_folder_path, 'production.out')
    if os.path.exists(production_out):
        # Read only the tail to avoid loading multi-GB files
        try:
            with open(production_out, 'rb') as fh:
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(0, size - 4096), 0)
                tail = fh.read().decode('utf-8', errors='replace')
            if 'Total wall time:' in tail or 'FINAL RESULTS' in tail:
                return "✅ Completed"
        except OSError:
            pass

    return "🔄 Running"


def get_md_mmgbsa_status(folder, mmgbsa_results_json):
    """
    Return an emoji string indicating MM-GBSA computation status for an MD assay.

    Priority:
      - DB column populated → completed (parsed and stored)
      - mmgbsa_results.dat exists → background run finished, pending parse
      - mmgbsa/ folder exists → running in background
      - otherwise → not computed
    """
    if mmgbsa_results_json is not None and str(mmgbsa_results_json).strip():
        return '✅ Completed'
    if folder and os.path.isdir(folder):
        mmgbsa_folder = os.path.join(folder, 'mmgbsa')
        if os.path.isdir(mmgbsa_folder):
            if os.path.exists(os.path.join(mmgbsa_folder, 'mmgbsa_results.dat')):
                return '✅ Completed'
            return '🔄 Running'
    return '⬜ Not computed'


def get_md_assay_registers(db_path):
    """
    Retrieve all entries from the 'md_assays' table in md_registers.db.

    Args:
        db_path (str): Path to the md_registers.db SQLite database.

    Returns:
        pd.DataFrame: DataFrame with all MD assay records plus a 'status' column, or None on error.
    """
    columns = ['assay_id', 'md_assay', 'description', 'assay_folder_path',
               'docking_assay_id', 'ligand_name', 'pose_id',
               'receptor_template_name', 'md_parameters']
    try:
        conn = sqlite3.connect(db_path)
        # mmgbsa_results column is added lazily by _store_mmgbsa_results(); handle absence gracefully
        existing_cols = [r[1] for r in conn.execute("PRAGMA table_info(md_assays)").fetchall()]
        fetch_cols = columns + (['mmgbsa_results'] if 'mmgbsa_results' in existing_cols else [])
        query = f"SELECT {', '.join(fetch_cols)} FROM md_assays"
        df = pd.read_sql_query(query, conn)
        conn.close()
    except Exception:
        return None

    if 'mmgbsa_results' not in df.columns:
        df['mmgbsa_results'] = None

    try:
        df['status'] = df['assay_folder_path'].apply(get_md_assay_status)
    except Exception:
        df['status'] = '⬜ Not started'

    try:
        df['mmgbsa'] = df.apply(
            lambda r: get_md_mmgbsa_status(r['assay_folder_path'], r['mmgbsa_results']),
            axis=1,
        )
    except Exception:
        df['mmgbsa'] = '⬜ Not computed'

    return df

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

def get_assay_engine(project_path, assay_name):
    """
    Return the docking engine string for the named assay, or None if not found.
    """
    registry_db = os.path.join(project_path, "docking", "docking_registers", "docking_assays.db")
    try:
        conn = sqlite3.connect(registry_db)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT docking_engine FROM docking_assays WHERE assay_name = ?", (assay_name,)
        )
        row = cursor.fetchone()
        conn.close()
        return row[0] if row else None
    except Exception:
        return None


def _load_moldock_class():
    """
    Import MolDock from the source tree that contains this file, bypassing
    whatever version of tidyscreen the Streamlit conda env has installed.
    streamlit_functions.py lives at <pkg>/misc/streamlit_functions.py and
    moldock.py lives at <pkg>/moldock/moldock.py.
    """
    import importlib.util
    src_pkg = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    moldock_file = os.path.join(src_pkg, "moldock", "moldock.py")
    spec = importlib.util.spec_from_file_location("_tidyscreen_moldock_src", moldock_file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.MolDock


def extract_poses_for_assay(project_name, project_path, assay_name, selection, score_column=None, overwrite=False):
    """
    Extract docked poses for *assay_name* using the native MolDock method.

    *selection* is the string code passed to the engine-specific extractor:
      AutoDockGPU — "1" most stable, "2" most populated, "3" both, "4" all
      AutoDock Vina — "1" most stable, "2" all

    *score_column*, when given, is used as the "lowest energy" criterion for
    selections "1" and "3" instead of prompting interactively (which cannot
    happen from this non-interactive GUI process).

    *overwrite*, when True, replaces a previous extraction's output directory
    without prompting. When False and a previous extraction's output directory
    already contains files, MolDock.extract_docked_poses() raises FileExistsError
    (there is no terminal here to confirm interactively) — that message is
    returned to the caller as-is so the GUI can offer to retry with overwrite.

    Returns (success: bool, message: str).
    """
    try:
        MolDock = _load_moldock_class()
        moldock = MolDock.from_path(project_name, project_path)
        output_dir = moldock.extract_docked_poses(
            assay=assay_name, selection=selection, score_column=score_column, overwrite=overwrite
        )
        if output_dir:
            return (True, f"Poses extracted to: {output_dir}")
        return (False, "Extraction returned no output directory — check console for details.")
    except Exception as e:
        return (False, str(e))


def get_available_score_columns(results_db_path):
    """
    Return the scoring columns in the assay's Results table that have at least one
    valid (non-NULL, non -1) value, in priority order: docking_score,
    mmgbsa_total_energy, mmgbsa_gas_energy. Mirrors MolDock._prompt_score_column(),
    so the GUI can offer the same choice as the interactive CLI prompt.
    """
    try:
        conn = sqlite3.connect(results_db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(Results)")
        existing = {row[1] for row in cursor.fetchall()}
        candidates = []
        for col in ['docking_score', 'mmgbsa_total_energy', 'mmgbsa_gas_energy']:
            if col in existing:
                cursor.execute(
                    f"SELECT COUNT(*) FROM Results WHERE {col} IS NOT NULL AND {col} != -1"
                )
                if cursor.fetchone()[0] > 0:
                    candidates.append(col)
        conn.close()
        return candidates
    except Exception:
        return []


def get_compounds_for_lig_names(project_path, assay_name, lig_names):
    """
    Retrieve smiles, name, and flag from the chemspace table associated with
    *assay_name* for the given *lig_names* list.

    Steps:
      1. Look up table_name for assay_name in docking_assays.db
      2. Query chemspace.db for rows whose name is in lig_names

    Returns (DataFrame, None) on success or (None, error_string) on failure.
    The DataFrame always has columns [smiles, name, flag]; flag is '' when the
    column is absent from the chemspace table.
    """
    if not lig_names:
        return None, "No ligand names provided."

    # ── Step 1: resolve chemspace table name from the docking registry ────────
    registry_db = os.path.join(
        project_path, "docking", "docking_registers", "docking_assays.db"
    )
    if not os.path.exists(registry_db):
        return None, f"Docking registry not found: {registry_db}"
    try:
        conn = sqlite3.connect(registry_db)
        row = conn.execute(
            "SELECT table_name FROM docking_assays WHERE assay_name = ?", (assay_name,)
        ).fetchone()
        conn.close()
    except Exception as e:
        return None, f"Error reading docking registry: {e}"
    if not row:
        return None, f"Assay '{assay_name}' not found in docking registry."
    table_name = row[0]

    # ── Step 2: query chemspace.db ────────────────────────────────────────────
    chemspace_db = os.path.join(
        project_path, "chemspace", "processed_data", "chemspace.db"
    )
    if not os.path.exists(chemspace_db):
        return None, f"ChemSpace database not found: {chemspace_db}"
    try:
        conn = sqlite3.connect(chemspace_db)
        # Discover available columns — flag may be absent in some table variants
        cursor = conn.cursor()
        cursor.execute(f'PRAGMA table_info("{table_name}")')
        available = {r[1].lower() for r in cursor.fetchall()}
        if not available:
            conn.close()
            return None, f"Table '{table_name}' not found in chemspace database."
        if "inchi_key" not in available:
            conn.close()
            return None, f"Table '{table_name}' has no 'inchi_key' column."
        select_cols = ["smiles", "name"]
        if "flag" in available:
            select_cols.append("flag")
        placeholders = ",".join("?" * len(lig_names))
        query = (
            f'SELECT {", ".join(select_cols)} FROM "{table_name}" '
            f"WHERE inchi_key IN ({placeholders})"
        )
        df = pd.read_sql_query(query, conn, params=list(lig_names))
        conn.close()
    except Exception as e:
        return None, f"Error querying chemspace table '{table_name}': {e}"

    if df.empty:
        return None, f"No compounds found in '{table_name}' matching the given LigNames."

    # Ensure flag column is always present
    if "flag" not in df.columns:
        df["flag"] = ""
    return df[["smiles", "name", "flag"]], None


def get_ligand_name_map(project_path, assay_name, lig_names):
    """
    Resolve a LigName (chemspace inchi_key) -> compound 'name' map for the
    given *lig_names*, using the chemspace table associated with *assay_name*.

    Returns {} if the docking registry, chemspace database, table, or the
    'inchi_key'/'name' columns cannot be resolved.
    """
    if not lig_names:
        return {}

    registry_db = os.path.join(
        project_path, "docking", "docking_registers", "docking_assays.db"
    )
    if not os.path.exists(registry_db):
        return {}
    try:
        conn = sqlite3.connect(registry_db)
        row = conn.execute(
            "SELECT table_name FROM docking_assays WHERE assay_name = ?", (assay_name,)
        ).fetchone()
        conn.close()
    except Exception:
        return {}
    if not row:
        return {}
    table_name = row[0]

    chemspace_db = os.path.join(
        project_path, "chemspace", "processed_data", "chemspace.db"
    )
    if not os.path.exists(chemspace_db):
        return {}
    try:
        conn = sqlite3.connect(chemspace_db)
        cursor = conn.cursor()
        cursor.execute(f'PRAGMA table_info("{table_name}")')
        available = {r[1].lower() for r in cursor.fetchall()}
        if "inchi_key" not in available or "name" not in available:
            conn.close()
            return {}
        placeholders = ",".join("?" * len(lig_names))
        df = pd.read_sql_query(
            f'SELECT inchi_key, name FROM "{table_name}" WHERE inchi_key IN ({placeholders})',
            conn,
            params=list(lig_names),
        )
        conn.close()
    except Exception:
        return {}
    return dict(zip(df["inchi_key"], df["name"]))


def add_ligand_name_column(df, project_path, assay_name, ligname_col="LigName"):
    """
    Insert a 'name' column right after *ligname_col*, mapping each LigName
    (chemspace inchi_key) to its compound name via the chemspace table linked
    to *assay_name*. Returns *df* unchanged if it is None/empty, lacks
    *ligname_col*, or the chemspace lookup cannot be resolved.
    """
    if df is None or df.empty or ligname_col not in df.columns:
        return df
    name_map = get_ligand_name_map(
        project_path, assay_name, df[ligname_col].dropna().unique().tolist()
    )
    if not name_map:
        return df
    df = df.copy()
    df.insert(df.columns.get_loc(ligname_col) + 1, "name", df[ligname_col].map(name_map))
    return df


def create_chemspace_subset(project_path, assay_name, inchi_keys, new_table_name):
    """
    Create a new table in chemspace.db as a subset of the compound table used for
    *assay_name*, copying ALL columns so the new table is ready for docking.

    Returns (True, success_message) or (False, error_message).
    """
    if not inchi_keys:
        return False, "No compounds to subset."
    if not new_table_name or not new_table_name.strip():
        return False, "Table name cannot be empty."
    new_table_name = new_table_name.strip().replace(" ", "_")

    # Resolve source table from docking registry
    registry_db = os.path.join(
        project_path, "docking", "docking_registers", "docking_assays.db"
    )
    if not os.path.exists(registry_db):
        return False, f"Docking registry not found: {registry_db}"
    try:
        conn = sqlite3.connect(registry_db)
        row = conn.execute(
            "SELECT table_name FROM docking_assays WHERE assay_name = ?", (assay_name,)
        ).fetchone()
        conn.close()
    except Exception as e:
        return False, f"Error reading docking registry: {e}"
    if not row:
        return False, f"Assay '{assay_name}' not found in docking registry."
    source_table = row[0]

    # Create subset in chemspace.db
    chemspace_db = os.path.join(
        project_path, "chemspace", "processed_data", "chemspace.db"
    )
    if not os.path.exists(chemspace_db):
        return False, f"ChemSpace database not found: {chemspace_db}"
    try:
        conn = sqlite3.connect(chemspace_db)
        cursor = conn.cursor()
        cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (new_table_name,),
        )
        if cursor.fetchone():
            conn.close()
            return False, f"Table '{new_table_name}' already exists in ChemSpace database."
        placeholders = ",".join("?" * len(inchi_keys))
        cursor.execute(
            f'CREATE TABLE "{new_table_name}" AS '
            f'SELECT * FROM "{source_table}" WHERE inchi_key IN ({placeholders})',
            list(inchi_keys),
        )
        count = cursor.execute(
            f'SELECT COUNT(*) FROM "{new_table_name}"'
        ).fetchone()[0]
        conn.commit()
        conn.close()
        return True, f"Created '{new_table_name}' with {count} compound(s)."
    except Exception as e:
        return False, f"Error creating subset table: {e}"


def get_rf_prediction_tables(results_db_path):
    """
    Return a list of model-name suffixes for every rf_predictions_* table
    present in the results DB (e.g. ['MyModel'] for 'rf_predictions_MyModel').
    Returns an empty list when the DB doesn't exist or has no such tables.
    """
    if not results_db_path or not os.path.exists(results_db_path):
        return []
    try:
        conn = sqlite3.connect(results_db_path)
        rows = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'rf_predictions_%'"
        ).fetchall()
        conn.close()
        prefix = "rf_predictions_"
        return [r[0][len(prefix):] for r in rows]
    except Exception:
        return []


def get_rf_classified_poses(results_db_path, model_suffix):
    """
    Join rf_predictions_{model_suffix} with the Results table on Pose_ID.

    Returns {"positive": DataFrame, "negative": DataFrame}.
    Positive poses are sorted by rf_prob_positive DESC; negative ASC.
    Columns: Pose_ID, LigName, run_number, docking_score, rf_prob_positive.
    Also includes 'mmgbsa_total_energy' and 'mmgbsa_gas_energy' when present,
    mirroring get_docking_results().
    """
    table = f"rf_predictions_{model_suffix}"
    empty = {"positive": pd.DataFrame(), "negative": pd.DataFrame()}
    if not os.path.exists(results_db_path):
        return empty
    try:
        conn = sqlite3.connect(results_db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(Results)")
        existing = {row[1] for row in cursor.fetchall()}
        optional_columns = [c for c in ('mmgbsa_total_energy', 'mmgbsa_gas_energy') if c in existing]
        r_columns = ", ".join(f"R.{c}" for c in ['Pose_ID', 'LigName', 'run_number', 'docking_score'] + optional_columns)
        df_pos = pd.read_sql_query(
            f"""
            SELECT {r_columns}, P.rf_prob_positive
            FROM Results AS R
            JOIN "{table}" AS P ON R.Pose_ID = P.Pose_ID
            WHERE P.rf_label = 1
            ORDER BY P.rf_prob_positive DESC
            """,
            conn,
        )
        df_neg = pd.read_sql_query(
            f"""
            SELECT {r_columns}, P.rf_prob_positive
            FROM Results AS R
            JOIN "{table}" AS P ON R.Pose_ID = P.Pose_ID
            WHERE P.rf_label = 0
            ORDER BY P.rf_prob_positive ASC
            """,
            conn,
        )
        conn.close()
        return {"positive": df_pos, "negative": df_neg}
    except Exception:
        return empty


def get_rf_ligand_summary(results_db_path, model_suffix):
    """
    Aggregate RF pose-level predictions per LigName.

    A ligand can have multiple docked poses (distinct run_number values), each
    classified independently by the RF model — a ligand is not guaranteed to
    get the same predicted label across all of its poses. This summarizes,
    per LigName, how many poses were predicted positive/negative.

    Returns a DataFrame with columns:
        LigName, n_poses, n_positive, n_negative, frac_positive,
        consensus_label, best_docking_score, max_rf_prob_positive
    sorted by frac_positive DESC, then n_poses DESC.
    consensus_label is 1 when frac_positive >= 0.5, else 0.
    best_docking_score is the minimum (most favorable) docking_score among
    the ligand's poses.
    Returns an empty DataFrame if the DB/table is missing or has no rows.
    """
    table = f"rf_predictions_{model_suffix}"
    if not os.path.exists(results_db_path):
        return pd.DataFrame()
    try:
        conn = sqlite3.connect(results_db_path)
        df = pd.read_sql_query(
            f"""
            SELECT R.LigName, R.run_number, R.docking_score,
                   P.rf_label, P.rf_prob_positive
            FROM Results AS R
            JOIN "{table}" AS P ON R.Pose_ID = P.Pose_ID
            """,
            conn,
        )
        conn.close()
    except Exception:
        return pd.DataFrame()

    if df.empty:
        return df

    summary = df.groupby("LigName").agg(
        n_poses=("rf_label", "size"),
        n_positive=("rf_label", "sum"),
        best_docking_score=("docking_score", "min"),
        max_rf_prob_positive=("rf_prob_positive", "max"),
    ).reset_index()
    summary["n_positive"] = summary["n_positive"].astype(int)
    summary["n_negative"] = summary["n_poses"] - summary["n_positive"]
    summary["frac_positive"] = summary["n_positive"] / summary["n_poses"]
    summary["consensus_label"] = (summary["frac_positive"] >= 0.5).astype(int)
    summary = summary[
        ["LigName", "n_poses", "n_positive", "n_negative", "frac_positive",
         "consensus_label", "best_docking_score", "max_rf_prob_positive"]
    ].sort_values(["frac_positive", "n_poses"], ascending=False).reset_index(drop=True)
    return summary


def get_rf_compound_summary(results_db_path, model_suffix, project_path=None, assay_name=None):
    """
    Aggregate the per-ligand RF summary further by chemical compound,
    grouping stereoisomers together.

    LigName values are InChIKeys; the first hyphen-delimited block encodes
    the molecular skeleton (connectivity) while the second block encodes
    stereochemistry/isotopes. Two LigNames sharing the first block are the
    same compound in different stereoisomeric forms. This groups the
    per-ligand summary (see get_rf_ligand_summary) by that shared prefix
    and reports averages across the grouped stereoisomers.

    If *project_path* and *assay_name* are given, a 'compound_name' column
    is also included, derived from each LigName's chemspace 'name' (see
    get_ligand_name_map) with the '_stereo_<n>' stereoisomer-enumeration
    suffix stripped (e.g. 'aspirin_stereo_2' -> 'aspirin'). All stereoisomers
    grouped under the same compound_key are expected to share this base
    name; the first non-null value found in the group is used.

    Returns a DataFrame with columns:
        compound_key, [compound_name,] n_representatives, avg_n_poses,
        avg_n_positive, avg_n_negative, avg_frac_positive,
        avg_consensus_label, avg_best_docking_score, avg_max_rf_prob_positive
    where n_representatives is the number of distinct stereoisomer LigNames
    grouped under compound_key, and the avg_* columns are the mean of the
    corresponding per-ligand summary column across those LigNames.
    sorted by avg_frac_positive DESC, then n_representatives DESC.
    Returns an empty DataFrame if no per-ligand summary is available.
    """
    lig_summary = get_rf_ligand_summary(results_db_path, model_suffix)
    if lig_summary.empty:
        return lig_summary

    df = lig_summary.copy()
    df["compound_key"] = df["LigName"].str.split("-").str[0]

    if project_path and assay_name:
        name_map = get_ligand_name_map(
            project_path, assay_name, df["LigName"].dropna().unique().tolist()
        )
        df["compound_name"] = df["LigName"].map(name_map).str.replace(
            r"_stereo_\d+$", "", regex=True
        )

    avg_cols = ["n_poses", "n_positive", "n_negative", "frac_positive",
                "consensus_label", "best_docking_score", "max_rf_prob_positive"]

    agg_kwargs = {"n_representatives": ("LigName", "nunique")}
    if "compound_name" in df.columns:
        agg_kwargs["compound_name"] = ("compound_name", "first")
    agg_kwargs.update({f"avg_{col}": (col, "mean") for col in avg_cols})

    compound_summary = df.groupby("compound_key").agg(**agg_kwargs).reset_index()

    ordered_cols = ["compound_key"]
    if "compound_name" in compound_summary.columns:
        ordered_cols.append("compound_name")
    ordered_cols.append("n_representatives")
    ordered_cols += [f"avg_{col}" for col in avg_cols]

    compound_summary = compound_summary[ordered_cols].sort_values(
        ["avg_frac_positive", "n_representatives"], ascending=False
    ).reset_index(drop=True)

    return compound_summary


def find_pose_pdb(results_db_path, ligname, run_number):
    """
    Search extracted pose directories adjacent to results_db_path for
    {ligname}_{run_number}.pdb.  Returns the full path on the first hit,
    or None if the file is not found in any directory.
    """
    results_dir = os.path.dirname(os.path.abspath(results_db_path))
    filename = f"{ligname}_{run_number}.pdb"
    for subdir in ("most_stable_poses", "most_populated_poses",
                   "most_populated_and_stable_poses", "all_poses"):
        candidate = os.path.join(results_dir, subdir, filename)
        if os.path.exists(candidate):
            return candidate
    return None


def resolve_ligpose_to_pdb(project_path, ligpose):
    """
    Resolve a 'ligpose' identifier (format '{pose_file_stem}_{assay_name}',
    produced by export_training_set_fingerprints_as_csv_bytes /
    MachineLearning.export_training_set_fingerprints_csv) back to an actual
    pose PDB file on disk.

    The merge is ambiguous in general (both the pose stem and the assay name
    may contain underscores), so this tries every known assay_name in the
    project as a suffix match, longest name first to minimise false matches.

    Returns (pdb_path, assay_name, ligname, run_number), or
    (None, None, None, None) if the identifier can't be resolved (unknown
    assay, malformed string, or the pose hasn't been extracted to disk yet).
    """
    registry_db = os.path.join(project_path, "docking", "docking_registers", "docking_assays.db")
    if not ligpose or not os.path.exists(registry_db):
        return None, None, None, None
    try:
        conn = sqlite3.connect(registry_db)
        assay_names = [r[0] for r in conn.execute("SELECT assay_name FROM docking_assays").fetchall()]
        conn.close()
    except Exception:
        return None, None, None, None

    for assay_name in sorted(assay_names, key=len, reverse=True):
        suffix = f"_{assay_name}"
        if not ligpose.endswith(suffix):
            continue
        pose_stem = ligpose[: -len(suffix)]
        parts = pose_stem.rsplit("_", 1)
        if len(parts) != 2:
            continue
        ligname, run_number_str = parts
        try:
            run_number = int(run_number_str)
        except ValueError:
            continue
        results_db_path = os.path.join(
            project_path, "docking", "docking_assays", assay_name, "results", f"{assay_name}.db"
        )
        pdb_path = find_pose_pdb(results_db_path, ligname, run_number)
        if pdb_path:
            return pdb_path, assay_name, ligname, run_number
    return None, None, None, None


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
        query_pdb_models = """
            SELECT m.file_id, m.pdb_model_name, m.project_name,
                   t.checked_pdb_path AS pdb_receptor_path,
                   m.filename, m.description, m.notes
            FROM pdb_models m
            LEFT JOIN pdb_templates t ON t.pdb_model_name = m.pdb_model_name
        """
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


def _rebase_pose_path(project_path: str, stored_path: str, assay_name: str, directory: str, pose_file: str) -> str:
    """Return stored_path if it exists on disk; otherwise reconstruct from project structure."""
    if os.path.exists(stored_path):
        return stored_path
    rebased = os.path.join(project_path, "docking", "docking_assays", assay_name, "results", directory, pose_file)
    return rebased


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
        df["pose_full_path"] = df.apply(
            lambda r: _rebase_pose_path(project_path, r["pose_full_path"], r["assay_name"], r["directory"], r["pose_file"]),
            axis=1,
        )
        return df
    except Exception:
        return None


def get_binders_registry_with_scores(project_path: str, binder_type: str, assay_name: str, results_db_path: str) -> "pd.DataFrame":
    """
    Like get_binders_registry(), but filtered to *assay_name* and enriched with
    docking_score (and mmgbsa_total_energy/mmgbsa_gas_energy when present) looked
    up from the assay's Results table.

    Binder registry rows only store pose_file/directory, not Pose_ID, so each
    pose_file ('{LigName}_{run_number}.pdb', the naming convention used by
    extract_docked_poses()/find_pose_pdb()) is parsed back into LigName/run_number
    and merged against Results on those columns. Rows whose pose_file can't be
    parsed, or that have no matching Results row, keep NaN for the score columns.

    Returns None if the registry itself is missing; an empty/unmodified DataFrame
    if there are no rows for this assay or the Results table can't be read.
    """
    df = get_binders_registry(project_path, binder_type)
    if df is None:
        return None
    if df.empty:
        return df

    df = df[df["assay_name"] == assay_name].reset_index(drop=True)
    if df.empty or not results_db_path or not os.path.exists(results_db_path):
        return df

    try:
        conn = sqlite3.connect(results_db_path)
        cursor = conn.cursor()
        cursor.execute("PRAGMA table_info(Results)")
        existing = {row[1] for row in cursor.fetchall()}
        score_columns = ['docking_score'] + [c for c in ('mmgbsa_total_energy', 'mmgbsa_gas_energy') if c in existing]
        results_df = pd.read_sql_query(
            f"SELECT LigName, run_number, {', '.join(score_columns)} FROM Results", conn
        )
        conn.close()
    except Exception:
        return df

    def _parse_pose_file(pose_file):
        stem = pose_file[:-4] if pose_file.lower().endswith(".pdb") else pose_file
        parts = stem.rsplit("_", 1)
        if len(parts) != 2:
            return None, None
        ligname, run_number_str = parts
        try:
            return ligname, int(run_number_str)
        except ValueError:
            return None, None

    parsed = df["pose_file"].apply(_parse_pose_file)
    df["LigName"] = [p[0] for p in parsed]
    df["run_number"] = [p[1] for p in parsed]

    df = df.merge(results_df, on=["LigName", "run_number"], how="left")
    return df


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


def move_binder_to_negative(project_path: str, assay_name: str, pose_file: str, directory: str, pose_full_path: str) -> str:
    """
    Move a pose from the positive binders registry to the negative binders registry.

    Removes the (assay_name, pose_file, directory) row from positive_binders.db first,
    then inserts it into negative_binders.db — removing first means
    save_negative_binder()'s own "already a positive binder" conflict check no
    longer sees the row, so the insert isn't blocked.

    Returns 'moved', 'not_found' (nothing to remove from positive), 'duplicate'
    (already present in negative_binders — still removed from positive), or
    'error:<message>'.
    """
    try:
        remove_status = remove_binder(project_path, "positive", assay_name, pose_file, directory)
        if remove_status == "not_found":
            return "not_found"
        if remove_status.startswith("error:"):
            return remove_status
        save_status = save_negative_binder(project_path, assay_name, pose_file, directory, pose_full_path)
        if save_status.startswith("error:"):
            return save_status
        return "moved" if save_status == "saved" else save_status
    except Exception as e:
        return f"error:{e}"


def move_binder_to_positive(project_path: str, assay_name: str, pose_file: str, directory: str, pose_full_path: str) -> str:
    """
    Move a pose from the negative binders registry to the positive binders registry.

    Mirror of move_binder_to_negative(): removes the (assay_name, pose_file,
    directory) row from negative_binders.db first, then inserts it into
    positive_binders.db — removing first means save_positive_binder()'s own
    "already a negative binder" conflict check no longer sees the row, so the
    insert isn't blocked.

    Returns 'moved', 'not_found' (nothing to remove from negative), 'duplicate'
    (already present in positive_binders — still removed from negative), or
    'error:<message>'.
    """
    try:
        remove_status = remove_binder(project_path, "negative", assay_name, pose_file, directory)
        if remove_status == "not_found":
            return "not_found"
        if remove_status.startswith("error:"):
            return remove_status
        save_status = save_positive_binder(project_path, assay_name, pose_file, directory, pose_full_path)
        if save_status.startswith("error:"):
            return save_status
        return "moved" if save_status == "saved" else save_status
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
        if df.empty:
            return None
        df["pose_full_path"] = df.apply(
            lambda r: _rebase_pose_path(project_path, r["pose_full_path"], r["assay_name"], r["directory"], r["pose_file"]),
            axis=1,
        )
        return df
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


def get_training_set_prolif_condition_options(project_path: str, training_set_id: str) -> "list[dict]":
    """
    Return the distinct ProLIF conditions for which fingerprints have been
    computed on a given training set snapshot.

    Args:
        project_path (str): Root path of the active project.
        training_set_id (str): ID of the training set snapshot.

    Returns:
        list[dict]: [{'id': <prolif_conditions_id>, 'description': <str or None>}, ...]
                    ordered by id. Empty list if none found.
    """
    db_path = os.path.join(project_path, "ml", "training_sets", "training_sets_snapshots.db")
    if not os.path.exists(db_path):
        return []
    try:
        conn = sqlite3.connect(db_path)
        rows = pd.read_sql_query(
            "SELECT DISTINCT prolif_conditions_id FROM training_set_fingerprints WHERE training_set_id = ? ORDER BY prolif_conditions_id",
            conn, params=(training_set_id,),
        )
        conn.close()
    except Exception:
        return []

    if rows.empty:
        return []

    ids = rows["prolif_conditions_id"].tolist()
    desc_map = {}
    params_db = os.path.join(project_path, "docking", "params", "params.db")
    try:
        if os.path.exists(params_db):
            conn = sqlite3.connect(params_db)
            conds = pd.read_sql_query("SELECT id, description FROM ProLIF_Conditions", conn)
            conn.close()
            desc_map = dict(zip(conds["id"], conds["description"]))
    except Exception:
        pass

    return [{"id": cid, "description": desc_map.get(cid)} for cid in ids]


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


@st.cache_data(show_spinner=False)
def get_table_columns(db_path: str, table_name: str, mtime=None) -> list:
    """
    Return the list of column names for a given table in a SQLite database.

    Cached per (db_path, table_name, mtime) so re-running the page (e.g. from
    toggling an unrelated checkbox) doesn't re-query the schema every time.

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


@st.cache_data(show_spinner=False)
def read_table_columns_as_dataframe(db_path: str, table_name: str, columns: list, limit: int = None, mtime=None) -> "pd.DataFrame":
    """
    Read selected columns from a SQLite table and return as a pandas DataFrame.

    Cached per (db_path, table_name, columns, limit, mtime): on large tables this
    query can be expensive, and Streamlit reruns the whole script on every widget
    interaction (e.g. toggling a column checkbox), so without caching the full
    table would be re-read from disk on every click.

    Args:
        db_path (str): Path to the SQLite database.
        table_name (str): Name of the table.
        columns (list): List of column names to include.
        limit (int, optional): If given, only read the first `limit` rows
            (pushed down as a SQL LIMIT) instead of loading the whole table.

    Returns:
        pd.DataFrame: DataFrame with the selected columns, or empty DataFrame on error.
    """
    try:
        col_expr = ", ".join(f"[{c}]" for c in columns)
        conn = sqlite3.connect(db_path)
        query = f"SELECT {col_expr} FROM [{table_name}]"
        if limit is not None:
            query += f" LIMIT {int(limit)}"
        df = pd.read_sql_query(query, conn)
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


def drop_chemspace_tables(db_path: str, table_names: list) -> dict:
    """
    Drop one or more tables from a SQLite database.

    Mirrors the SQL steps of ChemSpace.drop_table() (drop any indexes named
    after the table, then the table itself) without that method's
    input()/print() prompts, so it can be driven from the Streamlit UI with
    its own confirmation flow.

    Args:
        db_path (str): Path to the SQLite database.
        table_names (list[str]): Names of the tables to drop.

    Returns:
        dict[str, str]: {table_name: 'dropped'} on success, or
            {table_name: 'error:<message>'} for tables that failed.

    Connects with isolation_level=None (autocommit) so each DROP INDEX/DROP
    TABLE commits immediately instead of accumulating in one open
    transaction, and closes the connection via `finally` so a lock can't be
    left behind if this is interrupted partway through (see
    get_table_signatures for why a bare try/except isn't enough).
    """
    results = {}
    conn = None
    try:
        conn = sqlite3.connect(db_path, isolation_level=None, timeout=10)
        cursor = conn.cursor()
        for table_name in table_names:
            try:
                cursor.execute(
                    "SELECT name FROM sqlite_master WHERE type='index' AND (name LIKE ? OR name LIKE ?)",
                    (f"idx_{table_name}_%", f"index_{table_name}_%")
                )
                for (index_name,) in cursor.fetchall():
                    cursor.execute(f"DROP INDEX IF EXISTS [{index_name}]")
                cursor.execute(f"DROP TABLE IF EXISTS [{table_name}]")
                results[table_name] = "dropped"
            except Exception as e:
                results[table_name] = f"error:{e}"
    except Exception as e:
        for table_name in table_names:
            results.setdefault(table_name, f"error:{e}")
    finally:
        if conn is not None:
            conn.close()
    return results


def depict_table_to_images(db_path: str, table_name: str,
                           max_molecules: int = 25,
                           molecules_per_image: int = 25,
                           mol_image_size: tuple = (300, 300),
                           legend_cols: "list[str] | str | None" = None) -> list:
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
        legend_cols (list[str] | str): Column(s) to use as molecule label. When multiple
                          columns are given, the label is composed by joining their
                          values with '_'. Defaults to 'id' if present, otherwise the
                          first column.

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

    # Determine legend column(s): use caller's choice (str or list), else 'id', else first column
    if isinstance(legend_cols, str):
        legend_cols = [legend_cols]
    effective_legend_cols = [c for c in (legend_cols or []) if c in subset.columns]
    if not effective_legend_cols:
        if 'id' in subset.columns:
            effective_legend_cols = ['id']
        else:
            effective_legend_cols = [subset.columns[0]]

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
            legend_val = '_'.join(str(getattr(row, c, '')) for c in effective_legend_cols)[:40]
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


def export_training_set_fingerprints_as_csv_bytes(project_path: str, training_set_id: str, prolif_conditions_id: int) -> "bytes | None":
    """
    Build a CSV of the precomputed fingerprints for the given training set snapshot
    and ProLIF conditions, and return it as UTF-8 bytes suitable for st.download_button.

    A pose may have multiple fingerprints if computed under different ProLIF
    conditions; prolif_conditions_id selects which one to export so that rows
    from incompatible feature spaces are never mixed in the same CSV.

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
               WHERE training_set_id = ? AND prolif_conditions_id = ?
               ORDER BY binder_type DESC, assay_name, pose_file""",
            conn, params=(training_set_id, prolif_conditions_id)
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
            """SELECT m.filename, m.pdb_blob, m.pdb_model_name, t.checked_pdb_path
               FROM pdb_models m
               LEFT JOIN pdb_templates t ON t.pdb_model_name = m.pdb_model_name
               WHERE m.file_id = ?""",
            (int(file_id),)
        )
        row = cursor.fetchone()
        conn.close()

        if row is None:
            return False, f"No PDB model found for file_id '{file_id}'"

        filename, pdb_blob, pdb_model_name, checked_pdb_path = row

        # Stub receptors (created by import_receptor_model) have an empty blob and no
        # filename. Fall back to the receptor_checked.pdb stored on disk.
        if not pdb_blob:
            if checked_pdb_path and os.path.isfile(checked_pdb_path):
                fallback_name = filename or f"{pdb_model_name}.pdb"
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, fallback_name)
                import shutil
                shutil.copy2(checked_pdb_path, output_path)
                return True, output_path
            return False, (
                f"No PDB data available for stub receptor '{pdb_model_name}'. "
                "The receptor was imported without a PDB file."
            )

        fallback_name = filename or f"{pdb_model_name}.pdb"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, fallback_name)

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


def ensure_rf_trained_models_schema(db_path):
    """
    Create/migrate the rf_trained_models table so it carries everything needed
    to fully reconstruct the RF training "Results" tab for a stored model, not
    just the summary scalars (roc_auc/accuracy/macro_f1/cv_roc_mean/cv_roc_std).

    sklearn/numpy objects aren't natively storable in SQLite, so the confusion
    matrix, classification report, feature importances, per-fold CV scores and
    GridSearchCV best params are persisted as JSON text columns. Models saved
    before this migration will have NULL in these columns; callers should treat
    that as "detailed results unavailable for this model" rather than an error.
    """
    conn = sqlite3.connect(db_path)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS rf_trained_models (
            model_id                    INTEGER PRIMARY KEY AUTOINCREMENT,
            model_name                  TEXT NOT NULL,
            description                 TEXT,
            training_set_id             TEXT,
            roc_auc                     REAL,
            accuracy                    REAL,
            macro_f1                    REAL,
            cv_roc_mean                 REAL,
            cv_roc_std                  REAL,
            model_pkl                   BLOB NOT NULL,
            created_at                  TEXT NOT NULL,
            confusion_matrix_json       TEXT,
            classification_report_json  TEXT,
            feature_importances_json    TEXT,
            cv_scores_json              TEXT,
            best_params_json            TEXT,
            best_cv_score                REAL,
            test_predictions_json       TEXT
        )
    """)
    existing_cols = {row[1] for row in conn.execute("PRAGMA table_info(rf_trained_models)").fetchall()}
    migrations = {
        "confusion_matrix_json":      "TEXT",
        "classification_report_json": "TEXT",
        "feature_importances_json":   "TEXT",
        "cv_scores_json":             "TEXT",
        "best_params_json":           "TEXT",
        "best_cv_score":              "REAL",
        "test_predictions_json":      "TEXT",
    }
    for col, col_type in migrations.items():
        if col not in existing_cols:
            conn.execute(f"ALTER TABLE rf_trained_models ADD COLUMN {col} {col_type}")
    conn.commit()
    conn.close()
