from tidyscreen import tidyscreen
import sys
import os
import contextlib
import sqlite3
import json
import pandas as pd
from typing import Optional
from tqdm import tqdm

# Add the parent directory to path to import our local tidyscreen module
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir)

# Import from our local tidyscreen module
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject
from tidyscreen.actionlog.action_logger import log_all_public_methods


def _predict_pose_worker(task: dict) -> dict:
    """
    Module-level worker (required for ProcessPoolExecutor pickling) that restores
    a single docked pose, builds its AMBER complex, optionally minimizes it, and
    computes its ProLIF fingerprint.

    Runs in its own process, so it re-derives a MolDock instance via
    MolDock.from_path() rather than receiving a live instance. Ligand-specific
    tleap files (prepin/frcmod) must already be prepared (see
    MachineLearning._classify_poses_parallel) and are passed in read-only, since
    concurrent workers preparing the same ligand would race on its files.

    Each pose gets its own scratch subdirectory (temp_restored_poses/pose_<id>/)
    because MolDock's private helpers write fixed-name tleap/sander input files
    (complex.in, min.in, minimize.in, ...) into the directory they're given —
    fine when one pose is processed at a time, but colliding if two poses shared
    a directory concurrently.
    """
    import os
    import sys
    import shutil
    from moldock.moldock import MolDock

    pose_id = task['pose_id']
    ligname = task['ligname']
    verbose = task['verbose']

    devnull_fd = saved_out_fd = saved_err_fd = None
    if not verbose:
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        saved_out_fd = os.dup(1)
        saved_err_fd = os.dup(2)
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)

    pose_dir = None
    try:
        moldock = MolDock.from_path(task['project_name'], task['project_path'])

        shared_dir, restored_pdb = moldock._restore_single_docked_pose(
            task['results_db'], ligname, task['run_number']
        )
        pose_dir = os.path.join(shared_dir, f"pose_{pose_id}")
        work_dir = os.path.join(pose_dir, "work")
        os.makedirs(work_dir, exist_ok=True)

        prmtop_file, inpcrd_file, output_pdb_file = moldock._prepare_complex_prmtop_inpcrd(
            task['prepin_file'], task['frcmod_file'], task['assay_info'],
            work_dir, restored_pdb, pose_id, ligname=ligname
        )

        if task['minimize']:
            _, output_pdb_file = moldock._minimize_complex(
                prmtop_file, inpcrd_file, work_dir, ligname, pose_id, output_pdb_file
            )

        fps_df = moldock._compute_prolif_fingerprints3(
            task['prolif_params_dict'], output_pdb_file, task['receptor_file'],
            prmtop_file, inpcrd_file, ligname, pose_id
        )
        fps_df.columns = fps_df.columns.droplevel(0)
        if task['renumbering_dict']:
            fps_df = fps_df.rename(columns=task['renumbering_dict'], level=0)
        fps_df.columns = [f"{a}_{b}" for a, b in fps_df.columns]

        fp_row = fps_df.iloc[0].to_dict()
        return {
            'ok': True,
            'pose_id': pose_id,
            'lig_name': ligname,
            'docking_score': task['docking_score'],
            **{k: bool(v) for k, v in fp_row.items()},
        }
    except Exception as e:
        return {'ok': False, 'pose_id': pose_id, 'lig_name': ligname, 'error': str(e)}
    finally:
        if not verbose:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_out_fd, 1)
            os.dup2(saved_err_fd, 2)
            os.close(saved_out_fd)
            os.close(saved_err_fd)
            os.close(devnull_fd)
        if pose_dir and os.path.isdir(pose_dir):
            shutil.rmtree(pose_dir, ignore_errors=True)


def _run_background_rf_predictions(job_file: str) -> None:
    """
    Entry point for the detached background process launched by
    MachineLearning._launch_background_rf_predictions(). Runs as a separate
    `python -c` process with no connection to the interactive session that
    launched it, so everything needed is reloaded from the job spec file
    (plain, JSON-serializable data) rather than received as live objects.
    """
    import json
    import pickle

    with open(job_file) as f:
        job = json.load(f)

    project = ActivateProject(job['project_name'])
    ml = MachineLearning(project)

    from moldock.moldock import MolDock
    moldock = MolDock.from_path(job['project_name'], job['project_path'])

    conn = sqlite3.connect(job['rf_db_path'])
    row = conn.execute(
        "SELECT model_pkl FROM rf_trained_models WHERE model_id = ?", (job['model_id'],)
    ).fetchone()
    conn.close()
    if not row or not row[0]:
        print(f"[ERROR] Model blob missing for model_id {job['model_id']}")
        return
    model = pickle.loads(row[0])

    conn = sqlite3.connect(job['results_db'])
    poses_df = pd.read_sql_query(
        "SELECT Pose_ID, LigName, docking_score, run_number FROM Results", conn
    )
    conn.close()

    print(f"🧬 Background classification started for {len(poses_df)} poses "
          f"(job: {os.path.basename(job_file)})")

    fp_records, n_ok, n_failed = ml._classify_poses_parallel(
        moldock, poses_df, job['results_db'], job['assay_info'],
        job['prolif_params_dict'], job['receptor_file'], job['renumbering_dict'],
        job['minimize'], job['max_workers'], job['verbose']
    )

    print(f"\n🎉 Background classification complete: {n_ok} succeeded, {n_failed} failed.")

    ml._align_and_store_rf_predictions(model, fp_records, job['results_db'], job['predictions_table'])


@log_all_public_methods
class MachineLearning:
    """
    MachineLearning class for managing machine learning workflows within a project.
    Uses an ActivateProject object to access project information and database functionality.
    """

    def __init__(self, project_obj: ActivateProject):
        """
        Initialize MachineLearning with an ActivateProject object.

        Args:
            project_obj (ActivateProject): An instantiated ActivateProject object
        """
        if not isinstance(project_obj, ActivateProject):
            raise TypeError("MachineLearning requires an ActivateProject object")

        if not project_obj.project_exists():
            raise ValueError(f"Project '{project_obj.name}' not found. Please create the project first.")

        self.project = project_obj
        self.name = project_obj.name
        self.path = project_obj.path
        self.description = getattr(project_obj, 'description', None)
        self.id = getattr(project_obj, 'id', None)
        self.created_date = getattr(project_obj, 'created_date', None)

        # Set up chemspace database path within the project directory
        self.__chemspace_db = os.path.join(self.path, 'chemspace/processed_data', 'chemspace.db')

        # Set up ML results database path within the project directory
        self.__ml_db = os.path.join(self.path, 'ml', 'ml_results.db')

        # Set up ML models storage path within the project directory
        self.__ml_models_path = os.path.join(self.path, 'ml', 'models')

        # Set up ML results path within the project directory
        self.__ml_results_path = os.path.join(self.path, 'ml', 'results')

        # Set up ML training sets path within the project directory
        self.__training_sets_db = os.path.join(self.path, 'ml', 'training_sets', 'training_sets_snapshots.db')
        self.__docking_params_db = os.path.join(self.path, 'docking', 'params', 'params.db')
        self.__docking_assays_db = os.path.join(self.path, 'docking', 'docking_registers', 'docking_assays.db')
        self.__pdbs_db = os.path.join(self.path, 'docking', 'receptors', 'pdbs.db')

        # Ensure ML directories exist
        os.makedirs(self.__ml_models_path, exist_ok=True)
        os.makedirs(self.__ml_results_path, exist_ok=True)

    # -------------------------------------------------------------------------
    # Internal helpers
    # -------------------------------------------------------------------------

    def __get_connection(self, db_path: str) -> sqlite3.Connection:
        """Return a sqlite3 connection to the given database path."""
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Database not found: {db_path}")
        return sqlite3.connect(db_path)

    # -------------------------------------------------------------------------
    # Feature preparation
    # -------------------------------------------------------------------------

    def prepare_features(self, table_name: str, feature_columns: list, label_column: str) -> Optional[pd.DataFrame]:
        """
        Load features and labels from the chemspace database for model training.

        Args:
            table_name (str): Name of the table in the chemspace database.
            feature_columns (list): List of column names to use as features.
            label_column (str): Column name to use as label/target.

        Returns:
            pd.DataFrame or None: DataFrame with features and label column, or None on error.
        """
        try:
            conn = self.__get_connection(self.__chemspace_db)
            cols = ", ".join(feature_columns + [label_column])
            df = pd.read_sql_query(f"SELECT {cols} FROM {table_name}", conn)
            conn.close()
            print(f"\n[MachineLearning] Loaded {len(df)} rows from '{table_name}'.")
            return df
        except Exception as e:
            print(f"\n[MachineLearning] Error preparing features: {e}")
            return None

    # -------------------------------------------------------------------------
    # Model training
    # -------------------------------------------------------------------------

    def train_model(self, model, X_train, y_train, model_name: str) -> None:
        """
        Train a scikit-learn compatible model and save it to the models directory.

        Args:
            model: A scikit-learn compatible estimator.
            X_train: Training features.
            y_train: Training labels.
            model_name (str): Name used to save the model file (without extension).
        """
        try:
            import joblib
            print(f"\n[MachineLearning] Training model '{model_name}'...")
            model.fit(X_train, y_train)
            model_path = os.path.join(self.__ml_models_path, f"{model_name}.joblib")
            joblib.dump(model, model_path)
            print(f"[MachineLearning] Model saved to: {model_path}")
        except Exception as e:
            print(f"\n[MachineLearning] Error training model: {e}")

    # -------------------------------------------------------------------------
    # Prediction
    # -------------------------------------------------------------------------

    def predict(self, model_name: str, X) -> Optional[object]:
        """
        Load a saved model and run predictions.

        Args:
            model_name (str): Name of the saved model file (without extension).
            X: Input features for prediction.

        Returns:
            Predictions array or None on error.
        """
        try:
            import joblib
            model_path = os.path.join(self.__ml_models_path, f"{model_name}.joblib")
            if not os.path.exists(model_path):
                raise FileNotFoundError(f"Model not found: {model_path}")
            model = joblib.load(model_path)
            predictions = model.predict(X)
            print(f"\n[MachineLearning] Predictions generated using model '{model_name}'.")
            return predictions
        except Exception as e:
            print(f"\n[MachineLearning] Error during prediction: {e}")
            return None

    # -------------------------------------------------------------------------
    # Evaluation
    # -------------------------------------------------------------------------

    def evaluate(self, model_name: str, X_test, y_test) -> Optional[dict]:
        """
        Evaluate a saved model on test data and return common metrics.

        Args:
            model_name (str): Name of the saved model file (without extension).
            X_test: Test features.
            y_test: True labels.

        Returns:
            dict with evaluation metrics, or None on error.
        """
        try:
            from sklearn.metrics import accuracy_score, roc_auc_score, classification_report
            predictions = self.predict(model_name, X_test)
            if predictions is None:
                return None
            metrics = {
                "accuracy": accuracy_score(y_test, predictions),
                "report": classification_report(y_test, predictions),
            }
            try:
                metrics["roc_auc"] = roc_auc_score(y_test, predictions)
            except Exception:
                pass
            print(f"\n[MachineLearning] Evaluation results for '{model_name}':")
            for k, v in metrics.items():
                print(f"  {k}: {v}")
            return metrics
        except Exception as e:
            print(f"\n[MachineLearning] Error during evaluation: {e}")
            return None

    # -------------------------------------------------------------------------
    # Results persistence
    # -------------------------------------------------------------------------

    def save_results(self, df: pd.DataFrame, table_name: str) -> None:
        """
        Save a results DataFrame to the ML results database.

        Args:
            df (pd.DataFrame): DataFrame containing results to store.
            table_name (str): Target table name in the ML results database.
        """
        try:
            conn = sqlite3.connect(self.__ml_db)
            df.to_sql(table_name, conn, if_exists='replace', index=False)
            conn.close()
            print(f"\n[MachineLearning] Results saved to table '{table_name}' in '{self.__ml_db}'.")
        except Exception as e:
            print(f"\n[MachineLearning] Error saving results: {e}")

    def load_results(self, table_name: str) -> Optional[pd.DataFrame]:
        """
        Load results from the ML results database.

        Args:
            table_name (str): Table name to load from the ML results database.

        Returns:
            pd.DataFrame or None on error.
        """
        try:
            conn = self.__get_connection(self.__ml_db)
            df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
            conn.close()
            print(f"\n[MachineLearning] Loaded {len(df)} rows from table '{table_name}'.")
            return df
        except Exception as e:
            print(f"\n[MachineLearning] Error loading results: {e}")
            return None

    # -------------------------------------------------------------------------
    # ProLIF fingerprint computation on training set snapshots
    # -------------------------------------------------------------------------

    def compute_training_set_fingerprints(self, minimize=True, clean_files=True, verbose: bool = False):
        """
        Compute ProLIF fingerprints for all poses in a consolidated training set
        snapshot, following the same AMBER + ProLIF pipeline as
        MolDock.compute_fingerprints().

        The user is prompted to:
          1. Select a training set snapshot.
          2. Select a ProLIF conditions record.
          3. Confirm whether complex minimization should be performed.

        Results are stored in the 'training_set_fingerprints' table of
        training_sets_snapshots.db, keyed on
        (training_set_id, assay_name, pose_file, directory, prolif_conditions_id).

        The mapping from pose_file → LigName + run_number relies on the filename
        convention  ``{LigName}_{run_number}.pdb`` written by MolDock during pose
        extraction.  Note: the suffix is the per-ligand run_number, NOT the global
        Pose_ID from the Ringtail Results table.
        The assay_name column is used to look up assay metadata (receptor path,
        renumbering dict, results DB path) from the docking registers.

        Args:
            verbose (bool): If True, print full per-pose processing details,
                including output from the tleap/sander/ambpdb subprocess calls
                made while preparing each pose. If False (default), all of that
                output is suppressed and a tqdm progress counter is shown instead
                (same behaviour as make_rf_predictions_on_docking_assay()).
        """

        import shutil
        from datetime import datetime

        from moldock.moldock import MolDock

        # ── 1. Select training set ────────────────────────────────────────────
        snapshots = self._list_training_set_snapshots()
        if not snapshots:
            print("\n❌ No training set snapshots found. Consolidate a training set first.")
            return

        print(f"\n📦 AVAILABLE TRAINING SET SNAPSHOTS")
        print("=" * 70)
        for i, snap in enumerate(snapshots, 1):
            print(f"   {i}. {snap['training_set_id']}  "
                  f"(+{snap['n_positives']} / -{snap['n_negatives']})  "
                  f"created: {snap['created_at']}")
        print("=" * 70)

        while True:
            try:
                sel = input("\nSelect training set (number or 'c' to cancel): ").strip()
                if sel.lower() == 'c':
                    return
                idx = int(sel) - 1
                if 0 <= idx < len(snapshots):
                    selected_snap = snapshots[idx]
                    break
                print(f"❌ Enter a number between 1 and {len(snapshots)}")
            except ValueError:
                print("❌ Enter a valid number")

        training_set_id = selected_snap['training_set_id']
        print(f"\n✅ Selected training set: {training_set_id}")

        # ── 2. Select ProLIF conditions ───────────────────────────────────────
        prolif_params_dict, condition_selection = self._prompt_prolif_conditions()
        if prolif_params_dict is None:
            return

        # ── 3. Confirm minimization ───────────────────────────────────────────
        while True:
            ans = input("\nMinimize complex before fingerprint computation? (y/n): ").strip().lower()
            if ans in ('y', 'n'):
                minimize = ans == 'y'
                break
            print("❌ Please enter 'y' or 'n'")

        # ── 4. Load training set entries ──────────────────────────────────────
        conn = sqlite3.connect(self.__training_sets_db)
        entries_df = pd.read_sql_query(
            "SELECT * FROM training_set_entries WHERE training_set_id = ?",
            conn, params=(training_set_id,)
        )
        conn.close()

        if entries_df.empty:
            print(f"\n❌ No entries found for training set '{training_set_id}'")
            return

        print(f"\n📋 Loaded {len(entries_df)} entries "
              f"({(entries_df['binder_type'] == 'positive').sum()} positive, "
              f"{(entries_df['binder_type'] == 'negative').sum()} negative)")

        # ── 5. Load assay info keyed by assay_name ────────────────────────────
        assay_names = entries_df['assay_name'].unique().tolist()
        assay_info_by_name = self._load_assay_info_by_name(assay_names)
        if not assay_info_by_name:
            print("\n❌ Could not load assay information from the docking registry")
            return

        # ── 6. Load renumbering dicts keyed by assay_name ────────────────────
        renumbering_dict_by_name = {}
        for assay_name, assay_info in assay_info_by_name.items():
            template_name = assay_info.get('receptor_info', {}).get('template_name')
            if template_name:
                renumbering_dict_by_name[assay_name] = self._load_renumbering_dict(template_name)
            else:
                renumbering_dict_by_name[assay_name] = {}

        # ── 7. Ensure fingerprints table exists ───────────────────────────────
        self._ensure_fingerprints_table()

        # ── 8. Instantiate MolDock to access the AMBER + ProLIF pipeline ─────
        moldock = MolDock(self.project)

        output_dirs = []
        processed_ligands_by_assay = {}  # {assay_name: {ligname: (prepin_file, frcmod_file)}}
        n_computed = 0
        n_skipped = 0
        n_failed = 0

        # When quiet, tqdm writes through a duplicated stdout fd so the progress
        # counter keeps displaying even while fd 1/2 are redirected to devnull
        # below (needed to silence subprocess-inherited tleap/sander/ambpdb output).
        original_stdout = os.fdopen(os.dup(1), 'w') if not verbose else sys.stdout

        try:
            pbar = tqdm(
                entries_df.iterrows(), total=len(entries_df),
                desc="Computing fingerprints", file=original_stdout,
                dynamic_ncols=True, disable=verbose,
            )
            for _, row in pbar:
                assay_name = row['assay_name']
                pose_file  = row['pose_file']
                directory  = row['directory']
                binder_type = row['binder_type']
                label = 1 if binder_type == 'positive' else 0

                try:
                    # Skip if already computed for this conditions set
                    conn = sqlite3.connect(self.__training_sets_db)
                    existing = conn.execute(
                        """SELECT id FROM training_set_fingerprints
                           WHERE training_set_id=? AND assay_name=? AND pose_file=?
                             AND directory=? AND prolif_conditions_id=?""",
                        (training_set_id, assay_name, pose_file, directory, condition_selection)
                    ).fetchone()
                    conn.close()
                    if existing:
                        n_skipped += 1
                        if verbose:
                            print(f"  ⏭️  Already computed, skipping: {pose_file}")
                        continue

                    assay_info = assay_info_by_name.get(assay_name)
                    if not assay_info:
                        n_skipped += 1
                        if verbose:
                            print(f"  ❌ Assay '{assay_name}' not found in registry — skipping {pose_file}")
                        continue

                    # Parse LigName and run_number from pose_file ({LigName}_{run_number}.pdb).
                    # The suffix is the per-ligand run_number used during pose extraction,
                    # NOT the global Pose_ID from the Ringtail Results table.
                    stem = os.path.splitext(pose_file)[0]
                    parts = stem.rsplit('_', 1)
                    if len(parts) != 2:
                        n_skipped += 1
                        if verbose:
                            print(f"  ❌ Cannot parse pose_file '{pose_file}' — skipping")
                        continue
                    ligname, pose_id_str = parts
                    try:
                        pose_id = int(pose_id_str)  # this is the run_number
                    except ValueError:
                        n_skipped += 1
                        if verbose:
                            print(f"  ❌ Cannot parse run_number from '{pose_id_str}' — skipping {pose_file}")
                        continue

                    # Derive receptor_checked.pdb path from the pdbqt_file stored in assay info
                    receptor_pdbqt = assay_info.get('receptor_info', {}).get('pdbqt_file')
                    if not receptor_pdbqt:
                        n_skipped += 1
                        if verbose:
                            print(f"  ❌ No receptor pdbqt_file for assay '{assay_name}' — skipping {pose_file}")
                        continue
                    receptor_pdbqt = self._remap_project_path(receptor_pdbqt)
                    receptor_file = os.path.join(os.path.dirname(receptor_pdbqt), 'receptor_checked.pdb')

                    assay_folder_path = assay_info['assay_folder_path']
                    assay_id = assay_info['assay_id']
                    results_db = os.path.join(assay_folder_path, 'results', f'assay_{assay_id}.db')

                    renumbering_dict = renumbering_dict_by_name.get(assay_name, {})

                    if verbose:
                        print(f"\n  🔬 Processing {pose_file} ({binder_type}) — assay: {assay_name}")

                    suppress = contextlib.nullcontext() if verbose else self._suppress_output()
                    with suppress:
                        # Step A: Restore docked pose PDB from the assay results DB
                        output_dir, output_file = moldock._restore_single_docked_pose(results_db, ligname, pose_id)

                        if output_dir not in output_dirs:
                            output_dirs.append(output_dir)

                        # Step B: Prepare ligand AMBER files (once per unique ligname per assay)
                        processed_ligands_by_assay.setdefault(assay_name, {})
                        if ligname not in processed_ligands_by_assay[assay_name]:
                            prepin_file, frcmod_file = moldock._prepare_ligand_tleap_files(
                                ligname, assay_info, output_dir
                            )
                            processed_ligands_by_assay[assay_name][ligname] = (prepin_file, frcmod_file)
                        else:
                            prepin_file, frcmod_file = processed_ligands_by_assay[assay_name][ligname]

                        # Step C: Prepare complex prmtop and inpcrd files using tleap
                        prmtop_file, inpcrd_file, output_pdb_file = moldock._prepare_complex_prmtop_inpcrd(
                            prepin_file, frcmod_file, assay_info, output_dir, output_file, pose_id, ligname=ligname
                        )

                        # Step D: Optionally minimize complex (mirrors compute_fingerprints behaviour)
                        if minimize:
                            try:
                                if verbose:
                                    print("    MINIMIZING COMPLEX...")
                                min_rst_cpptraj_file, output_pdb_file = moldock._minimize_complex(
                                    prmtop_file, inpcrd_file, output_dir, ligname, pose_id, output_pdb_file
                                )
                                rst_cpptraj_file = min_rst_cpptraj_file
                            except Exception as e:
                                if verbose:
                                    print(f"    ⚠️  Minimization failed: {e}. Using original inpcrd.")
                                rst_cpptraj_file = inpcrd_file
                        else:
                            rst_cpptraj_file = inpcrd_file

                        # Step E: Compute ProLIF fingerprints via ambpdb + MDAnalysis + ProLIF
                        fps_df = moldock._compute_prolif_fingerprints3(
                            prolif_params_dict, output_pdb_file, receptor_file,
                            prmtop_file, inpcrd_file, ligname, pose_id
                        )

                        # Step F: Flatten MultiIndex columns and apply residue renumbering
                        fps_df.columns = fps_df.columns.droplevel(0)
                        if renumbering_dict:
                            fps_df = fps_df.rename(columns=renumbering_dict, level=0)
                        fps_df.columns = [f"{a}_{b}" for a, b in fps_df.columns]

                        # Step G: Serialise the single fingerprint row to JSON
                        fp_row = fps_df.iloc[0].to_dict()
                        fp_json = json.dumps({k: bool(v) for k, v in fp_row.items()})

                        # Step H: Store in training_set_fingerprints table
                        conn = sqlite3.connect(self.__training_sets_db)
                        conn.execute(
                            """INSERT OR REPLACE INTO training_set_fingerprints
                               (training_set_id, assay_name, pose_file, directory,
                                binder_type, label, prolif_conditions_id, minimized,
                                fingerprint_json, computed_at)
                               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                            (training_set_id, assay_name, pose_file, directory,
                             binder_type, label, condition_selection, int(minimize),
                             fp_json, datetime.now().isoformat())
                        )
                        conn.commit()
                        conn.close()

                    n_computed += 1
                    if verbose:
                        print(f"    ✅ Fingerprint stored for {pose_file}")

                except Exception as e:
                    n_failed += 1
                    if verbose:
                        print(f"    ❌ Error processing {pose_file}: {e} — skipping")
                    continue
                finally:
                    if not verbose:
                        pbar.set_postfix(computed=n_computed, skipped=n_skipped, failed=n_failed)

        finally:
            pbar.close()
            if not verbose:
                original_stdout.close()
            if clean_files:
                for d in output_dirs:
                    shutil.rmtree(d, ignore_errors=True)
                if verbose:
                    print("\n✅ Temporary files cleaned up.")

        print(f"\n🎉 Fingerprint computation complete for training set '{training_set_id}' "
              f"({n_computed} computed, {n_skipped} skipped, {n_failed} failed)")

    def load_training_set_features(self) -> tuple:
        """
        Interactively select a training set snapshot and a ProLIF conditions record,
        then return (X, y) ready for scikit-learn classifiers.

        X is a pd.DataFrame where each row is a pose fingerprint (binary int columns).
        y is a pd.Series of integer labels: 1 = positive binder, 0 = negative binder.

        Returns (None, None) if no fingerprints are found or the user cancels.
        """

        # ── Select training set ───────────────────────────────────────────────
        snapshots = self._list_training_set_snapshots()
        if not snapshots:
            print("\n❌ No training set snapshots found.")
            return None, None

        print(f"\n📦 AVAILABLE TRAINING SET SNAPSHOTS")
        print("=" * 70)
        for i, snap in enumerate(snapshots, 1):
            print(f"   {i}. {snap['training_set_id']}  "
                  f"(+{snap['n_positives']} / -{snap['n_negatives']})  "
                  f"created: {snap['created_at']}")
        print("=" * 70)

        while True:
            try:
                sel = input("\nSelect training set (number or 'c' to cancel): ").strip()
                if sel.lower() == 'c':
                    return None, None
                idx = int(sel) - 1
                if 0 <= idx < len(snapshots):
                    training_set_id = snapshots[idx]['training_set_id']
                    break
                print(f"❌ Enter a number between 1 and {len(snapshots)}")
            except ValueError:
                print("❌ Enter a valid number")

        # ── Select ProLIF conditions ──────────────────────────────────────────
        _, condition_selection = self._prompt_prolif_conditions()
        if condition_selection is None:
            return None, None

        # ── Load fingerprints ─────────────────────────────────────────────────
        if not os.path.exists(self.__training_sets_db):
            print("\n❌ Training sets database not found.")
            return None, None

        conn = sqlite3.connect(self.__training_sets_db)
        rows = pd.read_sql_query(
            """SELECT assay_name, pose_file, binder_type, label, fingerprint_json
               FROM training_set_fingerprints
               WHERE training_set_id = ? AND prolif_conditions_id = ?""",
            conn, params=(training_set_id, condition_selection)
        )
        conn.close()

        if rows.empty:
            print(f"\n❌ No fingerprints computed for training set '{training_set_id}' "
                  f"with conditions ID {condition_selection}. "
                  f"Run compute_training_set_fingerprints() first.")
            return None, None

        fp_dicts = [json.loads(fp) for fp in rows['fingerprint_json']]
        X = pd.DataFrame(fp_dicts).fillna(0).astype(int)
        y = rows['label'].reset_index(drop=True)

        print(f"\n[MachineLearning] Loaded {len(X)} fingerprints, "
              f"{X.shape[1]} features for training set '{training_set_id}' "
              f"(conditions ID {condition_selection})")
        print(f"  Label distribution: {y.value_counts().to_dict()}")
        return X, y

    def export_training_set_fingerprints_csv(self):
        """
        Export precomputed fingerprints for a training set snapshot to a CSV file.

        The user is prompted to select:
          1. A training set snapshot (filtered by training_set_id).
          2. A ProLIF conditions record.
          3. An output CSV file path.

        The exported CSV has the following column layout:
          - ``ligpose``  : ``{pose_file}__{assay_name}`` identifier for each pose.
          - <fp columns> : one binary integer column per interaction feature,
                           derived from ``fingerprint_json``.
          - ``label``    : 1 = positive binder, 0 = negative binder.
        """

        # ── Select training set ───────────────────────────────────────────────
        snapshots = self._list_training_set_snapshots()
        if not snapshots:
            print("\n❌ No training set snapshots found. Consolidate a training set first.")
            return

        print(f"\n📦 AVAILABLE TRAINING SET SNAPSHOTS")
        print("=" * 70)
        for i, snap in enumerate(snapshots, 1):
            print(f"   {i}. {snap['training_set_id']}  "
                  f"(+{snap['n_positives']} / -{snap['n_negatives']})  "
                  f"created: {snap['created_at']}")
        print("=" * 70)

        while True:
            try:
                sel = input("\nSelect training set (number or 'c' to cancel): ").strip()
                if sel.lower() == 'c':
                    return
                idx = int(sel) - 1
                if 0 <= idx < len(snapshots):
                    training_set_id = snapshots[idx]['training_set_id']
                    break
                print(f"❌ Enter a number between 1 and {len(snapshots)}")
            except ValueError:
                print("❌ Enter a valid number")

        # ── Select ProLIF conditions ──────────────────────────────────────────
        _, condition_selection = self._prompt_prolif_conditions()
        if condition_selection is None:
            return

        # ── Load fingerprint rows ─────────────────────────────────────────────
        if not os.path.exists(self.__training_sets_db):
            print("\n❌ Training sets database not found.")
            return

        conn = sqlite3.connect(self.__training_sets_db)
        rows = pd.read_sql_query(
            """SELECT assay_name, pose_file, label, fingerprint_json
               FROM training_set_fingerprints
               WHERE training_set_id = ? AND prolif_conditions_id = ?
               ORDER BY binder_type DESC, assay_name, pose_file""",
            conn, params=(training_set_id, condition_selection)
        )
        conn.close()

        if rows.empty:
            print(f"\n❌ No fingerprints found for training set '{training_set_id}' "
                  f"with conditions ID {condition_selection}. "
                  f"Run compute_training_set_fingerprints() first.")
            return

        # ── Build the output DataFrame ────────────────────────────────────────
        ligpose_col = (rows['pose_file'].str.replace('.pdb', '', regex=False) + '_' + rows['assay_name']).reset_index(drop=True)
        label_col   = rows['label'].reset_index(drop=True)

        fp_dicts = [json.loads(fp) for fp in rows['fingerprint_json']]
        fp_df = pd.DataFrame(fp_dicts).fillna(0).astype(int)
        fp_df = fp_df.reset_index(drop=True)

        out_df = pd.concat(
            [ligpose_col.rename('ligpose'), fp_df, label_col.rename('label')],
            axis=1
        )

        print(f"\n📊 {len(out_df)} poses  |  {fp_df.shape[1]} fingerprint features")
        print(f"   Label distribution: {label_col.value_counts().to_dict()}")

        # ── Prompt for output path ────────────────────────────────────────────
        safe_id = training_set_id.replace(' ', '_')
        default_csv = os.path.join(
            self.path, 'ml', 'training_sets',
            f'{safe_id}_cond{condition_selection}_fingerprints.csv'
        )
        path_input = input(
            f"\nOutput CSV path [default: {default_csv}]: "
        ).strip()
        output_path = path_input if path_input else default_csv

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        out_df.to_csv(output_path, index=False)
        print(f"\n✅ Fingerprints exported to: {output_path}")

    # -------------------------------------------------------------------------
    # Private helpers for fingerprint computation
    # -------------------------------------------------------------------------

    def _list_training_set_snapshots(self) -> list:
        """Return a list of snapshot dicts from training_sets_snapshots.db."""
        if not os.path.exists(self.__training_sets_db):
            return []
        try:
            conn = sqlite3.connect(self.__training_sets_db)
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='training_set_snapshots'")
            if not cursor.fetchone():
                conn.close()
                return []
            cursor.execute(
                "SELECT training_set_id, created_at, n_positives, n_negatives, notes "
                "FROM training_set_snapshots ORDER BY created_at ASC"
            )
            rows = cursor.fetchall()
            conn.close()
            return [
                {'training_set_id': r[0], 'created_at': r[1],
                 'n_positives': r[2], 'n_negatives': r[3], 'notes': r[4]}
                for r in rows
            ]
        except Exception as e:
            print(f"Error listing snapshots: {e}")
            return []

    def _prompt_prolif_conditions(self) -> tuple:
        """
        Non-interactively read ProLIF conditions by prompting the user to select
        from the ProLIF_Conditions table in the docking params DB.

        Returns (prolif_params_dict, conditions_id) or (None, None) on cancel/error.
        """
        if not os.path.exists(self.__docking_params_db):
            print(f"\n❌ Docking params database not found: {self.__docking_params_db}")
            return None, None

        try:
            conn = sqlite3.connect(self.__docking_params_db)
            cursor = conn.cursor()
            cursor.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='ProLIF_Conditions'"
            )
            if not cursor.fetchone():
                print("\n❌ No ProLIF_Conditions table found. "
                      "Create ProLIF conditions via MolDock.create_prolif_conditions() first.")
                conn.close()
                return None, None

            cursor.execute("SELECT id, description FROM ProLIF_Conditions ORDER BY id ASC")
            records = cursor.fetchall()
            if not records:
                print("\n❌ No ProLIF conditions registered in the database.")
                conn.close()
                return None, None

            print(f"\n🔬 AVAILABLE PROLIF CONDITIONS")
            print("=" * 70)
            for i, (rec_id, description) in enumerate(records, 1):
                print(f"   {i}. ID {rec_id}: {description}")
            print("=" * 70)

            while True:
                try:
                    sel = input("\nSelect conditions (number or 'c' to cancel): ").strip()
                    if sel.lower() == 'c':
                        conn.close()
                        return None, None
                    sel_idx = int(sel) - 1
                    if 0 <= sel_idx < len(records):
                        selected_id = records[sel_idx][0]
                        break
                    print(f"❌ Enter a number between 1 and {len(records)}")
                except ValueError:
                    print("❌ Enter a valid number")

            cursor.execute(
                "SELECT conditions FROM ProLIF_Conditions WHERE id = ?", (selected_id,)
            )
            result = cursor.fetchone()
            conn.close()

            if not result:
                print(f"\n❌ Could not retrieve conditions for ID {selected_id}")
                return None, None

            prolif_params_dict = json.loads(result[0])
            print(f"\n✅ Selected ProLIF conditions ID {selected_id}")
            return prolif_params_dict, selected_id

        except Exception as e:
            print(f"\n❌ Error loading ProLIF conditions: {e}")
            return None, None

    def _load_assay_info_by_name(self, assay_names: list) -> dict:
        """
        Load assay info dicts keyed by assay_name from the docking assays registry.

        Returns a dict {assay_name: assay_info_dict}.
        """
        if not os.path.exists(self.__docking_assays_db):
            print(f"\n❌ Docking assays DB not found: {self.__docking_assays_db}")
            return {}
        result = {}
        try:
            conn = sqlite3.connect(self.__docking_assays_db)
            placeholders = ','.join(['?' for _ in assay_names])
            cursor = conn.cursor()
            cursor.execute(
                f"""SELECT assay_id, assay_name, table_name, assay_folder_path, receptor_info
                    FROM docking_assays WHERE assay_name IN ({placeholders})""",
                assay_names
            )
            for row in cursor.fetchall():
                assay_id, assay_name, table_name, assay_folder_path, receptor_info_json = row
                try:
                    receptor_info = json.loads(receptor_info_json) if receptor_info_json else {}
                except Exception:
                    receptor_info = {}
                if receptor_info.get('pdbqt_file'):
                    receptor_info['pdbqt_file'] = self._remap_project_path(receptor_info['pdbqt_file'])
                result[assay_name] = {
                    'assay_id': assay_id,
                    'assay_name': assay_name,
                    'table_name': table_name,
                    'assay_folder_path': self._remap_project_path(assay_folder_path),
                    'receptor_info': receptor_info,
                }
            conn.close()
        except Exception as e:
            print(f"Error loading assay info: {e}")
        return result

    def _load_renumbering_dict(self, template_name: str) -> dict:
        """Load the residue renumbering dict for a receptor template from pdbs.db."""
        if not os.path.exists(self.__pdbs_db):
            return {}
        try:
            conn = sqlite3.connect(self.__pdbs_db)
            cursor = conn.cursor()
            cursor.execute(
                "SELECT renumbering_dict FROM pdb_templates WHERE pdb_template_name = ?",
                (template_name,)
            )
            row = cursor.fetchone()
            conn.close()
            if row and row[0]:
                return json.loads(row[0])
        except Exception as e:
            print(f"Error loading renumbering dict for '{template_name}': {e}")
        return {}

    def _remap_project_path(self, stored_path: str) -> str:
        """
        Remap a stored absolute path to the current project location.

        When a project is imported to a new location, paths stored in the
        SQLite databases still reference the original directory.  This method
        detects the first occurrence of a known top-level project subdirectory
        inside *stored_path* and replaces everything before it with self.path,
        leaving the relative suffix intact.
        """
        if not stored_path or os.path.exists(stored_path):
            return stored_path
        for marker in ('/docking/', '/chemspace/', '/ml/', '/dynamics/'):
            pos = stored_path.find(marker)
            if pos != -1:
                relative = stored_path[pos + 1:]  # strip the leading '/'
                return os.path.join(self.path, relative)
        return stored_path

    @contextlib.contextmanager
    def _suppress_output(self):
        """
        Redirect stdout/stderr to devnull at the OS file-descriptor level.

        Several MolDock steps (tleap, sander, ambpdb) shell out via subprocess
        with inherited stdio, so redirecting sys.stdout alone would not silence
        them — the fd itself must be swapped.
        """
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        saved_stdout_fd = os.dup(1)
        saved_stderr_fd = os.dup(2)
        sys.stdout.flush()
        sys.stderr.flush()
        try:
            os.dup2(devnull_fd, 1)
            os.dup2(devnull_fd, 2)
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_stdout_fd, 1)
            os.dup2(saved_stderr_fd, 2)
            os.close(saved_stdout_fd)
            os.close(saved_stderr_fd)
            os.close(devnull_fd)

    def _ensure_fingerprints_table(self):
        """Create the training_set_fingerprints table in the snapshots DB if absent.
        Also migrates existing tables that pre-date the 'minimized' column."""
        conn = sqlite3.connect(self.__training_sets_db)
        conn.execute(
            """CREATE TABLE IF NOT EXISTS training_set_fingerprints (
                id                   INTEGER PRIMARY KEY AUTOINCREMENT,
                training_set_id      TEXT    NOT NULL,
                assay_name           TEXT    NOT NULL,
                pose_file            TEXT    NOT NULL,
                directory            TEXT    NOT NULL,
                binder_type          TEXT    NOT NULL,
                label                INTEGER NOT NULL,
                prolif_conditions_id INTEGER NOT NULL,
                minimized            INTEGER NOT NULL DEFAULT 0,
                fingerprint_json     TEXT    NOT NULL,
                computed_at          TEXT    NOT NULL,
                UNIQUE(training_set_id, assay_name, pose_file, directory, prolif_conditions_id, minimized)
            )"""
        )
        # Migrate existing tables that lack the 'minimized' column
        existing_cols = [r[1] for r in conn.execute(
            "PRAGMA table_info(training_set_fingerprints)"
        ).fetchall()]
        if 'minimized' not in existing_cols:
            conn.execute(
                "ALTER TABLE training_set_fingerprints ADD COLUMN minimized INTEGER NOT NULL DEFAULT 0"
            )
        conn.commit()
        conn.close()

    # -------------------------------------------------------------------------
    # RF predictions on docking assays
    # -------------------------------------------------------------------------

    @staticmethod
    def _chunk_poses_by_ligand(poses_df: pd.DataFrame, chunk_size: int) -> list:
        """
        Group poses_df rows into chunks of roughly chunk_size poses, never
        splitting a single ligand's poses across two chunks (so each ligand's
        tleap prep only ever needs to run once, in whichever chunk it falls in).
        A ligand with more poses than chunk_size simply produces an oversized
        chunk on its own.
        """
        chunks = []
        current_groups = []
        current_count = 0
        for _, group in poses_df.groupby('LigName', sort=False):
            current_groups.append(group)
            current_count += len(group)
            if current_count >= chunk_size:
                chunks.append(pd.concat(current_groups))
                current_groups = []
                current_count = 0
        if current_groups:
            chunks.append(pd.concat(current_groups))
        return chunks

    def _classify_poses_parallel(self, moldock, poses_df, results_db, assay_info,
                                  prolif_params_dict, receptor_file, renumbering_dict,
                                  minimize, max_workers, verbose, chunk_size=None) -> tuple:
        """
        Classify docked poses in parallel using a single process pool kept alive
        across chunks.

        Large assays (hundreds of thousands to millions of poses) are processed
        in chunks (default ~4x max_workers poses per chunk, never splitting a
        ligand's poses across chunks) instead of preparing every ligand and
        restoring every pose up front. Preparing everything up front would leave
        every unique ligand's tleap files (~5 files each) plus every pose's
        restored/withH PDBs sitting in the flat temp_restored_poses/ directory
        for the whole run — millions of entries in one directory for a 1M-pose
        assay, which degrades filesystem performance badly. Instead, each chunk
        prepares only its own ligands, submits its own pose tasks to the shared
        pool, and wipes temp_restored_poses/ before the next chunk starts, so
        the directory only ever holds one chunk's worth of files.

        Ligand-specific tleap files (prepin/frcmod) are prepared once per unique
        ligand within a chunk, sequentially, so concurrent workers never race to
        write the same ligand's files. Each pose then runs its own restore/
        build-complex/minimize/fingerprint pipeline in an isolated scratch
        directory (see _predict_pose_worker) so fixed-name tleap/sander input
        files never collide between workers.

        Returns (fp_records, n_ok, n_failed).
        """
        import shutil
        from concurrent.futures import ProcessPoolExecutor, as_completed

        shared_output_dir = os.path.join(os.path.dirname(results_db), 'temp_restored_poses')

        if chunk_size is None:
            chunk_size = max(max_workers * 4, 1)

        chunks = self._chunk_poses_by_ligand(poses_df, chunk_size)

        fp_records = []
        n_ok = 0
        n_failed = 0
        n_chunks = len(chunks)
        chunks_done = 0

        pbar = tqdm(total=len(poses_df), desc="Classifying docked poses (parallel)", disable=verbose)
        pbar.set_postfix(classified=n_ok, failed=n_failed, chunks=f"{chunks_done}/{n_chunks}")
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                for chunk_idx, chunk_df in enumerate(chunks, 1):
                    os.makedirs(shared_output_dir, exist_ok=True)

                    unique_lignames = chunk_df['LigName'].unique().tolist()
                    ligand_files = {}
                    failed_ligands = []

                    if verbose:
                        print(f"\n🧪 Chunk {chunk_idx}/{len(chunks)}: preparing tleap files "
                              f"for {len(unique_lignames)} ligand(s), {len(chunk_df)} pose(s)...")

                    for ligname in unique_lignames:
                        suppress = contextlib.nullcontext() if verbose else self._suppress_output()
                        try:
                            with suppress:
                                result = moldock._prepare_ligand_tleap_files(
                                    ligname, assay_info, shared_output_dir
                                )
                        except SystemExit:
                            # _prepare_ligand_tleap_files calls sys.exit(1) if the ligand's
                            # SDF cannot be converted to PDB (e.g. missing from ChemSpace) —
                            # treat that as a per-ligand failure instead of killing the run.
                            result = None

                        if not result or result[0] is None or result[1] is None:
                            failed_ligands.append(ligname)
                            continue
                        ligand_files[ligname] = result

                    if failed_ligands:
                        print(f"\n⚠️  Failed to prepare tleap files for {len(failed_ligands)} "
                              f"ligand(s) in chunk {chunk_idx}/{len(chunks)} "
                              f"(antechamber/parmchk2/SDF error — rerun with verbose=True to "
                              f"see details): {', '.join(failed_ligands)}")
                        print("   Poses for these ligands will be skipped.")

                    tasks = []
                    for _, row in chunk_df.iterrows():
                        ligname = row['LigName']
                        if ligname not in ligand_files:
                            n_failed += 1
                            pbar.update(1)
                            pbar.set_postfix(classified=n_ok, failed=n_failed,
                                              chunks=f"{chunks_done}/{n_chunks}")
                            continue
                        prepin_file, frcmod_file = ligand_files[ligname]
                        tasks.append({
                            'project_path': self.path,
                            'project_name': self.name,
                            'results_db': results_db,
                            'ligname': ligname,
                            'run_number': row['run_number'],
                            'pose_id': row['Pose_ID'],
                            'docking_score': row['docking_score'],
                            'assay_info': assay_info,
                            'prolif_params_dict': prolif_params_dict,
                            'receptor_file': receptor_file,
                            'renumbering_dict': renumbering_dict,
                            'minimize': minimize,
                            'prepin_file': prepin_file,
                            'frcmod_file': frcmod_file,
                            'verbose': verbose,
                        })

                    futures = [executor.submit(_predict_pose_worker, task) for task in tasks]
                    for future in as_completed(futures):
                        result = future.result()
                        if result.get('ok'):
                            n_ok += 1
                            fp_records.append({k: v for k, v in result.items() if k != 'ok'})
                        else:
                            n_failed += 1
                            if verbose:
                                print(f"    ❌ Error processing Pose_ID {result['pose_id']} "
                                      f"({result['lig_name']}): {result.get('error')} — skipping")
                        pbar.update(1)
                        pbar.set_postfix(classified=n_ok, failed=n_failed,
                                          chunks=f"{chunks_done}/{n_chunks}")

                    # Wipe this chunk's ligand-prep + restored-pose scratch files before
                    # the next chunk starts, so temp_restored_poses/ never accumulates
                    # entries for the whole (potentially million-pose) run.
                    shutil.rmtree(shared_output_dir, ignore_errors=True)

                    chunks_done += 1
                    pbar.set_postfix(classified=n_ok, failed=n_failed,
                                      chunks=f"{chunks_done}/{n_chunks}")
        finally:
            pbar.close()
            shutil.rmtree(shared_output_dir, ignore_errors=True)
            if verbose:
                print("\n✅ Temporary files cleaned up.")

        return fp_records, n_ok, n_failed

    def _align_and_store_rf_predictions(self, model, fp_records: list,
                                         results_db: str, predictions_table: str) -> None:
        """
        Align computed pose fingerprints to a trained RF model's feature space,
        predict, and store the results in the assay's results DB.

        Shared by the foreground (sequential/parallel) and background
        classification paths of make_rf_predictions_on_docking_assay() so both
        write predictions identically.
        """
        if not fp_records:
            print("\n❌ No fingerprints were successfully computed. Aborting prediction.")
            return

        pred_df = pd.DataFrame(fp_records)
        meta_cols = ['pose_id', 'lig_name', 'docking_score']
        X = pred_df.drop(columns=meta_cols)

        expected_cols = list(model.feature_names_in_)
        X = X.reindex(columns=expected_cols, fill_value=0).fillna(0).astype(int)

        labels = model.predict(X)
        probs = model.predict_proba(X)[:, 1]

        print(f"\n📊 Predictions complete: {int((labels == 1).sum())} positives / "
              f"{int((labels == 0).sum())} negatives out of {len(labels)} poses")

        try:
            conn = sqlite3.connect(results_db)
            conn.execute(f"""
                CREATE TABLE IF NOT EXISTS {predictions_table} (
                    Pose_ID          INTEGER PRIMARY KEY,
                    rf_label         INTEGER NOT NULL,
                    rf_prob_positive REAL NOT NULL
                )
            """)
            conn.executemany(
                f"""INSERT OR REPLACE INTO {predictions_table}
                   (Pose_ID, rf_label, rf_prob_positive)
                   VALUES (?, ?, ?)""",
                [
                    (
                        int(pred_df.iloc[i]['pose_id']),
                        int(labels[i]),
                        float(probs[i]),
                    )
                    for i in range(len(pred_df))
                ]
            )
            conn.commit()
            conn.close()
            print(f"\n✅ Predictions written to '{predictions_table}' table in:\n   {results_db}")
        except Exception as e:
            print(f"\n❌ Error writing predictions to database: {e}")

    def _launch_background_rf_predictions(self, assay_info, model_meta, results_db,
                                           predictions_table, prolif_params_dict,
                                           receptor_file, renumbering_dict, minimize,
                                           max_workers, verbose) -> None:
        """
        Serialize everything needed to classify an assay's poses and store
        predictions, then launch it as a detached background process — a
        separate `python -c` process started with start_new_session=True so it
        survives this interactive session/terminal closing — and return
        immediately. See the module-level _run_background_rf_predictions() for
        the corresponding entry point.

        A plain subprocess is used (rather than multiprocessing.Process)
        because non-daemonic multiprocessing children are implicitly joined at
        interpreter exit, which would block the caller's script from actually
        returning until the classification finished — defeating the purpose.
        """
        import json
        import subprocess
        from datetime import datetime

        jobs_dir = os.path.join(self.path, 'ml', 'results')
        os.makedirs(jobs_dir, exist_ok=True)
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        job_stem = f"rf_prediction_job_{assay_info['assay_id']}_{model_meta['model_id']}_{stamp}"
        job_file = os.path.join(jobs_dir, f"{job_stem}.json")
        log_file = os.path.join(jobs_dir, f"{job_stem}.log")
        rf_db_path = os.path.join(self.__ml_models_path, 'rf_trained_models.db')

        job_spec = {
            'project_name': self.name,
            'project_path': self.path,
            'results_db': results_db,
            'predictions_table': predictions_table,
            'assay_info': assay_info,
            'prolif_params_dict': prolif_params_dict,
            'receptor_file': receptor_file,
            'renumbering_dict': renumbering_dict,
            'minimize': minimize,
            'model_id': model_meta['model_id'],
            'rf_db_path': rf_db_path,
            'max_workers': max_workers,
            'verbose': verbose,
        }
        with open(job_file, 'w') as f:
            json.dump(job_spec, f)

        cmd = [
            sys.executable, '-c',
            "from tidyscreen.ml.ml_functions import _run_background_rf_predictions; "
            f"_run_background_rf_predictions({job_file!r})"
        ]
        with open(log_file, 'w') as log_fh:
            proc = subprocess.Popen(
                cmd, stdout=log_fh, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, start_new_session=True, cwd=self.path
            )

        print(f"\n🚀 Classification launched in the background (PID {proc.pid}).")
        print(f"   Log file: {log_file}")
        print(f"   Predictions will be written to table '{predictions_table}' in:\n   {results_db}")

    def make_rf_predictions_on_docking_assay(self, verbose: bool = False):
        """
        Apply a trained RF model to all poses in a selected docking assay.

        The user is prompted to select a docking assay from the project registry.

        Args:
            verbose (bool): If True, print full per-pose processing details,
                including output from the tleap/sander/ambpdb subprocess calls
                made while preparing each pose. If False (default), all of that
                output is suppressed and a tqdm progress counter is shown instead.
        """

        if not os.path.exists(self.__docking_assays_db):
            print(f"\n❌ Docking assays database not found: {self.__docking_assays_db}")
            return

        try:
            conn = sqlite3.connect(self.__docking_assays_db)
            cursor = conn.cursor()
            cursor.execute(
                """SELECT assay_id, assay_name, table_name, assay_folder_path, receptor_info
                   FROM docking_assays
                   WHERE project_name = ?
                   ORDER BY assay_id ASC""",
                (self.name,)
            )
            rows = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"\n❌ Error reading docking assays: {e}")
            return

        if not rows:
            print("\n❌ No docking assays found for this project.")
            return

        print(f"\n🧬 DOCKING ASSAYS IN PROJECT '{self.name}'")
        print("=" * 70)
        assays = []
        for assay_id, assay_name, table_name, assay_folder_path, receptor_info_json in rows:
            try:
                receptor_info = json.loads(receptor_info_json) if receptor_info_json else {}
            except Exception:
                receptor_info = {}
            assay = {
                'assay_id': assay_id,
                'assay_name': assay_name,
                'table_name': table_name,
                'assay_folder_path': self._remap_project_path(assay_folder_path),
                'receptor_info': receptor_info,
            }
            assays.append(assay)
            print(f"   {len(assays)}. [{assay_id}] {assay_name}  (table: {table_name})")
        print("=" * 70)

        while True:
            try:
                sel = input("\nSelect docking assay (number or 'c' to cancel): ").strip()
                if sel.lower() == 'c':
                    print("❌ Operation cancelled.")
                    return
                idx = int(sel) - 1
                if 0 <= idx < len(assays):
                    selected_assay = assays[idx]
                    break
                print(f"❌ Enter a number between 1 and {len(assays)}")
            except ValueError:
                print("❌ Enter a valid number")

        print(f"\n✅ Selected assay: [{selected_assay['assay_id']}] {selected_assay['assay_name']}")

        # ── 2. Select RF model ────────────────────────────────────────────────
        model, model_meta = self._prompt_rf_model()
        if model is None:
            return

        print(f"\n✅ Selected model: [{model_meta['model_id']}] {model_meta['model_name']}")

        # ── 3. Load all poses from the assay Results table ────────────────────
        results_db = os.path.join(
            selected_assay['assay_folder_path'], 'results',
            f"assay_{selected_assay['assay_id']}.db"
        )

        if not os.path.isfile(results_db):
            print(f"\n❌ Results database not found: {results_db}")
            return

        try:
            conn = sqlite3.connect(results_db)
            tables = [r[0] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()]
            if 'Results' not in tables:
                print(f"\n❌ 'Results' table not found in {results_db}")
                conn.close()
                return
            poses_df = pd.read_sql_query(
                "SELECT Pose_ID, LigName, docking_score, run_number FROM Results",
                conn
            )
            conn.close()
        except Exception as e:
            print(f"\n❌ Error reading Results table: {e}")
            return

        if poses_df.empty:
            print("\n❌ No poses found in the Results table.")
            return

        print(f"\n📋 Loaded {len(poses_df)} poses from assay '{selected_assay['assay_name']}'")

        # ── 3b. Check if the predictions table already exists ─────────────────
        import re
        safe_model_name = re.sub(r'[^a-zA-Z0-9_]', '_', model_meta['model_name'])
        predictions_table = f"rf_predictions_{safe_model_name}"

        try:
            conn = sqlite3.connect(results_db)
            existing_tables = [r[0] for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()]
            conn.close()
        except Exception as e:
            print(f"\n❌ Error checking existing tables: {e}")
            return

        if predictions_table in existing_tables:
            print(f"\n⚠️  Table '{predictions_table}' already exists in the assay database.")
            while True:
                ans = input("Delete existing predictions and recompute? (y/n): ").strip().lower()
                if ans == 'y':
                    try:
                        conn = sqlite3.connect(results_db)
                        conn.execute(f"DROP TABLE {predictions_table}")
                        conn.commit()
                        conn.close()
                        print(f"✅ Table '{predictions_table}' deleted.")
                    except Exception as e:
                        print(f"\n❌ Failed to delete table: {e}")
                        return
                    break
                elif ans == 'n':
                    print("❌ Operation cancelled.")
                    return
                print("❌ Please enter 'y' or 'n'")

        # ── 4. Resolve ProLIF conditions and minimization from training metadata ─
        prolif_params_dict, condition_id, minimize = self._resolve_prolif_conditions_for_model(model_meta)
        if prolif_params_dict is None:
            return

        # ── 5. Resolve receptor and renumbering dict from assay metadata ──────
        receptor_pdbqt = self._remap_project_path(
            selected_assay['receptor_info'].get('pdbqt_file', '')
        )
        receptor_file = os.path.join(os.path.dirname(receptor_pdbqt), 'receptor_checked.pdb')

        template_name = selected_assay['receptor_info'].get('template_name')
        renumbering_dict = self._load_renumbering_dict(template_name) if template_name else {}

        # ── 6. For large assays, offer parallel classification ────────────────
        from moldock.moldock import MolDock

        moldock = MolDock(self.project)

        n_poses = len(poses_df)
        use_parallel = False
        max_workers = None

        if n_poses > 100:
            gpu_note = ""
            if minimize and moldock._check_gpu_available():
                gpu_note = (
                    " Note: minimization uses the GPU, so parallel workers will "
                    "contend for it — consider a small worker count."
                )
            while True:
                ans = input(
                    f"\nThis assay has {n_poses} poses.{gpu_note}\n"
                    f"Use parallel processing to classify poses faster? (y/n): "
                ).strip().lower()
                if ans in ('y', 'n'):
                    use_parallel = ans == 'y'
                    break
                print("❌ Please enter 'y' or 'n'")

            if use_parallel:
                cpu_total = os.cpu_count() or 4
                default_workers = min(cpu_total, n_poses)
                while True:
                    raw = input(f"Number of parallel workers [default {default_workers}]: ").strip()
                    if not raw:
                        max_workers = default_workers
                        break
                    try:
                        max_workers = max(1, int(raw))
                        break
                    except ValueError:
                        print("❌ Enter a valid integer")

                while True:
                    ans = input(
                        "\nRun this classification in the background? Control returns "
                        "immediately; predictions are written to the DB once the "
                        "background job finishes. (y/n): "
                    ).strip().lower()
                    if ans in ('y', 'n'):
                        run_in_background = ans == 'y'
                        break
                    print("❌ Please enter 'y' or 'n'")

                if run_in_background:
                    self._launch_background_rf_predictions(
                        selected_assay, model_meta, results_db, predictions_table,
                        prolif_params_dict, receptor_file, renumbering_dict,
                        minimize, max_workers, verbose
                    )
                    return

        # ── 7. Compute ProLIF fingerprints for each pose ───────────────────────
        if use_parallel:
            fp_records, n_ok, n_failed = self._classify_poses_parallel(
                moldock, poses_df, results_db, selected_assay, prolif_params_dict,
                receptor_file, renumbering_dict, minimize, max_workers, verbose
            )
        else:
            import shutil

            processed_ligands = {}   # {ligname: (prepin_file, frcmod_file)}
            fp_records = []           # [{pose_id, lig_name, docking_score, **fp_dict}, ...]
            output_dir = None
            n_ok = 0
            n_failed = 0

            # When quiet, tqdm writes through a duplicated stdout fd so the progress
            # counter keeps displaying even while fd 1/2 are redirected to devnull
            # below (needed to silence subprocess-inherited tleap/sander/ambpdb output).
            original_stdout = os.fdopen(os.dup(1), 'w') if not verbose else sys.stdout

            try:
                pbar = tqdm(
                    poses_df.iterrows(), total=len(poses_df),
                    desc="Classifying docked poses", file=original_stdout,
                    dynamic_ncols=True, disable=verbose,
                )
                for _, row in pbar:
                    pose_id       = row['Pose_ID']
                    run_number    = row['run_number']
                    ligname       = row['LigName']
                    docking_score = row['docking_score']

                    if verbose:
                        print(f"\n  🔬 Processing Pose_ID: {pose_id}  LigName: {ligname}")

                    suppress = contextlib.nullcontext() if verbose else self._suppress_output()
                    try:
                        with suppress:
                            output_dir, output_file = moldock._restore_single_docked_pose(
                                results_db, ligname, run_number
                            )

                            if ligname not in processed_ligands:
                                prepin_file, frcmod_file = moldock._prepare_ligand_tleap_files(
                                    ligname, selected_assay, output_dir
                                )
                                processed_ligands[ligname] = (prepin_file, frcmod_file)
                            else:
                                prepin_file, frcmod_file = processed_ligands[ligname]

                            prmtop_file, inpcrd_file, output_pdb_file = moldock._prepare_complex_prmtop_inpcrd(
                                prepin_file, frcmod_file, selected_assay, output_dir,
                                output_file, pose_id, ligname=ligname
                            )

                            if minimize:
                                try:
                                    if verbose:
                                        print("    MINIMIZING COMPLEX...")
                                    _, output_pdb_file = moldock._minimize_complex(
                                        prmtop_file, inpcrd_file, output_dir, ligname, pose_id, output_pdb_file
                                    )
                                except Exception as e:
                                    if verbose:
                                        print(f"    ⚠️  Minimization failed: {e}. Using original structure.")

                            fps_df = moldock._compute_prolif_fingerprints3(
                                prolif_params_dict, output_pdb_file, receptor_file,
                                prmtop_file, inpcrd_file, ligname, pose_id
                            )

                            fps_df.columns = fps_df.columns.droplevel(0)
                            if renumbering_dict:
                                fps_df = fps_df.rename(columns=renumbering_dict, level=0)
                            fps_df.columns = [f"{a}_{b}" for a, b in fps_df.columns]

                            fp_row = fps_df.iloc[0].to_dict()
                            fp_records.append({
                                'pose_id':       pose_id,
                                'lig_name':      ligname,
                                'docking_score': docking_score,
                                **{k: bool(v) for k, v in fp_row.items()},
                            })

                        n_ok += 1
                        if verbose:
                            print(f"    ✅ Fingerprint computed for Pose_ID {pose_id}")

                    except Exception as e:
                        n_failed += 1
                        if verbose:
                            print(f"    ❌ Error processing Pose_ID {pose_id} ({ligname}): {e} — skipping")
                        continue
                    finally:
                        if not verbose:
                            pbar.set_postfix(classified=n_ok, failed=n_failed)

            finally:
                pbar.close()
                if not verbose:
                    original_stdout.close()
                if output_dir:
                    shutil.rmtree(output_dir, ignore_errors=True)
                    if verbose:
                        print("\n✅ Temporary files cleaned up.")

        # ── 8. Align features to the model's training space, predict, and store ─
        self._align_and_store_rf_predictions(model, fp_records, results_db, predictions_table)

    def _resolve_prolif_conditions_for_model(self, model_meta: dict) -> tuple:
        """
        Recover the ProLIF conditions ID, params dict, and minimization flag used
        when the training fingerprints for the given RF model were computed.

        If the model has a training_set_id, queries training_set_fingerprints for the
        distinct (prolif_conditions_id, minimized) combinations for that snapshot:
          - exactly one  → auto-selected, no prompt
          - more than one → user picks from the list
        If training_set_id is None (CSV-trained model), falls back to interactive
        _prompt_prolif_conditions() and asks the user about minimization.

        Returns (prolif_params_dict, condition_id, minimize_bool)
        or (None, None, None) on cancel/error.
        """
        training_set_id = model_meta.get('training_set_id')

        if training_set_id:
            if not os.path.exists(self.__training_sets_db):
                print("\n❌ Training sets database not found.")
                return None, None, None

            try:
                conn = sqlite3.connect(self.__training_sets_db)
                rows = conn.execute(
                    "SELECT DISTINCT prolif_conditions_id, minimized "
                    "FROM training_set_fingerprints WHERE training_set_id = ?",
                    (training_set_id,)
                ).fetchall()
                conn.close()
            except Exception as e:
                print(f"\n❌ Error reading training set fingerprints: {e}")
                return None, None, None

            if not rows:
                print(f"\n⚠️  No fingerprints found for training set '{training_set_id}'. "
                      "Falling back to manual selection.")
                prolif_params_dict, condition_id = self._prompt_prolif_conditions()
                if prolif_params_dict is None:
                    return None, None, None
                while True:
                    ans = input("\nMinimize complex before fingerprint computation? (y/n): ").strip().lower()
                    if ans in ('y', 'n'):
                        return prolif_params_dict, condition_id, ans == 'y'
                    print("❌ Please enter 'y' or 'n'")

            if len(rows) == 1:
                selected_id, minimized_int = rows[0]
                minimize = bool(minimized_int)
                print(f"\n✅ Auto-selected ProLIF condition ID {selected_id}, "
                      f"minimized={'yes' if minimize else 'no'} "
                      f"(recovered from training set '{training_set_id}')")
            else:
                print(f"\n⚠️  Multiple featurization settings found for training set '{training_set_id}':")
                for i, (cid, mini) in enumerate(rows, 1):
                    print(f"   {i}. Condition ID {cid}  |  minimized={'yes' if mini else 'no'}")
                while True:
                    try:
                        sel = input("\nSelect setting to match training (number or 'c' to cancel): ").strip()
                        if sel.lower() == 'c':
                            return None, None, None
                        idx = int(sel) - 1
                        if 0 <= idx < len(rows):
                            selected_id, minimized_int = rows[idx]
                            minimize = bool(minimized_int)
                            break
                        print(f"❌ Enter a number between 1 and {len(rows)}")
                    except ValueError:
                        print("❌ Enter a valid number")
        else:
            print("\n⚠️  Model has no linked training set (CSV-trained). "
                  "Please select ProLIF conditions and minimization manually.")
            prolif_params_dict, condition_id = self._prompt_prolif_conditions()
            if prolif_params_dict is None:
                return None, None, None
            while True:
                ans = input("\nMinimize complex before fingerprint computation? (y/n): ").strip().lower()
                if ans in ('y', 'n'):
                    return prolif_params_dict, condition_id, ans == 'y'
                print("❌ Please enter 'y' or 'n'")

        # Load the prolif_params_dict for the resolved condition ID
        if not os.path.exists(self.__docking_params_db):
            print(f"\n❌ Docking params database not found: {self.__docking_params_db}")
            return None, None, None

        try:
            conn = sqlite3.connect(self.__docking_params_db)
            row = conn.execute(
                "SELECT conditions FROM ProLIF_Conditions WHERE id = ?", (selected_id,)
            ).fetchone()
            conn.close()
        except Exception as e:
            print(f"\n❌ Error reading ProLIF conditions: {e}")
            return None, None, None

        if not row:
            print(f"\n❌ ProLIF condition ID {selected_id} not found in params database.")
            return None, None, None

        return json.loads(row[0]), selected_id, minimize

    def _prompt_rf_model(self):
        """
        List trained RF models stored in the project and prompt the user to select one.

        Returns (model, metadata_dict) where model is the deserialized sklearn object
        and metadata_dict contains all non-blob columns.  Returns (None, None) on
        cancel or error.
        """
        import pickle

        rf_db_path = os.path.join(self.__ml_models_path, 'rf_trained_models.db')

        if not os.path.exists(rf_db_path):
            print(f"\n❌ No trained models found. Database not present: {rf_db_path}")
            return None, None

        try:
            conn = sqlite3.connect(rf_db_path)
            cursor = conn.cursor()
            cursor.execute(
                """SELECT model_id, model_name, description, training_set_id,
                          roc_auc, accuracy, macro_f1, cv_roc_mean, cv_roc_std, created_at
                   FROM rf_trained_models
                   ORDER BY model_id ASC"""
            )
            rows = cursor.fetchall()
            conn.close()
        except Exception as e:
            print(f"\n❌ Error reading RF models database: {e}")
            return None, None

        if not rows:
            print("\n❌ No trained RF models found in the database.")
            return None, None

        columns = ('model_id', 'model_name', 'description', 'training_set_id',
                   'roc_auc', 'accuracy', 'macro_f1', 'cv_roc_mean', 'cv_roc_std', 'created_at')

        print(f"\n🤖 TRAINED RF MODELS IN PROJECT '{self.name}'")
        print("=" * 70)
        models_meta = []
        for row in rows:
            meta = dict(zip(columns, row))
            models_meta.append(meta)
            cv_str = (f"{meta['cv_roc_mean']:.4f} ± {meta['cv_roc_std']:.4f}"
                      if meta['cv_roc_mean'] is not None else "—")
            print(f"   {len(models_meta)}. [{meta['model_id']}] {meta['model_name']}")
            print(f"        ROC-AUC: {meta['roc_auc']:.4f}  Accuracy: {meta['accuracy']:.4f}"
                  f"  Macro-F1: {meta['macro_f1']:.4f}  CV ROC-AUC: {cv_str}")
            print(f"        Training set: {meta['training_set_id'] or '—'}"
                  f"   Created: {meta['created_at']}")
            if meta['description']:
                print(f"        Description: {meta['description']}")
        print("=" * 70)

        while True:
            try:
                sel = input("\nSelect model (number or 'c' to cancel): ").strip()
                if sel.lower() == 'c':
                    print("❌ Operation cancelled.")
                    return None, None
                idx = int(sel) - 1
                if 0 <= idx < len(models_meta):
                    selected_meta = models_meta[idx]
                    break
                print(f"❌ Enter a number between 1 and {len(models_meta)}")
            except ValueError:
                print("❌ Enter a valid number")

        try:
            conn = sqlite3.connect(rf_db_path)
            cursor = conn.cursor()
            cursor.execute(
                "SELECT model_pkl FROM rf_trained_models WHERE model_id = ?",
                (selected_meta['model_id'],)
            )
            row = cursor.fetchone()
            conn.close()
            if not row or not row[0]:
                print(f"\n❌ Model blob missing for model_id {selected_meta['model_id']}")
                return None, None
            model = pickle.loads(row[0])
        except Exception as e:
            print(f"\n❌ Failed to load model: {e}")
            return None, None

        return model, selected_meta
