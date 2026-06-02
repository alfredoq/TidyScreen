from tidyscreen import tidyscreen
import sys
import os
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

    def compute_training_set_fingerprints(self, minimize=True, clean_files=True):
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

        try:
            for _, row in entries_df.iterrows():
                assay_name = row['assay_name']
                pose_file  = row['pose_file']
                directory  = row['directory']
                binder_type = row['binder_type']
                label = 1 if binder_type == 'positive' else 0

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
                    print(f"  ⏭️  Already computed, skipping: {pose_file}")
                    continue

                assay_info = assay_info_by_name.get(assay_name)
                if not assay_info:
                    print(f"  ❌ Assay '{assay_name}' not found in registry — skipping {pose_file}")
                    continue

                # Parse LigName and run_number from pose_file ({LigName}_{run_number}.pdb).
                # The suffix is the per-ligand run_number used during pose extraction,
                # NOT the global Pose_ID from the Ringtail Results table.
                stem = os.path.splitext(pose_file)[0]
                parts = stem.rsplit('_', 1)
                if len(parts) != 2:
                    print(f"  ❌ Cannot parse pose_file '{pose_file}' — skipping")
                    continue
                ligname, pose_id_str = parts
                try:
                    pose_id = int(pose_id_str)  # this is the run_number
                except ValueError:
                    print(f"  ❌ Cannot parse run_number from '{pose_id_str}' — skipping {pose_file}")
                    continue

                # Derive receptor_checked.pdb path from the pdbqt_file stored in assay info
                receptor_pdbqt = assay_info.get('receptor_info', {}).get('pdbqt_file')
                if not receptor_pdbqt:
                    print(f"  ❌ No receptor pdbqt_file for assay '{assay_name}' — skipping {pose_file}")
                    continue
                receptor_file = '/'.join(receptor_pdbqt.split('/')[:-1]) + '/receptor_checked.pdb'

                assay_folder_path = assay_info['assay_folder_path']
                assay_id = assay_info['assay_id']
                results_db = os.path.join(assay_folder_path, 'results', f'assay_{assay_id}.db')

                renumbering_dict = renumbering_dict_by_name.get(assay_name, {})

                print(f"\n  🔬 Processing {pose_file} ({binder_type}) — assay: {assay_name}")

                try:
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
                            print("    MINIMIZING COMPLEX...")
                            min_rst_cpptraj_file, output_pdb_file = moldock._minimize_complex(
                                prmtop_file, inpcrd_file, output_dir, ligname, pose_id, output_pdb_file
                            )
                            rst_cpptraj_file = min_rst_cpptraj_file
                        except Exception as e:
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
                            binder_type, label, prolif_conditions_id, fingerprint_json, computed_at)
                           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                        (training_set_id, assay_name, pose_file, directory,
                         binder_type, label, condition_selection, fp_json,
                         datetime.now().isoformat())
                    )
                    conn.commit()
                    conn.close()

                    print(f"    ✅ Fingerprint stored for {pose_file}")

                except Exception as e:
                    print(f"    ❌ Error processing {pose_file}: {e} — skipping")
                    continue

        finally:
            if clean_files:
                for d in output_dirs:
                    shutil.rmtree(d, ignore_errors=True)
                print("\n✅ Temporary files cleaned up.")

        print(f"\n🎉 Fingerprint computation complete for training set '{training_set_id}'")

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

        # ── Load fingerprint rows ─────────────────────────────────────────────
        if not os.path.exists(self.__training_sets_db):
            print("\n❌ Training sets database not found.")
            return

        conn = sqlite3.connect(self.__training_sets_db)
        rows = pd.read_sql_query(
            """SELECT assay_name, pose_file, label, fingerprint_json
               FROM training_set_fingerprints
               WHERE training_set_id = ?
               ORDER BY binder_type DESC, assay_name, pose_file""",
            conn, params=(training_set_id,)
        )
        conn.close()

        if rows.empty:
            print(f"\n❌ No fingerprints found for training set '{training_set_id}'. "
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
            f'{safe_id}_fingerprints.csv'
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
                result[assay_name] = {
                    'assay_id': assay_id,
                    'assay_name': assay_name,
                    'table_name': table_name,
                    'assay_folder_path': assay_folder_path,
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

    def _ensure_fingerprints_table(self):
        """Create the training_set_fingerprints table in the snapshots DB if absent."""
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
                fingerprint_json     TEXT    NOT NULL,
                computed_at          TEXT    NOT NULL,
                UNIQUE(training_set_id, assay_name, pose_file, directory, prolif_conditions_id)
            )"""
        )
        conn.commit()
        conn.close()
