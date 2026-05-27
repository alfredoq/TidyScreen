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
