"""
Random Forest classification pipeline.

Can be imported as a module (by the Streamlit app) or run directly:

    python rf_classifier.py data.csv --label label --drop ligpose --tune --save model.joblib
"""
from __future__ import annotations

import argparse
import logging

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

DEFAULT_PARAM_GRID: dict = {
    "n_estimators": [100, 200, 500],
    "max_depth": [None, 10, 20],
    "min_samples_split": [2, 5, 10],
    "min_samples_leaf": [1, 2, 4],
    "max_features": ["sqrt", "log2"],
}


# ── Data preparation ───────────────────────────────────────────────────────

def prepare_features(
    df: pd.DataFrame,
    label_col: str = "label",
    drop_cols: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (X, y) from a DataFrame, dropping the label and any extra columns."""
    to_drop = [label_col] + (drop_cols or [])
    to_drop = [c for c in to_drop if c in df.columns]
    X = df.drop(columns=to_drop)
    y = df[label_col]
    return X, y


def split_data(
    X: pd.DataFrame,
    y: pd.Series,
    test_size: float = 0.2,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
    return train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)


# ── Model training ─────────────────────────────────────────────────────────

def train_baseline(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_estimators: int = 200,
    max_depth: int | None = None,
    min_samples_split: int = 2,
    min_samples_leaf: int = 1,
    max_features: str = "sqrt",
    cv: int = 5,
    random_state: int = 42,
) -> tuple[RandomForestClassifier, np.ndarray]:
    """Fit a Random Forest with fixed hyperparameters and evaluate by cross-validation."""
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        min_samples_leaf=min_samples_leaf,
        max_features=max_features,
        random_state=random_state,
        n_jobs=-1,
    )
    rf.fit(X_train, y_train)
    cv_scores = cross_val_score(rf, X_train, y_train, cv=cv, scoring="roc_auc", n_jobs=-1)
    return rf, cv_scores


def tune_hyperparameters(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    param_grid: dict | None = None,
    cv: int = 5,
    scoring: str = "roc_auc",
    random_state: int = 42,
) -> tuple[RandomForestClassifier, dict, float]:
    """Run GridSearchCV and return (best_model, best_params, best_cv_score)."""
    grid = param_grid if param_grid is not None else DEFAULT_PARAM_GRID
    search = GridSearchCV(
        RandomForestClassifier(random_state=random_state),
        grid,
        cv=cv,
        scoring=scoring,
        n_jobs=-1,
        verbose=1,
    )
    search.fit(X_train, y_train)
    return search.best_estimator_, search.best_params_, search.best_score_


# ── Evaluation ─────────────────────────────────────────────────────────────

def evaluate_model(
    model: RandomForestClassifier,
    X_test: pd.DataFrame,
    y_test: pd.Series,
) -> dict:
    """Return a dict with confusion matrix, classification report, ROC-AUC, and predictions."""
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    return {
        "confusion_matrix": confusion_matrix(y_test, y_pred),
        "classification_report": classification_report(y_test, y_pred, output_dict=True),
        "roc_auc": roc_auc_score(y_test, y_prob),
        "y_pred": y_pred,
        "y_prob": y_prob,
    }


def feature_importances(
    model: RandomForestClassifier,
    feature_names: list[str],
    top_n: int = 20,
) -> pd.Series:
    return (
        pd.Series(model.feature_importances_, index=feature_names)
        .sort_values(ascending=False)
        .head(top_n)
    )


# ── Persistence ────────────────────────────────────────────────────────────

def save_model(model: RandomForestClassifier, path: str) -> None:
    joblib.dump(model, path)
    log.info("Model saved → %s", path)


def load_model(path: str) -> RandomForestClassifier:
    return joblib.load(path)


# ── Full pipeline (CLI / batch use) ───────────────────────────────────────

def run_pipeline(
    csv_path: str,
    label_col: str = "label",
    drop_cols: list[str] | None = None,
    test_size: float = 0.2,
    cv: int = 5,
    random_state: int = 42,
    tune: bool = False,
    param_grid: dict | None = None,
    top_n: int = 20,
    save_path: str | None = None,
) -> dict:
    df = pd.read_csv(csv_path)
    X, y = prepare_features(df, label_col=label_col, drop_cols=drop_cols)
    X_train, X_test, y_train, y_test = split_data(X, y, test_size=test_size, random_state=random_state)

    if tune:
        log.info("Running GridSearchCV…")
        model, best_params, best_cv_score = tune_hyperparameters(
            X_train, y_train, param_grid=param_grid, cv=cv, random_state=random_state
        )
        cv_scores = None
    else:
        log.info("Training baseline model…")
        model, cv_scores = train_baseline(
            X_train, y_train, cv=cv, random_state=random_state
        )
        best_params, best_cv_score = None, None

    eval_res = evaluate_model(model, X_test, y_test)
    fi = feature_importances(model, X.columns.tolist(), top_n=top_n)

    if save_path:
        save_model(model, save_path)

    return {
        "model": model,
        "cv_scores": cv_scores,
        "best_params": best_params,
        "best_cv_score": best_cv_score,
        "evaluation": eval_res,
        "feature_importances": fi,
    }


# ── CLI entry point ────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Random Forest classifier pipeline")
    parser.add_argument("csv", help="Path to input CSV file")
    parser.add_argument("--label", default="label", help="Label column name (default: label)")
    parser.add_argument("--drop", nargs="*", default=["ligpose"], metavar="COL",
                        help="Extra columns to drop before training")
    parser.add_argument("--test-size", type=float, default=0.2, metavar="FRAC")
    parser.add_argument("--cv", type=int, default=5, metavar="N", help="CV folds")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tune", action="store_true", help="Run GridSearchCV")
    parser.add_argument("--top-n", type=int, default=20, help="Top N features to print")
    parser.add_argument("--save", metavar="PATH", help="Save trained model to PATH")
    args = parser.parse_args()

    out = run_pipeline(
        csv_path=args.csv,
        label_col=args.label,
        drop_cols=args.drop,
        test_size=args.test_size,
        cv=args.cv,
        random_state=args.seed,
        tune=args.tune,
        top_n=args.top_n,
        save_path=args.save,
    )

    if out["cv_scores"] is not None:
        cv = out["cv_scores"]
        print(f"\nCross-validation ROC-AUC: {cv}")
        print(f"Mean ± SD: {cv.mean():.4f} ± {cv.std():.4f}")

    if out["best_params"] is not None:
        print(f"\nBest hyperparameters: {out['best_params']}")
        print(f"Best CV ROC-AUC: {out['best_cv_score']:.4f}")

    ev = out["evaluation"]
    print(f"\nTest ROC-AUC: {ev['roc_auc']:.4f}")
    print("\nConfusion Matrix:")
    print(ev["confusion_matrix"])
    print("\nClassification Report:")
    print(pd.DataFrame(ev["classification_report"]).T.to_string())
    print(f"\nTop {args.top_n} Feature Importances:")
    print(out["feature_importances"].to_string())
