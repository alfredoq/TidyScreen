import streamlit as st
from tidyscreen import tidyscreen
from tidyscreen.chemspace.chemspace import ChemSpace, _get_descriptor_funcs
from tidyscreen.moldock.moldock import MolDock
from tidyscreen.moldyn.moldyn import MolDyn
from tidyscreen.ml.ml_functions import MachineLearning
import io
import contextlib
import inspect
import sys
import os
import glob
import shutil
import json
import sqlite3
import py3Dmol
import streamlit_functions as st_funcs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rf_classifier import (
    evaluate_model,
    feature_importances,
    prepare_features,
    split_data,
    train_baseline,
    tune_hyperparameters,
)

tidyscreen_package_path = sys.argv[1]


def _toggle_ml_positive_viewer():
    key = "show_ml_positive_poses"
    st.session_state[key] = not st.session_state.get(key, False)


def _toggle_ml_negative_viewer():
    key = "show_ml_negative_poses"
    st.session_state[key] = not st.session_state.get(key, False)


def _toggle_show_projects():
    key = "show_projects"
    st.session_state[key] = not st.session_state.get(key, False)


def _toggle_ts_set_viewer():
    key = "show_ts_set_viewer"
    st.session_state[key] = not st.session_state.get(key, False)


def _toggle_ts_snapshot_inspector():
    key = "show_ts_snapshot_inspector"
    st.session_state[key] = not st.session_state.get(key, False)


def _clear_docking_assay_scoped_state():
    """
    Drop session_state entries scoped to the previously selected docking assay.
    Docking analysis widgets (pose browser, ProLIF filters, ligand/pose selectors,
    Classified Poses Viewer) key their selection lists off the active assay's data
    (e.g. LigName, pose_id). Left in place after switching assays, a stale value
    looked up via list.index() against the new assay's option list raises
    "'<value>' is not in list" — call this whenever selected_assay_name changes.
    """
    _exact_keys = {
        "selected_lig_name", "select_lig_name",
        "selected_pose_id", "select_pose_id",
        "reference_pdb_data",
        "confirm_delete_all_poses",
        "extract_poses_criteria", "extract_poses_score_column",
        "select_prolif_table", "select_prolif_pose",
        "show_all_prolif",
    }
    _prefixes = (
        "btn_poses_", "pose_idx_",
        "show_prolif_", "select_prolif_table_", "select_prolif_pose_",
        "prolif_filter_", "show_all_prolif_",
        "cpv_",
    )
    for _key in list(st.session_state.keys()):
        if _key in _exact_keys or _key.startswith(_prefixes):
            del st.session_state[_key]


# ── GridSearch hyperparameter input callbacks ──────────────────────────────
# Each callback parses the text input, appends new values to the accumulated
# list, then clears the input box.  Defined at module level so they are
# available before the page block executes.

def _gs_cb_n_est():
    raw = st.session_state.get("rf_gs_n_est_extra", "").strip()
    if raw:
        new_vals = [int(t.strip()) for t in raw.split(",") if t.strip().lstrip("-").isdigit()]
        if new_vals:
            acc = st.session_state.get("rf_gs_n_est_accum", [])
            st.session_state["rf_gs_n_est_accum"] = list(dict.fromkeys(acc + new_vals))
    st.session_state["rf_gs_n_est_extra"] = ""


def _gs_cb_depth():
    raw = st.session_state.get("rf_gs_depth_extra", "").strip()
    if raw:
        new_vals = []
        for token in raw.split(","):
            token = token.strip()
            if not token:
                continue
            if token.lower() == "none":
                new_vals.append("None")
            else:
                try:
                    new_vals.append(str(int(token)))
                except ValueError:
                    pass
        if new_vals:
            acc = st.session_state.get("rf_gs_depth_accum", [])
            st.session_state["rf_gs_depth_accum"] = list(dict.fromkeys(acc + new_vals))
    st.session_state["rf_gs_depth_extra"] = ""


def _gs_cb_split():
    raw = st.session_state.get("rf_gs_split_extra", "").strip()
    if raw:
        new_vals = [int(t.strip()) for t in raw.split(",") if t.strip().lstrip("-").isdigit()]
        if new_vals:
            acc = st.session_state.get("rf_gs_split_accum", [])
            st.session_state["rf_gs_split_accum"] = list(dict.fromkeys(acc + new_vals))
    st.session_state["rf_gs_split_extra"] = ""


def _gs_cb_leaf():
    raw = st.session_state.get("rf_gs_leaf_extra", "").strip()
    if raw:
        new_vals = [int(t.strip()) for t in raw.split(",") if t.strip().lstrip("-").isdigit()]
        if new_vals:
            acc = st.session_state.get("rf_gs_leaf_accum", [])
            st.session_state["rf_gs_leaf_accum"] = list(dict.fromkeys(acc + new_vals))
    st.session_state["rf_gs_leaf_extra"] = ""


def _gs_cb_feat():
    raw = st.session_state.get("rf_gs_feat_extra", "").strip()
    if raw:
        new_vals = [t.strip() for t in raw.split(",") if t.strip()]
        if new_vals:
            acc = st.session_state.get("rf_gs_feat_accum", [])
            st.session_state["rf_gs_feat_accum"] = list(dict.fromkeys(acc + new_vals))
    st.session_state["rf_gs_feat_extra"] = ""


_WORKFLOW_VAR_BY_CLASS_LABEL = {
    "TidyScreen": "tidyscreen",
    "ChemSpace": "active_project_cs",
    "MolDock": "active_project_moldock",
    "MolDyn": "active_project_moldyn",
    "MachineLearning": "active_project_ml",
}


def _build_workflow_script(project_name, selections):
    _lines = [
        "",
        "from tidyscreen import tidyscreen",
        "from tidyscreen.chemspace.chemspace import ChemSpace",
        "from tidyscreen.moldock.moldock import MolDock",
        "from tidyscreen.moldyn.moldyn import MolDyn",
        "from tidyscreen.ml.ml_functions import MachineLearning",
        "",
        f'active_project = tidyscreen.ActivateProject("{project_name}")',
        "active_project_cs = ChemSpace(active_project)",
        "active_project_moldock = MolDock(active_project)",
        "active_project_ml = MachineLearning(active_project)",
        "active_project_moldyn = MolDyn(active_project)",
        "",
    ]

    for _class_label, _entries in selections.items():
        if not _entries:
            continue
        _var_name = _WORKFLOW_VAR_BY_CLASS_LABEL[_class_label]
        for _method_name, _method_signature, _method_summary in _entries:
            if _method_summary:
                _lines.append(f"## {_method_summary}")
            _lines.append(f"#{_var_name}.{_method_name}()")
            _lines.append("")

    return "\n".join(_lines)


def _render_python_api_section(header, obj, class_label, search_lower, include_classes=False, exclude_names=None):
    _exclude_names = exclude_names or set()
    with st.expander(header, expanded=True):
        st.write(
            f"Public methods available in `{class_label}` (helper methods, i.e. names "
            "starting with `_`, are excluded). Check the methods to include in a "
            "generated workflow file (see 'Build Workflow File' below)."
        )
        if include_classes:
            _module_name = getattr(obj, "__name__", None)
            _methods = [
                (name, member)
                for name, member in inspect.getmembers(
                    obj, predicate=lambda m: inspect.isfunction(m) or inspect.isclass(m)
                )
                if not name.startswith("_")
                and name not in _exclude_names
                and getattr(member, "__module__", None) == _module_name
            ]
        else:
            _methods = [
                (name, member)
                for name, member in inspect.getmembers(obj, predicate=inspect.isfunction)
                if not name.startswith("_") and name not in _exclude_names
            ]

        _all_entries = []
        for _method_name, _method in sorted(_methods, key=lambda item: item[0]):
            try:
                _method_signature = str(inspect.signature(_method))
            except (TypeError, ValueError):
                _method_signature = "(...)"
            _method_doc = inspect.getdoc(_method)
            _method_summary = _method_doc.strip().split("\n")[0] if _method_doc else ""
            _all_entries.append((_method_name, _method_signature, _method_summary))

        _display_entries = _all_entries
        if search_lower:
            _display_entries = [
                entry for entry in _all_entries
                if search_lower in entry[0].lower() or search_lower in entry[2].lower()
            ]

        st.caption(f"{len(_display_entries)} method(s) found.")

        _sel_key = f"python_api_selected_{class_label}"
        _prior_selected = set(st.session_state.get(_sel_key, []))

        if not _display_entries:
            st.info("No methods match the current filter.")
        else:
            _df_api = pd.DataFrame(
                {
                    "Select": [_name in _prior_selected for _name, _, _ in _display_entries],
                    "Method": [f"{_name}{_sig}" for _name, _sig, _ in _display_entries],
                    "Description": [_summary for _, _, _summary in _display_entries],
                }
            )
            _edited_api = st.data_editor(
                _df_api,
                column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                disabled=["Method", "Description"],
                hide_index=True,
                use_container_width=True,
                key=f"python_api_data_editor_{class_label}_{search_lower}",
            )
            _visible_names = {_name for _name, _, _ in _display_entries}
            _newly_selected = {
                _display_entries[_i][0]
                for _i in range(len(_display_entries))
                if _edited_api.iloc[_i]["Select"]
            }
            _prior_selected = (_prior_selected - _visible_names) | _newly_selected
            st.session_state[_sel_key] = sorted(_prior_selected)

        _selected_now = set(st.session_state.get(_sel_key, []))
        return [entry for entry in _all_entries if entry[0] in _selected_now]


st.set_page_config(page_title="TidyScreen App", layout="wide")

# Sidebar navigation
_SIDEBAR_SEPARATOR = "──────────"
_SIDEBAR_SECTION_CHEMSPACE = "ChemSpace"
_SIDEBAR_SECTION_HELPERS = "Helpers"
_SIDEBAR_SECTION_ML = "Machine learning"
_SIDEBAR_SECTION_MOLDYN = "MolDyn"
_SIDEBAR_SECTION_MOLDOCK = "MolDock"
_SIDEBAR_NON_SELECTABLE_OPTIONS = {_SIDEBAR_SEPARATOR, _SIDEBAR_SECTION_CHEMSPACE, _SIDEBAR_SECTION_HELPERS, _SIDEBAR_SECTION_ML, _SIDEBAR_SECTION_MOLDYN, _SIDEBAR_SECTION_MOLDOCK}
_SIDEBAR_SECTION_HEADERS = {_SIDEBAR_SECTION_CHEMSPACE, _SIDEBAR_SECTION_HELPERS, _SIDEBAR_SECTION_ML, _SIDEBAR_SECTION_MOLDYN, _SIDEBAR_SECTION_MOLDOCK}


def _on_sidebar_page_change():
    if st.session_state.get("_sidebar_page_radio") in _SIDEBAR_NON_SELECTABLE_OPTIONS:
        st.session_state["_sidebar_page_radio"] = st.session_state.get("_sidebar_last_page", "Project selection")
    else:
        st.session_state["_sidebar_last_page"] = st.session_state["_sidebar_page_radio"]


_logo_path = os.path.join(tidyscreen_package_path, "images", "tidyscreen_logo.png")
if os.path.exists(_logo_path):
    st.sidebar.image(_logo_path, use_container_width=True)
_sidebar_page_options = ("Project selection", _SIDEBAR_SEPARATOR, _SIDEBAR_SECTION_CHEMSPACE, "ChemSpace Inspection", "ChemSpace Actions", _SIDEBAR_SEPARATOR, _SIDEBAR_SECTION_MOLDOCK, "Receptors", "MolDock assays", "ProLIF Conditions", "Docking analysis", _SIDEBAR_SEPARATOR, _SIDEBAR_SECTION_MOLDYN, "MolDyn assays", "MolDyn analysis", _SIDEBAR_SEPARATOR, _SIDEBAR_SECTION_ML, "ML features management", "RF model training", _SIDEBAR_SEPARATOR, _SIDEBAR_SECTION_HELPERS, "Mol Viewer", "Reaction Viewer", "Python API")
_separator_rules = "\n".join(
    f"""
    section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type({i + 1}) {{
        pointer-events: none;
    }}
    section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type({i + 1}) > div:first-child {{
        display: none;
    }}
    """
    + (
        f"""
    section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type({i + 1}) p {{
        font-size: 0.875rem;
    }}
    section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type({i + 1}) {{
        padding-left: 0;
        margin-left: 0;
    }}
    section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type({i + 1}) > div:last-child {{
        padding-left: 0;
        margin-left: 0;
    }}
    """
        if opt in _SIDEBAR_SECTION_HEADERS
        else ""
    )
    for i, opt in enumerate(_sidebar_page_options)
    if opt in _SIDEBAR_NON_SELECTABLE_OPTIONS
)
st.sidebar.markdown(f"<style>{_separator_rules}</style>", unsafe_allow_html=True)
page = st.sidebar.radio(
    "Project management",
    _sidebar_page_options,
    key="_sidebar_page_radio",
    on_change=_on_sidebar_page_change,
)

## Persistent sidebar info: active project and assay
st.sidebar.divider()
active_project = st.session_state.get("activate_project_selectbox") or st.session_state.get("selected_project")
active_assay = st.session_state.get("select_assay_name") or st.session_state.get("selected_assay_name")
st.sidebar.markdown("**Active project:**")
st.sidebar.info(active_project if active_project else "None selected")
st.sidebar.markdown("**Active docking assay:**")
st.sidebar.info(active_assay if active_assay else "None selected")

if page == "Project selection":
    st.title("Project selection")
    st.write("Welcome to the TidyScreen main page.")

    db_path = os.path.join(tidyscreen_package_path, "projects_db", "projects_database.db")
    df = st_funcs.read_database_as_dataframe(db_path, "projects")
    
    ## Create a button to show the projects DataFrame
    st.button(
        f"{'Hide' if st.session_state.get('show_projects', False) else 'Show'} Projects",
        on_click=_toggle_show_projects
    )
    if st.session_state.get("show_projects", False):
        st.dataframe(df)
    
    ## Create a selectbox to choose and activate a project
    if not df.empty:
        project_names = df["name"].tolist()
        # Use session_state to persist selection
        if "selected_project" not in st.session_state:
            st.session_state["selected_project"] = project_names[0] if project_names else ""
        selected_project = st.selectbox(
            "Select a project to activate:",
            project_names,
            key="activate_project_selectbox",
            index=project_names.index(st.session_state["selected_project"]) if st.session_state["selected_project"] in project_names else 0
        )
        # Update session_state if selection changes
        if selected_project != st.session_state["selected_project"]:
            st.session_state["selected_project"] = selected_project
        if selected_project:
            project_path = df.loc[df["name"] == selected_project, "path"].values[0]
            st.session_state["active_project_path"] = project_path
            st.success(f"Activated project: {selected_project}")
    else:
        st.warning("No projects found in the database.")

    ## --- Notes Management ---
    st.divider()
    _notes_project = st.session_state.get("selected_project")
    st.subheader(f"Notes — {_notes_project}" if _notes_project else "Notes")

    if "active_project_path" not in st.session_state:
        st.info("Activate a project above to manage notes.")
    else:
        _notes_db = os.path.join(st.session_state["active_project_path"], "project_notes.db")

        _nc = sqlite3.connect(_notes_db)
        try:
            _nc.execute("""
                CREATE TABLE IF NOT EXISTS notes (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    title TEXT NOT NULL,
                    content TEXT NOT NULL,
                    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
                )
            """)
            _nc.commit()
        finally:
            _nc.close()

        def _fetch_notes():
            _c = sqlite3.connect(_notes_db)
            try:
                return _c.execute(
                    "SELECT id, title, created_at FROM notes ORDER BY created_at ASC"
                ).fetchall()
            finally:
                _c.close()

        ## Add note form — processed before fetching list so new note appears immediately
        with st.expander("Add a New Note", expanded=False):
            with st.form("notes_add_form"):
                _add_title = st.text_input("Title")
                _add_content = st.text_area("Content")
                _add_submitted = st.form_submit_button("Add Note")
            if _add_submitted:
                if not _add_title.strip() or not _add_content.strip():
                    st.error("Title and content cannot be empty.")
                else:
                    _c = sqlite3.connect(_notes_db)
                    try:
                        _c.execute(
                            "INSERT INTO notes (title, content) VALUES (?, ?)",
                            (_add_title.strip(), _add_content.strip())
                        )
                        _c.commit()
                    finally:
                        _c.close()
                    st.rerun()

        ## Notes list
        _notes = _fetch_notes()

        if not _notes:
            st.info("No notes yet for this project.")
        else:
            with st.expander(f"Notes list ({len(_notes)})", expanded=False):
                st.dataframe(
                    [{"ID": n[0], "Title": n[1], "Created At": n[2]} for n in _notes],
                    use_container_width=True,
                    hide_index=True,
                )

            st.divider()

            ## View / Modify / Delete
            _note_labels = [f"[{n[0]}] {n[1]}" for n in _notes]
            _sel_label = st.selectbox("Select a note:", _note_labels, key="notes_select")
            _sel_id = int(_sel_label.split("]")[0].strip("["))

            _c = sqlite3.connect(_notes_db)
            try:
                _note = _c.execute(
                    "SELECT id, title, content, created_at FROM notes WHERE id = ?",
                    (_sel_id,)
                ).fetchone()
            finally:
                _c.close()

            if _note:
                st.markdown(f"**Created At:** {_note[3]}")

                with st.expander("View Content", expanded=False):
                    st.markdown(_note[2])

                with st.expander("Modify Note", expanded=False):
                    with st.form(f"notes_modify_form_{_sel_id}"):
                        _mod_title = st.text_input("Title", value=_note[1])
                        _mod_content = st.text_area("Content", value=_note[2])
                        _mod_submitted = st.form_submit_button("Save Changes")
                    if _mod_submitted:
                        if not _mod_title.strip() or not _mod_content.strip():
                            st.error("Title and content cannot be empty.")
                        else:
                            _c = sqlite3.connect(_notes_db)
                            try:
                                _c.execute(
                                    "UPDATE notes SET title = ?, content = ? WHERE id = ?",
                                    (_mod_title.strip(), _mod_content.strip(), _sel_id)
                                )
                                _c.commit()
                            finally:
                                _c.close()
                            st.rerun()

                _del_key = f"confirm_delete_note_{_sel_id}"
                if not st.session_state.get(_del_key):
                    if st.button("Delete Note", key=f"btn_delete_note_{_sel_id}"):
                        st.session_state[_del_key] = True
                        st.rerun()
                else:
                    st.warning(f"Delete note '{_note[1]}'?")
                    _dc1, _dc2 = st.columns(2)
                    with _dc1:
                        if st.button("Yes, delete", key=f"btn_confirm_del_note_{_sel_id}"):
                            _c = sqlite3.connect(_notes_db)
                            try:
                                _c.execute("DELETE FROM notes WHERE id = ?", (_sel_id,))
                                _c.commit()
                            finally:
                                _c.close()
                            st.session_state[_del_key] = False
                            st.rerun()
                    with _dc2:
                        if st.button("Cancel", key=f"btn_cancel_del_note_{_sel_id}"):
                            st.session_state[_del_key] = False
                            st.rerun()

elif page == "ChemSpace Inspection":
    st.title("ChemSpace Inspection")
    st.write("Welcome to the ChemSpace page.")
    
    db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
    _chemspace_table_sigs = st_funcs.get_table_signatures(db_path)
    df = st_funcs.get_tables_info(db_path, signatures=_chemspace_table_sigs)
    
    ## Create a button to show the project tables info DataFrame
    if "show_tables_info" not in st.session_state:
        st.session_state["show_tables_info"] = False
    if st.button(f"{'Hide' if st.session_state['show_tables_info'] else 'Show'} ChemSpace Tables Info", key="btn_chemspace_tables_info"):
        st.session_state["show_tables_info"] = not st.session_state["show_tables_info"]
        st.rerun()
    if st.session_state["show_tables_info"]:
        st.dataframe(df)

        if df is not None and not df.empty:
            _drop_table_names = df["table"].tolist()
            # Sanitize any previously-selected names that no longer exist (e.g.
            # after a drop) *before* the widget below reads/instantiates its
            # state — Streamlit forbids changing a widget's session_state value
            # after it has been instantiated in the same run.
            if "drop_tables_select" in st.session_state:
                st.session_state["drop_tables_select"] = [
                    t for t in st.session_state["drop_tables_select"] if t in _drop_table_names
                ]
            _drop_selected = st.multiselect(
                "Select table(s) to drop:",
                _drop_table_names,
                key="drop_tables_select"
            )
            _drop_key = "confirm_drop_tables"
            if not _drop_selected:
                st.button(
                    "🗑️ Drop Tables",
                    key="btn_drop_tables",
                    disabled=True,
                    help="Select at least one table above to enable deletion.",
                )
            elif not st.session_state.get(_drop_key):
                if st.button("🗑️ Drop Tables", key="btn_drop_tables"):
                    st.session_state[_drop_key] = True
                    st.rerun()
            else:
                st.warning(
                    f"⚠️ Delete {len(_drop_selected)} table(s): **{', '.join(_drop_selected)}**? "
                    "This action cannot be undone!"
                )
                _dtc1, _dtc2 = st.columns(2)
                with _dtc1:
                    if st.button("Yes, delete", key="btn_confirm_drop_tables"):
                        _drop_results = st_funcs.drop_chemspace_tables(db_path, _drop_selected)
                        _dropped = [t for t, r in _drop_results.items() if r == "dropped"]
                        _failed = {t: r for t, r in _drop_results.items() if r != "dropped"}
                        if _dropped:
                            st.success(f"✅ Dropped table(s): {', '.join(_dropped)}")
                        for _t, _r in _failed.items():
                            st.error(f"❌ Failed to drop '{_t}': {_r}")
                        st.session_state[_drop_key] = False
                        st.rerun()
                with _dtc2:
                    if st.button("Cancel", key="btn_cancel_drop_tables"):
                        st.session_state[_drop_key] = False
                        st.rerun()

    ## Display full content of a selected table
    st.divider()
    st.subheader("Display Table")

    if "show_display_table" not in st.session_state:
        st.session_state["show_display_table"] = False

    if df is not None and not df.empty:
        if st.button(f"{'Hide' if st.session_state['show_display_table'] else 'Show'} Display Table", key="btn_display_table"):
            st.session_state["show_display_table"] = not st.session_state["show_display_table"]
            st.rerun()

        if st.session_state["show_display_table"]:
            display_table_names = df["table"].tolist()
            display_selected_table = st.selectbox(
                "Select a table to display:",
                display_table_names,
                key="display_table_select"
            )
            display_columns = st_funcs.get_table_columns(db_path, display_selected_table, mtime=_chemspace_table_sigs.get(display_selected_table))
            if display_columns:
                st.markdown("**Select columns to display:**")
                selected_display_cols = [col for col in display_columns if st.checkbox(col, value=True, key=f"display_col_{display_selected_table}_{col}")]
                if selected_display_cols:
                    _display_row_limit = 10000
                    _table_total_rows = int(df.loc[df["table"] == display_selected_table, "rows"].iloc[0])
                    _read_limit = _display_row_limit if _table_total_rows > _display_row_limit else None
                    display_df = st_funcs.read_table_columns_as_dataframe(
                        db_path, display_selected_table, selected_display_cols,
                        limit=_read_limit, mtime=_chemspace_table_sigs.get(display_selected_table)
                    )
                    if display_df is not None and not display_df.empty:
                        if _read_limit is not None:
                            st.info(
                                f"Table '{display_selected_table}' has {_table_total_rows:,} rows; "
                                f"showing only the first {_display_row_limit:,}."
                            )
                        _display_df_sel = display_df.copy()
                        _display_df_sel.insert(0, "Select", False)
                        _edited_display = st.data_editor(
                            _display_df_sel,
                            column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                            disabled=[c for c in _display_df_sel.columns if c != "Select"],
                            hide_index=True,
                            use_container_width=True,
                            key=f"data_editor_display_{display_selected_table}",
                        )
                        _selected_display_rows = _edited_display[_edited_display["Select"]]
                        _n_sel_display = len(_selected_display_rows)

                        st.markdown("---")
                        st.markdown("**📦 Subset Table**")
                        _subset_col1, _subset_col2, _subset_col3 = st.columns([3, 1, 1])
                        with _subset_col1:
                            _new_table_name = st.text_input(
                                "New table name:",
                                placeholder="e.g. my_subset",
                                key=f"subset_table_name_{display_selected_table}",
                            )
                        with _subset_col2:
                            st.write("")
                            st.write("")
                            if st.button(
                                f"Create Subset selected columns ({_n_sel_display} rows)",
                                key=f"btn_subset_table_{display_selected_table}",
                                disabled=(_n_sel_display == 0 or not _new_table_name.strip()),
                            ):
                                _row_indices = list(_selected_display_rows.index)
                                result = st_funcs.create_subset_table(
                                    db_path, display_selected_table,
                                    _new_table_name.strip(), _row_indices,
                                    columns=selected_display_cols
                                )
                                if result == "created":
                                    st.success(f"✅ Subset table **'{_new_table_name.strip()}'** created with {_n_sel_display} row(s) and {len(selected_display_cols)} column(s).")
                                elif result == "duplicate":
                                    st.error(f"A table named **'{_new_table_name.strip()}'** already exists. Choose a different name.")
                                else:
                                    st.error(f"Could not create subset: {result}")
                        with _subset_col3:
                            st.write("")
                            st.write("")
                            if st.button(
                                f"Create Subset All columns ({_n_sel_display} rows)",
                                key=f"btn_subset_all_cols_{display_selected_table}",
                                disabled=(_n_sel_display == 0 or not _new_table_name.strip()),
                            ):
                                _row_indices = list(_selected_display_rows.index)
                                result = st_funcs.create_subset_table(
                                    db_path, display_selected_table,
                                    _new_table_name.strip(), _row_indices,
                                    columns=None
                                )
                                if result == "created":
                                    st.success(f"✅ Subset table **'{_new_table_name.strip()}'** created with {_n_sel_display} row(s) and all {len(display_columns)} column(s).")
                                elif result == "duplicate":
                                    st.error(f"A table named **'{_new_table_name.strip()}'** already exists. Choose a different name.")
                                else:
                                    st.error(f"Could not create subset: {result}")
                    else:
                        st.info(f"Table '{display_selected_table}' is empty or could not be read.")
                else:
                    st.warning("Please select at least one column.")
            else:
                st.info(f"Table '{display_selected_table}' is empty or could not be read.")
    else:
        st.info("No tables found in the ChemSpace database for the active project.")

    ## Depict a selected table inline
    st.divider()
    st.subheader("Depict Table")

    if "show_depiction" not in st.session_state:
        st.session_state["show_depiction"] = False
    if "depiction_images" not in st.session_state:
        st.session_state["depiction_images"] = []
    if "show_depict_options" not in st.session_state:
        st.session_state["show_depict_options"] = False

    if df is not None and not df.empty:
        if st.button(f"{'Hide' if st.session_state['show_depict_options'] else 'Show'} Depict Table", key="btn_depict_table"):
            st.session_state["show_depict_options"] = not st.session_state["show_depict_options"]
            if not st.session_state["show_depict_options"]:
                st.session_state["show_depiction"] = False
                st.session_state["depiction_images"] = []
            st.rerun()

        if st.session_state["show_depict_options"]:
            table_names = df["table"].tolist()
            selected_table = st.selectbox("Select a table to depict:", table_names, key="depict_table_select")

            ## Load columns for the selected table to let the user pick the label
            table_columns = st_funcs.get_table_columns(db_path, selected_table, mtime=_chemspace_table_sigs.get(selected_table))
            default_label = "id" if "id" in table_columns else (table_columns[0] if table_columns else None)
            default_labels = [default_label] if default_label else []

            label_cols = st.multiselect(
                "Column(s) to use as molecule label:",
                table_columns,
                default=default_labels,
                key="depict_label_cols",
                help="If multiple columns are selected, the label is composed by joining their values with '_'."
            ) if table_columns else []

            col1, col2, col3 = st.columns(3)
            with col1:
                max_mols = st.number_input("Max molecules (-1 for all):", min_value=-1, value=25, step=1, key="depict_max_mols")
            with col2:
                mols_per_image = st.number_input("Molecules per grid image:", min_value=1, max_value=100, value=25, step=1, key="depict_mols_per_image")
            with col3:
                mol_size = st.number_input("Molecule cell size (px):", min_value=100, max_value=600, value=300, step=50, key="depict_mol_size")

            button_label = "Hide Depictions" if st.session_state["show_depiction"] else "Depict Selected Table"
            if st.button(button_label, key="btn_depict_selected"):
                if st.session_state["show_depiction"]:
                    st.session_state["show_depiction"] = False
                    st.session_state["depiction_images"] = []
                    st.rerun()
                else:
                    with st.spinner(f"Generating depictions for table '{selected_table}'…"):
                        try:
                            images = st_funcs.depict_table_to_images(
                                db_path=db_path,
                                table_name=selected_table,
                                max_molecules=int(max_mols),
                                molecules_per_image=int(mols_per_image),
                                mol_image_size=(int(mol_size), int(mol_size)),
                                legend_cols=label_cols
                            )
                            if images:
                                st.session_state["depiction_images"] = images
                                st.session_state["show_depiction"] = True
                                st.rerun()
                            else:
                                st.warning(f"No valid molecules found in table '{selected_table}'.")
                        except Exception as e:
                            st.error(f"Depiction failed: {e}")

            if st.session_state["show_depiction"] and st.session_state["depiction_images"]:
                images = st.session_state["depiction_images"]
                st.success(f"Showing {len(images)} grid image(s) for '{selected_table}'.")
                for i, img in enumerate(images):
                    st.image(img, caption=f"Grid {i + 1}", use_container_width=True)
    else:
        st.info("No tables found in the ChemSpace database for the active project.")

    ## Export a selected table to CSV
    st.divider()
    st.subheader("Export Table")

    if "show_export_options" not in st.session_state:
        st.session_state["show_export_options"] = False
    if "export_save_path" not in st.session_state:
        st.session_state["export_save_path"] = ""

    if df is not None and not df.empty:
        if st.button(f"{'Hide' if st.session_state['show_export_options'] else 'Show'} Export Table", key="btn_export_table"):
            st.session_state["show_export_options"] = not st.session_state["show_export_options"]
            if not st.session_state["show_export_options"]:
                st.session_state["export_save_path"] = ""
                st.session_state["show_export_path_input"] = False
            st.rerun()

        if st.session_state["show_export_options"]:
            export_table_names = df["table"].tolist()
            export_selected_table = st.selectbox(
                "Select a table to export:",
                export_table_names,
                key="export_table_select"
            )

            export_columns = st_funcs.get_table_columns(db_path, export_selected_table, mtime=_chemspace_table_sigs.get(export_selected_table))
            if export_columns:
                st.markdown("**Select columns to export:**")
                selected_cols = [col for col in export_columns if st.checkbox(col, value=True, key=f"export_col_{export_selected_table}_{col}")]

                if st.button("Save Table", key="btn_save_table_csv"):
                    if not selected_cols:
                        st.warning("Please select at least one column.")
                    else:
                        st.session_state["export_save_path"] = st.session_state.get("export_save_path", "")
                        st.session_state["show_export_path_input"] = True

                if st.session_state.get("show_export_path_input"):
                    save_path = st.text_input(
                        "Enter full path for the CSV file (e.g. /home/user/output.csv):",
                        value=st.session_state["export_save_path"],
                        key="export_path_input"
                    )
                    if st.button("Confirm Save", key="btn_confirm_save_csv"):
                        if not save_path.strip():
                            st.warning("Please enter a valid file path.")
                        else:
                            try:
                                export_df = st_funcs.read_table_columns_as_dataframe(db_path, export_selected_table, selected_cols, mtime=_chemspace_table_sigs.get(export_selected_table))
                                export_df.to_csv(save_path.strip(), index=False)
                                st.success(f"Table '{export_selected_table}' saved to: {save_path.strip()}")
                                st.session_state["show_export_path_input"] = False
                                st.session_state["export_save_path"] = ""
                            except Exception as e:
                                st.error(f"Export failed: {e}")
            else:
                st.warning(f"No columns found for table '{export_selected_table}'.")
    else:
        st.info("No tables found in the ChemSpace database for the active project.")

elif page == "ChemSpace Actions":
    st.title("ChemSpace Actions")
    st.write("Welcome to the ChemSpace Actions page.")

    st.subheader("Chemical Filtering: functional group based")

    with st.expander("🧪 Available Chemical Filters"):
        chem_filters = tidyscreen.get_chemical_filters()
        if chem_filters:
            chem_filters_df = pd.DataFrame(chem_filters, columns=["ID", "Filter Name", "SMARTS Pattern"])
            st.dataframe(chem_filters_df, use_container_width=True, hide_index=True)
        else:
            st.info("No chemical filters found in database.")

        st.markdown("**Add Filter**")
        _acf_col1, _acf_col2 = st.columns(2)
        with _acf_col1:
            _acf_name = st.text_input("Filter name", key="acf_filter_name")
        with _acf_col2:
            _acf_smarts = st.text_input("SMARTS pattern", key="acf_smarts_pattern")

        _acf_overwrite_key = "acf_confirm_overwrite"
        if not st.session_state.get(_acf_overwrite_key):
            if st.button("Add Filter", key="btn_add_chem_filter"):
                _acf_result = tidyscreen.add_chemical_filter_entry(_acf_name, _acf_smarts, overwrite=False)
                if _acf_result["exists"]:
                    st.session_state[_acf_overwrite_key] = True
                    st.rerun()
                elif _acf_result["success"]:
                    st.success(f"✅ {_acf_result['message']}")
                    st.rerun()
                else:
                    st.error(f"❌ {_acf_result['message']}")
        else:
            st.warning(f"⚠️ Filter '{_acf_name}' already exists. Overwrite it?")
            _acf_oc1, _acf_oc2 = st.columns(2)
            with _acf_oc1:
                if st.button("Yes, overwrite", key="btn_confirm_overwrite_filter"):
                    _acf_result = tidyscreen.add_chemical_filter_entry(_acf_name, _acf_smarts, overwrite=True)
                    st.session_state[_acf_overwrite_key] = False
                    if _acf_result["success"]:
                        st.success(f"✅ {_acf_result['message']}")
                    else:
                        st.error(f"❌ {_acf_result['message']}")
                    st.rerun()
            with _acf_oc2:
                if st.button("Cancel", key="btn_cancel_overwrite_filter"):
                    st.session_state[_acf_overwrite_key] = False
                    st.rerun()

        st.markdown("**Modify Filter**")
        if not chem_filters:
            st.caption("No filters available to modify.")
        else:
            _mcf_names = [name for _id, name, _smarts in chem_filters]
            _mcf_selected_name = st.selectbox("Select filter to modify", _mcf_names, key="mcf_selected_filter")
            _mcf_smarts_by_name = {name: smarts for _id, name, smarts in chem_filters}

            _mcf_col1, _mcf_col2 = st.columns(2)
            with _mcf_col1:
                _mcf_new_name = st.text_input(
                    "Filter name", value=_mcf_selected_name, key=f"mcf_new_name_{_mcf_selected_name}"
                )
            with _mcf_col2:
                _mcf_new_smarts = st.text_input(
                    "SMARTS pattern",
                    value=_mcf_smarts_by_name.get(_mcf_selected_name, ""),
                    key=f"mcf_new_smarts_{_mcf_selected_name}",
                )

            _mcf_overwrite_key = "mcf_confirm_overwrite"
            if not st.session_state.get(_mcf_overwrite_key):
                if st.button("Save Changes", key="btn_modify_chem_filter"):
                    _mcf_result = tidyscreen.modify_chemical_filter_entry(
                        _mcf_selected_name, _mcf_new_name, _mcf_new_smarts, overwrite=False
                    )
                    if _mcf_result["exists"]:
                        st.session_state[_mcf_overwrite_key] = True
                        st.rerun()
                    elif _mcf_result["success"]:
                        st.success(f"✅ {_mcf_result['message']}")
                        st.rerun()
                    else:
                        st.error(f"❌ {_mcf_result['message']}")
            else:
                st.warning(f"⚠️ Filter '{_mcf_new_name}' already exists. Overwrite it?")
                _mcf_oc1, _mcf_oc2 = st.columns(2)
                with _mcf_oc1:
                    if st.button("Yes, overwrite", key="btn_confirm_overwrite_modify_filter"):
                        _mcf_result = tidyscreen.modify_chemical_filter_entry(
                            _mcf_selected_name, _mcf_new_name, _mcf_new_smarts, overwrite=True
                        )
                        st.session_state[_mcf_overwrite_key] = False
                        if _mcf_result["success"]:
                            st.success(f"✅ {_mcf_result['message']}")
                        else:
                            st.error(f"❌ {_mcf_result['message']}")
                        st.rerun()
                with _mcf_oc2:
                    if st.button("Cancel", key="btn_cancel_overwrite_modify_filter"):
                        st.session_state[_mcf_overwrite_key] = False
                        st.rerun()

    with st.expander("🛠️ Create Filtering Workflow"):
        if "cfw_filters" not in st.session_state:
            st.session_state["cfw_filters"] = {}

        all_filters = tidyscreen.get_chemical_filters()

        if not all_filters:
            st.info("No chemical filters found in database. Add some first.")
        else:
            filters_by_name = {name: (filter_id, smarts) for filter_id, name, smarts in all_filters}
            available_names = [name for name in filters_by_name if name not in st.session_state["cfw_filters"]]

            st.markdown("**Add a filter to the workflow**")
            col1, col2, col3 = st.columns([3, 1, 1])
            with col1:
                cfw_selected_name = st.selectbox(
                    "Filter", available_names, key="cfw_selected_filter"
                ) if available_names else None
            with col2:
                cfw_instances = st.number_input(
                    "Required instances", min_value=0, value=1, step=1, key="cfw_instances"
                )
            with col3:
                st.markdown("<br>", unsafe_allow_html=True)
                cfw_add_disabled = not available_names
            cfw_step_description = st.text_input(
                "Description (rationale for adding it, optional)", key="cfw_step_description"
            )
            if st.button("Add Filter", key="btn_cfw_add_filter", disabled=cfw_add_disabled):
                filter_id, smarts = filters_by_name[cfw_selected_name]
                if ChemSpace._validate_smarts_pattern(smarts):
                    st.session_state["cfw_filters"][cfw_selected_name] = {
                        "instances": int(cfw_instances),
                        "smarts": smarts,
                        "filter_id": filter_id,
                        "description": cfw_step_description.strip(),
                    }
                    st.rerun()
                else:
                    st.error(f"Invalid SMARTS pattern for filter '{cfw_selected_name}'")

        st.markdown("**Current workflow**")
        if st.session_state["cfw_filters"]:
            cfw_df = pd.DataFrame([
                {
                    "Filter Name": name,
                    "Required Instances": info["instances"],
                    "SMARTS": info["smarts"],
                    "Description": info.get("description") or "—",
                }
                for name, info in st.session_state["cfw_filters"].items()
            ])
            st.dataframe(cfw_df, use_container_width=True, hide_index=True)

            cfw_remove_name = st.selectbox(
                "Remove a filter", list(st.session_state["cfw_filters"].keys()), key="cfw_remove_filter"
            )
            if st.button("Remove Filter", key="btn_cfw_remove_filter"):
                del st.session_state["cfw_filters"][cfw_remove_name]
                st.rerun()
        else:
            st.info("No filters added to the workflow yet.")

        st.markdown("**Save workflow**")
        cfw_workflow_name = st.text_input("Workflow name", key="cfw_workflow_name")
        cfw_description = st.text_area("Description (optional)", key="cfw_description")
        cfw_overwrite = st.checkbox("Overwrite if a workflow with this name already exists", key="cfw_overwrite")

        cfw_save_disabled = not st.session_state["cfw_filters"] or not cfw_workflow_name.strip()
        if not cfw_workflow_name.strip():
            st.caption("Enter a workflow name to save.")

        if st.button("Save Workflow", key="btn_cfw_save_workflow", disabled=cfw_save_disabled):
            chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
            workflow_filters = {name: info["instances"] for name, info in st.session_state["cfw_filters"].items()}
            filter_descriptions = {name: info.get("description", "") for name, info in st.session_state["cfw_filters"].items()}
            result = ChemSpace.save_filtering_workflow(
                chemspace_db_path, cfw_workflow_name, workflow_filters,
                cfw_description or None, cfw_overwrite,
                filter_descriptions=filter_descriptions
            )
            if result["success"]:
                st.success(f"Saved filtering workflow '{result['workflow_name']}' "
                           f"({result['filter_count']} filters, {result['total_instances']} total instances).")
                st.session_state["cfw_filters"] = {}
                st.rerun()
            else:
                st.error(result["message"])

    with st.expander("📋 List Filtering Workflows"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        filtering_workflows = ChemSpace.get_filtering_workflows(chemspace_db_path)
        if filtering_workflows:
            filtering_workflows_df = pd.DataFrame(filtering_workflows)
            st.dataframe(filtering_workflows_df, use_container_width=True, hide_index=True)
        else:
            st.info("No saved filtering workflows found.")

        st.markdown("**Delete filtering workflow**")
        if not filtering_workflows:
            st.button("🗑️ Delete filtering workflow", key="btn_delete_fw", disabled=True,
                       help="No filtering workflows available to delete.")
        else:
            dfw_workflow_names = [w["workflow_name"] for w in filtering_workflows]
            dfw_selected_name = st.selectbox("Select filtering workflow", dfw_workflow_names, key="dfw_selected_workflow")

            _del_fw_key = f"confirm_delete_fw_{dfw_selected_name}"
            if not st.session_state.get(_del_fw_key):
                if st.button("🗑️ Delete filtering workflow", key="btn_delete_fw"):
                    st.session_state[_del_fw_key] = True
                    st.rerun()
            else:
                st.warning(f"Delete filtering workflow **{dfw_selected_name}**?")
                _dfw_c1, _dfw_c2 = st.columns(2)
                with _dfw_c1:
                    if st.button("Yes, delete", key=f"btn_confirm_delete_fw_{dfw_selected_name}"):
                        _dfw_result = ChemSpace.delete_filtering_workflow_entry(chemspace_db_path, dfw_selected_name)
                        st.session_state[_del_fw_key] = False
                        if _dfw_result["success"]:
                            st.success(f"✅ {_dfw_result['message']}")
                        else:
                            st.error(f"❌ {_dfw_result['message']}")
                        st.rerun()
                with _dfw_c2:
                    if st.button("Cancel", key=f"btn_cancel_delete_fw_{dfw_selected_name}"):
                        st.session_state[_del_fw_key] = False
                        st.rerun()

            st.markdown("**Export filtering workflow**")
            if st.button(
                f"{'Hide Export' if st.session_state.get('show_export_fw_form') else '📤 Export filtering workflow'}: {dfw_selected_name}",
                key="btn_export_fw",
            ):
                st.session_state["show_export_fw_form"] = not st.session_state.get("show_export_fw_form", False)
                st.rerun()

            if st.session_state.get("show_export_fw_form"):
                _default_exp_fw_path = os.path.join(
                    st.session_state["active_project_path"],
                    "chemspace", "misc",
                    f"{dfw_selected_name}.json",
                )
                _exp_fw_path = st.text_input(
                    "Output JSON path:",
                    value=_default_exp_fw_path,
                    key=f"export_fw_path_{dfw_selected_name}",
                )
                if st.button("💾 Save", key=f"btn_save_fw_{dfw_selected_name}"):
                    try:
                        _all_filters = tidyscreen.get_chemical_filters()
                        _smarts_by_name = {_name: _smarts for _id, _name, _smarts in _all_filters}
                        _selected_wf = next(w for w in filtering_workflows if w["workflow_name"] == dfw_selected_name)
                        _payload = [
                            {"filter_name": _fname, "smarts": _smarts_by_name.get(_fname, "")}
                            for _fname in _selected_wf["filters_dict"]
                        ]
                        _out = _exp_fw_path.strip()
                        os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
                        with open(_out, "w", encoding="utf-8") as _f:
                            json.dump(_payload, _f, indent=2, ensure_ascii=False)
                        st.success(f"✅ Exported to: {_out}")
                        st.session_state["show_export_fw_form"] = False
                        st.rerun()
                    except Exception as _e:
                        st.error(f"❌ Export failed: {_e}")

    with st.expander("✏️ Edit Filtering Workflow"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        efw_workflows = ChemSpace.get_filtering_workflows(chemspace_db_path)

        if not efw_workflows:
            st.info("No saved filtering workflows found. Create one first.")
        else:
            efw_workflow_names = [w["workflow_name"] for w in efw_workflows]
            efw_selected_name = st.selectbox("Select workflow to edit", efw_workflow_names, key="efw_selected_workflow")
            efw_selected_wf = next(w for w in efw_workflows if w["workflow_name"] == efw_selected_name)

            if st.session_state.get("efw_loaded_workflow") != efw_selected_name:
                st.session_state["efw_filters"] = dict(efw_selected_wf["filters_dict"])
                st.session_state["efw_filter_descriptions"] = dict(efw_selected_wf.get("filter_descriptions") or {})
                st.session_state["efw_loaded_workflow"] = efw_selected_name

            efw_all_filters = tidyscreen.get_chemical_filters()

            if not efw_all_filters:
                st.info("No chemical filters found in database. Add some first.")
            else:
                efw_filters_by_name = {name: (filter_id, smarts) for filter_id, name, smarts in efw_all_filters}

                st.markdown("**Add or update a filter**")
                efw_col1, efw_col2, efw_col3 = st.columns([3, 1, 1])
                with efw_col1:
                    efw_selected_filter = st.selectbox(
                        "Filter", list(efw_filters_by_name.keys()), key="efw_selected_filter"
                    )
                with efw_col2:
                    _efw_existing_instances = st.session_state["efw_filters"].get(efw_selected_filter, 1)
                    efw_instances = st.number_input(
                        "Required instances", min_value=0, value=int(_efw_existing_instances), step=1,
                        key=f"efw_instances_{efw_selected_name}_{efw_selected_filter}"
                    )
                with efw_col3:
                    st.markdown("<br>", unsafe_allow_html=True)
                _efw_existing_description = st.session_state["efw_filter_descriptions"].get(efw_selected_filter, "")
                efw_step_description = st.text_input(
                    "Description (rationale for adding it, optional)", value=_efw_existing_description,
                    key=f"efw_description_{efw_selected_name}_{efw_selected_filter}"
                )
                if st.button("Add / Update Filter", key="btn_efw_add_filter"):
                    _, efw_smarts = efw_filters_by_name[efw_selected_filter]
                    if ChemSpace._validate_smarts_pattern(efw_smarts):
                        st.session_state["efw_filters"][efw_selected_filter] = int(efw_instances)
                        st.session_state["efw_filter_descriptions"][efw_selected_filter] = efw_step_description.strip()
                        st.rerun()
                    else:
                        st.error(f"Invalid SMARTS pattern for filter '{efw_selected_filter}'")

            st.markdown("**Current workflow parameters**")
            if st.session_state["efw_filters"]:
                efw_df = pd.DataFrame([
                    {
                        "Filter Name": name,
                        "Required Instances": instances,
                        "Description": st.session_state["efw_filter_descriptions"].get(name) or "—",
                    }
                    for name, instances in st.session_state["efw_filters"].items()
                ])
                st.dataframe(efw_df, use_container_width=True, hide_index=True)

                efw_remove_name = st.selectbox(
                    "Remove a filter", list(st.session_state["efw_filters"].keys()), key="efw_remove_filter"
                )
                if st.button("Remove Filter", key="btn_efw_remove_filter"):
                    del st.session_state["efw_filters"][efw_remove_name]
                    st.session_state["efw_filter_descriptions"].pop(efw_remove_name, None)
                    st.rerun()
            else:
                st.warning("No filters remaining in this workflow. Add at least one before saving.")

            st.markdown("**Save changes**")
            efw_description = st.text_area(
                "Description", value=efw_selected_wf["description"] or "",
                key=f"efw_description_input_{efw_selected_name}"
            )

            efw_save_disabled = not st.session_state["efw_filters"]
            if st.button("Save Changes", key="btn_efw_save", disabled=efw_save_disabled):
                efw_filter_descriptions = {
                    name: st.session_state["efw_filter_descriptions"].get(name, "")
                    for name in st.session_state["efw_filters"]
                }
                result = ChemSpace.save_filtering_workflow(
                    chemspace_db_path, efw_selected_name, st.session_state["efw_filters"],
                    efw_description or None, overwrite=True,
                    filter_descriptions=efw_filter_descriptions
                )
                if result["success"]:
                    st.success(f"Saved filtering workflow '{result['workflow_name']}' "
                               f"({result['filter_count']} filters, {result['total_instances']} total instances).")
                else:
                    st.error(result["message"])

    with st.expander("🔬 Filter Using Workflow"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        fuw_tables = [t for t in ChemSpace.list_tables_at_path(chemspace_db_path) if t != "filtering_workflows"]
        fuw_workflows = [w["workflow_name"] for w in ChemSpace.get_filtering_workflows(chemspace_db_path)]

        if not fuw_tables:
            st.info("No compound tables found in the ChemSpace database.")
        elif not fuw_workflows:
            st.info("No saved filtering workflows found. Create one first.")
        else:
            fuw_table_name = st.selectbox("Compound table", fuw_tables, key="fuw_table_name")
            fuw_workflow_name = st.selectbox("Filtering workflow", fuw_workflows, key="fuw_workflow_name")
            fuw_save_results = st.checkbox("Save filtered results to a new table", key="fuw_save_results")

            fuw_result_table_name = ""
            if fuw_save_results:
                fuw_result_table_name = st.text_input("Result table name", key="fuw_result_table_name")

            fuw_run_disabled = fuw_save_results and not fuw_result_table_name.strip()
            if fuw_run_disabled:
                st.caption("Enter a result table name to save the filtered results.")

            if st.button("Run Filtering", key="btn_fuw_run", disabled=fuw_run_disabled):
                # ActivateProject(name) looks up the project in this conda env's own
                # projects_database.db, which is separate from (and can be out of sync
                # with) the one the project was actually created in. Build a bare
                # ActivateProject carrying only what ChemSpace needs, using the path
                # already resolved by the Project selection page, instead of re-querying.
                project_obj = tidyscreen.ActivateProject.__new__(tidyscreen.ActivateProject)
                project_obj.name = st.session_state["selected_project"]
                project_obj.path = st.session_state["active_project_path"]
                project_obj._project_exists = True
                chemspace_obj = ChemSpace(project_obj)

                fuw_progress_bar = st.progress(0.0, text="Starting filtering workflow...")

                def _fuw_update_progress(done, total):
                    fraction = min(done / total, 1.0) if total else 1.0
                    fuw_progress_bar.progress(fraction, text=f"Filtering... {done}/{total}")

                fuw_log = io.StringIO()
                with contextlib.redirect_stdout(fuw_log):
                    filtered_df = chemspace_obj.filter_using_workflow(
                        table_name=fuw_table_name,
                        workflow_name=fuw_workflow_name,
                        save_results=fuw_save_results,
                        result_table_name=fuw_result_table_name.strip() if fuw_save_results else None,
                        progress_callback=_fuw_update_progress,
                    )
                fuw_progress_bar.progress(1.0, text="Filtering complete.")

                if filtered_df is not None and not filtered_df.empty:
                    st.success(f"Filtering complete: {len(filtered_df):,} compounds passed the workflow.")
                    st.dataframe(filtered_df, use_container_width=True, hide_index=True)
                else:
                    st.warning("No compounds passed the filtering workflow.")

                with st.expander("Execution log", expanded=filtered_df is None or filtered_df.empty):
                    st.code(fuw_log.getvalue() or "(no output)")

    st.subheader("Chemical Filtering: physicochemical properties based")

    with st.expander("🔬 Available Physicochemical Filters"):
        apf_descriptor_names = list(_get_descriptor_funcs().keys())

        apf_search = st.text_input("Search descriptors", key="apf_search")
        if apf_search.strip():
            apf_search_lower = apf_search.strip().lower()
            apf_rows = [(idx, name) for idx, name in enumerate(apf_descriptor_names, 1) if apf_search_lower in name.lower()]
        else:
            apf_rows = list(enumerate(apf_descriptor_names, 1))

        apf_df = pd.DataFrame(apf_rows, columns=["ID", "Descriptor Name"])
        st.dataframe(apf_df, use_container_width=True, hide_index=True)

    with st.expander("🛠️ Create Physicochemical Filtering Workflow"):
        if "cpfw_filters" not in st.session_state:
            st.session_state["cpfw_filters"] = {}

        cpfw_descriptor_names = list(_get_descriptor_funcs().keys())
        cpfw_available_names = [name for name in cpfw_descriptor_names if name not in st.session_state["cpfw_filters"]]

        st.markdown("**Add a descriptor to the workflow**")
        cpfw_col1, cpfw_col2, cpfw_col3 = st.columns(3)
        with cpfw_col1:
            cpfw_selected_name = st.selectbox(
                "Descriptor", cpfw_available_names, key="cpfw_selected_descriptor"
            ) if cpfw_available_names else None
        with cpfw_col2:
            cpfw_min_bound = st.number_input(
                "Lower bound", value=-9999.0, key="cpfw_min_bound"
            )
        with cpfw_col3:
            cpfw_max_bound = st.number_input(
                "Upper bound", value=9999.0, key="cpfw_max_bound"
            )

        if st.button("Add Descriptor", key="btn_cpfw_add_descriptor", disabled=not cpfw_available_names):
            if cpfw_min_bound > cpfw_max_bound:
                st.error(f"Lower bound ({cpfw_min_bound}) cannot be greater than upper bound ({cpfw_max_bound}).")
            else:
                st.session_state["cpfw_filters"][cpfw_selected_name] = {
                    "min": float(cpfw_min_bound),
                    "max": float(cpfw_max_bound),
                }
                st.rerun()

        st.markdown("**Current workflow**")
        if st.session_state["cpfw_filters"]:
            cpfw_df = pd.DataFrame([
                {"Descriptor Name": name, "Lower Bound": bounds["min"], "Upper Bound": bounds["max"]}
                for name, bounds in st.session_state["cpfw_filters"].items()
            ])
            st.dataframe(cpfw_df, use_container_width=True, hide_index=True)

            cpfw_remove_name = st.selectbox(
                "Remove a descriptor", list(st.session_state["cpfw_filters"].keys()), key="cpfw_remove_descriptor"
            )
            if st.button("Remove Descriptor", key="btn_cpfw_remove_descriptor"):
                del st.session_state["cpfw_filters"][cpfw_remove_name]
                st.rerun()
        else:
            st.info("No descriptors added to the workflow yet.")

        st.markdown("**Save workflow**")
        cpfw_workflow_name = st.text_input("Workflow name", key="cpfw_workflow_name")
        cpfw_description = st.text_area("Description (optional)", key="cpfw_description")
        cpfw_overwrite = st.checkbox("Overwrite if a workflow with this name already exists", key="cpfw_overwrite")

        cpfw_save_disabled = not st.session_state["cpfw_filters"] or not cpfw_workflow_name.strip()
        if not cpfw_workflow_name.strip():
            st.caption("Enter a workflow name to save.")

        if st.button("Save Workflow", key="btn_cpfw_save_workflow", disabled=cpfw_save_disabled):
            chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
            result = ChemSpace.save_physicochemical_filtering_workflow(
                chemspace_db_path, cpfw_workflow_name, st.session_state["cpfw_filters"],
                cpfw_description or None, cpfw_overwrite
            )
            if result["success"]:
                st.success(
                    f"✅ Successfully saved physicochemical filtering workflow!\n\n"
                    f"- **Workflow name:** {result['workflow_name']}\n"
                    f"- **Descriptors saved:** {result['filter_count']}\n"
                    f"- **Created:** {result['creation_date']}\n"
                    f"- **Description:** {result['description']}"
                )
                st.session_state["cpfw_filters"] = {}
                st.session_state["cpfw_last_saved"] = result
            else:
                st.error(result["message"])

        if st.session_state.get("cpfw_last_saved") and not st.session_state["cpfw_filters"]:
            _last = st.session_state["cpfw_last_saved"]
            st.info(
                f"Last saved workflow: **{_last['workflow_name']}** "
                f"({_last['filter_count']} descriptors, created {_last['creation_date']})"
            )

    with st.expander("📋 List Physicochemical Filtering Workflows"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        physicochemical_workflows = ChemSpace.get_physicochemical_filtering_workflows(chemspace_db_path)
        if physicochemical_workflows:
            physicochemical_workflows_df = pd.DataFrame(physicochemical_workflows)
            st.dataframe(physicochemical_workflows_df, use_container_width=True, hide_index=True)
        else:
            st.info("No saved physicochemical filtering workflows found.")

        st.markdown("**Delete physicochemical filtering workflow**")
        if not physicochemical_workflows:
            st.button("🗑️ Delete physicochemical filtering workflow", key="btn_delete_pfw", disabled=True,
                       help="No physicochemical filtering workflows available to delete.")
        else:
            dpfw_workflow_names = [w["workflow_name"] for w in physicochemical_workflows]
            dpfw_selected_name = st.selectbox("Select physicochemical filtering workflow", dpfw_workflow_names, key="dpfw_selected_workflow")

            _del_pfw_key = f"confirm_delete_pfw_{dpfw_selected_name}"
            if not st.session_state.get(_del_pfw_key):
                if st.button("🗑️ Delete physicochemical filtering workflow", key="btn_delete_pfw"):
                    st.session_state[_del_pfw_key] = True
                    st.rerun()
            else:
                st.warning(f"Delete physicochemical filtering workflow **{dpfw_selected_name}**?")
                _dpfw_c1, _dpfw_c2 = st.columns(2)
                with _dpfw_c1:
                    if st.button("Yes, delete", key=f"btn_confirm_delete_pfw_{dpfw_selected_name}"):
                        _dpfw_result = ChemSpace.delete_physicochemical_filtering_workflow_entry(chemspace_db_path, dpfw_selected_name)
                        st.session_state[_del_pfw_key] = False
                        if _dpfw_result["success"]:
                            st.success(f"✅ {_dpfw_result['message']}")
                        else:
                            st.error(f"❌ {_dpfw_result['message']}")
                        st.rerun()
                with _dpfw_c2:
                    if st.button("Cancel", key=f"btn_cancel_delete_pfw_{dpfw_selected_name}"):
                        st.session_state[_del_pfw_key] = False
                        st.rerun()

            st.markdown("**Export physicochemical filtering workflow**")
            if st.button(
                f"{'Hide Export' if st.session_state.get('show_export_pfw_form') else '📤 Export physicochemical filtering workflow'}: {dpfw_selected_name}",
                key="btn_export_pfw",
            ):
                st.session_state["show_export_pfw_form"] = not st.session_state.get("show_export_pfw_form", False)
                st.rerun()

            if st.session_state.get("show_export_pfw_form"):
                _default_exp_pfw_path = os.path.join(
                    st.session_state["active_project_path"],
                    "chemspace", "misc",
                    f"{dpfw_selected_name}.json",
                )
                _exp_pfw_path = st.text_input(
                    "Output JSON path:",
                    value=_default_exp_pfw_path,
                    key=f"export_pfw_path_{dpfw_selected_name}",
                )
                if st.button("💾 Save", key=f"btn_save_pfw_{dpfw_selected_name}"):
                    try:
                        _selected_pfw = next(w for w in physicochemical_workflows if w["workflow_name"] == dpfw_selected_name)
                        _payload = [
                            {"descriptor_name": _dname, "min": _bounds["min"], "max": _bounds["max"]}
                            for _dname, _bounds in _selected_pfw["filters_dict"].items()
                        ]
                        _out = _exp_pfw_path.strip()
                        os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
                        with open(_out, "w", encoding="utf-8") as _f:
                            json.dump(_payload, _f, indent=2, ensure_ascii=False)
                        st.success(f"✅ Exported to: {_out}")
                        st.session_state["show_export_pfw_form"] = False
                        st.rerun()
                    except Exception as _e:
                        st.error(f"❌ Export failed: {_e}")

    with st.expander("✏️ Edit Physicochemical Filtering Workflow"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        epfw_workflows = ChemSpace.get_physicochemical_filtering_workflows(chemspace_db_path)

        if not epfw_workflows:
            st.info("No saved physicochemical filtering workflows found. Create one first.")
        else:
            epfw_workflow_names = [w["workflow_name"] for w in epfw_workflows]
            epfw_selected_name = st.selectbox("Select workflow to edit", epfw_workflow_names, key="epfw_selected_workflow")
            epfw_selected_wf = next(w for w in epfw_workflows if w["workflow_name"] == epfw_selected_name)

            if st.session_state.get("epfw_loaded_workflow") != epfw_selected_name:
                st.session_state["epfw_filters"] = dict(epfw_selected_wf["filters_dict"])
                st.session_state["epfw_loaded_workflow"] = epfw_selected_name

            epfw_descriptor_names = list(_get_descriptor_funcs().keys())

            st.markdown("**Add or update a descriptor**")
            epfw_col1, epfw_col2, epfw_col3 = st.columns(3)
            with epfw_col1:
                epfw_selected_descriptor = st.selectbox("Descriptor", epfw_descriptor_names, key="epfw_selected_descriptor")
            _epfw_existing_bounds = st.session_state["epfw_filters"].get(epfw_selected_descriptor, {})
            with epfw_col2:
                epfw_min_bound = st.number_input(
                    "Lower bound", value=float(_epfw_existing_bounds.get("min", -9999.0)),
                    key=f"epfw_min_bound_{epfw_selected_name}_{epfw_selected_descriptor}"
                )
            with epfw_col3:
                epfw_max_bound = st.number_input(
                    "Upper bound", value=float(_epfw_existing_bounds.get("max", 9999.0)),
                    key=f"epfw_max_bound_{epfw_selected_name}_{epfw_selected_descriptor}"
                )

            if st.button("Add / Update Descriptor", key="btn_epfw_add_descriptor"):
                if epfw_min_bound > epfw_max_bound:
                    st.error(f"Lower bound ({epfw_min_bound}) cannot be greater than upper bound ({epfw_max_bound}).")
                else:
                    st.session_state["epfw_filters"][epfw_selected_descriptor] = {
                        "min": float(epfw_min_bound),
                        "max": float(epfw_max_bound),
                    }
                    st.rerun()

            st.markdown("**Current workflow parameters**")
            if st.session_state["epfw_filters"]:
                epfw_df = pd.DataFrame([
                    {"Descriptor Name": name, "Lower Bound": bounds["min"], "Upper Bound": bounds["max"]}
                    for name, bounds in st.session_state["epfw_filters"].items()
                ])
                st.dataframe(epfw_df, use_container_width=True, hide_index=True)

                epfw_remove_name = st.selectbox(
                    "Remove a descriptor", list(st.session_state["epfw_filters"].keys()), key="epfw_remove_descriptor"
                )
                if st.button("Remove Descriptor", key="btn_epfw_remove_descriptor"):
                    del st.session_state["epfw_filters"][epfw_remove_name]
                    st.rerun()
            else:
                st.warning("No descriptors remaining in this workflow. Add at least one before saving.")

            st.markdown("**Save changes**")
            epfw_description = st.text_area(
                "Description", value=epfw_selected_wf["description"] or "",
                key=f"epfw_description_input_{epfw_selected_name}"
            )

            epfw_save_disabled = not st.session_state["epfw_filters"]
            if st.button("Save Changes", key="btn_epfw_save", disabled=epfw_save_disabled):
                result = ChemSpace.save_physicochemical_filtering_workflow(
                    chemspace_db_path, epfw_selected_name, st.session_state["epfw_filters"],
                    epfw_description or None, overwrite=True
                )
                if result["success"]:
                    st.success(
                        f"✅ Successfully updated physicochemical filtering workflow!\n\n"
                        f"- **Workflow name:** {result['workflow_name']}\n"
                        f"- **Descriptors saved:** {result['filter_count']}\n"
                        f"- **Updated:** {result['creation_date']}\n"
                        f"- **Description:** {result['description']}"
                    )
                else:
                    st.error(result["message"])

    with st.expander("🔬 Filter Using Physicochemical Workflow"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        fupw_tables = [t for t in ChemSpace.list_tables_at_path(chemspace_db_path)
                       if t not in ("filtering_workflows", "physicochemical_filtering_workflows")]
        fupw_workflows = [w["workflow_name"] for w in ChemSpace.get_physicochemical_filtering_workflows(chemspace_db_path)]

        if not fupw_tables:
            st.info("No compound tables found in the ChemSpace database.")
        elif not fupw_workflows:
            st.info("No saved physicochemical filtering workflows found. Create one first.")
        else:
            fupw_table_name = st.selectbox("Compound table", fupw_tables, key="fupw_table_name")
            fupw_workflow_name = st.selectbox("Physicochemical workflow", fupw_workflows, key="fupw_workflow_name")
            fupw_save_results = st.checkbox("Save filtered results to a new table", key="fupw_save_results")

            fupw_result_table_name = ""
            if fupw_save_results:
                fupw_result_table_name = st.text_input("Result table name", key="fupw_result_table_name")

            fupw_run_disabled = fupw_save_results and not fupw_result_table_name.strip()
            if fupw_run_disabled:
                st.caption("Enter a result table name to save the filtered results.")

            if st.button("Run Filtering", key="btn_fupw_run", disabled=fupw_run_disabled):
                project_obj = tidyscreen.ActivateProject.__new__(tidyscreen.ActivateProject)
                project_obj.name = st.session_state["selected_project"]
                project_obj.path = st.session_state["active_project_path"]
                project_obj._project_exists = True
                chemspace_obj = ChemSpace(project_obj)

                fupw_progress_bar = st.progress(0.0, text="Starting physicochemical filtering workflow...")

                def _fupw_update_progress(done, total):
                    fraction = min(done / total, 1.0) if total else 1.0
                    fupw_progress_bar.progress(fraction, text=f"Filtering... {done}/{total}")

                fupw_log = io.StringIO()
                with contextlib.redirect_stdout(fupw_log):
                    filtered_df = chemspace_obj.filter_using_physicochemical_workflow(
                        table_name=fupw_table_name,
                        workflow_name=fupw_workflow_name,
                        save_results=fupw_save_results,
                        result_table_name=fupw_result_table_name.strip() if fupw_save_results else None,
                        progress_callback=_fupw_update_progress,
                    )
                fupw_progress_bar.progress(1.0, text="Filtering complete.")

                if filtered_df is not None and not filtered_df.empty:
                    st.success(f"Filtering complete: {len(filtered_df):,} compounds passed the workflow.")
                    st.dataframe(filtered_df, use_container_width=True, hide_index=True)
                else:
                    st.warning("No compounds passed the filtering workflow.")

                with st.expander("Execution log", expanded=filtered_df is None or filtered_df.empty):
                    st.code(fupw_log.getvalue() or "(no output)")

    st.divider()
    st.subheader("Synthetic workflows")

    with st.expander("🧬 Available Chemical Reactions"):
        chem_reactions = tidyscreen.get_chemical_reactions()
        if chem_reactions:
            chem_reactions_df = pd.DataFrame(chem_reactions, columns=["ID", "Reaction Name", "SMARTS Pattern"])
            st.dataframe(chem_reactions_df, use_container_width=True, hide_index=True)
        else:
            st.info("No chemical reactions found in database.")

        st.markdown("**Add chemical reaction**")
        _acr_col1, _acr_col2 = st.columns(2)
        with _acr_col1:
            _acr_name = st.text_input("Reaction name", key="acr_reaction_name")
        with _acr_col2:
            _acr_smarts = st.text_input("SMARTS pattern", key="acr_smarts_pattern")

        _acr_overwrite_key = "acr_confirm_overwrite"
        if not st.session_state.get(_acr_overwrite_key):
            if st.button("Add Reaction", key="btn_add_chem_reaction"):
                _acr_result = tidyscreen.add_chemical_reaction_entry(_acr_name, _acr_smarts, overwrite=False)
                if _acr_result["exists"]:
                    st.session_state[_acr_overwrite_key] = True
                    st.rerun()
                elif _acr_result["success"]:
                    st.success(f"✅ {_acr_result['message']}")
                    st.rerun()
                else:
                    st.error(f"❌ {_acr_result['message']}")
        else:
            st.warning(f"⚠️ Reaction '{_acr_name}' already exists. Overwrite it?")
            _acr_oc1, _acr_oc2 = st.columns(2)
            with _acr_oc1:
                if st.button("Yes, overwrite", key="btn_confirm_overwrite_reaction"):
                    _acr_result = tidyscreen.add_chemical_reaction_entry(_acr_name, _acr_smarts, overwrite=True)
                    st.session_state[_acr_overwrite_key] = False
                    if _acr_result["success"]:
                        st.success(f"✅ {_acr_result['message']}")
                    else:
                        st.error(f"❌ {_acr_result['message']}")
                    st.rerun()
            with _acr_oc2:
                if st.button("Cancel", key="btn_cancel_overwrite_reaction"):
                    st.session_state[_acr_overwrite_key] = False
                    st.rerun()

        st.markdown("**Delete chemical reaction**")
        if not chem_reactions:
            st.button("🗑️ Delete Reaction", key="btn_delete_cr", disabled=True,
                       help="No chemical reactions available to delete.")
        else:
            dcr_reaction_names = [name for _, name, _ in chem_reactions]
            dcr_selected_name = st.selectbox("Reaction to delete", dcr_reaction_names, key="dcr_selected_reaction")

            _del_cr_key = f"confirm_delete_cr_{dcr_selected_name}"
            if not st.session_state.get(_del_cr_key):
                if st.button("🗑️ Delete Reaction", key="btn_delete_cr"):
                    st.session_state[_del_cr_key] = True
                    st.rerun()
            else:
                st.warning(f"Delete chemical reaction **{dcr_selected_name}**?")
                _dcr_c1, _dcr_c2 = st.columns(2)
                with _dcr_c1:
                    if st.button("Yes, delete", key=f"btn_confirm_delete_cr_{dcr_selected_name}"):
                        _dcr_result = tidyscreen.delete_chemical_reaction_entry(dcr_selected_name)
                        st.session_state[_del_cr_key] = False
                        if _dcr_result["success"]:
                            st.success(f"✅ {_dcr_result['message']}")
                        else:
                            st.error(f"❌ {_dcr_result['message']}")
                        st.rerun()
                with _dcr_c2:
                    if st.button("Cancel", key=f"btn_cancel_delete_cr_{dcr_selected_name}"):
                        st.session_state[_del_cr_key] = False
                        st.rerun()

    with st.expander("🛠️ Create Reactions Workflow"):
        if "crw_reactions" not in st.session_state:
            st.session_state["crw_reactions"] = {}

        all_reactions = tidyscreen.get_chemical_reactions()

        if not all_reactions:
            st.info("No chemical reactions found in database. Add some first.")
        else:
            reactions_by_name = {name: (reaction_id, smarts) for reaction_id, name, smarts in all_reactions}
            crw_available_names = [name for name in reactions_by_name if name not in st.session_state["crw_reactions"]]

            st.markdown("**Add a reaction to the workflow**")
            crw_col1, crw_col2 = st.columns([3, 1])
            with crw_col1:
                crw_selected_name = st.selectbox(
                    "Reaction", crw_available_names, key="crw_selected_reaction"
                ) if crw_available_names else None
            with crw_col2:
                st.markdown("<br>", unsafe_allow_html=True)
                if st.button("Add Reaction", key="btn_crw_add_reaction", disabled=not crw_available_names):
                    reaction_id, smarts = reactions_by_name[crw_selected_name]
                    if ChemSpace._validate_reaction_smarts_pattern(smarts):
                        st.session_state["crw_reactions"][crw_selected_name] = {
                            "type": "reaction",
                            "reaction_id": reaction_id,
                            "smarts": smarts,
                            "order": len(st.session_state["crw_reactions"]) + 1,
                        }
                        st.rerun()
                    else:
                        st.error(f"Invalid SMARTS pattern for reaction '{crw_selected_name}'")

        st.markdown("**Add a table-merge step**")
        st.caption("Combines two tables (e.g. a product table with an existing one) into one, "
                   "usable as input for later reaction or merge steps. The two tables are chosen "
                   "when the workflow is applied.")
        crw_merge_col1, crw_merge_col2 = st.columns([3, 1])
        with crw_merge_col1:
            crw_merge_label = st.text_input(
                "Merge step label (optional)", key="crw_merge_label"
            )
        with crw_merge_col2:
            st.markdown("<br>", unsafe_allow_html=True)
            if st.button("Add Merge Step", key="btn_crw_add_merge"):
                crw_order = len(st.session_state["crw_reactions"]) + 1
                crw_label = crw_merge_label.strip() or f"Merge tables #{crw_order}"
                crw_step_key = crw_label
                crw_suffix = 1
                while crw_step_key in st.session_state["crw_reactions"]:
                    crw_suffix += 1
                    crw_step_key = f"{crw_label} ({crw_suffix})"
                st.session_state["crw_reactions"][crw_step_key] = {
                    "type": "merge",
                    "order": crw_order,
                }
                st.rerun()

        st.markdown("**Current workflow**")
        if st.session_state["crw_reactions"]:
            crw_df = pd.DataFrame([
                {
                    "Order": info["order"],
                    "Step Name": name,
                    "Type": info.get("type", "reaction"),
                    "SMARTS": info.get("smarts", "—"),
                }
                for name, info in sorted(
                    st.session_state["crw_reactions"].items(), key=lambda kv: kv[1]["order"]
                )
            ])
            st.dataframe(crw_df, use_container_width=True, hide_index=True)

            crw_remove_name = st.selectbox(
                "Remove a step", list(st.session_state["crw_reactions"].keys()), key="crw_remove_reaction"
            )
            if st.button("Remove Step", key="btn_crw_remove_reaction"):
                del st.session_state["crw_reactions"][crw_remove_name]
                st.rerun()
        else:
            st.info("No steps added to the workflow yet.")

        st.markdown("**Save workflow**")
        crw_workflow_name = st.text_input("Workflow name", key="crw_workflow_name")
        crw_description = st.text_area("Description (optional)", key="crw_description")
        crw_overwrite = st.checkbox("Overwrite if a workflow with this name already exists", key="crw_overwrite")

        crw_save_disabled = not st.session_state["crw_reactions"] or not crw_workflow_name.strip()
        if not crw_workflow_name.strip():
            st.caption("Enter a workflow name to save.")

        if st.button("Save Workflow", key="btn_crw_save_workflow", disabled=crw_save_disabled):
            chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
            workflow_reactions = {}
            for name, info in st.session_state["crw_reactions"].items():
                if info.get("type") == "merge":
                    workflow_reactions[f"merge_{info['order']}"] = {
                        "type": "merge", "name": name, "order": info["order"],
                    }
                else:
                    workflow_reactions[str(info["reaction_id"])] = {
                        "type": "reaction", "name": name, "smarts": info["smarts"], "order": info["order"],
                    }
            result = ChemSpace.save_reaction_workflow(
                chemspace_db_path, crw_workflow_name, workflow_reactions,
                crw_description or None, crw_overwrite
            )
            if result["success"]:
                st.success(f"Saved reaction workflow '{result['workflow_name']}' "
                           f"({result['reaction_count']} steps).")
                st.session_state["crw_reactions"] = {}
                st.rerun()
            else:
                st.error(result["message"])

    with st.expander("📋 List Chemical Reactions Workflows"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        reaction_workflows = ChemSpace.get_reaction_workflows(chemspace_db_path)
        if reaction_workflows:
            reaction_workflows_df = pd.DataFrame(reaction_workflows)
            st.dataframe(reaction_workflows_df, use_container_width=True, hide_index=True)
        else:
            st.info("No saved reaction workflows found.")

        st.markdown("**Delete reactions workflow**")
        if not reaction_workflows:
            st.button("🗑️ Delete Workflow", key="btn_delete_rw", disabled=True,
                       help="No reaction workflows available to delete.")
        else:
            drw_workflow_names = [w["workflow_name"] for w in reaction_workflows]
            drw_selected_name = st.selectbox("Workflow to delete", drw_workflow_names, key="drw_selected_workflow")

            _del_rw_key = f"confirm_delete_rw_{drw_selected_name}"
            if not st.session_state.get(_del_rw_key):
                if st.button("🗑️ Delete Workflow", key="btn_delete_rw"):
                    st.session_state[_del_rw_key] = True
                    st.rerun()
            else:
                st.warning(f"Delete reaction workflow **{drw_selected_name}**?")
                _drw_c1, _drw_c2 = st.columns(2)
                with _drw_c1:
                    if st.button("Yes, delete", key=f"btn_confirm_delete_rw_{drw_selected_name}"):
                        _drw_result = ChemSpace.delete_reaction_workflow_entry(chemspace_db_path, drw_selected_name)
                        st.session_state[_del_rw_key] = False
                        if _drw_result["success"]:
                            st.success(f"✅ {_drw_result['message']}")
                        else:
                            st.error(f"❌ {_drw_result['message']}")
                        st.rerun()
                with _drw_c2:
                    if st.button("Cancel", key=f"btn_cancel_delete_rw_{drw_selected_name}"):
                        st.session_state[_del_rw_key] = False
                        st.rerun()

        st.markdown("**Duplicate reactions workflow**")
        if not reaction_workflows:
            st.button("📄 Duplicate Workflow", key="btn_duplicate_rw", disabled=True,
                       help="No reaction workflows available to duplicate.")
        else:
            durw_workflow_names = [w["workflow_name"] for w in reaction_workflows]
            durw_selected_name = st.selectbox("Workflow to duplicate", durw_workflow_names, key="durw_selected_workflow")
            durw_selected_wf = next(w for w in reaction_workflows if w["workflow_name"] == durw_selected_name)
            durw_new_name = st.text_input(
                "New workflow name", value=f"{durw_selected_name}_copy", key=f"durw_new_name_{durw_selected_name}"
            )

            durw_disabled = not durw_new_name.strip()
            if st.button("📄 Duplicate Workflow", key="btn_duplicate_rw", disabled=durw_disabled):
                durw_result = ChemSpace.save_reaction_workflow(
                    chemspace_db_path, durw_new_name.strip(), durw_selected_wf["reactions_dict"],
                    durw_selected_wf["description"], overwrite=False
                )
                if durw_result["success"]:
                    st.success(f"Duplicated workflow '{durw_selected_name}' as '{durw_result['workflow_name']}' "
                               f"({durw_result['reaction_count']} steps).")
                    st.rerun()
                else:
                    st.error(durw_result["message"])

    with st.expander("✏️ Edit Chemical Reaction Workflow"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        erw_workflows = ChemSpace.get_reaction_workflows(chemspace_db_path)

        if not erw_workflows:
            st.info("No saved reaction workflows found. Create one first.")
        else:
            erw_workflow_names = [w["workflow_name"] for w in erw_workflows]
            erw_selected_name = st.selectbox("Select workflow to edit", erw_workflow_names, key="erw_selected_workflow")
            erw_selected_wf = next(w for w in erw_workflows if w["workflow_name"] == erw_selected_name)

            if st.session_state.get("erw_loaded_workflow") != erw_selected_name:
                erw_loaded_reactions = {}
                for erw_step_key, erw_step_info in erw_selected_wf["reactions_dict"].items():
                    erw_step = dict(erw_step_info)
                    if erw_step.get("type") != "merge":
                        erw_step["reaction_id"] = erw_step_key
                    erw_loaded_reactions[erw_step_info["name"]] = erw_step
                st.session_state["erw_reactions"] = erw_loaded_reactions
                st.session_state["erw_loaded_workflow"] = erw_selected_name

            erw_all_reactions = tidyscreen.get_chemical_reactions()

            if not erw_all_reactions:
                st.info("No chemical reactions found in database. Add some first.")
            else:
                erw_reactions_by_name = {name: (reaction_id, smarts) for reaction_id, name, smarts in erw_all_reactions}
                erw_available_names = [name for name in erw_reactions_by_name if name not in st.session_state["erw_reactions"]]

                st.markdown("**Add a reaction to the workflow**")
                erw_col1, erw_col2 = st.columns([3, 1])
                with erw_col1:
                    erw_selected_reaction = st.selectbox(
                        "Reaction", erw_available_names, key="erw_selected_reaction"
                    ) if erw_available_names else None
                with erw_col2:
                    st.markdown("<br>", unsafe_allow_html=True)
                    if st.button("Add Reaction", key="btn_erw_add_reaction", disabled=not erw_available_names):
                        erw_reaction_id, erw_smarts = erw_reactions_by_name[erw_selected_reaction]
                        if ChemSpace._validate_reaction_smarts_pattern(erw_smarts):
                            st.session_state["erw_reactions"][erw_selected_reaction] = {
                                "type": "reaction",
                                "reaction_id": erw_reaction_id,
                                "smarts": erw_smarts,
                                "order": len(st.session_state["erw_reactions"]) + 1,
                            }
                            st.rerun()
                        else:
                            st.error(f"Invalid SMARTS pattern for reaction '{erw_selected_reaction}'")

            st.markdown("**Add a table-merge step**")
            st.caption("Combines two tables (e.g. a product table with an existing one) into one, "
                       "usable as input for later reaction or merge steps. The two tables are chosen "
                       "when the workflow is applied.")
            erw_merge_col1, erw_merge_col2 = st.columns([3, 1])
            with erw_merge_col1:
                erw_merge_label = st.text_input(
                    "Merge step label (optional)", key="erw_merge_label"
                )
            with erw_merge_col2:
                st.markdown("<br>", unsafe_allow_html=True)
                if st.button("Add Merge Step", key="btn_erw_add_merge"):
                    erw_order = len(st.session_state["erw_reactions"]) + 1
                    erw_label = erw_merge_label.strip() or f"Merge tables #{erw_order}"
                    erw_step_key = erw_label
                    erw_suffix = 1
                    while erw_step_key in st.session_state["erw_reactions"]:
                        erw_suffix += 1
                        erw_step_key = f"{erw_label} ({erw_suffix})"
                    st.session_state["erw_reactions"][erw_step_key] = {
                        "type": "merge",
                        "order": erw_order,
                    }
                    st.rerun()

            st.markdown("**Current workflow parameters**")
            if st.session_state["erw_reactions"]:
                erw_df = pd.DataFrame([
                    {
                        "Order": info["order"],
                        "Step Name": name,
                        "Type": info.get("type", "reaction"),
                        "SMARTS": info.get("smarts", "—"),
                    }
                    for name, info in sorted(
                        st.session_state["erw_reactions"].items(), key=lambda kv: kv[1]["order"]
                    )
                ])
                st.dataframe(erw_df, use_container_width=True, hide_index=True)

                erw_remove_name = st.selectbox(
                    "Remove a step", list(st.session_state["erw_reactions"].keys()), key="erw_remove_reaction"
                )
                if st.button("Remove Step", key="btn_erw_remove_reaction"):
                    del st.session_state["erw_reactions"][erw_remove_name]
                    st.rerun()
            else:
                st.warning("No steps remaining in this workflow. Add at least one before saving.")

            st.markdown("**Save changes**")
            erw_description = st.text_area(
                "Description", value=erw_selected_wf["description"] or "",
                key=f"erw_description_input_{erw_selected_name}"
            )

            erw_save_disabled = not st.session_state["erw_reactions"]
            if st.button("Save Changes", key="btn_erw_save", disabled=erw_save_disabled):
                erw_workflow_reactions = {}
                for erw_name, erw_info in st.session_state["erw_reactions"].items():
                    if erw_info.get("type") == "merge":
                        erw_workflow_reactions[f"merge_{erw_info['order']}"] = {
                            "type": "merge", "name": erw_name, "order": erw_info["order"],
                        }
                    else:
                        erw_workflow_reactions[str(erw_info["reaction_id"])] = {
                            "type": "reaction", "name": erw_name, "smarts": erw_info["smarts"], "order": erw_info["order"],
                        }
                result = ChemSpace.save_reaction_workflow(
                    chemspace_db_path, erw_selected_name, erw_workflow_reactions,
                    erw_description or None, overwrite=True
                )
                if result["success"]:
                    st.success(f"Saved reaction workflow '{result['workflow_name']}' "
                               f"({result['reaction_count']} steps).")
                else:
                    st.error(result["message"])

    with st.expander("▶️ Apply Reaction Workflow"):
        chemspace_db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
        arw_workflows = ChemSpace.get_reaction_workflows(chemspace_db_path)

        if not arw_workflows:
            st.info("No saved reaction workflows found. Create one first.")
        else:
            arw_workflow_names = [w["workflow_name"] for w in arw_workflows]
            arw_selected_name = st.selectbox("Reaction workflow", arw_workflow_names, key="arw_selected_workflow")

            if st.button("Load Workflow", key="btn_arw_load_workflow"):
                arw_selected_workflow = next(w for w in arw_workflows if w["workflow_name"] == arw_selected_name)
                arw_sorted_steps = sorted(
                    arw_selected_workflow["reactions_dict"].items(), key=lambda kv: kv[1]["order"]
                )
                st.session_state["arw_run"] = {
                    "workflow_name": arw_selected_workflow["workflow_name"],
                    "steps": arw_sorted_steps,
                    "step_idx": 0,
                    "step_products": {},
                }
                st.rerun()

            arw_run = st.session_state.get("arw_run")
            if arw_run:
                st.markdown(f"**Applying workflow: '{arw_run['workflow_name']}'**")

                # ActivateProject(name) looks up the project in this conda env's own
                # projects_database.db, which is separate from (and can be out of sync
                # with) the one the project was actually created in. Build a bare
                # ActivateProject carrying only what ChemSpace needs, using the path
                # already resolved by the Project selection page, instead of re-querying.
                arw_project_obj = tidyscreen.ActivateProject.__new__(tidyscreen.ActivateProject)
                arw_project_obj.name = st.session_state["selected_project"]
                arw_project_obj.path = st.session_state["active_project_path"]
                arw_project_obj._project_exists = True
                arw_chemspace_obj = ChemSpace(arw_project_obj)

                arw_total_steps = len(arw_run["steps"])
                arw_step_idx = arw_run["step_idx"]

                if arw_step_idx >= arw_total_steps:
                    st.success(f"Workflow complete: {arw_total_steps} step(s) executed.")
                    arw_summary_rows = [
                        {
                            "Step": step,
                            "Reaction": info["reaction_name"],
                            "Products": info["product_count"],
                            "Table": info["table_name"],
                        }
                        for step, info in arw_run["step_products"].items()
                    ]
                    if arw_summary_rows:
                        st.dataframe(pd.DataFrame(arw_summary_rows), use_container_width=True, hide_index=True)
                    if st.button("Reset", key="btn_arw_reset"):
                        del st.session_state["arw_run"]
                        st.rerun()
                else:
                    _arw_reaction_id, arw_reaction_raw = arw_run["steps"][arw_step_idx]
                    arw_step_type = arw_reaction_raw.get("type", "reaction")

                    if arw_step_type == "merge":
                        arw_reaction_info = {"name": arw_reaction_raw["name"]}
                        st.markdown(f"**Step {arw_step_idx + 1}/{arw_total_steps}: {arw_reaction_info['name']}** (merge)")
                        st.caption("Merge step: combines two tables into one.")
                    else:
                        arw_reaction_info = {"name": arw_reaction_raw["name"], "smarts": arw_reaction_raw["smarts"]}
                        arw_reaction_type = arw_chemspace_obj._analyze_single_reaction_type(arw_reaction_info["smarts"])

                        st.markdown(f"**Step {arw_step_idx + 1}/{arw_total_steps}: {arw_reaction_info['name']}** ({arw_reaction_type})")
                        st.code(arw_reaction_info["smarts"])

                    # Available reactant sources: original tables plus product tables saved by earlier steps
                    arw_original_tables = [
                        t for t in arw_chemspace_obj.get_all_tables()
                        if t not in ("filtering_workflows", "reaction_workflows")
                    ]
                    arw_source_options = []
                    arw_source_labels = {}
                    for t in arw_original_tables:
                        arw_count = arw_chemspace_obj.get_compound_count(table_name=t)
                        arw_source_options.append(t)
                        arw_source_labels[t] = f"{t} ({arw_count} compounds) [Original]"
                    for arw_prev_step, arw_info in arw_run["step_products"].items():
                        if arw_prev_step < arw_step_idx + 1:
                            t = arw_info["table_name"]
                            arw_source_options.append(t)
                            arw_source_labels[t] = (
                                f"{t} ({arw_info['product_count']} products from "
                                f"Step {arw_prev_step}: {arw_info['reaction_name']}) [Products]"
                            )

                    if not arw_source_options:
                        st.warning("No reactant tables available for this step.")
                    elif arw_step_type == "merge":
                        st.caption("Select the two tables to combine.")
                        arw_mcol1, arw_mcol2 = st.columns(2)
                        with arw_mcol1:
                            arw_table_a = st.selectbox(
                                "Table A", arw_source_options,
                                format_func=lambda t: arw_source_labels[t], key=f"arw_merge_a_{arw_step_idx}"
                            )
                        with arw_mcol2:
                            arw_table_b = st.selectbox(
                                "Table B", arw_source_options,
                                format_func=lambda t: arw_source_labels[t], key=f"arw_merge_b_{arw_step_idx}"
                            )

                        arw_save_products = st.checkbox(
                            "Save merged compounds to a new table (required for use in later steps)",
                            value=True, key=f"arw_save_{arw_step_idx}"
                        )
                        arw_output_table_name = ""
                        if arw_save_products:
                            arw_default_name = f"step{arw_step_idx + 1}_{arw_reaction_info['name']}"
                            arw_output_table_name = st.text_input(
                                "Output table name", value=arw_default_name, key=f"arw_output_name_{arw_step_idx}"
                            )

                        arw_run_disabled = arw_save_products and not arw_output_table_name.strip()
                        if st.button("Run Step", key=f"btn_arw_run_step_{arw_step_idx}", disabled=arw_run_disabled):
                            arw_log = io.StringIO()
                            with contextlib.redirect_stdout(arw_log):
                                arw_result = arw_chemspace_obj.apply_reaction_workflow_merge_step(
                                    workflow_name=arw_run["workflow_name"],
                                    step_num=arw_step_idx + 1,
                                    table_a=arw_table_a,
                                    table_b=arw_table_b,
                                    output_table_name=arw_output_table_name.strip() if arw_save_products else None,
                                )

                            if arw_result["success"] and arw_result["output_table"]:
                                arw_run["step_products"][arw_step_idx + 1] = {
                                    "table_name": arw_result["output_table"],
                                    "product_count": arw_result["products_generated"],
                                    "reaction_name": arw_reaction_info["name"],
                                }

                            if arw_result["success"]:
                                st.success(f"Combined {arw_result['products_generated']} compound row(s).")
                                if arw_result["output_table"]:
                                    st.info(f"Saved to table '{arw_result['output_table']}'.")
                                elif arw_save_products:
                                    st.warning("No compounds were generated; nothing was saved.")
                                if arw_result["products"]:
                                    st.dataframe(
                                        pd.DataFrame(arw_result["products"][:20]),
                                        use_container_width=True, hide_index=True
                                    )
                            else:
                                st.error(f"Step failed, workflow halted: {arw_result['message']}")

                            with st.expander("Execution log", expanded=not arw_result["success"]):
                                st.code(arw_log.getvalue() or "(no output)")

                            if arw_result["success"]:
                                arw_run["step_idx"] += 1
                            st.session_state["arw_run"] = arw_run
                    else:
                        if arw_reaction_type == "bimolecular":
                            st.caption("Bimolecular reaction: select a primary and a secondary reactant table "
                                       "(the same table may be used for both).")
                            arw_col1, arw_col2 = st.columns(2)
                            with arw_col1:
                                arw_primary = st.selectbox(
                                    "Primary reactant table", arw_source_options,
                                    format_func=lambda t: arw_source_labels[t], key=f"arw_primary_{arw_step_idx}"
                                )
                            with arw_col2:
                                arw_secondary = st.selectbox(
                                    "Secondary reactant table", arw_source_options,
                                    format_func=lambda t: arw_source_labels[t], key=f"arw_secondary_{arw_step_idx}"
                                )
                            arw_primary_tables = [arw_primary]
                            arw_secondary_table = arw_secondary
                        else:
                            st.caption("Unimolecular reaction: select one or more reactant tables "
                                       "(their compounds are combined as the input set).")
                            arw_primary_tables = st.multiselect(
                                "Reactant table(s)", arw_source_options,
                                format_func=lambda t: arw_source_labels[t], key=f"arw_sources_{arw_step_idx}"
                            )
                            arw_secondary_table = None

                        arw_save_products = st.checkbox(
                            "Save products to a new table (required for use in later steps)",
                            value=True, key=f"arw_save_{arw_step_idx}"
                        )
                        arw_output_table_name = ""
                        if arw_save_products:
                            arw_default_name = f"step{arw_step_idx + 1}_{arw_reaction_info['name']}"
                            arw_output_table_name = st.text_input(
                                "Output table name", value=arw_default_name, key=f"arw_output_name_{arw_step_idx}"
                            )

                        arw_run_disabled = (
                            not arw_primary_tables
                            or (arw_reaction_type == "bimolecular" and not arw_secondary_table)
                            or (arw_save_products and not arw_output_table_name.strip())
                        )
                        if st.button("Run Step", key=f"btn_arw_run_step_{arw_step_idx}", disabled=arw_run_disabled):
                            arw_log = io.StringIO()
                            with contextlib.redirect_stdout(arw_log):
                                arw_result = arw_chemspace_obj.apply_reaction_workflow_step(
                                    reaction_info=arw_reaction_info,
                                    workflow_name=arw_run["workflow_name"],
                                    step_num=arw_step_idx + 1,
                                    primary_tables=arw_primary_tables,
                                    secondary_table=arw_secondary_table,
                                    output_table_name=arw_output_table_name.strip() if arw_save_products else None,
                                )

                            if arw_result["success"] and arw_result["output_table"]:
                                arw_run["step_products"][arw_step_idx + 1] = {
                                    "table_name": arw_result["output_table"],
                                    "product_count": arw_result["products_generated"],
                                    "reaction_name": arw_reaction_info["name"],
                                }

                            if arw_result["success"]:
                                st.success(f"Generated {arw_result['products_generated']} product(s).")
                                if arw_result["output_table"]:
                                    st.info(f"Saved to table '{arw_result['output_table']}'.")
                                elif arw_save_products:
                                    st.warning("No products were generated; nothing was saved.")
                                if arw_result["products"]:
                                    st.dataframe(
                                        pd.DataFrame(arw_result["products"][:20]),
                                        use_container_width=True, hide_index=True
                                    )
                                if arw_result["ambiguous_reactants"]:
                                    arw_ambig_lines = []
                                    for arw_amb in arw_result["ambiguous_reactants"]:
                                        if "reactant2_id" in arw_amb:
                                            arw_ambig_lines.append(
                                                f"- '{arw_amb['reactant1_name']}' (id={arw_amb['reactant1_id']}, "
                                                f"smiles={arw_amb.get('reactant1_smiles', 'n/a')}) + "
                                                f"'{arw_amb['reactant2_name']}' (id={arw_amb['reactant2_id']}, "
                                                f"smiles={arw_amb.get('reactant2_smiles', 'n/a')}): "
                                                f"{arw_amb['product_possibilities']} product possibilities"
                                            )
                                        else:
                                            arw_ambig_lines.append(
                                                f"- '{arw_amb['reactant_name']}' (id={arw_amb['reactant_id']}, "
                                                f"smiles={arw_amb.get('reactant_smiles', 'n/a')}): "
                                                f"{arw_amb['product_possibilities']} product possibilities"
                                            )
                                    st.warning(
                                        f"{len(arw_result['ambiguous_reactants'])} reactant(s)/pair(s) matched the "
                                        f"reaction at more than one site (multiple matching functional groups):\n\n"
                                        + "\n".join(arw_ambig_lines)
                                    )
                            else:
                                st.error(f"Step failed, workflow halted: {arw_result['message']}")

                            with st.expander("Execution log", expanded=not arw_result["success"]):
                                st.code(arw_log.getvalue() or "(no output)")

                            if arw_result["success"]:
                                arw_run["step_idx"] += 1
                            st.session_state["arw_run"] = arw_run
                            st.rerun()

elif page == "Receptors":
    st.title("Receptors")
    st.write(f"Welcome to the Receptors page for project: {st.session_state.get('selected_project', 'Unknown')}.")

    pdbs_db_path = os.path.join(st.session_state["active_project_path"], "docking", "receptors", "pdbs.db")
    
    receptors_db_path = os.path.join(st.session_state["active_project_path"], "docking", "receptors", "receptors.db")

    try:
        pdbs_models = st_funcs.get_pdbs_info(pdbs_db_path)
    except Exception:
        pdbs_models = None

    try:
        df_templates_models = st_funcs.get_pdb_templates_info(pdbs_db_path)
    except Exception:
        df_templates_models = None
        
    try:
        df_receptors_models = st_funcs.get_docking_models_info(receptors_db_path)
    except Exception:
        df_receptors_models = None

    if pdbs_models is None or pdbs_models.empty:
        st.markdown(
            f"<span style='font-size: 24px; color: orange;'>No receptor PDBs for project: </span><span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
            unsafe_allow_html=True
        )
    else:
        ## Create a button to show/hide the receptor PDBs DataFrame
        if "show_pdbs" not in st.session_state:
            st.session_state["show_pdbs"] = False

        if st.button(f"{'Hide' if st.session_state['show_pdbs'] else 'Show'} Receptor PDBs Details"):
            st.session_state["show_pdbs"] = not st.session_state["show_pdbs"]

        if st.session_state["show_pdbs"]:
            st.dataframe(pdbs_models)

        ## Export PDB Model toggle button
        if "show_export_pdb" not in st.session_state:
            st.session_state["show_export_pdb"] = False

        if st.button(f"{'Hide' if st.session_state['show_export_pdb'] else 'Export'} PDB Model", key="toggle_export_pdb"):
            st.session_state["show_export_pdb"] = not st.session_state["show_export_pdb"]

        if st.session_state["show_export_pdb"]:
            model_options = pdbs_models["pdb_model_name"].tolist()
            selected_model_name = st.selectbox(
                "Select PDB model to export:",
                model_options,
                key="export_pdb_model_select"
            )
            default_output_dir = os.path.join(
                st.session_state["active_project_path"], "docking", "receptors"
            )
            export_output_dir = st.text_input(
                "Output directory:",
                value=default_output_dir,
                key="export_pdb_output_dir"
            )
            if st.button("Export PDB", key="export_pdb_button"):
                selected_row = pdbs_models[pdbs_models["pdb_model_name"] == selected_model_name].iloc[0]
                file_id = selected_row["file_id"]
                success, result = st_funcs.export_pdb_model(pdbs_db_path, file_id, export_output_dir.strip())
                if success:
                    st.success(f"PDB exported to: {result}")
                else:
                    st.error(f"Export failed: {result}")

    if df_templates_models is None or df_templates_models.empty:
        st.markdown(
            f"<span style='font-size: 24px; color: orange;'>No receptor templates for project: </span><span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
            unsafe_allow_html=True
        )
    else:
        ## Create a button to show/hide the receptor templates DataFrame
        if "show_templates" not in st.session_state:
            st.session_state["show_templates"] = False

        if st.button(f"{'Hide' if st.session_state['show_templates'] else 'Show'} Receptor Templates Details"):
            st.session_state["show_templates"] = not st.session_state["show_templates"]

        if st.session_state["show_templates"]:
            st.dataframe(df_templates_models)

    if df_receptors_models is None or df_receptors_models.empty:
        st.markdown(
            f"<span style='font-size: 24px; color: orange;'>No docking receptor models for project: </span><span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
            unsafe_allow_html=True
        )
    else:
        ## Create a button to show/hide the docking receptor models DataFrame
        if "show_receptors_models" not in st.session_state:
            st.session_state["show_receptors_models"] = False

        if st.button(f"{'Hide' if st.session_state['show_receptors_models'] else 'Show'} Docking Receptor Models Details"):
            st.session_state["show_receptors_models"] = not st.session_state["show_receptors_models"]

        if st.session_state["show_receptors_models"]:
            st.dataframe(df_receptors_models)

elif page == "MolDock assays":
    st.title("MolDock assays")
    st.write("Welcome to the MolDock page.")

    ## --- Docking Methods ---
    st.subheader("Docking Methods")

    docking_methods_db_path = os.path.join(
        st.session_state["active_project_path"], "docking", "docking_registers", "docking_methods.db"
    )

    try:
        df_methods = st_funcs.get_docking_methods(docking_methods_db_path)
    except Exception:
        df_methods = None

    if df_methods is None or df_methods.empty:
        st.markdown(
            f"<span style='font-size: 24px; color: orange;'>No docking methods registered for project: </span>"
            f"<span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
            unsafe_allow_html=True
        )
    else:
        with st.expander("Docking Methods Table", expanded=False):
            _df_dm = df_methods[["id", "method_name", "docking_engine", "description", "created_date"]].copy()
            _df_dm.insert(0, "Select", False)
            _edited_dm = st.data_editor(
                _df_dm,
                column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                disabled=[c for c in _df_dm.columns if c != "Select"],
                hide_index=True,
                use_container_width=True,
                key="data_editor_docking_methods",
            )
            _sel_dm_rows = _edited_dm[_edited_dm["Select"]]
            st.session_state["export_dm_selected_name"] = (
                _sel_dm_rows.iloc[0]["method_name"] if len(_sel_dm_rows) == 1 else None
            )

            ## Clear export form when selection changes
            _export_dm_name = st.session_state.get("export_dm_selected_name")
            if st.session_state.get("_last_export_dm_name") != _export_dm_name:
                st.session_state["show_export_dm_form"] = False
                st.session_state["_last_export_dm_name"] = _export_dm_name

            ## Export / Import / Delete buttons (indented)
            _sp, _col_dm_btns = st.columns([1, 9])
            with _col_dm_btns:
                _col_exp, _col_imp, _col_dm_del = st.columns(3)

                with _col_exp:
                    if not _export_dm_name:
                        st.button(
                            "📤 Export Method",
                            key="btn_export_dm",
                            disabled=True,
                            help="Check exactly one method in the table to enable export.",
                        )
                    else:
                        if st.button(
                            f"{'Hide Export' if st.session_state.get('show_export_dm_form') else '📤 Export'}: {_export_dm_name}",
                            key="btn_export_dm",
                        ):
                            st.session_state["show_export_dm_form"] = not st.session_state.get("show_export_dm_form", False)

                        if st.session_state.get("show_export_dm_form"):
                            _full_row = df_methods[df_methods["method_name"] == _export_dm_name].iloc[0]
                            _default_exp_path = os.path.join(
                                st.session_state["active_project_path"],
                                "docking", "docking_registers",
                                f"{_export_dm_name}.json",
                            )
                            _exp_path = st.text_input(
                                "Output JSON path:",
                                value=_default_exp_path,
                                key=f"export_dm_path_{_export_dm_name}",
                            )
                            if st.button("💾 Save", key=f"btn_save_dm_{_export_dm_name}"):
                                try:
                                    _params = json.loads(_full_row["parameters"]) if _full_row["parameters"] else {}
                                    _lp = json.loads(_full_row["ligand_prep_params"]) if _full_row["ligand_prep_params"] else {}
                                    _payload = {
                                        "method_name": str(_full_row["method_name"]),
                                        "docking_engine": str(_full_row["docking_engine"]),
                                        "description": str(_full_row["description"] or ""),
                                        "parameters": _params,
                                        "ligand_prep_params": _lp,
                                    }
                                    _out = _exp_path.strip()
                                    os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
                                    with open(_out, "w", encoding="utf-8") as _f:
                                        json.dump(_payload, _f, indent=2, ensure_ascii=False)
                                    st.success(f"✅ Exported to: {_out}")
                                    st.session_state["show_export_dm_form"] = False
                                except Exception as _e:
                                    st.error(f"❌ Export failed: {_e}")

                with _col_imp:
                    if st.button(
                        f"{'Hide Import' if st.session_state.get('show_import_dm_form') else '📥 Import Method'}",
                        key="btn_import_dm",
                    ):
                        st.session_state["show_import_dm_form"] = not st.session_state.get("show_import_dm_form", False)

                    if st.session_state.get("show_import_dm_form"):
                        _imp_path = st.text_input(
                            "JSON file path:",
                            placeholder="/path/to/method.json",
                            key="import_dm_path",
                        )
                        if st.button("📥 Load", key="btn_load_dm"):
                            if not _imp_path.strip():
                                st.error("Please enter a file path.")
                            elif not os.path.exists(_imp_path.strip()):
                                st.error(f"File not found: {_imp_path.strip()}")
                            else:
                                try:
                                    with open(_imp_path.strip(), "r", encoding="utf-8") as _f:
                                        _imp_data = json.load(_f)
                                    _required = {"method_name", "docking_engine", "description", "parameters", "ligand_prep_params"}
                                    _missing = _required - set(_imp_data.keys())
                                    if _missing:
                                        st.error(f"JSON missing required fields: {', '.join(sorted(_missing))}")
                                    else:
                                        _docking_dir = os.path.join(
                                            st.session_state["active_project_path"],
                                            "docking", "docking_registers",
                                        )
                                        os.makedirs(_docking_dir, exist_ok=True)
                                        _imp_db = os.path.join(_docking_dir, "docking_methods.db")
                                        _nc = sqlite3.connect(_imp_db)
                                        try:
                                            _nc.execute("""
                                                CREATE TABLE IF NOT EXISTS docking_methods (
                                                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                                                    method_name TEXT UNIQUE NOT NULL,
                                                    docking_engine TEXT NOT NULL,
                                                    description TEXT,
                                                    parameters TEXT,
                                                    ligand_prep_params TEXT,
                                                    created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                                                )
                                            """)
                                            if _nc.execute(
                                                "SELECT COUNT(*) FROM docking_methods WHERE method_name = ?",
                                                (_imp_data["method_name"],),
                                            ).fetchone()[0] > 0:
                                                st.warning(f"Method '{_imp_data['method_name']}' already exists in this project.")
                                            else:
                                                _nc.execute("""
                                                    INSERT INTO docking_methods
                                                        (method_name, docking_engine, description, parameters, ligand_prep_params)
                                                    VALUES (?, ?, ?, ?, ?)
                                                """, (
                                                    _imp_data["method_name"],
                                                    _imp_data["docking_engine"],
                                                    _imp_data["description"],
                                                    json.dumps(_imp_data["parameters"], indent=2),
                                                    json.dumps(_imp_data["ligand_prep_params"], indent=2),
                                                ))
                                                _nc.commit()
                                                st.success(f"✅ Method '{_imp_data['method_name']}' imported successfully.")
                                                st.session_state["show_import_dm_form"] = False
                                                st.rerun()
                                        finally:
                                            _nc.close()
                                except json.JSONDecodeError as _e:
                                    st.error(f"Invalid JSON: {_e}")
                                except Exception as _e:
                                    st.error(f"Import failed: {_e}")

                with _col_dm_del:
                    _del_dm_key = f"confirm_delete_dm_{_export_dm_name}"
                    if not _export_dm_name:
                        st.button(
                            "🗑️ Delete Method",
                            key="btn_delete_dm",
                            disabled=True,
                            help="Check exactly one method in the table to enable deletion.",
                        )
                    elif not st.session_state.get(_del_dm_key):
                        if st.button("🗑️ Delete Method", key="btn_delete_dm"):
                            st.session_state[_del_dm_key] = True
                            st.rerun()
                    else:
                        st.warning(f"Delete **{_export_dm_name}**?")
                        _ddc1, _ddc2 = st.columns(2)
                        with _ddc1:
                            if st.button("Yes, delete", key=f"btn_confirm_delete_dm_{_export_dm_name}"):
                                try:
                                    _nc = sqlite3.connect(docking_methods_db_path)
                                    try:
                                        _nc.execute("DELETE FROM docking_methods WHERE method_name = ?", (_export_dm_name,))
                                        _nc.commit()
                                    finally:
                                        _nc.close()
                                    st.session_state[_del_dm_key] = False
                                    st.session_state["export_dm_selected_name"] = None
                                    st.rerun()
                                except Exception as _e:
                                    st.error(f"❌ Delete failed: {_e}")
                        with _ddc2:
                            if st.button("Cancel", key=f"btn_cancel_delete_dm_{_export_dm_name}"):
                                st.session_state[_del_dm_key] = False
                                st.rerun()

            ## Inspector — nested inside the same expander as the table above
            _export_dm_name = st.session_state.get("export_dm_selected_name")
            _sp, _col_dm_insp = st.columns([1, 9])
            with _col_dm_insp:
                method_names = df_methods["method_name"].dropna().unique().tolist()

                if st.session_state.get("selected_docking_method") not in method_names:
                    st.session_state["selected_docking_method"] = method_names[0] if method_names else None

                selected_method = st.selectbox(
                    "Select a docking method to inspect:",
                    method_names,
                    key="select_docking_method",
                    index=method_names.index(st.session_state["selected_docking_method"]) if st.session_state.get("selected_docking_method") in method_names else 0,
                )
                if selected_method != st.session_state.get("selected_docking_method"):
                    st.session_state["selected_docking_method"] = selected_method

                method_row = df_methods[df_methods["method_name"] == selected_method].iloc[0]

                st.markdown(f"### Method: `{method_row['method_name']}`")
                st.markdown(f"**Docking Engine:** {method_row['docking_engine']}")
                st.markdown(f"**Description:** {method_row['description'] or 'N/A'}")
                st.markdown(f"**Created:** {method_row['created_date']}")
                with st.expander("Docking Parameters", expanded=False):
                    st.code(method_row["parameters"] or "N/A", language="json")
                with st.expander("Ligand Preparation Parameters", expanded=False):
                    st.code(method_row["ligand_prep_params"] or "N/A", language="json")

    st.divider()

    ## --- Docking Assays ---
    pdbs_db_path = os.path.join(st.session_state["active_project_path"], "docking", "receptors", "pdbs.db")
    receptors_db_path = os.path.join(st.session_state["active_project_path"], "docking", "receptors", "receptors.db")
    docking_registries_db_path = os.path.join(st.session_state["active_project_path"], "docking", "docking_registers", "docking_assays.db")

    try:
        df = st_funcs.get_docking_assay_registers(docking_registries_db_path)
    except Exception:
        df = None

    if df is None or df.empty:
        st.markdown(
            f"<span style='font-size: 24px; color: orange;'>No docking assays for project: </span><span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
            unsafe_allow_html=True
        )
        
    else:
        with st.expander("Docking Assays Details", expanded=False):
            st.dataframe(df)

            # Selection box for assay_name
            if "assay_name" in df.columns:
                assay_names = df["assay_name"].dropna().unique().tolist()
                if "selected_assay_name" not in st.session_state and assay_names:
                    st.session_state["selected_assay_name"] = assay_names[0]
                selected_assay = st.selectbox(
                    "Select a docking assay:",
                    assay_names,
                    key="select_assay_name",
                    index=assay_names.index(st.session_state.get("selected_assay_name", assay_names[0])) if assay_names else 0
                )
                if selected_assay != st.session_state.get("selected_assay_name", None):
                    _clear_docking_assay_scoped_state()
                    st.session_state["selected_assay_name"] = selected_assay
                st.success(f"Selected assay: {selected_assay}")

        st.write(f"Selected Docking Assay: {st.session_state.get('selected_assay_name', 'None')}")


elif page == "MolDyn assays":
    st.title("MolDyn assays")
    st.write("Welcome to the MolDyn page.")

    if "active_project_path" not in st.session_state:
        st.warning("No active project. Please select a project first.")
    else:
        project_path = st.session_state["active_project_path"]

        ## --- MD Methods Registries ---
        st.subheader("MD Methods Registries")

        md_methods_db_path = os.path.join(project_path, "dynamics", "md_registers", "md_methods.db")

        try:
            df_md_methods = st_funcs.get_md_method_registers(md_methods_db_path)
        except Exception:
            df_md_methods = None

        if df_md_methods is None or df_md_methods.empty:
            st.markdown(
                f"<span style='font-size: 24px; color: orange;'>No MD methods registered for project: </span>"
                f"<span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
                unsafe_allow_html=True
            )
        else:
            with st.expander("MD Methods Table", expanded=False):
                _df_mdm = df_md_methods[["method_id", "method_name", "engine", "description"]].copy()
                _df_mdm.insert(0, "Select", False)
                _edited_mdm = st.data_editor(
                    _df_mdm,
                    column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                    disabled=[c for c in _df_mdm.columns if c != "Select"],
                    hide_index=True,
                    use_container_width=True,
                    key="data_editor_md_methods",
                )
                _sel_mdm_rows = _edited_mdm[_edited_mdm["Select"]]
                st.session_state["export_mdm_selected_name"] = (
                    _sel_mdm_rows.iloc[0]["method_name"] if len(_sel_mdm_rows) == 1 else None
                )

                ## Clear export form when selection changes
                _export_mdm_name = st.session_state.get("export_mdm_selected_name")
                if st.session_state.get("_last_export_mdm_name") != _export_mdm_name:
                    st.session_state["show_export_mdm_form"] = False
                    st.session_state["_last_export_mdm_name"] = _export_mdm_name

                ## Export / Import / Delete buttons (indented)
                _sp, _col_mdm_btns = st.columns([1, 9])
                with _col_mdm_btns:
                    _col_mdm_exp, _col_mdm_imp, _col_mdm_del = st.columns(3)

                    with _col_mdm_exp:
                        if not _export_mdm_name:
                            st.button(
                                "📤 Export Method",
                                key="btn_export_mdm",
                                disabled=True,
                                help="Check exactly one method in the table to enable export.",
                            )
                        else:
                            if st.button(
                                f"{'Hide Export' if st.session_state.get('show_export_mdm_form') else '📤 Export'}: {_export_mdm_name}",
                                key="btn_export_mdm",
                            ):
                                st.session_state["show_export_mdm_form"] = not st.session_state.get("show_export_mdm_form", False)

                            if st.session_state.get("show_export_mdm_form"):
                                _full_mdm_row = df_md_methods[df_md_methods["method_name"] == _export_mdm_name].iloc[0]
                                _default_mdm_exp_path = os.path.join(
                                    project_path, "dynamics", "md_registers",
                                    f"{_export_mdm_name}.json",
                                )
                                _mdm_exp_path = st.text_input(
                                    "Output JSON path:",
                                    value=_default_mdm_exp_path,
                                    key=f"export_mdm_path_{_export_mdm_name}",
                                )
                                if st.button("💾 Save", key=f"btn_save_mdm_{_export_mdm_name}"):
                                    try:
                                        _mdm_params = json.loads(_full_mdm_row["parameters"]) if _full_mdm_row["parameters"] else {}
                                        _mdm_payload = {
                                            "method_name": str(_full_mdm_row["method_name"]),
                                            "engine": str(_full_mdm_row["engine"]),
                                            "description": str(_full_mdm_row["description"] or ""),
                                            "parameters": _mdm_params,
                                        }
                                        _mdm_out = _mdm_exp_path.strip()
                                        os.makedirs(os.path.dirname(os.path.abspath(_mdm_out)), exist_ok=True)
                                        with open(_mdm_out, "w", encoding="utf-8") as _f:
                                            json.dump(_mdm_payload, _f, indent=2, ensure_ascii=False)
                                        st.success(f"✅ Exported to: {_mdm_out}")
                                        st.session_state["show_export_mdm_form"] = False
                                    except Exception as _e:
                                        st.error(f"❌ Export failed: {_e}")

                    with _col_mdm_imp:
                        if st.button(
                            f"{'Hide Import' if st.session_state.get('show_import_mdm_form') else '📥 Import Method'}",
                            key="btn_import_mdm",
                        ):
                            st.session_state["show_import_mdm_form"] = not st.session_state.get("show_import_mdm_form", False)

                        if st.session_state.get("show_import_mdm_form"):
                            _mdm_imp_path = st.text_input(
                                "JSON file path:",
                                placeholder="/path/to/md_method.json",
                                key="import_mdm_path",
                            )
                            if st.button("📥 Load", key="btn_load_mdm"):
                                if not _mdm_imp_path.strip():
                                    st.error("Please enter a file path.")
                                elif not os.path.exists(_mdm_imp_path.strip()):
                                    st.error(f"File not found: {_mdm_imp_path.strip()}")
                                else:
                                    try:
                                        with open(_mdm_imp_path.strip(), "r", encoding="utf-8") as _f:
                                            _mdm_imp_data = json.load(_f)
                                        _mdm_required = {"method_name", "engine", "description", "parameters"}
                                        _mdm_missing = _mdm_required - set(_mdm_imp_data.keys())
                                        if _mdm_missing:
                                            st.error(f"JSON missing required fields: {', '.join(sorted(_mdm_missing))}")
                                        else:
                                            _mdm_dir = os.path.join(project_path, "dynamics", "md_registers")
                                            os.makedirs(_mdm_dir, exist_ok=True)
                                            _mdm_imp_db = os.path.join(_mdm_dir, "md_methods.db")
                                            _nc = sqlite3.connect(_mdm_imp_db)
                                            try:
                                                _nc.execute("""
                                                    CREATE TABLE IF NOT EXISTS md_methods (
                                                        method_id INTEGER PRIMARY KEY AUTOINCREMENT,
                                                        method_name TEXT UNIQUE NOT NULL,
                                                        description TEXT,
                                                        engine TEXT,
                                                        parameters TEXT
                                                    )
                                                """)
                                                if _nc.execute(
                                                    "SELECT COUNT(*) FROM md_methods WHERE method_name = ?",
                                                    (_mdm_imp_data["method_name"],),
                                                ).fetchone()[0] > 0:
                                                    st.warning(f"Method '{_mdm_imp_data['method_name']}' already exists in this project.")
                                                else:
                                                    _nc.execute(
                                                        "INSERT INTO md_methods (method_name, engine, description, parameters) VALUES (?, ?, ?, ?)",
                                                        (
                                                            _mdm_imp_data["method_name"],
                                                            _mdm_imp_data["engine"],
                                                            _mdm_imp_data["description"],
                                                            json.dumps(_mdm_imp_data["parameters"], indent=2),
                                                        )
                                                    )
                                                    _nc.commit()
                                                    st.success(f"✅ Method '{_mdm_imp_data['method_name']}' imported successfully.")
                                                    st.session_state["show_import_mdm_form"] = False
                                                    st.rerun()
                                            finally:
                                                _nc.close()
                                    except json.JSONDecodeError as _e:
                                        st.error(f"Invalid JSON: {_e}")
                                    except Exception as _e:
                                        st.error(f"Import failed: {_e}")

                    with _col_mdm_del:
                        _del_mdm_key = f"confirm_delete_mdm_{_export_mdm_name}"
                        if not _export_mdm_name:
                            st.button(
                                "🗑️ Delete Method",
                                key="btn_delete_mdm",
                                disabled=True,
                                help="Check exactly one method in the table to enable deletion.",
                            )
                        elif not st.session_state.get(_del_mdm_key):
                            if st.button("🗑️ Delete Method", key="btn_delete_mdm"):
                                st.session_state[_del_mdm_key] = True
                                st.rerun()
                        else:
                            st.warning(f"Delete **{_export_mdm_name}**?")
                            _dmc1, _dmc2 = st.columns(2)
                            with _dmc1:
                                if st.button("Yes, delete", key=f"btn_confirm_delete_mdm_{_export_mdm_name}"):
                                    try:
                                        _nc = sqlite3.connect(md_methods_db_path)
                                        try:
                                            _nc.execute("DELETE FROM md_methods WHERE method_name = ?", (_export_mdm_name,))
                                            _nc.commit()
                                        finally:
                                            _nc.close()
                                        st.session_state[_del_mdm_key] = False
                                        st.session_state["export_mdm_selected_name"] = None
                                        st.rerun()
                                    except Exception as _e:
                                        st.error(f"❌ Delete failed: {_e}")
                            with _dmc2:
                                if st.button("Cancel", key=f"btn_cancel_delete_mdm_{_export_mdm_name}"):
                                    st.session_state[_del_mdm_key] = False
                                    st.rerun()

                ## Inspector — nested inside the same expander as the table above
                _sp, _col_mdm_insp = st.columns([1, 9])
                with _col_mdm_insp:
                    _mdm_names = df_md_methods["method_name"].dropna().unique().tolist()

                    if st.session_state.get("selected_md_method") not in _mdm_names:
                        st.session_state["selected_md_method"] = _mdm_names[0] if _mdm_names else None

                    _sel_md_method = st.selectbox(
                        "Select an MD method to inspect:",
                        _mdm_names,
                        key="select_md_method",
                        index=_mdm_names.index(st.session_state["selected_md_method"]) if st.session_state.get("selected_md_method") in _mdm_names else 0,
                    )
                    if _sel_md_method != st.session_state.get("selected_md_method"):
                        st.session_state["selected_md_method"] = _sel_md_method

                    _mdm_row = df_md_methods[df_md_methods["method_name"] == _sel_md_method].iloc[0]

                    st.markdown(f"### Method: `{_mdm_row['method_name']}`")
                    st.markdown(f"**Engine:** {_mdm_row['engine']}")
                    st.markdown(f"**Description:** {_mdm_row['description'] or 'N/A'}")
                    with st.expander("Parameters", expanded=False):
                        st.code(_mdm_row["parameters"] or "N/A", language="json")

        st.divider()

        ## --- MD Assay Registries ---
        st.subheader("MD Assay Registries")

        md_registers_db_path = os.path.join(project_path, "dynamics", "md_registers", "md_registers.db")

        try:
            df_md_assays = st_funcs.get_md_assay_registers(md_registers_db_path)
        except Exception:
            df_md_assays = None

        if df_md_assays is None or df_md_assays.empty:
            st.markdown(
                f"<span style='font-size: 24px; color: orange;'>No MD assays registered for project: </span>"
                f"<span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
                unsafe_allow_html=True
            )
        else:
            # Add receptor_template_name if not already present (query directly from DB)
            if 'receptor_template_name' not in df_md_assays.columns:
                try:
                    _conn = sqlite3.connect(md_registers_db_path)
                    _tpl = pd.read_sql_query(
                        "SELECT assay_id, receptor_template_name FROM md_assays", _conn
                    )
                    _conn.close()
                    df_md_assays = df_md_assays.merge(_tpl, on='assay_id', how='left')
                except Exception:
                    df_md_assays['receptor_template_name'] = None

            # Add status column if not already present
            if 'status' not in df_md_assays.columns:
                def _md_status(folder):
                    if not folder or not os.path.isdir(folder):
                        return '⬜ Not started'
                    out_files = ['min1.out', 'min2.out', 'heating.out', 'equilibration.out', 'production.out']
                    if not any(os.path.exists(os.path.join(folder, f)) for f in out_files):
                        return '⬜ Not started'
                    prod = os.path.join(folder, 'production.out')
                    if os.path.exists(prod):
                        try:
                            with open(prod, 'rb') as _fh:
                                _fh.seek(0, 2)
                                _fh.seek(max(0, _fh.tell() - 4096), 0)
                                _tail = _fh.read().decode('utf-8', errors='replace')
                                if 'Total wall time:' in _tail or 'FINAL RESULTS' in _tail:
                                    return '✅ Completed'
                        except OSError:
                            pass
                    return '🔄 Running'
                df_md_assays['status'] = df_md_assays['assay_folder_path'].apply(_md_status)

            # Add mmgbsa column if not already present (inline fallback, no st_funcs dependency)
            if 'mmgbsa' not in df_md_assays.columns:
                def _mmgbsa_status(row):
                    mmgbsa_json = row.get('mmgbsa_results')
                    folder = row.get('assay_folder_path')
                    if mmgbsa_json is not None and str(mmgbsa_json).strip():
                        return '✅ Completed'
                    if folder and os.path.isdir(str(folder)):
                        mmgbsa_dir = os.path.join(folder, 'mmgbsa')
                        if os.path.isdir(mmgbsa_dir):
                            if os.path.exists(os.path.join(mmgbsa_dir, 'mmgbsa_results.dat')):
                                return '✅ Completed'
                            return '🔄 Running'
                    return '⬜ Not computed'
                df_md_assays['mmgbsa'] = df_md_assays.apply(_mmgbsa_status, axis=1)

            if "show_md_assays" not in st.session_state:
                st.session_state["show_md_assays"] = False

            if st.button(f"{'Hide' if st.session_state['show_md_assays'] else 'Show'} MD Assay Registries"):
                st.session_state["show_md_assays"] = not st.session_state["show_md_assays"]

            if st.session_state["show_md_assays"]:
                _display_cols = ['assay_id', 'md_assay', 'description', 'status', 'mmgbsa',
                                 'receptor_template_name', 'docking_assay_id',
                                 'ligand_name', 'pose_id', 'assay_folder_path']
                _display_cols = [c for c in _display_cols if c in df_md_assays.columns]
                st.dataframe(
                    df_md_assays[_display_cols],
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "status": st.column_config.TextColumn("Status", width="small"),
                        "mmgbsa": st.column_config.TextColumn("MMGBSA", width="small"),
                        "receptor_template_name": st.column_config.TextColumn("PDB Template"),
                    },
                )

                _sp, _col_desc_edit = st.columns([1, 9])
                with _col_desc_edit:
                    with st.expander("✏️ Edit Assay Description", expanded=False):
                        _assay_options = df_md_assays['md_assay'].tolist()
                        _selected_assay_name = st.selectbox(
                            "Select assay to edit description:",
                            options=_assay_options,
                            key="md_assay_desc_select",
                        )
                        if _selected_assay_name:
                            _sel_row = df_md_assays[df_md_assays['md_assay'] == _selected_assay_name].iloc[0]
                            _current_desc = _sel_row.get('description') or ''
                            _new_desc = st.text_area(
                                "Description:",
                                value=_current_desc,
                                key=f"md_assay_desc_text_{_selected_assay_name}",
                            )
                            if st.button("💾 Save description", key="btn_save_md_assay_desc"):
                                try:
                                    _conn = sqlite3.connect(md_registers_db_path)
                                    _conn.execute(
                                        "UPDATE md_assays SET description = ? WHERE md_assay = ?",
                                        (_new_desc, _selected_assay_name),
                                    )
                                    _conn.commit()
                                    _conn.close()
                                    st.success(f"Description updated for '{_selected_assay_name}'.")
                                    st.rerun()
                                except Exception as _e:
                                    st.error(f"Failed to update description: {_e}")

elif page == "MolDyn analysis":
    st.title("MolDyn analysis")

    if "active_project_path" not in st.session_state:
        st.warning("No active project. Please select a project first.")
    else:
        project_path = st.session_state["active_project_path"]
        md_registers_db_path = os.path.join(project_path, "dynamics", "md_registers", "md_registers.db")

        try:
            df_all_assays = st_funcs.get_md_assay_registers(md_registers_db_path)
        except Exception:
            df_all_assays = None

        if df_all_assays is None or df_all_assays.empty:
            st.markdown(
                f"<span style='font-size: 24px; color: orange;'>No MD assays registered for project: </span>"
                f"<span style='font-size: 30px; color: green; font-weight: bold;'><b>{st.session_state.get('selected_project', 'Unknown')}</b></span>",
                unsafe_allow_html=True
            )
        else:
            df_completed = df_all_assays[df_all_assays.get('status', pd.Series(dtype=str)) == '✅ Completed']
            if 'status' in df_all_assays.columns:
                df_completed = df_all_assays[df_all_assays['status'] == '✅ Completed']
            else:
                df_completed = df_all_assays

            if df_completed.empty:
                st.info("No completed MD assays found for the active project.")
            else:
                _completed_names = df_completed['md_assay'].tolist()
                _selected_md_assay = st.selectbox(
                    "Select completed MD assay:",
                    options=_completed_names,
                    key="moldyn_analysis_assay_select",
                )

                if _selected_md_assay:
                    _assay_row = df_completed[df_completed['md_assay'] == _selected_md_assay].iloc[0]
                    _assay_folder = _assay_row.get('assay_folder_path', '')
                    _ligand_name = _assay_row.get('ligand_name', '')
                    _pose_id = _assay_row.get('pose_id', '')
                    _description = _assay_row.get('description', '') or ''

                    if _description:
                        st.markdown(f"**Description:** {_description}")
                    st.markdown(f"**Assay folder:** `{_assay_folder}`")
                    if _ligand_name:
                        st.markdown(f"**Ligand:** {_ligand_name} &nbsp; **Pose:** {_pose_id}", unsafe_allow_html=True)

                    st.divider()
                    st.subheader("Analysis")

                    _sp, _col_btns = st.columns([1, 9])
                    with _col_btns:
                        with st.expander("📊 QC on trajectory", expanded=False):
                            _prod_out = os.path.join(_assay_folder, "production.out")
                            if not os.path.exists(_prod_out):
                                st.error(f"production.out not found: {_prod_out}")
                            else:
                                import re as _re
                                _times, _etots, _ektots, _eptots, _temps = [], [], [], [], []
                                try:
                                    with open(_prod_out, "r") as _fh:
                                        _lines = _fh.readlines()
                                    for _i, _line in enumerate(_lines):
                                        if "A V E R A G E S" in _line and "S T E P S" in _line:
                                            break
                                        if "NSTEP" in _line and "TIME(PS)" in _line:
                                            _m_time = _re.search(r"TIME\(PS\)\s*=\s*([\-\d.]+)", _line)
                                            _m_temp = _re.search(r"TEMP\(K\)\s*=\s*([\-\d.]+)", _line)
                                            if _m_time and _m_temp and _i + 1 < len(_lines):
                                                _eline = _lines[_i + 1]
                                                _m_etot = _re.search(r"Etot\s*=\s*([\-\d.]+)", _eline)
                                                _m_ektot = _re.search(r"EKtot\s*=\s*([\-\d.]+)", _eline)
                                                _m_ep = _re.search(r"EPtot\s*=\s*([\-\d.]+)", _eline)
                                                if _m_etot and _m_ektot and _m_ep:
                                                    _times.append(float(_m_time.group(1)))
                                                    _temps.append(float(_m_temp.group(1)))
                                                    _etots.append(float(_m_etot.group(1)))
                                                    _ektots.append(float(_m_ektot.group(1)))
                                                    _eptots.append(float(_m_ep.group(1)))
                                except Exception as _e:
                                    st.error(f"Error reading production.out: {_e}")

                                if _times:
                                    import matplotlib.pyplot as _plt

                                    def _mean(v): return sum(v) / len(v)

                                    _plots = [
                                        (_eptots, "EPtot (kcal/mol)", "steelblue"),
                                        (_etots,  "Etot (kcal/mol)",  "seagreen"),
                                        (_ektots, "EKtot (kcal/mol)", "darkorange"),
                                        (_temps,  "TEMP (K)",         "mediumpurple"),
                                    ]
                                    _fig, _axes = _plt.subplots(2, 2, figsize=(14, 8))
                                    for _ax, (_vals, _ylabel, _color) in zip(_axes.flat, _plots):
                                        _ax.plot(_times, _vals, linewidth=0.8, color=_color)
                                        _ax.axhline(_mean(_vals), color="tomato", linestyle="--",
                                                    linewidth=1, label=f"Mean: {_mean(_vals):.2f}")
                                        _ax.set_xlabel("Time (ps)")
                                        _ax.set_ylabel(_ylabel)
                                        _ax.legend(fontsize=8)
                                    _fig.suptitle(f"Trajectory QC — {_selected_md_assay}", fontsize=12)
                                    _plt.tight_layout()
                                    st.pyplot(_fig)
                                    _plt.close(_fig)
                                    st.caption(f"Frames parsed: {len(_times)}")
                                else:
                                    st.warning("No energy data found in production.out.")

                        with st.expander("📐 Structural analysis (RMSD)", expanded=False):
                            _prmtop = os.path.join(_assay_folder, "complex.prmtop")
                            _traj   = os.path.join(_assay_folder, "production.nc")
                            _inpcrd = os.path.join(_assay_folder, "complex.inpcrd")
                            _cpptraj_bin = "/home/fredy/Programas/Amber26/ambertools26/bin/cpptraj"

                            for _req, _label in [(_prmtop, "complex.prmtop"), (_traj, "production.nc"), (_inpcrd, "complex.inpcrd")]:
                                if not os.path.exists(_req):
                                    st.error(f"Required file not found: {_label}")
                                    st.stop()

                            _mask = st.text_input(
                                "cpptraj atom mask for RMSD (e.g. @CA, :1-100@CA, :MOL):",
                                value=":1-100",
                                key="rmsd_mask_input",
                            )

                            # Check mask atom/residue count (cached by assay + mask)
                            _mask_cache_key = f"rmsd_mask_info_{_assay_folder}_{_mask}"
                            if _mask_cache_key not in st.session_state:
                                import subprocess as _sp_mod2, tempfile as _tf2
                                _check_script = (
                                    f"parmstrip !{_mask}\n"
                                    f"parminfo\n"
                                    f"quit\n"
                                )
                                try:
                                    with _tf2.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as _tmp2:
                                        _tmp2.write(_check_script)
                                        _tmp2_path = _tmp2.name
                                    _chk = _sp_mod2.run(
                                        [_cpptraj_bin, _prmtop, "-i", _tmp2_path],
                                        capture_output=True, text=True,
                                    )
                                    import os as _os2; _os2.unlink(_tmp2_path)
                                    _out = _chk.stdout + _chk.stderr
                                    import re as _re2
                                    _m_atoms = _re2.search(r"Stripped parm:.*?(\d+)\s+atoms,\s*(\d+)\s+res", _out)
                                    if _m_atoms:
                                        st.session_state[_mask_cache_key] = (int(_m_atoms.group(1)), int(_m_atoms.group(2)))
                                    else:
                                        st.session_state[_mask_cache_key] = None
                                except Exception:
                                    st.session_state[_mask_cache_key] = None

                            _mask_info = st.session_state.get(_mask_cache_key)
                            if _mask_info:
                                st.caption(f"Mask matches **{_mask_info[0]} atoms** across **{_mask_info[1]} residues**")
                            elif _mask_info is None and _mask:
                                st.warning("Mask matched no atoms — check the syntax.")

                            def _toggle_rmsd_visibility():
                                st.session_state["rmsd_plot_visible"] = not st.session_state.get("rmsd_plot_visible", True)

                            _rmsd_btn_col, _rmsd_hide_col = st.columns([2, 2])
                            with _rmsd_btn_col:
                                _run_rmsd = st.button("▶️ Run RMSD", key="btn_run_rmsd")
                            with _rmsd_hide_col:
                                _rmsd_visible = st.session_state.get("rmsd_plot_visible", True)
                                _toggle_label = "🙈 Hide RMSD plot" if _rmsd_visible else "👁 Show RMSD plot"
                                st.button(_toggle_label, key="btn_toggle_rmsd", on_click=_toggle_rmsd_visibility)

                            if _run_rmsd:
                                import subprocess as _sp_mod
                                import tempfile as _tf

                                _rmsd_out = os.path.join(_assay_folder, "rmsd_analysis.dat")
                                _cpptraj_in = (
                                    f"parm {_prmtop}\n"
                                    f"trajin {_traj}\n"
                                    f"reference {_inpcrd} parm {_prmtop}\n"
                                    f"rms reference {_mask} out {_rmsd_out} mass\n"
                                    f"run\n"
                                    f"quit\n"
                                )
                                with _tf.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as _tmp:
                                    _tmp.write(_cpptraj_in)
                                    _tmp_path = _tmp.name

                                try:
                                    _result = _sp_mod.run(
                                        [_cpptraj_bin, "-i", _tmp_path],
                                        capture_output=True, text=True,
                                    )
                                    if _result.returncode != 0:
                                        st.error(f"cpptraj error:\n{_result.stderr}")
                                        st.session_state.pop("rmsd_plot_data", None)
                                    else:
                                        _frames, _rmsds = [], []
                                        with open(_rmsd_out) as _fh:
                                            for _line in _fh:
                                                _line = _line.strip()
                                                if _line.startswith("#") or not _line:
                                                    continue
                                                _parts = _line.split()
                                                if len(_parts) >= 2:
                                                    _frames.append(float(_parts[0]))
                                                    _rmsds.append(float(_parts[1]))
                                        st.session_state["rmsd_plot_data"] = {
                                            "frames": _frames, "rmsds": _rmsds,
                                            "mask": _mask, "assay": _selected_md_assay,
                                        }
                                        st.session_state["rmsd_plot_visible"] = True
                                except Exception as _e:
                                    st.error(f"Failed to run cpptraj: {_e}")
                                    st.session_state.pop("rmsd_plot_data", None)
                                finally:
                                    import os as _os
                                    _os.unlink(_tmp_path)

                            _rmsd_data = st.session_state.get("rmsd_plot_data")
                            if _rmsd_data and _rmsd_data.get("assay") == _selected_md_assay and st.session_state.get("rmsd_plot_visible", True):
                                import matplotlib.pyplot as _plt
                                _frames = _rmsd_data["frames"]
                                _rmsds  = _rmsd_data["rmsds"]
                                if _frames:
                                    _mean_rmsd = sum(_rmsds) / len(_rmsds)
                                    _fig, _ax = _plt.subplots(figsize=(10, 4))
                                    _ax.plot(_frames, _rmsds, linewidth=0.8, color="steelblue")
                                    _ax.axhline(_mean_rmsd, color="tomato", linestyle="--",
                                                linewidth=1, label=f"Mean: {_mean_rmsd:.3f} Å")
                                    _ax.set_xlabel("Frame")
                                    _ax.set_ylabel("RMSD (Å)")
                                    _ax.set_title(
                                        f"RMSD — {_selected_md_assay}  |  mask: {_rmsd_data['mask']}"
                                    )
                                    _ax.legend(fontsize=9)
                                    _plt.tight_layout()
                                    st.pyplot(_fig)
                                    _plt.close(_fig)
                                    st.caption(
                                        f"Frames: {len(_frames)}  |  "
                                        f"Min: {min(_rmsds):.3f}  |  "
                                        f"Max: {max(_rmsds):.3f}  |  "
                                        f"Mean: {_mean_rmsd:.3f} Å"
                                    )

                        with st.expander("📈 Structural analysis (RMSF)", expanded=False):
                            _prmtop = os.path.join(_assay_folder, "complex.prmtop")
                            _traj   = os.path.join(_assay_folder, "production.nc")
                            _cpptraj_bin = "/home/fredy/Programas/Amber26/ambertools26/bin/cpptraj"

                            for _req, _label in [(_prmtop, "complex.prmtop"), (_traj, "production.nc")]:
                                if not os.path.exists(_req):
                                    st.error(f"Required file not found: {_label}")
                                    st.stop()

                            _rmsf_mask = st.text_input(
                                "cpptraj atom mask for RMSF (e.g. :1-100, !:WAT):",
                                value=":1-100",
                                key="rmsf_mask_input",
                            )

                            # Mask atom/residue count (reuse same cache pattern as RMSD)
                            _rmsf_cache_key = f"rmsf_mask_info_{_assay_folder}_{_rmsf_mask}"
                            if _rmsf_cache_key not in st.session_state:
                                import subprocess as _sp_rmsf, tempfile as _tf_rmsf
                                _rmsf_check = f"parmstrip !{_rmsf_mask}\nparminfo\nquit\n"
                                try:
                                    with _tf_rmsf.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as _tmp_r:
                                        _tmp_r.write(_rmsf_check)
                                        _tmp_r_path = _tmp_r.name
                                    _chk_r = _sp_rmsf.run(
                                        [_cpptraj_bin, _prmtop, "-i", _tmp_r_path],
                                        capture_output=True, text=True,
                                    )
                                    import os as _os_r; _os_r.unlink(_tmp_r_path)
                                    import re as _re_r
                                    _m_r = _re_r.search(r"Stripped parm:.*?(\d+)\s+atoms,\s*(\d+)\s+res", _chk_r.stdout + _chk_r.stderr)
                                    st.session_state[_rmsf_cache_key] = (int(_m_r.group(1)), int(_m_r.group(2))) if _m_r else None
                                except Exception:
                                    st.session_state[_rmsf_cache_key] = None

                            _rmsf_info = st.session_state.get(_rmsf_cache_key)
                            if _rmsf_info:
                                st.caption(f"Mask matches **{_rmsf_info[0]} atoms** across **{_rmsf_info[1]} residues**")
                            elif _rmsf_info is None and _rmsf_mask:
                                st.warning("Mask matched no atoms — check the syntax.")

                            def _toggle_rmsf_visibility():
                                st.session_state["rmsf_plot_visible"] = not st.session_state.get("rmsf_plot_visible", True)

                            _rmsf_btn_col, _rmsf_hide_col = st.columns([2, 2])
                            with _rmsf_btn_col:
                                _run_rmsf = st.button("▶️ Run RMSF", key="btn_run_rmsf")
                            with _rmsf_hide_col:
                                _rmsf_visible = st.session_state.get("rmsf_plot_visible", True)
                                st.button(
                                    "🙈 Hide RMSF plot" if _rmsf_visible else "👁 Show RMSF plot",
                                    key="btn_toggle_rmsf",
                                    on_click=_toggle_rmsf_visibility,
                                )

                            if _run_rmsf:
                                import subprocess as _sp_rmsf2, tempfile as _tf_rmsf2
                                _rmsf_out = os.path.join(_assay_folder, "rmsf_analysis.dat")
                                _rmsf_script = (
                                    f"parm {_prmtop}\n"
                                    f"trajin {_traj}\n"
                                    f"atomicfluct out {_rmsf_out} {_rmsf_mask} byres\n"
                                    f"run\n"
                                    f"quit\n"
                                )
                                with _tf_rmsf2.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as _tmp_r2:
                                    _tmp_r2.write(_rmsf_script)
                                    _tmp_r2_path = _tmp_r2.name
                                try:
                                    _res_r = _sp_rmsf2.run(
                                        [_cpptraj_bin, "-i", _tmp_r2_path],
                                        capture_output=True, text=True,
                                    )
                                    if _res_r.returncode != 0:
                                        st.error(f"cpptraj error:\n{_res_r.stderr}")
                                        st.session_state.pop("rmsf_plot_data", None)
                                    else:
                                        _res_nums, _rmsf_vals = [], []
                                        with open(_rmsf_out) as _fh_r:
                                            for _line in _fh_r:
                                                _line = _line.strip()
                                                if _line.startswith("#") or not _line:
                                                    continue
                                                _parts = _line.split()
                                                if len(_parts) >= 2:
                                                    _res_nums.append(float(_parts[0]))
                                                    _rmsf_vals.append(float(_parts[1]))
                                        st.session_state["rmsf_plot_data"] = {
                                            "residues": _res_nums, "rmsf": _rmsf_vals,
                                            "mask": _rmsf_mask, "assay": _selected_md_assay,
                                        }
                                        st.session_state["rmsf_plot_visible"] = True
                                except Exception as _e:
                                    st.error(f"Failed to run cpptraj: {_e}")
                                    st.session_state.pop("rmsf_plot_data", None)
                                finally:
                                    import os as _os_r2; _os_r2.unlink(_tmp_r2_path)

                            _rmsf_data = st.session_state.get("rmsf_plot_data")
                            if _rmsf_data and _rmsf_data.get("assay") == _selected_md_assay and st.session_state.get("rmsf_plot_visible", True):
                                import matplotlib.pyplot as _plt
                                _res_nums = _rmsf_data["residues"]
                                _rmsf_vals = _rmsf_data["rmsf"]
                                if _res_nums:
                                    _fig, _ax = _plt.subplots(figsize=(10, 4))
                                    _ax.plot(_res_nums, _rmsf_vals, linewidth=0.8, color="seagreen")
                                    _mean_rmsf = sum(_rmsf_vals) / len(_rmsf_vals)
                                    _ax.axhline(_mean_rmsf, color="tomato", linestyle="--",
                                                linewidth=1, label=f"Mean: {_mean_rmsf:.3f} Å")
                                    _ax.set_xlabel("Residue")
                                    _ax.set_ylabel("RMSF (Å)")
                                    _ax.set_title(
                                        f"RMSF — {_selected_md_assay}  |  mask: {_rmsf_data['mask']}"
                                    )
                                    _ax.legend(fontsize=9)
                                    _plt.tight_layout()
                                    st.pyplot(_fig)
                                    _plt.close(_fig)
                                    st.caption(
                                        f"Residues: {len(_res_nums)}  |  "
                                        f"Min: {min(_rmsf_vals):.3f}  |  "
                                        f"Max: {max(_rmsf_vals):.3f}  |  "
                                        f"Mean: {_mean_rmsf:.3f} Å"
                                    )

                        with st.expander("🔬 Display MMGBSA results", expanded=False):
                            import json as _json, re as _re_mmgbsa

                            # --- source 1: DB JSON blob ---
                            _mmgbsa_json = _assay_row.get('mmgbsa_results')
                            _mmgbsa_results = None
                            _mmgbsa_params  = None
                            if _mmgbsa_json:
                                try:
                                    _payload = _json.loads(_mmgbsa_json)
                                    _mmgbsa_results = _payload.get('results', {})
                                    _mmgbsa_params  = _payload.get('parameters', {})
                                except Exception:
                                    pass

                            # --- source 2: mmgbsa_results.dat file ---
                            if not _mmgbsa_results:
                                _dat_path = os.path.join(_assay_folder, 'mmgbsa', 'mmgbsa_results.dat')
                                if os.path.exists(_dat_path):
                                    try:
                                        with open(_dat_path) as _fh_m:
                                            _dat_content = _fh_m.read()
                                        _comp_pat = _re_mmgbsa.compile(
                                            r'^\s*(VDWAALS|EEL|EGB|ESURF|EPOL|ENPOLAR|EDISPER|EPB|ENPB|ECAVITY'
                                            r'|DELTA G gas|DELTA G solv|DELTA TOTAL)\s+(-?\d+\.\d+)\s+(\d+\.\d+)',
                                            _re_mmgbsa.MULTILINE,
                                        )
                                        _mmgbsa_results = {}
                                        for _m in _comp_pat.finditer(_dat_content):
                                            _mmgbsa_results[_m.group(1).strip()] = {
                                                'average': float(_m.group(2)),
                                                'std_dev': float(_m.group(3)),
                                            }
                                        _bind = _re_mmgbsa.search(
                                            r'DELTA G binding\s*=\s*(-?\d+\.\d+)\s+\+/-\s+(\d+\.\d+)',
                                            _dat_content,
                                        )
                                        if _bind:
                                            _mmgbsa_results['DELTA_G_binding'] = {
                                                'average': float(_bind.group(1)),
                                                'std_dev': float(_bind.group(2)),
                                            }
                                        elif 'DELTA TOTAL' in _mmgbsa_results:
                                            _mmgbsa_results['DELTA_G_binding'] = _mmgbsa_results['DELTA TOTAL']
                                        if not _mmgbsa_results:
                                            _mmgbsa_results = None
                                    except Exception:
                                        _mmgbsa_results = None

                            if not _mmgbsa_results:
                                st.info("No MMGBSA results found for this assay.")
                            else:
                                # Parameters summary
                                if _mmgbsa_params:
                                    st.caption(
                                        f"GB model: igb={_mmgbsa_params.get('igb', '?')}  |  "
                                        f"Salt: {_mmgbsa_params.get('saltcon', '?')} M  |  "
                                        f"Ligand mask: {_mmgbsa_params.get('ligand_mask', '?')}  |  "
                                        f"Frames: {_mmgbsa_params.get('startframe', '?')}–"
                                        f"{_mmgbsa_params.get('endframe', '?')} "
                                        f"(every {_mmgbsa_params.get('interval', '?')})"
                                    )

                                # Energy components table
                                _component_order = [
                                    'VDWAALS', 'EEL', 'EGB', 'ESURF', 'EPOL',
                                    'DELTA G gas', 'DELTA G solv', 'DELTA TOTAL',
                                ]
                                _rows = []
                                for _comp in _component_order:
                                    if _comp in _mmgbsa_results:
                                        _v = _mmgbsa_results[_comp]
                                        _rows.append({
                                            'Component': _comp,
                                            'Average (kcal/mol)': f"{_v['average']:.3f}",
                                            'Std Dev':            f"{_v['std_dev']:.3f}",
                                        })
                                if _rows:
                                    st.dataframe(
                                        pd.DataFrame(_rows),
                                        use_container_width=True,
                                        hide_index=True,
                                    )

                                # Binding free energy highlight
                                _dg = _mmgbsa_results.get('DELTA_G_binding')
                                if _dg:
                                    st.markdown(
                                        f"**ΔG binding = {_dg['average']:.3f} ± {_dg['std_dev']:.3f} kcal/mol**"
                                    )

                        with st.expander("⚙️ Execute custom analysis", expanded=False):
                            _scripts_db = os.path.join(
                                project_path, "dynamics", "md_registers", "md_analysis_scripts.db"
                            )
                            _cpptraj_bin = "/home/fredy/Programas/Amber26/ambertools26/bin/cpptraj"
                            _prmtop = os.path.join(_assay_folder, "complex.prmtop")
                            _traj   = os.path.join(_assay_folder, "production.nc")
                            _inpcrd = os.path.join(_assay_folder, "complex.inpcrd")

                            # Ensure scripts table exists
                            _conn_s = sqlite3.connect(_scripts_db)
                            _conn_s.execute("""
                                CREATE TABLE IF NOT EXISTS analysis_scripts (
                                    script_id   INTEGER PRIMARY KEY AUTOINCREMENT,
                                    script_name TEXT UNIQUE NOT NULL,
                                    script_type TEXT NOT NULL DEFAULT 'cpptraj',
                                    script_code TEXT NOT NULL,
                                    description TEXT,
                                    created_date TEXT
                                )
                            """)
                            _conn_s.commit()
                            _existing = pd.read_sql_query(
                                "SELECT script_id, script_name, script_type, script_code, description FROM analysis_scripts ORDER BY script_name",
                                _conn_s,
                            )
                            _conn_s.close()

                            # Script selector
                            _script_options = ["-- New script --"] + _existing['script_name'].tolist()
                            _sel_script_name = st.selectbox(
                                "Load existing script:",
                                options=_script_options,
                                key="custom_script_select",
                            )

                            # Pre-fill fields when an existing script is selected
                            if _sel_script_name != "-- New script --" and not _existing.empty:
                                _sel_row_s = _existing[_existing['script_name'] == _sel_script_name].iloc[0]
                                _default_name  = _sel_row_s['script_name']
                                _default_type  = _sel_row_s['script_type']
                                _default_code  = _sel_row_s['script_code']
                                _default_desc  = _sel_row_s['description'] or ''
                            else:
                                _default_name  = ''
                                _default_type  = 'cpptraj'
                                _default_code  = ''
                                _default_desc  = ''

                            # When selection changes to an existing script, sync type and code
                            _prev_select_key = "custom_script_prev_select"
                            if (
                                _sel_script_name != "-- New script --"
                                and st.session_state.get(_prev_select_key) != _sel_script_name
                            ):
                                st.session_state["custom_script_name"] = _default_name
                                st.session_state["custom_script_type"] = _default_type
                                st.session_state["custom_script_code"] = _default_code
                            st.session_state[_prev_select_key] = _sel_script_name

                            _script_name_in = st.text_input(
                                "Script name:", value=_default_name, key="custom_script_name"
                            )
                            _script_type_in = st.radio(
                                "Script type:", options=["cpptraj", "python"],
                                index=0 if _default_type == "cpptraj" else 1,
                                horizontal=True, key="custom_script_type",
                            )

                            _cpptraj_hint = (
                                f"# Available variables (auto-substituted):\n"
                                f"# {{prmtop}}       = {_prmtop}\n"
                                f"# {{traj}}         = {_traj}\n"
                                f"# {{inpcrd}}       = {_inpcrd}\n"
                                f"# {{assay_folder}} = {_assay_folder}\n\n"
                                f"parm {{prmtop}}\n"
                                f"trajin {{traj}}\n"
                                f"# ... add cpptraj commands ...\n"
                                f"run\n"
                                f"quit\n"
                            )
                            _python_hint = (
                                "# Auto-substituted variables (read-only, pre-bound):\n"
                                f"# prmtop      = '{_prmtop}'\n"
                                f"# traj        = '{_traj}'\n"
                                f"# inpcrd      = '{_inpcrd}'\n"
                                f"# assay_folder = '{_assay_folder}'\n"
                                f"# assay_name  = '{_selected_md_assay}'\n"
                                "# st, pd, os  — Streamlit, pandas, os modules\n\n"
                                "import MDAnalysis as mda\n"
                                "import matplotlib.pyplot as plt\n"
                                "from MDAnalysis.analysis import rms\n\n"
                                "# Load topology and trajectory\n"
                                "u = mda.Universe(prmtop, traj)\n\n"
                                "# Select atoms\n"
                                "sel = u.select_atoms('protein and name CA')\n\n"
                                "# Compute RMSD relative to first frame\n"
                                "R = rms.RMSD(sel, select='protein and name CA')\n"
                                "R.run()\n\n"
                                "# Plot\n"
                                "fig, ax = plt.subplots(figsize=(10, 4))\n"
                                "ax.plot(R.results.rmsd[:, 1], R.results.rmsd[:, 2],\n"
                                "        linewidth=0.8, color='steelblue')\n"
                                "ax.set_xlabel('Frame')\n"
                                "ax.set_ylabel('RMSD (Å)')\n"
                                "ax.set_title(f'RMSD — {assay_name}')\n"
                                "st.pyplot(fig)\n"
                                "plt.close(fig)\n"
                            )

                            # When the type radio changes on a new script, reset the code box
                            _prev_type_key = "custom_script_prev_type"
                            _is_new_script = (_sel_script_name == "-- New script --")
                            if _is_new_script and st.session_state.get(_prev_type_key) != _script_type_in:
                                st.session_state["custom_script_code"] = (
                                    _cpptraj_hint if _script_type_in == "cpptraj" else _python_hint
                                )
                            st.session_state[_prev_type_key] = _script_type_in

                            _hint = _cpptraj_hint if _script_type_in == "cpptraj" else _python_hint
                            _script_code_in = st.text_area(
                                "Script code:",
                                value=_default_code if _default_code else _hint,
                                height=220,
                                key="custom_script_code",
                            )
                            _script_desc_in = st.text_input(
                                "Description (optional):", value=_default_desc, key="custom_script_desc"
                            )

                            # Action buttons
                            _cs_col_save, _cs_col_del, _cs_col_run = st.columns(3)

                            with _cs_col_save:
                                if st.button("💾 Save script", key="btn_save_custom_script"):
                                    if not _script_name_in.strip():
                                        st.warning("Enter a script name before saving.")
                                    else:
                                        try:
                                            _conn_s2 = sqlite3.connect(_scripts_db)
                                            _conn_s2.execute("""
                                                INSERT INTO analysis_scripts
                                                    (script_name, script_type, script_code, description, created_date)
                                                VALUES (?, ?, ?, ?, datetime('now'))
                                                ON CONFLICT(script_name) DO UPDATE SET
                                                    script_type  = excluded.script_type,
                                                    script_code  = excluded.script_code,
                                                    description  = excluded.description,
                                                    created_date = datetime('now')
                                            """, (_script_name_in.strip(), _script_type_in,
                                                  _script_code_in, _script_desc_in.strip()))
                                            _conn_s2.commit()
                                            _conn_s2.close()
                                            st.success(f"Script '{_script_name_in.strip()}' saved.")
                                            st.rerun()
                                        except Exception as _e:
                                            st.error(f"Error saving script: {_e}")

                            with _cs_col_del:
                                _del_confirmed = st.checkbox(
                                    "Confirm deletion", key="chk_del_custom_script"
                                )
                                if st.button(
                                    "🗑️ Delete script", key="btn_del_custom_script",
                                    disabled=not _del_confirmed,
                                ):
                                    _del_name = _script_name_in.strip()
                                    if not _del_name:
                                        st.warning("Enter the script name to delete.")
                                    else:
                                        try:
                                            _conn_s3 = sqlite3.connect(_scripts_db)
                                            _cur_del = _conn_s3.execute(
                                                "DELETE FROM analysis_scripts WHERE script_name = ?",
                                                (_del_name,)
                                            )
                                            _conn_s3.commit()
                                            _conn_s3.close()
                                            if _cur_del.rowcount:
                                                st.success(f"Script '{_del_name}' deleted.")
                                            else:
                                                st.warning(f"No script named '{_del_name}' found.")
                                            st.rerun()
                                        except Exception as _e:
                                            st.error(f"Error deleting script: {_e}")

                            with _cs_col_run:
                                if st.button("▶️ Run script", key="btn_run_custom_script"):
                                    if not _script_code_in.strip():
                                        st.warning("No code to run.")
                                    else:
                                        if _script_type_in == "cpptraj":
                                            import tempfile as _tf_cs, subprocess as _sp_cs
                                            _resolved = _script_code_in.format(
                                                prmtop=_prmtop, traj=_traj,
                                                inpcrd=_inpcrd, assay_folder=_assay_folder,
                                            )
                                            with _tf_cs.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as _tmp_cs:
                                                _tmp_cs.write(_resolved)
                                                _tmp_cs_path = _tmp_cs.name
                                            try:
                                                _res_cs = _sp_cs.run(
                                                    [_cpptraj_bin, "-i", _tmp_cs_path],
                                                    capture_output=True, text=True,
                                                )
                                                _out_cs = _res_cs.stdout + ("\n" + _res_cs.stderr if _res_cs.stderr else "")
                                                st.session_state["custom_script_output"] = _out_cs
                                                st.session_state["custom_script_rc"] = _res_cs.returncode
                                            except Exception as _e:
                                                st.session_state["custom_script_output"] = str(_e)
                                                st.session_state["custom_script_rc"] = -1
                                            finally:
                                                import os as _os_cs; _os_cs.unlink(_tmp_cs_path)
                                        else:
                                            import io as _io_cs
                                            _stdout_cap = _io_cs.StringIO()
                                            _exec_ns = {
                                                "prmtop": _prmtop, "traj": _traj,
                                                "inpcrd": _inpcrd, "assay_folder": _assay_folder,
                                                "assay_name": _selected_md_assay,
                                                "st": st, "pd": pd, "os": os,
                                            }
                                            try:
                                                exec(_script_code_in, _exec_ns)
                                                st.session_state["custom_script_output"] = "✅ Python script executed."
                                                st.session_state["custom_script_rc"] = 0
                                            except Exception as _e:
                                                st.session_state["custom_script_output"] = f"❌ {type(_e).__name__}: {_e}"
                                                st.session_state["custom_script_rc"] = -1

                            # Output display
                            _cs_out = st.session_state.get("custom_script_output")
                            if _cs_out:
                                _cs_rc = st.session_state.get("custom_script_rc", 0)
                                if _cs_rc != 0:
                                    st.error("Script exited with errors:")
                                st.code(_cs_out, language="bash")

elif page == "ProLIF Conditions":
    st.title("ProLIF Conditions")

    if "active_project_path" not in st.session_state:
        st.warning("No active project. Please select a project first.")
    else:
        project_path = st.session_state["active_project_path"]

        st.subheader("Registered ProLIF Conditions")

        _prolif_records = st_funcs.get_prolif_conditions_records(project_path)

        if not _prolif_records:
            st.info(
                "No ProLIF conditions registered yet. "
                "Create conditions using `MolDock.create_prolif_conditions()`."
            )
        else:
            ## Summary table: ID + Description
            st.dataframe(
                [{"ID": r["id"], "Description": r["description"]} for r in _prolif_records],
                use_container_width=True,
                hide_index=True,
            )

            st.divider()

            ## Selector
            _descriptions = [r["description"] for r in _prolif_records]
            _selected_desc = st.selectbox(
                "Select a condition to inspect:",
                _descriptions,
                key="prolif_conditions_inspect_select",
            )

            if _selected_desc:
                _selected_record = next(r for r in _prolif_records if r["description"] == _selected_desc)
                st.markdown(f"#### 📌 {_selected_desc}  *(ID {_selected_record['id']})*")
                st.json(_selected_record["conditions"], expanded=False)

                st.divider()
                _default_export_path = os.path.join(
                    project_path, 'docking', 'params', 'prolif_conditions.json'
                )
                _export_path_input = st.text_input(
                    "Export path:",
                    value=_default_export_path,
                    key=f"prolif_cond_export_path_{_selected_record['id']}",
                )
                if st.button("📤 Export ProLIF Conditions", key=f"btn_export_prolif_cond_{_selected_record['id']}"):
                    _result = st_funcs.export_prolif_conditions_to_file(
                        project_path, _selected_record['id'], _export_path_input.strip()
                    )
                    if _result == "exported":
                        st.success(f"✅ ProLIF conditions exported to: {_export_path_input.strip()}")
                    else:
                        st.error(f"❌ {_result[6:]}")


elif page == "Docking analysis":
    _analysis_col_title, _analysis_col_info = st.columns([2, 3])
    with _analysis_col_title:
        st.title("Docking analysis")
    with _analysis_col_info:
        st.markdown(
            f"<div style='padding-top: 14px;'>"
            f"<span style='color: grey;'>Project: </span>"
            f"<span style='color: green; font-weight: bold;'>{active_project or 'None'}</span>"
            f"&nbsp;&nbsp;|&nbsp;&nbsp;"
            f"<span style='color: grey;'>Assay: </span>"
            f"<span style='color: orange; font-weight: bold;'>{active_assay or 'None'}</span>"
            f"</div>",
            unsafe_allow_html=True
        )
    st.write("Welcome to the Analysis page.")
            
    if "selected_assay_name" not in st.session_state:
        st.warning("Select a docking assay")
    else:
        results_db_path = os.path.join(st.session_state["active_project_path"], "docking", "docking_assays", st.session_state["selected_assay_name"], "results", f"{st.session_state['selected_assay_name']}.db")

        df_results = st_funcs.get_docking_results(results_db_path)
        df_results = st_funcs.add_ligand_name_column(
            df_results, st.session_state["active_project_path"], st.session_state["selected_assay_name"]
        )

        df_mmpbsa_poses_results = st_funcs.get_mmpbsa_results(results_db_path)

        if df_results is None or df_results.empty:
            st.warning("No docking results")
        else:
            with st.expander(f"📊 {st.session_state['selected_assay_name']} Docking Results", expanded=False):
                st.write(df_results)

                if df_mmpbsa_poses_results is not None and not df_mmpbsa_poses_results.empty:
                    _sp, _col = st.columns([1, 9])
                    with _col:
                        if st.button("Update MMGBSA data", key="btn_update_mmgbsa_results"):
                            n_updated, errs = st_funcs.populate_results_with_mmgbsa_energies(results_db_path)
                            if errs:
                                for msg in errs:
                                    st.warning(msg)
                            if n_updated > 0:
                                st.success(f"Updated mmgbsa_total_energy and mmgbsa_gas_energy for {n_updated} pose(s) in the Results table.")
                            else:
                                st.info("No poses were updated.")

            ## Single "View Poses" collapsible menu — always shown; extraction controls
            ## are nested inside when no poses have been extracted yet.
            extracted_poses = st_funcs.get_extracted_poses_info(results_db_path)
            any_active = any(e["active"] for e in extracted_poses)

            with st.expander("🔍 View Poses", expanded=False):
                if not any_active:
                    _sp, _ep_col = st.columns([1, 9])
                    with _ep_col:
                        _engine = st_funcs.get_assay_engine(
                            st.session_state["active_project_path"],
                            st.session_state["selected_assay_name"]
                        )
                        if _engine == "AutoDockGPU":
                            _criteria_labels = [
                                "1 - Most stable poses",
                                "2 - Most populated poses",
                                "3 - Most stable + most populated",
                                "4 - All poses",
                            ]
                        else:
                            _criteria_labels = [
                                "1 - Most stable poses",
                                "2 - All poses",
                            ]
                        _ep_inner_left, _ep_inner_mid, _ep_inner_right = st.columns([3, 3, 2])
                        with _ep_inner_left:
                            _criteria = st.selectbox(
                                "Extraction criteria:",
                                _criteria_labels,
                                key="extract_poses_criteria",
                            )

                        ## "Most stable" (1) and "Most stable + most populated" (3) rank by a
                        ## numeric score column; offer a choice when more than one is available,
                        ## mirroring MolDock._prompt_score_column() from the interactive CLI.
                        _needs_score_col = _criteria.startswith("1") or _criteria.startswith("3")
                        _score_col_labels = {
                            'docking_score': 'Docking score',
                            'mmgbsa_total_energy': 'MMGBSA total energy',
                            'mmgbsa_gas_energy': 'MMGBSA gas energy',
                        }
                        _selected_score_col = None
                        with _ep_inner_mid:
                            if _needs_score_col:
                                _score_col_candidates = st_funcs.get_available_score_columns(results_db_path)
                                if len(_score_col_candidates) > 1:
                                    _selected_score_col = st.selectbox(
                                        "Scoring column:",
                                        _score_col_candidates,
                                        format_func=lambda c: _score_col_labels.get(c, c),
                                        key="extract_poses_score_column",
                                    )

                        with _ep_inner_right:
                            st.write("")
                            st.write("")
                            if st.button("Extract Poses", key="btn_extract_poses"):
                                _selection = _criteria[0]  # leading digit "1", "2", "3", or "4"
                                with st.spinner("Extracting poses..."):
                                    _ok, _msg = st_funcs.extract_poses_for_assay(
                                        st.session_state["selected_project"],
                                        st.session_state["active_project_path"],
                                        st.session_state["selected_assay_name"],
                                        _selection,
                                        score_column=_selected_score_col,
                                    )
                                if _ok:
                                    st.success(_msg)
                                    st.rerun()
                                else:
                                    st.error(_msg)
                else:
                    ## Reference PDB uploader (shared across all pose folders)
                    _sp, _col = st.columns([1, 9])
                    with _col:
                        ref_file = st.file_uploader(
                            "Load a reference PDB for superimposition (optional):",
                            type=["pdb"],
                            key="reference_pdb_uploader"
                        )
                        if ref_file is not None:
                            st.session_state["reference_pdb_data"] = ref_file.read().decode("utf-8")
                            st.success(f"Reference loaded: {ref_file.name}")
                        elif "reference_pdb_data" not in st.session_state:
                            st.session_state["reference_pdb_data"] = None

                    ## Delete all extracted poses (PDB files + folders)
                    _sp, _col = st.columns([1, 9])
                    with _col:
                        _del_poses_key = "confirm_delete_all_poses"
                        if not st.session_state.get(_del_poses_key):
                            if st.button("🗑️ Delete all poses", key="btn_delete_all_poses"):
                                st.session_state[_del_poses_key] = True
                                st.rerun()
                        else:
                            st.warning("Delete all extracted poses and their folders? This cannot be undone.")
                            _dpc1, _dpc2 = st.columns(2)
                            with _dpc1:
                                if st.button("Yes, delete all poses", key="btn_confirm_delete_all_poses"):
                                    for entry in extracted_poses:
                                        if os.path.isdir(entry["path"]):
                                            shutil.rmtree(entry["path"])
                                    st.session_state[_del_poses_key] = False
                                    st.success("All extracted poses deleted.")
                                    st.rerun()
                            with _dpc2:
                                if st.button("Cancel", key="btn_cancel_delete_all_poses"):
                                    st.session_state[_del_poses_key] = False
                                    st.rerun()

                    for entry in extracted_poses:
                        label = f"{entry['directory']} ({entry['count']} PDB files)" if entry["active"] else entry["directory"]
                        key = f"btn_poses_{entry['directory']}"
                        if entry["active"]:
                            if key not in st.session_state:
                                st.session_state[key] = False
                            _sp, _col = st.columns([1, 9])
                            with _col:
                                if st.button(f"{'Hide' if st.session_state[key] else 'Show'} {label}", key=key + "_btn"):
                                    st.session_state[key] = not st.session_state[key]
                                    st.rerun()
                            if st.session_state[key]:
                                st.divider()
                                filename_to_pose_id = st_funcs.get_filename_to_pose_id_map(results_db_path)
                                pdb_files = sorted(
                                    glob.glob(os.path.join(entry["path"], "*.pdb")),
                                    key=lambda f: filename_to_pose_id.get(os.path.basename(f), {}).get("pose_id", float("inf"))
                                )
                                pdb_names = [os.path.basename(f) for f in pdb_files]
                                idx_key = f"pose_idx_{entry['directory']}"
                                if idx_key not in st.session_state:
                                    st.session_state[idx_key] = 0
                                ## Clamp index to valid range (e.g. folder contents changed)
                                st.session_state[idx_key] = max(0, min(st.session_state[idx_key], len(pdb_names) - 1))

                                ## Selectbox with NO key — driven entirely by idx_key via index param
                                _sp, _col = st.columns([2, 8])
                                with _col:
                                    selected_file = st.selectbox(
                                        f"Select a pose file from {entry['directory']}:",
                                        pdb_names,
                                        index=st.session_state[idx_key],
                                    )
                                ## Keep index in sync when user picks from the selectbox directly
                                if selected_file:
                                    st.session_state[idx_key] = pdb_names.index(selected_file)

                                if selected_file:
                                    full_path = os.path.join(entry["path"], selected_file)
                                    with open(full_path, "r") as f:
                                        pdb_data = f.read()
                                    view = py3Dmol.view(width=800, height=500)
                                    ## Add selected pose (model 0) — stick + spectrum cartoon
                                    view.addModel(pdb_data, "pdb")
                                    view.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                                    ## Add reference structure (model 1) if loaded — CPK spheres
                                    if st.session_state.get("reference_pdb_data"):
                                        view.addModel(st.session_state["reference_pdb_data"], "pdb")
                                        view.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                                    view.zoomTo()

                                    ## Navigation prev/next arrow buttons flanking the viewer
                                    current_idx = st.session_state[idx_key]
                                    has_ref = st.session_state.get("reference_pdb_data")
                                    prev_col, viewer_col, next_col = st.columns([0.5, 8, 0.5 if not has_ref else 1])
                                    with prev_col:
                                        st.write("")  ## vertical alignment spacer
                                        st.write("")
                                        st.write("")
                                        if st.button("◀", key=f"prev_pose_{entry['directory']}", disabled=(current_idx == 0)):
                                            st.session_state[idx_key] = current_idx - 1
                                            st.rerun()
                                    with viewer_col:
                                        st.components.v1.html(view.write_html(), height=520)
                                    with next_col:
                                        st.write("")
                                        st.write("")
                                        st.write("")
                                        if st.button("▶", key=f"next_pose_{entry['directory']}", disabled=(current_idx >= len(pdb_names) - 1)):
                                            st.session_state[idx_key] = current_idx + 1
                                            st.rerun()
                                        if has_ref:
                                            st.write("")
                                            flag_pos_key = f"flag_pos_{entry['directory']}_{selected_file}"
                                            if st.button("✅ Flag positive binding pose", key=flag_pos_key):
                                                result = st_funcs.save_positive_binder(
                                                    project_path=st.session_state["active_project_path"],
                                                    assay_name=st.session_state.get("selected_assay_name", "unknown"),
                                                    pose_file=selected_file,
                                                    directory=entry["directory"],
                                                    pose_full_path=full_path,
                                                )
                                                if result == "saved":
                                                    st.success(f"Saved: {selected_file}")
                                                elif result == "duplicate":
                                                    st.warning("Already flagged as positive.")
                                                elif result == "conflict":
                                                    st.error("Already flagged as negative binder.")
                                                else:
                                                    st.error(result)
                                            flag_neg_key = f"flag_neg_{entry['directory']}_{selected_file}"
                                            if st.button("❌ Flag negative binding pose", key=flag_neg_key):
                                                result = st_funcs.save_negative_binder(
                                                    project_path=st.session_state["active_project_path"],
                                                    assay_name=st.session_state.get("selected_assay_name", "unknown"),
                                                    pose_file=selected_file,
                                                    directory=entry["directory"],
                                                    pose_full_path=full_path,
                                                )
                                                if result == "saved":
                                                    st.success(f"Saved: {selected_file}")
                                                elif result == "duplicate":
                                                    st.warning("Already flagged as negative.")
                                                elif result == "conflict":
                                                    st.error("Already flagged as positive binder.")
                                                else:
                                                    st.error(result)

                                    ## pose counter label
                                    pose_info = filename_to_pose_id.get(selected_file, {})
                                    pose_id_label = pose_info.get("pose_id", "?")
                                    lig_label = pose_info["ligname"] if pose_info else "?"
                                    st.caption(f"Pose ID: {pose_id_label}  |  {lig_label}  ({current_idx + 1} of {len(pdb_names)})")

                                    ## VMD script creation
                                    with st.expander("🎬 Create VMD Script", expanded=False):
                                        _default_vmd_path = os.path.join(
                                            os.path.expanduser("~"), "Desktop",
                                            os.path.splitext(selected_file)[0] + "_vmd.tcl"
                                        )
                                        _vmd_path = st.text_input(
                                            "Save script to:",
                                            value=_default_vmd_path,
                                            key=f"vmd_path_{entry['directory']}_{selected_file}",
                                        )
                                        if st.button("💾 Save VMD Script", key=f"btn_save_vmd_{entry['directory']}_{selected_file}"):
                                            try:
                                                _ref_pdb_path = None
                                                if st.session_state.get("reference_pdb_data"):
                                                    _ref_pdb_path = os.path.join(
                                                        os.path.dirname(os.path.abspath(_vmd_path.strip())),
                                                        "reference.pdb"
                                                    )
                                                    with open(_ref_pdb_path, "w") as _rf:
                                                        _rf.write(st.session_state["reference_pdb_data"])
                                                _script = st_funcs.generate_vmd_script(
                                                    os.path.abspath(full_path), _ref_pdb_path
                                                )
                                                _vmd_out = _vmd_path.strip()
                                                os.makedirs(os.path.dirname(os.path.abspath(_vmd_out)), exist_ok=True)
                                                with open(_vmd_out, "w") as _sf:
                                                    _sf.write(_script)
                                                st.success(f"VMD script saved to: {_vmd_out}")
                                                if _ref_pdb_path:
                                                    st.info(f"Reference PDB saved alongside: {_ref_pdb_path}")
                                            except Exception as _e:
                                                st.error(f"Could not save VMD script: {_e}")

                                ## ProLIF Fingerprints for this pose directory — always follows
                                ## the pose currently selected above (no separate pose picker),
                                ## so the FPS shown can never drift from the extracted pose set,
                                ## regardless of which score column was used to extract it.
                                dir_name = entry['directory']
                                dir_pose_ids = [
                                    filename_to_pose_id[f]["pose_id"]
                                    for f in pdb_names if f in filename_to_pose_id
                                ]
                                prolif_tables_dir = st_funcs.get_prolif_tables_by_pose_ids(results_db_path, dir_pose_ids)
                                prolif_dir_key = f"show_prolif_{dir_name}"
                                if not prolif_tables_dir:
                                    _sp, _col = st.columns([2, 8])
                                    with _col:
                                        st.button("Show ProLIF Fingerprints", disabled=True, key=f"btn_prolif_disabled_{dir_name}")
                                else:
                                    if prolif_dir_key not in st.session_state:
                                        st.session_state[prolif_dir_key] = False
                                    _dir_label = {
                                        "most_stable_poses": "most_stable",
                                        "most_populated_poses": "most_populated",
                                        "most_populated_and_stable_poses": "most_populated_and_stable",
                                        "all_poses": "all_poses",
                                    }.get(dir_name, dir_name)
                                    _sp, _col = st.columns([2, 8])
                                    with _col:
                                        if st.button(
                                            f"{'Hide' if st.session_state[prolif_dir_key] else 'Show'} ProLIF FPS - {_dir_label}",
                                            key=f"btn_prolif_{dir_name}"
                                        ):
                                            st.session_state[prolif_dir_key] = not st.session_state[prolif_dir_key]
                                            st.rerun()
                                    if st.session_state[prolif_dir_key]:
                                        selected_table_dir = st.selectbox(
                                            "Select a ProLIF condition table:",
                                            prolif_tables_dir,
                                            key=f"select_prolif_table_{dir_name}"
                                        )
                                        pose_label_map = st_funcs.get_pose_labels_for_pose_ids(results_db_path, dir_pose_ids)
                                        selected_prolif_pose_dir = filename_to_pose_id.get(selected_file, {}).get("pose_id")
                                        st.caption(
                                            f"Fingerprints for the pose selected above: "
                                            f"{pose_label_map.get(selected_prolif_pose_dir, selected_file)}"
                                        )
                                        fps_dir_df = st_funcs.get_prolif_fingerprint_for_pose(results_db_path, selected_table_dir, selected_prolif_pose_dir)
                                        reserved_dir = {"pose_id"}
                                        all_columns_dir = st_funcs.get_prolif_all_column_names_by_pose_ids(results_db_path, selected_table_dir, dir_pose_ids)
                                        all_itypes_dir = sorted(set(
                                            col.split("_")[-1] for col in all_columns_dir if col not in reserved_dir
                                        ))
                                        if fps_dir_df is not None and not fps_dir_df.empty:
                                            st.markdown("**Filter by interaction type:**")
                                            filter_cols_dir = st.columns(min(len(all_itypes_dir), 6))
                                            selected_types_dir = []
                                            for i, itype in enumerate(all_itypes_dir):
                                                ck_key = f"prolif_filter_{dir_name}_{itype}"
                                                if ck_key not in st.session_state:
                                                    st.session_state[ck_key] = True
                                                if filter_cols_dir[i % len(filter_cols_dir)].checkbox(itype, key=ck_key):
                                                    selected_types_dir.append(itype)
                                            keep_dir = [c for c in fps_dir_df.columns
                                                        if c in reserved_dir or c.split("_")[-1] in selected_types_dir]
                                            st.dataframe(fps_dir_df[keep_dir], use_container_width=True)
                                        else:
                                            st.warning(f"No fingerprint data found for pose {selected_file} in table {selected_table_dir}.")
                                            selected_types_dir = []

                                        show_all_dir_key = f"show_all_prolif_{dir_name}"
                                        if show_all_dir_key not in st.session_state:
                                            st.session_state[show_all_dir_key] = False
                                        _sp, _col = st.columns([3, 7])
                                        with _col:
                                            if st.button(
                                                f"{'Hide' if st.session_state[show_all_dir_key] else 'Show'} all fingerprints for {selected_table_dir}",
                                                key=f"btn_all_prolif_{dir_name}"
                                            ):
                                                st.session_state[show_all_dir_key] = not st.session_state[show_all_dir_key]
                                                st.rerun()
                                        if st.session_state[show_all_dir_key]:
                                            all_fps_dir_df = st_funcs.get_all_prolif_fingerprints_by_pose_ids(results_db_path, selected_table_dir, dir_pose_ids)
                                            if all_fps_dir_df is not None and not all_fps_dir_df.empty:
                                                all_fps_dir_df.insert(1, "pose_name", all_fps_dir_df["pose_id"].map(lambda pid: pose_label_map.get(pid, str(pid))))
                                                all_fps_dir_df = all_fps_dir_df.sort_values("pose_name").reset_index(drop=True)
                                                keep_all_dir = [c for c in all_fps_dir_df.columns
                                                                if c in reserved_dir | {"pose_name"} or c.split("_")[-1] in selected_types_dir]
                                                st.dataframe(all_fps_dir_df[keep_all_dir], use_container_width=True)
                                            else:
                                                st.warning("Could not retrieve fingerprints for all poses.")
                                st.divider()
                        else:
                            _sp, _col = st.columns([1, 9])
                            with _col:
                                st.button(label, disabled=True, key=key + "_btn")

            if df_mmpbsa_poses_results is None or df_mmpbsa_poses_results.empty:
                st.button(f"No MMPBSA results available for {st.session_state['selected_assay_name']}", disabled=True)
            else:
                with st.expander(f"📈 {st.session_state['selected_assay_name']} MMPBSA Results", expanded=False):
                    # Create a selection box for MMPBSA poses.
                    poses = sorted(df_mmpbsa_poses_results['pose_id'].dropna().tolist())
                    if poses:
                        pose_label_map = st_funcs.get_pose_labels_for_pose_ids(results_db_path, poses)

                        default_pose = st.session_state.get('selected_pose_id', poses[0])
                        if default_pose not in poses:
                            default_pose = poses[0]

                        _sp, _col = st.columns([1, 9])
                        with _col:
                            selected_pose = st.selectbox(
                                f"Select a docked pose: {st.session_state['selected_assay_name']}",
                                poses,
                                key="select_pose_id",
                                index=poses.index(default_pose),
                                format_func=lambda pid: pose_label_map.get(pid, str(pid))
                            )

                        # Collapse sub-panels and update tracking when the selection changes
                        if selected_pose != st.session_state.get('selected_pose_id', None):
                            st.session_state['selected_pose_id'] = selected_pose
                            st.session_state['show_mmpbsa_summary_data'] = False
                            st.session_state['show_mmpbsa_decomposition_data'] = False

                        selected_pose_label = pose_label_map.get(selected_pose, str(selected_pose))
                        _sp, _col = st.columns([1, 9])
                        with _col:
                            st.success(f"Selected Pose: {selected_pose_label}")
                        
                        # Using the selected pose, retrieve de MMPBSA data from the database
                        
                        mmpbsa_data = st_funcs.get_mmpbsa_data_for_pose(results_db_path, selected_pose)
                        
                        # Display the retrieved MMPBSA data
                        if mmpbsa_data is not None and not mmpbsa_data.empty:
                            _sp, _col = st.columns([1, 9])
                            with _col:
                                st.write(mmpbsa_data)
                        
                        
                        ## Add a button to show/hide the MMPBSA data for the selected pose
                        if "show_mmpbsa_summary_data" not in st.session_state:
                            st.session_state["show_mmpbsa_summary_data"] = False
                        
                        _sp, _col = st.columns([1, 9])
                        with _col:
                            if st.button(f"{'Hide' if st.session_state['show_mmpbsa_summary_data'] else 'Show'} MMPBSA Summary Data for {selected_pose_label}"):
                                st.session_state["show_mmpbsa_summary_data"] = not st.session_state["show_mmpbsa_summary_data"]
                                st.rerun()
                        if st.session_state["show_mmpbsa_summary_data"]:
                            # Show sum of the energy components for the selected pose
                            if mmpbsa_data is not None and not mmpbsa_data.empty:
                                total_energy = mmpbsa_data['total'].sum()
                                gas_energy = mmpbsa_data['gas'].sum()
                                ele_energy = mmpbsa_data['ele'].sum()
                                vdw_energy = mmpbsa_data['vdw'].sum()
                                polar_solvation = mmpbsa_data['polar_solvation'].sum()
                                nonpolar_solvation = mmpbsa_data['nonpolar_solvation'].sum()
                                _sp, _col = st.columns([2, 8])
                                with _col:
                                    st.write("MMPBSA Energy Components and Total Energy:")
                                    st.write(f"Total Energy: {total_energy:.2f} kcal/mol")
                                    st.write(f"Gas Energy: {gas_energy:.2f} kcal/mol")
                                    st.write(f"Electrostatic Energy: {ele_energy:.2f} kcal/mol")
                                    st.write(f"Van der Waals Energy: {vdw_energy:.2f} kcal/mol")
                                    st.write(f"Polar Solvation Energy: {polar_solvation:.2f} kcal/mol")
                                    st.write(f"Nonpolar Solvation Energy: {nonpolar_solvation:.2f} kcal/mol")
                        
                        ## Add a button to show per residue decomposition data for the selected pose
                        if "show_mmpbsa_decomposition_data" not in st.session_state:
                            st.session_state["show_mmpbsa_decomposition_data"] = False

                        _sp, _col = st.columns([1, 9])
                        with _col:
                            if st.button(
                                f"{'Hide' if st.session_state['show_mmpbsa_decomposition_data'] else 'Show'} "
                                f"MMPBSA Per-Residue Decomposition Data for {selected_pose_label}"
                            ):
                                st.session_state["show_mmpbsa_decomposition_data"] = not st.session_state["show_mmpbsa_decomposition_data"]
                                st.rerun()

                        # IMPORTANT: keep this outside the button block
                        if st.session_state["show_mmpbsa_decomposition_data"]:
                            if "energy_threshold" not in st.session_state:
                                st.session_state["energy_threshold"] = -1.0

                            _sp, _col = st.columns([2, 8])
                            with _col:
                                st.session_state["energy_threshold"] = st.number_input(
                                    f"Set Energy Threshold for Per-Residue Decomposition Data for Pose {selected_pose} (kcal/mol):",
                                    value=st.session_state["energy_threshold"],
                                    step=0.1,
                                    format="%.2f"
                                )

                            # create a button to show the electrostatic energy profile for the selected pose using the MMPBSA data and the energy threshold
                            if "show_electrostatic_profile" not in st.session_state:
                                st.session_state["show_electrostatic_profile"] = False

                            # Indent to the right (tab-like) using a spacer column
                            _spacer, col_profile = st.columns([2, 8])
                            with col_profile:
                                if st.button(
                                    f"{'Hide' if st.session_state['show_electrostatic_profile'] else 'Show'} "
                                    f"Electrostatic Energy Profile for Pose {selected_pose}",
                                    key=f"btn_show_ele_profile_{selected_pose}"
                                ):
                                    st.session_state["show_electrostatic_profile"] = not st.session_state["show_electrostatic_profile"]
                                    st.rerun()

                            if st.session_state["show_electrostatic_profile"]:
                                df_electrostatic = st_funcs.get_electrostatic_profile_for_pose(
                                    mmpbsa_data, st.session_state["energy_threshold"]
                                )

                                st.write(
                                    f"Electrostatic Energy Profile for Pose {selected_pose} "
                                    f"(Threshold: {st.session_state['energy_threshold']} kcal/mol):"
                                )
                                st.dataframe(df_electrostatic)

                                ele_fig = st_funcs.create_mmpbsa_component_plot(
                                    df_electrostatic,
                                    x_col="residue",
                                    y_col="ele",
                                    title=f"Electrostatic Energy Profile for Pose {selected_pose}",
                                    xlabel="Residue",
                                    ylabel="Electrostatic Energy (kcal/mol)"
                                )
                                st.pyplot(ele_fig, use_container_width=False, clear_figure=True)
                                
                            # create a button to show the vdw energy profile for the selected pose using the MMPBSA data and the energy threshold
                            if "show_vdw_profile" not in st.session_state:
                                st.session_state["show_vdw_profile"] = False

                            # Indent to the right (tab-like) using a spacer column
                            _spacer, col_profile = st.columns([2, 8])
                            with col_profile:
                                if st.button(
                                    f"{'Hide' if st.session_state['show_vdw_profile'] else 'Show'} "
                                    f"Van der Waals Energy Profile for Pose {selected_pose}",
                                    key=f"btn_show_vdw_profile_{selected_pose}"
                                ):
                                    st.session_state["show_vdw_profile"] = not st.session_state["show_vdw_profile"]
                                    st.rerun()

                            if st.session_state["show_vdw_profile"]:
                                df_vdw = st_funcs.get_vdw_profile_for_pose(
                                    mmpbsa_data, st.session_state["energy_threshold"]
                                )

                                st.write(
                                    f"Van der Waals Energy Profile for Pose {selected_pose} "
                                    f"(Threshold: {st.session_state['energy_threshold']} kcal/mol):"
                                )
                                st.dataframe(df_vdw)

                                vdw_fig = st_funcs.create_mmpbsa_component_plot(
                                    df_vdw,
                                    x_col="residue",
                                    y_col="vdw",
                                    title=f"Van der Waals Energy Profile for Pose {selected_pose}",
                                    xlabel="Residue",
                                    ylabel="Van der Waals Energy (kcal/mol)"
                                )
                                st.pyplot(vdw_fig, use_container_width=False, clear_figure=True)
                                
                            
                            # create a button to show the polar_solvation energy profile for the selected pose using the MMPBSA data and the energy threshold
                            if "show_polar_solvation_profile" not in st.session_state:
                                st.session_state["show_polar_solvation_profile"] = False

                            # Indent to the right (tab-like) using a spacer column
                            _spacer, col_profile = st.columns([2, 8])
                            with col_profile:
                                if st.button(
                                    f"{'Hide' if st.session_state['show_polar_solvation_profile'] else 'Show'} "
                                    f"Polar Solvation Energy Profile for Pose {selected_pose}",
                                    key=f"btn_show_polar_solvation_profile_{selected_pose}"
                                ):
                                    st.session_state["show_polar_solvation_profile"] = not st.session_state["show_polar_solvation_profile"]
                                    st.rerun()

                            if st.session_state["show_polar_solvation_profile"]:
                                df_polar_solvation = st_funcs.get_polsolv_profile_for_pose(
                                    mmpbsa_data, st.session_state["energy_threshold"]
                                )

                                st.write(
                                    f"Polar Solvation Energy Profile for Pose {selected_pose} "
                                    f"(Threshold: {st.session_state['energy_threshold']} kcal/mol):"
                                )
                                st.dataframe(df_polar_solvation)

                                polar_solvation_fig = st_funcs.create_mmpbsa_component_plot(
                                    df_polar_solvation,
                                    x_col="residue",
                                    y_col="polar_solvation",
                                    title=f"Polar Solvation Energy Profile for Pose {selected_pose}",
                                    xlabel="Residue",
                                    ylabel="Polar Solvation Energy (kcal/mol)"
                                )
                                st.pyplot(polar_solvation_fig, use_container_width=False, clear_figure=True)
                                
                            
                            # create a button to show the non_polar_solvation energy profile for the selected pose using the MMPBSA data and the energy threshold
                            if "show_non_polar_solvation_profile" not in st.session_state:
                                st.session_state["show_non_polar_solvation_profile"] = False

                            # Indent to the right (tab-like) using a spacer column
                            _spacer, col_profile = st.columns([2, 8])
                            with col_profile:
                                if st.button(
                                    f"{'Hide' if st.session_state['show_non_polar_solvation_profile'] else 'Show'} "
                                    f"Non-Polar Solvation Energy Profile for Pose {selected_pose}",
                                    key=f"btn_show_non_polar_solvation_profile_{selected_pose}"
                                ):
                                    st.session_state["show_non_polar_solvation_profile"] = not st.session_state["show_non_polar_solvation_profile"]
                                    st.rerun()

                            if st.session_state["show_non_polar_solvation_profile"]:
                                df_non_polar_solvation = st_funcs.get_nonpolsolv_profile_for_pose(
                                    mmpbsa_data, st.session_state["energy_threshold"]
                                )

                                st.write(
                                    f"Non-Polar Solvation Energy Profile for Pose {selected_pose} "
                                    f"(Threshold: {st.session_state['energy_threshold']} kcal/mol):"
                                )
                                st.dataframe(df_non_polar_solvation)

                                non_polar_solvation_fig = st_funcs.create_mmpbsa_component_plot(
                                    df_non_polar_solvation,
                                    x_col="residue",
                                    y_col="non_polar_solvation",
                                    title=f"Non-Polar Solvation Energy Profile for Pose {selected_pose}",
                                    xlabel="Residue",
                                    ylabel="Non-Polar Solvation Energy (kcal/mol)"
                                )
                                st.pyplot(non_polar_solvation_fig, use_container_width=False, clear_figure=True)
                                
                                
                                
                            # create a button to show the gas energy profile for the selected pose using the MMPBSA data and the energy threshold
                            if "show_gas_profile" not in st.session_state:
                                st.session_state["show_gas_profile"] = False

                            # Indent to the right (tab-like) using a spacer column
                            _spacer, col_profile = st.columns([2, 8])
                            with col_profile:
                                if st.button(
                                    f"{'Hide' if st.session_state['show_gas_profile'] else 'Show'} "
                                    f"Gas Energy Profile for Pose {selected_pose}",
                                    key=f"btn_show_gas_profile_{selected_pose}"
                                ):
                                    st.session_state["show_gas_profile"] = not st.session_state["show_gas_profile"]
                                    st.rerun()

                            if st.session_state["show_gas_profile"]:
                                df_gas = st_funcs.get_gas_profile_for_pose(
                                    mmpbsa_data, st.session_state["energy_threshold"]
                                )

                                st.write(
                                    f"Gas Energy Profile for Pose {selected_pose} "
                                    f"(Threshold: {st.session_state['energy_threshold']} kcal/mol):"
                                )
                                st.dataframe(df_gas)

                                gas_fig = st_funcs.create_mmpbsa_component_plot(
                                    df_gas,
                                    x_col="residue",
                                    y_col="gas",
                                    title=f"Gas Energy Profile for Pose {selected_pose}",
                                    xlabel="Residue",
                                    ylabel="Gas Energy (kcal/mol)"
                                )
                                st.pyplot(gas_fig, use_container_width=False, clear_figure=True)
                                
                                
                            # create a button to show the total energy profile for the selected pose using the MMPBSA data and the energy threshold
                            if "show_total_profile" not in st.session_state:
                                st.session_state["show_total_profile"] = False

                            # Indent to the right (tab-like) using a spacer column
                            _spacer, col_profile = st.columns([2, 8])
                            with col_profile:
                                if st.button(
                                    f"{'Hide' if st.session_state['show_total_profile'] else 'Show'} "
                                    f"Total Energy Profile for Pose {selected_pose}",
                                    key=f"btn_show_total_profile_{selected_pose}"
                                ):
                                    st.session_state["show_total_profile"] = not st.session_state["show_total_profile"]
                                    st.rerun()

                            if st.session_state["show_total_profile"]:
                                df_total = st_funcs.get_total_profile_for_pose(
                                    mmpbsa_data, st.session_state["energy_threshold"]
                                )

                                st.write(
                                    f"Total Energy Profile for Pose {selected_pose} "
                                    f"(Threshold: {st.session_state['energy_threshold']} kcal/mol):"
                                )
                                st.dataframe(df_total)

                                total_fig = st_funcs.create_mmpbsa_component_plot(
                                    df_total,
                                    x_col="residue",
                                    y_col="total",
                                    title=f"Total Energy Profile for Pose {selected_pose}",
                                    xlabel="Residue",
                                    ylabel="Total Energy (kcal/mol)"
                                )
                                st.pyplot(total_fig, use_container_width=False, clear_figure=True)
                        
                        
                        
                    else:
                        _sp, _col = st.columns([1, 9])
                        with _col:
                            st.warning("No valid pose IDs found in MMPBSA results")

            # ProLIF fingerprints button
            prolif_tables = st_funcs.get_prolif_tables(results_db_path)
            if not prolif_tables:
                st.button("Show ProLIF Fingerprints", disabled=True)
            else:
                with st.expander("🧬 ProLIF Fingerprints", expanded=False):
                    _sp, _col = st.columns([1, 9])
                    with _col:
                        selected_table = st.selectbox(
                            "Select a ProLIF condition table:",
                            prolif_tables,
                            key="select_prolif_table"
                        )
                        prolif_poses_df = st_funcs.get_prolif_poses(results_db_path, selected_table)
                        if prolif_poses_df is not None and not prolif_poses_df.empty:
                            pose_options = prolif_poses_df["pose_id"].tolist()
                            selected_prolif_pose = st.selectbox(
                                "Select a pose to display fingerprints:",
                                pose_options,
                                key="select_prolif_pose"
                            )
                            fps_df = st_funcs.get_prolif_fingerprint_for_pose(results_db_path, selected_table, selected_prolif_pose)

                            if fps_df is not None and not fps_df.empty:
                                ## Build interaction-type filter checkboxes from the union of ALL poses' columns
                                reserved = {"pose_id"}
                                all_columns = st_funcs.get_prolif_all_column_names(results_db_path, selected_table)
                                all_interaction_types = sorted(set(
                                    col.split("_")[-1] for col in all_columns if col not in reserved
                                ))
                                st.markdown("**Filter by interaction type:**")
                                filter_cols = st.columns(min(len(all_interaction_types), 6))
                                selected_types = []
                                for i, itype in enumerate(all_interaction_types):
                                    ck_key = f"prolif_filter_{itype}"
                                    if ck_key not in st.session_state:
                                        st.session_state[ck_key] = True
                                    checked = filter_cols[i % len(filter_cols)].checkbox(itype, key=ck_key)
                                    if checked:
                                        selected_types.append(itype)

                                def filter_fps(df):
                                    keep = [c for c in df.columns
                                            if c in reserved or c.split("_")[-1] in selected_types]
                                    return df[keep]

                                st.dataframe(filter_fps(fps_df), use_container_width=True)
                            else:
                                st.warning(f"No fingerprint data found for pose {selected_prolif_pose}.")

                            ## Button to concatenate and display all poses for the selected condition
                            if "show_all_prolif" not in st.session_state:
                                st.session_state["show_all_prolif"] = False

                            _sp2, _col2 = st.columns([1, 9])
                            with _col2:
                                if st.button(f"{'Hide' if st.session_state['show_all_prolif'] else 'Show'} all fingerprints for {selected_table}"):
                                    st.session_state["show_all_prolif"] = not st.session_state["show_all_prolif"]
                                    st.rerun()

                            if st.session_state["show_all_prolif"]:
                                all_fps_df = st_funcs.get_all_prolif_fingerprints(results_db_path, selected_table)
                                if all_fps_df is not None and not all_fps_df.empty:
                                    st.dataframe(filter_fps(all_fps_df), use_container_width=True)
                                else:
                                    st.warning("Could not retrieve fingerprints for all poses.")
                        else:
                            st.warning("No poses found in the selected ProLIF table.")

            with st.expander("📊 Docking Score Histogram", expanded=False):
                _sp, _col = st.columns([1, 9])
                with _col:
                    # Create a selectbox for unique LigName values if present
                    if 'LigName' in df_results.columns:
                        lig_names = df_results['LigName'].dropna().unique().tolist()
                        if st.session_state.get('selected_lig_name') not in lig_names:
                            st.session_state['selected_lig_name'] = lig_names[0] if lig_names else None
                        selected_lig = st.selectbox(
                            f"Select a Ligand (LigName) in Assay: {st.session_state['selected_assay_name']}",
                            lig_names,
                            key="select_lig_name",
                            index=lig_names.index(st.session_state['selected_lig_name']) if lig_names else 0
                        )
                        # Always update session state and plot on selection change
                        if selected_lig != st.session_state.get('selected_lig_name', None):
                            st.session_state['selected_lig_name'] = selected_lig
                        st.success(f"Selected Ligand: {selected_lig}")
                        histogram = st_funcs.construct_hist_for_ligand(df_results, selected_lig)
                        st.pyplot(histogram, use_container_width=False, clear_figure=True)

            ## Classified Poses Viewer
            with st.expander("🏷️ Classified Poses Viewer", expanded=False):
                _cpv_project_path = st.session_state.get("active_project_path")
                _cpv_assay = st.session_state.get("selected_assay_name")
                _cpv_rf_models = st_funcs.get_rf_prediction_tables(results_db_path)

                ## Reference PDB uploader — shared across both columns
                _sp, _cpv_ref_col = st.columns([1, 9])
                with _cpv_ref_col:
                    _cpv_ref_file = st.file_uploader(
                        "Load a reference PDB for superimposition (optional):",
                        type=["pdb"],
                        key="cpv_reference_pdb_uploader",
                    )
                    if _cpv_ref_file is not None:
                        st.session_state["cpv_reference_pdb_data"] = _cpv_ref_file.read().decode("utf-8")
                        st.success(f"Reference loaded: {_cpv_ref_file.name}")
                    elif "cpv_reference_pdb_data" not in st.session_state:
                        st.session_state["cpv_reference_pdb_data"] = None

                if _cpv_rf_models:
                    ## ── RF predictions path ──────────────────────────────────────────
                    if len(_cpv_rf_models) > 1:
                        _cpv_model = st.selectbox(
                            "RF model:", _cpv_rf_models, key="cpv_model_select"
                        )
                    else:
                        _cpv_model = _cpv_rf_models[0]
                        st.caption(f"RF model: **{_cpv_model}**")

                    _cpv_classified = st_funcs.get_rf_classified_poses(results_db_path, _cpv_model)
                    _cpv_df_pos = _cpv_classified["positive"]
                    _cpv_df_neg = _cpv_classified["negative"]
                    _cpv_df_pos = st_funcs.add_ligand_name_column(_cpv_df_pos, _cpv_project_path, _cpv_assay)
                    _cpv_df_neg = st_funcs.add_ligand_name_column(_cpv_df_neg, _cpv_project_path, _cpv_assay)

                    _cpv_ligname_filter = st.text_input(
                        "Filter by LigName (substring):", key="cpv_ligname_filter"
                    )
                    if _cpv_ligname_filter:
                        _cpv_df_pos = _cpv_df_pos[
                            _cpv_df_pos["LigName"].str.contains(_cpv_ligname_filter, case=False, na=False)
                        ].reset_index(drop=True)
                        _cpv_df_neg = _cpv_df_neg[
                            _cpv_df_neg["LigName"].str.contains(_cpv_ligname_filter, case=False, na=False)
                        ].reset_index(drop=True)

                    _cpv_pred_filter_cols = list(_cpv_df_pos.columns) if not _cpv_df_pos.empty else list(_cpv_df_neg.columns)
                    if _cpv_pred_filter_cols:
                        st.caption("Filter Positive/Negative Predictions by column (case-insensitive substring match):")
                        _cpv_pred_col_filter_boxes = st.columns(len(_cpv_pred_filter_cols))
                        _cpv_pred_col_filter_values = {}
                        for _cpv_pred_fbox, _cpv_pred_fcol in zip(_cpv_pred_col_filter_boxes, _cpv_pred_filter_cols):
                            with _cpv_pred_fbox:
                                _cpv_pred_col_filter_values[_cpv_pred_fcol] = st.text_input(
                                    _cpv_pred_fcol, key=f"cpv_pred_colfilter_{_cpv_pred_fcol}"
                                )
                        for _cpv_pred_fcol, _cpv_pred_fval in _cpv_pred_col_filter_values.items():
                            if _cpv_pred_fval:
                                if _cpv_pred_fcol in _cpv_df_pos.columns:
                                    _cpv_df_pos = _cpv_df_pos[
                                        _cpv_df_pos[_cpv_pred_fcol].astype(str).str.contains(
                                            _cpv_pred_fval, case=False, na=False
                                        )
                                    ].reset_index(drop=True)
                                if _cpv_pred_fcol in _cpv_df_neg.columns:
                                    _cpv_df_neg = _cpv_df_neg[
                                        _cpv_df_neg[_cpv_pred_fcol].astype(str).str.contains(
                                            _cpv_pred_fval, case=False, na=False
                                        )
                                    ].reset_index(drop=True)

                    _cpv_has_pos = not _cpv_df_pos.empty
                    _cpv_has_neg = not _cpv_df_neg.empty

                    if not _cpv_has_pos and not _cpv_has_neg:
                        st.info("No classified poses found for the selected model.")
                    else:
                        _cpv_col_pos, _cpv_col_neg = st.columns(2)

                        with _cpv_col_pos:
                            st.markdown(f"### ✅ Positive Predictions ({len(_cpv_df_pos)})")
                            if not _cpv_has_pos:
                                st.info("No positive predictions.")
                            else:
                                _cpv_pos_labels = [
                                    f"{r['LigName']}_{r['run_number']}  (prob: {r['rf_prob_positive']:.2f})"
                                    for _, r in _cpv_df_pos.iterrows()
                                ]
                                if (st.session_state.get("cpv_pos_selectbox") not in _cpv_pos_labels):
                                    st.session_state["cpv_pos_selectbox"] = _cpv_pos_labels[0]

                                _cpv_pos_cur_idx = _cpv_pos_labels.index(st.session_state["cpv_pos_selectbox"])
                                _cpv_pos_prev_col, _cpv_pos_mid_col, _cpv_pos_next_col = st.columns([1, 4, 1])
                                with _cpv_pos_prev_col:
                                    if st.button("◀ Prev", key="cpv_pos_prev_btn", disabled=(_cpv_pos_cur_idx == 0)):
                                        st.session_state["cpv_pos_selectbox"] = _cpv_pos_labels[_cpv_pos_cur_idx - 1]
                                        st.rerun()
                                with _cpv_pos_next_col:
                                    if st.button("Next ▶", key="cpv_pos_next_btn", disabled=(_cpv_pos_cur_idx == len(_cpv_pos_labels) - 1)):
                                        st.session_state["cpv_pos_selectbox"] = _cpv_pos_labels[_cpv_pos_cur_idx + 1]
                                        st.rerun()
                                with _cpv_pos_mid_col:
                                    st.caption(f"Pose {_cpv_pos_cur_idx + 1} of {len(_cpv_pos_labels)}")

                                _cpv_pos_sel = st.selectbox(
                                    "Select a positive pose:",
                                    _cpv_pos_labels,
                                    key="cpv_pos_selectbox",
                                )
                                _cpv_pos_row = _cpv_df_pos.iloc[_cpv_pos_labels.index(_cpv_pos_sel)]
                                _cpv_pos_pdb_path = st_funcs.find_pose_pdb(
                                    results_db_path, _cpv_pos_row["LigName"], _cpv_pos_row["run_number"]
                                )
                                if _cpv_pos_pdb_path:
                                    with open(_cpv_pos_pdb_path, "r") as _f:
                                        _cpv_pos_pdb = _f.read()
                                    _cpv_pos_view = py3Dmol.view(width=480, height=420)
                                    _cpv_pos_view.addModel(_cpv_pos_pdb, "pdb")
                                    _cpv_pos_view.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                                    if st.session_state.get("cpv_reference_pdb_data"):
                                        _cpv_pos_view.addModel(st.session_state["cpv_reference_pdb_data"], "pdb")
                                        _cpv_pos_view.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                                    _cpv_pos_view.zoomTo()
                                    st.components.v1.html(_cpv_pos_view.write_html(), height=430)
                                    st.caption(
                                        f"Pose ID: {_cpv_pos_row['Pose_ID']}  |  "
                                        f"Score: {_cpv_pos_row['docking_score']:.2f}  |  "
                                        f"Prob: {_cpv_pos_row['rf_prob_positive']:.3f}"
                                    )
                                else:
                                    st.warning(
                                        f"PDB not found for {_cpv_pos_row['LigName']}_{_cpv_pos_row['run_number']}. "
                                        "Extract poses first."
                                    )
                                st.divider()
                                st.dataframe(_cpv_df_pos, use_container_width=True)
                                _cpv_pos_lig_names = _cpv_df_pos["LigName"].unique().tolist()
                                _cpv_pos_compounds, _cpv_pos_err = st_funcs.get_compounds_for_lig_names(
                                    _cpv_project_path, _cpv_assay, _cpv_pos_lig_names
                                )
                                _cpv_btn_pos, _cpv_sub_pos = st.columns(2)
                                with _cpv_btn_pos:
                                    if _cpv_pos_compounds is not None:
                                        st.download_button(
                                            "📥 Export compounds as CSV",
                                            data=_cpv_pos_compounds.to_csv(index=False),
                                            file_name=f"{_cpv_assay}_positive_compounds.csv",
                                            mime="text/csv",
                                            key="cpv_export_pos",
                                        )
                                    else:
                                        st.warning(f"Export unavailable: {_cpv_pos_err}")
                                with _cpv_sub_pos:
                                    if st.button("🗂️ Subset binders", key="cpv_subset_btn_pos"):
                                        st.session_state["cpv_show_subset_pos"] = not st.session_state.get("cpv_show_subset_pos", False)
                                if st.session_state.get("cpv_show_subset_pos"):
                                    with st.form(key="cpv_subset_form_pos"):
                                        _cpv_new_tbl_pos = st.text_input("New table name:", key="cpv_subset_name_pos")
                                        if st.form_submit_button("Create subset"):
                                            _ok, _msg = st_funcs.create_chemspace_subset(
                                                _cpv_project_path, _cpv_assay,
                                                _cpv_pos_lig_names, _cpv_new_tbl_pos
                                            )
                                            if _ok:
                                                st.success(_msg)
                                                st.session_state["cpv_show_subset_pos"] = False
                                            else:
                                                st.error(_msg)

                        with _cpv_col_neg:
                            st.markdown(f"### ❌ Negative Predictions ({len(_cpv_df_neg)})")
                            if not _cpv_has_neg:
                                st.info("No negative predictions.")
                            else:
                                _cpv_neg_labels = [
                                    f"{r['LigName']}_{r['run_number']}  (prob: {r['rf_prob_positive']:.2f})"
                                    for _, r in _cpv_df_neg.iterrows()
                                ]
                                if (st.session_state.get("cpv_neg_selectbox") not in _cpv_neg_labels):
                                    st.session_state["cpv_neg_selectbox"] = _cpv_neg_labels[0]

                                _cpv_neg_cur_idx = _cpv_neg_labels.index(st.session_state["cpv_neg_selectbox"])
                                _cpv_neg_prev_col, _cpv_neg_mid_col, _cpv_neg_next_col = st.columns([1, 4, 1])
                                with _cpv_neg_prev_col:
                                    if st.button("◀ Prev", key="cpv_neg_prev_btn", disabled=(_cpv_neg_cur_idx == 0)):
                                        st.session_state["cpv_neg_selectbox"] = _cpv_neg_labels[_cpv_neg_cur_idx - 1]
                                        st.rerun()
                                with _cpv_neg_next_col:
                                    if st.button("Next ▶", key="cpv_neg_next_btn", disabled=(_cpv_neg_cur_idx == len(_cpv_neg_labels) - 1)):
                                        st.session_state["cpv_neg_selectbox"] = _cpv_neg_labels[_cpv_neg_cur_idx + 1]
                                        st.rerun()
                                with _cpv_neg_mid_col:
                                    st.caption(f"Pose {_cpv_neg_cur_idx + 1} of {len(_cpv_neg_labels)}")

                                _cpv_neg_sel = st.selectbox(
                                    "Select a negative pose:",
                                    _cpv_neg_labels,
                                    key="cpv_neg_selectbox",
                                )
                                _cpv_neg_row = _cpv_df_neg.iloc[_cpv_neg_labels.index(_cpv_neg_sel)]
                                _cpv_neg_pdb_path = st_funcs.find_pose_pdb(
                                    results_db_path, _cpv_neg_row["LigName"], _cpv_neg_row["run_number"]
                                )
                                if _cpv_neg_pdb_path:
                                    with open(_cpv_neg_pdb_path, "r") as _f:
                                        _cpv_neg_pdb = _f.read()
                                    _cpv_neg_view = py3Dmol.view(width=480, height=420)
                                    _cpv_neg_view.addModel(_cpv_neg_pdb, "pdb")
                                    _cpv_neg_view.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                                    if st.session_state.get("cpv_reference_pdb_data"):
                                        _cpv_neg_view.addModel(st.session_state["cpv_reference_pdb_data"], "pdb")
                                        _cpv_neg_view.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                                    _cpv_neg_view.zoomTo()
                                    st.components.v1.html(_cpv_neg_view.write_html(), height=430)
                                    st.caption(
                                        f"Pose ID: {_cpv_neg_row['Pose_ID']}  |  "
                                        f"Score: {_cpv_neg_row['docking_score']:.2f}  |  "
                                        f"Prob: {_cpv_neg_row['rf_prob_positive']:.3f}"
                                    )
                                else:
                                    st.warning(
                                        f"PDB not found for {_cpv_neg_row['LigName']}_{_cpv_neg_row['run_number']}. "
                                        "Extract poses first."
                                    )
                                st.divider()
                                st.dataframe(_cpv_df_neg, use_container_width=True)
                                _cpv_neg_lig_names = _cpv_df_neg["LigName"].unique().tolist()
                                _cpv_neg_compounds, _cpv_neg_err = st_funcs.get_compounds_for_lig_names(
                                    _cpv_project_path, _cpv_assay, _cpv_neg_lig_names
                                )
                                _cpv_btn_neg, _cpv_sub_neg = st.columns(2)
                                with _cpv_btn_neg:
                                    if _cpv_neg_compounds is not None:
                                        st.download_button(
                                            "📥 Export compounds as CSV",
                                            data=_cpv_neg_compounds.to_csv(index=False),
                                            file_name=f"{_cpv_assay}_negative_compounds.csv",
                                            mime="text/csv",
                                            key="cpv_export_neg",
                                        )
                                    else:
                                        st.warning(f"Export unavailable: {_cpv_neg_err}")
                                with _cpv_sub_neg:
                                    if st.button("🗂️ Subset binders", key="cpv_subset_btn_neg"):
                                        st.session_state["cpv_show_subset_neg"] = not st.session_state.get("cpv_show_subset_neg", False)
                                if st.session_state.get("cpv_show_subset_neg"):
                                    with st.form(key="cpv_subset_form_neg"):
                                        _cpv_new_tbl_neg = st.text_input("New table name:", key="cpv_subset_name_neg")
                                        if st.form_submit_button("Create subset"):
                                            _ok, _msg = st_funcs.create_chemspace_subset(
                                                _cpv_project_path, _cpv_assay,
                                                _cpv_neg_lig_names, _cpv_new_tbl_neg
                                            )
                                            if _ok:
                                                st.success(_msg)
                                                st.session_state["cpv_show_subset_neg"] = False
                                            else:
                                                st.error(_msg)

                        st.divider()
                        st.markdown("### 🧬 Per-Ligand Predicted Label Summary")
                        st.caption(
                            "A ligand's poses (run_number) can each get a different "
                            "predicted label — this shows how consistently each ligand "
                            "was classified across its poses."
                        )
                        _cpv_lig_summary = st_funcs.get_rf_ligand_summary(results_db_path, _cpv_model)
                        if _cpv_lig_summary.empty:
                            st.info("No per-ligand summary available.")
                        else:
                            _cpv_lig_min_positive = st.slider(
                                "Minimum positive poses per ligand:", 0,
                                int(_cpv_lig_summary["n_poses"].max()), 0,
                                key="cpv_lig_min_positive",
                            )
                            _cpv_lig_filtered = _cpv_lig_summary[
                                _cpv_lig_summary["n_positive"] >= _cpv_lig_min_positive
                            ].reset_index(drop=True)
                            _cpv_lig_filtered = st_funcs.add_ligand_name_column(
                                _cpv_lig_filtered, _cpv_project_path, _cpv_assay
                            )
                            if _cpv_ligname_filter:
                                _cpv_lig_filtered = _cpv_lig_filtered[
                                    _cpv_lig_filtered["LigName"].str.contains(_cpv_ligname_filter, case=False, na=False)
                                ].reset_index(drop=True)

                            _cpv_lig_col_filter_boxes = st.columns(len(_cpv_lig_filtered.columns))
                            _cpv_lig_col_filter_values = {}
                            for _cpv_lig_fbox, _cpv_lig_fcol in zip(_cpv_lig_col_filter_boxes, _cpv_lig_filtered.columns):
                                with _cpv_lig_fbox:
                                    _cpv_lig_col_filter_values[_cpv_lig_fcol] = st.text_input(
                                        _cpv_lig_fcol, key=f"cpv_lig_colfilter_{_cpv_lig_fcol}"
                                    )
                            for _cpv_lig_fcol, _cpv_lig_fval in _cpv_lig_col_filter_values.items():
                                if _cpv_lig_fval:
                                    _cpv_lig_filtered = _cpv_lig_filtered[
                                        _cpv_lig_filtered[_cpv_lig_fcol].astype(str).str.contains(
                                            _cpv_lig_fval, case=False, na=False
                                        )
                                    ].reset_index(drop=True)

                            st.dataframe(_cpv_lig_filtered, use_container_width=True)
                            st.download_button(
                                "📥 Export per-ligand summary as CSV",
                                data=_cpv_lig_filtered.to_csv(index=False),
                                file_name=f"{_cpv_assay}_{_cpv_model}_ligand_summary.csv",
                                mime="text/csv",
                                key="cpv_export_ligand_summary",
                            )

                        st.divider()
                        st.markdown("### 🧪 Per-Compound Predicted Label Summary")
                        st.caption(
                            "Groups stereoisomers of the same compound together (LigNames "
                            "sharing the first InChIKey block, before the first dash) and "
                            "reports averaged per-ligand statistics across the group."
                        )
                        _cpv_cmpd_summary = st_funcs.get_rf_compound_summary(
                            results_db_path, _cpv_model, _cpv_project_path, _cpv_assay
                        )
                        if _cpv_cmpd_summary.empty:
                            st.info("No per-compound summary available.")
                        else:
                            _cpv_cmpd_filtered = _cpv_cmpd_summary
                            if _cpv_ligname_filter:
                                _cpv_cmpd_filtered = _cpv_cmpd_filtered[
                                    _cpv_cmpd_filtered["compound_key"].str.contains(_cpv_ligname_filter, case=False, na=False)
                                ].reset_index(drop=True)
                            st.dataframe(_cpv_cmpd_filtered, use_container_width=True)
                            st.bar_chart(
                                _cpv_cmpd_filtered.set_index("compound_key")["avg_frac_positive"].head(30)
                            )
                            st.download_button(
                                "📥 Export per-compound summary as CSV",
                                data=_cpv_cmpd_filtered.to_csv(index=False),
                                file_name=f"{_cpv_assay}_{_cpv_model}_compound_summary.csv",
                                mime="text/csv",
                                key="cpv_export_compound_summary",
                            )

                else:
                    ## ── Fallback: manual positive/negative binder flags ───────────
                    _cpv_df_pos = st_funcs.get_binders_registry(_cpv_project_path, "positive")
                    _cpv_df_neg = st_funcs.get_binders_registry(_cpv_project_path, "negative")

                    if _cpv_df_pos is not None and not _cpv_df_pos.empty:
                        _cpv_df_pos = _cpv_df_pos[_cpv_df_pos["assay_name"] == _cpv_assay].reset_index(drop=True)
                    if _cpv_df_neg is not None and not _cpv_df_neg.empty:
                        _cpv_df_neg = _cpv_df_neg[_cpv_df_neg["assay_name"] == _cpv_assay].reset_index(drop=True)

                    _cpv_has_pos = _cpv_df_pos is not None and not _cpv_df_pos.empty
                    _cpv_has_neg = _cpv_df_neg is not None and not _cpv_df_neg.empty

                    if not _cpv_has_pos and not _cpv_has_neg:
                        st.info(f"No classified poses for assay '{_cpv_assay}'.")
                    else:
                        _cpv_col_pos, _cpv_col_neg = st.columns(2)

                        with _cpv_col_pos:
                            st.markdown(f"### ✅ Positive Binders ({len(_cpv_df_pos) if _cpv_has_pos else 0})")
                            if not _cpv_has_pos:
                                st.info("No positive binders for this assay.")
                            else:
                                _cpv_pos_labels = [
                                    f"{r['pose_file']}  [{r['directory']}]"
                                    for _, r in _cpv_df_pos.iterrows()
                                ]
                                _cpv_pos_sel = st.selectbox(
                                    "Select a positive pose:", _cpv_pos_labels, key="cpv_pos_selectbox"
                                )
                                _cpv_pos_row = _cpv_df_pos.iloc[_cpv_pos_labels.index(_cpv_pos_sel)]
                                _cpv_pos_path = _cpv_pos_row["pose_full_path"]
                                if os.path.exists(_cpv_pos_path):
                                    with open(_cpv_pos_path, "r") as _f:
                                        _cpv_pos_pdb = _f.read()
                                    _cpv_pos_view = py3Dmol.view(width=480, height=420)
                                    _cpv_pos_view.addModel(_cpv_pos_pdb, "pdb")
                                    _cpv_pos_view.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                                    if st.session_state.get("cpv_reference_pdb_data"):
                                        _cpv_pos_view.addModel(st.session_state["cpv_reference_pdb_data"], "pdb")
                                        _cpv_pos_view.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                                    _cpv_pos_view.zoomTo()
                                    st.components.v1.html(_cpv_pos_view.write_html(), height=430)
                                    st.caption(
                                        f"{_cpv_pos_row['pose_file']}  |  {_cpv_pos_row['directory']}  |  "
                                        f"flagged: {_cpv_pos_row['flagged_at'][:10]}"
                                    )
                                else:
                                    st.warning(f"PDB file not found: {_cpv_pos_path}")
                                if st.button("Remove from positive binders", key="cpv_btn_remove_pos"):
                                    _cpv_rem = st_funcs.remove_binder(
                                        project_path=_cpv_project_path,
                                        binder_type="positive",
                                        assay_name=_cpv_pos_row["assay_name"],
                                        pose_file=_cpv_pos_row["pose_file"],
                                        directory=_cpv_pos_row["directory"],
                                    )
                                    if _cpv_rem == "removed":
                                        st.success(f"Removed: {_cpv_pos_row['pose_file']}")
                                        st.rerun()
                                    else:
                                        st.error(f"Could not remove: {_cpv_rem}")

                        with _cpv_col_neg:
                            st.markdown(f"### ❌ Negative Binders ({len(_cpv_df_neg) if _cpv_has_neg else 0})")
                            if not _cpv_has_neg:
                                st.info("No negative binders for this assay.")
                            else:
                                _cpv_neg_labels = [
                                    f"{r['pose_file']}  [{r['directory']}]"
                                    for _, r in _cpv_df_neg.iterrows()
                                ]
                                _cpv_neg_sel = st.selectbox(
                                    "Select a negative pose:", _cpv_neg_labels, key="cpv_neg_selectbox"
                                )
                                _cpv_neg_row = _cpv_df_neg.iloc[_cpv_neg_labels.index(_cpv_neg_sel)]
                                _cpv_neg_path = _cpv_neg_row["pose_full_path"]
                                if os.path.exists(_cpv_neg_path):
                                    with open(_cpv_neg_path, "r") as _f:
                                        _cpv_neg_pdb = _f.read()
                                    _cpv_neg_view = py3Dmol.view(width=480, height=420)
                                    _cpv_neg_view.addModel(_cpv_neg_pdb, "pdb")
                                    _cpv_neg_view.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                                    if st.session_state.get("cpv_reference_pdb_data"):
                                        _cpv_neg_view.addModel(st.session_state["cpv_reference_pdb_data"], "pdb")
                                        _cpv_neg_view.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                                    _cpv_neg_view.zoomTo()
                                    st.components.v1.html(_cpv_neg_view.write_html(), height=430)
                                    st.caption(
                                        f"{_cpv_neg_row['pose_file']}  |  {_cpv_neg_row['directory']}  |  "
                                        f"flagged: {_cpv_neg_row['flagged_at'][:10]}"
                                    )
                                else:
                                    st.warning(f"PDB file not found: {_cpv_neg_path}")
                                if st.button("Remove from negative binders", key="cpv_btn_remove_neg"):
                                    _cpv_rem = st_funcs.remove_binder(
                                        project_path=_cpv_project_path,
                                        binder_type="negative",
                                        assay_name=_cpv_neg_row["assay_name"],
                                        pose_file=_cpv_neg_row["pose_file"],
                                        directory=_cpv_neg_row["directory"],
                                    )
                                    if _cpv_rem == "removed":
                                        st.success(f"Removed: {_cpv_neg_row['pose_file']}")
                                        st.rerun()
                                    else:
                                        st.error(f"Could not remove: {_cpv_rem}")


elif page == "ML features management":
    st.title("Machine Learning")

    if "active_project_path" not in st.session_state:
        st.warning("No active project. Please select a project first.")
    else:
        project_path = st.session_state["active_project_path"]

        with st.expander("📚 Training Set Registries", expanded=False):
            if "ml_all_tables_expanded" not in st.session_state:
                st.session_state["ml_all_tables_expanded"] = False
            if st.button(
                "📂 Expand All Tables" if not st.session_state["ml_all_tables_expanded"] else "📁 Collapse All Tables",
                key="btn_ml_toggle_all_tables",
            ):
                st.session_state["ml_all_tables_expanded"] = not st.session_state["ml_all_tables_expanded"]
                st.rerun()

            _ts_indent, _ts_content = st.columns([0.05, 0.95])
            with _ts_content:
                ## --- Positive binders ---
                st.markdown("### ✅ Positive Binders")
                df_pos = st_funcs.get_binders_registry(project_path, "positive")
                if df_pos is None or df_pos.empty:
                    st.info("No positive binders registered yet.")
                else:
                    st.success(f"{len(df_pos)} positive binder(s) registered.")
                    with st.expander("📋 Positive Binders Table", expanded=st.session_state.get("ml_all_tables_expanded", False)):
                        _df_pos_display = df_pos.copy()
                        _df_pos_display.insert(0, "Select", False)
                        _edited_pos = st.data_editor(
                            _df_pos_display,
                            column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                            disabled=[c for c in _df_pos_display.columns if c != "Select"],
                            hide_index=True,
                            use_container_width=True,
                            key="data_editor_pos_binders",
                        )
                    _selected_pos = _edited_pos[_edited_pos["Select"]]

                    ## --- View / Clear / Delete Selected Positive Poses ---
                    pos_viewer_key = "show_ml_positive_poses"
                    if pos_viewer_key not in st.session_state:
                        st.session_state[pos_viewer_key] = False
                    _col_pos_view, _col_pos_clear, _col_pos_delete, _col_pos_move = st.columns([3, 2, 2, 3])
                    with _col_pos_view:
                        st.button(
                            f"{'Hide' if st.session_state[pos_viewer_key] else 'View'} Positive Poses",
                            key="btn_ml_positive_poses",
                            on_click=_toggle_ml_positive_viewer
                        )
                    with _col_pos_clear:
                        if not st.session_state.get("confirm_clear_positive"):
                            if st.button("🗑️ Clear Positive Binders", key="btn_clear_positive_binders"):
                                st.session_state["confirm_clear_positive"] = True
                                st.rerun()
                        else:
                            st.warning("⚠️ This will delete **all** positive binders. Are you sure?")
                            _ccpos1, _ccpos2 = st.columns(2)
                            with _ccpos1:
                                if st.button("✅ Yes, clear all", key="btn_confirm_clear_positive"):
                                    result = st_funcs.clear_binders(project_path, "positive")
                                    st.session_state["confirm_clear_positive"] = False
                                    if result == "cleared":
                                        st.session_state[pos_viewer_key] = False
                                        st.rerun()
                                    else:
                                        st.error(f"Could not clear positive binders: {result}")
                            with _ccpos2:
                                if st.button("❌ Cancel", key="btn_cancel_clear_positive"):
                                    st.session_state["confirm_clear_positive"] = False
                                    st.rerun()
                    with _col_pos_delete:
                        _n_sel_pos = len(_selected_pos)
                        if not st.session_state.get("confirm_delete_selected_positive"):
                            if st.button(
                                f"🗑️ Delete Selected ({_n_sel_pos})",
                                key="btn_delete_selected_positive",
                                disabled=(_n_sel_pos == 0),
                            ):
                                st.session_state["confirm_delete_selected_positive"] = True
                                st.rerun()
                        else:
                            st.warning(f"⚠️ Delete **{_n_sel_pos}** selected row(s)? Are you sure?")
                            _cdpos1, _cdpos2 = st.columns(2)
                            with _cdpos1:
                                if st.button("✅ Yes, delete", key="btn_confirm_delete_selected_positive"):
                                    for _, _row in _selected_pos.iterrows():
                                        st_funcs.remove_binder(
                                            project_path=project_path,
                                            binder_type="positive",
                                            assay_name=_row["assay_name"],
                                            pose_file=_row["pose_file"],
                                            directory=_row["directory"],
                                        )
                                    st.session_state["confirm_delete_selected_positive"] = False
                                    st.rerun()
                            with _cdpos2:
                                if st.button("❌ Cancel", key="btn_cancel_delete_selected_positive"):
                                    st.session_state["confirm_delete_selected_positive"] = False
                                    st.rerun()
                    with _col_pos_move:
                        _n_sel_pos_move = len(_selected_pos)
                        if not st.session_state.get("confirm_move_selected_positive"):
                            if st.button(
                                f"➡️ Send to Negative Binder Table ({_n_sel_pos_move})",
                                key="btn_move_selected_positive",
                                disabled=(_n_sel_pos_move == 0),
                            ):
                                st.session_state["confirm_move_selected_positive"] = True
                                st.rerun()
                        else:
                            st.warning(f"⚠️ Move **{_n_sel_pos_move}** selected row(s) to the Negative Binders table?")
                            _cmpos1, _cmpos2 = st.columns(2)
                            with _cmpos1:
                                if st.button("✅ Yes, move", key="btn_confirm_move_selected_positive"):
                                    _move_errors = []
                                    for _, _row in _selected_pos.iterrows():
                                        _mv_result = st_funcs.move_binder_to_negative(
                                            project_path=project_path,
                                            assay_name=_row["assay_name"],
                                            pose_file=_row["pose_file"],
                                            directory=_row["directory"],
                                            pose_full_path=_row["pose_full_path"],
                                        )
                                        if _mv_result not in ("moved", "duplicate"):
                                            _move_errors.append(f"{_row['pose_file']}: {_mv_result}")
                                    st.session_state["confirm_move_selected_positive"] = False
                                    if _move_errors:
                                        st.error("Some rows could not be moved: " + "; ".join(_move_errors))
                                    st.rerun()
                            with _cmpos2:
                                if st.button("❌ Cancel", key="btn_cancel_move_selected_positive"):
                                    st.session_state["confirm_move_selected_positive"] = False
                                    st.rerun()

                    if st.session_state[pos_viewer_key]:
                        ref_file_pos = st.file_uploader(
                            "Load a reference PDB for superimposition (optional):",
                            type=["pdb"],
                            key="ml_pos_ref_pdb_uploader"
                        )
                        if ref_file_pos is not None:
                            st.session_state["ml_pos_ref_pdb_data"] = ref_file_pos.read().decode("utf-8")
                            st.success(f"Reference loaded: {ref_file_pos.name}")
                        elif "ml_pos_ref_pdb_data" not in st.session_state:
                            st.session_state["ml_pos_ref_pdb_data"] = None

                        pose_paths_pos = df_pos["pose_full_path"].tolist()
                        pose_labels_pos = [
                            f"{row['assay_name']} / {row['pose_file']}"
                            for _, row in df_pos.iterrows()
                        ]
                        idx_key_pos = "ml_pos_pose_idx"
                        if idx_key_pos not in st.session_state:
                            st.session_state[idx_key_pos] = 0
                        st.session_state[idx_key_pos] = max(0, min(st.session_state[idx_key_pos], len(pose_labels_pos) - 1))

                        selected_pos = st.selectbox(
                            "Select a positive binder pose:",
                            pose_labels_pos,
                            index=st.session_state[idx_key_pos],
                        )
                        if selected_pos:
                            st.session_state[idx_key_pos] = pose_labels_pos.index(selected_pos)

                        full_path_pos = pose_paths_pos[st.session_state[idx_key_pos]]
                        if os.path.exists(full_path_pos):
                            with open(full_path_pos, "r") as f:
                                pdb_data_pos = f.read()
                            view_pos = py3Dmol.view(width=800, height=500)
                            view_pos.addModel(pdb_data_pos, "pdb")
                            view_pos.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                            if st.session_state.get("ml_pos_ref_pdb_data"):
                                view_pos.addModel(st.session_state["ml_pos_ref_pdb_data"], "pdb")
                                view_pos.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                            view_pos.zoomTo()

                            cur_pos = st.session_state[idx_key_pos]
                            has_ref_pos = st.session_state.get("ml_pos_ref_pdb_data")
                            prev_col, viewer_col, next_col = st.columns([0.5, 8, 0.5 if not has_ref_pos else 1])
                            with prev_col:
                                st.write(""); st.write(""); st.write("")
                                if st.button("◀", key="ml_pos_prev", disabled=(cur_pos == 0)):
                                    st.session_state[idx_key_pos] = cur_pos - 1
                                    st.rerun()
                            with viewer_col:
                                st.components.v1.html(view_pos.write_html(), height=520)
                            with next_col:
                                st.write(""); st.write(""); st.write("")
                                if st.button("▶", key="ml_pos_next", disabled=(cur_pos >= len(pose_labels_pos) - 1)):
                                    st.session_state[idx_key_pos] = cur_pos + 1
                                    st.rerun()
                            st.caption(f"Pose {cur_pos + 1} of {len(pose_labels_pos)}")

                            ## VMD script creation
                            with st.expander("🎬 Create VMD Script", expanded=False):
                                _default_vmd_path_pos = os.path.join(
                                    os.path.expanduser("~"), "Desktop",
                                    os.path.splitext(os.path.basename(full_path_pos))[0] + "_vmd.tcl"
                                )
                                _vmd_path_pos = st.text_input(
                                    "Save script to:",
                                    value=_default_vmd_path_pos,
                                    key=f"vmd_path_ml_pos_{cur_pos}",
                                )
                                if st.button("💾 Save VMD Script", key=f"btn_save_vmd_ml_pos_{cur_pos}"):
                                    try:
                                        _ref_pdb_path_pos = None
                                        if st.session_state.get("ml_pos_ref_pdb_data"):
                                            _ref_pdb_path_pos = os.path.join(
                                                os.path.dirname(os.path.abspath(_vmd_path_pos.strip())),
                                                "reference.pdb"
                                            )
                                            with open(_ref_pdb_path_pos, "w") as _rf:
                                                _rf.write(st.session_state["ml_pos_ref_pdb_data"])
                                        _script_pos = st_funcs.generate_vmd_script(
                                            os.path.abspath(full_path_pos), _ref_pdb_path_pos
                                        )
                                        _vmd_out_pos = _vmd_path_pos.strip()
                                        os.makedirs(os.path.dirname(os.path.abspath(_vmd_out_pos)), exist_ok=True)
                                        with open(_vmd_out_pos, "w") as _sf:
                                            _sf.write(_script_pos)
                                        st.success(f"VMD script saved to: {_vmd_out_pos}")
                                        if _ref_pdb_path_pos:
                                            st.info(f"Reference PDB saved alongside: {_ref_pdb_path_pos}")
                                    except Exception as _e:
                                        st.error(f"Could not save VMD script: {_e}")

                            ## Remove current pose from positive registry
                            current_row_pos = df_pos.iloc[st.session_state[idx_key_pos]]
                            if st.button("🗑️ Remove pose from set", key=f"ml_pos_remove_{st.session_state[idx_key_pos]}"):
                                result = st_funcs.remove_binder(
                                    project_path=project_path,
                                    binder_type="positive",
                                    assay_name=current_row_pos["assay_name"],
                                    pose_file=current_row_pos["pose_file"],
                                    directory=current_row_pos["directory"],
                                )
                                if result == "removed":
                                    st.session_state[idx_key_pos] = max(0, st.session_state[idx_key_pos] - 1)
                                    st.rerun()
                                else:
                                    st.error(f"Could not remove pose: {result}")
                        else:
                            st.warning(f"Pose file not found on disk: {full_path_pos}")

                st.divider()

                ## --- Negative binders ---
                st.markdown("### ❌ Negative Binders")
                df_neg = st_funcs.get_binders_registry(project_path, "negative")
                if df_neg is None or df_neg.empty:
                    st.info("No negative binders registered yet.")
                else:
                    st.success(f"{len(df_neg)} negative binder(s) registered.")
                    with st.expander("📋 Negative Binders Table", expanded=st.session_state.get("ml_all_tables_expanded", False)):
                        _df_neg_display = df_neg.copy()
                        _df_neg_display.insert(0, "Select", False)
                        _edited_neg = st.data_editor(
                            _df_neg_display,
                            column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                            disabled=[c for c in _df_neg_display.columns if c != "Select"],
                            hide_index=True,
                            use_container_width=True,
                            key="data_editor_neg_binders",
                        )
                    _selected_neg = _edited_neg[_edited_neg["Select"]]

                    ## --- View / Clear / Delete Selected Negative Poses ---
                    neg_viewer_key = "show_ml_negative_poses"
                    if neg_viewer_key not in st.session_state:
                        st.session_state[neg_viewer_key] = False
                    _col_neg_view, _col_neg_clear, _col_neg_delete, _col_neg_move = st.columns([3, 2, 2, 3])
                    with _col_neg_view:
                        st.button(
                            f"{'Hide' if st.session_state[neg_viewer_key] else 'View'} Negative Poses",
                            key="btn_ml_negative_poses",
                            on_click=_toggle_ml_negative_viewer
                        )
                    with _col_neg_clear:
                        if not st.session_state.get("confirm_clear_negative"):
                            if st.button("🗑️ Clear Negative Binders", key="btn_clear_negative_binders"):
                                st.session_state["confirm_clear_negative"] = True
                                st.rerun()
                        else:
                            st.warning("⚠️ This will delete **all** negative binders. Are you sure?")
                            _ccneg1, _ccneg2 = st.columns(2)
                            with _ccneg1:
                                if st.button("✅ Yes, clear all", key="btn_confirm_clear_negative"):
                                    result = st_funcs.clear_binders(project_path, "negative")
                                    st.session_state["confirm_clear_negative"] = False
                                    if result == "cleared":
                                        st.session_state[neg_viewer_key] = False
                                        st.rerun()
                                    else:
                                        st.error(f"Could not clear negative binders: {result}")
                            with _ccneg2:
                                if st.button("❌ Cancel", key="btn_cancel_clear_negative"):
                                    st.session_state["confirm_clear_negative"] = False
                                    st.rerun()
                    with _col_neg_delete:
                        _n_sel_neg = len(_selected_neg)
                        if not st.session_state.get("confirm_delete_selected_negative"):
                            if st.button(
                                f"🗑️ Delete Selected ({_n_sel_neg})",
                                key="btn_delete_selected_negative",
                                disabled=(_n_sel_neg == 0),
                            ):
                                st.session_state["confirm_delete_selected_negative"] = True
                                st.rerun()
                        else:
                            st.warning(f"⚠️ Delete **{_n_sel_neg}** selected row(s)? Are you sure?")
                            _cdneg1, _cdneg2 = st.columns(2)
                            with _cdneg1:
                                if st.button("✅ Yes, delete", key="btn_confirm_delete_selected_negative"):
                                    for _, _row in _selected_neg.iterrows():
                                        st_funcs.remove_binder(
                                            project_path=project_path,
                                            binder_type="negative",
                                            assay_name=_row["assay_name"],
                                            pose_file=_row["pose_file"],
                                            directory=_row["directory"],
                                        )
                                    st.session_state["confirm_delete_selected_negative"] = False
                                    st.rerun()
                            with _cdneg2:
                                if st.button("❌ Cancel", key="btn_cancel_delete_selected_negative"):
                                    st.session_state["confirm_delete_selected_negative"] = False
                                    st.rerun()
                    with _col_neg_move:
                        _n_sel_neg_move = len(_selected_neg)
                        if not st.session_state.get("confirm_move_selected_negative"):
                            if st.button(
                                f"➡️ Send to Positive Binder Table ({_n_sel_neg_move})",
                                key="btn_move_selected_negative",
                                disabled=(_n_sel_neg_move == 0),
                            ):
                                st.session_state["confirm_move_selected_negative"] = True
                                st.rerun()
                        else:
                            st.warning(f"⚠️ Move **{_n_sel_neg_move}** selected row(s) to the Positive Binders table?")
                            _cmneg1, _cmneg2 = st.columns(2)
                            with _cmneg1:
                                if st.button("✅ Yes, move", key="btn_confirm_move_selected_negative"):
                                    _move_errors_neg = []
                                    for _, _row in _selected_neg.iterrows():
                                        _mv_result = st_funcs.move_binder_to_positive(
                                            project_path=project_path,
                                            assay_name=_row["assay_name"],
                                            pose_file=_row["pose_file"],
                                            directory=_row["directory"],
                                            pose_full_path=_row["pose_full_path"],
                                        )
                                        if _mv_result not in ("moved", "duplicate"):
                                            _move_errors_neg.append(f"{_row['pose_file']}: {_mv_result}")
                                    st.session_state["confirm_move_selected_negative"] = False
                                    if _move_errors_neg:
                                        st.error("Some rows could not be moved: " + "; ".join(_move_errors_neg))
                                    st.rerun()
                            with _cmneg2:
                                if st.button("❌ Cancel", key="btn_cancel_move_selected_negative"):
                                    st.session_state["confirm_move_selected_negative"] = False
                                    st.rerun()

                    if st.session_state[neg_viewer_key]:
                        ref_file_neg = st.file_uploader(
                            "Load a reference PDB for superimposition (optional):",
                            type=["pdb"],
                            key="ml_neg_ref_pdb_uploader"
                        )
                        if ref_file_neg is not None:
                            st.session_state["ml_neg_ref_pdb_data"] = ref_file_neg.read().decode("utf-8")
                            st.success(f"Reference loaded: {ref_file_neg.name}")
                        elif "ml_neg_ref_pdb_data" not in st.session_state:
                            st.session_state["ml_neg_ref_pdb_data"] = None

                        pose_paths_neg = df_neg["pose_full_path"].tolist()
                        pose_labels_neg = [
                            f"{row['assay_name']} / {row['pose_file']}"
                            for _, row in df_neg.iterrows()
                        ]
                        idx_key_neg = "ml_neg_pose_idx"
                        if idx_key_neg not in st.session_state:
                            st.session_state[idx_key_neg] = 0
                        st.session_state[idx_key_neg] = max(0, min(st.session_state[idx_key_neg], len(pose_labels_neg) - 1))

                        selected_neg = st.selectbox(
                            "Select a negative binder pose:",
                            pose_labels_neg,
                            index=st.session_state[idx_key_neg],
                        )
                        if selected_neg:
                            st.session_state[idx_key_neg] = pose_labels_neg.index(selected_neg)

                        full_path_neg = pose_paths_neg[st.session_state[idx_key_neg]]
                        if os.path.exists(full_path_neg):
                            with open(full_path_neg, "r") as f:
                                pdb_data_neg = f.read()
                            view_neg = py3Dmol.view(width=800, height=500)
                            view_neg.addModel(pdb_data_neg, "pdb")
                            view_neg.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                            if st.session_state.get("ml_neg_ref_pdb_data"):
                                view_neg.addModel(st.session_state["ml_neg_ref_pdb_data"], "pdb")
                                view_neg.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                            view_neg.zoomTo()

                            cur_neg = st.session_state[idx_key_neg]
                            has_ref_neg = st.session_state.get("ml_neg_ref_pdb_data")
                            prev_col, viewer_col, next_col = st.columns([0.5, 8, 0.5 if not has_ref_neg else 1])
                            with prev_col:
                                st.write(""); st.write(""); st.write("")
                                if st.button("◀", key="ml_neg_prev", disabled=(cur_neg == 0)):
                                    st.session_state[idx_key_neg] = cur_neg - 1
                                    st.rerun()
                            with viewer_col:
                                st.components.v1.html(view_neg.write_html(), height=520)
                            with next_col:
                                st.write(""); st.write(""); st.write("")
                                if st.button("▶", key="ml_neg_next", disabled=(cur_neg >= len(pose_labels_neg) - 1)):
                                    st.session_state[idx_key_neg] = cur_neg + 1
                                    st.rerun()
                            st.caption(f"Pose {cur_neg + 1} of {len(pose_labels_neg)}")

                            ## VMD script creation
                            with st.expander("🎬 Create VMD Script", expanded=False):
                                _default_vmd_path_neg = os.path.join(
                                    os.path.expanduser("~"), "Desktop",
                                    os.path.splitext(os.path.basename(full_path_neg))[0] + "_vmd.tcl"
                                )
                                _vmd_path_neg = st.text_input(
                                    "Save script to:",
                                    value=_default_vmd_path_neg,
                                    key=f"vmd_path_ml_neg_{cur_neg}",
                                )
                                if st.button("💾 Save VMD Script", key=f"btn_save_vmd_ml_neg_{cur_neg}"):
                                    try:
                                        _ref_pdb_path_neg = None
                                        if st.session_state.get("ml_neg_ref_pdb_data"):
                                            _ref_pdb_path_neg = os.path.join(
                                                os.path.dirname(os.path.abspath(_vmd_path_neg.strip())),
                                                "reference.pdb"
                                            )
                                            with open(_ref_pdb_path_neg, "w") as _rf:
                                                _rf.write(st.session_state["ml_neg_ref_pdb_data"])
                                        _script_neg = st_funcs.generate_vmd_script(
                                            os.path.abspath(full_path_neg), _ref_pdb_path_neg
                                        )
                                        _vmd_out_neg = _vmd_path_neg.strip()
                                        os.makedirs(os.path.dirname(os.path.abspath(_vmd_out_neg)), exist_ok=True)
                                        with open(_vmd_out_neg, "w") as _sf:
                                            _sf.write(_script_neg)
                                        st.success(f"VMD script saved to: {_vmd_out_neg}")
                                        if _ref_pdb_path_neg:
                                            st.info(f"Reference PDB saved alongside: {_ref_pdb_path_neg}")
                                    except Exception as _e:
                                        st.error(f"Could not save VMD script: {_e}")

                            ## Remove current pose from negative registry
                            current_row_neg = df_neg.iloc[st.session_state[idx_key_neg]]
                            if st.button("🗑️ Remove pose from set", key=f"ml_neg_remove_{st.session_state[idx_key_neg]}"):
                                result = st_funcs.remove_binder(
                                    project_path=project_path,
                                    binder_type="negative",
                                    assay_name=current_row_neg["assay_name"],
                                    pose_file=current_row_neg["pose_file"],
                                    directory=current_row_neg["directory"],
                                )
                                if result == "removed":
                                    st.session_state[idx_key_neg] = max(0, st.session_state[idx_key_neg] - 1)
                                    st.rerun()
                                else:
                                    st.error(f"Could not remove pose: {result}")
                        else:
                            st.warning(f"Pose file not found on disk: {full_path_neg}")

        st.divider()

        ## --- Consolidate training set ---
        with st.expander("📦 Consolidate Training Set", expanded=False):
            st.markdown(
                "Create a named snapshot of the current positive and negative binders registries "
                "for use in model training."
            )

            with st.form("consolidate_ts_form"):
                ts_id = st.text_input("Training Set ID", placeholder="e.g. ts_run1_2026")
                ts_notes = st.text_area("Notes (optional)", "")
                submitted = st.form_submit_button("Consolidate")

            if submitted:
                if not ts_id.strip():
                    st.error("Please provide a Training Set ID.")
                else:
                    n_pos = len(df_pos) if df_pos is not None and not df_pos.empty else 0
                    n_neg = len(df_neg) if df_neg is not None and not df_neg.empty else 0
                    result = st_funcs.consolidate_training_set(project_path, ts_id.strip(), ts_notes.strip())
                    if result == "saved":
                        st.success(f"✅ Training set **'{ts_id.strip()}'** consolidated — {n_pos} positive(s), {n_neg} negative(s).")
                    elif result == "duplicate":
                        st.warning(f"⚠️ A training set with ID **'{ts_id.strip()}'** already exists. Choose a different ID.")
                    else:
                        st.error(f"❌ {result}")

        st.divider()
        with st.expander("🗄️ Existing Training Set Snapshots", expanded=False):
            df_snaps = st_funcs.get_training_set_snapshots(project_path)
            if df_snaps is None or df_snaps.empty:
                st.info("No training sets consolidated yet.")
            else:
                _fp_status = st_funcs.get_training_set_fingerprint_status(project_path)
                _prolif_conditions_used = st_funcs.get_training_set_prolif_conditions_used(project_path)
                _df_snaps_display = df_snaps.copy()
                _notes_loc = _df_snaps_display.columns.get_loc("notes") + 1
                _df_snaps_display.insert(
                    _notes_loc,
                    "fingerprints",
                    _df_snaps_display["training_set_id"].map(lambda ts: _fp_status.get(ts, "❌ None")),
                )
                _df_snaps_display.insert(
                    _notes_loc + 1,
                    "prolif_conditions",
                    _df_snaps_display["training_set_id"].map(lambda ts: _prolif_conditions_used.get(ts, "—")),
                )
                _df_snaps_display.insert(0, "Select", False)
                with st.expander("📋 Training Set Snapshots Table", expanded=st.session_state.get("ml_all_tables_expanded", False)):
                    _edited_snaps = st.data_editor(
                        _df_snaps_display,
                        column_config={"Select": st.column_config.CheckboxColumn("Select", default=False)},
                        disabled=[c for c in _df_snaps_display.columns if c != "Select"],
                        hide_index=True,
                        use_container_width=True,
                        key="data_editor_snaps",
                    )
                _selected_snaps = _edited_snaps[_edited_snaps["Select"]]

                ## --- View set poses / Delete Selected ---
                snap_ids = df_snaps["training_set_id"].tolist()
                selected_ts = st.selectbox("Select a training set to inspect:", snap_ids, key="ts_viewer_select")

                ts_viewer_key = "show_ts_set_viewer"
                if ts_viewer_key not in st.session_state:
                    st.session_state[ts_viewer_key] = False
                ts_inspector_key = "show_ts_snapshot_inspector"
                if ts_inspector_key not in st.session_state:
                    st.session_state[ts_inspector_key] = False
                _col_ts_view, _col_ts_restore, _col_ts_inspect, _col_ts_delete = st.columns([3, 3, 3, 2])
                with _col_ts_view:
                    st.button(
                        f"{'Hide' if st.session_state[ts_viewer_key] else 'View'} Set Poses",
                        key="btn_ts_set_viewer",
                        on_click=_toggle_ts_set_viewer
                    )
                with _col_ts_restore:
                    _restore_ts = _selected_snaps["training_set_id"].iloc[0] if len(_selected_snaps) == 1 else selected_ts
                    if st.button("↩️ Restore to Binders", key="btn_restore_to_binders"):
                        _restore_result = st_funcs.restore_binders_from_snapshot(project_path, _restore_ts)
                        if isinstance(_restore_result, str) and _restore_result.startswith("error:"):
                            st.error(f"❌ {_restore_result[6:]}")
                        else:
                            _pos_total = _restore_result['pos_inserted'] + _restore_result['pos_skipped']
                            _neg_total = _restore_result['neg_inserted'] + _restore_result['neg_skipped']
                            st.success(
                                f"Restored from **{_restore_ts}**: "
                                f"✅ {_pos_total} positive "
                                f"({_restore_result['pos_inserted']} new, {_restore_result['pos_skipped']} already existed), "
                                f"❌ {_neg_total} negative "
                                f"({_restore_result['neg_inserted']} new, {_restore_result['neg_skipped']} already existed)."
                            )
                            st.rerun()
                with _col_ts_inspect:
                    st.button(
                        f"{'Hide' if st.session_state[ts_inspector_key] else 'Inspect'} Set Snapshot",
                        key="btn_ts_snapshot_inspector",
                        on_click=_toggle_ts_snapshot_inspector
                    )
                with _col_ts_delete:
                    _n_sel_snaps = len(_selected_snaps)
                    if not st.session_state.get("confirm_delete_selected_snaps"):
                        if st.button(
                            f"🗑️ Delete Selected ({_n_sel_snaps})",
                            key="btn_delete_selected_snaps",
                            disabled=(_n_sel_snaps == 0),
                        ):
                            st.session_state["confirm_delete_selected_snaps"] = True
                            st.rerun()
                    else:
                        st.warning(f"⚠️ Delete **{_n_sel_snaps}** snapshot(s) and all associated data? Are you sure?")
                        _cdsnap1, _cdsnap2 = st.columns(2)
                        with _cdsnap1:
                            if st.button("✅ Yes, delete", key="btn_confirm_delete_selected_snaps"):
                                _ids_to_delete = _selected_snaps["training_set_id"].tolist()
                                result = st_funcs.delete_training_set_snapshots(project_path, _ids_to_delete)
                                st.session_state["confirm_delete_selected_snaps"] = False
                                if result == "deleted":
                                    st.session_state[ts_viewer_key] = False
                                    st.rerun()
                                else:
                                    st.error(f"Could not delete snapshots: {result}")
                        with _cdsnap2:
                            if st.button("❌ Cancel", key="btn_cancel_delete_selected_snaps"):
                                st.session_state["confirm_delete_selected_snaps"] = False
                                st.rerun()

                ## --- Export Fingerprints CSV ---
                if len(_selected_snaps) == 0:
                    st.button(
                        "⬇️ Export Fingerprints CSV",
                        key="btn_export_fps_csv",
                        disabled=True,
                        help="Select a training set snapshot using the checkbox to enable export.",
                    )
                else:
                    _export_ts_id = _selected_snaps["training_set_id"].iloc[0]
                    _cond_options = st_funcs.get_training_set_prolif_condition_options(project_path, _export_ts_id)
                    if not _cond_options:
                        st.button(
                            "⬇️ Export Fingerprints CSV",
                            key="btn_export_fps_csv",
                            disabled=True,
                            help=f"No fingerprints computed for '{_export_ts_id}' yet. Run compute_training_set_fingerprints() first.",
                        )
                    else:
                        _cond_labels = {
                            (f"ID {c['id']}: {c['description']}" if c['description'] else f"ID {c['id']}"): c['id']
                            for c in _cond_options
                        }
                        _cond_label_sel = st.selectbox(
                            "ProLIF conditions", list(_cond_labels.keys()),
                            key=f"export_fps_cond_{_export_ts_id}",
                        )
                        _export_cond_id = _cond_labels[_cond_label_sel]
                        _csv_bytes = st_funcs.export_training_set_fingerprints_as_csv_bytes(project_path, _export_ts_id, _export_cond_id)
                        _safe_ts_id = _export_ts_id.replace(' ', '_')
                        _default_csv_path = os.path.join(
                            project_path, 'ml', 'training_sets',
                            f'{_safe_ts_id}_cond{_export_cond_id}_fingerprints.csv'
                        )
                        _export_path = st.text_input(
                            "Output CSV path:",
                            value=_default_csv_path,
                            key=f"export_fps_csv_path_{_export_ts_id}_{_export_cond_id}",
                        )
                        if st.button("⬇️ Export Fingerprints CSV", key="btn_export_fps_csv"):
                            try:
                                _out = _export_path.strip()
                                os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
                                with open(_out, 'wb') as _f:
                                    _f.write(_csv_bytes)
                                st.success(f"Fingerprints exported to: {_out}")
                            except Exception as _e:
                                st.error(f"Export failed: {_e}")

                ## --- Snapshot inspector panel ---
                if st.session_state[ts_inspector_key] and selected_ts:
                    df_inspect = st_funcs.get_training_set_entries(project_path, selected_ts)
                    st.markdown(f"##### 🔍 Snapshot: **{selected_ts}**")
                    if df_inspect is None or df_inspect.empty:
                        st.info("No entries found for this snapshot.")
                    else:
                        df_pos_inspect = df_inspect[df_inspect["binder_type"] == "positive"].drop(columns=["binder_type"], errors="ignore")
                        df_neg_inspect = df_inspect[df_inspect["binder_type"] == "negative"].drop(columns=["binder_type"], errors="ignore")
                        _icol1, _icol2 = st.columns(2)
                        with _icol1:
                            st.markdown(f"**✅ Positive binders ({len(df_pos_inspect)})**")
                            if df_pos_inspect.empty:
                                st.info("None.")
                            else:
                                st.dataframe(df_pos_inspect, use_container_width=True, hide_index=True)
                        with _icol2:
                            st.markdown(f"**❌ Negative binders ({len(df_neg_inspect)})**")
                            if df_neg_inspect.empty:
                                st.info("None.")
                            else:
                                st.dataframe(df_neg_inspect, use_container_width=True, hide_index=True)

                if st.session_state[ts_viewer_key] and selected_ts:
                    df_ts_entries = st_funcs.get_training_set_entries(project_path, selected_ts)
                    if df_ts_entries is None or df_ts_entries.empty:
                        st.info("No entries found for this training set.")
                    else:
                        ## Reference PDB uploader
                        ref_file_ts = st.file_uploader(
                            "Load a reference PDB for superimposition (optional):",
                            type=["pdb"],
                            key="ml_ts_ref_pdb_uploader"
                        )
                        if ref_file_ts is not None:
                            st.session_state["ml_ts_ref_pdb_data"] = ref_file_ts.read().decode("utf-8")
                            st.success(f"Reference loaded: {ref_file_ts.name}")
                        elif "ml_ts_ref_pdb_data" not in st.session_state:
                            st.session_state["ml_ts_ref_pdb_data"] = None

                        ## Build pose labels with binder type prefix
                        BINDER_ICON = {"positive": "✅ POSITIVE", "negative": "❌ NEGATIVE"}
                        ts_pose_labels = [
                            f"{BINDER_ICON.get(row['binder_type'], row['binder_type'])} — {row['assay_name']} / {row['pose_file']}"
                            for _, row in df_ts_entries.iterrows()
                        ]
                        ts_pose_paths = df_ts_entries["pose_full_path"].tolist()

                        _select_key_ts = "ml_ts_pose_select"
                        _idx_key_ts = "ml_ts_pose_idx"

                        # Initialise / clamp integer index
                        if _idx_key_ts not in st.session_state:
                            st.session_state[_idx_key_ts] = 0
                        if st.session_state[_idx_key_ts] >= len(ts_pose_labels):
                            st.session_state[_idx_key_ts] = 0

                        # Initialise selectbox only if absent or stale
                        if _select_key_ts not in st.session_state or st.session_state[_select_key_ts] not in ts_pose_labels:
                            st.session_state[_select_key_ts] = ts_pose_labels[st.session_state[_idx_key_ts]]

                        def _ts_go_prev():
                            if st.session_state[_idx_key_ts] > 0:
                                st.session_state[_idx_key_ts] -= 1
                                st.session_state[_select_key_ts] = ts_pose_labels[st.session_state[_idx_key_ts]]

                        def _ts_go_next():
                            if st.session_state[_idx_key_ts] < len(ts_pose_labels) - 1:
                                st.session_state[_idx_key_ts] += 1
                                st.session_state[_select_key_ts] = ts_pose_labels[st.session_state[_idx_key_ts]]

                        selected_ts_pose = st.selectbox(
                            "Select a pose:",
                            ts_pose_labels,
                            key=_select_key_ts,
                        )

                        # Sync integer index when user picks via dropdown
                        _dropdown_idx = ts_pose_labels.index(selected_ts_pose)
                        if _dropdown_idx != st.session_state[_idx_key_ts]:
                            st.session_state[_idx_key_ts] = _dropdown_idx

                        cur_ts = st.session_state[_idx_key_ts]
                        current_ts_row = df_ts_entries.iloc[cur_ts]
                        full_path_ts = ts_pose_paths[cur_ts]

                        ## Binder type badge
                        if current_ts_row["binder_type"] == "positive":
                            st.success("✅ Positive binder")
                        else:
                            st.error("❌ Negative binder")

                        if os.path.exists(full_path_ts):
                            with open(full_path_ts, "r") as f:
                                pdb_data_ts = f.read()
                            view_ts = py3Dmol.view(width=800, height=500)
                            view_ts.addModel(pdb_data_ts, "pdb")
                            view_ts.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                            if st.session_state.get("ml_ts_ref_pdb_data"):
                                view_ts.addModel(st.session_state["ml_ts_ref_pdb_data"], "pdb")
                                view_ts.setStyle({"model": 1}, {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6}, "stick": {"colorscheme": "elementColors", "opacity": 0.6}})
                            view_ts.zoomTo()

                            has_ref_ts = st.session_state.get("ml_ts_ref_pdb_data")

                            _ts_fps = st_funcs.get_training_set_fingerprints_for_pose(
                                project_path, selected_ts,
                                current_ts_row["assay_name"], current_ts_row["pose_file"],
                                current_ts_row["directory"]
                            )

                            # Pre-compute frequency data and apply threshold — shared by both columns
                            _ts_fp_meta = {}  # keyed by prolif_conditions_id
                            if _ts_fps:
                                for _fp_entry in _ts_fps:
                                    _cid = _fp_entry["prolif_conditions_id"]
                                    _freq_df_full = st_funcs.get_training_set_interaction_frequencies(
                                        project_path, selected_ts, _cid
                                    )
                                    _all_ints = st_funcs.get_all_training_set_interactions(
                                        project_path, selected_ts, _cid
                                    )
                                    _thresh_key = f"ml_ts_fp_thresh_{selected_ts}_{_cid}"
                                    if _thresh_key not in st.session_state:
                                        st.session_state[_thresh_key] = 0.0
                                    _ts_fp_meta[_cid] = {
                                        "freq_df_full": _freq_df_full,
                                        "all_interactions": _all_ints,
                                        "thresh_key": _thresh_key,
                                    }

                            if _ts_fps:
                                _viewer_area, _fp_area = st.columns([3, 2])
                            else:
                                _viewer_area = st.container()
                                _fp_area = None

                            with _viewer_area:
                                prev_col, viewer_col, next_col = st.columns([0.5, 8, 0.5 if not has_ref_ts else 1])
                                with prev_col:
                                    st.write(""); st.write(""); st.write("")
                                    st.button("◀", key="ml_ts_prev", disabled=(cur_ts == 0), on_click=_ts_go_prev)
                                with viewer_col:
                                    st.components.v1.html(view_ts.write_html(), height=520)
                                with next_col:
                                    st.write(""); st.write(""); st.write("")
                                    st.button("▶", key="ml_ts_next", disabled=(cur_ts >= len(ts_pose_labels) - 1), on_click=_ts_go_next)
                                st.caption(f"Pose {cur_ts + 1} of {len(ts_pose_labels)} — {current_ts_row['assay_name']} / {current_ts_row['pose_file']}")

                                ## Delete pose from Positive / Negative Binders registries
                                st.divider()
                                _df_pos_reg = st_funcs.get_binders_registry(project_path, "positive")
                                _df_neg_reg = st_funcs.get_binders_registry(project_path, "negative")

                                def _pose_in_registry(df, row):
                                    if df is None or df.empty:
                                        return False
                                    return (
                                        (df["assay_name"] == row["assay_name"]) &
                                        (df["pose_file"] == row["pose_file"]) &
                                        (df["directory"] == row["directory"])
                                    ).any()

                                _in_pos = _pose_in_registry(_df_pos_reg, current_ts_row)
                                _in_neg = _pose_in_registry(_df_neg_reg, current_ts_row)

                                if not _in_pos and not _in_neg:
                                    st.info(
                                        f"ℹ️ **{current_ts_row['pose_file']}** is not registered "
                                        f"in Positive or Negative Binders."
                                    )
                                else:
                                    _found_labels = []
                                    if _in_pos:
                                        _found_labels.append("✅ Positive Binders")
                                    if _in_neg:
                                        _found_labels.append("❌ Negative Binders")

                                    # Reflag target: only defined when pose is in exactly one registry
                                    if _in_pos and not _in_neg:
                                        _reflag_target = "negative"
                                        _reflag_label = "❌ Negative"
                                    elif _in_neg and not _in_pos:
                                        _reflag_target = "positive"
                                        _reflag_label = "✅ Positive"
                                    else:
                                        _reflag_target = None
                                        _reflag_label = None

                                    _del_binders_key = f"confirm_del_binders_{selected_ts}_{cur_ts}"
                                    _reflag_key = f"confirm_reflag_binders_{selected_ts}_{cur_ts}"

                                    _btn_col_del, _btn_col_ref = st.columns(2)

                                    with _btn_col_del:
                                        if not st.session_state.get(_del_binders_key):
                                            if st.button(
                                                f"🗑️ Delete from {' & '.join(_found_labels)}",
                                                key=f"btn_del_binders_{selected_ts}_{cur_ts}",
                                            ):
                                                st.session_state[_del_binders_key] = True
                                                st.session_state[_reflag_key] = False
                                                st.rerun()
                                        else:
                                            st.warning(
                                                f"Delete **{current_ts_row['pose_file']}** from "
                                                f"{', '.join(_found_labels)}?"
                                            )
                                            _dbc1, _dbc2 = st.columns(2)
                                            with _dbc1:
                                                if st.button("Yes, delete", key=f"btn_confirm_del_binders_{selected_ts}_{cur_ts}"):
                                                    if _in_pos:
                                                        st_funcs.remove_binder(
                                                            project_path=project_path,
                                                            binder_type="positive",
                                                            assay_name=current_ts_row["assay_name"],
                                                            pose_file=current_ts_row["pose_file"],
                                                            directory=current_ts_row["directory"],
                                                        )
                                                    if _in_neg:
                                                        st_funcs.remove_binder(
                                                            project_path=project_path,
                                                            binder_type="negative",
                                                            assay_name=current_ts_row["assay_name"],
                                                            pose_file=current_ts_row["pose_file"],
                                                            directory=current_ts_row["directory"],
                                                        )
                                                    st.session_state[_del_binders_key] = False
                                                    st.rerun()
                                            with _dbc2:
                                                if st.button("Cancel", key=f"btn_cancel_del_binders_{selected_ts}_{cur_ts}"):
                                                    st.session_state[_del_binders_key] = False
                                                    st.rerun()

                                    with _btn_col_ref:
                                        if _reflag_target is None:
                                            st.button(
                                                "🔄 Reflag",
                                                key=f"btn_reflag_{selected_ts}_{cur_ts}",
                                                disabled=True,
                                                help="Pose is in both registries — remove from one first.",
                                            )
                                        elif not st.session_state.get(_reflag_key):
                                            if st.button(
                                                f"🔄 Reflag as {_reflag_label}",
                                                key=f"btn_reflag_{selected_ts}_{cur_ts}",
                                            ):
                                                st.session_state[_reflag_key] = True
                                                st.session_state[_del_binders_key] = False
                                                st.rerun()
                                        else:
                                            st.warning(
                                                f"Move **{current_ts_row['pose_file']}** to {_reflag_label}?"
                                            )
                                            _rfc1, _rfc2 = st.columns(2)
                                            with _rfc1:
                                                if st.button("Yes, reflag", key=f"btn_confirm_reflag_{selected_ts}_{cur_ts}"):
                                                    _src = "positive" if _in_pos else "negative"
                                                    st_funcs.remove_binder(
                                                        project_path=project_path,
                                                        binder_type=_src,
                                                        assay_name=current_ts_row["assay_name"],
                                                        pose_file=current_ts_row["pose_file"],
                                                        directory=current_ts_row["directory"],
                                                    )
                                                    if _reflag_target == "positive":
                                                        st_funcs.save_positive_binder(
                                                            project_path=project_path,
                                                            assay_name=current_ts_row["assay_name"],
                                                            pose_file=current_ts_row["pose_file"],
                                                            directory=current_ts_row["directory"],
                                                            pose_full_path=current_ts_row["pose_full_path"],
                                                        )
                                                    else:
                                                        st_funcs.save_negative_binder(
                                                            project_path=project_path,
                                                            assay_name=current_ts_row["assay_name"],
                                                            pose_file=current_ts_row["pose_file"],
                                                            directory=current_ts_row["directory"],
                                                            pose_full_path=current_ts_row["pose_full_path"],
                                                        )
                                                    st.session_state[_reflag_key] = False
                                                    st.rerun()
                                            with _rfc2:
                                                if st.button("Cancel", key=f"btn_cancel_reflag_{selected_ts}_{cur_ts}"):
                                                    st.session_state[_reflag_key] = False
                                                    st.rerun()
                                st.divider()

                                ## VMD script creation
                                with st.expander("🎬 Create VMD Script", expanded=False):
                                    _default_vmd_path_ts = os.path.join(
                                        os.path.expanduser("~"), "Desktop",
                                        os.path.splitext(os.path.basename(full_path_ts))[0] + "_vmd.tcl"
                                    )
                                    _vmd_path_ts = st.text_input(
                                        "Save script to:",
                                        value=_default_vmd_path_ts,
                                        key=f"vmd_path_ml_ts_{cur_ts}",
                                    )
                                    if st.button("💾 Save VMD Script", key=f"btn_save_vmd_ml_ts_{cur_ts}"):
                                        try:
                                            _ref_pdb_path_ts = None
                                            if st.session_state.get("ml_ts_ref_pdb_data"):
                                                _ref_pdb_path_ts = os.path.join(
                                                    os.path.dirname(os.path.abspath(_vmd_path_ts.strip())),
                                                    "reference.pdb"
                                                )
                                                with open(_ref_pdb_path_ts, "w") as _rf:
                                                    _rf.write(st.session_state["ml_ts_ref_pdb_data"])
                                            _script_ts = st_funcs.generate_vmd_script(
                                                os.path.abspath(full_path_ts), _ref_pdb_path_ts
                                            )
                                            _vmd_out_ts = _vmd_path_ts.strip()
                                            os.makedirs(os.path.dirname(os.path.abspath(_vmd_out_ts)), exist_ok=True)
                                            with open(_vmd_out_ts, "w") as _sf:
                                                _sf.write(_script_ts)
                                            st.success(f"VMD script saved to: {_vmd_out_ts}")
                                            if _ref_pdb_path_ts:
                                                st.info(f"Reference PDB saved alongside: {_ref_pdb_path_ts}")
                                        except Exception as _e:
                                            st.error(f"Could not save VMD script: {_e}")

                                ## Interaction frequency table with threshold filter
                                if _ts_fps:
                                    for _fp_entry in _ts_fps:
                                        _cond_id = _fp_entry["prolif_conditions_id"]
                                        _meta = _ts_fp_meta[_cond_id]
                                        _freq_df_full = _meta["freq_df_full"]
                                        if _freq_df_full is not None:
                                            _thresh = st.number_input(
                                                "Min. frequency (%) to display:",
                                                min_value=0.0, max_value=100.0, step=5.0,
                                                key=_meta["thresh_key"],
                                            )
                                            _freq_df_shown = _freq_df_full[_freq_df_full["Frequency (%)"] >= _thresh]
                                            with st.expander(
                                                f"📊 Interaction frequencies — Conditions ID {_cond_id} "
                                                f"({len(_freq_df_shown)}/{len(_freq_df_full)} interactions, {len(df_ts_entries)} poses)",
                                                expanded=False
                                            ):
                                                st.dataframe(_freq_df_shown, use_container_width=True, hide_index=True)

                            if _fp_area is not None:
                                with _fp_area:
                                    st.markdown("##### 🔬 ProLIF Fingerprint")
                                    for _fp_entry in _ts_fps:
                                        _cond_id = _fp_entry["prolif_conditions_id"]
                                        _fp_dict = _fp_entry["fingerprint"]
                                        _active = [k for k, v in _fp_dict.items() if v]
                                        _meta = _ts_fp_meta[_cond_id]

                                        # Filter all_interactions by current frequency threshold
                                        _thresh = st.session_state.get(_meta["thresh_key"], 0.0)
                                        _freq_df_full = _meta["freq_df_full"]
                                        if _freq_df_full is not None and _thresh > 0:
                                            _above = set(_freq_df_full.loc[_freq_df_full["Frequency (%)"] >= _thresh, "Interaction"])
                                            _all_interactions = [i for i in _meta["all_interactions"] if i in _above]
                                        else:
                                            _all_interactions = _meta["all_interactions"]

                                        # Per-conditions filter state key
                                        _filter_key = f"ml_ts_fp_filter_{selected_ts}_{_cond_id}"
                                        if _filter_key not in st.session_state:
                                            st.session_state[_filter_key] = set(_all_interactions)

                                        with st.expander(f"Conditions ID {_cond_id} — filter interactions", expanded=False):
                                            _col_a, _col_b = st.columns(2)
                                            _all_selected = st.session_state[_filter_key] == set(_all_interactions)
                                            _none_selected = len(st.session_state[_filter_key]) == 0
                                            if _col_a.button("✅ All" if _all_selected else "☐ All",
                                                             key=f"ml_ts_fp_all_{selected_ts}_{_cond_id}",
                                                             type="primary" if _all_selected else "secondary"):
                                                st.session_state[_filter_key] = set(_all_interactions)
                                                for _iname in _all_interactions:
                                                    st.session_state[f"ml_ts_fp_cb_{selected_ts}_{_cond_id}_{_iname}"] = True
                                                st.rerun()
                                            if _col_b.button("✅ None" if _none_selected else "☐ None",
                                                             key=f"ml_ts_fp_none_{selected_ts}_{_cond_id}",
                                                             type="primary" if _none_selected else "secondary"):
                                                st.session_state[_filter_key] = set()
                                                for _iname in _all_interactions:
                                                    st.session_state[f"ml_ts_fp_cb_{selected_ts}_{_cond_id}_{_iname}"] = False
                                                st.rerun()
                                            for _iname in _all_interactions:
                                                _cb_key = f"ml_ts_fp_cb_{selected_ts}_{_cond_id}_{_iname}"
                                                _checked = st.checkbox(
                                                    _iname,
                                                    value=(_iname in st.session_state[_filter_key]),
                                                    key=_cb_key,
                                                )
                                                if _checked:
                                                    st.session_state[_filter_key].add(_iname)
                                                else:
                                                    st.session_state[_filter_key].discard(_iname)

                                        _visible = [k for k in _active if k in st.session_state[_filter_key]]
                                        _n_selected = len(st.session_state[_filter_key])
                                        st.markdown(f"**{len(_visible)} shown** ({_n_selected} selected / {len(_all_interactions)} total)")
                                        if _visible:
                                            st.dataframe(
                                                {"Interaction": _visible},
                                                use_container_width=True,
                                                hide_index=True,
                                            )
                                        elif _active:
                                            st.info("All active interactions are filtered out.")
                                        else:
                                            st.info("No interactions detected for this pose.")

                        else:
                            st.warning(f"Pose file not found on disk: {full_path_ts}")

elif page == "RF model training":
    st.title("RF model training")

    # ── Session state defaults ─────────────────────────────────────────────
    for _k in ("rf_df", "rf_model", "rf_eval_results", "rf_fi", "rf_cv_scores",
               "rf_best_params", "rf_best_cv_score", "rf_accuracy", "rf_macro_f1",
               "rf_test_predictions"):
        if _k not in st.session_state:
            st.session_state[_k] = None

    # ── Helpers for parsing custom hyperparameter values ──────────────────
    def _rf_parse_ints(text):
        result = []
        for token in text.split(","):
            token = token.strip()
            if token:
                try:
                    result.append(int(token))
                except ValueError:
                    st.warning(f"Ignored invalid integer value: '{token}'")
        return result

    def _rf_parse_depth(text):
        result = []
        for token in text.split(","):
            token = token.strip()
            if not token:
                continue
            if token.lower() == "none":
                result.append("None")
            else:
                try:
                    result.append(str(int(token)))
                except ValueError:
                    st.warning(f"Ignored invalid max_depth value: '{token}'")
        return result

    def _rf_parse_strs(text):
        return [t.strip() for t in text.split(",") if t.strip()]

    def _rf_combine(selected, extra):
        return list(dict.fromkeys(selected + extra))

    # ── Tabs ──────────────────────────────────────────────────────────────
    rf_tab_data, rf_tab_train, rf_tab_results, rf_tab_stored = st.tabs(["Data", "Training", "Results", "Stored models"])

    # ══ Tab 1 — Data ══════════════════════════════════════════════════════
    with rf_tab_data:
        st.subheader("Load Dataset")
        source = st.radio(
            "Data source",
            ["TidyScreen training set", "CSV file"],
            horizontal=True,
            key="rf_data_source",
        )

        if source == "TidyScreen training set":
            project_path = st.session_state.get("active_project_path")
            if not project_path:
                st.warning("No active project. Activate a project in the TidyScreen page first.")
            else:
                df_snaps = st_funcs.get_training_set_snapshots(project_path)
                if df_snaps is None or df_snaps.empty:
                    st.warning("No training set snapshots found for this project.")
                else:
                    snap_ids = df_snaps["training_set_id"].tolist()
                    selected_snap = st.selectbox("Training set snapshot", snap_ids, key="rf_snap_select")
                    fp_status = st_funcs.get_training_set_fingerprint_status(project_path)
                    st.caption(f"Fingerprint status: {fp_status.get(selected_snap, '❌ None')}")
                    _rf_cond_options = st_funcs.get_training_set_prolif_condition_options(project_path, selected_snap)
                    if not _rf_cond_options:
                        st.info("No fingerprints computed for this snapshot. Compute fingerprints first in the ML features management page.")
                    else:
                        _rf_cond_labels = {
                            (f"ID {c['id']}: {c['description']}" if c['description'] else f"ID {c['id']}"): c['id']
                            for c in _rf_cond_options
                        }
                        _rf_cond_label_sel = st.selectbox(
                            "ProLIF conditions", list(_rf_cond_labels.keys()),
                            key=f"rf_cond_select_{selected_snap}",
                        )
                        _rf_cond_id = _rf_cond_labels[_rf_cond_label_sel]
                        if st.button("Load fingerprints", key="rf_load_btn"):
                            csv_bytes = st_funcs.export_training_set_fingerprints_as_csv_bytes(project_path, selected_snap, _rf_cond_id)
                            if csv_bytes is None:
                                st.error("No fingerprints found for this snapshot/conditions combination.")
                            else:
                                st.session_state.rf_df = pd.read_csv(io.BytesIO(csv_bytes))
                                st.success(f"Loaded fingerprints from snapshot '{selected_snap}' (ProLIF conditions ID {_rf_cond_id}).")

        else:
            file_path = st.text_input("CSV file path", placeholder="/path/to/dataset.csv", key="rf_csv_path")
            uploaded = st.file_uploader("… or upload CSV", type="csv", key="rf_csv_upload")
            if file_path:
                try:
                    st.session_state.rf_df = pd.read_csv(file_path)
                except FileNotFoundError:
                    st.error(f"File not found: {file_path}")
                except Exception as exc:
                    st.error(f"Could not read file: {exc}")
            elif uploaded is not None:
                st.session_state.rf_df = pd.read_csv(uploaded)

        if st.session_state.rf_df is not None:
            _df = st.session_state.rf_df
            st.success(f"Dataset loaded — {_df.shape[0]:,} rows × {_df.shape[1]} columns")
            st.markdown("---")
            st.subheader("Dataset Overview")
            st.dataframe(_df, use_container_width=True)
            st.markdown("---")
            ov1, ov2 = st.columns(2)
            with ov1:
                st.markdown("**Class counts**")
                _label_preview = "label" if "label" in _df.columns else _df.columns[-1]
                counts = _df[_label_preview].value_counts().reset_index()
                counts.columns = ["Class", "Count"]
                st.bar_chart(counts.set_index("Class")["Count"])
                st.dataframe(counts, use_container_width=True, hide_index=True)
            with ov2:
                st.markdown("**Column info**")
                _info = pd.DataFrame({
                    "Column":   _df.columns,
                    "dtype":    _df.dtypes.astype(str).values,
                    "non-null": _df.notna().sum().values,
                    "null":     _df.isna().sum().values,
                })
                st.dataframe(_info, use_container_width=True, hide_index=True)

    # ══ Tab 2 — Training ══════════════════════════════════════════════════
    with rf_tab_train:
        _df_train = st.session_state.rf_df
        if _df_train is None:
            st.info("Please load a dataset in the **Data** tab first.")
        else:
            all_cols = _df_train.columns.tolist()

            st.subheader("Configuration")
            cfg1, cfg2 = st.columns(2)
            with cfg1:
                st.markdown("**Dataset columns**")
                _default_label_idx = all_cols.index("label") if "label" in all_cols else 0
                label_col = st.selectbox("Label column", all_cols, index=_default_label_idx, key="rf_label_col")
                _candidate_drop = [c for c in all_cols if c != label_col]
                _default_drop = [c for c in ("ligpose",) if c in _candidate_drop]
                drop_cols = st.multiselect("Extra columns to drop", _candidate_drop, default=_default_drop, key="rf_drop_cols")
            with cfg2:
                st.markdown("**Train / Test Split**")
                test_size = st.slider("Test set fraction", 0.10, 0.40, 0.20, 0.05, format="%.2f", key="rf_test_size")
                random_state = int(st.number_input("Random state", 0, 9999, 42, step=1, key="rf_random_state"))
                st.markdown("**Cross-Validation**")
                cv_folds = st.slider("Number of folds", 3, 10, 5, key="rf_cv_folds")
                st.markdown("**Feature Importances**")
                top_n = st.slider("Top N features", 5, 50, 20, key="rf_top_n")

            st.markdown("**Model Parameters**")
            use_gridsearch = st.checkbox("Hyperparameter tuning (GridSearchCV)", value=False, key="rf_use_gridsearch")

            if not use_gridsearch:
                bp1, bp2 = st.columns(2)
                with bp1:
                    n_estimators      = int(st.number_input("n_estimators",      min_value=10,  max_value=5000, value=200, step=10, key="rf_n_est"))
                    min_samples_split = int(st.number_input("min_samples_split", min_value=2,   max_value=50,   value=2,   step=1,  key="rf_min_split"))
                    min_samples_leaf  = int(st.number_input("min_samples_leaf",  min_value=1,   max_value=50,   value=1,   step=1,  key="rf_min_leaf"))
                with bp2:
                    max_features    = st.selectbox("max_features", ["sqrt", "log2"], key="rf_max_features")
                    unlimited_depth = st.checkbox("max_depth = None (unlimited)", value=True, key="rf_unlimited_depth")
                    max_depth = None if unlimited_depth else int(st.number_input("max_depth", min_value=1, max_value=200, value=10, step=1, key="rf_max_depth"))

            if use_gridsearch:
                with st.expander("Hyperparameter Grid", expanded=True):
                    st.caption("Type custom values (comma-separated) and press Enter — they are added as selected chips and the box is cleared.")
                    h1, h2 = st.columns(2)
                    with h1:
                        # ── n_estimators ──────────────────────────────────
                        _n_est_extra = st.session_state.get("rf_gs_n_est_accum", [])
                        _n_est_opts  = sorted(set([100, 200, 500] + _n_est_extra))
                        _n_est_prev  = sorted(set(v for v in st.session_state.get("rf_gs_n_est", [100, 200, 500]) if v in _n_est_opts))
                        st.session_state["rf_gs_n_est"] = sorted(set(_n_est_prev + _n_est_extra))
                        _n_est_sel   = st.multiselect("n_estimators", _n_est_opts, key="rf_gs_n_est")
                        st.text_input("Add n_estimators", placeholder="e.g. 300, 1000", key="rf_gs_n_est_extra", on_change=_gs_cb_n_est)

                        # ── max_depth (numeric sort, "None" last) ─────────
                        _depth_extra = st.session_state.get("rf_gs_depth_accum", [])
                        _depth_all   = list(dict.fromkeys(["None", "10", "20"] + _depth_extra))
                        _depth_opts  = sorted((v for v in _depth_all if v != "None"), key=int) + (["None"] if "None" in _depth_all else [])
                        _depth_prev  = [v for v in st.session_state.get("rf_gs_depth", ["None", "10", "20"]) if v in _depth_opts]
                        _depth_merged = list(dict.fromkeys(_depth_prev + _depth_extra))
                        st.session_state["rf_gs_depth"] = sorted((v for v in _depth_merged if v != "None"), key=int) + (["None"] if "None" in _depth_merged else [])
                        _depth_sel   = st.multiselect("max_depth", _depth_opts, key="rf_gs_depth")
                        st.text_input("Add max_depth", placeholder="e.g. 5, 30, None", key="rf_gs_depth_extra", on_change=_gs_cb_depth)

                        # ── min_samples_split ─────────────────────────────
                        _split_extra = st.session_state.get("rf_gs_split_accum", [])
                        _split_opts  = sorted(set([2, 5, 10] + _split_extra))
                        _split_prev  = sorted(set(v for v in st.session_state.get("rf_gs_split", [2, 5, 10]) if v in _split_opts))
                        st.session_state["rf_gs_split"] = sorted(set(_split_prev + _split_extra))
                        _split_sel   = st.multiselect("min_samples_split", _split_opts, key="rf_gs_split")
                        st.text_input("Add min_samples_split", placeholder="e.g. 3, 7", key="rf_gs_split_extra", on_change=_gs_cb_split)

                    with h2:
                        # ── min_samples_leaf ──────────────────────────────
                        _leaf_extra  = st.session_state.get("rf_gs_leaf_accum", [])
                        _leaf_opts   = sorted(set([1, 2, 4] + _leaf_extra))
                        _leaf_prev   = sorted(set(v for v in st.session_state.get("rf_gs_leaf", [1, 2, 4]) if v in _leaf_opts))
                        st.session_state["rf_gs_leaf"] = sorted(set(_leaf_prev + _leaf_extra))
                        _leaf_sel    = st.multiselect("min_samples_leaf", _leaf_opts, key="rf_gs_leaf")
                        st.text_input("Add min_samples_leaf", placeholder="e.g. 3, 6", key="rf_gs_leaf_extra", on_change=_gs_cb_leaf)

                        # ── max_features (alphabetical) ───────────────────
                        _feat_extra  = st.session_state.get("rf_gs_feat_accum", [])
                        _feat_opts   = sorted(set(["sqrt", "log2"] + _feat_extra))
                        _feat_prev   = sorted(set(v for v in st.session_state.get("rf_gs_feat", ["sqrt", "log2"]) if v in _feat_opts))
                        st.session_state["rf_gs_feat"] = sorted(set(_feat_prev + _feat_extra))
                        _feat_sel    = st.multiselect("max_features", _feat_opts, key="rf_gs_feat")
                        st.text_input("Add max_features", placeholder="e.g. 0.5, 0.8", key="rf_gs_feat_extra", on_change=_gs_cb_feat)

            st.markdown("---")
            if st.button("Train Model", type="primary", key="rf_train_btn"):
                X, y = prepare_features(_df_train, label_col=label_col, drop_cols=drop_cols)
                X_train, X_test, y_train, y_test = split_data(X, y, test_size=test_size, random_state=random_state)

                if use_gridsearch:
                    _depth_vals = [None if v == "None" else int(v) for v in _depth_sel]
                    param_grid = {
                        "n_estimators":      _n_est_sel,
                        "max_depth":         _depth_vals,
                        "min_samples_split": _split_sel,
                        "min_samples_leaf":  _leaf_sel,
                        "max_features":      _feat_sel,
                    }
                    with st.spinner("Running GridSearchCV — this may take several minutes…"):
                        model, best_params, best_cv_score = tune_hyperparameters(
                            X_train, y_train, param_grid=param_grid, cv=cv_folds, random_state=random_state,
                        )
                    st.session_state.rf_cv_scores     = None
                    st.session_state.rf_best_params   = best_params
                    st.session_state.rf_best_cv_score = best_cv_score
                else:
                    with st.spinner("Training baseline model…"):
                        model, cv_scores = train_baseline(
                            X_train, y_train,
                            n_estimators=n_estimators,
                            max_depth=max_depth,
                            min_samples_split=min_samples_split,
                            min_samples_leaf=min_samples_leaf,
                            max_features=max_features,
                            cv=cv_folds,
                            random_state=random_state,
                        )
                    st.session_state.rf_cv_scores     = cv_scores
                    st.session_state.rf_best_params   = None
                    st.session_state.rf_best_cv_score = None

                st.session_state.rf_eval_results = evaluate_model(model, X_test, y_test)
                st.session_state.rf_fi           = feature_importances(model, X.columns.tolist(), top_n=top_n)
                st.session_state.rf_model        = model
                _trained_report = st.session_state.rf_eval_results["classification_report"]
                st.session_state.rf_accuracy = _trained_report["accuracy"]
                st.session_state.rf_macro_f1 = _trained_report["macro avg"]["f1-score"]

                # Per-test-pose predictions, so the Results tab can show which
                # poses fall in each confusion-matrix quadrant (TP/FP/FN/TN).
                # X_test keeps the original DataFrame index, so identity columns
                # (whatever the user dropped as non-feature columns) can be
                # looked back up from _df_train and stay row-aligned with
                # y_test/y_pred/y_prob (all in X_test's row order).
                _id_cols = [c for c in drop_cols if c in _df_train.columns]
                _test_df = (
                    _df_train.loc[X_test.index, _id_cols].reset_index(drop=True)
                    if _id_cols else pd.DataFrame(index=range(len(X_test)))
                )
                _test_df["true_label"] = y_test.to_numpy()
                _test_df["pred_label"] = st.session_state.rf_eval_results["y_pred"]
                _test_df["pred_prob"]  = st.session_state.rf_eval_results["y_prob"]
                _test_df["quadrant"] = np.select(
                    [
                        (_test_df["true_label"] == 1) & (_test_df["pred_label"] == 1),
                        (_test_df["true_label"] == 0) & (_test_df["pred_label"] == 1),
                        (_test_df["true_label"] == 1) & (_test_df["pred_label"] == 0),
                    ],
                    ["TP", "FP", "FN"],
                    default="TN",
                )
                st.session_state.rf_test_predictions = _test_df

                st.success("Model trained successfully! See the **Results** tab.")

    # ══ Tab 3 — Results ═══════════════════════════════════════════════════
    with rf_tab_results:
        if st.session_state.rf_model is None:
            st.info("Train a model in the **Training** tab first.")
        else:
            ev     = st.session_state.rf_eval_results
            fi     = st.session_state.rf_fi
            report = ev.get("classification_report")
            _accuracy = report["accuracy"] if report else st.session_state.get("rf_accuracy")
            _macro_f1 = report["macro avg"]["f1-score"] if report else st.session_state.get("rf_macro_f1")

            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Test ROC-AUC", f"{ev['roc_auc']:.4f}")
            m2.metric("Accuracy",     f"{_accuracy:.4f}" if _accuracy is not None else "—")
            m3.metric("Macro F1",     f"{_macro_f1:.4f}" if _macro_f1 is not None else "—")
            if st.session_state.rf_cv_scores is not None:
                _cv = st.session_state.rf_cv_scores
                m4.metric("Mean CV ROC-AUC", f"{_cv.mean():.4f} ± {_cv.std():.4f}")
            elif st.session_state.rf_best_cv_score is not None:
                m4.metric("Best CV ROC-AUC", f"{st.session_state.rf_best_cv_score:.4f}")

            st.markdown("---")

            if st.session_state.rf_cv_scores is not None:
                st.subheader("Cross-Validation ROC-AUC per Fold")
                cv_arr = st.session_state.rf_cv_scores
                cv_df  = pd.DataFrame({
                    "Fold":    [f"Fold {i + 1}" for i in range(len(cv_arr))],
                    "ROC-AUC": cv_arr,
                })
                fig, ax = plt.subplots(figsize=(7, 3))
                ax.bar(cv_df["Fold"], cv_df["ROC-AUC"], color="steelblue", edgecolor="white")
                ax.axhline(cv_arr.mean(), color="crimson", linestyle="--", linewidth=1.5,
                           label=f"Mean = {cv_arr.mean():.4f}")
                ax.set_ylim(max(0, cv_arr.min() - 0.05), min(1.0, cv_arr.max() + 0.05))
                ax.set_ylabel("ROC-AUC")
                ax.legend()
                ax.spines[["top", "right"]].set_visible(False)
                st.pyplot(fig)
                plt.close(fig)

            if st.session_state.rf_best_params is not None:
                st.subheader("Best Hyperparameters (GridSearchCV)")
                bp_df = pd.DataFrame(st.session_state.rf_best_params.items(), columns=["Parameter", "Value"])
                st.dataframe(bp_df, use_container_width=True, hide_index=True)

            st.markdown("---")

            col_cm, col_fi = st.columns(2)
            with col_cm:
                st.subheader("Confusion Matrix")
                cm = ev.get("confusion_matrix")
                if cm is None:
                    st.info("Not available for this stored model (saved before detailed results were tracked).")
                else:
                    fig, ax = plt.subplots(figsize=(4, 3))
                    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=False, linewidths=0.5, ax=ax)
                    ax.set_xlabel("Predicted label")
                    ax.set_ylabel("True label")
                    fig.tight_layout()
                    st.pyplot(fig)
                    plt.close(fig)
            with col_fi:
                if fi is None or fi.empty:
                    st.subheader("Feature Importances")
                    st.info("Not available for this stored model (saved before detailed results were tracked).")
                else:
                    st.subheader(f"Top {len(fi)} Feature Importances")
                    fig, ax = plt.subplots(figsize=(5, max(3, len(fi) * 0.28)))
                    colors = plt.cm.viridis_r(np.linspace(0.2, 0.85, len(fi)))
                    ax.barh(fi.index[::-1], fi.values[::-1], color=colors[::-1])
                    ax.set_xlabel("Mean decrease in impurity")
                    ax.spines[["top", "right"]].set_visible(False)
                    fig.tight_layout()
                    st.pyplot(fig)
                    plt.close(fig)

            with st.expander("🔍 Poses by Confusion Matrix Quadrant", expanded=False):
                _tp_df = st.session_state.get("rf_test_predictions")
                if _tp_df is None or _tp_df.empty:
                    st.info("Not available for this stored model (saved before detailed results were tracked).")
                else:
                    _quad_labels = {
                        "TP": "✅ True Positive",
                        "FP": "❌ False Positive",
                        "FN": "❌ False Negative",
                        "TN": "✅ True Negative",
                    }
                    _quad_counts = _tp_df["quadrant"].value_counts()
                    _qc1, _qc2, _qc3, _qc4 = st.columns(4)
                    for _qcol, _q in zip((_qc1, _qc2, _qc3, _qc4), ("TP", "FP", "FN", "TN")):
                        _qcol.metric(_quad_labels[_q], int(_quad_counts.get(_q, 0)))

                    _sel_quad = st.radio(
                        "Show test-set poses classified as:",
                        ["TP", "FP", "FN", "TN"],
                        format_func=lambda q: _quad_labels[q],
                        horizontal=True,
                        key="rf_quadrant_select",
                    )
                    _quad_df = _tp_df[_tp_df["quadrant"] == _sel_quad].reset_index(drop=True)
                    if _quad_df.empty:
                        st.info(f"No poses in the {_quad_labels[_sel_quad]} quadrant.")
                    else:
                        st.dataframe(_quad_df, use_container_width=True)
                        st.download_button(
                            f"Download {_sel_quad} poses (.csv)",
                            data=_quad_df.to_csv(index=False).encode(),
                            file_name=f"{_sel_quad.lower()}_poses.csv",
                            mime="text/csv",
                            key=f"rf_dl_quad_{_sel_quad}",
                        )

                        st.markdown("---")
                        st.subheader("Pose Viewer")
                        if "ligpose" not in _quad_df.columns:
                            st.caption(
                                "Pose viewer unavailable — no 'ligpose' identity column found "
                                "(only present when training data was loaded from a TidyScreen "
                                "training set snapshot)."
                            )
                        else:
                            _rfq_ref_sp, _rfq_ref_col = st.columns([1, 9])
                            with _rfq_ref_col:
                                _rfq_ref_file = st.file_uploader(
                                    "Load a reference PDB for superimposition (optional):",
                                    type=["pdb"],
                                    key="rfq_reference_pdb_uploader",
                                )
                                if _rfq_ref_file is not None:
                                    st.session_state["rfq_reference_pdb_data"] = _rfq_ref_file.read().decode("utf-8")
                                    st.success(f"Reference loaded: {_rfq_ref_file.name}")
                                elif "rfq_reference_pdb_data" not in st.session_state:
                                    st.session_state["rfq_reference_pdb_data"] = None

                            _rfq_project_path = st.session_state.get("active_project_path")
                            if not _rfq_project_path:
                                st.warning("No active project. Activate a project in the TidyScreen page first.")
                            else:
                                _rfq_labels = [
                                    f"{r['ligpose']}  (true={int(r['true_label'])}, "
                                    f"pred={int(r['pred_label'])}, prob={r['pred_prob']:.2f})"
                                    for _, r in _quad_df.iterrows()
                                ]
                                _rfq_key = f"rfq_pose_selectbox_{_sel_quad}"
                                if st.session_state.get(_rfq_key) not in _rfq_labels:
                                    st.session_state[_rfq_key] = _rfq_labels[0]

                                _rfq_cur_idx = _rfq_labels.index(st.session_state[_rfq_key])
                                _rfq_prev_col, _rfq_mid_col, _rfq_next_col = st.columns([1, 4, 1])
                                with _rfq_prev_col:
                                    if st.button("◀ Prev", key=f"rfq_prev_btn_{_sel_quad}",
                                                 disabled=(_rfq_cur_idx == 0)):
                                        st.session_state[_rfq_key] = _rfq_labels[_rfq_cur_idx - 1]
                                        st.rerun()
                                with _rfq_next_col:
                                    if st.button("Next ▶", key=f"rfq_next_btn_{_sel_quad}",
                                                 disabled=(_rfq_cur_idx == len(_rfq_labels) - 1)):
                                        st.session_state[_rfq_key] = _rfq_labels[_rfq_cur_idx + 1]
                                        st.rerun()
                                with _rfq_mid_col:
                                    st.caption(f"Pose {_rfq_cur_idx + 1} of {len(_rfq_labels)}")

                                _rfq_sel_label = st.selectbox(
                                    "Select a pose:", _rfq_labels, key=_rfq_key,
                                )
                                _rfq_row = _quad_df.iloc[_rfq_labels.index(_rfq_sel_label)]

                                _rfq_pdb_path, _rfq_assay, _rfq_ligname, _rfq_run = st_funcs.resolve_ligpose_to_pdb(
                                    _rfq_project_path, _rfq_row["ligpose"]
                                )
                                if not _rfq_pdb_path:
                                    st.warning(
                                        f"PDB not found for pose '{_rfq_row['ligpose']}'. "
                                        "Extract poses first."
                                    )
                                else:
                                    with open(_rfq_pdb_path, "r") as _f:
                                        _rfq_pdb = _f.read()
                                    _rfq_view = py3Dmol.view(width=800, height=500)
                                    _rfq_view.addModel(_rfq_pdb, "pdb")
                                    _rfq_view.setStyle({"model": 0}, {"stick": {}, "cartoon": {"color": "spectrum"}})
                                    if st.session_state.get("rfq_reference_pdb_data"):
                                        _rfq_view.addModel(st.session_state["rfq_reference_pdb_data"], "pdb")
                                        _rfq_view.setStyle(
                                            {"model": 1},
                                            {"sphere": {"colorscheme": "elementColors", "scale": 0.3, "opacity": 0.6},
                                             "stick": {"colorscheme": "elementColors", "opacity": 0.6}},
                                        )
                                    _rfq_view.zoomTo()
                                    st.components.v1.html(_rfq_view.write_html(), height=510)
                                    st.caption(
                                        f"Assay: {_rfq_assay}  |  Ligand: {_rfq_ligname}  |  Run: {_rfq_run}  |  "
                                        f"True: {int(_rfq_row['true_label'])}  |  "
                                        f"Pred: {int(_rfq_row['pred_label'])}  |  "
                                        f"Prob: {_rfq_row['pred_prob']:.3f}"
                                    )

            st.subheader("Classification Report")
            if report is None:
                st.info("Not available for this stored model (saved before detailed results were tracked).")
            else:
                _rows = [{"class": cls, **vals} for cls, vals in report.items() if isinstance(vals, dict)]
                rep_df = (
                    pd.DataFrame(_rows)
                    .set_index("class")
                    .rename(columns={"precision": "Precision", "recall": "Recall",
                                     "f1-score": "F1-score", "support": "Support"})
                )
                rep_df["Support"] = rep_df["Support"].astype(int)
                st.dataframe(
                    rep_df.style.format({"Precision": "{:.4f}", "Recall": "{:.4f}", "F1-score": "{:.4f}"}),
                    use_container_width=True,
                )

            st.subheader("Feature Importance Table")
            if fi is None or fi.empty:
                st.info("Not available for this stored model (saved before detailed results were tracked).")
            else:
                fi_df = fi.reset_index()
                fi_df.columns = ["Feature", "Importance"]
                fi_df.insert(0, "Rank", range(1, len(fi_df) + 1))
                st.dataframe(
                    fi_df.style.format({"Importance": "{:.6f}"}),
                    use_container_width=True, hide_index=True, height=400,
                )
                st.download_button(
                    "Download feature importances (.csv)",
                    data=fi_df.to_csv(index=False).encode(),
                    file_name="feature_importances.csv",
                    mime="text/csv",
                    key="rf_dl_fi",
                )

            st.markdown("---")
            st.subheader("Export Trained Model")
            _buf = io.BytesIO()
            joblib.dump(st.session_state.rf_model, _buf)
            st.download_button(
                label="Download model (.joblib)",
                data=_buf.getvalue(),
                file_name="rf_model.joblib",
                mime="application/octet-stream",
                key="rf_dl_model",
            )

            st.markdown("---")
            st.subheader("Save Model to Project Database")
            _save_project_path = st.session_state.get("active_project_path")
            if not _save_project_path:
                st.warning("No active project. Activate a project in the TidyScreen page first.")
            else:
                _rf_snap_ref = st.session_state.get("rf_snap_select")
                if _rf_snap_ref:
                    st.caption(f"Training set snapshot linked to this model: **{_rf_snap_ref}**")
                else:
                    st.caption("Training set snapshot: *none* (model was trained from a CSV file)")

                _sv1, _sv2 = st.columns(2)
                with _sv1:
                    _rf_save_name = st.text_input(
                        "Model name", placeholder="e.g. rf_v1_baseline", key="rf_save_name"
                    )
                with _sv2:
                    _rf_save_desc = st.text_input(
                        "Description (optional)", placeholder="Short note about this model", key="rf_save_desc"
                    )

                if st.button("Save model (.pkl) to project DB", type="primary", key="rf_save_db_btn"):
                    if not _rf_save_name.strip():
                        st.error("Please enter a model name before saving.")
                    else:
                        try:
                            import pickle, datetime
                            _pkl_bytes = pickle.dumps(st.session_state.rf_model)
                            _ev = st.session_state.rf_eval_results
                            _rpt = _ev["classification_report"]
                            _roc = float(_ev["roc_auc"])
                            _acc = float(_rpt["accuracy"])
                            _f1  = float(_rpt["macro avg"]["f1-score"])
                            _cv_mean = float(st.session_state.rf_cv_scores.mean()) if st.session_state.rf_cv_scores is not None else None
                            _cv_std  = float(st.session_state.rf_cv_scores.std())  if st.session_state.rf_cv_scores is not None else None

                            # Full Results-tab data, serialised as JSON so it can be
                            # reconstructed later via the "Stored models" load button.
                            _cm_json = json.dumps(_ev["confusion_matrix"].tolist())
                            _rpt_json = json.dumps(_rpt)
                            _fi_json = json.dumps(st.session_state.rf_fi.to_dict())
                            _cv_json = (
                                json.dumps(st.session_state.rf_cv_scores.tolist())
                                if st.session_state.rf_cv_scores is not None else None
                            )
                            _bp_json = (
                                json.dumps(st.session_state.rf_best_params)
                                if st.session_state.rf_best_params is not None else None
                            )
                            _best_cv = (
                                float(st.session_state.rf_best_cv_score)
                                if st.session_state.rf_best_cv_score is not None else None
                            )
                            _tp_json = (
                                st.session_state.rf_test_predictions.to_json(orient="records")
                                if st.session_state.get("rf_test_predictions") is not None else None
                            )

                            _models_dir = os.path.join(_save_project_path, "ml", "models")
                            os.makedirs(_models_dir, exist_ok=True)
                            _db_path = os.path.join(_models_dir, "rf_trained_models.db")

                            st_funcs.ensure_rf_trained_models_schema(_db_path)

                            _conn = sqlite3.connect(_db_path)
                            _cur  = _conn.cursor()
                            _cur.execute(
                                """INSERT INTO rf_trained_models
                                   (model_name, description, training_set_id, roc_auc, accuracy, macro_f1,
                                    cv_roc_mean, cv_roc_std, model_pkl, created_at,
                                    confusion_matrix_json, classification_report_json,
                                    feature_importances_json, cv_scores_json,
                                    best_params_json, best_cv_score, test_predictions_json)
                                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                                (
                                    _rf_save_name.strip(),
                                    _rf_save_desc.strip(),
                                    _rf_snap_ref,
                                    _roc, _acc, _f1,
                                    _cv_mean, _cv_std,
                                    _pkl_bytes,
                                    datetime.datetime.now().isoformat(timespec="seconds"),
                                    _cm_json, _rpt_json, _fi_json, _cv_json,
                                    _bp_json, _best_cv, _tp_json,
                                ),
                            )
                            _conn.commit()
                            _new_id = _cur.lastrowid
                            _conn.close()
                            st.success(
                                f"Model **'{_rf_save_name.strip()}'** saved as ID {_new_id} "
                                f"in `{os.path.relpath(_db_path, _save_project_path)}`."
                            )
                        except Exception as _exc:
                            st.error(f"Failed to save model: {_exc}")

    # ══ Tab 4 — Stored models ═════════════════════════════════════════════
    with rf_tab_stored:
        st.subheader("Stored RF Models")

        if "rf_stored_load_message" in st.session_state:
            _msg_kind, _msg_text = st.session_state.pop("rf_stored_load_message")
            getattr(st, _msg_kind)(_msg_text)

        _stored_project_path = st.session_state.get("active_project_path")
        if not _stored_project_path:
            st.warning("No active project. Activate a project in the TidyScreen page first.")
        else:
            _stored_db_path = os.path.join(_stored_project_path, "ml", "models", "rf_trained_models.db")
            if not os.path.exists(_stored_db_path):
                st.info("No stored models found for this project yet. Train and save a model first.")
            else:
                try:
                    st_funcs.ensure_rf_trained_models_schema(_stored_db_path)
                    _sc = sqlite3.connect(_stored_db_path)
                    _sm_df = pd.read_sql_query(
                        """SELECT model_id, model_name, description, training_set_id,
                                  roc_auc, accuracy, macro_f1, cv_roc_mean, cv_roc_std,
                                  created_at, confusion_matrix_json
                           FROM rf_trained_models
                           ORDER BY model_id DESC""",
                        _sc,
                    )
                    _sc.close()

                    if _sm_df.empty:
                        st.info("No models stored yet.")
                    else:
                        st.dataframe(
                            _sm_df.drop(columns=["confusion_matrix_json"]).style.format({
                                "roc_auc":    "{:.4f}",
                                "accuracy":   "{:.4f}",
                                "macro_f1":   "{:.4f}",
                                "cv_roc_mean": lambda v: f"{v:.4f}" if pd.notna(v) else "—",
                                "cv_roc_std":  lambda v: f"{v:.4f}" if pd.notna(v) else "—",
                            }),
                            use_container_width=True,
                            hide_index=True,
                        )

                        st.markdown("---")
                        st.subheader("Load a Stored Model")
                        st.caption(
                            "Loads a stored model's data into the **Results** tab, exactly as it "
                            "appeared right after training."
                        )
                        _sm_labels = [
                            f"[{r.model_id}] {r.model_name}  ({r.created_at})"
                            + ("" if pd.notna(r.confusion_matrix_json) else "  — summary only")
                            for r in _sm_df.itertuples()
                        ]
                        _sm_sel_label = st.selectbox(
                            "Select a model to load:", _sm_labels, key="rf_stored_load_select",
                        )
                        _sm_sel_id = int(_sm_df.iloc[_sm_labels.index(_sm_sel_label)]["model_id"])
                        _sm_sel_name = str(_sm_df.iloc[_sm_labels.index(_sm_sel_label)]["model_name"])

                        _sm_load_col, _sm_del_col = st.columns(2)
                        with _sm_load_col:
                            _sm_load_clicked = st.button("📂 Load into Results tab", key="rf_stored_load_btn")
                        with _sm_del_col:
                            _sm_del_clicked = st.button("🗑️ Delete Model", key="rf_stored_delete_btn")

                        if _sm_load_clicked:
                            try:
                                _lc = sqlite3.connect(_stored_db_path)
                                _lrow = _lc.execute(
                                    """SELECT model_name, roc_auc, accuracy, macro_f1,
                                              model_pkl, confusion_matrix_json,
                                              classification_report_json, feature_importances_json,
                                              cv_scores_json, best_params_json, best_cv_score,
                                              test_predictions_json
                                       FROM rf_trained_models WHERE model_id = ?""",
                                    (_sm_sel_id,),
                                ).fetchone()
                                _lc.close()

                                if _lrow is None:
                                    st.session_state["rf_stored_load_message"] = (
                                        "error", "Selected model could not be found (it may have been deleted)."
                                    )
                                else:
                                    (_l_name, _l_roc, _l_acc, _l_f1, _l_pkl,
                                     _l_cm_json, _l_rep_json, _l_fi_json,
                                     _l_cv_json, _l_bp_json, _l_best_cv, _l_tp_json) = _lrow

                                    import pickle as _pickle
                                    _l_model = _pickle.loads(_l_pkl)
                                    _l_cm = np.array(json.loads(_l_cm_json)) if _l_cm_json else None
                                    _l_report = json.loads(_l_rep_json) if _l_rep_json else None
                                    _l_fi_dict = json.loads(_l_fi_json) if _l_fi_json else None
                                    _l_fi = pd.Series(_l_fi_dict, dtype=float) if _l_fi_dict else pd.Series(dtype=float)
                                    _l_cv_scores = np.array(json.loads(_l_cv_json)) if _l_cv_json else None
                                    _l_best_params = json.loads(_l_bp_json) if _l_bp_json else None
                                    _l_test_preds = (
                                        pd.read_json(io.StringIO(_l_tp_json), orient="records")
                                        if _l_tp_json else None
                                    )

                                    st.session_state.rf_model = _l_model
                                    st.session_state.rf_eval_results = {
                                        "confusion_matrix": _l_cm,
                                        "classification_report": _l_report,
                                        "roc_auc": _l_roc,
                                    }
                                    st.session_state.rf_fi = _l_fi
                                    st.session_state.rf_cv_scores = _l_cv_scores
                                    st.session_state.rf_best_params = _l_best_params
                                    st.session_state.rf_best_cv_score = (
                                        float(_l_best_cv) if _l_best_cv is not None else None
                                    )
                                    st.session_state.rf_accuracy = _l_acc
                                    st.session_state.rf_macro_f1 = _l_f1
                                    st.session_state.rf_test_predictions = _l_test_preds

                                    if _l_cm_json is None:
                                        st.session_state["rf_stored_load_message"] = (
                                            "warning",
                                            f"Loaded **{_l_name}**, but it was saved before detailed "
                                            "results (confusion matrix, classification report, feature "
                                            "importances, per-fold CV scores) were tracked — only the "
                                            "summary metrics are available. Retrain and re-save to "
                                            "capture full details.",
                                        )
                                    else:
                                        st.session_state["rf_stored_load_message"] = (
                                            "success",
                                            f"Loaded **{_l_name}**. Switch to the **Results** tab to view it.",
                                        )
                                st.rerun()
                            except Exception as _exc:
                                st.error(f"Failed to load model: {_exc}")

                        _sm_del_confirm_key = f"confirm_delete_rf_model_{_sm_sel_id}"
                        if _sm_del_clicked:
                            st.session_state[_sm_del_confirm_key] = True
                            st.rerun()

                        if st.session_state.get(_sm_del_confirm_key):
                            st.warning(f"Delete model **{_sm_sel_name}** (ID {_sm_sel_id})? This cannot be undone.")
                            _sm_del_yes_col, _sm_del_no_col = st.columns(2)
                            with _sm_del_yes_col:
                                if st.button("Yes, delete", key=f"btn_confirm_delete_rf_{_sm_sel_id}"):
                                    try:
                                        _dc = sqlite3.connect(_stored_db_path)
                                        _dc.execute("DELETE FROM rf_trained_models WHERE model_id = ?", (_sm_sel_id,))
                                        _dc.commit()
                                        _dc.close()
                                        st.session_state[_sm_del_confirm_key] = False
                                        st.session_state["rf_stored_load_message"] = (
                                            "success", f"Deleted model **{_sm_sel_name}** (ID {_sm_sel_id})."
                                        )
                                        st.rerun()
                                    except Exception as _exc:
                                        st.error(f"Failed to delete model: {_exc}")
                            with _sm_del_no_col:
                                if st.button("Cancel", key=f"btn_cancel_delete_rf_{_sm_sel_id}"):
                                    st.session_state[_sm_del_confirm_key] = False
                                    st.rerun()

                except Exception as _exc:
                    st.error(f"Failed to load stored models: {_exc}")

elif page == "Mol Viewer":
    st.title("Molecular Viewer")

    def _mv_draw_molecule(smiles, smarts=None):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, 0, None
            highlight_atoms, highlight_bonds, match_count = [], [], 0
            if smarts:
                query = Chem.MolFromSmarts(smarts)
                if query is None:
                    return None, 0, "Invalid SMARTS pattern."
                matches = mol.GetSubstructMatches(query)
                match_count = len(matches)
                if match_count > 0:
                    atom_set, bond_set = set(), set()
                    for match in matches:
                        atom_set.update(match)
                        for qbond in query.GetBonds():
                            bond = mol.GetBondBetweenAtoms(match[qbond.GetBeginAtomIdx()], match[qbond.GetEndAtomIdx()])
                            if bond is not None:
                                bond_set.add(bond.GetIdx())
                    highlight_atoms = sorted(atom_set)
                    highlight_bonds = sorted(bond_set)
            img = Draw.MolToImage(mol, highlightAtoms=highlight_atoms, highlightBonds=highlight_bonds)
            return img, match_count, None
        except Exception as e:
            return None, 0, f"Error parsing input: {e}"

    def _mv_draw_3d(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            mol = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol) != 0:
                return "Embedding failed: RDKit could not generate 3D coordinates for this molecule."
            return Chem.MolToPDBBlock(mol)
        except Exception as e:
            return f"3D generation error: {e}"

    smiles_input = st.text_input("Enter a SMILES string:", key="mv_smiles")
    smarts_input = st.text_input("Optional SMARTS to highlight:", key="mv_smarts")

    _mv_descriptor_funcs = _get_descriptor_funcs()
    _mv_descriptor_names = list(_mv_descriptor_funcs.keys())
    mv_descriptor_name = st.selectbox("Select an RDKit descriptor to compute:", _mv_descriptor_names, key="mv_descriptor")

    if smiles_input:
        col1, col2 = st.columns(2)
        img, match_count, draw_error = _mv_draw_molecule(smiles_input, smarts_input.strip())
        pdb_block = _mv_draw_3d(smiles_input)
        with col1:
            st.subheader("2D Structure")
            if draw_error:
                st.error(draw_error)
            if img is not None:
                st.image(img, use_container_width=True)
                if smarts_input.strip():
                    if match_count > 0:
                        st.success(f"SMARTS matches found: {match_count}")
                    else:
                        st.warning("No SMARTS match found in this molecule.")

                mv_mol = Chem.MolFromSmiles(smiles_input)
                mv_descriptor_func = _mv_descriptor_funcs.get(mv_descriptor_name)
                if mv_mol is not None and mv_descriptor_func is not None:
                    try:
                        mv_descriptor_value = mv_descriptor_func(mv_mol)
                        st.metric(mv_descriptor_name, f"{mv_descriptor_value:.4g}")
                    except Exception as _mv_desc_error:
                        st.error(f"Error computing descriptor '{mv_descriptor_name}': {_mv_desc_error}")
            else:
                st.error("Invalid SMILES string.")
        with col2:
            st.subheader("3D Structure")
            if pdb_block is None:
                st.error("Unable to render 3D molecule. Try a different SMILES string.")
            elif isinstance(pdb_block, str) and (pdb_block.startswith("Embedding failed") or pdb_block.startswith("3D generation error")):
                st.error(pdb_block)
            else:
                view = py3Dmol.view(width=400, height=400)
                view.addModel(pdb_block, 'pdb')
                view.setStyle({'stick': {}})
                view.zoomTo()
                st.components.v1.html(view.write_html(), height=420)
    else:
        st.info("Enter a SMILES string above to visualize the molecule.")

elif page == "Reaction Viewer":
    st.title("Reaction Viewer")
    st.write(
        "Simulate a SMARTS-based reaction. The underlying mechanism (RDKit "
        "ReactionFromSmarts + RunReactants) is homologous to "
        "ChemSpace.apply_reaction_workflow(), applied here to one manually "
        "entered reactant (or two, for a bimolecular scheme) instead of a "
        "whole table of compounds."
    )

    rxv_scheme = st.radio(
        "Reaction scheme:",
        ["Unimolecular (1 reactant)", "Bimolecular (2 reactants)"],
        key="rxv_scheme",
    )
    rxv_is_bimolecular = rxv_scheme.startswith("Bimolecular")

    rxv_col1, rxv_col2 = st.columns(2)
    with rxv_col1:
        rxv_reactant1 = st.text_input("Reactant 1 SMILES:", key="rxv_reactant1")
    with rxv_col2:
        rxv_reactant2 = st.text_input(
            "Reactant 2 SMILES:",
            key="rxv_reactant2",
            disabled=not rxv_is_bimolecular,
            help="Only used for the bimolecular scheme.",
        )

    rxv_smarts = st.text_area(
        "Reaction SMARTS:",
        key="rxv_smarts",
        placeholder="e.g. [C:1](=[O:2])[OH].[N:3]>>[C:1](=[O:2])[N:3]",
    )

    if st.button("Run Reaction", key="btn_run_reaction"):
        rxv_errors = []
        if not rxv_reactant1.strip():
            rxv_errors.append("Reactant 1 SMILES is required.")
        if rxv_is_bimolecular and not rxv_reactant2.strip():
            rxv_errors.append("Reactant 2 SMILES is required for a bimolecular reaction.")
        if not rxv_smarts.strip():
            rxv_errors.append("Reaction SMARTS is required.")

        for _rxv_err in rxv_errors:
            st.error(f"❌ {_rxv_err}")

        if not rxv_errors:
            rxv_mol1 = Chem.MolFromSmiles(rxv_reactant1.strip())
            if rxv_mol1 is None:
                st.error(f"❌ Invalid SMILES for reactant 1: '{rxv_reactant1.strip()}'")

            rxv_mol2 = None
            if rxv_is_bimolecular:
                rxv_mol2 = Chem.MolFromSmiles(rxv_reactant2.strip())
                if rxv_mol2 is None:
                    st.error(f"❌ Invalid SMILES for reactant 2: '{rxv_reactant2.strip()}'")

            if rxv_mol1 is not None and (not rxv_is_bimolecular or rxv_mol2 is not None):
                # ReactionFromSmarts doesn't just return None on bad input — some
                # malformed SMARTS (e.g. multi-step "a>>b>>c") raise instead.
                rxv_rxn, rxv_smarts_error = None, None
                try:
                    rxv_rxn = AllChem.ReactionFromSmarts(rxv_smarts.strip())
                except Exception as e:
                    rxv_smarts_error = str(e)

                if rxv_rxn is None:
                    st.error(f"❌ Invalid reaction SMARTS: {rxv_smarts_error or repr(rxv_smarts.strip())}")
                else:
                    rxv_expected = rxv_rxn.GetNumReactantTemplates()
                    rxv_provided = 2 if rxv_is_bimolecular else 1
                    if rxv_expected != rxv_provided:
                        st.warning(
                            f"⚠️ Reaction SMARTS expects {rxv_expected} reactant(s), but "
                            f"{rxv_provided} were provided. Attempting to run anyway."
                        )

                    try:
                        rxv_reactants = (rxv_mol1, rxv_mol2) if rxv_is_bimolecular else (rxv_mol1,)
                        rxv_results = rxv_rxn.RunReactants(rxv_reactants)
                    except Exception as e:
                        st.error(f"❌ Error running reaction: {e}")
                        rxv_results = ()

                    rxv_product_smiles = []
                    for rxv_product_set in rxv_results:
                        for rxv_product_mol in rxv_product_set:
                            try:
                                Chem.SanitizeMol(rxv_product_mol)
                                rxv_product_smiles.append(Chem.MolToSmiles(rxv_product_mol))
                            except Exception:
                                continue

                    # De-duplicate while preserving order (a reaction can match the
                    # same site more than once and yield identical products)
                    rxv_seen = set()
                    rxv_unique_products = []
                    for _smi in rxv_product_smiles:
                        if _smi not in rxv_seen:
                            rxv_seen.add(_smi)
                            rxv_unique_products.append(_smi)

                    st.divider()
                    st.subheader("Reaction Scheme")

                    rxv_n_reactants = 2 if rxv_is_bimolecular else 1
                    if rxv_unique_products:
                        st.success(f"✅ Reaction successful: {len(rxv_unique_products)} distinct product(s) generated.")
                        rxv_col_weights = [2] * rxv_n_reactants + [1] + [2] * len(rxv_unique_products)
                    else:
                        st.warning("⚠️ No valid product was generated. Check the SMARTS pattern and reactant(s).")
                        rxv_col_weights = [2] * rxv_n_reactants

                    # Single row: reactant(s) -> arrow -> product(s), so the whole
                    # scheme reads left to right like a reaction diagram.
                    rxv_scheme_cols = st.columns(rxv_col_weights)
                    rxv_col_i = 0

                    with rxv_scheme_cols[rxv_col_i]:
                        st.image(Draw.MolToImage(rxv_mol1), caption="Reactant 1", use_container_width=True)
                    rxv_col_i += 1

                    if rxv_is_bimolecular:
                        with rxv_scheme_cols[rxv_col_i]:
                            st.image(Draw.MolToImage(rxv_mol2), caption="Reactant 2", use_container_width=True)
                        rxv_col_i += 1

                    if rxv_unique_products:
                        with rxv_scheme_cols[rxv_col_i]:
                            st.markdown(
                                "<div style='text-align:center; font-size:2em; padding-top:70px;'>&#8594;</div>",
                                unsafe_allow_html=True,
                            )
                        rxv_col_i += 1

                        for _i, _smi in enumerate(rxv_unique_products, 1):
                            _prod_mol = Chem.MolFromSmiles(_smi)
                            with rxv_scheme_cols[rxv_col_i]:
                                if _prod_mol is not None:
                                    st.image(Draw.MolToImage(_prod_mol), caption=f"Product {_i}", use_container_width=True)
                                else:
                                    st.error(f"Could not render product {_i}")
                                st.code(_smi, language=None)
                            rxv_col_i += 1

elif page == "Python API":
    st.title("Python API")

    _api_search = st.text_input(
        "Search methods:",
        key="python_api_search",
        placeholder="Filter by method name or description...",
    )
    _api_search_lower = _api_search.strip().lower()

    _selected_tidyscreen = _render_python_api_section(
        "TidyScreen Methods",
        tidyscreen,
        "TidyScreen",
        _api_search_lower,
        include_classes=True,
        exclude_names={"ActivateProject"},
    )

    st.divider()

    _selected_chemspace = _render_python_api_section("ChemSpace Methods", ChemSpace, "ChemSpace", _api_search_lower)

    st.divider()

    _selected_moldock = _render_python_api_section("MolDock Methods", MolDock, "MolDock", _api_search_lower)

    st.divider()

    _selected_moldyn = _render_python_api_section("MolDyn Methods", MolDyn, "MolDyn", _api_search_lower)

    st.divider()

    _selected_ml = _render_python_api_section("ML Methods", MachineLearning, "MachineLearning", _api_search_lower)

    st.divider()
    st.subheader("Build Workflow File")
    st.write(
        "Checked methods are chained (in commented form, ready to be filled in and "
        "uncommented) into a standalone workflow script, following the same pattern "
        "as a hand-written chaining script."
    )

    _api_selections = {
        "TidyScreen": _selected_tidyscreen,
        "ChemSpace": _selected_chemspace,
        "MolDock": _selected_moldock,
        "MolDyn": _selected_moldyn,
        "MachineLearning": _selected_ml,
    }
    _api_total_selected = sum(len(v) for v in _api_selections.values())
    st.caption(f"{_api_total_selected} method(s) selected across all sections.")

    if _api_total_selected == 0:
        st.info("Check methods above to include them in the generated workflow file.")
    else:
        _workflow_project_name = st.text_input(
            "Project name (used in tidyscreen.ActivateProject(...)):",
            value=st.session_state.get("selected_project", ""),
            key="python_api_workflow_project_name",
        )
        _workflow_code = _build_workflow_script(
            _workflow_project_name.strip() or "project_name", _api_selections
        )
        st.code(_workflow_code, language="python")

        _default_workflow_path = os.path.join(
            st.session_state.get("active_project_path", os.getcwd()),
            f"{(_workflow_project_name.strip() or 'workflow')}_workflow.py",
        )
        _workflow_out_path = st.text_input(
            "Output .py path:",
            value=_default_workflow_path,
            key="python_api_workflow_output_path",
        )
        if st.button("💾 Save Workflow File", key="btn_save_python_api_workflow"):
            try:
                _out = _workflow_out_path.strip()
                os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
                with open(_out, "w", encoding="utf-8") as _f:
                    _f.write(_workflow_code)
                st.success(f"✅ Workflow saved to: {_out}")
            except Exception as _e:
                st.error(f"❌ Save failed: {_e}")
