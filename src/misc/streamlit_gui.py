import streamlit as st
from tidyscreen import tidyscreen
import io
import sys
import os
import glob
import json
import sqlite3
import py3Dmol
import streamlit_functions as st_funcs

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


st.set_page_config(page_title="TidyScreen App", layout="wide")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio(
    "Go to",
    ("TidyScreen", "ChemSpace", "Receptors", "MolDock", "MolDyn", "ProLIF Conditions", "Analysis", "ML features management", "RF model training")
)

## Persistent sidebar info: active project and assay
st.sidebar.divider()
active_project = st.session_state.get("activate_project_selectbox") or st.session_state.get("selected_project")
active_assay = st.session_state.get("select_assay_name") or st.session_state.get("selected_assay_name")
st.sidebar.markdown("**Active project:**")
st.sidebar.info(active_project if active_project else "None selected")
st.sidebar.markdown("**Active docking assay:**")
st.sidebar.info(active_assay if active_assay else "None selected")

if page == "TidyScreen":
    st.title("TidyScreen")
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

elif page == "ChemSpace":
    st.title("ChemSpace")
    st.write("Welcome to the ChemSpace page.")
    
    db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
    df = st_funcs.get_tables_info(db_path)
    
    ## Create a button to show the project tables info DataFrame
    if "show_tables_info" not in st.session_state:
        st.session_state["show_tables_info"] = False
    if st.button(f"{'Hide' if st.session_state['show_tables_info'] else 'Show'} ChemSpace Tables Info", key="btn_chemspace_tables_info"):
        st.session_state["show_tables_info"] = not st.session_state["show_tables_info"]
    if st.session_state["show_tables_info"]:
        st.dataframe(df)

    ## Display full content of a selected table
    st.divider()
    st.subheader("Display Table")

    if "show_display_table" not in st.session_state:
        st.session_state["show_display_table"] = False

    if df is not None and not df.empty:
        if st.button(f"{'Hide' if st.session_state['show_display_table'] else 'Display'} Table", key="btn_display_table"):
            st.session_state["show_display_table"] = not st.session_state["show_display_table"]

        if st.session_state["show_display_table"]:
            display_table_names = df["table"].tolist()
            display_selected_table = st.selectbox(
                "Select a table to display:",
                display_table_names,
                key="display_table_select"
            )
            display_columns = st_funcs.get_table_columns(db_path, display_selected_table)
            if display_columns:
                st.markdown("**Select columns to display:**")
                selected_display_cols = [col for col in display_columns if st.checkbox(col, value=True, key=f"display_col_{display_selected_table}_{col}")]
                if selected_display_cols:
                    display_df = st_funcs.read_table_columns_as_dataframe(db_path, display_selected_table, selected_display_cols)
                    if display_df is not None and not display_df.empty:
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
                        _subset_col1, _subset_col2 = st.columns([3, 1])
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
                                f"Create Subset ({_n_sel_display} rows)",
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
                                    st.success(f"✅ Subset table **'{_new_table_name.strip()}'** created with {_n_sel_display} row(s).")
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
        if st.button(f"{'Hide' if st.session_state['show_depict_options'] else 'Depict'} Table", key="btn_depict_table"):
            st.session_state["show_depict_options"] = not st.session_state["show_depict_options"]
            if not st.session_state["show_depict_options"]:
                st.session_state["show_depiction"] = False
                st.session_state["depiction_images"] = []

        if st.session_state["show_depict_options"]:
            table_names = df["table"].tolist()
            selected_table = st.selectbox("Select a table to depict:", table_names, key="depict_table_select")

            ## Load columns for the selected table to let the user pick the label
            table_columns = st_funcs.get_table_columns(db_path, selected_table)
            default_label = "id" if "id" in table_columns else (table_columns[0] if table_columns else None)
            default_idx = table_columns.index(default_label) if default_label in table_columns else 0

            label_col = st.selectbox(
                "Column to use as molecule label:",
                table_columns,
                index=default_idx,
                key="depict_label_col"
            ) if table_columns else None

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
                else:
                    with st.spinner(f"Generating depictions for table '{selected_table}'…"):
                        try:
                            images = st_funcs.depict_table_to_images(
                                db_path=db_path,
                                table_name=selected_table,
                                max_molecules=int(max_mols),
                                molecules_per_image=int(mols_per_image),
                                mol_image_size=(int(mol_size), int(mol_size)),
                                legend_col=label_col
                            )
                            if images:
                                st.session_state["depiction_images"] = images
                                st.session_state["show_depiction"] = True
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
        if st.button(f"{'Hide' if st.session_state['show_export_options'] else 'Export'} Table", key="btn_export_table"):
            st.session_state["show_export_options"] = not st.session_state["show_export_options"]
            if not st.session_state["show_export_options"]:
                st.session_state["export_save_path"] = ""
                st.session_state["show_export_path_input"] = False

        if st.session_state["show_export_options"]:
            export_table_names = df["table"].tolist()
            export_selected_table = st.selectbox(
                "Select a table to export:",
                export_table_names,
                key="export_table_select"
            )

            export_columns = st_funcs.get_table_columns(db_path, export_selected_table)
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
                                export_df = st_funcs.read_table_columns_as_dataframe(db_path, export_selected_table, selected_cols)
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

elif page == "MolDock":
    st.title("MolDock")
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

        ## Inspector — always visible, collapsed expander above does not hide it
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
        ## Create a button to show/hide the docking assays DataFrame
        if "show_assays" not in st.session_state:
            st.session_state["show_assays"] = False

        if st.button(f"{'Hide' if st.session_state['show_assays'] else 'Show'} Docking Assays Details"):
            st.session_state["show_assays"] = not st.session_state["show_assays"]

        if st.session_state["show_assays"]:
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
                st.session_state["selected_assay_name"] = selected_assay
            st.success(f"Selected assay: {selected_assay}")

        st.write(f"Selected Docking Assay: {st.session_state.get('selected_assay_name', 'None')}")


elif page == "MolDyn":
    st.title("MolDyn")
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
            ## Show/hide summary table
            if "show_md_methods" not in st.session_state:
                st.session_state["show_md_methods"] = False

            if st.button(f"{'Hide' if st.session_state['show_md_methods'] else 'Show'} MD Methods Table"):
                st.session_state["show_md_methods"] = not st.session_state["show_md_methods"]
                if not st.session_state["show_md_methods"]:
                    st.session_state["export_mdm_selected_name"] = None
                st.rerun()

            if st.session_state["show_md_methods"]:
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

            ## Select a method to inspect details (visible only when table is shown)
            if st.session_state.get("show_md_methods"):
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
            if "show_md_assays" not in st.session_state:
                st.session_state["show_md_assays"] = False

            if st.button(f"{'Hide' if st.session_state['show_md_assays'] else 'Show'} MD Assay Registries"):
                st.session_state["show_md_assays"] = not st.session_state["show_md_assays"]

            if st.session_state["show_md_assays"]:
                st.dataframe(df_md_assays, use_container_width=True, hide_index=True)

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


elif page == "Analysis":
    st.title("Analysis")
    st.write("Welcome to the Analysis page.")
            
    if "selected_assay_name" not in st.session_state:
        st.warning("Select a docking assay")
    else:
        results_db_path = os.path.join(st.session_state["active_project_path"], "docking", "docking_assays", st.session_state["selected_assay_name"], "results", f"{st.session_state['selected_assay_name']}.db")

        df_results = st_funcs.get_docking_results(results_db_path)

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

            ## Single "View Poses" button that expands to show the four folder buttons
            extracted_poses = st_funcs.get_extracted_poses_info(results_db_path)
            any_active = any(e["active"] for e in extracted_poses)

            if not any_active:
                st.button("View Poses", disabled=True)
            else:
                with st.expander("🔍 View Poses", expanded=False):
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

                                ## ProLIF Fingerprints for this pose directory
                                dir_name = entry['directory']
                                dir_pose_ids = st_funcs.get_pose_ids_for_directory(results_db_path, dir_name)
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
                                        prolif_poses_dir_df = st_funcs.get_prolif_poses_by_pose_ids(results_db_path, selected_table_dir, dir_pose_ids)
                                        if prolif_poses_dir_df is not None and not prolif_poses_dir_df.empty:
                                            pose_options_dir = prolif_poses_dir_df["pose_id"].tolist()
                                            pose_label_map = st_funcs.get_pose_labels_for_pose_ids(results_db_path, pose_options_dir)
                                            pose_options_dir = sorted(pose_options_dir, key=lambda pid: pose_label_map.get(pid, str(pid)))
                                            selected_prolif_pose_dir = st.selectbox(
                                                "Select a pose to display fingerprints:",
                                                pose_options_dir,
                                                format_func=lambda pid: pose_label_map.get(pid, str(pid)),
                                                key=f"select_prolif_pose_{dir_name}"
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
                                                st.warning(f"No fingerprint data found for pose {selected_prolif_pose_dir}.")
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
                                        else:
                                            st.warning(f"No ProLIF poses found for {dir_name}.")
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
                        if 'selected_lig_name' not in st.session_state and lig_names:
                            st.session_state['selected_lig_name'] = lig_names[0]
                        selected_lig = st.selectbox(
                            f"Select a Ligand (LigName) in Assay: {st.session_state['selected_assay_name']}",
                            lig_names,
                            key="select_lig_name",
                            index=lig_names.index(st.session_state.get('selected_lig_name', lig_names[0])) if lig_names else 0
                        )
                        # Always update session state and plot on selection change
                        if selected_lig != st.session_state.get('selected_lig_name', None):
                            st.session_state['selected_lig_name'] = selected_lig
                        st.success(f"Selected Ligand: {selected_lig}")
                        histogram = st_funcs.construct_hist_for_ligand(df_results, selected_lig)
                        st.pyplot(histogram, use_container_width=False, clear_figure=True)
                
                
    
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
                    _col_pos_view, _col_pos_clear, _col_pos_delete = st.columns([3, 2, 2])
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
                    _col_neg_view, _col_neg_clear, _col_neg_delete = st.columns([3, 2, 2])
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
                    _csv_bytes = st_funcs.export_training_set_fingerprints_as_csv_bytes(project_path, _export_ts_id)
                    _safe_ts_id = _export_ts_id.replace(' ', '_')
                    if _csv_bytes is None:
                        st.button(
                            "⬇️ Export Fingerprints CSV",
                            key="btn_export_fps_csv",
                            disabled=True,
                            help=f"No fingerprints computed for '{_export_ts_id}' yet. Run compute_training_set_fingerprints() first.",
                        )
                    else:
                        _default_csv_path = os.path.join(
                            project_path, 'ml', 'training_sets',
                            f'{_safe_ts_id}_fingerprints.csv'
                        )
                        _export_path = st.text_input(
                            "Output CSV path:",
                            value=_default_csv_path,
                            key=f"export_fps_csv_path_{_export_ts_id}",
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
    st.info("RF model training functionality coming soon.")
