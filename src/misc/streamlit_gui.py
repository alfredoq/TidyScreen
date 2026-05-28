import streamlit as st
from tidyscreen import tidyscreen
import io
import sys
import os
import glob
import py3Dmol
import streamlit_functions as st_funcs

tidyscreen_package_path = sys.argv[1]


st.set_page_config(page_title="TidyScreen App", layout="wide")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio(
    "Go to",
    ("TidyScreen", "ChemSpace", "Receptors", "Docking Methods", "MolDock", "Analysis", "ML")
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
    if st.button("Show Projects"):
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

elif page == "ChemSpace":
    st.title("ChemSpace")
    st.write("Welcome to the ChemSpace page.")
    
    db_path = os.path.join(st.session_state["active_project_path"], "chemspace", "processed_data", "chemspace.db")
    df = st_funcs.get_tables_info(db_path)
    
    ## Create a button to show the project tables info DataFrame
    if st.button("Show ChemSpace Tables Info"):
        st.dataframe(df)

    ## Display full content of a selected table
    st.divider()
    st.subheader("Display Table")

    if "show_display_table" not in st.session_state:
        st.session_state["show_display_table"] = False

    if df is not None and not df.empty:
        if st.button("Display Table"):
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
                        st.dataframe(display_df, use_container_width=True)
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
        if st.button("Depict Table"):
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
            if st.button(button_label):
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
        if st.button("Export Table"):
            st.session_state["show_export_options"] = not st.session_state["show_export_options"]
            if not st.session_state["show_export_options"]:
                st.session_state["export_save_path"] = ""

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

                if st.button("Save Table"):
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
                    if st.button("Confirm Save"):
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

        if st.button("Export PDB Model", key="toggle_export_pdb"):
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

elif page == "Docking Methods":
    st.title("Docking Methods")
    st.write(f"Docking methods registered for project: {st.session_state.get('selected_project', 'Unknown')}.")

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
        ## Show/hide summary table
        if "show_docking_methods" not in st.session_state:
            st.session_state["show_docking_methods"] = False

        if st.button(f"{'Hide' if st.session_state['show_docking_methods'] else 'Show'} Docking Methods Table"):
            st.session_state["show_docking_methods"] = not st.session_state["show_docking_methods"]

        if st.session_state["show_docking_methods"]:
            st.dataframe(df_methods[["id", "method_name", "docking_engine", "description", "created_date"]])

        ## Select a method to inspect details
        method_names = df_methods["method_name"].dropna().unique().tolist()
        if "selected_docking_method" not in st.session_state and method_names:
            st.session_state["selected_docking_method"] = method_names[0]

        selected_method = st.selectbox(
            "Select a docking method to inspect:",
            method_names,
            key="select_docking_method",
            index=method_names.index(st.session_state.get("selected_docking_method", method_names[0])) if method_names else 0
        )
        if selected_method != st.session_state.get("selected_docking_method"):
            st.session_state["selected_docking_method"] = selected_method

        method_row = df_methods[df_methods["method_name"] == selected_method].iloc[0]

        st.markdown(f"### Method: `{method_row['method_name']}`")
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"**Docking Engine:** {method_row['docking_engine']}")
            st.markdown(f"**Description:** {method_row['description'] or 'N/A'}")
            st.markdown(f"**Created:** {method_row['created_date']}")
        with col2:
            st.markdown("**Docking Parameters:**")
            st.code(method_row["parameters"] or "N/A", language="json")
            st.markdown("**Ligand Preparation Parameters:**")
            st.code(method_row["ligand_prep_params"] or "N/A", language="json")

elif page == "MolDock":
    st.title("MolDock")
    st.write("Welcome to the MolDock page.")

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
            ## Button to show/hide docking results DataFrame
            if "show_results" not in st.session_state:
                st.session_state["show_results"] = False

            if st.button(f"{'Hide' if st.session_state['show_results'] else 'Show'} {st.session_state['selected_assay_name']} Docking Results"):
                st.session_state["show_results"] = not st.session_state["show_results"]

            if st.session_state["show_results"]:
                st.write(df_results)

            ## Single "View Poses" button that expands to show the four folder buttons
            extracted_poses = st_funcs.get_extracted_poses_info(results_db_path)
            any_active = any(e["active"] for e in extracted_poses)

            if not any_active:
                st.button("View Poses", disabled=True)
            else:
                if "show_view_poses" not in st.session_state:
                    st.session_state["show_view_poses"] = False

                if st.button(f"{'Hide' if st.session_state['show_view_poses'] else 'View'} Poses"):
                    st.session_state["show_view_poses"] = not st.session_state["show_view_poses"]

                if st.session_state["show_view_poses"]:
                    ## Reference PDB uploader (shared across all pose folders)
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
                            if st.button(f"{'Hide' if st.session_state[key] else 'Show'} {label}", key=key + "_btn"):
                                st.session_state[key] = not st.session_state[key]
                            if st.session_state[key]:
                                pdb_files = sorted(glob.glob(os.path.join(entry["path"], "*.pdb")))
                                pdb_names = [os.path.basename(f) for f in pdb_files]
                                idx_key = f"pose_idx_{entry['directory']}"
                                if idx_key not in st.session_state:
                                    st.session_state[idx_key] = 0
                                ## Clamp index to valid range (e.g. folder contents changed)
                                st.session_state[idx_key] = max(0, min(st.session_state[idx_key], len(pdb_names) - 1))

                                ## Selectbox with NO key — driven entirely by idx_key via index param
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
                                                    st.warning("Already flagged.")
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
                                                    st.warning("Already flagged.")
                                                else:
                                                    st.error(result)

                                    ## pose counter label
                                    st.caption(f"Pose {current_idx + 1} of {len(pdb_names)}")
                        else:
                            st.button(label, disabled=True, key=key + "_btn")

            if df_mmpbsa_poses_results is None or df_mmpbsa_poses_results.empty:
                st.button(f"No MMPBSA results available for {st.session_state['selected_assay_name']}", disabled=True)
            else:
                ## Button to show/hide MMPBSA results DataFrame
                if "show_mmpbsa_results" not in st.session_state:
                    st.session_state["show_mmpbsa_results"] = False

                if st.button(f"{'Hide' if st.session_state['show_mmpbsa_results'] else 'Show'} {st.session_state['selected_assay_name']} MMPBSA Results"):
                    st.session_state["show_mmpbsa_results"] = not st.session_state["show_mmpbsa_results"]

                if st.session_state["show_mmpbsa_results"]:
                    # Create a selection box for MMPBSA poses.
                    poses = df_mmpbsa_poses_results['pose_id'].dropna().tolist()
                    if poses:
                        default_pose = st.session_state.get('selected_pose_id', poses[0])
                        if default_pose not in poses:
                            default_pose = poses[0]

                        selected_pose = st.selectbox(
                            f"Select a docked pose: {st.session_state['selected_assay_name']}",
                            poses,
                            key="select_pose_id",
                            index=poses.index(default_pose)
                        )

                        # Always update session state and plot on selection change
                        if selected_pose != st.session_state.get('selected_pose_id', None):
                            st.session_state['selected_pose_id'] = selected_pose

                        st.success(f"Selected Pose: {selected_pose}")
                        
                        # Using the selected pose, retrieve de MMPBSA data from the database
                        
                        mmpbsa_data = st_funcs.get_mmpbsa_data_for_pose(results_db_path, selected_pose)
                        
                        # Display the retrieved MMPBSA data
                        if mmpbsa_data is not None and not mmpbsa_data.empty:
                            st.write(mmpbsa_data)
                        
                        
                        ## Add a button to show/hide the MMPBSA data for the selected pose
                        if "show_mmpbsa_summary_data" not in st.session_state:
                            st.session_state["show_mmpbsa_summary_data"] = False
                        
                        if st.button(f"{'Hide' if st.session_state['show_mmpbsa_summary_data'] else 'Show'} MMPBSA Summary Data for Pose {selected_pose}"):
                            st.session_state["show_mmpbsa_summary_data"] = not st.session_state["show_mmpbsa_summary_data"]
                        if st.session_state["show_mmpbsa_summary_data"]:
                            # Show sum of the energy components for the selected pose
                            if mmpbsa_data is not None and not mmpbsa_data.empty:
                                total_energy = mmpbsa_data['total'].sum()
                                gas_energy = mmpbsa_data['gas'].sum()
                                ele_energy = mmpbsa_data['ele'].sum()
                                vdw_energy = mmpbsa_data['vdw'].sum()
                                polar_solvation = mmpbsa_data['polar_solvation'].sum()
                                nonpolar_solvation = mmpbsa_data['nonpolar_solvation'].sum()
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

                        if st.button(
                            f"{'Hide' if st.session_state['show_mmpbsa_decomposition_data'] else 'Show'} "
                            f"MMPBSA Per-Residue Decomposition Data for Pose {selected_pose}"
                            ):
                            
                            st.session_state["show_mmpbsa_decomposition_data"] = not st.session_state["show_mmpbsa_decomposition_data"]

                        # IMPORTANT: keep this outside the button block
                        if st.session_state["show_mmpbsa_decomposition_data"]:
                            if "energy_threshold" not in st.session_state:
                                st.session_state["energy_threshold"] = -1.0

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
                            _spacer, col_profile = st.columns([1, 5])
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
                            _spacer, col_profile = st.columns([1, 5])
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
                            _spacer, col_profile = st.columns([1, 5])
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
                            _spacer, col_profile = st.columns([1, 5])
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
                            _spacer, col_profile = st.columns([1, 5])
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
                            _spacer, col_profile = st.columns([1, 5])
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
                        st.warning("No valid pose IDs found in MMPBSA results")

            # ProLIF fingerprints button
            prolif_tables = st_funcs.get_prolif_tables(results_db_path)
            if not prolif_tables:
                st.button("Show ProLIF Fingerprints", disabled=True)
            else:
                if "show_prolif" not in st.session_state:
                    st.session_state["show_prolif"] = False

                if st.button(f"{'Hide' if st.session_state['show_prolif'] else 'Show'} ProLIF Fingerprints"):
                    st.session_state["show_prolif"] = not st.session_state["show_prolif"]

                if st.session_state["show_prolif"]:
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

                        if st.button(f"{'Hide' if st.session_state['show_all_prolif'] else 'Show'} all fingerprints for {selected_table}"):
                            st.session_state["show_all_prolif"] = not st.session_state["show_all_prolif"]

                        if st.session_state["show_all_prolif"]:
                            all_fps_df = st_funcs.get_all_prolif_fingerprints(results_db_path, selected_table)
                            if all_fps_df is not None and not all_fps_df.empty:
                                st.dataframe(filter_fps(all_fps_df), use_container_width=True)
                            else:
                                st.warning("Could not retrieve fingerprints for all poses.")
                    else:
                        st.warning("No poses found in the selected ProLIF table.")

            # Create a button to show/hide the histogram of docking scores for the selected ligand in the selected assay
            if st.button("Show/Hide Histogram"):
                 st.session_state["show_histogram"] = not st.session_state.get("show_histogram", False)
            if st.session_state.get("show_histogram", False):

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
                
                
    
elif page == "ML":
    st.title("Machine Learning")

    if "active_project_path" not in st.session_state:
        st.warning("No active project. Please select a project first.")
    else:
        project_path = st.session_state["active_project_path"]

        st.subheader("Training Set Registries")

        ## --- Positive binders ---
        st.markdown("### ✅ Positive Binders")
        df_pos = st_funcs.get_binders_registry(project_path, "positive")
        if df_pos is None or df_pos.empty:
            st.info("No positive binders registered yet.")
        else:
            st.success(f"{len(df_pos)} positive binder(s) registered.")
            st.dataframe(df_pos, use_container_width=True)

        st.divider()

        ## --- Negative binders ---
        st.markdown("### ❌ Negative Binders")
        df_neg = st_funcs.get_binders_registry(project_path, "negative")
        if df_neg is None or df_neg.empty:
            st.info("No negative binders registered yet.")
        else:
            st.success(f"{len(df_neg)} negative binder(s) registered.")
            st.dataframe(df_neg, use_container_width=True)
