import streamlit as st
from tidyscreen import tidyscreen
import io
import sys
import os
import streamlit_functions as st_funcs

tidyscreen_package_path = sys.argv[1]


st.set_page_config(page_title="TidyScreen App", layout="wide")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio(
    "Go to",
    ("TidyScreen", "ChemSpace", "Receptors", "MolDock", "Analysis")
)

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

            if df_mmpbsa_poses_results is None or df_mmpbsa_poses_results.empty:
                st.warning("No MMPBSA results")
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
                        

                        
                        
                    else:
                        st.warning("No valid pose IDs found in MMPBSA results")

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
                
                
    