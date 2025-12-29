import streamlit as st
#import tidyscreen
import io
import sys

st.set_page_config(page_title="TidyScreen App", layout="wide")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio(
    "Go to",
    ("TidyScreen", "ChemSpace", "MolDock", "Analysis")
)

if page == "TidyScreen":
    st.title("TidyScreen")
    st.write("Welcome to the TidyScreen main page.")

    # if st.button("Show Projects"):
    #     buffer = io.StringIO()
    #     sys_stdout = sys.stdout
    #     sys.stdout = buffer
    #     try:
    #         tidyscreen.projects()
       
    #     finally:
    #         sys.stdout = sys_stdout
    #     output = buffer.getvalue()
    #     st.text(output)

elif page == "ChemSpace":
    st.title("ChemSpace")
    st.write("Welcome to the ChemSpace page.")

elif page == "MolDock":
    st.title("MolDock")
    st.write("Welcome to the MolDock page.")

elif page == "Analysis":
    st.title("Analysis")
    st.write("Welcome to the Analysis page.")