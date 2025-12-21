#!/bin/bash

# Ask user for environment name
read -p "Enter the name for the conda environment (default: tidyscreen): " ENV_NAME
ENV_NAME=${ENV_NAME:-tidyscreen}

# Check if the main TidyScreen environment already exists
if conda env list | grep -q "^$ENV_NAME "; then
    echo "Warning: Conda environment '$ENV_NAME' already exists."
    read -p "Do you want to replace it? (y/N): " REPLACE_ENV
    REPLACE_ENV=${REPLACE_ENV:-n}
    
    if [[ "$REPLACE_ENV" =~ ^[Yy]$ ]]; then
        echo "Removing existing environment: $ENV_NAME"
        echo "y" | conda env remove -n $ENV_NAME -y
    else
        echo "Installation cancelled. Environment '$ENV_NAME' already exists."
        exit 1
    fi
fi

# Check if the auxiliaty adt environment already exists
if conda env list | grep -q "adt"; then
    echo "Warning: Conda environment 'adt' already exists."
    read -p "Do you want to replace it? (y/N): " REPLACE_ENV
    REPLACE_ENV=${REPLACE_ENV:-n}
    
    if [[ "$REPLACE_ENV" =~ ^[Yy]$ ]]; then
        echo "Removing existing environment: adt"
        echo "y" | conda env remove -n adt -y
    else
        echo "Installation cancelled. Environment 'adt' already exists."
        exit 1
    fi
fi

echo "Creating conda environment: $ENV_NAME"
echo "y" | conda create -n $ENV_NAME python=3.12

echo "Installing TidyScreen"

conda run -n $ENV_NAME pip install git+https://github.com/alfredoq/TidyScreen

echo "Installing Dependencies"

echo "y" | conda install -n $ENV_NAME -c conda-forge ambertools==23.6 espaloma espaloma_charge chemicalite visidata vmd-python vina pdbfixer

conda run -n $ENV_NAME pip install git+https://github.com/forlilab/Meeko@develop

echo "y" | conda install -n $ENV_NAME -c bioconda autodock autogrid

echo "y" | conda install -n $ENV_NAME -c conda-forge redis-server

echo "Installing Ersilia Hub"

# Try to delete if already exists
if [ -d "/tmp/ersilia" ]; then
    rm -rf /tmp/ersilia
    echo "folder deleted"
fi

cd /tmp && git clone https://github.com/ersilia-os/ersilia.git

cd ersilia && conda run -n $ENV_NAME pip install -e . 

# Execute BentoML for the environment to retrieve the custom installation

conda run -n $ENV_NAME bentoml

# A second environment is needed for certain Autodock Tools functionalities requiring Python 2.7

echo "y" | conda create -n adt python=2.7

echo "y" | conda install -n adt -c insilichem autodocktools-prepare
