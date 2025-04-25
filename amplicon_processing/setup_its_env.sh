#!/bin/bash

### ----- SETUP CONDA ENVIRONMENT -----
# Create new conda environment
ENV_NAME="its_pipeline_env"
echo "Creating conda environment: $ENV_NAME"

# Create environment with python 3.8
conda create -n $ENV_NAME python=3.8 -y

# Activate the environment
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Add required channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install packages from software_list.txt
while read -r package; do
    if [ -n "$package" ]; then  # Skip empty lines
        echo "Installing $package..."
        conda install -y $package
    fi
done < "$(dirname "$0")/software_list.txt"

# Verify installation
echo "Verifying package installations..."
conda list

# Print success message
echo "Conda environment setup complete. You can now run the ITS pipeline."
echo "To activate the environment in the future, use:"
echo "conda activate $ENV_NAME" 