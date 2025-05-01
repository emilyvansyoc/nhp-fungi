#!/bin/bash

# create python 2 environment to run Andrew's scripts in

# run step-by-step, not as a bash file

# create new conda environment
conda create -n py2 python=2.7

# activate
conda activate py2

# test python install and version
python --version

# install dependencies from Andrew's scripts
conda install cogent -y
conda install pandas -y
conda install scipy -y
conda install statsmodels -y
conda install seaborn -y
conda install biopython -y

# install ete2 from archive
# first have to install pip
python -m ensurepip --upgrade # this is from the internet, no idea what it does
pip install --upgrade ete2
