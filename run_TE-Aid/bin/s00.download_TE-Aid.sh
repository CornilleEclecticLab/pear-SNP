#!/usr/bin/env bash

git clone https://github.com/clemgoub/TE-Aid.git

cd TE-Aid

# Set the conda channel priority as default, "flexible" 
conda config --describe channel_priority
conda config --set channel_priority flexible

# Install TE-Aid
mamba env create -f TE_AID.yml

# Set the conda channel priority as "strict", which is recommand by bioconda
conda config --set channel_priority strict

