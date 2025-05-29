#!/bin/bash

mkdir -p ../database

cd ../database

wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz

gunzip uniref90.fasta.gz
