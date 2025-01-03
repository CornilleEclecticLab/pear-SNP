#!/usr/bin/bash

#SBATCH -e s00.download_database.%J.err
#SBATCH -o s00.download_database.%J.out
#SBATCH --mem=120G



module load eggnog-mapper/2.1.12

# Specific hmmer DB and full Diamond and full MMseqs2 DB
download_eggnog_data.py  -f -H -d 33090 --dbname Viridiplantae --data_dir ../database -y

# Specific Diamond or MMseqs2 DB
create_dbs.py --dbname Viridiplantae --taxids 33090 --data_dir ../database -y
