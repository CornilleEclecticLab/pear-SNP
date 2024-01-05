#!/bin/bash
# By Xilong CHEN
# Create date: 2023-10-24
# Contact: chen_xilong@outlook.com
#SBATCH -J s02.CLUMPAK
#SBATCH -o s02.CLUMPAK.out
#SBATCH -e s02.CLUMPAK.err

module purge
module load bioinfo/CLUMPAK-v1.1

cd /work/xchen/XCHEN_apple_genome/analysis_2022_Mar/run_ECRNA_FastStructure/output/s02.CLUMPAK

#Create symbolic link (or copy) mcl, distruct, fonts and CLUMPP directory in your working directory (adjust paramfile parameter if needed)
ln -s /usr/local/bioinfo/src/CLUMPAK/CLUMPAK-v1.1/26_03_2015_CLUMPAK/CLUMPAK/CLUMPP .
ln -s /usr/local/bioinfo/src/CLUMPAK/CLUMPAK-v1.1/26_03_2015_CLUMPAK/CLUMPAK/mcl .
ln -s /usr/local/bioinfo/src/CLUMPAK/CLUMPAK-v1.1/26_03_2015_CLUMPAK/CLUMPAK/distruct .
ln -s /usr/local/bioinfo/src/CLUMPAK/CLUMPAK-v1.1/26_03_2015_CLUMPAK/CLUMPAK/fonts .
ln ../s01.generate_fast/input_clumpak.zip
CLUMPAK.pl \
    --id 85018752 \
    --dir ./output_clumpak \
    --file input_clumpak.zip \
    --inputtype admixture

# copy summary pdf and zip

