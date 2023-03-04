# create a snpEff database using a gff3 and genomic DNA fasta file... (note, the chromosome names must match in the 2 files)

# reference: # https://www.biostars.org/p/50963/


mkdir snpEff
cd snpEff

# download and unpress the snpEff
# 2022/11/2 version 5.1 is latest
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip 
gunzip snpEff_latest_core.zip 


# download and unpress the fasta and annotation
cd snpEff


## wget https://www.rosaceae.org/rosaceae_downloads/Prunus_persica/Prunus_persica-genome.v2.0.a1/genes/Prunus_persica_v2.0.a1.gene.gff3.gz 
## wget https://www.rosaceae.org/rosaceae_downloads/Prunus_persica/Prunus_persica-genome.v2.0.a1/assembly/Prunus_persica_v2.0.a1_scaffolds.fasta.gz

# set env
DBNAME=Pyr.cuiguan

# link or copy files to data directory
mkdir data
mkdir data/$DBNAME
mv ./GWHBAOS00000000.gff ./data/$DBNAME/genes.gff
mv ./GWHBAOS00000000.genome.fasta ./data/$DBNAME/sequences.fa


#Edit snpEff.config and insert your specific database information:
echo "$DBNAME.genome : $DBNAME" >> snpEff.config

# Build the database
# Note: Use -noCheckCds and -noCheckProtein , as the cds and proteins IDs didn't match with fasta.
/tools/java/jdk-15.0.1/bin/java -jar snpEff.jar build -gff3 -v -noCheckCds -noCheckProtein  $DBNAME

