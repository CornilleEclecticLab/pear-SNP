# wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar


# Clone repository
git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git SIFT4G_Create_Genomic_DB.NOGIT.

# Cp def file
cp SIFT4G_Create_Genomic_DB.def ./SIFT4G_Create_Genomic_DB.NOGIT. 

# Build sif container
cd ./SIFT4G_Create_Genomic_DB.NOGIT.

apptainer build sift4g_db.sif SIFT4G_Create_Genomic_DB.def
