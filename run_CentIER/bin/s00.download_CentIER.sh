# Create centier enviroment and install ltr_retriever
conda create -n centier -c bioconda -c conda-forge ltr_retriever=3.0.0
conda activate centier
conda install python=3.8
pip install pyfastx==2.1.0 numpy pandas scipy

# install gt
cd ~/software

# if you have not download gt, please use this link to get it.
wget https://github.com/genometools/genometools/archive/refs/tags/v1.6.5.tar.gz

tar -zxf v1.6.5.tar.gz
cd genometools_1.6.5
make -j4 cairo=no


# install CentIER
cd ~/software
# if you have a dictoinary name software to install software. 
git clone  https://github.com/simon19891216/CentIER.git

echo 'export PATH=$PATH:~/software/genometools_1.6.5/bin' >> s00.load_software.sh
echo 'export PATH=$PATH:~/software/CentIER' >> s00.load_software.sh
