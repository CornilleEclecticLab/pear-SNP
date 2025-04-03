# create a new environment and install needed pacakages from conda
conda create -y -n NLGenomeSweeper -c bioconda -c conda-forge \
    python=3.6 blast muscle=3.8.1551 samtools bedtools hmmer transdecoder openjdk
conda activate NLGenomeSweeper

# install interproscan with Panther
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.45-80.0/interproscan-5.45-80.0-64-bit.tar.gz
tar -xzf interproscan-5.45-80.0-64-bit.tar.gz
ln -s $(pwd)/interproscan-5.45-80.0/interproscan.sh $(dirname $(which python))/interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-14.1.tar.gz
tar -xzf panther-data-14.1.tar.gz -C interproscan-5.45-80.0/data/

# download NLGenomeSweeper and install
git clone https://github.com/ntoda03/NLGenomeSweeper.git
ln -s $(pwd)/NLGenomeSweeper/NLGenomeSweeper $(dirname $(which python))/NLGenomeSweeper
NLGenomeSweeper -h
