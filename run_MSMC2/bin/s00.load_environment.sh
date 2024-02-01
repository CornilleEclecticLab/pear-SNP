source ~/.bashrc
mamba activate msmc2
export PATH=/shared/ifbstor1/projects/pear_snp3/pear/run_MSMC2_infer_population/bin/msmc-tools:$PATH
module load bcftools/1.14


# Set work path
WORK_BIN=$(pwd)
WORK_DIR=$(dirname ${WORK_BIN})
REF_PREFIX="GWHBAOS00000000.genome.genmap"
CHR="chr76"


# Set run name for output file
RUN_NAME="chr76_A301_P11-1"