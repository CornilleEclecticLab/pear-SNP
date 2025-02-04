module load gcc/11.2.0
module load conda
source activate Fbranch

WORK_DIR=/shared/projects/pear_snp3/pear/run_Dsuite_genome/
export PATH="$WORK_DIR/bin/Dsuite.NOGIT./Build/":$PATH
export PATH="$WORK_DIR/bin/Dsuite.NOGIT./utils/":$PATH
