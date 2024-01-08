BATCH="branch15.bathch16.pear_Dec2023_noclone"

BASE_PREFIX="pear_Dec2023_noclone.Combine_chr.geno20_maf005.anno.syno.thin8k"

PROJECT_DIR="/shared/ifbstor1/projects/pear_snp3/pear"

########################
# The following variables will be called by other scripts
########################
WORK_DIR="${PROJECT_DIR}/run_FastStructure"

BED_PREFIX="/${PROJECT_DIR}/run_Population_genetic_filter/output/\
${BATCH}/\
${BASE_PREFIX}"

OUTPUT_DIR="${WORK_DIR}/output/${BATCH}"

# the minimum and maximum K to run fastStructure
MIN_K=2
MAX_K=15

LOAD_PERL5="module load perl/5.26.2"

CLUMPAK_DIR="${WORK_DIR}/bin/CLUMPAK/26_03_2015_CLUMPAK/CLUMPAK.NOGIT."

# make sure the dependent perl5 lib is installed in the following path via cpanm or other methods
export PERL5LIB="/shared/ifbstor1/home/ynie/perl5/lib/perl5":$PERL5LIB
