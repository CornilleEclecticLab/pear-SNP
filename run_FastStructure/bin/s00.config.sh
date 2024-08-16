##### All the samples
# s01 have been run
# BATCH="branch15.batch21.pear_Dec2023.noclone"
# BASE_PREFIX="pear_Dec2023.noclone.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 run
# s02 run
# BATCH="branch15.batch24.pear_Dec2023.noClone"
# BASE_PREFIX="pear_Dec2023.noClone.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 run
# s02 has been run
# BATCH="branch15.batch30.noclone.removeConflictss"
# BASE_PREFIX="pear_Dec2023.noClone.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 has been run
# s02 has been run 
# BATCH="branch15.batch42.full"
# BASE_PREFIX="pear_Dec2023.full.Combine_chr.geno20_maf005.anno.syno.thin8k"

##### Estern samples
# s01 run
# s02 run
# BATCH="branch15.batch22.pear_Dec2023.Asia"
# BASE_PREFIX="pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 run
# s02 run
# BATCH="branch15.batch25.pear_Dec2023.Asia"
# BASE_PREFIX="pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 has been run
# s02 has been run
# BATCH="branch15.batch26.pear_Dec2023.Asia"
# BASE_PREFIX="pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 has been run
# s02 has been run
# s03 has been run
# BATCH="branch15.batch31.Asia.removeConflicts"
# BASE_PREFIX="pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k"

##### Western samples
# s01 has been run
# s02 has been run
BATCH="branch17.batch45.Europe"
BASE_PREFIX="pear_Jul2024_ref_comm.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k"


# s01 has been run
# s02 has been run
# BATCH="branch15.batch43.pear_Dec2023.Europe"
# BASE_PREFIX="pear_Dec2023.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 have been run
# s02 run
# BATCH="branch15.batch23.pear_Dec2023.Europe"
# BASE_PREFIX="pear_Dec2023.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k"

# s01 have been run
# s02 have been run
# BATCH="branch15.batch27.pear_Dec2023.Europe"
# BASE_PREFIX="pear_Dec2023.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k"


# BATCH="branch15.batch18.pear_Dec2023_noclone_no34Minority_no8conflict"
# BASE_PREFIX="pear_Dec2023_noclone_no34Minority_no8conflict.Combine_chr.geno20_maf005.anno.syno.thin8k"

# Cultivar
# s01 have been run
# BATCH="branch15.batch28.Cultivar"
# BASE_PREFIX="pear_Dec2023.Cultivar.Combine_chr.geno20_maf005.anno.syno.thin8k"

# Wild
# s01 have been run
# s02 have been run
# BATCH="branch15.batch29.Wild"
# BASE_PREFIX="pear_Dec2023.Wild.Combine_chr.geno20_maf005.anno.syno.thin8k"


PROJECT_DIR="/shared/ifbstor1/projects/pear_snp3/pear"

########################################################################
# The following variables will be called by other scripts, e.g. s01,s02
########################################################################
WORK_DIR="${PROJECT_DIR}/run_FastStructure"

BED_PREFIX="/${PROJECT_DIR}/run_Population_genetic_filter_refer_on_communis/output/\
${BATCH}/\
${BASE_PREFIX}"

OUTPUT_DIR="${WORK_DIR}/output/${BATCH}"

# the minimum and maximum K to run fastStructure
MIN_K=2
MAX_K=16

LOAD_GHOSTSCRIPT='source activate ghostscript'

LOAD_PERL5="module load perl/5.26.2"

LOAD_R="module load r/4.2.3"

CLUMPAK_DIR="${WORK_DIR}/bin/CLUMPAK.NOGIT./26_03_2015_CLUMPAK/CLUMPAK"

# make sure the dependent perl5 lib is installed in the following path via cpanm or other methods
export PERL5LIB="/shared/ifbstor1/home/ynie/perl5/lib/perl5":$PERL5LIB
