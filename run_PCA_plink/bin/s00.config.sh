WORKDIR="/shared/ifbstor1/projects/pear_snp3/pear/run_PCA_plink/"

LOAD_PLINK="module load plink/1.90b6.18"

branch="branch15.batch37.Cultivar_pyri_CN"
set="pear_Dec2023.Cultivar_pyri_CN.Combine_chr.geno20_maf005.anno.syno.thin8k"

OUTDIR="$WORKDIR/output/${branch}"
INPUT="/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter/output/${branch}/${set}"
