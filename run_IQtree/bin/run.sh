#!/bin/bash

#SBATCH -J job
#SBATCH -o job.out
#SBATCH -e job.err
module purge

module load bioinfo/iqtree-2.0.6

# This -m MFP is wrong
# iqtree -s ./pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.min4.phy -st DNA -m MFP -bb 1000 -bnni -alrt 1000 -nt 6 -o "CRR180425,CRR180437,CRR180442,CRR180443,CRR180449"

# Convert file format from vcf to phy
# python3 vcf2phylip.py -r -i pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.vcf.gz


# Without -alrt 1000, it would be faster
iqtree -s ./pear.Combine_Chr.geno20_maf005.anno.syno.thin8k.min4.phy -st DNA -m GTR+ASC -bb 1000 -bnni -alrt 1000 -nt 10


