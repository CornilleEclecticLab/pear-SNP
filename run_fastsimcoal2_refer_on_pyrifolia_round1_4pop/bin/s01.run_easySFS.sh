#!/bin/bash

#SBATCH -J s01.easySFS.sh
#SBATCH -o s01.easySFS.sh.%J.out
#SBATCH -e s01.easySFS.sh.%J.err
#SBATCH -c 8
#SBATCH --mem=64G

# @File     :   s01.run_easySFS.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi modified from Xilong's script
# @Time(Europe/Paris):   2025/09/04
#Description:
	# 

. s00.config.As.sh
. s00.load_conda.sh
. activate easySFS
. s00.load_r.sh


OUTPUT_PREP_VCF=$OUTPUT
OUTPUT=../output


mkdir -p $OUTPUT/easySFS
# prepare a pure vcf file
gunzip -c $OUTPUT_PREP_VCF/east.M5F1.syn_inter_LDprune.vcf.gz \
    >$OUTPUT/easySFS/east.M5F1.syn_inter_LDprune.vcf

# prepare a pop/group file
sed -e 's/pyri_JP/pyriJP/g' \
    -e 's/Sand_CN-SE/SandCNSE/g' \
    -e 's/Sand_CN-SW/SandCNSW/g' \
    -e 's/Sand_CN/SandCN/g' \
    ../input/s01.sample_tab_population.eastern4pop_betu_pash_Sand_CN_Sand_CN-SW.txt \
    >$OUTPUT/easySFS/east.group.set.txt

# easySFS part, here the vcf is not the RADvcf, so use -a to consider the full sites
# Preview
# jubail py3 with pandas, scipy, Aphid env
# IFB easySFS conda env

easySFS.py \
    -i $OUTPUT/easySFS/east.M5F1.syn_inter_LDprune.vcf \
    -p $OUTPUT/easySFS/east.group.set.txt \
    --order betu,pash,SandCN,SandCNSW \
    -a \
    --preview \
    >$OUTPUT/easySFS/preview.output.txt
# find the candidate project from output of preview
## convert output for R plot
perl ./scripts/easySFSpreview_output_format.pl \
    $OUTPUT/easySFS/preview.output.txt \
    >$OUTPUT/easySFS/preview.output.Rdata.txt
## plot them and output the max value
Rscript ./scripts/easySFSpreview_plot.R \
    >$OUTPUT/easySFS/preview.output.proj.result.txt
## the result is
# betu, pash, SandCN, Sand_CNSW
# 16,34,54,20


# Here the projection should also be used for tpl and est file
# easySFS running

## the multiSFS take too much time and memory for big pop size, here I
## remove this function

python3 ./scripts/easySFS_withoutMultiSFS.py \
    -i $OUTPUT/easySFS/east.M5F1.syn_inter_LDprune.vcf \
    -p $OUTPUT/easySFS/east.group.set.txt \
    --order betu,pash,SandCN,SandCNSW \
    -a -f \
    --proj 16,34,54,20 \
    -o $OUTPUT/easySFS \
    --prefix east