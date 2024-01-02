#!/usr/bin/env bash
# 2023-11-30 
for i in $(cat s01.input.whole_pear_filter_passed_sites_vcf.file_list.txt) ; 
    do echo $i > s01.input.$(echo $i | cut -d "." -f4).path.txt; 
done


for i in s01.input.GWHBAOS00000*.path.txt; do
    CHR_ID=$(echo "${i}" | cut -d "." -f3)
    JOB_ID=$(sbatch --mem=32G s01.branch15.prepare_test_vcf.sh "$(cat "${i}")" "pear_Oct2023_whole.${CHR_ID}" "$(cat "${i}")" | awk '{print $4}')
    echo -e "${CHR_ID}\t${JOB_ID}" >> s01.job_id.txt
done
