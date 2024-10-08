**Please always edit this note.md file on `master` branch!!**
on other branch, do 
`git checkout master -- NOTE.md`

# Branch pre
this branch performed on genotoul and portable (Rscript)

## Conclusion
for the project of Zhang et al 2021, we calculate the kinship set the cutoff at 0.354, then remain
251, and *`211`* individuals with sequencing depth > 10x.

# Branch QC
## Two bad files the sequences are not paired. **REMOVE them**
- SAMN12691580 
- SNAM12691720, 

## Three files are missing. **YES, They are not accessable on ENA or NCBI**
- SRR7135561
- SRR7135583
- SRR7135514  

## Some sample re-performed fastp for severl times as they didn'd pass 'Per Base Sequence Content'.
### Three samples from Teng dataset run fastp for two times, s03 or s05, s05 using input s05.1.fastp_cutoff.dic.txt and the same input fastq with s03, as there were warning on 'Per Base Sequence Content'.
Following are the output of s03, they run fastp for one time, the input fastq files of these are performed s05
- P11-8.s03.old  
- P13-6.s03.old  
- P6-1.s03.old 

### Two samples from Teng dataset run fastp for three times,s03 or s05 or s05, s05 using input s05.2.fastp_cutoff.dic.txt
- P11-8
- P6-1


## Three Europe Wild samples we sequenced removed for low quality
During the Further filter steps:
A326, A352, A794

# 2023-03-17 I checked the GATK SNP calling steps *.err and *.out files:
Did follwing steps on later record accounts on different clusters:
```bash
    cat s01.*/*.out
    grep failed s01.*/*.out
    tail -n 2 s01.*/*.err | awk 'NR%4==2{print}'
    tail -n 2 s01.*/*.err | awk 'NR%4==2{print}' | grep -v done
```
For secondary dir:
```bash
    grep failed s01.*/*/*.out
```

Check the job status using seff  
using hash to renew the job id for every script as the last running
```bash
    cat s01.*/*.txt |sort -n -k2 |awk '{a[$1]=$2}END{for (i in a) print i,a[i]}'|cut -d' ' -f2 |xargs -I% seff %| grep State
    cat s01.*/*.txt |sort -n -k2 |awk '{a[$1]=$2}END{for (i in a) print i,a[i]}'|cut -d' ' -f2 |xargs -I% seff %| grep State | grep -v COMPLETED
```

using this command to list all the dirs and secondary dir (if exit):
```bash
    find s01* -type d
```
ensure there is no secondary dir:
```
    find s01*/* -type d
```

ruilin@genoutoul:  
    [-] s01.branch.bam_to_gvcf_and_statics.ECLECTIC  
    [-] s01.fastq_to_gvcf_and_statics.pear.Teng  
    [-] s01.fastq_to_gvcf_and_statics.pear.Wu2018EC   
    [-] s01.fastq_to_gvcf_and_statics.pear.Wu2018EW  
ynie@genoutoul:  
    [-] s01.fastq_to_gvcf_and_statics.pear.Li2021  
    [-] s01.fastq_to_gvcf_and_statics.pear.Li2021_part2  
    [-] s01.fastq_to_gvcf_and_statics.pear.Li2021_part2_2  
    [-] s01.fastq_to_gvcf_and_statics.pear.Loquat  
    [-] s01.fastq_to_gvcf_and_statics.pear.Wu2018AC  
    [-] s01.fastq_to_gvcf_and_statics.pear.Wu2018AW  
    [-] s01.fastq_to_gvcf_and_statics.pear.Wu2018_part2  
    [-] s01.fastq_to_gvcf_and_statics.pear.Zhang2021  
    [-] s01.fastq_to_gvcf_and_statics.pear.Zhang2021_b2_p1  
    [-] s01.fastq_to_gvcf_and_statics.pear.Zhang2021_b3  
ynie@IFB:   
    [-] s01.fastq_to_gvcf_and_statics.pear.Zhang2021_b2_p2  
    [x] s01.fastq_to_gvcf_and_statics.pear.Zhang2021_b2_p2/DISCARD.run_with_problem_tmp  
    [x] s01.fastq_to_gvcf_and_statics.pear.Zhang2021_b2_p2/DISCARD.run1.gatk_variant_calling_stoped_as_Java_env_memory_is_not_enough  
    [-] s01.fastq_to_gvcf_and_statics.pear.Zhang2021_b2_p3c  

As you see, the *.err files in `DISCARD.run_with_problem_tmp` and `DISCARD.run1.gatk_variant_calling_stoped_as_Java_env_memory_is_not_enough` could not pass the check steps, showing the steps are effective. 

# On 2024-10-07
Due to the IFB shoutdown maintain until 2024-10-11 18:00. some data transfermed to genobioinfo for analysing, You should transfer back if changed.
[ ] pear/run_Stacks/input
[ ] pear/run_Stacks/output
[ ] pear/run_RAiSD_refer_on_communis/input
[ ] pear/run_RAiSD_refer_on_communis/output
[ ] pear/run_Population_genetic_filter/output/branch18.batch47.full_noAdmix_base_branch15_batch42_unlink
[ ] pear/run_Population_genetic_filter/output/branch15.batch46.noAdmix
[ ] pear/run_OmegaPlus_refer_on_communis/input
[ ] pear/run_OmegaPlus_refer_on_communis/output
[ ] pear/run_SweeD_refer_on_communis/output
[ ] pear/run_SweeD_refer_on_communis/output 

Regenerated results on:
[ ] pear/run_IQtree_sampling/output
