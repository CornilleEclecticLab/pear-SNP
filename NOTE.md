**Please always edit this note.md file on `master` branch!!**
on other branch, do 
`git checkout master -- NOTE.md`

# Branch pre
this branch performed on genotoul and portable (Rscript)

## Conclusion
for the project of Zhang et al 2021, we calculate the kinship set the cutoff at 0.354, then remain
251, and *`211`* individuals with sequencing depth > 10x.

# QC
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

