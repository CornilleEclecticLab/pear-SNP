
#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     :   s06.variant_invariant_concat.py
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/02/03 15:18:53
#Description:
    # 


import datetime
import os
import getpass
import sys

start_time = datetime.datetime.now()
print("{0:=^79}".format(' Start '))



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>                Here is the start fo setting part.                  >>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Set path
user = getpass.getuser()
work_dir = '/work' + user + '/pear/GATK_variant_calling'
group = 'pear.Loquat'


## NO need to change >>
prefix = 's06.variant_invariant_concat' + group

bin_dir = work_dir + '/bin'
input_dir = work_dir + '/input'
output_dir = work_dir + '/output'

scripts_dir = bin_dir + '/' + prefix
input5_dir = output_dir + 's05.hard_filter_chr_vcf.' + group
output6_dir = output_dir + '/' + prefix
## NO need to change Change <<



# Set input file names

## The path of a text file record the paths of fastq file, in which, the line starts with # would be ignored, every line for one record(fastq file), and the fastq files of the same sample should be in a same dir with their sample name.
fastq_list = input_dir + "/input.clean_data_list.txt"

## Reference genome (fasta) path
ref_genome = input_dir + '/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.genome.fasta'

# The path of a txt file records the scaffolds IDs in the reference genome. If there are too many scaffolds, you could save some of them in another file, and write it's path in the following file. In this file, {input_dir} could be recognized as the variant.
chromosomes = input_dir + '/all_scaffolds.list'



# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<                Here is the end fo setting part.                    <<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



# Ask if run the scripts
run = ''
print("Welcome!")
while run == '':
    run = input("Please chose 1 or 2 to continue: \n\t 1) generate the scripts and run them (submit the jobs using sbatch)\
                \n\t 2) only generate the scripts, without running them.\n")
    
    if run != '1' and run !='2':
        print("Unrecognized input, please try again.\n")
        run = ''



# Create dir for saving scripts and output
os.system('mkdir -p ' + scripts_dir)
os.system('mkdir -p ' + output6_dir)


# Change the working dir
os.chdir(scripts_dir)



# Create the g.vcf.gz file list for every scaffold
## Read chromosomes list
with open(chromosomes,'r') as fo:
    chromosome = fo.read().strip().format(input_dir=input_dir).split()


for chr in chromosome:
    chr_base = os.path.splitext(os.path.basename(chr))[0]
    
    ## Create the bash script
    sh = scripts_dir + '/' + prefix + '.' + chr_base +".sh"
    with open(sh,'w') as fo: 
        # Shebang
        content_header = "#!/usr/bin/env bash \n\
            # bash file generated" + __file__ +"\n"
        
        fo.write(content_header)
        
        ## SBATCH settings
        content = '#SBATCH -J s01.' + chr_base + '\n#SBATCH -o ' + \
            prefix + '.' + chr_base + '.out \n' + '#SBATCH -e ' + prefix + '.' + chr_base + '.err \n'
        fo.write(content)

        ## Load softwares
        content = '''
module purge

module load bioinfo/gatk-4.1.9.0
module load bioinfo/bcftools-1.14

'''
        fo.write(content)
        
        ## bcftools filter, .vcf.gz
        content = '''

bcftools filter -S . -e 'FMT/DP<5 | FMT/RGQ<20 | FMT/DP>100' \\
    {input5_dir}/{group}.{chr_base}.rmQFI.combine.filter_passed_sites.vcf.gz \\
    | bcftools filter -i 'ALT="."' \\
    | bcftools filter -e 'F_MISSING > 0.2' -O z4 \\
    -o {output6_dir}/{group}.{chr_base}.filtered_pixy_invariant.vcf.gz
    
tabix -p vcf {output6_dir}/{group}.{chr_base}.filtered_pixy_invariant.vcf.gz


bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20 | FMT/DP>100' \\
    {input5_dir}/{group}.{chr_base}.rmQFI.combine.filter_passed_sites.vcf.gz \\
    | bcftools filter --SnpGap 10 \\
    | bcftools view -m2 -M2 -v snps \\
    | bcftools filter -e 'F_MISSING > 0.2' -O z4 \\
    -o {output6_dir}/{group}.{chr_base}.filtered_pixy_variant.vcf.gz

tabix -p vcf {output6_dir}/{group}.{chr_base}.filtered_pixy_variant.vcf.gz


bcftools concat \\
--allow-overlaps \\
{output6_dir}/{group}.{chr_base}.filtered_pixy_invariant.vcf.gz \\
{output6_dir}/{group}.{chr_base}.filtered_pixy_variant.vcf.gz \\
-O z4 \\
-o {output6_dir}/{group}.{chr_base}.filtered_pixy_concat.vcf.gz

tabix -p vcf {output6_dir}/{group}.{chr_base}.filtered_pixy_concat.vcf.gz

'''.format(ref_genome=ref_genome,input5_dir=input5_dir,chr_base=chr_base,output6_dir=output6_dir,group=group)
        fo.write(content)
        
    ## Submit this script as a job
    if run == '1':
        command = "sbatch -c 2 --mem=5G "+sh
        submit = os.popen(command, 'r')
        job_id = submit.read().strip().split()[-1]
        print("\n\nInformation: Dealing with "+ group + ' ' + chr_base)
        print("\nInformation: A job has been submitted: \n\t" + command +'\n')



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
