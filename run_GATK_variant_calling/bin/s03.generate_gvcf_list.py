#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     :   s03.generate_gvcf_list.py
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023-01-05 18:00:27
#Description:


import datetime,os
start_time = datetime.datetime.now()
print("{0:=^80}".format(' Start '))



# >>>>><<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>             Here is the start fo setting part.                     >>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Set path
work_dir = '/data/atipe-workspace/ynie/pear/run_GATK_variant_calling'
groups      = 'pear.Li2021 pear.Teng pear.Wu2018AC pear.Wu2018EC pear.Zhang2021'.split()


## NO need to change >>
input_dir   = work_dir+'/input'
output_dir  = work_dir+'/output'
prefix = 's03.generate_gvcf_list'


output3_dir = output_dir+'/'+prefix

## NO need to change Change <<



# The path of a txt file records the scaffolds IDs in the reference genome. If there are too many scaffolds, you could save some of them in another file, and write it's path in the following file. In this file,  {input_dir} could be recognized as the variant.
chromosomes = input_dir+'/all_scaffolds.list'



# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<                   Here is the end fo setting part.                 <<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        


# Create dir for saving scripts and output
os.system('mkdir -p '+output3_dir)



# Create the g.vcf.gz file list for every scaffold
## Read chromosomes list
with open(chromosomes,'r') as fo:
    chromosome = fo.read().strip().format(input_dir=input_dir).split()

## For every chromosome
for chr in chromosome:
    chr_base = os.path.splitext(os.path.basename(chr))[0]

    ## Create the chromosomes list file
    gvcf_list = output3_dir+'/groups.'+chr_base+".list"
    with open(gvcf_list, 'w') as fo:
        
        ## For every group
        for group in groups:
            
            # The path of a text file record the paths of fastq file, in which, the line starts with # would be ignored, every line for one record(fastq file), and the fastq files of the same sample should be in a same dir with their sample name.
            fastq_list = input_dir+"/"+group+".input.clean_data_list.txt"

            ## Read samples list
            with open(fastq_list, 'r') as fr:

                dic = {}
                for line in fr:
                    line = line.strip()

                    if line.startswith('#'):
                        continue

                    path = line
                    file_name = os.path.basename(path)
                    dir_name = os.path.basename(os.path.dirname(path))
                    sample_name = dir_name

                    dic[sample_name] = dic.get(sample_name, [])
                    dic[sample_name].append(path)

            ## For every sample
            for sample, path in dic.items():
                output1_dir = output_dir+'/s01.fastq_to_gvcf_and_statics.'+group
                fo.write(output1_dir+'/'+group+'.' +
                            sample+'.'+chr_base+'.g.vcf.gz\n')




end_time = datetime.datetime.now()
print('')
print(' END '.center(80,'='))
print(str(end_time-start_time).center(80))
