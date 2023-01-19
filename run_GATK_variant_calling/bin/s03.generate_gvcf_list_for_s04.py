#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     :   s03.generate_gvcf_list_for_s04.py
# @Version  :   2.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023-01-05 18:00:27
#Description:
    # v2.0.0 
    # 2023-01-15 21:32:23
    # Compatible with extra g.vcf file path
    


import datetime,os,warnings
start_time = datetime.datetime.now()
print("{0:=^80}".format(' Start '))



# >>>>><<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>             Here is the start fo setting part.                     >>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Set path
work_dir = '/data/atipe-workspace/ynie/pear/run_GATK_variant_calling'
groups   = 'pear.Li2021 pear.Teng pear.Wu2018AC pear.Wu2018EC pear.Zhang2021'.split()
groups   = 'pear.Li2021 pear.Teng pear.Wu2018AC pear.Wu2018EC'.split()  

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



# Read the allowed sample ID list.
allow_samples=[]
with open(input_dir+"/s03.input.allow_sample_ID_list.txt",'r') as fr:
    for line in fr:
        line = line.strip()
        if line.startswith("#"):
            continue
        elif line != '':
            allow_samples.append(line)
            
            
    if len(allow_samples)  != len(set(allow_samples)):
        warnings.warn("[ERROR]: Duplicated samples in :" +
                     input_dir+"/s03.input.allow_sample_ID_list.txt")
        exit()
        


# Create the g.vcf.gz file list for every scaffold
## Read chromosomes list
with open(chromosomes,'r') as fo:
    chromosome = fo.read().strip().format(input_dir=input_dir).split()
    # chromosome = map(lambda chr: os.path.splitext(os.path.basename(chr))[0], chromosome)


## For every chromosome
for chr in chromosome:
    chr_base = os.path.splitext(os.path.basename(chr))[0]
    pass_sample = []

    ## Create the chromosomes list file
    gvcf_list = output3_dir+'/pear.'+chr_base+".gvcf.list"
    with open(gvcf_list, 'w') as fo:
        
        
        ## For every line in g.vcf file list
        s03input_gvcf_list=input_dir+'/s03.input.gvcf_list.txt'
        with open(s03input_gvcf_list,'r') as fr:
            for line in fr:
                line = line.strip()
                if line.startswith("#") or line == '':
                    continue
                try:
                    gvcf_sample = os.path.basename(line).split('.')[-5]
                    gvcf_chr    = os.path.basename(line).split('.')[-4]
                except:
                    print(line)
                
                if gvcf_sample in allow_samples and gvcf_chr == chr_base:
                    fo.write(line+'\n')
                    pass_sample.append(gvcf_sample)
    
    
    if len(pass_sample) != len(allow_samples):
        warnings.warn("[ERROR] Unequal sample amount in pass list from allow list!")
        print("Interrupted at chromosome: "+chr_base)
        exit()

             

end_time = datetime.datetime.now()
print('')
print(' END '.center(80,'='))
print(str(end_time-start_time).center(80))
