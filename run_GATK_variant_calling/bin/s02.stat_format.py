#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     :   s02.stat_format.py
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2022/12/30 11:40:38
#Description:
    # 

import datetime, os, getpass
start_time = datetime.datetime.now()
print("{0:=^40}".format(' Start '))



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>Here is the start fo setting part.>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Set path
work_dir 	= '/work/ynie/pear/run_GATK_variant_calling'
group       = 'pear.Wu2018_part2'


## NO need to change >>
input_dir   = work_dir+'/input'
output_dir  = work_dir+'/output'
prefix      = 's02.stat_format.'+group
scripts_dir = work_dir+'/bin/'+prefix
output1_dir = output_dir+'/s01.fastq_to_gvcf_and_statics.'+group
output2_dir = output_dir+'/'+prefix

## NO need to change Change <<



# Set input file names

## The path of a text file record the paths of fastq file, in which, the line starts with # would be ignored, every line for one record(fastq file), and the fastq files of the same sample should be in a same dir with their sample name.
fastq_list = input_dir+"/"+group+".input.clean_data_list.txt"



# The path of a txt file records the scaffolds IDs in the reference genome. If there are too many scaffolds, you could save some of them in another file, and write it's path in the following file. In this file,  {input_dir} could be recognized as the variant.
chromosomes = input_dir+'/all_scaffolds.list'



# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<Here is the end fo setting part.<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



# Read paths of fastq files
with open(fastq_list,'r') as fo:
    
    dic = {}
    for line in fo:
        line = line.strip()
        
        if line.startswith('#'):
            continue
        
        path = line
        file_name = os.path.basename(path)
        dir_name = os.path.basename(os.path.dirname(path))
        sample_name = dir_name
        
        dic[sample_name] = dic.get(sample_name,[])
        dic[sample_name].append(path)
        


# Create dir for saving scripts and output
os.system('mkdir -p '+output2_dir)



# Create statistical table
with open(output2_dir+'/'+prefix+".txt",'w') as fo:
    header = "Group\tSample\tAvg_depth\tBreadth_coverage\tMarkdu_avg_depth\tMarkdu_breadth_coverage\n"
    fo.write(header)
    
    for sample, path in dic.items():
        
        with open(output1_dir+"/"+group+"."+sample+".sorted.avg_depth.txt",'r') as txt:
            depth = txt.read().strip()
        
        with open(output1_dir+"/"+group+"."+sample+".sorted.breadth_coverage.txt", 'r') as txt:
            coverage = txt.read().strip()
            
        with open(output1_dir+"/"+group+"."+sample+".sorted.markdu.avg_depth.txt", 'r') as txt:
            markdu_depth = txt.read().strip()

        with open(output1_dir+"/"+group+"."+sample+".sorted.markdu.breadth_coverage.txt", 'r') as txt:
            markdu_coverage = txt.read().strip()
        
        fo.write(group+'\t'+sample+'\t'+depth+'\t'+coverage+'\t'+markdu_depth+'\t'+markdu_coverage+'\n')




end_time = datetime.datetime.now()
print('')
print(' END '.center(40,'='))
print(str(end_time-start_time).center(40))
