#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     :   s01.raw_qc.py
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2022/12/09 12:19:21
#Description:
    # 

import datetime
start_time = datetime.datetime.now()
print("{0:=^40}".format(' Start '))



import os, warnings


work_dir = "/data/atipe-workspace/ynie/pear/run_QC"

data_dic = {}

with open(work_dir+"/input/raw_data.list",'r') as data_list:
    for line in data_list:
        line = line.strip()
        
        if line.startswith("#"):
            continue
        
        path = line
        fastq_file_name = os.path.basename(path)
        suffix = str(fastq_file_name).split("_")[-1]
        
        # if not suffix.startswith('1') and not suffix.startswith('2') or '_' not in fastq_file_name:
        #     warnings.warn("WARNING: Unexpected file name: "+fastq_file_name)
            
        dirname = os.path.dirname(path)
        samplename = os.path.basename(dirname)
        data_dic[samplename] = data_dic.get(samplename,[])
        data_dic[samplename].append([dirname,fastq_file_name])
        
        
for sample, ls in data_dic.items():
    if len(ls)%2==1:
        warnings.warn(\
"WARNING: this sample has odd number of files, \
it may not sequenced by double-end sequencing: " + sample)


content_head = '''#!/usr/bin/env bash
#Description:
    # this file is created by ../s01.raw_qc.py

'''

shdir=work_dir+'/bin/s01.raw_qc'
if not os.path.exists(shdir):
    os.system("mkdir "+shdir)


for sample,ls in data_dic.items():
    with open("s01.raw_qc/"+sample+".s01.raw_qc.sh",'w') as fo:
        outdir=work_dir+'/output/s01.raw_qc/'+sample
        content = content_head + 'mkdir -p '+outdir+' \n'
        fo.write(content)
        for i in ls:
            dirname = i[0]
            filename = i[1]
            content = 'fastqc -t 2 -f fastq -o ' + outdir +' ' +dirname+'/'+filename+' &&\\\n'
            fo.write(content)
        fo.write("echo done")



end_time = datetime.datetime.now()
print('')
print(' END '.center(40,'='))
print(str(end_time-start_time).center(40))
