#!/usr/bin/env bash

# @File     :   bin.stat.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2022/12/07 12:36:00
#Description:
	# 

seqkit version >> info.software_version.txt

seqkit stat -a ../input/hap/Ppy362_CCS3.HiFiasm.*.fasta.gz > ../output/stat.txt
seqkit stat -a ../input/primary_alternatif/*.fasta.gz >> ../output/stat.txt
