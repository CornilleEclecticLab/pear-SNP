#!/usr/bin/env bash

# @File     :   s06.1.prepear_data_for_stats_clean_reads_length.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/04/12 15:24:56
#Description:
    # 


# mkdir -p ../input/s06.input
# ln -s ../output/s03.fastp_and_qc/* ../input/s06.input/

# run Multi QC
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

multiqc -n /data/atipe-workspace/ynie/pear/run_QC/input/s06.input/s06_input.html \
/data/atipe-workspace/ynie/pear/run_QC/input/s06.input/*/*/ 
