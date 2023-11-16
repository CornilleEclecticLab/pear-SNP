#!/usr/bin/env bash

# @File     :   s00.general_command.seff_check_slurm_jobs_status.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/10/24 12:05:29
#Description:
    #


# Usage
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <slurm_job_id.txt>"
    exit 1
fi

# Main
while read -r line; do
    ID=$(awk '{print $NF}' <<< "$line")
    NF=$(awk '{print NF}' <<< "$line")

    if [[ $NF -eq 2 ]]; then
        seff $ID 2>&1 | grep State | awk -v line="$line" '{print line "\t" $0}'
    elif [[ $NF -eq 1 ]]; then
        seff $ID 2>&1 | grep State | awk '{print "Job_ID: '$line'\t" $0}'
    fi
done < "$1"