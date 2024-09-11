#!/usr/bin/env bash

# @File     :   s01.run_PopLDdecay.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/07/10 17:19:00
#Description:
    # 

# Load config
source s00.config.sh

# Start
start_time = $(date +%s)
"${PopLDdecayBin}/PopLDdecay" \
    -InVCF "${VCF}" \
    -OutStat "${OUTPUTDIR}/${PopListPrefix}" \
    -SubPop "${PopList}" \
    -OutType 2 &&\
echo "Done!"
echo "Finished at $(date)"
echo "Finished in $(( $(date +%s) - $start_time )) seconds"
