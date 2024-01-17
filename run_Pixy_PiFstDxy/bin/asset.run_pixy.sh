#!/usr/bin/env bash

# @File     :   building.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2023/05/04 22:53:00
#Description:
	# 



# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/local/bioinfo/src/Miniconda/Miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/local/bioinfo/src/Miniconda/Miniconda3/etc/profile.d/conda.sh" ]; then
        . "/usr/local/bioinfo/src/Miniconda/Miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/usr/local/bioinfo/src/Miniconda/Miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate pixy

pixy --stats pi fst dxy \
--vcf test.vcf.gz \
--populations Ag1000_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 4 \
--output_folder output \
--output_prefix pixy_output \
--bypass_invariant_check yes


