#!/usr/bin/env bash

#SBATCH -J s06.dot_plotter.sh
#SBATCH -o s06.dot_plotter.sh.%J.out
#SBATCH -e s06.dot_plotter.sh.%J.err

# @File     :   s06.dot_plotter.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/02/10 16:22:14
#Description:
    # 

source s00.load_MCScanX.sh
source s00.load_java.sh
source s00_config.py

mkdir -p $s06_dir

for prefix in comm_pyri pyri_comm ; do
    java "${downstream_analyses}/dot_plotter.java" \
        -g ${s04_dir}/${prefix}.gff \
        -s ${s04_dir}/${prefix}.collinearity \
        -c ../input/dot_plot.ctl \
        -o ${s06_dir}/${prefix}.dot_plot.png

    java "${downstream_analyses}/dual_synteny_plotter.java" \
        -g ${s04_dir}/${prefix}.gff \
        -s ${s04_dir}/${prefix}.collinearity \
        -c ../input/dot_plot.ctl \
        -o ${s06_dir}/${prefix}.dual_synteny_plot.png
done