#!/usr/bin/env bash

#SBATCH -J s02.run_Astral.sh
#SBATCH -o s02.run_Astral.sh.%J.out
#SBATCH -e s02.run_Astral.sh.%J.err
#SBATCH -c 1

# @File     :   s02.run_Astral.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2024/10/04 10:35:53
#Description:
    # 

module load java-jdk

mkdir -p ../output/s02.run_Astral

cat ../output/s01.sampling_and_run_IQtree/*.treefile \
    > ../output/s02.run_Astral/concatenated.treefile

ls ../output/s01.sampling_and_run_IQtree/*.ufboot \
    > ../output/s02.run_Astral/ufboot.list

java -jar Astral.NOGIT./astral.5.7.8.jar \
    -i ../output/s02.run_Astral/concatenated.treefile \
    -b ../output/s02.run_Astral/ufboot.list \
    -o ../output/s02.run_Astral/astral.species_boot.tree

tail -n 1 ../output/s02.run_Astral/astral.species_boot.tree \
    > ../output/s02.run_Astral/astral.species.tree
