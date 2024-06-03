#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s05.split_time_estimation.py
# @Version  : 1.2.1
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023-07-06 19:32:47
# @Description:
#    This script is used to generate sub-scripts for the "smc++ split" command to estimate the split time.

#    Version 1.0.1 2023-07-07 00:35:19
#    Fix a bug in the "smc++ split" command, that is, wrong input file name (pop2.pop1 rather than pop1.pop2).

#    Version 1.0.2: 2023-07-09 10:37:26
#    Shortened the slurm job file name.
#    Fixed bugs.

#    Version 1.0.3: 2023-07-11 16:19:23
#    Didn't use for chr in chr_list to generate the same sub-script for each population pair.

#    Version 1.0.4: 2023-07-13 17:28:26
#    Book more cpu cores for the "smc++ split" command.

#    Version 1.1.0: 2023-08-06 16:53:57
#    Try to infer split time pop2 from pop1, and use less smc files.

#    Version 1.2.0: 2023-08-07 22:45:11
#    Use both pop21 and pop12 smc files (medium)

#    Version 1.2.1: 2024-05-21 18:13:50
#    Don't need load singularity module.


import datetime
import os
import sys
import textwrap
from warnings import warn
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.2.1"
s01_script_basename = "s01.vcf2smc_by_chr"
s02_script_basename = "s02.estimate_by_population"
s04_script_basename = "s04.vcf2smc_to_prepare_for_split"
script_basename  = 's05.split_time_estimation'
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = work_dir+'/input/'
output_dir = work_dir+'/output/'+script_basename
sub_script_dir = work_dir+'/bin/'+script_basename
s01_output_dir = work_dir+'/output/'+s01_script_basename
s02_output_dir = work_dir+'/output/'+s02_script_basename
s04_output_dir = work_dir+'/output/'+s04_script_basename



vcf = input_dir+"/s01.input.vcf.gz"
chromosome_list = input_dir+"/s01.scaffolds_list.txt"
individual_population_list = input_dir+"/s01.individuals_and_populations_list.txt" # Format: individual_name population_name
population_pair_list = input_dir+"/s04.population_pair_list.txt" # Format: population1 population2
load_singularity = '# module load system/singularity-3.7.3 # singularity is installed and callable on the cluster without loading a module'
cpu_cores = "20"



os.system("mkdir -p "+sub_script_dir)
os.system("mkdir -p "+output_dir)



def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)



## Read the chromosome list
with open(chromosome_list,'r') as fr:
    chr_list = fr.read().strip().split('\n')
    info=f'Reading the "chromosomes list": {",".join(chr_list)}'
    print(wrap(info))



## Read the individuals and populations list
population_dict = {}
with open(individual_population_list,'r') as fr:
    for line in fr:
        if line.startswith("#") or line.strip()=='':
            continue 
        line = line.strip().split()
        individual = line[0]
        population = line[1]
        info=f'Reading the "individuals_and_populations_list": {individual} {population}'
        print(wrap(info))

        population_dict[population] = population_dict.get(population,list())
        population_dict[population].append(individual)



## Read the population pair list
if os.path.exists(population_pair_list):
    ls = list()
    with open(population_pair_list,'r') as fr:
        for line in fr:
            if line.startswith("#") or line.strip()=='':
                continue

            line = line.strip().split()
            pop1 = line[0]
            pop2 = line[1]

            if pop1 == pop2:
                warn(wrap(f'Population "{pop1}" is paired with itself.'))
                continue
            if (pop1,pop2) in ls or (pop2,pop1) in ls:
                warn(wrap(f'Population pair "{pop1} {pop2}" is duplicated.'))
                continue

            ls.append((pop1,pop2))
            info=f'Reading the "population_pair_list": {pop1} {pop2}'
            print(wrap(info))
    population_pair_list = ls


elif not os.path.exists(population_pair_list):
    warn(f'Did not find the "{population_pair_list}" file, it will define the populations pairs for SMC++ with this format: \n    Population1  Population2.\n You could run {s04_script_basename}.py to generate this file.')
    exit(1)



print(wrap(f'Population pair list: {",".join([str(i) for i in population_pair_list])}'))

if len(population_pair_list) == 0:
    warn(wrap(f'No population pair was found in the "{population_pair_list}" file. If you do not want to use the "{population_pair_list}" file, you could remove it.'))
    exit(1)



def individual_list_to_string(individual_list):
    return ",".join(individual_list)


def get_smc_list_for_pop_by_chr(output_dir, pops):
    return " \\\n        ".join([f"{output_dir}/{pop}.{chr}.smc.gz" for pop in pops for chr in chr_list])


for (pop1,pop2) in population_pair_list:
    sub_script_basename = f"split.{pop1}.{pop2}.pair"
    sub_script = f"{sub_script_dir}/{sub_script_basename}.sh"

    inds1 = individual_list_to_string(population_dict[pop1])
    inds2 = individual_list_to_string(population_dict[pop2])


    with open(sub_script,"w") as fo:
        content = f'''#!/usr/bin/env bash

# This sub-script was generated by {script_basename}.py V{version} at {start_time}


#SBATCH -J {sub_script_basename}
#SBATCH -o {sub_script_basename}.%J.out
#SBATCH -e {sub_script_basename}.%J.err
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem 120G


{load_singularity}


singularity run -B  {work_dir}:{work_dir} \\
    {work_dir}/bin/smcpp.sif split \\
        -o {output_dir}/{pop1}.{pop2}/ \\
        --cores {cpu_cores} \\
        {s02_output_dir}/{pop1}/model.final.json \\
        {s02_output_dir}/{pop2}/model.final.json \\
        {get_smc_list_for_pop_by_chr(s01_output_dir, [pop1,pop2])} \\
        {get_smc_list_for_pop_by_chr(s04_output_dir, [f'{pop1}.{pop2}',f'{pop2}.{pop1}'])} 
'''
        fo.write(content)

# Run pop2 pop1 split
    sub_script_basename = f"split.{pop2}.{pop1}.pair"
    sub_script = f"{sub_script_dir}/{sub_script_basename}.sh"

    inds1 = individual_list_to_string(population_dict[pop1])
    inds2 = individual_list_to_string(population_dict[pop2])


    with open(sub_script,"w") as fo:
        content = f'''#!/usr/bin/env bash

# This sub-script was generated by {script_basename}.py V{version} at {start_time}


#SBATCH -J {sub_script_basename}
#SBATCH -o {sub_script_basename}.%J.out
#SBATCH -e {sub_script_basename}.%J.err
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem 120G


{load_singularity}


singularity run -B  {work_dir}:{work_dir} \\
    {work_dir}/bin/smcpp.sif split \\
        -o {output_dir}/{pop2}.{pop1}/ \\
        --cores {cpu_cores} \\
        {s02_output_dir}/{pop2}/model.final.json \\
        {s02_output_dir}/{pop1}/model.final.json \\
        {get_smc_list_for_pop_by_chr(s01_output_dir, [pop1,pop2])} \\
        {get_smc_list_for_pop_by_chr(s04_output_dir, [f'{pop1}.{pop2}',f'{pop2}.{pop1}'])} 
'''
        fo.write(content)


# Run pop2 po1 with less smc
    sub_script_basename = f"split.{pop2}.{pop1}.less.pair"
    sub_script = f"{sub_script_dir}/{sub_script_basename}.sh"

    inds1 = individual_list_to_string(population_dict[pop1])
    inds2 = individual_list_to_string(population_dict[pop2])


    with open(sub_script,"w") as fo:
        content = f'''#!/usr/bin/env bash

# This sub-script was generated by {script_basename}.py V{version} at {start_time}


#SBATCH -J {sub_script_basename}
#SBATCH -o {sub_script_basename}.%J.out
#SBATCH -e {sub_script_basename}.%J.err
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem 120G


{load_singularity}


singularity run -B  {work_dir}:{work_dir} \\
    {work_dir}/bin/smcpp.sif split \\
        -o {output_dir}/{pop2}.{pop1}.less/ \\
        --cores {cpu_cores} \\
        {s02_output_dir}/{pop2}/model.final.json \\
        {s02_output_dir}/{pop1}/model.final.json \\
        {get_smc_list_for_pop_by_chr(s04_output_dir, [f'{pop2}.{pop1}'])} 
'''
        fo.write(content)


# Run pop1 po2 with less smc
    sub_script_basename = f"split.{pop1}.{pop2}.less.pair"
    sub_script = f"{sub_script_dir}/{sub_script_basename}.sh"

    inds1 = individual_list_to_string(population_dict[pop1])
    inds2 = individual_list_to_string(population_dict[pop2])


    with open(sub_script,"w") as fo:
        content = f'''#!/usr/bin/env bash

# This sub-script was generated by {script_basename}.py V{version} at {start_time}


#SBATCH -J {sub_script_basename}
#SBATCH -o {sub_script_basename}.%J.out
#SBATCH -e {sub_script_basename}.%J.err
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem 120G


{load_singularity}


singularity run -B  {work_dir}:{work_dir} \\
    {work_dir}/bin/smcpp.sif split \\
        -o {output_dir}/{pop1}.{pop2}.less/ \\
        --cores {cpu_cores} \\
        {s02_output_dir}/{pop1}/model.final.json \\
        {s02_output_dir}/{pop2}/model.final.json \\
        {get_smc_list_for_pop_by_chr(s04_output_dir, [f'{pop1}.{pop2}'])} 
'''
        fo.write(content)




# Run pop2 po1 with medium smc
    sub_script_basename = f"split.{pop2}.{pop1}.medium.pair"
    sub_script = f"{sub_script_dir}/{sub_script_basename}.sh"

    inds1 = individual_list_to_string(population_dict[pop1])
    inds2 = individual_list_to_string(population_dict[pop2])


    with open(sub_script,"w") as fo:
        content = f'''#!/usr/bin/env bash

# This sub-script was generated by {script_basename}.py V{version} at {start_time}


#SBATCH -J {sub_script_basename}
#SBATCH -o {sub_script_basename}.%J.out
#SBATCH -e {sub_script_basename}.%J.err
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem 120G


{load_singularity}


singularity run -B  {work_dir}:{work_dir} \\
    {work_dir}/bin/smcpp.sif split \\
        -o {output_dir}/{pop2}.{pop1}.medium/ \\
        --cores {cpu_cores} \\
        {s02_output_dir}/{pop2}/model.final.json \\
        {s02_output_dir}/{pop1}/model.final.json \\
        {get_smc_list_for_pop_by_chr(s04_output_dir, [f'{pop1}.{pop2}',f'{pop2}.{pop1}'])} 
'''
        fo.write(content)


# Run pop1 po2 with medium smc
    sub_script_basename = f"split.{pop1}.{pop2}.medium.pair"
    sub_script = f"{sub_script_dir}/{sub_script_basename}.sh"

    inds1 = individual_list_to_string(population_dict[pop1])
    inds2 = individual_list_to_string(population_dict[pop2])


    with open(sub_script,"w") as fo:
        content = f'''#!/usr/bin/env bash

# This sub-script was generated by {script_basename}.py V{version} at {start_time}


#SBATCH -J {sub_script_basename}
#SBATCH -o {sub_script_basename}.%J.out
#SBATCH -e {sub_script_basename}.%J.err
#SBATCH --cpus-per-task={cpu_cores}
#SBATCH --mem 120G


{load_singularity}


singularity run -B  {work_dir}:{work_dir} \\
    {work_dir}/bin/smcpp.sif split \\
        -o {output_dir}/{pop1}.{pop2}.medium/ \\
        --cores {cpu_cores} \\
        {s02_output_dir}/{pop1}/model.final.json \\
        {s02_output_dir}/{pop2}/model.final.json \\
        {get_smc_list_for_pop_by_chr(s04_output_dir, [f'{pop1}.{pop2}',f'{pop2}.{pop1}'])} 
'''
        fo.write(content)

# Another way to write the content of the sub-script
'''
        # {s02_output_dir}/{pop1}.*.smc.gz \\
        # {s02_output_dir}/{pop2}.*.smc.gz \\
        # {s04_output_dir}/{pop1}.{pop2}.*.smc.gz \\
        # {s04_output_dir}/{pop2}.{pop1}.*.smc.gz \\
'''



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
