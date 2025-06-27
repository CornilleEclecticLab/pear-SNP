#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s09.assess_simulate_stat.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/11/12 17:26:38
# @Description:
#    

import datetime
import sys
import textwrap
import os

import allel
import glob
import argparse

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output')
sub_output_dir = os.path.join(output_dir,script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)


# Ensure output directory exists
os.makedirs(sub_output_dir, exist_ok=True)

# Simulation directory
parser = argparse.ArgumentParser(
    description="Process genotype data and calculate statistics.")
parser.add_argument(
    '-f', '--folder', type=str, 
    default=os.path.join(output_dir, 's08.simulate_by_population_piecewise'),
    help="The folder of VCF files to be calculated.")
args = parser.parse_args()




# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Function to read vcf
def read_vcf(file):
    return allel.read_vcf(file)


# Function to calculate pi, ho, he
def popStat(file):
    vcf_header = allel.read_vcf_headers(file)
    for header in vcf_header:
        for line in header:
            if line.startswith("##contig"):
                contig_length = int(line.split('length=')[-1].replace('>', '')) # Only one contig in the simulation vcf file

    callset = read_vcf(file)
    
    pos = callset['variants/POS']
    gt = allel.GenotypeArray(callset['calldata/GT'])
    ac = gt.count_alleles()
    af = ac.to_frequencies()

    pi = allel.sequence_diversity(pos, ac, start = 1, stop = contig_length)
    
    theta_hat_w = allel.watterson_theta(
        pos, ac, start=1, stop=contig_length)
    ho = allel.heterozygosity_observed(gt)
    mean_ho = ho.mean()
    
    he = allel.heterozygosity_expected(af,ploidy=2)
    mean_he = he.mean()
    
    fa = allel.inbreeding_coefficient(gt)
    mean_fa = fa.mean()

    return [pi, theta_hat_w, mean_ho, mean_he, mean_fa]


# Function to write results to csv
def write_stats_to_csv(results, file):
    with open(file, 'w') as fw:
        fw.write('population\ttheat_hat_w\tpi\tho\the\tfa\n')
        for pop, stats in results.items():
            stats = '\t'.join(map(str, stats))
            fw.write(f"{pop}\t{stats}\n")

def main():
    # results = {}
    # for pop, files in simulation_files.items():
        # print(f'Processing population: {pop}')
        # result = popStat(files[0])
        # results[pop] = result 
    # write_stats_to_csv(results, os.path.join(
        # sub_output_dir, f"{script_basename}.csv"))
    glob_vcf = glob.glob(args.folder + '/*norm.vcf.gz')
    print(wrap79(f'{len(glob_vcf)} vcf files found in {args.folder}'))
    for file_path in glob_vcf:
        dic = {}
        file_name = os.path.basename(file_path)
        print(f'Processing {file_name}')
        dic[file_name] = popStat(file_path)
        write_stats_to_csv(dic, os.path.join(sub_output_dir, f"{file_name}.stat.csv"))    
        

main()

end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))