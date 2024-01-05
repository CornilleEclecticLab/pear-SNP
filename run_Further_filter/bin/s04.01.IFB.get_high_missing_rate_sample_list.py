#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.01.IFB.get_high_missing_rate_sample_list.py
# @Version  : 1.0.1
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/12/11 14:44:12
# @Description:
#

import argparse
import datetime
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')

# Set the missing rate threshold value (%)
# The value 40% is ensured as the best by testing from 20%, 30%, 40% and 50%.
threshold = 40

version = "1.0.1"


parser = argparse.ArgumentParser(
    prog='Get High Missing Rate Sample',
    description=(
        'This script reads the bcftools statstic file and generate a list of '
        'samples with a high missing rate to be removed from the SNP dataset.'
    )
)

parser.add_argument(
    "stats_file",
    help=(
        'The bcftools stats file. Please note the bcftools statistic file '
        'should be generated with the parameter "-s -", as the "PSC" filed in '
        'the stats file is necessary.'
    )
)

args = parser.parse_args()


high_missing_list = []
output_file = open(args.stats_file + '.missing_rate.list.txt', 'w')

with open(args.stats_file, 'r') as file:
    lines = file.readlines()
    SNP_num = 0
    for line in lines:
        if line.strip().startswith('SN') and 'number of SNPs' in line:
            SNP_num = line.strip().split('\t')[3]
            # print(line)
            print('SNP number: ', SNP_num)
            break

    for line in lines:
        if line.strip().startswith('PSC'):
            columns = line.strip().split('\t')
            sample_id, missing = columns[2], columns[13]
            missing_rate = float(missing) / float(SNP_num)
            percentage = "{:.2%}".format(missing_rate)
            output_file.write(sample_id + '\t' + percentage + '\n')

            if missing_rate * 100 >= float(threshold):
                high_missing_list.append(sample_id)
                print(sample_id, percentage)

output_file.close()


output_file = args.stats_file + '.higher' + \
    str(threshold) + '_missing_rate_sample.list.txt'

with open(output_file, 'w') as file:
    for s in high_missing_list:
        file.write(s + '\n')

end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time - start_time).center(79))
