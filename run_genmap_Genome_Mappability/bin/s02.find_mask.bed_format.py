#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.find_mask.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/07/18 21:35:05
# @Description:
#   This script is used to find the masked regions (low score) in the genome based on the mappability score.
# @Update: v1.1.0 2024-04-17 16:57:30
#   1. Output the pass regions (high score) in the genome.
#   2. Fix the bug that first window pass region is not print.
#   3. Fix the bug that the last window pass region is repeated.
#   4. Fix the issue that joint lines of genmap results using empty string, normally the result is single line for each chromosome.

# @Update: v2.0.0 2024-05-22 14:54:17
#    1. Update the output filename.

# @Update: v3.0.0 2024-05-28 21:54:23
#    1. Negative score is priority to be masked, i.e. check the final average score of the region that supposed to pass.

# @Update: v3.0.1 2024-05-30 00:49:01
#    1. Fix an issue that some continues regions are not merged.

# @Update: v3.0.2 2024-06-07 16:42:10
#    1. Fix the bug that happened in the first window not masked.

import datetime
import textwrap
import sys
import os
# import pandas as pd

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "3.0.2"
script_basename = "s02.find_mask"
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output')

# Parameters
step = 50000
window = 100000   # step * 2
threshold = 0.9   # window average mappability score


def warp(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Set default parameters
genmap_out_base = os.path.join(output_dir, "GWHBAOS00000000.genome.genmap")
genmap_txt = genmap_out_base+".txt"
genmap_size = genmap_out_base+".chrom.sizes"
mask = genmap_out_base+".mask.bed"
mask_with_score = genmap_out_base+".mask_with_score.bed"
mask_pass = genmap_out_base+".pass.bed"
mask_pass_with_score = genmap_out_base+".pass_with_score.bed"


if len(sys.argv) != 4:
    message = f'''Usage: python3 {sys.argv[0]} [ <.genmap_chrom.sizes> <.genmap.txt> <output_mask.bed> ]'''
    print(warp(message))

    if len(sys.argv) == 1:
        message = f'''Now running with default parameters: {sys.argv[0]} {genmap_size} {genmap_txt} {mask}'''
        print(warp(message))
    else:
        sys.exit(1)

elif len(sys.argv) == 4:
    genmap_size = sys.argv[1]
    genmap_txt = sys.argv[2]
    mask = sys.argv[3]

# Load genome size use pandas, or not. The file is small, so it's not necessary to use pandas.
# chr_size = pd.read_csv(genmap_size, sep='\t', header=None, names=['chr','size'])
chr_size = {}
with open(genmap_size, 'r') as fo:
    lines = fo.readlines()
    [chr_size.setdefault(line.split('\t')[0], line.split('\t')[1])
     for line in lines]


# Function: Read genmap results
def read_file(file):
    dic_score = {}

    with open(file, 'r') as fo:
        content = fo.read().strip().split('>')
        for sequences in content[1:]:
            lines = sequences.strip().split('\n')
            id = lines[0].split()[0]
            seq = ' '.join(lines[1:])
            dic_score[id] = [float(score) for score in seq.split()]
            # If use pandas, the code is:
            # dic_score[id] = pd.Series([float(score) for score in seq.split()])

    return dic_score


dic_score = read_file(genmap_txt)

# Function: get average score for a region


def get_average_score(score_list, start, end):
    return sum(score_list[start:end]) / (end-start)
    # Here if use pandas, the code is: but is slower than the above code.
    # return score_list[start:end].mean()


# Output
pass_region = open(mask_pass, 'w')
pass_with_score = open(mask_pass_with_score, 'w')
mask = open(mask, 'w')
with open(mask_with_score, 'w') as fo:
    for id, scores in dic_score.items():
        # If read the genome size using pandas:
        # length  = int(chr_size.loc[chr_size['chr'] == id, 'size'].values[0])
        length = int(chr_size[id])

        # # Xilong's strategy
        # for m in range(0,length-window+1, window):
        #     start1 = m
        #     end1 = m+window
        # write into file1

        # for n in range(step,length-window+1, window):
        #     start2 = n
        #     end2 = n+window
        # write into file2

        # sort and merge file1 and file2

        # # My strategy
        # for i in range(0, max(1,length-window+1), step):
        #     start = i
        #     end = min(i+window, length)
        #     # try:
        #     #     average_score = sum(scores[start:end]) / (end-start)
        #     # except ZeroDivisionError:
        #     #     average_score = 0

        #     # if average_score < thresholds:
        #     #     fo.write(f'{id}\t{start+1}\t{end}\t{average_score}\n')

        # Flowing strategy considers the overlap between two windows
        m_start = m_end = 0  # Masked region
        p_start = p_end = 0  # Passed region (unmask, or high score)

        for i in range(0, length, step):
            start = i
            end = min(i+window, length)
            average_score = get_average_score(scores, start, end)

            if average_score < threshold:
                p_end = min(start, p_end)
                if m_start == m_end == 0:
                    m_start = start
                    m_end = end
                    if p_start != p_end:
                        p_average_score = get_average_score(
                            scores, p_start, p_end)
                        if p_average_score < threshold:
                            m_start = min(start, p_start)
                            m_end = max(end, p_end)
                            # raise ValueError(
                            #     'BUG:The average score of the passed region is lower than the threshold.')
                        elif p_average_score >= threshold:
                            pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
                            pass_with_score.write(
                            f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
                            print(
                                f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')
                elif start <= m_end:
                    # m_start = m_start
                    m_end = end

                # elif start > m_end:
                #     m_average_score = get_average_score(scores, m_start, m_end)
                #     mask.write(f'{id}\t{m_start}\t{m_end}\n')
                #     fo.write(f'{id}\t{m_start}\t{m_end}\t{m_average_score}\n')
                #     p_average_score = get_average_score(scores, p_start, p_end)
                #     pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
                #     pass_with_score.write(
                #                     f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
                #     print(f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')
                #     if p_average_score < threshold:
                #         raise ValueError(
                #             'BUG[2]:The average score of the passed region is lower than the threshold.')
                elif start > m_end:
                    p_average_score = get_average_score(scores, p_start, p_end)
                    to_write_result = True
                    if p_average_score < threshold:
                        m_end = max(end, p_end)
                        to_write_result = False
                    elif p_average_score >= threshold:
                        pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
                        pass_with_score.write(
                            f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
                        print(
                            f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')
                        m_average_score = get_average_score(
                            scores, m_start, m_end)
                        mask.write(f'{id}\t{m_start}\t{m_end}\n')
                        fo.write(
                            f'{id}\t{m_start}\t{m_end}\t{m_average_score}\n')

                        m_start = start
                        m_end = end

            elif average_score >= threshold:
                p_start = max(p_start, m_end)
                p_end = end

            if end == length:
                break

        if p_start != p_end and m_end != end:
            p_average_score = get_average_score(scores, p_start, p_end)
            if p_average_score < threshold:
                m_end = max(m_end, p_end)
            elif p_average_score >= threshold:
                pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
                pass_with_score.write(
                    f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
                print(
                    f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')

        if m_start != m_end:
            m_average_score = get_average_score(scores, m_start, m_end)
            mask.write(f'{id}\t{m_start}\t{m_end}\n')
            fo.write(f'{id}\t{m_start}\t{m_end}\t{m_average_score}\n')

mask.close()
pass_with_score.close()
pass_region.close()


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
