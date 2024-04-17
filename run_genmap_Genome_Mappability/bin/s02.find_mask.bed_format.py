#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.find_mask.py
# @Version  : 1.1.0
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



import datetime
import textwrap
import sys
import os
import pandas as pd

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "1.1.0"
script_basename = "s02.find_mask"
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output')

# Parameters
step = 50000
window = 100000 # step * 2
thresholds = 0.9 # 90% of the window should be mapped to the genome


def warp(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Set default parameters
genmap_out_base = os.path.join(output_dir, "GWHBAOS00000000.genome.genmap")
genmap_txt = genmap_out_base+".txt"
genmap_size = genmap_out_base+".chrom.sizes"
mask = genmap_out_base+".mask.bed"
mask_with_score = genmap_out_base+".mask_with_score.bed"


if len(sys.argv) !=4:
    message=f'''Usage: python3 {sys.argv[0]} [ <.genmap_chrom.sizes> <.genmap.txt> <output_mask.bed> ]'''
    print(warp(message))

    if len(sys.argv) == 1:
        message=f'''Now running with default parameters: {sys.argv[0]} {genmap_size} {genmap_txt} {mask}'''
        print(warp(message))
    else:
        sys.exit(1)

elif len(sys.argv) == 4:
    genmap_size = sys.argv[1]
    genmap_txt = sys.argv[2]
    mask = sys.argv[3]


chr_size = pd.read_csv(genmap_size, sep='\t', header=None, names=['chr','size'])

def read_file(file): # read a genmap.txt file
    dic_score = {}

    with open(file,'r') as fo:
        content=fo.read().strip().split('>')
        for sequences in content[1:]:
            lines=sequences.strip().split('\n')
            id = lines[0].split()[0]
            seq = ' '.join(lines[1:])
            dic_score[id] = [float(score) for score in seq.split()]

    return dic_score

dic_score = read_file(genmap_txt)

def get_average_score(score_list, start, end):
    average_score = sum(score_list[start:end]) / (end-start)
    return average_score


pass_region = open(mask+'.pass.bed','w')
pass_with_score = open(mask+'.pass_with_score.bed','w')
mask = open(mask, 'w')
with open(mask_with_score,'w') as fo:
    for id, scores in dic_score.items():

        length  = int(chr_size.loc[chr_size['chr'] == id, 'size'].values[0])

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
        # for i in range(0,max(1,length-window+1), step):
        for i in range(0,length, step):
            start = i
            end = min(i+window, length)
            average_score = get_average_score(scores, start, end)

            if average_score < thresholds:
                p_end = min(start, p_end)
                if m_start == m_end == 0:
                    m_start = start
                    m_end = end

                    if p_start != p_end:
                        p_average_score = get_average_score(scores, p_start, p_end)
                        pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
                        pass_with_score.write(f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
                        print(f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')

                elif start <= m_end:
                    m_start = m_start
                    m_end = end

                elif start > m_end:
                    m_average_score = get_average_score(scores, m_start, m_end)
                    mask.write(f'{id}\t{m_start}\t{m_end}\n')
                    fo.write(f'{id}\t{m_start}\t{m_end}\t{m_average_score}\n')

                    p_average_score = get_average_score(scores, p_start, p_end)
                    pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
                    pass_with_score.write(
                                    f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
                    print(f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')

                    m_start = start
                    m_end = end

            elif average_score >= thresholds:
                    p_start = max(p_start, m_end) if m_end != 0 else p_start
                    p_end = end

            if end == length:
                break

        if m_start != m_end:
            m_average_score = get_average_score(scores, m_start, m_end)
            mask.write(f'{id}\t{m_start}\t{m_end}\n')
            fo.write(f'{id}\t{m_start}\t{m_end}\t{m_average_score}\n')


        if p_start != p_end and m_end != end:
            p_average_score = get_average_score(scores, p_start, p_end)
            pass_region.write(f'{id}\t{p_start}\t{p_end}\n')
            pass_with_score.write(f'{id}\t{p_start}\t{p_end}\t{p_average_score}\n')
            print(f'Passed region:  {id}\t{p_start}\t{p_end}\t{p_average_score}')

mask.close()
pass_with_score.close()
pass_region.close()


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
