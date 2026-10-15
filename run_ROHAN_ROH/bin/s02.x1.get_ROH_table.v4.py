#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.x1.get_short_long_ROH.py
# @Author   : NIE Yuqi
# @Time(Europe/Paris): 2026/06/05 18:02:09
# @Description:
#  Read *.mid.hmmp.gz from ROHan output,
#  merge contiguous ROH windows into runs,
# then classify  each run as short or long ROH in bed format.


import sys
import gzip

try:
    hmmp_file = sys.argv[1]
except:
    print(f'Usage: python {sys.argv[0]} <hmmp_file> [<short_length(bp) default:1Mb> <p_threshold default:0.99>]')
    sys.exit(1)

try:
    short_threshold = int(sys.argv[2])
    p_threshold = float(sys.argv[3])
except:
    short_threshold = 1000_000 # 1Mb
    p_threshold = 0.99


with gzip.open(hmmp_file, 'rt') as f, \
    open(hmmp_file+'.long.bed', 'w') as fl, \
    open(hmmp_file+'.short.bed', 'w') as fs, \
    open(hmmp_file+'.summary.tsv', 'w') as fsum:

    def flush(run):
        if run is None:
            return
        chrom, start, end, proh, pnonroh = run
        length = end - start
        if length <= short_threshold:
            long_or_short_type = "short"
            fs.write(f"{chrom}\t{start}\t{end}\t{proh}\t{pnonroh}\n")
        else:
            long_or_short_type = "long"
            fl.write(f"{chrom}\t{start}\t{end}\t{proh}\t{pnonroh}\n")
        fsum.write(f"{chrom}\t{start}\t{end}\t{proh}\t{pnonroh}\t{length}\t{long_or_short_type}\n")

    cur = None

    header = None
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            if header is None:
                header = line  #  "#CHROM  BEGIN   END     p[ROH]  p[nonROH]"
                fsum.write(header+"\tLength\tLong_or_short_ROH_type\n")
                if header!="#CHROM\tBEGIN\tEND\tp[ROH]\tp[nonROH]":
                    print(header)
                    raise ValueError("Header of the input file is unexpeceted.")
            else:
                raise ValueError("Multiple headers found.")
            continue

        chrom,start,end,proh,pnonroh = line.split()

        # NA for un-assessed windows (no data / low mappability).
        if proh == "NA" or pnonroh == "NA":
            flush(cur); cur = None
            continue

        start = int(start)
        end = int(end)
        proh = float(proh)
        pnonroh = float(pnonroh)

        # instead of writing each ROH window on its own, grow
        # a run across adjacent ROH windows and only write it on flush.
        if not (proh > p_threshold and proh > pnonroh):
            flush(cur); cur = None
            continue

        if cur is not None and chrom == cur[0] and start == cur[2]:
            cur[2] = end
            cur[3] = min(cur[3], proh)
            cur[4] = max(cur[4], pnonroh)
        else:
            flush(cur)
            cur = [chrom, start, end, proh, pnonroh]

    flush(cur)