#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @Description:
#    Sum s12.x1 per-chromosome counts per sample, and write one burden table per arm

import glob
import os
import sys

out_dir = '../output/s12.downsample_AD'
ref_genomes = ['PPY', 'PCOM']
n_chr = 17

pop = {}
for ref in ref_genomes:
    with open(f'../input/s01.sample_tab_population_{ref}.txt') as f:
        for line in f:
            sid, p = line.split()[0:2]
            pop[sid] = p

files = sorted(glob.glob(os.path.join(out_dir, 'counts', '*.counts.tsv')))
if len(files) != n_chr * len(ref_genomes):
    sys.exit(f'Expected {n_chr * len(ref_genomes)} count files, found {len(files)}')

# totals[arm][sample] = [del_all, del_ho, del_he, syn_site]
totals = {}
for fp in files:
    with open(fp) as f:
        next(f)
        for line in f:
            arm, _, sid, *vals = line.rstrip('\n').split('\t')
            t = totals.setdefault(arm, {}).setdefault(sid, [0, 0, 0, 0])
            for i, v in enumerate(vals):
                t[i] += int(v)

header = ('Sample_ID\tPopulation\tAllelic_state\tDeleterious_allele_count\t'
          'Deleterious_per_synonymous_SNP\tRegion_type\tSynonymous_SNP_count\n')
for arm, per_sample in totals.items():
    out = os.path.join(out_dir, f'{arm}.synonymous.pop_mutation_burden.Rdata.txt')
    with open(out, 'w') as fo:
        fo.write(header)
        for sid, (d_all, d_ho, d_he, syn) in per_sample.items():
            for vt, de in (('All', d_all), ('Ho', d_ho), ('He', d_he)):
                ratio = de / syn if de and syn else 0.0
                fo.write(f'{sid}\t{pop[sid]}\t{vt}\t{de}\t{ratio}\tall_masked\t{syn}\n')
    print(f'{arm}: {len(per_sample)} samples -> {out}')

# Check
s06 = {}
with open('../output/s06.dm_sweep/region.synonymous.pop_mutation_burden.Rdata.txt') as f:
    next(f)
    for line in f:
        sid, _, vt, de, _, region, syn = line.rstrip('\n').split('\t')
        if region == 'all_masked':
            s06[(sid, vt)] = (int(de), int(syn))
n_diff = 0
for sid, (d_all, d_ho, d_he, syn) in totals['orig'].items():
    for vt, de in (('All', d_all), ('Ho', d_ho), ('He', d_he)):
        if s06.get((sid, vt)) != (de, syn):
            n_diff += 1
            if n_diff <= 5:
                print(f'DIFF {sid} {vt}: s12 {(de, syn)} vs s06 {s06.get((sid, vt))}')
print(f'orig vs s06 all_masked: {n_diff} differing rows out of {3 * len(totals["orig"])}')
