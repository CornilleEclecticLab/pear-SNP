#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s05_genetic_burden_by_roh.py
# @Author   : NIE Yuqi
# @Time(Europe/Paris): 2026/07/30 14:10:24
# @Description:
#    



import datetime
import sys
import textwrap
import os
import glob
import gzip
import bisect
from s05_config import *

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
# os.makedirs(sub_output_dir, exist_ok=True)


REGIONS = ('ROH', 'nonROH')
ZYGOSITIES = ('Allele', 'Homozygous', 'Heterozygous')
CLASSES = ('DELETERIOUS', 'TOLERATED', 'SYNONYMOUS')


def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' ' * 4)


# 1. ROH runs -> per-sample, per-chrom searchable index                       #
def build_roh_index(summary_files):
    """
    Parse per-sample *.summary.tsv (merged ROH runs, short + long) into a fast
    membership index.

    Returns: {sample_id: {chrom: (starts, ends)}}, {sample_id: roh_span_bp}
      starts/ends are parallel, sorted-by-start lists of non-overlapping runs,
      enabling O(log n) lookup with bisect. v3's merge step guarantees runs are
      non-overlapping and adjacent-merged, so a single bisect is sufficient.
    """
    roh_index = {}
    roh_span_bp = {}

    for filepath in summary_files:
        # "<sample>.mid.hmmp.gz.summary.tsv".
        sample_id = os.path.basename(filepath)
        for suffix in ('.mid.hmmp.gz.summary.tsv', '.summary.tsv'):
            if sample_id.endswith(suffix):
                sample_id = sample_id[:-len(suffix)]
                break

        by_chrom = {}          # {chrom: [(start, end), ...]}
        span = 0
        with open(filepath, 'r') as fi:
            for line in fi:
                if line.startswith('#') or not line.strip():
                    continue
                # columns: CHROM BEGIN END p[ROH] p[nonROH] Length type
                parts = line.rstrip('\n').split('\t')
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                by_chrom.setdefault(chrom, []).append((start, end))
                span += end - start

        # Sort each chrom's runs and split into parallel start/end lists.
        sample_index = {}
        for chrom, runs in by_chrom.items():
            runs.sort()
            starts = [r[0] for r in runs]
            ends = [r[1] for r in runs]
            # confirm non-overlap.
            for i in range(1, len(starts)):
                if starts[i] < ends[i - 1]:
                    raise ValueError(
                        f'Overlapping ROH runs for {sample_id} on {chrom}: '
                        f'[{starts[i-1]},{ends[i-1]}) & [{starts[i]},{ends[i]}). '
                        f'Expected merged, non-overlapping runs.')
            sample_index[chrom] = (starts, ends)

        if sample_id in roh_index:
            raise ValueError(f'Duplicate ROH summary for sample {sample_id}.')
        roh_index[sample_id] = sample_index
        roh_span_bp[sample_id] = span

    return roh_index, roh_span_bp


def in_roh(sample_index, chrom, pos_1based):
    """True if 1-based VCF position falls inside one of the sample's ROH runs."""
    chrom_runs = sample_index.get(chrom)
    if not chrom_runs:
        return False
    starts, ends = chrom_runs
    p0 = pos_1based - 1                         # to 0-based coordinate
    i = bisect.bisect_right(starts, p0) - 1     # rightmost run with start <= p0
    return i >= 0 and p0 < ends[i]              # inside [start, end)?


# SIFT4G annotation reader
def read_sift4g_annotation(file_path):
    annotations = {}
    with open(file_path, 'r') as file:
        header = None
        for line in file:
            line = line.strip()
            parts = line.split('\t')
            if header is None and line.startswith('CHROM'):
                header = parts
                continue

            chr = parts[0]
            pos = parts[1]
            ref = parts[2]
            alt = parts[3]
            variant_type = parts[8]
            sift_score = parts[12]
            prediction = parts[16]

            if variant_type not in ['NONSYNONYMOUS', 'SYNONYMOUS']:
                continue
            if variant_type == 'NONSYNONYMOUS' and sift_score == 'NA':
                continue
            if variant_type == 'SYNONYMOUS':
                prediction = 'SYNONYMOUS'
                sift_score = 0.0

            annotations.setdefault(chr, {})
            annotations[chr][pos] = {
                'variant_type': variant_type,
                'sift_score': float(sift_score),
                'prediction': prediction,
                'ref': ref,
                'alt': alt,
            }
    return annotations


# Genotype parsing
def parse_vcf_genotype(gt_index, genotype):
    gt = genotype.split(':')[gt_index]
    gt = gt.replace('/', '|').replace('\\', '|').split('|')

    if len(gt) != 2:
        raise ValueError(f'Invalid genotype format: {genotype}')

    sum1_m = sum(1 for g in gt if g == '1')
    sum0_m = sum(1 for g in gt if g == '0')

    if '.' in gt and sum1_m + sum0_m == 0:
        return None
    elif '.' in gt and sum1_m + sum0_m == 1:
        raise ValueError(f'Invalid genotype format with missing data: {genotype}')

    if sum1_m == 2:
        return 'Homozygous'
    elif sum1_m == 1:
        return 'Heterozygous'
    elif sum1_m == 0 and sum0_m == 2:
        return 'Non-variant'
    else:
        raise ValueError(f'Invalid genotype format: {genotype}')


def new_sample_counter():
    """
    Per-sample counter, now indexed by REGION first.
      counts[region][zygosity][class]   -> alt-allele / genotype counts
      sites[region][class]              -> # genotyped (non-missing) annotated
                                           sites, i.e. the normalisation
                                           denominator per region.
    """
    return {
        'counts': {r: {z: {c: 0 for c in CLASSES} for z in ZYGOSITIES}
                   for r in REGIONS},
        'sites': {r: {c: 0 for c in CLASSES} for r in REGIONS},
    }


# VCF parsing, now ROH-aware
def parse_vcf(file_path, sift_annotations, roh_index):
    variants = {}

    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt') as file:
        header = None
        for line in file:
            line = line.strip()
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line.lstrip('#').strip().split('\t')
                continue

            parts = line.split('\t')
            chr, pos, id, ref, alt, qual, filter, info, format, *genotypes = parts

            if pos not in sift_annotations.get(chr, {}):
                continue

            ann = sift_annotations[chr][pos]
            prediction = ann['prediction']

            if ref != ann['ref'] or alt != ann['alt']:
                raise ValueError(
                    f'Variant {chr}:{pos} has different ref/alt in VCF and '
                    f'SIFT4G annotation: VCF({ref}/{alt}) vs '
                    f'SIFT4G({ann["ref"]}/{ann["alt"]})')

            gt_index = format.split(':').index('GT')
            pos_int = int(pos)

            for i, sample_id in enumerate(header[9:]):
                gt = parse_vcf_genotype(gt_index, genotypes[i])
                if gt is None:                     # missing -> not a callable site
                    continue

                variants.setdefault(sample_id, new_sample_counter())

                # Per-individual ROH membership.
                sidx = roh_index.get(sample_id)
                if sidx is None:
                    raise ValueError(
                        f'No ROH summary found for sample {sample_id}; cannot '
                        f'assign ROH/non-ROH. Check roh_summary_dir_dic.')
                region = 'ROH' if in_roh(sidx, chr, pos_int) else 'nonROH'

                # Normalisation denominator: this is a genotyped annotated site
                # for this individual in this region, regardless of zygosity.
                variants[sample_id]['sites'][region][prediction] += 1

                c = variants[sample_id]['counts'][region]
                if gt == 'Homozygous':
                    c['Allele'][prediction] += 2
                    c['Homozygous'][prediction] += 2
                elif gt == 'Heterozygous':
                    c['Allele'][prediction] += 1
                    c['Heterozygous'][prediction] += 1
                # 'Non-variant' contributes only to the site denominator.

    return variants


# Sample/pop map (mirrors s04 structure)
def read_individual_pop_map(file_path):
    individual_pop_map = {}
    pop_individual = {}
    with open(file_path, 'r') as file:
        for line in file:
            ind, pop = line.strip().split()[0:2]
            individual_pop_map[ind] = pop
            pop_individual.setdefault(pop, []).append(ind)
    return individual_pop_map, pop_individual


def merge_sample_counters(dst, src):
    """Accumulate src (one VCF's parse result) into dst across chrom/pop VCFs."""
    for region in REGIONS:
        for cls in CLASSES:
            dst['sites'][region][cls] += src['sites'][region][cls]
        for zyg in ZYGOSITIES:
            for cls in CLASSES:
                dst['counts'][region][zyg][cls] += src['counts'][region][zyg][cls]


def process_vcf_files(ref_genome, sift_annotations, roh_index):
    burden_lines = []
    denom_lines = []
    sample_variants = {}

    ind_pop_dic, pop_ind_dic = read_individual_pop_map(
        individuals_map_file_dic[ref_genome])

    with open(chr_list_file_dic[ref_genome], 'r') as file:
        chr_list = [line.strip() for line in file if line.strip()]

    processed = {chr: {pop: False for pop in pop_ind_dic} for chr in chr_list}
    total_vcf_num = len(chr_list) * len(pop_ind_dic)

    vcf_file_list = glob.glob(
        os.path.join(vcf_file_dir_dic[ref_genome], '*.vcf.gz'))

    n_vcf = 0
    for vcf_file in vcf_file_list:
        pop_vcf = os.path.basename(vcf_file).split('.')[-3]
        chr_vcf = os.path.basename(vcf_file).split('.')[-4]

        if pop_vcf not in pop_ind_dic or chr_vcf not in chr_list:
            continue
        if processed[chr_vcf][pop_vcf]:
            raise ValueError(f'VCF for {chr_vcf}/{pop_vcf} already processed.')

        n_vcf += 1
        print(f'Processing {n_vcf}/{total_vcf_num}: {vcf_file}')
        s_variant = parse_vcf(vcf_file, sift_annotations, roh_index)

        for sample_id, counter in s_variant.items():
            if pop_vcf != ind_pop_dic[sample_id]:
                raise ValueError(
                    f'Sample {sample_id} in {vcf_file} != population {pop_vcf}.')
            if sample_id not in sample_variants:
                sample_variants[sample_id] = counter
            else:
                merge_sample_counters(sample_variants[sample_id], counter)

        processed[chr_vcf][pop_vcf] = True

    for chr, pop_dict in processed.items():
        for pop, done in pop_dict.items():
            if not done:
                raise ValueError(f'No VCF found for {chr} in population {pop}.')

    print(f'{len(sample_variants)} samples aggregated.\n')

    for sample_id, v in sample_variants.items():
        pop = ind_pop_dic[sample_id]
        for region in REGIONS:
            for zyg in ZYGOSITIES:
                cc = v['counts'][region][zyg]
                burden_lines.append(
                    f"{sample_id}\t{pop}\t{region}\t{zyg}\t"
                    f"{cc['DELETERIOUS']}\t{cc['TOLERATED']}\t{cc['SYNONYMOUS']}\n")
            sc = v['sites'][region]
            denom_lines.append(
                f"{sample_id}\t{pop}\t{region}\t"
                f"{sc['DELETERIOUS']}\t{sc['TOLERATED']}\t{sc['SYNONYMOUS']}\t"
                f"{roh_span_bp_global.get(sample_id, 'NA')}\n")

    return burden_lines, denom_lines


def main():
    global roh_span_bp_global

    burden_header = ('Individual\tPopulation\tRegion\tType\t'
                     'Deleterious_allele_count\tTolerated_allele_count\t'
                     'Synonymous_allele_count\n')
    denom_header = ('Individual\tPopulation\tRegion\t'
                    'Genotyped_deleterious_sites\tGenotyped_tolerated_sites\t'
                    'Genotyped_synonymous_sites\tROH_span_bp\n')

    concat_burden = open(
        os.path.join(sub_output_dir, 'pop_mutation_burden_by_ROH.txt'), 'w')
    concat_denom = open(
        os.path.join(sub_output_dir, 'pop_ROH_denominators.txt'), 'w')
    concat_burden.write(burden_header)
    concat_denom.write(denom_header)

    for ref_genome in ref_genomes:
        print(f'Processing {ref_genome}...\n')

        sift_annotations = read_sift4g_annotation(
            os.path.join(sift4g_dir_dic[ref_genome],
                         'merged.SIFTannotations.xls.cleaned'))
        print(f'Read {sum(len(v) for v in sift_annotations.values())} '
              f'SIFT4G annotations for {ref_genome}.\n')

        summary_files = glob.glob(
            os.path.join(roh_summary_dir_dic[ref_genome], '*', '*.summary.tsv'))
        if not summary_files:
            raise FileNotFoundError(
                f'No *.summary.tsv under {roh_summary_dir_dic[ref_genome]}. '
                f'Run s02.x1.get_ROH_table.py first.')
        roh_index, roh_span_bp_global = build_roh_index(summary_files)
        print(f'Loaded ROH runs for {len(roh_index)} samples.\n')

        burden_lines, denom_lines = process_vcf_files(
            ref_genome, sift_annotations, roh_index)

        with open(os.path.join(sub_output_dir,
                               ref_genome + '_mutation_burden_by_ROH.txt'), 'w') as f:
            f.write(burden_header)
            f.writelines(burden_lines)
        with open(os.path.join(sub_output_dir,
                               ref_genome + '_ROH_denominators.txt'), 'w') as f:
            f.write(denom_header)
            f.writelines(denom_lines)

        concat_burden.writelines(burden_lines)
        concat_denom.writelines(denom_lines)

    concat_burden.close()
    concat_denom.close()


if __name__ == "__main__":
    require_paths()
    global roh_span_bp_global
    roh_span_bp_global = {}
    main()

end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time - start_time).center(79))
