#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s06.x1.dm_extract_sweep.py
# @Author   : NIE Yuqi (based on CHEN Xilong's methodology)
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/06/06 12:36:25
# @Description:
#

from s04_config import *
from s06_config import *
import datetime
import sys
import textwrap
import os
import glob
import argparse
import subprocess
try:
    from intervaltree import IntervalTree
    USE_INTERVALTREE = True
except ImportError:
    USE_INTERVALTREE = False
    print("TIP: Active or install 'intervaltree' to speed up interval operations.")

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir, 'input')
output_dir = os.path.join(work_dir, 'output')

temp_dir = os.path.join(work_dir, 'temp')

# Ensure output directory exists
os.makedirs(temp_dir, exist_ok=True)


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Python3 version of dm_extract_sweep.pl")
    parser.add_argument("output_dir", help="Output directory (e.g. ../output/s06.dm_sweep)")
    parser.add_argument(
        "region_type", help="Region type (e.g. selection, control, all_masked)")
    parser.add_argument("bed_files", nargs="*", help="One or more BED files named by *.*.pop.*.bed")
    return parser.parse_args()


def intersect_and_calc_cds_length(cds_bed, region_bed, tmp_bedfile):
    """
    Use bedtools intersect to calculate the intersection of cds_bed and region_bed, output to tmp_bedfile.
    Then use awk to calculate the total length of all intersected segments, and return this length (int).
    """
    # bedtools intersect
    load_bedtools = 'source s00.load_bedtools.sh'
    if args.region_type == 'control':
        cmd_int = f"bedtools intersect -nonamecheck -v -a {cds_bed} -b {region_bed}"
    else:
        cmd_int = f"bedtools intersect -nonamecheck -wa -a {cds_bed} -b {region_bed}"
    with open(tmp_bedfile, "w") as outf:
        subprocess.run(f"{load_bedtools} && {cmd_int}", stdout=outf, shell=True, check=True, executable="/bin/bash")
    # awk to calculate length
    # awk '{sum += $3-$2}END{print sum}' tmp_bedfile
    cmd_awk = ["awk", "{sum += $3-$2}END{print sum}", tmp_bedfile]
    result = subprocess.run(cmd_awk, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    total_len_str = result.stdout.strip()
    try:
        return int(total_len_str)
    except ValueError:
        return 0


def build_bed_intervals(tmp_bedfile):
    """
    Read tmp_bedfile, each line has 3 columns: chr, start, end.
    Return a dict: bed_dic = { "Chr01": [(start1, end1), (start2,end2), ...], ... }
                 or tree_dic = {"Chr01: IntervalTree, ...}
    """

    bed_dic = {}
    tree_dic = {}

    with open(tmp_bedfile) as fin:
        for line in fin:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            s = int(cols[1])
            e = int(cols[2])
            bed_dic.setdefault(chrom, []).append((s, e))
            if USE_INTERVALTREE:
                tree_dic.setdefault(chrom, IntervalTree())
                tree_dic[chrom].addi(s, e)
    if USE_INTERVALTREE:
        return tree_dic
    return bed_dic


def position_in_intervals(pos, intervals):
    # pos is 1 based index from VCF
    # intervals is a list of tuples (start, end) or an IntervalTree.
    # BED rule: 0-based start, 1-based end.
    
    zpos = int(pos) - 1  # Convert to 0-based index
    
    if USE_INTERVALTREE:
        return bool(intervals[zpos]) # the intervals include the lower bound but not the upper bound [start, end)

    else:
        for (s, e) in intervals:
            if zpos >= int(s) and zpos < int(e):
                return True
        return False


class VCFParser:
    """
    Parse VCF files, including:
      - parse_genotype(): parse non-variant, homozygous, heterozygous from GT field.
      - parse_vcf():parse a VCF file, according to sift_annotations, and bed-region.
    """

    def __init__(self, sift_annotations: dict):
        # sift_annotations is returned by read_sift4g_annotation()
        self.sift_annotations = sift_annotations

    def parse_genotype(self, gt_index: int, genotype: str): # -> str or None
        """
        Split the genotype string ("0/1:..."), locate the GT field by gt_index,
        determine if it is 'Homozygous', 'Heterozygous', 'Non-variant', or return None (if missing ./. ).
        """
        # Extract GT field
        gt = genotype.split(':')[gt_index]
        alleles = gt.replace('/', '|').replace('\\', '|').split('|')

        if len(alleles) != 2:
            raise ValueError(f'Invalid genotype format: {genotype}')

        # Count occurrences of 0 and 1
        sum1 = sum(1 for a in alleles if a == '1')
        sum0 = sum(1 for a in alleles if a == '0')

        # In case of missing
        if '.' in alleles and sum1 + sum0 == 0:
            return None
        # If only one missing site, it maybe a format error
        if '.' in alleles and sum1 + sum0 == 1:
            raise ValueError(
                f'Invalid genotype format with missing data: {genotype}')

        if sum1 == 2:
            return 'Homozygous'
        elif sum1 == 1:
            return 'Heterozygous'
        elif sum1 == 0 and sum0 == 2:
            return 'Non-variant'
        else:
            raise ValueError(f'Invalid genotype format: {genotype}')

    def read_sift4g_annotation(self, chrom, pos, ref, alt):
        # If there is no corresponding chrom:pos in sift_annotations, skip
        if chrom not in self.sift_annotations or pos not in self.sift_annotations[chrom]:
            return None

        ann = self.sift_annotations[chrom][pos]
        prediction = ann['prediction']  # 'DELETERIOUS' or 'TOLERATED'

        # Check if ref/alt match
        if ref != ann['ref'] or alt != ann['alt']:
            raise ValueError(
                f'Variant {chrom}:{pos} in VCF({ref}/{alt}) and SIFT4G annotation ({ann["ref"]}/{ann["alt"]}) do not match.'
            )
        return ann

    def parse_vcf(self, file_path: str, bed_region: dict) -> dict:
        """
        Iterate through each line in a single VCF file. If there is a corresponding (Chr, Pos) record in sift_annotations,
        then for each sample, count the DEL / TOL number of Allele / Homozygous / Heterozygous.

        Return value:
          variants = {
             sample_id1: {
               'Allele': {'DELETERIOUS': int, 'TOLERATED': int},
               'Homozygous': {...},
               'Heterozygous': {...}
             },
             sample_id2: { ... }
          }
        """
        variants = {}
        with open(file_path, 'r') as fin:
            for line in fin:
                line = line.strip()
                # Ignore meta header
                if line.startswith('##'):
                    continue
                # Process header line to get sample ID column names
                if line.startswith('#CHROM'):
                    header = line.lstrip('#').split('\t')
                    sample_ids = header[9:]  # From the 10th column, each column is a sample
                    continue

                parts = line.split('\t')
                chrom = parts[0]
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                format_field = parts[8]
                genotypes = parts[9:]

                # check if the position is in the bed_region
                if chrom not in bed_region.keys():
                    continue

                if not position_in_intervals(int(pos), bed_region[chrom]):
                    continue

                # Read SIFT4G annotation for this variant
                ann = self.read_sift4g_annotation(chrom, pos, ref, alt)
                if ann is None:
                    continue
                prediction = ann['prediction']  # 'DELETERIOUS' or 'TOLERATED'

                # Locate which column GT is in the FORMAT field
                fmt_cols = format_field.split(':')
                gt_index = fmt_cols.index('GT')

                # For each sample, parse genotype and accumulate
                for idx, sample_id in enumerate(sample_ids):
                    gt = self.parse_genotype(gt_index, genotypes[idx])
                    # If gt is None (missing), do nothing
                    if gt is None:
                        continue

                    # Initialize sample statistics structure
                    variants.setdefault(sample_id, {
                        'Allele': {'DELETERIOUS': 0, 'TOLERATED': 0},
                        'Homozygous': {'DELETERIOUS': 0, 'TOLERATED': 0},
                        'Heterozygous': {'DELETERIOUS': 0, 'TOLERATED': 0}
                    })

                    if gt == 'Homozygous':
                        # Homozygous alt -> Alt allele counts as 2
                        variants[sample_id]['Allele'][prediction] += 2
                        variants[sample_id]['Homozygous'][prediction] += 2
                    elif gt == 'Heterozygous':
                        variants[sample_id]['Allele'][prediction] += 1
                        variants[sample_id]['Heterozygous'][prediction] += 1
                    # "Non-variant" does not need to accumulate

        return variants


# Function
def read_sift4g_annotation(file_path):
    # Arg file_path is the path to the SIFT4G annotation file (.xls or .xls.cleaned)
    # Return the following dictionary
    annotations = {} # {chrom: {pos:
    #           {'type': str, 'score': float, 'prediction': str}}}
    with open(file_path, 'r') as file:
        header = None
        for line in file:
            line = line.strip()
            parts = line.split('\t')
            if header is None and line.startswith('CHROM'):
                header = parts
                continue

            chrom = parts[0]
            pos = parts[1]
            ref = parts[2]
            alt = parts[3]
            variant_type = parts[8]  # e.g. SYNONYMOUS, NONSYNONYMOUS, START-LOST, STOP-GAIN, STOP-LOSS
            sift_score = parts[12]
            prediction = parts[16]

            if variant_type != 'NONSYNONYMOUS' or sift_score == 'NA':
                continue

            annotations.setdefault(chrom, {})
            annotations[chrom][pos] = {
                'variant_type' : variant_type,
                'sift_score' : float(sift_score),
                'prediction' : prediction,
                'ref' : ref,
                'alt' : alt
            }
    return annotations


# Read individuals pop map
def read_individual_pop_map(file_path):
    # The file is expected to have two columns wihtout header: individual_id and population_id
    individual_pop_map = {}  # {ind_id: pop_id}
    pop_individual = {}  # {pop_id: [ind_id1, ind_id2, ...]}
    with open(file_path, 'r') as file:
        for line in file:
            ind, pop = line.strip().split()[0:2]
            individual_pop_map[ind] = pop
            pop_individual[pop] = pop_individual.get(pop, [])
            pop_individual[pop].append(ind)

    return individual_pop_map, pop_individual

# Function to process VCF files and check for missing populations
def process_vcf_files(ref_genome, sift_annotations, cds_length_dic, bed_region_per_pop):

    # Return str lines to be written to output file
    write_lines = []
    sample_stats = {}

    ind_pop_dic, pop_ind_dic = read_individual_pop_map(
        individuals_map_file_dic[ref_genome])

    chr_list_file = chr_list_file_dic[ref_genome]

    with open(chr_list_file, 'r') as file:
        chr_list = [line.strip() for line in file if line.strip()]

    processed_vcf_dic = {chrom: {pop: False for pop in pop_ind_dic.keys()} for chrom in chr_list}
    total_vcf_num = len(chr_list) * len(pop_ind_dic.keys())

    vcf_parser = VCFParser(sift_annotations)

    vcf_file_dir = vcf_file_dir_dic[ref_genome]
    vcf_file_list = glob.glob(os.path.join(vcf_file_dir, '*.vcf'))
    n_vcf = 0

    for vcf_file in vcf_file_list:
        if not vcf_file.endswith('.vcf'):
            continue
        basename = os.path.basename(vcf_file)
        pop_vcf = basename.split('.')[-2]   # VCF filename format: *.chr.pop.vcf
        chr_vcf = basename.split('.')[-3]

        if pop_vcf not in pop_ind_dic.keys() or chr_vcf not in chr_list:
            # print(
            #     f'Skipping {vcf_file} due to mismatch in population {pop_vcf} or chromosome {chr_vcf}.\n')
            continue

        if processed_vcf_dic[chr_vcf][pop_vcf]:
            raise ValueError(
                f'VCF file {vcf_file} for chromosome {chr_vcf} and population {pop_vcf} has already been processed.')

        n_vcf += 1
        print(f'Processing {n_vcf}/{total_vcf_num} VCF file: {vcf_file} for population {pop_vcf} and chromosome {chr_vcf}...\n')
        
        if pop_vcf not in bed_region_per_pop or chr_vcf not in bed_region_per_pop[pop_vcf]:
            processed_vcf_dic[chr_vcf][pop_vcf] = True
            continue
        
        s_variant= vcf_parser.parse_vcf(vcf_file, bed_region=bed_region_per_pop[pop_vcf])
        print(f'Parsed {len(s_variant)} samples from VCF file {vcf_file}.\n')

        for sample_id, variants in s_variant.items():
            if pop_vcf != ind_pop_dic[sample_id]:
                raise ValueError(
                    f'Sample {sample_id} in VCF file {vcf_file} does not match population {pop_vcf}.')

            if sample_id not in sample_stats:
                sample_stats[sample_id] = variants
                continue
            sample_stats[sample_id]['Allele']['DELETERIOUS'] += variants['Allele']['DELETERIOUS']
            sample_stats[sample_id]['Homozygous']['DELETERIOUS'] += variants['Homozygous']['DELETERIOUS']
            sample_stats[sample_id]['Heterozygous']['DELETERIOUS'] += variants['Heterozygous']['DELETERIOUS']
            sample_stats[sample_id]['Allele']['TOLERATED'] += variants['Allele']['TOLERATED']
            sample_stats[sample_id]['Homozygous']['TOLERATED'] += variants['Homozygous']['TOLERATED']
            sample_stats[sample_id]['Heterozygous']['TOLERATED'] += variants['Heterozygous']['TOLERATED']

        processed_vcf_dic[chr_vcf][pop_vcf] = True

    for chrom, pop_dict in processed_vcf_dic.items():
        for pop, processed in pop_dict.items():
            if not processed:
                raise ValueError(f'No VCF file found for {chrom} in population {pop}.')

    print(len(sample_stats), 'samples found in all VCF files for all populations and chromosomes.\n')
    for sample_id, variants in sample_stats.items():
        pop = ind_pop_dic[sample_id]
        cds_len = cds_length_dic[pop]
        if cds_len == 0:
            raise ValueError(
                f'Population {pop} CDS length is 0, please check the BED files.')
        for var_type, vt in (['Allele', 'All'],  ['Homozygous', 'Ho'], ['Heterozygous', 'He']):
            de_nu = variants[var_type]['DELETERIOUS']
            de_per_cds = de_nu / cds_len if de_nu else 0.0
            # tol_nu = variants[var_type]['TOLERATED']
            # tol_per_cds = tol_nu / cds_len if tol_nu else 0.0
            write_lines.append(
               f"{sample_id}\t{pop}\t{vt}\t{de_nu}\t{de_per_cds}\t{args.region_type}\n")

    return write_lines

def main():

    # header = f'Sample_ID\tPopulation\tAllele_Count\tHomozygous_Count\tHeterozygous_Count\n'
    header  = f'ind\tpop\ttype\tde_nu\tde_per_cds\tregion_type\n'
    with open(os.path.join(args.output_dir, args.region_type+'.pop_mutation_burden.txt'), 'w') as concat_file:
        concat_file.write(header)
        for ref_genome in ref_genomes:
            cds_length_dic = {}
            bed_region_per_pop = {}  # {pop: {chrom: [(start1, end1), (start2, end2), ...]}

            cds_bed = cds_masked_bed_dic[ref_genome]

            ind_pop_dic, pop_ind_dic = read_individual_pop_map(
                individuals_map_file_dic[ref_genome])

            if args.region_type != 'all_masked':
                for region_bed in args.bed_files:
                    if not region_bed.endswith('.bed'):
                        continue

                    pop = os.path.basename(region_bed).split('.')[2]

                    if pop not in pop_ind_dic:
                        continue

                    tmp_bedfile = os.path.join(temp_dir, f'{ref_genome}_{pop}.{args.region_type}.type_vs_cds.bed')

                    cds_length_dic[pop] = intersect_and_calc_cds_length(cds_bed, region_bed, tmp_bedfile)
                    
                    bed_region_per_pop[pop] = build_bed_intervals(tmp_bedfile)
                    print(f'Population {pop} has {len(bed_region_per_pop[pop])} chromosomes in bed file on {ref_genome}.\n')
            elif args.region_type == 'all_masked':
                region_bed = cds_bed  # Use cds_bed as the region_bed
                for pop in pop_ind_dic.keys():
                    tmp_bedfile = os.path.join(temp_dir, f'{ref_genome}_{pop}.{args.region_type}.type_vs_cds.bed')
                    cds_length_dic[pop] = intersect_and_calc_cds_length(cds_bed, region_bed, tmp_bedfile)
                    bed_region_per_pop[pop] = build_bed_intervals(tmp_bedfile)
                    print(f'Population {pop} has {len(bed_region_per_pop[pop])} chromosomes in bed file on {ref_genome}.\n')

            with open(os.path.join(args.output_dir, args.region_type+'.'+ref_genome+'_mutation_burden.txt'), 'w') as out_file:
                out_file.write(header)

                print(f'Processing {ref_genome}...\n')
                sift_annotations = read_sift4g_annotation(
                    os.path.join(sift4g_dir_dic[ref_genome], 'merged.SIFTannotations.xls.cleaned'))

                sift4g_annotation_num = sum(
                    len(sift_annotations[chrom]) for chrom in sift_annotations.keys())

                print(f'Read {sift4g_annotation_num} SIFT4G annotations for {ref_genome} completed.\n')

                write_lines = process_vcf_files(
                    ref_genome, sift_annotations, cds_length_dic, bed_region_per_pop)
                out_file.writelines(write_lines)
                concat_file.writelines(write_lines)


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    main()


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
