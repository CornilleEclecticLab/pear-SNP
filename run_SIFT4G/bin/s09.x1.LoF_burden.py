#!/usr/bin/env python3
# _*_ coding: utf-8 _*_


from s04_config import *
from s06_config import *
import datetime
import sys
import os
import glob
import argparse
import subprocess
import gzip

try:
    from intervaltree import IntervalTree
    USE_INTERVALTREE = True
except ImportError:
    USE_INTERVALTREE = False
    print("TIP: Active or install 'intervaltree' to speed up interval operations.")

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')


version = "1.0.0"
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
work_dir = os.path.dirname(script_path)
output_dir = os.path.join(work_dir, 'output')
temp_dir = os.path.join(work_dir, 'temp')
os.makedirs(temp_dir, exist_ok=True)

lof_dir_dic = {
    ref: os.path.join(output_dir, 's08.calculate_LoF', ref) for ref in ref_genomes
}

LOF = 'LOF'
SYN = 'SYNONYMOUS'

# Canonical splice donor/acceptor sites lie 1-2 bp outside the exon, i.e. outside
# cds_masked.bed; pad each CDS interval so they are counted.
SPLICE_PAD = 2

# High-confidence LoF classes (LOFTEE-like core classes; SNPs only, so no frameshift)
HC_TERMS = {'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'}


def parse_args():
    parser = argparse.ArgumentParser(
        description="LoF burden per individual (denominator: synonymous SNP)")
    parser.add_argument(
        "output_dir", help="Output directory (e.g. ../output/s09.LoF_burden/consensus)")
    parser.add_argument(
        "region_type", help="Region type (selection, control, all_masked)")
    parser.add_argument(
        "lof_set", choices=['consensus', 'union', 'hc'],
        help="consensus = VEP & snpEff, union = VEP | snpEff, hc = high-confidence (see header)")
    parser.add_argument("bed_files", nargs="*",
                        help="One or more BED files named by *.*.pop.*.bed")
    return parser.parse_args()


def intersect_and_calc_cds_length(cds_bed, region_bed, tmp_bedfile):
    load_bedtools = 'source s00.load_bedtools.sh'
    if args.region_type == 'control':
        cmd_int = f"bedtools intersect -nonamecheck -v -a {cds_bed} -b {region_bed}"
    else:
        cmd_int = f"bedtools intersect -nonamecheck -wa -a {cds_bed} -b {region_bed}"
    with open(tmp_bedfile, "w") as outf:
        subprocess.run(f"{load_bedtools} && {cmd_int}", stdout=outf,
                       shell=True, check=True, executable="/bin/bash")
    total_len = 0
    with open(tmp_bedfile) as fin:
        for line in fin:
            cols = line.split("\t")
            if len(cols) >= 3:
                total_len += int(cols[2]) - int(cols[1])
    return total_len


def pad_cds_bed(cds_bed, padded_bed):
    cmd = (f"source s00.load_bedtools.sh && "
           f"awk -v p={SPLICE_PAD} 'BEGIN{{OFS=\"\\t\"}} {{s=$2-p; if (s<0) s=0; print $1, s, $3+p}}' {cds_bed} "
           f"| bedtools sort | bedtools merge > {padded_bed}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    return padded_bed


def build_bed_intervals(tmp_bedfile):
    bed_dic = {}
    tree_dic = {}
    with open(tmp_bedfile) as fin:
        for line in fin:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            chrom, s, e = cols[0], int(cols[1]), int(cols[2])
            bed_dic.setdefault(chrom, []).append((s, e))
            if USE_INTERVALTREE:
                tree_dic.setdefault(chrom, IntervalTree())
                tree_dic[chrom].addi(s, e)
    return tree_dic if USE_INTERVALTREE else bed_dic


def position_in_intervals(pos, intervals):
    # pos: 1-based VCF; BED: [start, end) 0-based
    zpos = int(pos) - 1
    if USE_INTERVALTREE:
        return bool(intervals[zpos])
    for (s, e) in intervals:
        if s <= zpos < e:
            return True
    return False


def load_site_classes(ref_genome, lof_set):
    # Return {chrom: 
    #          {pos: (ref, alt, 'LOF' | 'SYNONYMOUS')
    #           }
    #         }
    # LoF wins over synonymous
    sites = {}

    syn_file = os.path.join(lof_dir_dic[ref_genome], 'merged.synonymous_consensus.tsv')
    with open(syn_file) as fin:
        for line in fin:
            chrom, pos, ref, alt = line.rstrip('\n').split('\t')[0:4]
            sites.setdefault(chrom, {})[pos] = (ref, alt, SYN)

    lof_file = os.path.join(lof_dir_dic[ref_genome], 'merged.LoF_merged.tsv')
    with open(lof_file) as fin:
        header = fin.readline().rstrip('\n').split('\t')
        i_vep, i_snp = header.index('IS_VEP'), header.index('IS_SNPEFF')
        i_tag, i_csq = header.index('SNPEFF_LOF_TAG'), header.index('VEP_CONSEQUENCE')
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            is_vep, is_snp = parts[i_vep] == '1', parts[i_snp] == '1'
            if lof_set == 'consensus':
                keep = is_vep and is_snp
            elif lof_set == 'union':
                keep = is_vep or is_snp
            else:  # hc
                keep = (is_vep and is_snp and parts[i_tag] != '.'
                        and any(t in HC_TERMS for t in parts[i_csq].replace(',', '&').split('&')))
            if keep:
                sites.setdefault(parts[0], {})[parts[1]] = (parts[2], parts[3], LOF)
    return sites


def parse_genotype(gt_index, genotype):
    gt = genotype.split(':')[gt_index]
    alleles = gt.replace('/', '|').replace('\\', '|').split('|')
    if len(alleles) != 2:
        raise ValueError(f'Invalid genotype format: {genotype}')
    sum1 = sum(1 for a in alleles if a == '1')
    sum0 = sum(1 for a in alleles if a == '0')
    if '.' in alleles and sum1 + sum0 == 0:
        return None
    if '.' in alleles and sum1 + sum0 == 1:
        raise ValueError(f'Invalid genotype format with missing data: {genotype}')
    if sum1 == 2:
        return 'Homozygous'
    if sum1 == 1:
        return 'Heterozygous'
    if sum0 == 2:
        return 'Non-variant'
    raise ValueError(f'Invalid genotype format: {genotype}')


def new_stats():
    return {k: {LOF: 0, SYN: 0} for k in ('Total_site', 'Allele', 'Homozygous', 'Heterozygous')}


def parse_vcf(file_path, site_classes, bed_region):
    variants = {}
    with gzip.open(file_path, 'rt') as fin:
        for line in fin:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                sample_ids = line.rstrip('\n').lstrip('#').split('\t')[9:]
                continue

            parts = line.rstrip('\n').split('\t')
            chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]

            if chrom not in site_classes or pos not in site_classes[chrom]:
                continue
            if chrom not in bed_region or not position_in_intervals(pos, bed_region[chrom]):
                continue

            ann_ref, ann_alt, cls = site_classes[chrom][pos]
            if ref != ann_ref or alt != ann_alt:
                raise ValueError(
                    f'Variant {chrom}:{pos} in VCF({ref}/{alt}) and LoF annotation ({ann_ref}/{ann_alt}) do not match.')

            gt_index = parts[8].split(':').index('GT')
            for idx, sample_id in enumerate(sample_ids):
                gt = parse_genotype(gt_index, parts[9 + idx])
                if gt is None or gt == 'Non-variant':
                    continue
                st = variants.setdefault(sample_id, new_stats())
                st['Total_site'][cls] += 1
                if gt == 'Homozygous':
                    st['Allele'][cls] += 2
                    st['Homozygous'][cls] += 2
                else:
                    st['Allele'][cls] += 1
                    st['Heterozygous'][cls] += 1
    return variants


def read_individual_pop_map(file_path):
    individual_pop_map = {}
    pop_individual = {}
    with open(file_path) as fin:
        for line in fin:
            ind, pop = line.strip().split()[0:2]
            individual_pop_map[ind] = pop
            pop_individual.setdefault(pop, []).append(ind)
    return individual_pop_map, pop_individual


def process_vcf_files(ref_genome, site_classes, cds_length_dic, bed_region_per_pop):
    write_lines = []
    write_lines_summary = []
    sample_stats = {}

    ind_pop_dic, pop_ind_dic = read_individual_pop_map(individuals_map_file_dic[ref_genome])
    with open(chr_list_file_dic[ref_genome]) as fin:
        chr_list = [line.strip() for line in fin if line.strip()]

    processed = {c: {p: False for p in pop_ind_dic} for c in chr_list}
    total_vcf_num = len(chr_list) * len(pop_ind_dic)
    n_vcf = 0

    for vcf_file in glob.glob(os.path.join(vcf_file_dir_dic[ref_genome], '*.vcf.gz')):
        basename = os.path.basename(vcf_file)
        # VCF filename format: *.chr.pop.vcf.gz
        pop_vcf = basename.split('.')[-3]
        chr_vcf = basename.split('.')[-4]
        if pop_vcf not in pop_ind_dic or chr_vcf not in chr_list:
            continue
        if processed[chr_vcf][pop_vcf]:
            raise ValueError(f'VCF file {vcf_file} for {chr_vcf} / {pop_vcf} has already been processed.')

        n_vcf += 1
        print(f'Processing {n_vcf}/{total_vcf_num} VCF file: {vcf_file}')

        if pop_vcf not in bed_region_per_pop or chr_vcf not in bed_region_per_pop[pop_vcf]:
            processed[chr_vcf][pop_vcf] = True
            continue

        for sample_id, st in parse_vcf(vcf_file, site_classes, bed_region_per_pop[pop_vcf]).items():
            if pop_vcf != ind_pop_dic[sample_id]:
                raise ValueError(f'Sample {sample_id} in {vcf_file} does not match population {pop_vcf}.')
            if sample_id not in sample_stats:
                sample_stats[sample_id] = st
                continue
            for k in st:
                for cls in (LOF, SYN):
                    sample_stats[sample_id][k][cls] += st[k][cls]

        processed[chr_vcf][pop_vcf] = True

    for chrom, pop_dict in processed.items():
        for pop, done in pop_dict.items():
            if not done:
                raise ValueError(f'No VCF file found for {chrom} in population {pop}.')

    print(f'{len(sample_stats)} samples found for {ref_genome}.')
    for sample_id, st in sample_stats.items():
        pop = ind_pop_dic[sample_id]
        cds_len = cds_length_dic[pop]
        if cds_len == 0:
            raise ValueError(f'Population {pop} CDS length is 0, please check the BED files.')

        denominator = st['Total_site'][SYN]
        for var_type, vt in (('Allele', 'All'), ('Homozygous', 'Ho'), ('Heterozygous', 'He')):
            lof_nu = st[var_type][LOF]
            lof_per_syn = lof_nu / denominator if lof_nu and denominator else 0.0
            write_lines.append(
                f"{sample_id}\t{pop}\t{vt}\t{lof_nu}\t{lof_per_syn}\t{args.region_type}\t{denominator}\n")

        lof_site = st['Total_site'][LOF]
        if lof_site != st['Homozygous'][LOF] // 2 + st['Heterozygous'][LOF]:
            raise ValueError(f"DEVELOPMENTAL BUG: Sample {sample_id} {pop} LoF count mismatch")

        write_lines_summary.append(
            f"{sample_id}\t{pop}\t{args.region_type}\t{cds_len}\t"
            f"{st['Allele'][LOF]}\t{st['Homozygous'][LOF]}\t{st['Heterozygous'][LOF]}\t{lof_site}\t"
            f"{st['Allele'][SYN]}\t{st['Total_site'][SYN]}\n")

    return write_lines, write_lines_summary


def main():
    header = 'Sample_ID\tPopulation\tAllelic_state\tLoF_allele_count\tLoF_per_synonymous_site\tRegion_type\tSynonymous_site_count\n'
    header_summary = 'Sample_ID\tPopulation\tRegion_type\tCDS_length\tLoF_allele_count\tLoF_homozygous_count\tLoF_heter_count\tLoF_site_count\tSynonymous_allele_count\tSynonymous_site_count\n'

    out_file = os.path.join(args.output_dir, f'{args.region_type}.synonymous.pop_mutation_burden.txt')
    out_summary = os.path.join(args.output_dir, f'{args.region_type}.synonymous.pop_mutation_burden.summary.txt')

    with open(out_file, 'w') as fout, open(out_summary, 'w') as fsum:
        fout.write(header)
        fsum.write(header_summary)

        for ref_genome in ref_genomes:
            cds_length_dic = {}
            bed_region_per_pop = {}
            cds_bed = pad_cds_bed(
                cds_masked_bed_dic[ref_genome],
                os.path.join(temp_dir, f'{ref_genome}.s09_LoF_{args.lof_set}.{args.region_type}.cds_masked.pad{SPLICE_PAD}.bed'))
            ind_pop_dic, pop_ind_dic = read_individual_pop_map(individuals_map_file_dic[ref_genome])

            if args.region_type in ('selection', 'control'):
                region_beds = [(os.path.basename(b).split('.')[2], b)
                               for b in args.bed_files if b.endswith('.bed')]
            elif args.region_type == 'all_masked':
                region_beds = [(pop, cds_bed) for pop in pop_ind_dic]
            else:
                raise ValueError(f'Invalid region type: {args.region_type}.')

            for pop, region_bed in region_beds:
                if pop not in pop_ind_dic:
                    continue
                tmp_bedfile = os.path.join(
                    temp_dir, f'{ref_genome}_{pop}.s09_LoF_{args.lof_set}.{args.region_type}.type_vs_cds.bed')
                cds_length_dic[pop] = intersect_and_calc_cds_length(cds_bed, region_bed, tmp_bedfile)
                bed_region_per_pop[pop] = build_bed_intervals(tmp_bedfile)

            print(f'Processing {ref_genome}...')
            site_classes = load_site_classes(ref_genome, args.lof_set)
            n_lof = sum(1 for c in site_classes.values() for v in c.values() if v[2] == LOF)
            n_syn = sum(1 for c in site_classes.values() for v in c.values() if v[2] == SYN)
            print(f'{ref_genome}: {n_lof} LoF ({args.lof_set}) and {n_syn} synonymous sites loaded.')

            write_lines, write_lines_summary = process_vcf_files(
                ref_genome, site_classes, cds_length_dic, bed_region_per_pop)
            fout.writelines(write_lines)
            fsum.writelines(write_lines_summary)


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    main()


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time - start_time).center(79))
