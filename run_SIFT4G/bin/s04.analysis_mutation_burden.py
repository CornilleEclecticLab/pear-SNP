#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s04.analysis_mutation_burden.py
# @Author   : NIE Yuqi (based on CHEN Xilong's methodology)
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/06/03 17:17:00
# @Description:
#    

from s04_config import *
import datetime
import sys
import textwrap
import os
import glob
import warnings

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


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Function
def read_sift4g_annotation(file_path):
    # Arg file_path is the path to the SIFT4G annotation file (.xls or .xls.cleaned)
    # Return the following dictionary
    annotations = {} # {chr: {pos:
                     #           {'type': str, 'score': float, 'prediction': str}}}
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
            variant_type = parts[8]  # e.g. SYNONYMOUS, NONSYNONYMOUS, START-LOST, STOP-GAIN, STOP-LOSS
            sift_score = parts[12]
            prediction = parts[16]
            
            if variant_type != 'NONSYNONYMOUS' or sift_score == 'NA':
                continue

            annotations.setdefault(chr, {})
            annotations[chr][pos] = {
                'variant_type' : variant_type,
                'sift_score' : float(sift_score),
                'prediction' : prediction,
                'ref' : ref,
                'alt' : alt
            }
    return annotations


# Funtion to parse genotype in VCF
def parse_vcf_genotype(gt_index, genotype): 
    gt = genotype.split(':')[gt_index]
    gt = gt.replace('/', '|').replace('\\', '|').split('|')

    if len(gt) != 2:
        raise ValueError(f'Invalid genotype format: {genotype}')
    
    sum1_m = sum([1 for g in gt if g == '1'])
    sum0_m = sum([1 for g in gt if g == '0'])
    
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


# Function to parse VCF
def parse_vcf(file_path, sift_annotations):
    # Arg file_path is the path to the VCF file
    # sift_annotations is the dictionary returned by func read_sift4g_annotation
    # return following dictionary
    variants = {}

    file = open(file_path, 'r')
    
    # process each variant
    for line in file:
        line = line.strip()
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            header = line.lstrip('#').strip().split('\t')
            continue
            
        parts = line.split('\t')
        chr, pos, id, ref, alt, qual, filter, info, format, *genotypes = parts
            
        if pos not in sift_annotations[chr].keys():
            continue
            
        prediction = sift_annotations[chr][pos]['prediction']
            
        # Chekc for ref/alt mismatch
        if ref != sift_annotations[chr][pos]['ref'] or alt != sift_annotations[chr][pos]['alt']:
            raise ValueError(
                f'Variant {chr}:{pos} has different ref/alt in VCF and SIFT4G annotation: '
                f'VCF({ref}/{alt}) vs SIFT4G({sift_annotations[chr][pos]["ref"]}/{sift_annotations[chr][pos]["alt"]})')

        gt_index = format.split(':').index('GT')

        # process each genotype
        for i, sample_id in enumerate(header[9:]):
            genotype = genotypes[i]
            gt = parse_vcf_genotype(gt_index, genotype)

            variants.setdefault(sample_id, {'Allele': {'DELETERIOUS': 0, 'TOLERATED': 0},
                                            'Homozygous': {'DELETERIOUS': 0, 'TOLERATED': 0}, 
                                            'Heterozygous': {'DELETERIOUS': 0, 'TOLERATED': 0}})

            if gt == 'Homozygous':
                variants[sample_id]['Allele'][prediction] += 2     # Alt allele counts as 2
                variants[sample_id]['Homozygous'][prediction] += 2
            elif gt == 'Heterozygous':
                variants[sample_id]['Allele'][prediction] += 1
                variants[sample_id]['Heterozygous'][prediction] += 1

    return variants

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
def process_vcf_files(ref_genome, sift_annotations):

    # Return str lines to be written to output file
    write_lines = []
    sample_variants = {}
    
    ind_pop_dic, pop_ind_dic = read_individual_pop_map(
        individuals_map_file_dic[ref_genome])

    chr_list_file = chr_list_file_dic[ref_genome]

    with open(chr_list_file, 'r') as file:
        chr_list = [line.strip() for line in file if line.strip()]

    processed_vcf_dic = {chr: {pop: False for pop in pop_ind_dic.keys()} for chr in chr_list}
    total_vcf_num = len(chr_list) * len(pop_ind_dic.keys())

    vcf_file_dir = vcf_file_dir_dic[ref_genome]
    vcf_file_list = glob.glob(os.path.join(vcf_file_dir, '*.vcf'))
    
    n_vcf = 0
    for vcf_file in vcf_file_list:
        if not vcf_file.endswith('.vcf'):
            continue
        pop_vcf = os.path.basename(vcf_file).split('.')[-2]   # VCF filename format: *.chr.pop.vcf
        chr_vcf = os.path.basename(vcf_file).split('.')[-3]

        if pop_vcf not in pop_ind_dic.keys() or chr_vcf not in chr_list:
            # print(
            #     f'Skipping {vcf_file} due to mismatch in population {pop_vcf} or chromosome {chr_vcf}.\n')
            continue
        
        if not processed_vcf_dic[chr_vcf][pop_vcf]:
            n_vcf += 1
            print(f'Processing {n_vcf}/{total_vcf_num} VCF file: {vcf_file} for  population {pop_vcf} and chromosome {chr_vcf}...\n')
            s_variant= parse_vcf(vcf_file, sift_annotations)
            print(f'Parsed {len(s_variant)} samples from VCF file {vcf_file}.\n')
        else:
            raise ValueError(
                f'VCF file {vcf_file} for chromosome {chr_vcf} and population {pop_vcf} has already been processed.')
        
        for sample_id in s_variant.keys():
            if pop_vcf != ind_pop_dic[sample_id]:
                raise ValueError(
                    f'Sample {sample_id} in VCF file {vcf_file} does not match population {pop_vcf}.')

        for sample_id, variants in s_variant.items():
            if sample_id not in sample_variants:
                sample_variants[sample_id] = variants
                continue
            sample_variants[sample_id]['Allele']['DELETERIOUS'] += variants['Allele']['DELETERIOUS']
            sample_variants[sample_id]['Homozygous']['DELETERIOUS'] += variants['Homozygous']['DELETERIOUS']
            sample_variants[sample_id]['Heterozygous']['DELETERIOUS'] += variants['Heterozygous']['DELETERIOUS']
            sample_variants[sample_id]['Allele']['TOLERATED'] += variants['Allele']['TOLERATED']
            sample_variants[sample_id]['Homozygous']['TOLERATED'] += variants['Homozygous']['TOLERATED']
            sample_variants[sample_id]['Heterozygous']['TOLERATED'] += variants['Heterozygous']['TOLERATED']

        processed_vcf_dic[chr_vcf][pop_vcf] = True

    for chr, pop_dict in processed_vcf_dic.items():
        for pop, processed in pop_dict.items():
            if not processed:
                raise ValueError(f'No VCF file found for {chr} in population {pop}.')
    
    print(len(sample_variants), 'samples found in all VCF files for all populations and chromosomes.\n')
    for sample_id, variants in sample_variants.items():
        pop = ind_pop_dic[sample_id]
        for var_type, vt in (['Allele','All'],  ['Homozygous', 'Ho'], ['Heterozygous','He']):
           write_lines.append(
               f"{sample_id}\t{pop}\t{vt}\t{variants[var_type]['DELETERIOUS']}\t{variants[var_type]['TOLERATED']}\n")

    return write_lines

def main():
    # header = f'Sample_ID\tPopulation\tAllele_Count\tHomozygous_Count\tHeterozygous_Count\n'
    header  = f'ind\tpop\ttype\tde_nu\ttol_nu\n'
    with open(os.path.join(sub_output_dir, 'pop_mutation_burden.txt'), 'w') as concat_file:
        concat_file.write(header)
        for ref_genome in ref_genomes:
            with open(os.path.join(sub_output_dir, ref_genome+'_mutation_burden.txt'), 'w') as out_file:
                out_file.write(header)
                
                print(f'Processing {ref_genome}...\n')
                sift_annotations = read_sift4g_annotation(
                    os.path.join(sift4g_dir_dic[ref_genome], 'merged.SIFTannotations.xls.cleaned'))        

                sift4g_annotation_num = sum(
                    len(sift_annotations[chr]) for chr in sift_annotations.keys())

                print(f'Read {sift4g_annotation_num} SIFT4G annotations for {ref_genome} completed.\n')

                write_lines = process_vcf_files(ref_genome, sift_annotations)
                out_file.writelines(write_lines)
                concat_file.writelines(write_lines)
        
if __name__ == "__main__":
    main()


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))