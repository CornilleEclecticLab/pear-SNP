#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : cal_TE_load.py
# @Author   : NIE Yuqi
# @Time(Europe/Paris): 2026/06/30 13:31:01
# @Description:
#

import datetime
import gzip
import os
import sys
import textwrap
from collections import Counter, defaultdict


import s00_config as cfg 



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))


# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Function get sample_population dic
def get_sample_pop_dic(sample_pop_file):
    with open(sample_pop_file, 'r') as f:
        sample_pop_dic = {}
        for line in f:
            sample_id, population = line.strip().split()
            sample_pop_dic[sample_id] = population
    return sample_pop_dic


# Function to initial dict for counting 
def new_node():
    # Initial node with 0 count
    return {'All': 0, 'Ho': 0, 'He': 0, 'n_called_alleles': 0, 'n_sites': 0}
def add_node(dst, src):
    for k in ('All', 'Ho', 'He', 'n_called_alleles', 'n_sites'):
        dst[k] += src[k]
 

# Function to load and retrieve the classification from REPET classif file
def load_classification(classification_files):
    dic = {}
    for file in classification_files:
        with open(file, 'r') as f:
            header = False
            for line in f:
                if header is False:
                    if line.strip().startswith('Seq_name'):
                        header = line.strip().split()
                        continue
                    else:
                        print(line.strip())
                        raise ValueError(f"Header not found in {file}")
                
                #Seq_name        length  strand  confused        class   order   Wcode   sFamily CI      coding  struct  other
                line_list = line.strip().split()
                for n, col in enumerate(header):
                    if col == 'Seq_name':
                        seq_name = line_list[n]
                        dic[seq_name] = dic.get(seq_name, {})
                    
                    if col in ['class', 'order', 'Wcode', 'sFamily', 'CI']:
                        col_val = line_list[n].strip().split('|')[0]
                    else:
                        col_val = line_list[n]
                    dic[seq_name][col] = col_val
    return dic


# Functions to parse vcf
# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
class vcf():
    def __init__(self, sample_pop_dic, classification_dic):
        self.sample_pop_dic = sample_pop_dic
        self.classification_dic = classification_dic

    # Function get classification
    def get_classification(self, seq_list):
        wcode_ls = [self.classification_dic.get(seq, {}).get('Wcode') for seq in seq_list]
        wcode_counts = Counter(wcode_ls)    
        top_wcode, top_count = wcode_counts.most_common(1)[0]

        for seq in seq_list:
            if self.classification_dic.get(seq, {}).get('Wcode') == top_wcode:
                if_confused_classification = self.classification_dic.get(seq, {}).get('confused')
                _class = self.classification_dic.get(seq, {}).get('class')
                _order = self.classification_dic.get(seq, {}).get('order')
                Wcode = self.classification_dic.get(seq, {}).get('Wcode')
                _sFamily = self.classification_dic.get(seq, {}).get('sFamily')

                if top_count > len(wcode_ls) / 2:
                    break
                else:
                    confused_set = set([ self.classification_dic.get(seq, {}).get('confused') for seq in seq_list ])
                    if len(confused_set) == 1:
                        if_confused_classification = confused_set.pop()
                    else:
                        if_confused_classification = "True"
                    _class += '_ambigous_classif'
                    _order += '_ambigous_classif'
                    Wcode = '_ambigous_classif'
                    _sFamily = '_ambigous_classif'
                    break

        return if_confused_classification, _class, _order, Wcode, _sFamily


    # Function to parse vcf INFO
    def parse_info(self, info_field: str) -> dict:
        # SVTYPE=DNA;MEI=PCOM_TEdenovoGr-B-G7562-Map8_reversed;MEPRED=PASS;
        if not info_field:
            return {}
        
        return { 
            parts[0]:parts[1] if len(parts) == 2 else True
            for term in info_field.strip().split(';')
            if term and (parts := term.split('=', 1))
        }

    def parse_gt(self, gt_field, pol):
        """Allele-level counting;
          present_alleles  : observed TE-present alleles (ins:'1' count as 1, :'0' count as 1)
          n_called_alleles : observed (non-missing) alleles, 0/1/2
          state            : 'Homozygous'|'Heterozygous'|'No_TE'|None
          All = Ho + He
        """
        # 0/0, 0/1, 0/. 1/1, 1/., 1/. 
        if '/' in gt_field:
            alleles = gt_field.split('/')
        elif '|' in gt_field:
            alleles = gt_field.split('|')
        else:
            raise ValueError(f"GT field {gt_field} is not in the expected format")
        
        if len(alleles) != 2:
            raise ValueError(f"GT field {gt_field} is not in the expected format, expected 2 alleles, got {len(alleles)}")

        sum1 = sum(1 for a in alleles if a=='1')
        sum0 = sum(1 for a in alleles if a=='0')
        n_called = sum0 + sum1
        if n_called != sum(1 for a in alleles if a != '.'):
            raise ValueError(f"DEVELOPMENT Error: Invalid GT: {gt_field}")

        if n_called == 0:
            return 0, 0, 0, None

        n_present = sum1 if pol == 'ins' else sum0

        if n_called == 2:
            state = ('Ho' if n_present == 2      # Homozygous
                    else 'He' if n_present == 1  # Heterozygous
                    else 'No_TE')
            n_site = 1
        elif n_called == 1:
            state = ('He' if n_present == 1
                    else 'No_TE')
            n_site = 1
        else: 
            raise ValueError(f"DEVELOPMENT Error: Biallelic GT? : {gt_field}")
            
        return n_present, n_called, n_site, state


# Filter vcf on FILTER field
    def filter_vcf(self, filter_field):
        ##FILTER=<ID=LC,Description="Low confidence">                                                                                                      >
        ##FILTER=<ID=NU,Description="Not unique">                                                                                                          >
        ##FILTER=<ID=S,Description="Spanning read num is outlier">                                                                                         >
        ##FILTER=<ID=F,Description="Shorter than 50-bp">                                                                                                   >
        ##FILTER=<ID=D,Description="Relative depth of breakpoint is outlier">                                                                              >
        ##FILTER=<ID=G,Description="Outliers during genotyping">                                                                                           >
        ##FILTER=<ID=R,Description="No discordant read stat available">                                                                                    >
        ##FILTER=<ID=Y,Description="Variants on chrY.
        
        filter_list = filter_field.split(';')
        if filter_list[0] == '.':
            raise ValueError(f"Filter field in the VCF file is '.', unexpected in this ")
        if 'PASS' in filter_list:
            return True
        elif 'LC' in filter_list or \
            'NU' in filter_list or \
            'S' in filter_list or \
            'F' in filter_list or \
            'D' in filter_list or \
            'G' in filter_list or \
            'R' in filter_list or \
            'Y' in filter_list:
            return False
        elif 'M' in filter_list:
            return True
        else:
            return False

    # Function distingwish ins and  (polarization)
    def polar_te(self, id_field):
        if '_ins_' in id_field:
            return 'ins'
        elif '_abs_' in id_field:
            return 'abs'
        else:
            raise ValueError(f"TE type in {id_field} is not in the expected format")

    # Function to parse file name as sample_ID
    def get_sample_ID(self, filename):
        if len(filename.split('.')) == 1:
            sample_id = filename
        else:
            sample_id = filename.split('.')[2]
        if sample_id in self.sample_pop_dic:
            return sample_id
        else:
            raise ValueError(f"Sample ID {sample_id} in {filename} not found in sample population dictionary")

    # Function to parse vcf
    def parse_vcf(self, vcf_file, vcf_dic):
        
        if vcf_file.endswith('vcf.gz'):
            fh = gzip.open(vcf_file, 'rt')
        elif vcf_file.endswith('vcf'):
            fh = open(vcf_file, 'r')
        else:
            raise ValueError(f"Not a .vcf or .vcf.gz file: {vcf_file}")
        
        for line in fh:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
                header = line[1:].strip().split('\t')
                samples = header[9:]
                sample_ids = [self.get_sample_ID(sample) for sample in samples]
                sample_vs_id_dic = {sample:sample_ids[i] for i, sample in enumerate(samples)}
                continue
            lines = line.strip().split('\t')
            _chrom = lines[header.index('CHROM')]
            _pos = lines[header.index('POS')]
            _id = lines[header.index('ID')]
            _ref = lines[header.index('REF')]
            _alt = lines[header.index('ALT')]
            _qual = lines[header.index('QUAL')]
            _filter = lines[header.index('FILTER')]
            _info = lines[header.index('INFO')]
            _format = lines[header.index('FORMAT')]

            if not self.filter_vcf(_filter):
                continue

            if ',' in _alt:
                print('>>> DEVELOPMENT Error: Multiple ALT:')
                print(lines)
                
            pol = self.polar_te(_id)  # abs / ins
            te_id = self.parse_info(_info).get('MEI', None)

            if te_id is not None:
                te_id_ls = te_id.split('|')
                if_confused_classification, _class, _order, Wcode, _sFamily = self.get_classification(te_id_ls)
                classification = f"{_class}/{_order}/{_sFamily}|{Wcode}"
            else:
                classification = 'Unknown'
            
            gt_index = _format.split(':').index('GT')

            for sample, genotype in zip(samples, lines[9:]):
                sample_id = sample_vs_id_dic[sample]
                n_present, n_called, n_site, allelic_state = self.parse_gt(genotype.split(':')[gt_index], pol)
                
                if n_called == 0:
                    continue
                node = vcf_dic[sample_id][pol][classification]
                node['n_called_alleles'] += n_called
                node['n_sites'] += n_site
                node['All'] += n_present

                if allelic_state == 'Ho':
                    node['Ho'] += 2
                elif allelic_state == 'He':
                    node['He'] += 1

        fh.close()
        return vcf_dic

# main
def main():
    sample_pop_dic = get_sample_pop_dic(cfg.sample_pop_file)
    output_path = os.path.join(cfg.output_dir, 'TE_load.per_individual.txt')

    vcf_dic = defaultdict(lambda: defaultdict(lambda: defaultdict(new_node)))
    classfi_dic = load_classification(cfg.classif_files)

    for vcf_file in cfg.vcfs:
        print(f'Parsing {vcf_file}')
        vcf(sample_pop_dic, classfi_dic).parse_vcf(vcf_file, vcf_dic)

    with open(output_path,'w') as f:
        f.write("Sample_ID\tPopulation\tPolarization\tTE_classification\tAllelic_state\tPresent_TE_allele_count\tNumber_of_called_alleles\tNumber_of_sites\n")

        def write_row(s, pol, classification, node):
            pop = sample_pop_dic[s]
            for t in ('All', 'Ho', 'He'):
                f.write(f"{s}\t{pop}\t{pol}\t{classification}\t{t}\t{node[t]}\t{node['n_called_alleles']}\t{node['n_sites']}\n")

        for s in sorted(vcf_dic):
            grand = new_node()
            for pol in sorted(vcf_dic[s]):
                pol_total = new_node()
                for classification, node in sorted(vcf_dic[s][pol].items()):
                    write_row(s, pol, classification, node)    # If no need to wirte by each TE_classification, comment this line 
                    add_node(pol_total, node)
                write_row(s, pol, 'All_classification', pol_total)
                add_node(grand, pol_total)
            write_row(s, 'All_present', 'All_classification', grand)


if __name__ == '__main__':
    start_time = datetime.datetime.now()
    print(f'{" Start ":=^79}')

    main()

    end_time = datetime.datetime.now()
    print('')
    print(' END '.center(79,'='))
    print(str(end_time-start_time).center(79))