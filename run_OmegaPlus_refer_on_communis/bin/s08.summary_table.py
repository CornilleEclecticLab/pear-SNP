#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s08.summary_table.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2024/12/29 14:59:31
# @Description:
#    

import datetime
import sys
import textwrap
import os

from s00_config_summary_table import *
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
# os.makedirs(sub_script_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


# Load positive selection genes
pos_genes = {}
with open(positive_selection_output_list) as f:
    files = f.read().strip().split('\n')

common_and_specific_pop_files = {}
for file in files:
    with open(file) as f:
        for l in f:
            if l.startswith('#'):
                continue
            splits = l.strip().split('\t')
            gene = splits[0]
            method = splits[1]
            try:
                pop = splits[2]
            except IndexError:
                pop = os.path.basename(file).split('.')[2]
                if method=='RAiSD' and pop == 'by_win':
                    pop = os.path.basename(file).split('.')[3]
            pos_genes[gene] = pos_genes.get(gene, {})   # pos_genes {gene:{}}
            pos_genes[gene][pop] = pos_genes[gene].get(pop, [])  # pos_genes {gene:{pop:[]}}
            if method in pos_genes[gene][pop]:
                print(f'Warning: {gene} {pop} {method} already in {positive_selection_output_list}')
                #raise ValueError(f'Warning: {gene} {pop} {method} already in {positive_selection_output_list}')
            pos_genes[gene][pop].append(method)                  # pos_genes {gene:
                                                                 #                {pop1:[method1, method2, ...]}
                                                                 #                 pop2:[method1, ...]}
                                                                 # }
for pop in set(pop for pop_methods in pos_genes.values() for pop in pop_methods):
    common_and_specific_pop_files[pop] = open(os.path.join(
        sub_output_dir, f'common_and_specific.positive_selection_genes.{pop}.txt'), 'w')
    # Gene, method, pop
    for gene, pop_methods in pos_genes.items():
        for pop_key in pop_methods.keys():
            if pop_key == pop:
                common_and_specific_pop_files[pop].write(gene+'\t'+pop_key+'\t'+','.join(pop_methods[pop_key])+'\n')
    
# Load chr ID and accession map
chr_ID_accession_map = {}
if os.path.exists(chrID_map):
    with open(chrID_map) as f:
        for l in f:
            splits = l.strip().split()
            chr_ID_accession_map[splits[0]] = splits[1]

# Load gene list from gff file
gff_dic = {}
with open(gff_ori) as f:
    for l in f:
        if l.startswith('#'):
            continue
        splits = l.strip().split('\t')
        chr_accession = splits[0]
        type = splits[2]
        if type != 'gene':
            continue
        start = splits[3]
        end = splits[4]
        attributions = splits[8].replace(" ",'').split(';')
        dic = {i.split('=')[0]:i.split('=')[1] for i in attributions if i != ''}
        if 'Name' in dic:
            name = dic['Name']
        elif 'Accession' in dic:
            name = dic['Accession']
        else:
            raise ValueError(f'No Name or Accession in {gff_ori} line {l}')
        gff_dic[name] = {'chr_accession': chr_accession, 'start':start, 'end':end}

# Load gene collinearity information
collinearity_dic = {}
with open(collinearity) as f:
    for line in f:
        gene1, gene2 = line.strip().split()
        collinearity_dic[gene1] = collinearity_dic.get(gene1, [])
        collinearity_dic[gene1].append(gene2)
        collinearity_dic[gene2] = collinearity_dic.get(gene2, [])
        collinearity_dic[gene2].append(gene1)
        

# Load gene function from eggNOG
eggNOG_dic = {}
with open(eggNOG) as f:
    for l in f:
        if l.startswith('##'):
            continue
        elif l.startswith('#'):
            header = l.strip()[1:].split('\t')
            continue
        
        splits = l.strip().split('\t')
        
        query = splits[header.index('query')]
        if query.startswith('GWH') and "Gene=" in query:
            gene = query.split('Gene=')[1].split('\\t')[0]
        elif query.startswith('pycom'):
            gene = query.split('\\t')[0]
        else:
            raise ValueError(f'Unknown query format in emapper.annotations: {query}')

        name = splits[header.index('Preferred_name')]
        
        pfam = splits[header.index('PFAMs')]
        
        description = splits[header.index('Description')]
        
        kegg_ko = splits[header.index('KEGG_ko')]
        
        kegg_pathways = splits[header.index('KEGG_Pathway')]

        
        go = splits[header.index('GOs')]
        
        eggNOG_dic[gene] = {'name':name,
                            'pfam':pfam,
                            'description':description,
                            'kegg_ko':kegg_ko,
                            'kegg_pathways':kegg_pathways,
                            'go':go}


# Laod blast annotation
blast_annotation_dic = {}
with open(blast_annotation) as f:
    lines = f.readlines() 
    blast_header = '\t'.join(lines[0].strip().split('\t')[1:])
    blast_annotation_dic = {line.strip().split('\t')[0]:'\t'.join(line.strip().split('\t')[1:]) for line in lines[1:]}


# Write summary table
## Expected columns: Gene_ID, Chr, Start, End, Detected_methods, common_vs_specific, Detedted_populations, Gene_name Gene_description
with open(os.path.join(sub_output_dir,'positive_selection_summary_table.tsv'),'w') as f:
    header = '\t'.join(['Gene_ID', 'Chr_accession', 'Chr_ID', 'Start', 'End', 'Detected_methods',
                        'Collinarity_with_the_other_reference_genome', 'Common_vs_specific_in_population',
                        'Detected_populations', 'Gene_name', 'PFAM', 'Gene_description', blast_header,
                        'GO', 'KEGG_Pathway', 'Reference_articles'])+'\n'
    f.write(header)
    
    specific_pop_files = {}
    for pop in set(pop for pop_methods in pos_genes.values() for pop in pop_methods):
        specific_pop_files[pop] = open(os.path.join(sub_output_dir, f'specific.positive_selection_genes_list.{pop}.txt'), 'w')

    specific_summary = open(os.path.join(sub_output_dir, 'specific.positive_selection_genes.summary.txt'), 'w')
    specific_summary.write(header)
    
    gene_list = sorted(pos_genes.keys())
    for gene in gene_list:
        if gene in gff_dic:
            chr_accession = gff_dic[gene]['chr_accession']
            if chr_accession.startswith("Chr") or chr_accession.startswith("chr"):
                chr = chr_accession
            elif chr_accession in chr_ID_accession_map:
                chr = chr_ID_accession_map[chr_accession]
            else:
                raise ValueError(f'Unknown chr accession: {chr_accession}')
            
            start = gff_dic[gene]['start']
            end = gff_dic[gene]['end']
        else:
            raise ValueError(f'Gene {gene} not found in gff file {gff_ori}')
        
        if gene in collinearity_dic:
            collinearity = ','.join(collinearity_dic[gene])
        else:
            collinearity = 'Not_in_collinearity'
        detected_methods = ','.join([','.join(methods)
                                    for methods in pos_genes[gene].values()])
        # remove duplicated methods
        detected_methods = ','.join(list(set(detected_methods.split(','))))

        if len(pos_genes[gene]) > 1:
            common_vs_specific = 'common'
        elif len(pos_genes[gene]) == 1:
            common_vs_specific = 'specific'
        else:
            raise ValueError(f'No population detected for gene {gene}')
        
        detected_populations = ','.join(list(pos_genes[gene].keys()))
        
        gene_name = eggNOG_dic.get(gene,{}).get('name','NA')
        
        pfam = eggNOG_dic.get(gene,{}).get('pfam','NA')
        
        go = eggNOG_dic.get(gene,{}).get('go','NA')
        
        kegg_pathways = eggNOG_dic.get(gene,{}).get('kegg_pathways','NA')
        kegg_pathways = kegg_pathways.split(',')# eg. - | ko03008, map03008
        kegg_pathways = ','.join(
            [kp for kp in kegg_pathways if not kp.startswith('map')])
        
        gene_description = eggNOG_dic.get(gene,{}).get('description','NA')
        
        blast_annotation = blast_annotation_dic.get(gene, str('NA\t'*len(blast_header.split('\t'))).strip())
        
        write_line = f'{gene}\t{chr_accession}\t{chr}\t{start}\t{end}\t{detected_methods}\t{collinearity}\t{common_vs_specific}\t{detected_populations}\t{gene_name}\t{pfam}\t{gene_description}\t{blast_annotation}\t{go}\t{kegg_pathways}\t\n'

        f.write(write_line)

        if common_vs_specific == 'specific':
            specific_pop = list(pos_genes[gene].keys())[0]
            specific_pop_files[specific_pop].write(gene+'\n')
            specific_summary.write(write_line)

# Close files
for pop_file in common_and_specific_pop_files.values():
    pop_file.close()
for pop_file in specific_pop_files.values():
    pop_file.close()
specific_summary.close()


end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))
