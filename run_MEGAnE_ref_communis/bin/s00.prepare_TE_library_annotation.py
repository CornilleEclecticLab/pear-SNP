#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s00.prepare_TE_library_annotation.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/06/18 12:32:38
# @Description:
#    

import datetime
import sys
import textwrap
import os
import re
start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')



# Define the input files
consensus_file = os.path.join(
    input_dir, 'Pyrus_cat_refTE_denovoLibTEs_filtered_MCL.fa')

classif_file = os.path.join(
    input_dir, 'Pyrus_cat_refTE_denovoLibTEs_filtered.classif')

tesorter_file = os.path.join(
    work_dir,
    '../run_TEsorter/output/Pyrus_cat_refTE_denovoLibTEs_filtered_MCL.fa.rexdb-plant.cls.tsv')


# Define the output files
output_file = os.path.join(
    input_dir,
    'Pyrus_cat_refTE_denovoLibTEs_filtered_MCL.annotation.fa')
output_file_RM = os.path.join(
    input_dir,
    'Pyrus_cat_refTE_denovoLibTEs_filtered_MCL.RepeatMasker.fa')


# Funtion to read TEsorter classification
def load_te_sorter_classification(te_consensus_sorter):
    header = None
    dic = {}
    with open(te_consensus_sorter, 'r') as f:
        for line in f:
            splits = line.strip().split('\t')
            if line.startswith('#') and header is None:
                header = splits
                continue
            if header is None:
                raise ValueError(f'Header not found in {te_consensus_sorter}')
            seq_id = '_'.join(splits[header.index('#TE')].split('_')[1:])
            
            superfamily = splits[header.index('Superfamily')]
            
            # Merge superfamily names with PASTEC
            if superfamily == 'PIF_Harbinger':
                superfamily = superfamily.replace('_', '-')

            elif superfamily == 'EnSpm_CACTA':
                superfamily = superfamily.replace('EnSpm_', '')
           
            if seq_id not in dic:
                dic[seq_id] = {
                    'order': splits[header.index('Order')],
                    'superfamily': superfamily,
                    'clade': splits[header.index('Clade')]
                }
            else:
                raise ValueError(
                    f'Duplicate seq_id {seq_id} in {te_consensus_sorter} line {line}')
    return dic


# Funtion to parse coding column in PASTEC output
def parse_coding_column(coding_str):
    """
    Parse the coding colum from PASTEC output, e.g.:
    
    coding=(TE_BLRtx: Gypsy-8_PX-I:ClassI:LTR:Gypsy: 14.65%; TE_BLRx: Gypsy-8_PX-I_1p:ClassI:LTR:Gypsy: 27.74%)

    Return the best matched TE category based on the highest percentage.
    """
    if not coding_str.startswith("coding="):
        raise ValueError("The coding string must start with 'coding='.")

    # Get the content inside the parentheses ()
    content = coding_str[len("coding=("):-1]

    if content == 'NA':
        return 'NA'

    # Split by semicolon; or comma,
    entries = [e.strip() for e in re.split(r'[;,]', content) if e.strip()]

    best_cat = None
    best_pct = -1.0

    for entry in entries:
        # Split by colon:
        parts = [p.strip() for p in entry.split(':')]
        if len(parts) < 3:
            continue

        # Get percentage value
        pct_str = parts[-1].rstrip('%').strip()
        try:
            pct = float(pct_str)
        except ValueError:
            continue

        # index of "Class"
        class_idx = next((i for i, p in enumerate(
            parts) if p.startswith('Class')), None)
        if class_idx is None or class_idx >= len(parts) - 1:
            continue

        # Class:order:superfamily
        category = ':'.join(parts[class_idx:-1])

        # fix best category
        if pct > best_pct:
            best_pct = pct
            best_cat = category

    if best_cat is None and 'profiles' in content:
        return 'NA'
    if best_cat is None:
        raise ValueError("No valid TE categories found in the coding string.")
    return best_cat


# Funtion to read TE families
def load_te_families(te_consensus_classif_file, tesort_db):
    te_families = {}
    # Read TE families from the consensus classification file
    with open(te_consensus_classif_file, 'r') as f:
        header = None
        for line in f:
            splits = line.strip().split('\t')
            if header is None:
                header = splits
                continue
            seq_name = splits[header.index('Seq_name')]
            length = splits[header.index('length')]
            confused = splits[header.index('confused')]
            if confused == 'False' or confused == 'True':
                pass
            else:
                raise ValueError(f'Unknown confused value {confused} in {te_consensus_classif_file} line {line}')

            class_ = splits[header.index('class')] # can be 'I', 'II', 'Unclassified', or 'NA', or 'I|II' etc.

            # Consensus class in 'class' values I or II if consensus classified as TE, value NA if consensus classified as not TE,
            # if class_ == 'NA':
            #     continue

            order = splits[header.index('order')]
            superfamily = splits[header.index('sFamily')]
            coding = splits[header.index('coding')]
            seq_name = splits[header.index('Seq_name')]
            
            # WARNING: if consensus is confused, the 'class', 'order', 'Wcode', 'sFamily' and 'CI' fields will contain all information separated by |.

            if confused == "True":
                ci_ls = [ int(i) for i in  splits[header.index('CI')].split('|') ]
                max_ci = max(ci_ls)
                for i, ci in enumerate(ci_ls):
                    if ci == max_ci:
                        class_ = splits[header.index('class')].split('|')[i]
                        order = splits[header.index('order')].split('|')[i]
                        superfamily = splits[header.index('sFamily')].split('|')[i]
                        break

            coding_category = 'NA'
            # if class_ == 'I' or class_ == 'II' or class_ == 'Unclassified' or class_ == 'NA':
            coding_category = parse_coding_column(coding)

            # Check class with coding category
            if coding_category != 'NA':
                if not coding_category.startswith('Class'):
                    raise ValueError(
                        f'Coding category {coding_category} does not start with Class in {te_consensus_classif_file} line {line}')

                class_coding = coding_category.split(':')[0].lstrip('Class')

                if class_ != 'Unclassified' or class_ != 'NA':
                    if class_coding != class_:
                        raise ValueError(
                            f'Class {class_} does not match coding category {coding_category} in {te_consensus_classif_file} line {line}')

                if class_ == 'Unclassified' and class_coding in ('I', 'II'):
                    class_ = class_coding

            # Compare order
            if coding_category != 'NA':
                order_coding = coding_category.split(':')[1]
                sfamily_coding = coding_category.split(':')[2]
                if order_coding == '?':
                    pass
                elif order_coding != '?' and order == 'Unclassified':
                    # print(f'XXXXXXXXXXX{seq_name} Order {order} != {order_coding}, {line}')
                    order = order_coding                                # Assign order from coding category
                elif order_coding != '?' and order != 'Unclassified':
                    if order_coding!= order:
                        print(f'XXXXXXXXXXX Order != coding column: {order} <= {order_coding}, {line}')
                        order = order_coding                            # Assign order from coding category

                if sfamily_coding == '?':
                    pass
                elif sfamily_coding != '?' and superfamily == 'NA':
                    if order_coding == order:                           # Only assign superfamily if order matches
                        superfamily = sfamily_coding                    # Assign superfamily from coding category

            # Improve class, order and superfamily classification using TEsort
            classI_orders = ('LTR', 'SINE', 'LINE', 'PLE', 'DIRS')
            classII_orders = ('TIR', 'Crypton', 'Helitron', 'Maverick')

            if order == 'Unclassified'or order == 'NA':
                if tesort_db.get(seq_name, None) is not None:
                    
                    order_tesort = tesort_db[seq_name]['order']
                    if (class_ =='I' and order_tesort in classI_orders) or \
                        (class_ =='II' and order_tesort in classII_orders):
                        print(">>> Order", {order}, "<=", {order_tesort}, seq_name, "is classified by TEsorter as", tesort_db[seq_name])
                        order = order_tesort
                    elif class_ == 'Unclassified':
                        if order_tesort in classI_orders:
                            print(">>> Order", {order}, "<=", {order_tesort}, seq_name, "is classified by TEsorter as", tesort_db[seq_name])
                            order = order_tesort
                            class_ = 'I'  # Assign class I if order is in class I orders
                        elif order_tesort in classII_orders:
                            print(">>> Order", {order}, "<=", {order_tesort}, seq_name, "is classified by TEsorter as", tesort_db[seq_name])
                            order = order_tesort
                            class_ = 'II'

                        # raise ValueError(
                        #     f'Class {class_} and order {order_tesort} for {seq_name} may mismatch in {te_consensus_classif_file} and te_consensus_sorter')

            if superfamily == 'NA' or superfamily == 'Unclassified':
                if tesort_db.get(seq_name, None) is not None:
                    superfamily_tesort = tesort_db[seq_name]['superfamily']
                    if order == tesort_db[seq_name]['order'] and superfamily_tesort != 'unknown':
                        print(">>>>>> Superfamily", superfamily,"<=", superfamily_tesort, seq_name, "is classified by TEsorter as", tesort_db[seq_name])
                        superfamily = superfamily_tesort


            family_hierarchy = {
                'class': class_,
                'order': order,
                'sFamily': superfamily,
                'length': length
            }

            te_families[seq_name] = family_hierarchy
    return te_families



# main
te_families = load_te_families(classif_file, load_te_sorter_classification(tesorter_file))
print(type(te_families))

with open(consensus_file, 'r') as f_consensus, \
     open(output_file, 'w') as f_output, \
     open(output_file_RM, 'w') as f_output_rm:
         for line in f_consensus:
            line = line.strip()
            if line.startswith('>'):
                # e.g. >MCLCluster449Mb2_PCOM_TEdenovoGr-B-G1-Map20
                fa_id = line.strip().split()[0][1:] 
                
                # e.g. PCOM_TEdenovoGr-B-G1-Map20
                short_id = '_'.join(fa_id.split('_')[1:])
    
                d = te_families.get(short_id)
                
                classification = f"{d.get('class', 'NA')}/{d.get('order', 'NA')}/{d.get('sFamily', 'NA')}"
                
                id_classification = f"{short_id}\t{classification}\t."
                print(f'>{id_classification}')
                
                f_output.write(f'>{id_classification}\n')
                f_output_rm.write(f'>{short_id}#{classification}\n')
            else:
                f_output.write(line + '\n')
                f_output_rm.write(line + '\n')

# Function to wrap text to 79 characters
def wrap79(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)






end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))