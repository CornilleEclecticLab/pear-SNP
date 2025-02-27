
#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @File     : s02.prepare_GO_KEGG.py
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2025/02/11 11:44:37
# @Description:
#    

import datetime
import sys
import textwrap
import os
import glob

from s00_config import *

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


anno_files = glob.glob(f"{output_dir}/*.annotations")

# Read reference GO terms
ref_go_hash = {}
with open(ref_go_f, "r") as in_ref_go:
    for line in in_ref_go:
        if line.startswith("!"):
            continue
        arr = line.split("\t")
        go_term = arr[4]
        ref_go_hash[go_term] = 1

# Write reference GO lists to output files
with open(os.path.join(sub_output_dir, "ref.go.list.txt"), "w") as ou_rgo:
    for go in sorted(ref_go_hash.keys()):
        ou_rgo.write(f"{go}\n")


# Read reference KEGG terms
ref_kegg_hash = {}
with open(ref_kegg_f, "r") as in_ref_kegg, \
        open(os.path.join(sub_output_dir, "KEGG_TERM2NAME.table.txt"), "w") as ou_kegg_name:

    # Write header to output file
    ou_kegg_name.write("ko\tKEGG_pathway\n")

    for line in in_ref_kegg:
        line = line.strip()
        if line == "":
            continue
        if line[0] not in '0123456789':
            continue
        if not line.startswith("0"):
            print(f"WARNING: Unexpected line not starts with 0: {line}")

        arr = line.split("  ")
        if len(arr[0]) != 5:
            raise ValueError(f"KO ID is not 5 characters: {line}")

        ko = "ko" + arr[0]
        ref_kegg_hash[ko] = 1
        print(f"{ko} :: {arr[1]}")
        ou_kegg_name.write(f"{ko}\t{arr[1]}\n")

# Write reference KEGG lists to output files
with open(os.path.join(sub_output_dir, "ref.kegg.list.txt"), "w") as ou_rkegg:
    for kegg in sorted(ref_kegg_hash.keys()):
        ou_rkegg.write(f"{kegg}\n")


# Read annotation file and process GO terms
for anno_file in anno_files:
    basename = os.path.basename(anno_file)
    print(f"Processing {basename} ...")
    prefix = basename.replace(".emapper.annotations", "")
    sp_go_hash = {}
    sp_kegg_hash = {}
    with open(anno_file, "r") as in1, \
            open(os.path.join(sub_output_dir, f"{prefix}.GO_GID.table.txt"), "w") as ou_gene_id,\
            open(os.path.join(sub_output_dir, f"{prefix}.GO_GO_term.table.txt"), "w") as ou_go_term,\
            open(os.path.join(sub_output_dir, f"{prefix}.KEGG_TERM2GENE.table.txt"), "w") as ou_kegg_ko:

        # Write headers to output files
        ou_gene_id.write("GID\tSYMBOL\tGENENAME\n")
        ou_go_term.write("GID\tGO\tEVIDENCE\n")
        ou_kegg_ko.write("ko\tGID\n")

        for line in in1:
            line = line.strip()
            if line.startswith("#"):
                continue
            arr = line.split("\t")

            if "Gene=" in arr[0]:
                gene_id = arr[0].split("Gene=")[1].split("\\t")[0]
            else:
                gene_id = arr[0].split("\\t")[0]

            gene_discription = arr[7]
            gene_name = arr[8]
            go_line = arr[9]
            kegg_line = arr[12]

            if kegg_line != "-":
                kos = kegg_line.split(",")
                n_ko, n_tko = 0, 0
                for ko in kos:
                    if "map" in ko:
                        continue
                    n_tko += 1
                    sp_kegg_hash[ko] = 1
                    if ko in ref_kegg_hash:
                        n_ko += 1
                        ou_kegg_ko.write(f"{ko}\t{gene_id}\n")
                print(f"{gene_id}:\ttotal_ko={n_tko}\tannotation_ko={n_ko}")

            if go_line != "-":
                gos = go_line.split(",")

                ############################################
                # TODO: change the same gene_id three times?
                ############################################
                ou_gene_id.write(f"{gene_id}\t{gene_id}\t{gene_id}\n")
                n_go = 0
                for go_term in gos:
                    sp_go_hash[go_term] = 1
                    if go_term in ref_go_hash:
                        n_go += 1

                        ou_go_term.write(f"{gene_id}\t{go_term}\tegg\n")
                print(f"{gene_id}:\ttotal_go={len(gos)}\tannotation_go={n_go}")

    # Write GO lists to output files
    with open(os.path.join(sub_output_dir, f"{prefix}.go.list.txt"), "w") as ou_sgo:
        for go in sorted(sp_go_hash.keys()):
            ou_sgo.write(f"{go}\n")

    # Write KEGG lists to output files
    with open(os.path.join(sub_output_dir, f"{prefix}.kegg.list.txt"), "w") as ou_skegg:
        for kegg in sorted(sp_kegg_hash.keys()):
            ou_skegg.write(f"{kegg}\n")


end_time = datetime.datetime.now()
print('')
print(' END '.center(79, '='))
print(str(end_time-start_time).center(79))
