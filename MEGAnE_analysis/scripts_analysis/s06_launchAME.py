#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create upstream FASTA files for TFBS enrichment analysis.

From MEGAnE_summary_{SPECIES}.tsv:
  - Selected genes WITH TIPs
  - Selected genes WITHOUT TIPs
  - Non-selected genes

Also count genes with / without TIPs per condition.
"""

import os
import pandas as pd
import subprocess

# ======================================================
# PARAMETERS
# ======================================================

SPECIES = "PCOM"   # "PPY" or "PCOM"

BASE_DIR = "../output"
TFBS_DIR = f"{BASE_DIR}/{SPECIES}/TFBS"
os.makedirs(TFBS_DIR, exist_ok=True)

SUMMARY_TSV = f"{BASE_DIR}/{SPECIES}/MEGAnE_summary_{SPECIES}.tsv"

GENOME_FASTA = (
    f"../data/{SPECIES}/assembly/GWHBAOS00000000.genome.fasta.reformat"
    if SPECIES == "PPY"
    else f"../data/{SPECIES}/assembly/PyrusCommunis_BartlettDHv2.0.fasta"
)

UPSTREAM_LEN = 2000

MOTIF_DB = "../data/motifs/JASPAR2022_plants.meme" # wget https://jaspar.elixir.no/download/data/2022/CORE/JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt -O JASPAR2022_plants.meme

# Output files
BED_FILES = {
    "sel_with_TIP": f"{TFBS_DIR}/selected_with_TIP.bed",
    "sel_no_TIP": f"{TFBS_DIR}/selected_no_TIP.bed",
    "non_selected": f"{TFBS_DIR}/non_selected.bed",
}

FASTA_FILES = {
    k: v.replace(".bed", ".fa")
    for k, v in BED_FILES.items()
}

COUNT_FILE = f"{TFBS_DIR}/gene_counts.txt"

AME_DIR = f"{TFBS_DIR}/AME"
os.makedirs(AME_DIR, exist_ok=True)

# ======================================================
# FUNCTIONS
# ======================================================

def load_and_expand_genes(df):
    """
    Extract gene-level information and explode multi-gene fields.
    Remove TIPs without associated genes.
    """

    gdf = df[
        [
            "CHROM",
            "gene_Gene_ID",
            "gene_Start",
            "gene_End",
            "gene_Strand",
            "positive_selection",
        ]
    ].copy()

    # Remove TIPs without genes
    gdf = gdf[gdf["gene_Gene_ID"] != "NA"]

    # Split multi-gene fields
    gdf = gdf.assign(
        gene_Gene_ID=gdf["gene_Gene_ID"].str.split(","),
        gene_Start=gdf["gene_Start"].str.split(","),
        gene_End=gdf["gene_End"].str.split(","),
        gene_Strand=gdf["gene_Strand"].str.split(","),
        positive_selection=gdf["positive_selection"].str.split(","),
    ).explode(
        ["gene_Gene_ID", "gene_Start", "gene_End", "gene_Strand", "positive_selection"]
    )

    # Drop invalid rows
    gdf = gdf.dropna(
        subset=["gene_Gene_ID", "gene_Start", "gene_End", "gene_Strand"]
    )

    gdf["gene_Start"] = gdf["gene_Start"].astype(int)
    gdf["gene_End"] = gdf["gene_End"].astype(int)

    return gdf


def compute_upstream_coordinates(gdf):
    """
    Compute strand-aware upstream regions (2 kb).
    """

    def upstream(row):
        if row["gene_Strand"] == "+":
            start = max(1, row["gene_Start"] - UPSTREAM_LEN)
            end = row["gene_Start"]
        else:
            start = row["gene_End"]
            end = row["gene_End"] + UPSTREAM_LEN
        return pd.Series([start, end])

    gdf[["up_start", "up_end"]] = gdf.apply(upstream, axis=1)
    return gdf


def write_bed(df, out_bed):
    """
    Write BED file (0-based start).
    """

    bed = pd.DataFrame({
        "chrom": df["CHROM"],
        "start": df["up_start"] - 1,
        "end": df["up_end"],
        "name": df["gene_Gene_ID"],
        "score": 0,
        "strand": df["gene_Strand"],
    })

    bed = bed.drop_duplicates(subset=["name"])
    bed.to_csv(out_bed, sep="\t", header=False, index=False)


def bed_to_fasta(bed, fasta):
    """
    Extract sequences using bedtools getfasta.
    """

    cmd = (
        f"module load all gencore/3 && module load bedtools/2.31.1 && "
        f"bedtools getfasta -fi {GENOME_FASTA} -bed {bed} -s -nameOnly > {fasta}"
    )
    subprocess.run(cmd, shell=True, check=True)

def run_ame(target_fa, background_fa, out_dir, label):
    """
    Run AME (MEME suite) for TFBS enrichment.
    """

    os.makedirs(out_dir, exist_ok=True)

    cmd = (
        f"module load all gencore/3 && "
        f"module load meme/5.5.8 && "
        f"ame "
        f"--verbose 3 "
        f"--oc {out_dir} "
        f"--scoring avg "
        f"--method fisher "
        f"--hit-lo-fraction 0.25 "
        f"--evalue-report-threshold 10.0 "
        f"--control {background_fa} "
        f"{target_fa} "
        f"{MOTIF_DB}"
    )
    print(f"🔥 Running AME: {label}")
    subprocess.run(cmd, shell=True, check=True)

# ======================================================
# MAIN
# ======================================================

def main():

    print("🔹 Loading summary...")
    df = pd.read_csv(SUMMARY_TSV, sep="\t")

    print("🔹 Extracting gene-level table...")
    gdf = load_and_expand_genes(df)
    gdf = compute_upstream_coordinates(gdf)

    # -------------------------------
    # Gene sets
    # -------------------------------

    sel = gdf[gdf["positive_selection"].str.lower() == "true"]
    non_sel = gdf[gdf["positive_selection"].str.lower() == "false"]

    genes_with_TIP = set(sel["gene_Gene_ID"])
    all_selected_genes = set(sel["gene_Gene_ID"])

    # Selected genes WITHOUT TIPs
    sel_no_tip = pd.read_csv(
        f"../data/{SPECIES}/genes/positive_selection_summary_table.East.tsv"
        if SPECIES == "PPY"
        else f"../data/{SPECIES}/genes/positive_selection_summary_table.West.upstream1kb.2025-09-02.txt",
        sep="\t"
    )

    sel_no_tip = sel_no_tip[
        ~sel_no_tip["Gene_ID"].isin(genes_with_TIP)
    ]

    # -------------------------------
    # Write BED + FASTA
    # -------------------------------

    print("🔹 Writing FASTA files...")

    write_bed(sel, BED_FILES["sel_with_TIP"])
    bed_to_fasta(BED_FILES["sel_with_TIP"], FASTA_FILES["sel_with_TIP"])

    write_bed(non_sel, BED_FILES["non_selected"])
    bed_to_fasta(BED_FILES["non_selected"], FASTA_FILES["non_selected"])

    # Selected genes without TIPs (reuse upstream from table)
    sel_no_tip_bed = sel_no_tip.rename(
        columns={
            "Chr_accession": "chrom",
            "Start_upstream2Kb": "start",
            "End_upstream2Kb": "end",
            "Strand": "strand",
            "Gene_ID": "name",
        }
    )

    sel_no_tip_bed["start"] -= 1
    sel_no_tip_bed[["chrom", "start", "end", "name", "strand"]] \
        .assign(score=0) \
        .to_csv(BED_FILES["sel_no_TIP"], sep="\t", header=False, index=False)

    bed_to_fasta(BED_FILES["sel_no_TIP"], FASTA_FILES["sel_no_TIP"])

    # -------------------------------
    # 🔥 COUNT GENES
    # -------------------------------

    with open(COUNT_FILE, "w") as fh:
        fh.write("Gene counts per condition\n\n")
        fh.write(f"Selected genes WITH TIPs : {len(genes_with_TIP)}\n")
        fh.write(f"Selected genes WITHOUT TIPs : {len(sel_no_tip)}\n")
        fh.write(f"Non-selected genes : {non_sel['gene_Gene_ID'].nunique()}\n")

    print("✅ FASTA files written in:", TFBS_DIR)
    print("🔥 Gene counts written →", COUNT_FILE)

    # -------------------------------
    # 🔥 RUN AME
    # -------------------------------

    run_ame(
        FASTA_FILES["sel_with_TIP"],
        FASTA_FILES["non_selected"],
        f"{AME_DIR}/selected_with_TIP_vs_non_selected",
        "Selected WITH TIP vs Non-selected"
    )

    run_ame(
        FASTA_FILES["sel_no_TIP"],
        FASTA_FILES["non_selected"],
        f"{AME_DIR}/selected_no_TIP_vs_non_selected",
        "Selected WITHOUT TIP vs Non-selected"
    )

if __name__ == "__main__":
    main()
