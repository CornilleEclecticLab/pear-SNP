#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build MEGAnE_summary_[PCOM/PPY].tsv

For each transposable element (either MEI or MEA):
  - Extract CHROM, POS, ID
  - List samples carrying the TE (duplicated if homozygous)
  - Associate samples with populations (via individuals_table)
  - Include TE classification info (class, order, family)
  - If the TE is located within ≤2kb upstream of any gene (strand-aware):
      Add gene information (Gene_ID, Start, End, Strand, Symbol, positive_selection)
"""

import os
import gzip
import warnings
import pandas as pd
import numpy as np
import re
from collections import Counter
import pyranges as pr
import subprocess
from goatools.obo_parser import GODag
import matplotlib.pyplot as plt
import seaborn as sns


warnings.filterwarnings("ignore")

# ------------------- File paths -------------------

#SPECIES = 'PCOM'
SPECIES = 'PPY'

VCF_MEI = f"../output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEI_jointcall.vcf.gz"
VCF_MEA = f"../output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEA_jointcall.vcf.gz"

if SPECIES == 'PCOM':
    CLASSIF = f"../data/{SPECIES}/TE_annotation_URGI/{SPECIES}_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
    GENE_TSV_SELECTED = f"../data/{SPECIES}/genes/positive_selection_summary_table.West.upstream1kb.2025-09-02.txt" #2,156 genes Source: /scratch/yn2515/pear_snp_tip/selection_gene_list_2025-09-02/positive_selection_summary_table.West.upstream1kb.2025-09-02.txt
    GENES_GFF = f"../data/{SPECIES}/genes/PyrusCommunis_BartlettDHv2.0.gff" # Source: /scratch/yn2515/pear_snp_tip/assembly/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.gff
    GENES_GFF_CLEAN = f"../data/{SPECIES}/genes/PyrusCommunis_BartlettDHv2.0-rmIntrons.gff"
    POP_FILE = f"../data/{SPECIES}/s01.individuals_table.txt"
    RHO_FILE = f"../data/{SPECIES}/Rho/FastEPRR_step4.comm_Dessert.Chrs.withChr.cat.bed"
    GENOME_FILE = f"../data/{SPECIES}/assembly/PyrusCommunis_BartlettDHv2.0.fasta"
    TES_GFF = f"../data/{SPECIES}/TE_annotation_URGI/PCOM_TE_annotation.gff3"
elif SPECIES == 'PPY':
    CLASSIF = f"../data/{SPECIES}/TE_annotation_URGI/{SPECIES}3_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
    GENE_TSV_SELECTED = f"../data/{SPECIES}/genes/positive_selection_summary_table.East.tsv" # 3,527 genes, Source: /scratch/yn2515/pear_snp_tip/selection_gene_list_2025-09-02/positive_selection_summary_table.East.tsv
    GENES_GFF = f"../data/{SPECIES}/genes/GWHBAOS00000000.gff" # Source: /scratch/yn2515/pear_snp_tip/assembly/Pyrus_pyrifolia_Cuiguan_Gao2021/GWHBAOS00000000.gff
    GENES_GFF_CLEAN = f"../data/{SPECIES}/genes/GWHBAOS00000000.gff"
    POP_FILE = f"../data/{SPECIES}/s01.individuals_table.txt3"
    RHO_FILE = f"../data/{SPECIES}/Rho/FastEPRR_step4.pyri_JP.Chrs.withChr.cat.bed"
    GENOME_FILE = f"../data/{SPECIES}/assembly/GWHBAOS00000000.genome.fasta.reformat"
    TES_GFF = f"../data/{SPECIES}/TE_annotation_URGI/PPY_TE_annotation.gff3.renameChr"

OUT_DIR = f"../output/{SPECIES}"
OUT_TSV = f"../output/{SPECIES}/MEGAnE_summary_{SPECIES}.tsv"

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

# ================================================================
# 1️⃣ LOAD METADATA
# ================================================================

def load_classification(path):
    """
    Load TE classification information.

    Returns dictionaries mapping TE IDs to:
      - class
      - order
      - sFamily
    """
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("Unclassified")
    if "sFamily" in df.columns:
        df["sFamily"] = df["sFamily"].replace("NA", "Unclassified")
    class_dict = dict(zip(df["Seq_name"], df["class"]))
    order_dict = dict(zip(df["Seq_name"], df["order"]))
    fam_dict = dict(zip(df["Seq_name"], df.get("sFamily", ["Unclassified"] * len(df))))
    return class_dict, order_dict, fam_dict

def load_populations(path):
    """
    Load individuals table and return a mapping:
        sample_id (from VCF) → population name
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    if "#ID_vcf" in df.columns:
        df = df.rename(columns={"#ID_vcf": "ID_vcf"})
    df["Population"] = df["Population"].fillna("Unknown")
    return df.set_index("ID_vcf")["Population"].to_dict()

# ================================================================
# 2️⃣ GENE INFO HANDLING
# ================================================================

def filter_gene_gff(gff_in, gff_out):
    """
    Keep only valid gene features from a GFF file.
    """

    df = pd.read_csv(
        gff_in,
        sep="\t",
        comment="#",
        header=None,
        dtype=str
    )

    df = df[df[2] == "gene"]

    df[3] = df[3].astype(int)
    df[4] = df[4].astype(int)

    df = df[df[3] <= df[4]]

    df.to_csv(gff_out, sep="\t", header=False, index=False)

    print(f"[OK] Gene-only GFF written → {gff_out}")

def load_genes_selected( tsv, SPECIES ):
    """
    Load genes with their precomputed 2kb upstream coordinates
    from the provided TSV (strand-aware).

    Required columns:
      Gene_ID,
      Chr_ID,
      Start,
      End,
      Strand,
      Start_upstream2Kb,
      End_upstream2Kb,
      Symbol (optional),
      GO_ID,
      Ath_description
    """
    
    if SPECIES == 'PCOM':
        chr_name = "Chr_ID"
    elif SPECIES == 'PPY':
        chr_name = "Chr_accession"
        
    df = pd.read_csv(tsv, sep="\t", dtype=str).fillna("NA")
    genes = []
    for _, r in df.iterrows():
        try:
            g = {
                "Gene_ID": r["Gene_ID"],
                "Chr": r[chr_name],
                "Start": int(r["Start"]),
                "End": int(r["End"]),
                "Strand": r["Strand"],
                "up1": int(r["Start_upstream2Kb"]),
                "up2": int(r["End_upstream2Kb"]),
                "Symbol": r.get("Symbol", "NA"),
                "positive_selection": "true",
                "GO_ID": r.get("GO_ID", "NA"),
                "Ath_description": r.get("Ath_description", "NA")
            }
            #Ensure coordinates are in correct order (up1 < up2)
            if g["up1"] > g["up2"]:
                g["up1"], g["up2"] = g["up2"], g["up1"]
            genes.append(g)
        except Exception:
            pass
    return genes

def load_genes_nonselected(gff, selected_ids):
    genes = []
    with open(gff) as fh:
        for l in fh:
            if l.startswith("#") or not l.strip():
                continue
            c = l.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            start, end = int(c[3]), int(c[4])
            strand = c[6]
            gid = ""
            for p in c[8].split(";"):
                if p.startswith("ID="):
                    v = p.split("=", 1)[1]
                    gid = v.split("gene:", 1)[1] if v.startswith("gene:") else v
            if not gid or gid in selected_ids:
                continue
            if strand == "+":
                u1 = max(1, start - 2000)
                u2 = start
            else:
                u1 = end
                u2 = end + 2000
            genes.append({
                "Gene_ID": gid,
                "Chr": c[0], "Start": start, "End": end,
                "Strand": strand, "up1": u1, "up2": u2,
                "Symbol": "NA",
                "positive_selection": "false",
                "GO_ID": "NA",
                "Ath_description": "NA"
            })
    return genes

def build_upstream_index(genes):
    """
    Build a simple index of genes grouped by chromosome.
    Each entry contains the list of genes and their upstream intervals.
    """
    idx = {}
    for g in genes:
        idx.setdefault(g["Chr"], []).append(g)
    return idx


def find_upstream_genes(chr_, pos, idx):
    """
    Return a list of genes whose upstream region (2kb)
    overlaps with the given genomic position.
    """
    if chr_ not in idx:
        return []
    return [g for g in idx[chr_] if g["up1"] <= pos <= g["up2"]]

# ================================================================
# 3️⃣ PARSE VCF FILES
# ================================================================

def parse_vcf(path, te_type, pop_dict, cdict, odict, fdict, upstream_idx):
    """
    Parse a compressed VCF file and extract all relevant TE information.

    For each variant:
      - Extract CHROM, POS, ID
      - Determine TE consensus IDs and their classification
      - Determine samples and populations carrying the TE
      - Check for upstream genes within ≤2kb
    """
    rows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            # Skip metadata lines
            if line.startswith("##"):
                continue

            # Extract sample names from header
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                raw_samples = header[9:]
                # Clean sample names: keep only the third part if formatted like "X.Y.Z"
                samples = [s.split(".")[2] if len(s.split(".")) >= 3 else s for s in raw_samples]
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, TE_id, INFO, FORMAT = parts[0], int(parts[1]), parts[2], parts[7], parts[8]
            genos = parts[9:]
            if "GT" not in FORMAT.split(":"): # extract Genotype info
                continue
            gt_idx = FORMAT.split(":").index("GT")

            # ---- Extract END position (from INFO, e.g. "0END=110545") ----
            end = pos  # default value if not found
            for kv in INFO.split(";"):
                if kv.startswith("0END="):
                    try:
                        end = int(kv.split("=")[1])
                    except ValueError:
                        end = pos
                    break

            # ---- Extract SVLEN ----
            svlen = "NA"
            for kv in INFO.split(";"):
                if kv.startswith("SVLEN="):
                    try:
                        svlen = int(kv.split("=")[1])
                    except:
                        svlen = "NA"
                    break

            # ---- Extract TE IDs from INFO ----
            mei_ids = []
            for kv in INFO.split(";"):
                if kv.startswith("MEI="):
                    mei_ids = kv.split("=")[1].split("|")
                    break

            # ---- Classification fields ----
            cons = ",".join(mei_ids) if mei_ids else "NA"
            nte = len(mei_ids)
            cls = ",".join([cdict.get(i, "Unclassified") for i in mei_ids])
            odr = ",".join([odict.get(i, "Unclassified") for i in mei_ids])
            fam = ",".join([fdict.get(i, "Unclassified") for i in mei_ids])

            # ---- Genotype parsing (sample presence/absence) ----
            s_list, p_list = [], []
            for s_clean, g in zip(samples, genos):
                g_info = g.split(":")
                if gt_idx >= len(g_info):
                    continue
                gt = g_info[gt_idx]
                if gt in (".", "./.", ".|."):
                    continue

                # Convert genotype to allele dosage
                vals = [int(a) for a in gt.replace("|", "/").split("/") if a != "."]
                dosage = sum(vals) if te_type == "MEI" else (2 - sum(vals))

                # Duplicate sample entries for homozygotes
                if dosage > 0:
                    for _ in range(int(dosage)):
                        s_list.append(s_clean)
                        p_list.append(pop_dict.get(s_clean, "Unknown"))

            # ---- Check if the TE position falls within any gene's 2kb upstream region ----
            up_genes = find_upstream_genes(chrom, pos, upstream_idx)
            if up_genes:
                g_ids = [g["Gene_ID"] for g in up_genes]
                starts = [str(g["Start"]) for g in up_genes]
                ends = [str(g["End"]) for g in up_genes]
                strands = [g["Strand"] for g in up_genes]
                symbols = [g["Symbol"] for g in up_genes]
                #ps = [str(g["positive_selection"]).lower() for g in up_genes]
                ps = [g["positive_selection"] for g in up_genes]
                GO_ID = [g["GO_ID"] for g in up_genes]
                Ath_description = [g["Ath_description"] for g in up_genes]
            else:
                g_ids = starts = ends = strands = symbols = ps = GO_ID = Ath_description = ["NA"]

            # ---- Build the result row ----
            rows.append({
                "CHROM": chrom,
                "POS": pos,
                "END": end,
                "END_tmp": pos + 1,
                "SVLEN": svlen,
                "ID": TE_id,
                "samples": ",".join(s_list) if s_list else "NA",
                "populations": ",".join(p_list) if p_list else "NA",
                "nb_TEs": nte,
                "consensus": cons,
                "class": cls,
                "order": odr,
                "family": fam,
                "gene_Gene_ID": ",".join(map(str, g_ids)),
                "gene_Start": ",".join(map(str, starts)),
                "gene_End": ",".join(map(str, ends)),
                "gene_Strand": ",".join(map(str, strands)),
                "gene_Symbol": ",".join(map(str, symbols)),
                "positive_selection": ",".join(map(str, ps)),
                "GO_ID": ",".join(map(str, GO_ID)),
                "Ath_description": ",".join(map(str, Ath_description)),
            })
    return rows


def _normalize_token(tok: str) -> str:
    t = str(tok).strip()
    if "|" in t:
        t = t.split("|", 1)[0].strip()  # 1) garder avant le premier pipe
    if t in {"", "NA", "NaN", "nan", "N/A", "."}:
        t = "Unclassified"              # 2) normaliser
    return t
    
def _choose_dominant(cell: str) -> str:
    if not isinstance(cell, str):
        cell = "" if pd.isna(cell) else str(cell)
    raw_tokens = [t.strip() for t in cell.split(",") if t.strip() != ""]
    if not raw_tokens:
        return "Unclassified"
    tokens = [_normalize_token(t) for t in raw_tokens]

    # 3) préférer la valeur non-Unclassified majoritaire
    non_uncl = [t for t in tokens if t != "Unclassified"]
    if non_uncl:
        counts = Counter(non_uncl)
        max_count = max(counts.values())
        # stabilité : on garde la première parmi les plus fréquentes
        for t in tokens:
            if t != "Unclassified" and counts[t] == max_count:
                return t
    return "Unclassified"
    
def _replace_family_NA_tokens(s: str) -> str:
    # remplace tout token 'NA' dans family par 'Unclassified' (prend en compte ',' et '|')
    if not isinstance(s, str):
        s = "" if pd.isna(s) else str(s)
    parts = re.split(r'([,\|])', s)  # on garde les délimiteurs
    out = []
    for p in parts:
        if p in {",", "|"}:
            out.append(p)
        else:
            tok = p.strip()
            out.append("Unclassified" if tok == "NA" else p)
    return "".join(out)


def add_rho_v1(df):
    
    #print(df)
    df_rho = pd.read_csv(RHO_FILE, sep="\t")
    #print(df_rho)

    # Créer des objets PyRanges
    gr1 = pr.PyRanges(df.rename(columns={'CHROM': 'Chromosome', 'POS': 'Start', 'END_tmp': 'End'}))
    gr2 = pr.PyRanges(df_rho.rename(columns={'#Chr': 'Chromosome'}))
    
    # Intersection : trouve les lignes de df1 qui tombent dans une fenêtre de df2
    #merged = gr1.join(gr2)
    merged = gr1.join(gr2, how="left")
    
    #merge
    df_merged = merged.as_df()

    df_merged = df_merged.rename(columns={'Chromosome': 'CHROM', 'Start': 'POS', 'End': 'END_tmp', 'Start_b': 'Start_Rho', 'End_b': 'End_Rho'})

    #print(df_merged.head())
    
    return df_merged

def add_rho(df):

    df_rho = pd.read_csv(RHO_FILE, sep="\t")

    gr1 = pr.PyRanges(df.rename(columns={
        'CHROM': 'Chromosome',
        'POS': 'Start',
        'END_tmp': 'End'
    }))
    gr2 = pr.PyRanges(df_rho.rename(columns={'#Chr': 'Chromosome'}))

    merged = gr1.join(gr2, how="left")
    df_new = merged.as_df()

    df_new = df_new.rename(columns={
        'Chromosome': 'CHROM',
        'Start': 'POS',
        'End': 'END_tmp',
        'Start_b': 'Start_Rho',
        'End_b': 'End_Rho'
    })

    # 🔑 ICI on ajoute AUSSI Rho
    df = df.merge(
        df_new[['POS', 'ID', 'Start_Rho', 'End_Rho', 'Rho']],
        on=['POS', 'ID'],
        how='left'
    )

    return df




def add_TEdensity_v1(df, gff_file, genome_file, OUT_DIR, window_size):
    
    print(df)
    
    out_prefix="TE_density"

    # Créer les fichiers temporaires
    windows_file = f"{OUT_DIR}/{out_prefix}_100kb.bed"
    out_file = f"{OUT_DIR}/{out_prefix}_100kb_cov.bed"

    # Commande bedtools makewindows
    faidx_cmd = f"module load all gencore/2 && module load samtools/1.9 && samtools faidx {genome_file}"
    subprocess.run(faidx_cmd, shell=True, check=True)
    
    # Commande bedtools makewindows
    makewindows_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools makewindows -g {genome_file}.fai -w {window_size} > {windows_file}"
    subprocess.run(makewindows_cmd, shell=True, check=True)

    # Commande bedtools coverage
    coverage_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools coverage -a {windows_file} -b {gff_file} > {out_file}"
    subprocess.run(coverage_cmd, shell=True, check=True)

    print(f"[OK] Densité TE calculée → {out_file}")

    # Charger le résultat dans pandas si tu veux l’ajouter à df :
    df_te_cov = pd.read_csv(out_file, sep="\t", header=None)
    df_te_cov.columns = ["Chromosome", "Start", "End", "NB_TE_OVERLAP", "BASES_COVERED", "WINDOW_SIZE", "TE_density"]
    print(df_te_cov)
    
    # Créer des objets PyRanges
    gr1 = pr.PyRanges(df.rename(columns={'CHROM': 'Chromosome', 'POS': 'Start', 'END_tmp': 'End'}))
    gr2 = pr.PyRanges(df_te_cov)
    # Intersection : trouve les lignes de df1 qui tombent dans une fenêtre de df2
    #merged = gr1.join(gr2)
    merged = gr1.join(gr2, how="left")
    #merge
    df_merged = merged.as_df()
    df_merged = df_merged.rename(columns={'Chromosome': 'CHROM', 'Start': 'POS', 'End': 'END_tmp', 'Start_b': 'Start_TEcov', 'End_b': 'End_TEcov'})
    print(df_merged)
    
    return df_merged


def add_TEdensity(df, gff_file, genome_file, OUT_DIR, window_size):

    print(df)
    
    out_prefix="TE_density"

    # Créer les fichiers temporaires
    windows_file = f"{OUT_DIR}/{out_prefix}_100kb.bed"
    out_file = f"{OUT_DIR}/{out_prefix}_100kb_cov.bed"

    # Commande bedtools makewindows
    faidx_cmd = f"module load all gencore/2 && module load samtools/1.9 && samtools faidx {genome_file}"
    subprocess.run(faidx_cmd, shell=True, check=True)
    
    # Commande bedtools makewindows
    makewindows_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools makewindows -g {genome_file}.fai -w {window_size} > {windows_file}"
    subprocess.run(makewindows_cmd, shell=True, check=True)

    # Commande bedtools coverage
    coverage_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools coverage -a {windows_file} -b {gff_file} > {out_file}"
    subprocess.run(coverage_cmd, shell=True, check=True)

    print(f"[OK] Densité TE calculée → {out_file}")
    
    df_te_cov = pd.read_csv(out_file, sep="\t", header=None)
    df_te_cov.columns = [
        "Chromosome", "Start", "End",
        "NB_TE_OVERLAP", "BASES_COVERED",
        "WINDOW_SIZE", "TE_density"
    ]

    gr1 = pr.PyRanges(df.rename(columns={
        'CHROM': 'Chromosome',
        'POS': 'Start',
        'END_tmp': 'End'
    }))
    gr2 = pr.PyRanges(df_te_cov)

    merged = gr1.join(gr2, how="left")
    df_new = merged.as_df()

    df_new = df_new.rename(columns={
        'Chromosome': 'CHROM',
        'Start': 'POS',
        'End': 'END_tmp',
        'Start_b': 'Start_TEcov',
        'End_b': 'End_TEcov'
    })

    df = df.merge(
        df_new[['POS', 'ID',
                'NB_TE_OVERLAP', 'BASES_COVERED',
                'WINDOW_SIZE', 'TE_density',
                'Start_TEcov', 'End_TEcov']],
        on=['POS', 'ID'],
        how='left'
    )

    return df



def add_GO_level0( df ):
    """
    Add a column with GO Level-0 (root) categories for each gene.
    Level 0 corresponds to the three main GO domains:
    'biological process', 'molecular function', or 'cellular component'.
    """

    # Load GO structure
    go_dag = GODag("go-basic.obo")

    # Function to get the Level-0 ancestors (root categories)
    def get_level0_terms(go_list):
        if pd.isna(go_list):
            return []
        gos = [g.strip() for g in go_list.split(",") if g.startswith("GO:")]
        level0 = set()
        for go_id in gos:
            if go_id not in go_dag:
                continue
            go_term = go_dag[go_id]

            # Climb up to ancestors until we reach level 0 (root)
            ancestors = go_term.get_all_parents()
            for a in ancestors:
                if a in go_dag and go_dag[a].level == 0:
                    level0.add(go_dag[a].name)

            # If itself is level 0
            if go_term.level == 0:
                level0.add(go_term.name)

        return list(level0)

    # Apply to GO column
    df["GO_Level0"] = df["GO_ID"].apply(get_level0_terms)
    return df


def classify_go_terms_old(go_terms):
    # Define categories
    developmental_keywords = [
        "develop", "signal", "morphogen", "differentiation", "growth",
        "pattern", "embryo", "organogenesis", "morphogenesis",
        "cell fate", "cell communication", "receptor", "transduction",
        "cell migration", "cell proliferation"
    ]
    immunity_keywords = [
        "immune", "defense", "response to", "inflammatory", "antigen",
        "pathogen", "infection", "cytokine", "chemokine", "leukocyte",
        "macrophage", "innate", "complement", "interferon"
    ]
    go_terms = [t.lower() for t in go_terms if isinstance(t, str)]
    
    # Category flags
    has_dev = any(any(k in t for k in developmental_keywords) for t in go_terms)
    has_imm = any(any(k in t for k in immunity_keywords) for t in go_terms)
    
    if has_imm and has_dev:
        return "Ambiguous"
    elif has_imm:
        return "Immunity"
    elif has_dev:
        return "Developmental & Signaling"
    else:
        return "Other"
        

def classify_go_terms(go_terms):
    """
    Classify GO terms into broad biological categories.
    Returns 'Ambiguous' if multiple categories apply.
    """
    if not go_terms:
        return "Other"

    go_terms = [t.lower() for t in go_terms if isinstance(t, str)]

    # --- Keyword groups ---
    developmental_keywords = [
        "develop", "signal", "morphogen", "differentiation", "growth",
        "pattern", "embryo", "organogenesis", "morphogenesis",
        "cell fate", "cell communication", "receptor", "transduction",
        "cell migration", "cell proliferation"
    ]
    immunity_keywords = [
        "immune", "inflammatory", "antigen", #"defense", 
        "pathogen", "infection", "cytokine", "chemokine", "leukocyte",
        "macrophage", "innate", "complement", "interferon", "disease"
    ]
    # metabolism_keywords = [
    #     "metabolic", "biosynthesis", "catabolic", "oxidation", "glycolysis", "anabolism"
    # ]
    # cellular_keywords = [
    #     "cell cycle", "cell division", "cellular process",
    #     "cell organization", "cell communication", "cell death"
    # ]
    response_keywords = [
        "response to", "stress", "stimulus", "detoxification"
    ]
    # function_keywords = [
    #     "binding", "catalytic", "transporter", "enzyme", "regulator", "ligase", "hydrolase"
    # ]
    # structure_keywords = [
    #     "organelle", "membrane", "complex", "cell part", "cytoskeleton"
    # ]

    # --- Detection flags ---
    categories = []
    if any(any(k in t for k in immunity_keywords) for t in go_terms):
        categories.append("Immunity")
    #if any(any(k in t for k in developmental_keywords) for t in go_terms):
    #    categories.append("Developmental & Signaling")
    # if any(any(k in t for k in metabolism_keywords) for t in go_terms):
    #     categories.append("Metabolism")
    #if any(any(k in t for k in response_keywords) for t in go_terms):
    #    categories.append("Response to Stimulus")
    # if any(any(k in t for k in cellular_keywords) for t in go_terms):
    #     categories.append("Cellular Process")
    # if any(any(k in t for k in function_keywords) for t in go_terms):
    #     categories.append("Molecular Function")
    # if any(any(k in t for k in structure_keywords) for t in go_terms):
    #     categories.append("Cell Structure")

    # --- Decision logic ---
    if len(set(categories)) == 0:
        return "Other"
    elif len(set(categories)) == 1:
        return categories[0]
    else:
        return "Ambiguous"


def add_GO_level2( df ):

    # Load GO structure
    go_dag = GODag("go-basic.obo")

    # Function to get level-2 ancestors
    def get_level2_terms(go_list):
        if pd.isna(go_list):
            return []
        gos = [g.strip() for g in go_list.split(",") if g.startswith("GO:")]
        level2 = set()
        for go_id in gos:
            if go_id not in go_dag:
                continue
            go_term = go_dag[go_id]
            # Climb up to ancestors until level 2
            ancestors = [a for a in go_term.get_all_parents() if a in go_dag]
            for a in ancestors:
                if go_dag[a].level == 2:
                    level2.add(go_dag[a].name)
            # Include itself if already level 2
            if go_term.level == 2:
                level2.add(go_term.name)
        return list(level2)
    
    # Apply to your GO column
    df["GO_Level2"] = df["GO_ID"].apply(get_level2_terms)
    df["GO_Level2_Category"] = df["GO_Level2"].apply(classify_go_terms)

    
    return df
    
    
# def is_NBS(df):
#     df['is_NBS'] = df['gene_Symbol'].apply(
#         lambda x: "NBS" if "NBS" in str(x) else "Other"
#     )
#     return df

# def is_Disease(df):
#     df['Ath_description_disease'] = df['Ath_description'].astype(str).str.contains(
#         "disease", case=False, na=False
#     ).map({True: "Disease", False: "Other"})
#     return df

def flag_immunity(df):
    # List of immunity-related keywords
    immunity_keywords = [
        "immune", "pathogen", "infection", "disease",
        #"cytokine", "chemokine", "leukocyte", "macrophage", "interferon", "inflammatory"
    ]

    # Condition 1: "disease" found in Ath_description (case-insensitive)
    cond_description = df['Ath_description'].astype(str).str.contains(
        r"disease|pathogen",
        case=False,
        na=False
    )

    # Condition 2: gene_symbol contains "NBS" (case-insensitive)
    cond_NBS = df['gene_Symbol'].astype(str).str.contains(
        "nbs", case=False, na=False
    )

    # Condition 3: GO_Level2 contains any of the immunity keywords
    cond_GO = df["GO_Level2"].astype(str).str.contains(
        "|".join(immunity_keywords),
        case=False,
        na=False
    )

    # Final immunity flag: True if ANY of the conditions is True
    df["Immunity"] = cond_description | cond_NBS | cond_GO

    return df

def add_gene_density_v1(df, genes_gff, genome_file, OUT_DIR, window_size):
    
    print(df)
    
    out_prefix = "gene_density"

    # Temporary files
    windows_file = f"{OUT_DIR}/{out_prefix}_100kb.bed"
    out_file = f"{OUT_DIR}/{out_prefix}_100kb_cov.bed"

    # Genome index
    faidx_cmd = f"module load all gencore/2 && module load samtools/1.9 && samtools faidx {genome_file}"
    subprocess.run(faidx_cmd, shell=True, check=True)
    
    # Genomic windows
    makewindows_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools makewindows -g {genome_file}.fai -w {window_size} > {windows_file}"
    subprocess.run(makewindows_cmd, shell=True, check=True)

    # Gene coverage
    coverage_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools coverage -a {windows_file} -b {genes_gff} > {out_file}"
    subprocess.run(coverage_cmd, shell=True, check=True)

    print(f"[OK] Gene density computed → {out_file}")

    # Load coverage results
    df_gene_cov = pd.read_csv(out_file, sep="\t", header=None)
    df_gene_cov.columns = [
        "Chromosome",
        "Start",
        "End",
        "NB_GENE_OVERLAP",
        "BASES_COVERED",
        "WINDOW_SIZE",
        "gene_density"
    ]
    print(df_gene_cov)
    
    # PyRanges join (same logic as TE density)
    gr1 = pr.PyRanges(df.rename(columns={
        "CHROM": "Chromosome",
        "POS": "Start",
        "END_tmp": "End"
    }))
    gr2 = pr.PyRanges(df_gene_cov)

    #merged = gr1.join(gr2)
    merged = gr1.join(gr2, how="left")

    df_merged = merged.as_df()
    df_merged = df_merged.rename(columns={
        "Chromosome": "CHROM",
        "Start": "POS",
        "End": "END_tmp",
        "Start_b": "Start_GENcov",
        "End_b": "End_GENcov"
    })

    print(df_merged)
    
    return df_merged

def add_gene_density(df, genes_gff, genome_file, OUT_DIR, window_size):

    print(df)
    
    out_prefix = "gene_density"

    # Temporary files
    windows_file = f"{OUT_DIR}/{out_prefix}_100kb.bed"
    out_file = f"{OUT_DIR}/{out_prefix}_100kb_cov.bed"

    # Genome index
    faidx_cmd = f"module load all gencore/2 && module load samtools/1.9 && samtools faidx {genome_file}"
    subprocess.run(faidx_cmd, shell=True, check=True)
    
    # Genomic windows
    makewindows_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools makewindows -g {genome_file}.fai -w {window_size} > {windows_file}"
    subprocess.run(makewindows_cmd, shell=True, check=True)

    # Gene coverage
    coverage_cmd = f"module load all gencore/3 && module load bedtools/2.31.1 && bedtools coverage -a {windows_file} -b {genes_gff} > {out_file}"
    subprocess.run(coverage_cmd, shell=True, check=True)

    print(f"[OK] Gene density computed → {out_file}")
    
    df_gene_cov = pd.read_csv(out_file, sep="\t", header=None)
    df_gene_cov.columns = [
        "Chromosome", "Start", "End",
        "NB_GENE_OVERLAP", "BASES_COVERED",
        "WINDOW_SIZE", "gene_density"
    ]

    gr1 = pr.PyRanges(df.rename(columns={
        'CHROM': 'Chromosome',
        'POS': 'Start',
        'END_tmp': 'End'
    }))
    gr2 = pr.PyRanges(df_gene_cov)

    merged = gr1.join(gr2, how="left")
    df_new = merged.as_df()

    df_new = df_new.rename(columns={
        'Chromosome': 'CHROM',
        'Start': 'POS',
        'End': 'END_tmp',
        'Start_b': 'Start_GENcov',
        'End_b': 'End_GENcov'
    })

    df = df.merge(
        df_new[['POS', 'ID',
                'NB_GENE_OVERLAP', 'BASES_COVERED',
                'WINDOW_SIZE', 'gene_density',
                'Start_GENcov', 'End_GENcov']],
        on=['POS', 'ID'],
        how='left'
    )

    return df




# ================================================================
# 4️⃣ MAIN EXECUTION PIPELINE
# ================================================================

def main():

    if not os.path.exists(GENES_GFF_CLEAN):
        filter_gene_gff(GENES_GFF, GENES_GFF_CLEAN)

    print("🔹 Loading metadata...")
    cdict, odict, fdict = load_classification(CLASSIF)
    pop_dict = load_populations(POP_FILE)

    print("🔹 Loading genes...")
    sel_genes = load_genes_selected( GENE_TSV_SELECTED, SPECIES )
    sel_ids = {g["Gene_ID"] for g in sel_genes}
    non_sel = load_genes_nonselected(GENES_GFF, sel_ids)
    all_genes = sel_genes + non_sel
    upstream_idx = build_upstream_index(all_genes)
    
    print("🔹 Parsing VCF files...")
    rows_mei = parse_vcf(VCF_MEI, "MEI", pop_dict, cdict, odict, fdict, upstream_idx)
    rows_mea = parse_vcf(VCF_MEA, "MEA", pop_dict, cdict, odict, fdict, upstream_idx)

    # Combine MEI + MEA results and sort by genomic position
    df = pd.DataFrame(rows_mei + rows_mea).sort_values(["CHROM", "POS"])
    
    # 4) Clean family name (NA becomes Unclassified)
    df["family"] = df["family"].apply(_replace_family_NA_tokens)
    
    # 5) Create 3 new columns with basic TE classif curation
    df["class_curated"]  = df["class"].apply(_choose_dominant)
    df["order_curated"]  = df["order"].apply(_choose_dominant)
    df["family_curated"] = df["family"].apply(_choose_dominant)
    
    # Add Rho info (basic recombination rate)
    df = add_rho(df)
    
    # Add TE density
    window_size = 50_000
    df = add_TEdensity(df, TES_GFF, GENOME_FILE, OUT_DIR, window_size)
    # Add gene density
    if SPECIES == 'PCOM':
        df = add_gene_density(df, GENES_GFF_CLEAN, GENOME_FILE, OUT_DIR, window_size)
    else:
        df = add_gene_density(df, GENES_GFF, GENOME_FILE, OUT_DIR, window_size)

    # Fetch GO ID into level 1
    df = add_GO_level0( df )
    
    # Reduce GO ID into level 2
    df = add_GO_level2( df )
    
    # Highligh 'NBS' locus (related to immunity)
    # df = is_NBS( df )
    
    # Highligh 'disease' locus based on 'Ath_description'
    #df = is_Disease( df )
    
    # Highligh immunity function
    df = flag_immunity( df )
    
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"✅ Summary written: {OUT_TSV}")
    print(f"Total variants: {len(df)}")
    
    # to check there is no redundant lines
    assert df.shape[0] == df[["POS", "ID"]].drop_duplicates().shape[0]

# ================================================================
# 5️⃣ ENTRY POINT
# ================================================================

if __name__ == "__main__":
    main()
