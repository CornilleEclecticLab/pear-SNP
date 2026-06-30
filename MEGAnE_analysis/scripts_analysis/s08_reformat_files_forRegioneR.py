#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline helper for regioneR inputs.

Tasks:
 1. TSV (positive selection summary) -> BED (genes of interest, 0-based half-open)
 2. FASTA -> chrom.sizes (contig sizes file)
 3. MEGAnE summary table -> BED files per specific population (10 kb windows)
 4. GFF genome annotation -> BED (all genes, the "universe")
 5. BED (genes of interest) + Strand from TSV + chrom.sizes -> BED (2kbp upstream)
"""

import os
import re
import sys
import pandas as pd
from typing import Dict, Iterable, Optional, Tuple

# ============================================================
# Utils
# ============================================================
def parse_chrom_sizes(sizes_path: str) -> Dict[str, int]:
    """Read a chrom.sizes file (2 columns: name, length) into a dict."""
    sizes: Dict[str, int] = {}
    with open(sizes_path) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, L = line.rstrip("\n").split("\t")[:2]
            sizes[chrom] = int(L)
    return sizes


# ============================================================
# 1) TSV -> BED (genes of interest)  ——  1-based (TSV)  ->  0-based (BED)
# ============================================================
def tsv_to_bed(tsv_in: str, bed_out: str) -> None:
    """
    Convert a TSV with columns: Chr_ID, Start, End, Gene_ID, Strand
    into a BED (0-based half-open) with 4 columns:
        chr   start0   end0   Gene_ID

    TSV is expected to be 1-based inclusive (gene tables usual convention).
    We convert: start0 = Start - 1 ; end0 = End.
    """
    df = pd.read_csv(tsv_in, sep="\t")

    needed = ["Chr_ID", "Start", "End", "Gene_ID", "Strand"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {tsv_in}: {missing}")

    # Convert to true BED (0-based half-open)
    bed_df = df[["Chr_ID", "Start", "End", "Gene_ID"]].copy()
    bed_df["Start"] = pd.to_numeric(bed_df["Start"], errors="raise").astype(int) - 1
    bed_df["End"]   = pd.to_numeric(bed_df["End"],   errors="raise").astype(int)

    bed_df.to_csv(bed_out, sep="\t", index=False, header=False)
    print(f"[OK] BED (genes of interest; 0-based): {bed_out} ({len(bed_df)} rows)")


# ============================================================
# 2) FASTA -> chrom.sizes
# ============================================================
def fasta_sizes(fasta_path: str):
    """Return a list of (seqname, length) extracted directly from FASTA."""
    sizes = []
    with open(fasta_path, "r") as f:
        name = None
        length = 0
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    sizes.append((name, length))
                name = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            sizes.append((name, length))
    return sizes


def write_sizes(fasta_in: str, sizes_out: str) -> None:
    """Write a chrom.sizes file: 2 columns [name, length]."""
    sizes = fasta_sizes(fasta_in)
    with open(sizes_out, "w") as out:
        for n, L in sizes:
            out.write(f"{n}\t{L}\n")
    print(f"[OK] chrom.sizes: {sizes_out} ({len(sizes)} sequences)")


# ============================================================
# 3) MEGAnE summary -> BED per population
# ============================================================
def clean_name(pop: str) -> str:
    """Sanitize a population name for filenames."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(pop))


def unique_population(pop_str: str):
    """
    Return the population if all listed populations in the string are identical.
    Example: "Dessert" -> "Dessert"; "Dessert,Perry" -> None.
    """
    if pd.isna(pop_str) or str(pop_str).strip() == "":
        return None
    parts = [p.strip() for p in str(pop_str).split(",") if p.strip() != ""]
    if not parts:
        return None
    return parts[0] if len(set(parts)) == 1 else None


def megane_to_beds(mega_in: str) -> int:
    """
    Read MEGAnE_summary_PCOM.tsv and create one BED per specific population,
    plus one global BED file combining all populations.

    Input columns: CHROM, POS, END, ID, populations
    Output: 
        <mega_in>.<population>.bed  → one per specific population
        <mega_in>.ALL.bed           → all populations combined
    Columns in each BED:
        chr   start   end   ID
    """
    df = pd.read_csv(mega_in, sep="\t")

    needed = ["CHROM", "POS", "END", "ID", "populations"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {mega_in}: {missing}")

    df["unique_pop"] = df["populations"].apply(unique_population)
    spec_df = df.dropna(subset=["unique_pop"]).copy()

    #chr
    spec_df["chr"] = spec_df["CHROM"].astype(str)
    
    #start
    spec_df["start"] = pd.to_numeric(spec_df["POS"], errors="coerce").astype("Int64")
    spec_df = spec_df.dropna(subset=["start"]).copy()
    spec_df["start"] = spec_df["start"].astype(int)
    
    #end
    spec_df["end"] = pd.to_numeric(spec_df["END"], errors="coerce").astype("Int64")
    spec_df = spec_df.dropna(subset=["end"]).copy()
    spec_df["end"] = spec_df["end"].astype(int)
    
    #ID
    spec_df["name"] = spec_df["ID"].astype(str)

    created = 0
    # ---- One BED per specific population ----
    for pop, g in spec_df.groupby("unique_pop", dropna=True):
        suffix = clean_name(pop)
        out_path = f"{mega_in}.{suffix}.bed"
        g[["chr", "start", "end", "name"]].to_csv(out_path, sep="\t", index=False, header=False)
        print(f"[OK] BED (population '{pop}'): {out_path} ({len(g)} loci)")
        created += 1
        
    # ---- Add a global BED file with all populations combined ----
    all_bed_path = f"{mega_in}.ALL.bed"
    all_df = spec_df[["chr", "start", "end", "name"]].drop_duplicates()
    all_df.to_csv(all_bed_path, sep="\t", index=False, header=False)
    print(f"[OK] Global BED (all populations): {all_bed_path} ({len(all_df)} loci)")
    created += 1

    if created == 1:
        print("[INFO] Only global BED file created (no specific populations found).")
    else:
        print(f"[OK] {created} BED file(s) created in total.")
    return created


# ============================================================
# 4) GFF -> BED (universe of all genes)
# ============================================================
def gff_to_bed(gff_in: str, bed_out: str) -> None:
    """
    Convert a GFF3 genome annotation to BED (all genes).
    Keep only 'gene' features, convert coordinates to 0-based half-open.
    BED columns: chr   start   end   gene_id
    """
    with open(gff_in) as fin, open(bed_out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = fields
            if feature != "gene":
                continue

            bed_start = int(start) - 1  # GFF is 1-based inclusive
            bed_end = int(end)          # BED is 0-based half-open

            gene_id = "NA"
            for attr in attrs.split(";"):
                if attr.startswith("ID=gene:"):
                    gene_id = attr.replace("ID=gene:", "")
                    break
                elif attr.startswith("ID="):
                    gene_id = attr.replace("ID=", "")
                    break
                elif attr.startswith("Name="):
                    gene_id = attr.replace("Name=", "")
                    break

            fout.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_id}\n")
    print(f"[OK] BED (all genes): {bed_out}")


# ============================================================
# 5) BED + Strand (from TSV) + chrom.sizes -> BED (2kbp upstream)
# ============================================================
def load_strand_map_from_tsv(tsv_in: str) -> Dict[str, str]:
    """Return {Gene_ID -> Strand} from TSV."""
    df = pd.read_csv(tsv_in, sep="\t", usecols=["Gene_ID", "Strand"])
    # In case of duplicates, check consistency; keep the first and warn if conflicting
    d: Dict[str, str] = {}
    for gid, s in df.itertuples(index=False, name=None):
        s = str(s).strip()
        if s not in {"+", "-"}:
            continue
        if gid in d and d[gid] != s:
            print(f"[WARN] Conflicting strand for {gid}: '{d[gid]}' vs '{s}'. Keeping first.", file=sys.stderr)
            continue
        d[str(gid)] = s
    return d


def upstream_2kb_from_bed(bed_in: str, tsv_in: str, sizes_path: str, out_bed: str, upstream: int = 2000) -> None:
    """
    Build a BED with 2kbp upstream of each gene using:
      - positions from BED (already 0-based half-open),
      - strand from TSV,
      - chromosome lengths from chrom.sizes.

    Rules (ensure one boundary equals the gene boundary):
      '+' strand: [start0 - 2000, start0)
      '-' strand: [end0, end0 + 2000)
    """
    chrom_sizes = parse_chrom_sizes(sizes_path)
    strand_map = load_strand_map_from_tsv(tsv_in)

    bed = pd.read_csv(bed_in, sep="\t", header=None, names=["chr", "start0", "end0", "Gene_ID"])
    out_rows = []
    skipped = 0

    for _, r in bed.iterrows():
        chrom = str(r["chr"])
        start0 = int(r["start0"])
        end0 = int(r["end0"])
        gid = str(r["Gene_ID"])

        if chrom not in chrom_sizes:
            skipped += 1
            print(f"[WARN] Chrom '{chrom}' not in sizes; skipping {gid}.", file=sys.stderr)
            continue

        strand = strand_map.get(gid)
        if strand not in {"+", "-"}:
            skipped += 1
            print(f"[WARN] Strand missing/invalid for {gid}; skipping.", file=sys.stderr)
            continue

        chr_len = chrom_sizes[chrom]

        if strand == "+":
            up_start = max(0, start0 - upstream)
            up_end   = end0
        else:  # strand == "-"
            up_start = end0
            up_end   = min(chr_len, end0 + upstream)

        if up_start < up_end:
            out_rows.append((chrom, up_start, up_end, gid))
        else:
            skipped += 1
            print(f"[WARN] Invalid upstream for {gid} ({chrom}:{up_start}-{up_end}); skipping.", file=sys.stderr)

    if out_rows:
        with open(out_bed, "w") as fout:
            for chrom, s, e, name in out_rows:
                fout.write(f"{chrom}\t{s}\t{e}\t{name}\n")
        print(f"[OK] BED (2kbp upstream): {out_bed} ({len(out_rows)} rows, skipped={skipped})")
    else:
        print("[INFO] No upstream intervals produced.")


# ============================================================
# Main
# ============================================================
def main():
    
    # 1) TSV with genes of interest -> BED (make true 0-based)
    tsv_in = "../data/PCOM/genes/positive_selection_summary_table.West.tsv"
    bed_out = tsv_in + ".bed"
    tsv_to_bed(tsv_in, bed_out)

    # 2) FASTA -> chrom.sizes
    fasta_in = "../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat"
    sizes_out = fasta_in + ".sizes"
    write_sizes(fasta_in, sizes_out)

    # 3) MEGAnE summary -> BED per population
    mega_in = "../output/PCOM/MEGAnE_summary_PCOM.tsv"
    megane_to_beds(mega_in)

    # 4) GFF -> BED (universe of all genes)
    gff_in = "../data/PCOM/genes/PyrusCommunis_BartlettDHv2.0.gff"
    gff_out = gff_in + ".bed"
    gff_to_bed(gff_in, gff_out)

    # 5) from BED (0-based) + Strand (TSV) + chrom.sizes -> 2kb upstream BED
    upstream_bed_out = tsv_in + ".2kbp_upstream.bed"
    upstream_2kb_from_bed(bed_out, tsv_in, sizes_out, upstream_bed_out, upstream=2000)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
