#!/usr/bin/env python3

"""
REPET2RM_fasta.py

Convert a REPET consensus FASTA and a classification table
into a RepeatMasker-style FASTA library.

Final header format:
    >FamilyName\tClass/Subclass\tSource

- FamilyName comes from the consensus FASTA header (e.g. '>name')
- Class/Subclass comes from the classification table
- Source is "REPET"
"""

import argparse
import csv
import gzip
import io
import re
from typing import Dict, Tuple, Optional, List

# Canonical top-level classes
CANON_CLASSES = [
    "LTR", "LINE", "SINE", "DNA", "RC", "Satellite", "Simple_repeat", "Unknown"
]

def open_maybe_gzip(path: str, mode: str = "rt"):
    """Open plain text or gzipped file automatically."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def normalize_name(raw: str) -> str:
    """Extract sequence name from FASTA header (remove '>', cut at first space)."""
    name = raw.strip()
    if name.startswith(">"):
        name = name[1:]
    # only keep first token before any whitespace
    return name.split()[0]

def reduce_to_two_levels(cls_str: str) -> Tuple[str, str]:
    """Reduce a classification string to 'Class/Subclass'."""
    if not cls_str:
        return ("Unknown", "Unknown")

    # Normalize separators to "/"
    s = re.sub(r"[|]+", "/", cls_str.strip())
    s = re.sub(r"\s+", "/", s)
    tokens_raw = [t for t in s.split("/") if t]
    tokens_low = [t.lower() for t in tokens_raw]

    # Map synonyms
    mapped_tokens: List[str] = []
    for t in tokens_low:
        if t in ("classi", "retrotransposon", "non-ltr"):
            continue
        if t in ("classii", "tir", "mite", "helitron", "hat", "tc1-mariner", "tc-mar", "tc-1", "pif-harb"):
            mapped_tokens.append("dna")
        else:
            mapped_tokens.append(t)

    # Find main class
    cls_idx = -1
    cls_out = "Unknown"
    for i, t in enumerate(mapped_tokens):
        if t in ("ltr", "line", "sine", "dna", "rc", "satellite", "simple_repeat", "simple", "unknown"):
            cls_out = {
                "ltr": "LTR", "line": "LINE", "sine": "SINE", "dna": "DNA",
                "rc": "RC", "satellite": "Satellite",
                "simple_repeat": "Simple_repeat", "simple": "Simple_repeat",
                "unknown": "Unknown",
            }[t]
            cls_idx = i
            break

    # Pick subclass (token after class)
    sub_out = "Unknown"
    if cls_idx >= 0:
        for j in range(cls_idx + 1, len(tokens_raw)):
            cand = tokens_raw[j]
            if cand and cand.lower() not in ("unknown", "other", "na"):
                sub_out = cand.replace(" ", "_")
                break

        # Special case for DNA subclasses
        if cls_out == "DNA":
            dna_map = {
                "tc1-mariner": "TcMar", "tc-mar": "TcMar", "tc1": "TcMar",
                "hat": "hAT", "hat-ac": "hAT-Ac",
                "pif-harb": "PIF-Harb", "tir": "TIR", "mite": "MITE",
            }
            key = sub_out.lower()
            if key in dna_map:
                sub_out = dna_map[key]

    return (cls_out, sub_out)

def read_classification_table(path: str) -> Dict[str, str]:
    """
    Read classification table and return mapping:
        consensus_name -> "Class/Subclass"
    Supports TSV/CSV/whitespace with or without header.
    """
    with open_maybe_gzip(path, "rt") as fh:
        text = fh.read()

    # Clean newlines and skip comments
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = [ln for ln in text.split("\n") if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return {}

    # Try delimiter detection
    sample = "\n".join(lines[:20])
    dialect = None
    for delim in ["\t", ",", ";"]:
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=delim)
            break
        except Exception:
            continue

    rows: List[List[str]] = []
    headers: Optional[List[str]] = None

    if dialect:
        reader = csv.reader(io.StringIO("\n".join(lines)), dialect)
        first = next(reader)
        # Detect header
        header_like = any(k.lower() in (
            "name","family","seqname","model","element","consensus",
            "te","id","classification","classif","class","subclass",
            "superfamily","order"
        ) for k in first)
        if header_like:
            headers = [h.strip() for h in first]
        else:
            rows.append(first)
        for row in reader:
            if row:
                rows.append([c.strip() for c in row])
    else:
        for idx, ln in enumerate(lines):
            parts = re.split(r"\s+", ln.strip())
            if idx == 0 and any(k.lower() in (
                "name","family","seqname","model","element","consensus","te",
                "id","classification","classif","class","subclass","superfamily","order"
            ) for k in parts):
                headers = parts
            else:
                rows.append(parts)

    def col_index(names: List[str], preferred: List[str]) -> Optional[int]:
        low = [n.lower() for n in names]
        for pref in preferred:
            if pref in low:
                return low.index(pref)
        return None

    name_idx = 0
    class_idx = None
    subclass_idx = None
    combined_idx = None

    if headers:
        name_idx = col_index(headers, ["name","family","seqname","model","element","consensus","te","id"]) or 0
        subclass_idx = col_index(headers, ["subclass","superfamily","order"])
        class_idx = col_index(headers, ["class","te_class"])
        combined_idx = col_index(headers, ["classification","classif","taxonomy","annot","assignation"])

    mapping: Dict[str, str] = {}
    for row in rows:
        if not row:
            continue
        nm = normalize_name(row[name_idx])
        cls, sub = "Unknown", "Unknown"
        if combined_idx is not None and combined_idx < len(row):
            cls, sub = reduce_to_two_levels(row[combined_idx])
        else:
            cls_str = row[class_idx] if (class_idx is not None and class_idx < len(row)) else ""
            sub_str = row[subclass_idx] if (subclass_idx is not None and subclass_idx < len(row)) else ""
            if cls_str and not sub_str:
                cls, sub = reduce_to_two_levels(cls_str)
            elif sub_str and not cls_str:
                cls, sub = reduce_to_two_levels(sub_str)
            elif cls_str or sub_str:
                cls = cls_str.strip() if cls_str else "Unknown"
                sub = sub_str.strip() if sub_str else "Unknown"

        if cls not in CANON_CLASSES:
            cls, sub = reduce_to_two_levels(f"{cls}/{sub}")

        mapping[nm] = f"{cls}/{sub}"

    return mapping

def write_rm_fasta(consensus_fa: str, class_map: Dict[str, str], out_fa: str) -> Tuple[int, int]:
    """
    Rewrite FASTA headers into RepeatMasker format:
        >Name\tClass/Subclass\tREPET
    """
    n_total = 0
    n_mapped = 0
    with open_maybe_gzip(consensus_fa, "rt") as fin, open(out_fa, "wt") as fout:
        current_name = None
        seq_buf = []

        def flush():
            nonlocal n_total, n_mapped, current_name, seq_buf
            if current_name is None:
                return
            n_total += 1
            nm = normalize_name(current_name)
            cls = class_map.get(nm, "Unknown/Unknown")
            if nm in class_map:
                n_mapped += 1
            fout.write(f">{nm}\t{cls}\tREPET\n")
            if seq_buf:
                fout.write("".join(seq_buf) + "\n")
            current_name = None
            seq_buf = []

        for line in fin:
            if line.startswith(">"):
                flush()
                current_name = line.strip()
            else:
                if current_name:
                    seq_buf.append(line.strip())
        flush()

    return n_total, n_mapped

def main():
    ap = argparse.ArgumentParser(description="Convert REPET consensus + classification to RepeatMasker-style FASTA headers.")
    ap.add_argument("--consensus", "-c", required=True, help="REPET consensus FASTA (can be .gz)")
    ap.add_argument("--classification", "-t", required=True, help="Classification table (TSV/CSV/whitespace; can be .gz)")
    ap.add_argument("--out", "-o", required=True, help="Output FASTA with RepeatMasker-style headers")
    args = ap.parse_args()

    class_map = read_classification_table(args.classification)
    n_total, n_mapped = write_rm_fasta(args.consensus, class_map, args.out)
    print(f"Wrote {args.out}")
    print(f"Total records: {n_total}")
    print(f"Classified records: {n_mapped} ({(n_mapped/n_total*100.0 if n_total else 0):.1f}%)")

if __name__ == "__main__":
    main()
