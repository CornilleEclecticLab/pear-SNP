#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import argparse

def read_lengths(path):
    """
    Read sequence lengths from either:
      - .fai (samtools index) or
      - TSV: seq<TAB>length
    Returns: dict {seq: length}
    """
    lens = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[1].isdigit():
                lens[cols[0]] = int(cols[1])
    return lens

def read_classif(path):
    """
    Read REPET classification table (TSV).
    Expected columns include at least:
      Seq_name, length, class, order, Wcode, sFamily
    Returns: dict {Seq_name: row_dict}
    """
    data = {}
    with open(path, "r", encoding="utf-8") as fh:
        header = None
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if i == 0:
                header = cols
                continue
            row = {header[j]: (cols[j] if j < len(cols) else "") for j in range(len(header))}
            key = row.get("Seq_name")
            if key:
                data[key] = row
    return data

def parse_attrs(s):
    """
    Parse GFF attribute string into a dict.
    Example:
      ID=...;Target=NAME 2220 2322;TargetLength=2542;AlignIdentity=78.57
    """
    out = {}
    for part in s.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v.strip()
        else:
            out[part] = ""
    return out

def parse_target(val):
    """
    Parse 'Target=NAME start end' attribute.
    Returns: (name, begin, end) with begin <= end when possible.
    """
    if not val:
        return None, None, None
    toks = val.split()
    if len(toks) < 3:
        return toks[0], None, None
    name = toks[0]
    try:
        tbeg = int(toks[1])
        tend = int(toks[2])
    except Exception:
        return name, None, None
    if tbeg is not None and tend is not None and tbeg > tend:
        tbeg, tend = tend, tbeg
    return name, tbeg, tend

# Map Wicker codes to RepeatMasker-like class/family.
WCODE_MAP = {
    # LTR retrotransposons
    "RLG": "LTR/Gypsy",
    "RLC": "LTR/Copia",
    "RLB": "LTR/Bel-Pao",
    "RLE": "LTR/ERV",
    "RLX": "LTR/Unknown",
    # LINEs
    "RIL": "LINE/L1",
    "RIR": "LINE/RTE",
    "RIX": "LINE/Unknown",
    # SINEs
    "RSX": "SINE/Unknown",
    # DNA TIR/others
    "DTA": "DNA/hAT",
    "DTC": "DNA/CACTA",
    "DTM": "DNA/MuDR",
    "DTT": "DNA/Tc1-Mariner",
    "DTP": "DNA/PIF-Harbinger",
    "DTH": "DNA/Helitron",
    "DTX": "DNA/Unknown",
    "DHX": "DNA/Unknown",
}

def norm_na(x):
    """Normalize NA-like strings to empty."""
    if x is None:
        return ""
    x = x.strip()
    return "" if x in {"NA", "na", ".", ""} else x

def classfam_from_classif(row):
    """
    Build 'class/family' using the classification row.
    Priority:
      1) sFamily + order/class
      2) Wcode
      3) order/class fallbacks
    """
    if not row:
        return None

    sfamily = norm_na(row.get("sFamily"))
    order = norm_na(row.get("order"))
    rclass = norm_na(row.get("class"))  # e.g., "I", "II", "Unclassified"
    wcode = norm_na(row.get("Wcode"))

    # 1) If we have a concrete subfamily name, prefer it.
    if sfamily:
        if order in {"LTR", "LINE", "SINE"}:
            return f"{order}/{sfamily}"
        # Otherwise map by class: II => DNA
        top = "DNA" if rclass in {"II", "ClassII", "2"} else ("LTR" if order == "LTR" else "")
        if top:
            return f"{top}/{sfamily}"
        # Last resort with sFamily
        return f"Unspecified/{sfamily}"

    # 2) Wcode mapping
    if wcode and wcode in WCODE_MAP:
        return WCODE_MAP[wcode]

    # 3) Fallbacks from order/class
    if order in {"LTR", "LINE", "SINE"}:
        return f"{order}/Unknown"
    if rclass in {"II", "ClassII", "2"}:
        return "DNA/Unknown"

    return "Unspecified"

CLASS_RE = re.compile(r'(ClassI{1,2}):([^:]+):([^:|()\s]+)')

def classfam_from_description(desc):
    """
    Infer 'class/family' from TargetDescription when classification is missing.
    Examples:
      'ClassII:TIR:Tc1-Mariner:?: ...' -> 'DNA/Tc1-Mariner'
      'ClassI:LTR:Gypsy:?: ...'        -> 'LTR/Gypsy'
    """
    if not desc:
        return None
    m = CLASS_RE.search(desc)
    if m:
        cl, group, fam = m.groups()
        if cl.startswith("ClassII"):
            return f"DNA/{fam}"
        if group in ("LTR", "LINE", "SINE"):
            return f"{group}/{fam}"
        return f"DNA/{fam}"
    m2 = re.search(r'Wcode:([A-Z]{3})', desc)
    if m2 and m2.group(1) in WCODE_MAP:
        return WCODE_MAP[m2.group(1)]
    return None

# --- NEW: normalize class/family for MEGAnE ---
def norm_classfam(cf: str) -> str:
    """
    Normalize to allowed top classes and ensure 'Class/Subclass' with safe tokens.
    """
    if not cf or cf.lower().startswith("unspecified"):
        return "Unknown/Unknown"
    cf = cf.strip().replace(" ", "_").replace("|", "_")
    if "/" not in cf:
        cf = f"{cf}/Unknown"
    top, sub = cf.split("/", 1)
    allowed = {"LTR","LINE","SINE","DNA","RC","Satellite","Simple_repeat","Unknown"}
    if top not in allowed:
        # Map common variants
        m = {"Simple":"Simple_repeat","simple":"Simple_repeat","simple_repeat":"Simple_repeat",
             "unknown":"Unknown","classii":"DNA","classi":"LTR"}
        top2 = m.get(top, top)
        if top2 not in allowed:
            top2 = "Unknown"
        top = top2
    sub = sub if sub else "Unknown"
    return f"{top}/{sub}"

def main():
    parser = argparse.ArgumentParser(
        description="Convert REPET TEannot GFF to RepeatMasker .out (approximate)."
    )
    parser.add_argument("gff", help="REPET GFF file")
    parser.add_argument("-l", "--lengths", help="Sequence lengths (.fai or TSV seq\\tlen)", default=None)
    parser.add_argument("-c", "--classif", help="REPET classification TSV", default=None)
    parser.add_argument("-o", "--out", help="Output .out (default: stdout)", default="-")
    args = parser.parse_args()

    # Optional: genome/contig lengths for '(left)' in query
    lens = read_lengths(args.lengths) if args.lengths else {}

    # Optional: classification table to refine class/family and consensus length
    classif = read_classif(args.classif) if args.classif else {}

    outfh = sys.stdout if args.out == "-" else open(args.out, "w", encoding="utf-8")

    # RepeatMasker-like header
    header = (
        "    SW  perc perc perc  query                      position in query        matching                                repeat          position in repeat\n"
        " score  div. del. ins.  sequence                   begin   end    (left)   repeat                                  class/family    begin   end    (left)   ID\n"
    )
    print(header, file=outfh)

    rid = 0
    with open(args.gff, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs_s = cols
            if ftype != "match":
                # Emit one .out line per 'match' only; 'match_part' are internal fragments
                continue

            # Query coordinates (1-based inclusive)
            try:
                qbeg = int(start)
                qend = int(end)
                if qbeg > qend:
                    qbeg, qend = qend, qbeg
            except Exception:
                continue

            attrs = parse_attrs(attrs_s)
            tname, tbeg, tend = parse_target(attrs.get("Target", ""))
            target_desc = attrs.get("TargetDescription", "")
            # Prefer TargetLength from GFF; else take 'length' from classification
            tlen = None
            if attrs.get("TargetLength"):
                try:
                    tlen = int(attrs["TargetLength"])
                except Exception:
                    tlen = None
            if tlen is None and tname and tname in classif:
                try:
                    tlen = int(classif[tname].get("length")) if classif[tname].get("length") else None
                except Exception:
                    tlen = None
            # Sometimes names differ by a trailing "_reversed" in classification
            if tlen is None and tname and (tname + "_reversed") in classif:
                try:
                    tlen = int(classif[tname + "_reversed"].get("length"))
                except Exception:
                    tlen = None

            # Identity → approximate divergence
            ident = None
            if attrs.get("AlignIdentity"):
                try:
                    ident = float(attrs["AlignIdentity"])
                except Exception:
                    ident = None
            elif attrs.get("Identity"):
                try:
                    ident = float(attrs["Identity"])
                except Exception:
                    ident = None
            div = (100.0 - ident) if ident is not None else 0.0

            # Unknown deletions/insertions → 0.0
            perc_del = 0.0
            perc_ins = 0.0

            # Smith–Waterman score is not provided by REPET → 0 (placeholder)
            sw = 0

            # RepeatMasker strand notation: '+' or 'C' (complement)
            rm_strand = '+' if strand == '+' else 'C'

            # Determine class/family (with normalization)
            classfam = None
            if tname and tname in classif:
                classfam = classfam_from_classif(classif[tname])
            elif tname and (tname + "_reversed") in classif:
                classfam = classfam_from_classif(classif[tname + "_reversed"])
            if not classfam:
                classfam = classfam_from_description(target_desc) or "Unknown/Unknown"
            classfam = norm_classfam(classfam)

            # Query '(left)': needs query length if provided
            left_q = "."
            qlen = lens.get(seqid)
            if qlen is not None:
                if rm_strand == '+':
                    left_q = str(max(0, qlen - qend))
                else:
                    left_q = str(max(0, qbeg - 1))

            # Repeat positions and '(left)' within the consensus
            rb, re_ = (tbeg or 0), (tend or 0)
            if rb and re_ and rb > re_:
                rb, re_ = re_, rb
            if tlen is not None and re_:
                if rm_strand == '+':
                    left_r = str(max(0, tlen - re_))
                else:
                    left_r = str(max(0, rb - 1))
            else:
                left_r = "."

            rid += 1

            # -------- FIXED OUTPUT FORMATTING (explicit spaces between fields) --------
            out_line = (
                f"{sw:7d}"
                f" {div:6.1f} {perc_del:5.1f} {perc_ins:5.1f}  "
                f"{seqid:<27s}"
                f" {qbeg:9d} {qend:9d}  ({left_q:>5})  "
                f"{rm_strand:1s}  "
                f"{(tname or 'NA'):<41s}"
                f" {classfam:<16s}"
                f" {rb:7d} {re_:8d}  ({left_r:>5})"
                f" {rid:8d}"
            )
            print(out_line, file=outfh)

    if outfh is not sys.stdout:
        outfh.close()

if __name__ == "__main__":
    main()


