#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MEGAnE — Figures runners:
- run_afs_chi2_summary(): AFS on top, χ² residuals below, per population
- run_chi2_tables(): write single TSV with all contingency tables + TSV with χ² stats
- run_me_classorder_heatmap(): single Class and Order heatmap (Fixed vs Polymorphic)

"""

import os, gzip, warnings, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import gridspec
from scipy.stats import chi2_contingency
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import patheffects as pe
import gzip
import subprocess
from Bio import Phylo
from PIL import Image, ImageDraw, ImageFont
from adjustText import adjust_text
import allel
import matplotlib.lines as mlines


warnings.filterwarnings("ignore", category=UserWarning)
#sns.set_style("whitegrid")


# ============================== Loaders ==================================

def load_classification(classif_file):
    df = pd.read_csv(classif_file, sep="\t", dtype=str)
    for col in ["Seq_name", "class", "order"]:
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {classif_file}")
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
        df.loc[df[c].isin(["", "nan", "None", "NaN"]), c] = np.nan
    df["class"]  = df["class"].fillna("NA")
    df["order"]  = df["order"].fillna("NA")
    if "sFamily" in df.columns:
        df["sFamily"] = df["sFamily"].fillna("NA")
    else:
        df["sFamily"] = "NA"
    class_dict   = dict(zip(df["Seq_name"], df["class"]))
    order_dict   = dict(zip(df["Seq_name"], df["order"]))
    sfamily_dict = dict(zip(df["Seq_name"], df["sFamily"]))
    return df, class_dict, order_dict, sfamily_dict

def load_population_table(pop_file):
    df = pd.read_csv(pop_file, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    if "#ID_vcf" in df.columns:
        df = df.rename(columns={"#ID_vcf": "ID_vcf"})
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
        df.loc[df[c].isin(["", "nan", "None", "NaN"]), c] = np.nan
    if "Population" not in df.columns:
        raise ValueError("Missing 'Population' column")
    df["Population"] = df["Population"].fillna("Unknown")
    return df.set_index("ID_vcf")

def _load_all_inputs( VCF_MEI, VCF_MEA, POP_FILE, CLASSIF_FILE ):
    _, class_dict, order_dict, sfamily_dict = load_classification(CLASSIF_FILE)
    df_pop = load_population_table(POP_FILE)
    G_mei, var_ids_mei = parse_vcf_genotypes_with_ids(VCF_MEI, class_dict, order_dict, sfamily_dict)
    G_mea, var_ids_mea = parse_vcf_genotypes_with_ids(VCF_MEA, class_dict, order_dict, sfamily_dict)
    return class_dict, order_dict, df_pop, G_mei, var_ids_mei, G_mea, var_ids_mea, sfamily_dict

# ============================== VCF parsing ==============================

def parse_vcf_genotypes_with_ids(vcf_file, class_dict, order_dict, sfamily_dict):
    """Return G (samples×variants) and var_ids (list of ME IDs per variant)."""
    samples, variant_genos, var_ids = [], [], []
    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                hdr = line.strip().split("\t")
                samples = hdr[9:]; continue
            cols = line.strip().split("\t")
            if len(cols) < 10: continue
            info = cols[7]; fmt = cols[8]; genotypes = cols[9:]
            mei_ids = []
            for kv in info.split(";"):
                if kv.startswith("MEI="):
                    mei_ids = kv.split("=")[1].split("|"); break
            var_ids.append(mei_ids if mei_ids else [])
            ff = fmt.split(":")
            if "GT" not in ff: continue
            gt_idx = ff.index("GT")
            row = []
            for g in genotypes:
                parts = g.split(":")
                if gt_idx >= len(parts): row.append(np.nan); continue
                gt = parts[gt_idx]
                if gt in (".","./.",".|."): row.append(np.nan); continue
            
                al = gt.replace("|","/").split("/") # normalize separators
                vals = [int(a) for a in al if a!="."]
                row.append(sum(vals) if len(vals)>=2 else np.nan) # version 1 where I count 2 for homozygous insertions (1/1)
                
                #### TO RM 
                #row.append(1 if "1" in al else 0) version where I count just 1 for homozygous insertions (1/1)

                # #version 3: handle partial missing genotypes, 1/. or 0/. are misising
                # if al == ["1", "."] or al == [".", "1"]:
                #     row.append(np.nan)
                #     continue
                # elif al == ["0", "."] or al == [".", "0"]:
                #     row.append(np.nan)
                #     continue
                # # standard complete cases
                # vals = [int(a) for a in al if a != "."]
                # row.append(sum(vals) if len(vals) == 2 else np.nan)
                #### TO RM 

            variant_genos.append(row)
    idx = [s.split(".")[2] if len(s.split("."))>=3 else s for s in samples]
    if not variant_genos:
        return pd.DataFrame(index=idx), var_ids
    G = pd.DataFrame(variant_genos).T
    G.columns = [f"var{i}" for i in range(G.shape[1])]
    G.index = idx
    return G, var_ids

def compute_maf_and_missing(G_mei, G_mea, maf_threshold=0.05, missing_threshold=0.10):
    """
    Filter MEI and MEA genotype matrices based on MAF and missing rate thresholds.

    Parameters
    ----------
    G_mei, G_mea : pandas.DataFrame
        Genotype matrices (samples × variants) with values 0, 1, 2, or NaN.
    maf_threshold : float
        Minimum minor allele frequency to retain a locus.
    missing_threshold : float
        Maximum allowed missing rate per locus.

    Returns
    -------
    G_mei_filt, G_mea_filt : pandas.DataFrame
        Filtered genotype matrices.
    """

    def maf_and_missing(G):
        # Compute allele frequency per locus
        allele_sums = G.sum(axis=0, skipna=True)
        allele_counts = 2 * G.notna().sum(axis=0)
        af = allele_sums / allele_counts
        maf = af.clip(upper=1 - af)  # minor allele frequency
        missing_rate = 1 - (G.notna().sum(axis=0) / G.shape[0])
        return maf, missing_rate

    # Compute metrics
    maf_mei, miss_mei = maf_and_missing(G_mei)
    maf_mea, miss_mea = maf_and_missing(G_mea)

    # Apply filtering
    keep_mei = (maf_mei > maf_threshold) & (miss_mei < missing_threshold)
    keep_mea = (maf_mea > maf_threshold) & (miss_mea < missing_threshold)

    # Filter genotype matrices
    G_mei_filt = G_mei.loc[:, keep_mei]
    G_mea_filt = G_mea.loc[:, keep_mea]

    # Print summary
    print(f"Filtered MEI: kept {keep_mei.sum()} / {len(keep_mei)} loci "
          f"(MAF > {maf_threshold}, missing < {missing_threshold})")
    print(f"Filtered MEA: kept {keep_mea.sum()} / {len(keep_mea)} loci "
          f"(MAF > {maf_threshold}, missing < {missing_threshold})")

    return G_mei_filt, G_mea_filt


# ===================== Presence / populations helpers ====================

def _dosage_te(G_values, type_name):
    """MEI: dosage = GT;  MEA: dosage = 2 - GT (MEA GT=0 means TE presence)."""
    g = np.array(G_values, dtype=float)
    return g if type_name == "MEI" else (2.0 - g)

def variant_population_sets(G, type_name, df_pop):
    """For each variant, set of populations with ≥1 sample carrying TE (dosage > 0)."""
    if G is None or G.empty:
        return []
    d = _dosage_te(G.values, type_name)
    present = d > 0
    pops = df_pop["Population"].reindex(G.index.astype(str)).fillna("Unknown").values
    out = []
    for j in range(G.shape[1]):
        out.append(set(pops[present[:, j]]) if present[:, j].any() else set())
    return out

def unique_sample_count_in_pop(pop, G_mei, G_mea, df_pop):
    """Unique sample count in 'pop' present in MEI and/or MEA matrices."""
    pop_ids = set(df_pop.index[df_pop["Population"] == pop].astype(str))
    mei_ids = set(G_mei.index.astype(str)) if not G_mei.empty else set()
    mea_ids = set(G_mea.index.astype(str)) if not G_mea.empty else set()
    return len((mei_ids | mea_ids) & pop_ids)

# ==================== Canonical class/order mapping =======================

CLASS_I_ORDERS  = {"LTR","LINE","SINE","DIRS","PLE","TRIM","LARD","GYPSY","COPIA","ERV","BEL"}
CLASS_II_ORDERS = {"TIR","HELITRON","MITE","HARBINGER","HAT","TC1","MARINER","TC1/MARINER","MERLIN",
                   "CHAPAEV","TRANSIB","MUDR","PIGGYBAC","SOLA","ZISUPTON","UCON","EULOR","CHOMPY",
                   "REP","HOPPER","MER"}
UNCLASSIFIED_TOKENS = {"UNCLASSIFIED","UNKNOWN","OTHER","OTHERS"}

def canonicalize_order_label(order_name: str) -> str:
    if not order_name:
        return "NA"
    o = str(order_name).strip(); u = o.upper()
    if u in UNCLASSIFIED_TOKENS: return "Unclassified"
    if u in {"GYPSY","COPIA","ERV","BEL","OTHER LTRS","OTHER LTR","LTR"}: return "LTR"
    if u in {"TC1","MARINER","TC1/MARINER","TC1-MARINER"}: return "Tc1/Mariner"
    if u == "HAT": return "hAT"
    keep = {"LTR","LINE","SINE","DIRS","PLE","TRIM","LARD","TIR","MITE"}
    return u if u in keep else o

def canonical_class_from_order(order_name: str) -> str:
    if not order_name: return "Unclassified"
    u = str(order_name).strip().upper()
    if u in UNCLASSIFIED_TOKENS: return "Unclassified"
    if (u in CLASS_I_ORDERS) or (u in {"GYPSY","COPIA","ERV","BEL","OTHER LTRS","OTHER LTR"}): return "Class I"
    if u in CLASS_II_ORDERS: return "Class II"
    return "Unclassified"

def _order_sort_key(o: str) -> tuple:
    s = "" if o is None else str(o)
    return (s.upper() == "NA", s.lower())


























# ===================== AFS + χ² (summary + TSVs) =========================

def _fmt_edge_prec(x: float, ndp: int) -> str:
    s = f"{x:.{ndp}f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "-0.") else s

def make_af_bin_labels(nbins: int) -> list[str]:
    if nbins <= 0: return []
    step = 1.0 / nbins
    ndp = max(1, min(6, int(math.ceil(-math.log10(step))) + 1))
    edges = np.linspace(0.0, 1.0, nbins + 1)
    edges = np.round(edges, decimals=ndp); edges[-1] = 1.0
    labels = []
    for i in range(nbins):
        left  = _fmt_edge_prec(edges[i], ndp)
        right = _fmt_edge_prec(edges[i + 1], ndp)
        closing = ")" if i < nbins - 1 else "]"
        labels.append(f"[{left}, {right}{closing}")
    if len(set(labels)) != len(labels):
        ndp2 = min(6, ndp + 2); labels = []
        for i in range(nbins):
            left  = _fmt_edge_prec(edges[i], ndp2)
            right = _fmt_edge_prec(edges[i + 1], ndp2)
            closing = ")" if i < nbins - 1 else "]"
            labels.append(f"[{left}, {right}{closing}")
    return labels

def build_contingency_af_en(af_vals, te_labels, nbins=10):
    bins = np.linspace(0.0, 1.0, nbins + 1)
    af = np.asarray(af_vals, float); labs = np.asarray(te_labels, object)
    ok = np.isfinite(af); af = np.clip(af[ok], 0.0, 1.0); labs = labs[ok]
    vocab = {"LTR","Class I non-LTR","Class II","Unclassified"}
    labs = np.array([x if x in vocab else "Unclassified" for x in labs], dtype=object)
    bin_idx = np.digitize(af, bins, right=False) - 1; bin_idx = np.clip(bin_idx, 0, nbins - 1)
    row_order = ["LTR","Class I non-LTR","Class II","Unclassified"]
    dfc = pd.DataFrame({"Allele frequency": bin_idx, "TE category": labs})
    ct = pd.crosstab(dfc["TE category"], dfc["Allele frequency"]).reindex(
        index=row_order, columns=list(range(nbins)), fill_value=0
    )
    af_ticklabels = make_af_bin_labels(nbins)
    ct.columns = af_ticklabels
    ct.index.name = "TE category"; ct.columns.name = "Allele frequency"
    return ct, af_ticklabels

def chi2_with_residuals_safe(ct_full: pd.DataFrame):
    ct_full = ct_full.copy(); ct_full[:] = ct_full.values.astype(float)
    row_mask = (ct_full.sum(axis=1) > 0); col_mask = (ct_full.sum(axis=0) > 0)
    ct = ct_full.loc[row_mask, col_mask]
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        E_full = pd.DataFrame(0.0, index=ct_full.index, columns=ct_full.columns)
        Z_full = pd.DataFrame(0.0, index=ct_full.index, columns=ct_full.columns)
        return 0.0, 1.0, 0, E_full, Z_full
    chi2, p, dof, expected = chi2_contingency(ct.values, correction=False)
    E = pd.DataFrame(expected, index=ct.index, columns=ct.columns)
    with np.errstate(divide="ignore", invalid="ignore"):
        Z = (ct - E) / np.sqrt(E)
    E_full = pd.DataFrame(0.0, index=ct_full.index, columns=ct_full.columns)
    Z_full = pd.DataFrame(0.0, index=ct_full.index, columns=ct_full.columns)
    E_full.loc[ct.index, ct.columns] = E
    Z_full.loc[ct.index, ct.columns] = Z
    Z_full = Z_full.replace([np.inf, -np.inf], 0.0).fillna(0.0)
    return chi2, p, dof, E_full, Z_full

def make_label_palette(labels):
    base_palette = {
        "LTR": (0.20, 0.60, 0.80),
        "Class I non-LTR": (0.90, 0.60, 0.00),
        "Class II": (0.30, 0.70, 0.30),
        "Unclassified": (0.60, 0.60, 0.60),
    }
    order = ["LTR","Class I non-LTR","Class II","Unclassified"]
    labels_present = [l for l in order if l in set(labels)]
    return {lab: base_palette[lab] for lab in labels_present}, labels_present

def counts_per_bin_label(af_vals, lab_vals, bins, label_order):
    if af_vals.size == 0 or lab_vals is None or len(lab_vals) == 0:
        return np.zeros((len(bins) - 1, len(label_order)), dtype=int)
    bin_idx = np.digitize(af_vals, bins, right=False) - 1
    bin_idx = np.clip(bin_idx, 0, len(bins) - 2)
    dfc = pd.DataFrame({"bin": bin_idx, "label": lab_vals})
    counts = dfc.value_counts().reset_index(name="count")
    mat = counts.pivot(index="bin", columns="label", values="count") \
               .reindex(columns=label_order).fillna(0).astype(int)
    mat = mat.reindex(index=pd.Index(range(len(bins) - 1), name="bin"), fill_value=0)
    return mat.values

def _kde_gaussian_1d(x, grid, bw=None):
    x = np.asarray(x, dtype=float); x = x[np.isfinite(x)]
    if x.size == 0: return np.zeros_like(grid)
    if bw is None:
        if x.size > 1:
            std = np.std(x, ddof=1)
            iqr = np.subtract(*np.percentile(x, [75, 25]))
            sigma = max(std, iqr / 1.349, 1e-6)
        else:
            sigma = 1e-6
        bw = max(1.06 * sigma * (x.size ** (-1/5)), 1e-3)
    u = (grid[None,:] - x[:,None]) / bw
    dens = np.exp(-0.5 * u**2) / math.sqrt(2*math.pi)
    return dens.sum(axis=0) / (x.size * bw)

def overlay_density_curves(ax, af_vals, lab_vals, label_palette, nbins):
    if af_vals.size == 0 or lab_vals is None or len(lab_vals) == 0: return
    af = np.asarray(af_vals, dtype=float); labs = np.asarray(lab_vals, dtype=object)
    valid = np.isfinite(af); af = np.clip(af[valid], 0.0, 1.0); labs = labs[valid]
    xgrid = np.linspace(0.0, 1.0, 512); bin_width = 1.0 / nbins
    for lab, col in label_palette.items():
        sel = af[labs == lab]
        if sel.size == 0: continue
        dens = _kde_gaussian_1d(sel, xgrid, bw=None)
        y = dens * sel.size * bin_width
        ax.plot(xgrid, y, lw=1.8, color=col, alpha=0.95, solid_capstyle="round")

def grouped_afs(ax, af_vals, labels=None, label_palette=None,
                show_xlabel=False, show_ylabel=False, nbins=10):
    bins = np.linspace(0.0, 1.0, nbins + 1); ax.set_xlim(0.0, 1.0)
    if labels is None or (label_palette is None):
        af_np = np.asarray(af_vals, dtype=float); af_np = af_np[~np.isnan(af_np)]
        if af_np.size == 0: ax.text(0.5, 0.5, "no data", ha="center", va="center"); ax.set_axis_off(); return
        af_np = np.clip(af_np, 0.0, 1.0)
        counts, _ = np.histogram(af_np, bins=bins)
        lefts = bins[:-1] + 0.05 * (bins[1:] - bins[:-1]); widths = 0.90 * (bins[1:] - bins[:-1])
        ax.bar(lefts, counts, width=widths, align='edge', color="grey", edgecolor="black", linewidth=0.3)
    else:
        af_np = np.asarray(af_vals, dtype=float); lab_arr = np.array(labels, dtype=object)
        n = min(len(af_np), len(lab_arr))
        if n == 0: ax.text(0.5, 0.5, "no data", ha="center", va="center"); ax.set_axis_off(); return
        af_np = af_np[:n]; lab_arr = lab_arr[:n]
        valid = ~np.isnan(af_np); af_np = np.clip(af_np[valid], 0.0, 1.0); lab_arr = lab_arr[valid]
        lab_arr = np.array([x if x in {"LTR","Class I non-LTR","Class II","Unclassified"} else "Unclassified" for x in lab_arr], dtype=object)
        label_order = list(label_palette.keys()); counts_mat = counts_per_bin_label(af_np, lab_arr, bins, label_order)
        nbins_local = len(bins) - 1; nlab = len(label_order)
        for b in range(nbins_local):
            left = bins[b]; bw = bins[b + 1] - bins[b]; inner_left = left + 0.05 * bw; group_w = 0.90 * bw
            bar_w = group_w / max(nlab, 1)
            for i, lab in enumerate(label_order):
                x = inner_left + i * bar_w; y = counts_mat[b, i]
                ax.bar(x, y, width=bar_w, align='edge', color=label_palette.get(lab, (0.7, 0.7, 0.7)),
                       edgecolor="black", linewidth=0.3)
        overlay_density_curves(ax, af_np, lab_arr, label_palette, nbins)
    if show_xlabel:
        ax.set_xlabel("Allele frequency (TE)", fontsize=9); ax.tick_params(axis='x', labelsize=8)
    else:
        ax.set_xlabel(""); ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("Count", fontsize=9); ax.tick_params(axis='y', labelsize=8)
    else:
        ax.set_ylabel(""); ax.set_yticklabels([])


def run_chi2_tables():
    """Compute χ² contingency tables/stats for top pops and write ONE TSV per artifact."""
    class_dict, order_dict, df_pop, G_mei, var_ids_mei, G_mea, var_ids_mea, sfamily_dict = _load_all_inputs( VCF_MEI, VCF_MEA, POP_FILE, CLASSIF_FILE )

    # simplified TE labels (4-class) for AFS/χ²
    def label_for_variant(ids):
        if not ids: return "Unclassified"
        def get_cls(i):  return str(class_dict.get(i, "Unclassified")).split("|")[0]
        def get_odr(i):  return str(order_dict.get(i, "NA")).split("|")[0]
        orders  = [get_odr(i) for i in ids]; classes = [get_cls(i) for i in ids]
        if any(odr == "LTR" for odr in orders): return "LTR"
        if any(cls == "II" for cls in classes): return "Class II"
        if any(cls == "I" for cls in classes):  return "Class I non-LTR"
        return "Unclassified"

    labels_mei = np.array([label_for_variant(ids) for ids in var_ids_mei], dtype=object) if len(var_ids_mei) else np.array([], object)
    labels_mea = np.array([label_for_variant(ids) for ids in var_ids_mea], dtype=object) if len(var_ids_mea) else np.array([], object)

    pops = sorted(df_pop["Population"].dropna().unique().tolist())
    contingency_blocks = []; stats_rows = []; af_labels_any = None

    for pop in pops:
        mask_mei = (df_pop["Population"].reindex(G_mei.index.astype(str)).fillna("Unknown").values == pop) if not G_mei.empty else np.array([], bool)
        mask_mea = (df_pop["Population"].reindex(G_mea.index.astype(str)).fillna("Unknown").values == pop) if not G_mea.empty else np.array([], bool)
        idx_mei_pop = np.where((_dosage_te(G_mei.values, "MEI")[mask_mei, :] > 0).any(axis=0))[0] if not G_mei.empty else np.array([], dtype=int)
        idx_mea_pop = np.where((_dosage_te(G_mea.values, "MEA")[mask_mea, :] > 0).any(axis=0))[0] if not G_mea.empty else np.array([], dtype=int)
        def _af(G, type_name, mask, cols_idx):
            if G.empty or cols_idx.size == 0 or not np.any(mask): return np.array([])
            d = _dosage_te(G.values[:, cols_idx], type_name)[mask, :]
            m = ~np.isnan(d); n_called = m.sum(axis=0); d[~m] = 0.0; te = d.sum(axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                af = np.where(n_called > 0, te / (2.0 * n_called), np.nan)
            return np.clip(af[~np.isnan(af)], 0.0, 1.0)
        af_mei = _af(G_mei, "MEI", mask_mei, idx_mei_pop); af_mea = _af(G_mea, "MEA", mask_mea, idx_mea_pop*0 + idx_mea_pop)
        af_all = np.concatenate([af_mei, af_mea])
        labs_all = np.concatenate([labels_mei[idx_mei_pop] if idx_mei_pop.size else np.array([], object),
                                   labels_mea[idx_mea_pop] if idx_mea_pop.size else np.array([], object)])
        if af_all.size == 0 or labs_all.size == 0: continue
        ct, af_ticklabels = build_contingency_af_en(af_all, labs_all, nbins=NBINS_CHI2)
        chi2, p, dof, E, Z = chi2_with_residuals_safe(ct)
        if af_labels_any is None: af_labels_any = af_ticklabels
        ct_block = ct.copy(); ct_block.insert(0, "Population", pop); ct_block.insert(1, "TE category", ct_block.index)
        ct_block.index = pd.RangeIndex(len(ct_block)); contingency_blocks.append(ct_block)
        stats_rows.append({"Population": pop, "Chi2": chi2, "dof": dof, "p_value": p})

    if contingency_blocks:
        all_ct = pd.concat(contingency_blocks, ignore_index=True)
        cols = ["Population", "TE category"] + (af_labels_any or [])
        all_ct = all_ct[cols]
        all_ct.to_csv(os.path.join(TSV_OUTDIR, "MEGAnE_chi2_contingency_all_pops.tsv"), sep="\t", index=False)
    if stats_rows:
        stats_df = pd.DataFrame(stats_rows)[["Population","Chi2","dof","p_value"]]
        stats_df.to_csv(os.path.join(TSV_OUTDIR, "MEGAnE_chi2_stats.tsv"), sep="\t", index=False)
    print(f"📝 Wrote χ² TSVs to: {TSV_OUTDIR}")


def plot_afs_chi2_summary_2x2(pop_to_aflabs, pop_to_Z, af_ticklabels, pops_order, outdir,
                              filename="MEGAnE_AFS_plus_Chi2.png",
                              nbins_for_afs=10,
                              show_cell_grid=False):
    pops = [p for p in pops_order if (p in pop_to_aflabs) and (p in pop_to_Z)]
    if not pops: print("⚠️ No populations for summary."); return
    vmax = max(abs(pop_to_Z[p].values).max() for p in pops); vmax = max(vmax, 1.0); vmin = -vmax
    if SPECIES == 'PCOM_old':
        fig = plt.figure(figsize=(20, 13.5))
        outer = gridspec.GridSpec(3, 2, figure=fig, wspace=0.18, hspace=0.33)
    elif SPECIES == 'PPY' or SPECIES == 'PCOM':
        fig = plt.figure(figsize=(20, 23))
        outer = gridspec.GridSpec(4, 2, figure=fig, wspace=0.18, hspace=0.37)
    ims = []
    all_cats = set()
    for p in pops:
        af_all, labs, n_samples = pop_to_aflabs[p]
        if labs is not None: all_cats |= set(labs.tolist())
    cat_palette, _ = make_label_palette(list(all_cats))
    legend_handles = [Patch(facecolor=cat_palette[c], edgecolor="black", linewidth=0.3, label=c)
                      for c in cat_palette.keys()]
    letters = ["a","b","c","d","e","f","g","h"]
    for k, pop in enumerate(pops):
        sub = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[k],
                                               height_ratios=[1.9, 1.5], hspace=0.15)
        ax_top = fig.add_subplot(sub[0, 0]); ax_bot = fig.add_subplot(sub[1, 0])
        af_all, labs_all, n_samples = pop_to_aflabs[pop]
        if labs_all is None or af_all.size == 0:
            ax_top.text(0.5,0.5,"no AFS data",ha="center",va="center"); ax_top.set_axis_off()
        else:
            grouped_afs(ax_top, af_all, labels=labs_all, label_palette=cat_palette,
                        show_xlabel=False, show_ylabel=True, nbins=nbins_for_afs)
        ax_top.text(0.5, 1.10, f"{pop} (N={n_samples})", transform=ax_top.transAxes, ha="center", va="bottom", fontsize=11)
        ax_top.text(-0.12, 1.08, f"{letters[k]}.", transform=ax_top.transAxes,
                    ha="right", va="bottom", fontsize=18, fontweight="bold")
        Z = pop_to_Z[pop]
        im = ax_bot.imshow(Z.values, aspect="auto", interpolation="nearest", cmap="coolwarm", vmin=vmin, vmax=vmax)
        ims.append(im)
        ax_bot.set_xlabel("Allele frequency", fontsize=9); ax_bot.set_ylabel("TE category", fontsize=9)
        nb = Z.shape[1]; ax_bot.set_xticks(range(nb))
        ax_bot.set_xticklabels(af_ticklabels, rotation=45, ha="right", fontsize=8)
        ax_bot.set_yticks(range(len(Z.index))); ax_bot.set_yticklabels(list(Z.index), fontsize=8)
        if show_cell_grid:
            ax_bot.set_xticks(np.arange(-.5, nb, 1), minor=True)
            ax_bot.set_yticks(np.arange(-.5, len(Z.index), 1), minor=True)
            ax_bot.grid(which="minor", color="white", linewidth=0.5)
    cax = fig.add_axes([0.915, 0.22, 0.012, 0.46])
    cb  = fig.colorbar(ims[0], cax=cax)
    cb.set_label("Standardized residual", rotation=90, labelpad=8)
    cb.ax.tick_params(labelsize=8)
    fig.legend(legend_handles, [h.get_label() for h in legend_handles], title="TE category",
               loc="upper left", bbox_to_anchor=(0.905, 0.78), frameon=True,
               fontsize=8, title_fontsize=9, borderaxespad=0.0)
    plt.tight_layout(rect=[0.03, 0.05, 0.89, 0.97])
    outpath = os.path.join(outdir, filename); plt.savefig(outpath, dpi=350, bbox_inches="tight")
    plt.close(fig); print(f"🧭 Summary saved: {outpath}")


def run_afs_chi2_summary():
    """Build the 2×2 summary figure (AFS + χ² residuals)."""
    class_dict, order_dict, df_pop, G_mei, var_ids_mei, G_mea, var_ids_mea, sfamily_dict = _load_all_inputs( VCF_MEI, VCF_MEA, POP_FILE, CLASSIF_FILE )

    print(df_pop)
    
    def label_for_variant(ids):
        if not ids: return "NA"
        def get_cls(i):  return str(class_dict.get(i, "NA")).split("|")[0]
        def get_odr(i):  return str(order_dict.get(i, "NA")).split("|")[0]
        orders  = [get_odr(i) for i in ids]; classes = [get_cls(i) for i in ids]
        if any(odr == "LTR" for odr in orders): return "LTR"
        if any(cls == "II" for cls in classes): return "Class II"
        if any(cls == "I" for cls in classes):  return "Class I non-LTR"
        return "Unclassified"

    labels_mei = np.array([label_for_variant(ids) for ids in var_ids_mei], dtype=object) if len(var_ids_mei) else np.array([], object)
    labels_mea = np.array([label_for_variant(ids) for ids in var_ids_mea], dtype=object) if len(var_ids_mea) else np.array([], object)

    pops = sorted(df_pop["Population"].dropna().unique().tolist())
    pop_to_Z = {}; pop_to_aflabs = {}; af_labels_any = None

    for pop in pops:
        mask_mei = (df_pop["Population"].reindex(G_mei.index.astype(str)).fillna("Unknown").values == pop) if not G_mei.empty else np.array([], bool)
        mask_mea = (df_pop["Population"].reindex(G_mea.index.astype(str)).fillna("Unknown").values == pop) if not G_mea.empty else np.array([], bool)
        
        # nb of samples for this pop
        n_samples = mask_mei.sum() + mask_mea.sum()

        idx_mei_pop = np.where((_dosage_te(G_mei.values, "MEI")[mask_mei, :] > 0).any(axis=0))[0] if not G_mei.empty else np.array([], dtype=int)
        idx_mea_pop = np.where((_dosage_te(G_mea.values, "MEA")[mask_mea, :] > 0).any(axis=0))[0] if not G_mea.empty else np.array([], dtype=int)
        def _af(G, type_name, mask, cols_idx):
            if G.empty or cols_idx.size == 0 or not np.any(mask): return np.array([])
            d = _dosage_te(G.values[:, cols_idx], type_name)[mask, :]
            m = ~np.isnan(d); n_called = m.sum(axis=0); d[~m] = 0.0; te = d.sum(axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                af = np.where(n_called > 0, te / (2.0 * n_called), np.nan)
            return np.clip(af[~np.isnan(af)], 0.0, 1.0)
        af_mei = _af(G_mei, "MEI", mask_mei, idx_mei_pop)
        af_mea = _af(G_mea, "MEA", mask_mea, idx_mea_pop)
        af_all = np.concatenate([af_mei, af_mea])
        labs_all = np.concatenate([labels_mei[idx_mei_pop] if idx_mei_pop.size else np.array([], object),
                                   labels_mea[idx_mea_pop] if idx_mea_pop.size else np.array([], object)])
        if af_all.size == 0 or labs_all.size == 0: continue
        pop_to_aflabs[pop] = (af_all, labs_all, n_samples)
        ct, af_ticklabels = build_contingency_af_en(af_all, labs_all, nbins=NBINS_CHI2)
        chi2, p, dof, E, Z = chi2_with_residuals_safe(ct)
        pop_to_Z[pop] = Z
        if af_labels_any is None: af_labels_any = af_ticklabels
        
    merge_pairs = [
        ("cauc", "pyra", "cauc+pyra"),
        ("comm_Dessert", "comm_Perry", "comm_Dessert+comm_Perry"),
    ]

    for popA, popB, new_name in merge_pairs:
        # Build merged masks (union of individuals from both populations)
        mask_mei_A = (df_pop["Population"].reindex(G_mei.index.astype(str)).fillna("Unknown").values == popA) if not G_mei.empty else np.array([], bool)
        mask_mei_B = (df_pop["Population"].reindex(G_mei.index.astype(str)).fillna("Unknown").values == popB) if not G_mei.empty else np.array([], bool)
        mask_mei   = mask_mei_A | mask_mei_B
    
        mask_mea_A = (df_pop["Population"].reindex(G_mea.index.astype(str)).fillna("Unknown").values == popA) if not G_mea.empty else np.array([], bool)
        mask_mea_B = (df_pop["Population"].reindex(G_mea.index.astype(str)).fillna("Unknown").values == popB) if not G_mea.empty else np.array([], bool)
        mask_mea   = mask_mea_A | mask_mea_B
    
        # Number of individuals in merged population
        n_samples = mask_mei.sum() + mask_mea.sum()
    
        # Identify loci present in the merged population
        idx_mei_pop = np.where((_dosage_te(G_mei.values, "MEI")[mask_mei, :] > 0).any(axis=0))[0] if not G_mei.empty else np.array([], dtype=int)
        idx_mea_pop = np.where((_dosage_te(G_mea.values, "MEA")[mask_mea, :] > 0).any(axis=0))[0] if not G_mea.empty else np.array([], dtype=int)
    
        # Recompute allele frequencies using merged masks
        af_mei = _af(G_mei, "MEI", mask_mei, idx_mei_pop)
        af_mea = _af(G_mea, "MEA", mask_mea, idx_mea_pop)
    
        # Concatenate allele frequencies and labels
        af_all = np.concatenate([af_mei, af_mea])
        labs_all = np.concatenate([
            labels_mei[idx_mei_pop] if idx_mei_pop.size else np.array([], object),
            labels_mea[idx_mea_pop] if idx_mea_pop.size else np.array([], object)
        ])
    
        # Skip if no data
        if af_all.size == 0 or labs_all.size == 0:
            continue
    
        # Store merged population in dictionaries
        pop_to_aflabs[new_name] = (af_all, labs_all, n_samples)
        ct, af_ticklabels = build_contingency_af_en(af_all, labs_all, nbins=NBINS_CHI2)
        chi2, p, dof, E, Z = chi2_with_residuals_safe(ct)
        pop_to_Z[new_name] = Z
    
        # Add merged population to the list (keep original ones too)
        if new_name not in pops:
            pops.append(new_name)
    
    # Finally call the plotting function
    if pop_to_Z and pop_to_aflabs and af_labels_any:
        plot_afs_chi2_summary_2x2(pop_to_aflabs, pop_to_Z, af_labels_any, pops, FIG_OUTDIR)




































############################################################################
# ===================== Fixed / Polymorphic by Class→Order =================

def compute_me_fixed_polymorphic_counts_by_class_order(
    G_mei, var_ids_mei, G_mea, var_ids_mea,
    class_dict, order_dict, sfamily_dict,
    df_pop
):
    """
    Create count tables with columns MultiIndex (TE class, TE order):
      - fixed_df      : ME present only in one population (private)
      - polymorphic_df: ME present in ≥2 populations (counted in each carrying pop)
    """

    def _infer_base_class(ids):
        if not ids:
            return "NA"
        cls = [str(class_dict.get(i, "NA")).split("|")[0] for i in ids]
        if any(c == "I" for c in cls):
            return "Class I"
        if any(c == "II" for c in cls):
            return "Class II"
        return "NA"

    #### TO RM
    # def _infer_raw_order(ids):
    #     if not ids:
    #         return "NA"
    #     os_ = [str(order_dict.get(i, "NA")).split("|")[0] for i in ids]
    #     if any(o == "LTR" for o in os_):
    #         return "LTR"
    #     for o in os_:
    #         if o and o != "NA":
    #             return o
    #     return "NA"
    #### TO RM

    def _infer_raw_order(ids):
        if not ids:
            return "NA"
        # Fetch order and family
        os_ = [str(order_dict.get(i, "NA")).split("|")[0].upper() for i in ids]
        sf_ = [str(sfamily_dict.get(i, "NA")).upper() for i in ids]
        # Split LTRs into Copia, Gypsy and other
        if any("LTR" in o for o in os_):
            if any("COPIA" in s for s in sf_):
                return "LTR - Copia"
            elif any("GYPSY" in s for s in sf_):
                return "LTR - Gypsy"
            else:
                return "LTR - Other"
        # Other orders wich are not LTR
        for o in os_:
            if o and o != "NA":
                return o
        return "NA"

    base_cls_mei = [_infer_base_class(ids) for ids in var_ids_mei] if len(var_ids_mei) else []
    raw_ord_mei  = [_infer_raw_order(ids)  for ids in var_ids_mei] if len(var_ids_mei) else []
    pops_mei     = variant_population_sets(G_mei, "MEI", df_pop) if G_mei is not None else []

    base_cls_mea = [_infer_base_class(ids) for ids in var_ids_mea] if len(var_ids_mea) else []
    raw_ord_mea  = [_infer_raw_order(ids)  for ids in var_ids_mea] if len(var_ids_mea) else []
    pops_mea     = variant_population_sets(G_mea, "MEA", df_pop) if G_mea is not None else []

    classes_all, orders_all = [], []
    for cls0, ord0 in zip(base_cls_mei + base_cls_mea, raw_ord_mei + raw_ord_mea):
        if not ord0:
            ord0 = "NA"
        u = str(ord0).strip().upper()

        if u in UNCLASSIFIED_TOKENS:
            # Cas explicite "Unclassified"
            classes_all.append("Unclassified")
            orders_all.append("Unclassified")

        elif u == "NA":
            if cls0 in {"Class I", "Class II"}:
                # Classe connue mais ordre manquant → créer "NA" dans ce bloc
                classes_all.append(cls0)
                orders_all.append("NA")
            else:
                # Ni classe ni ordre connus → vrai Unclassified
                classes_all.append("Unclassified")
                orders_all.append("Unclassified")
        else:
            # Cas normal : ordre reconnu
            #classes_all.append(canonical_class_from_order(u))
            #orders_all.append(canonicalize_order_label(ord0))
            if any(x in u for x in ["LTR - COPIA", "LTR - GYPSY", "LTR - OTHER"]):
                classes_all.append("Class I")
            else:
                classes_all.append(canonical_class_from_order(u))
            orders_all.append(canonicalize_order_label(ord0))
    
    popsets_all = pops_mei + pops_mea

    pops_seen = set()
    for s in popsets_all:
        pops_seen.update(s)
    pop_list = sorted(pops_seen)

    class_rank = {"Class I": 0, "Class II": 1, "Unclassified": 2}
    pairs_sorted = sorted(
        set(zip(classes_all, orders_all)),
        key=lambda t: (class_rank.get(t[0], 99),) + _order_sort_key(t[1])
    )
    cols = pd.MultiIndex.from_tuples(pairs_sorted, names=["TE class", "TE order"])

    fixed_df       = pd.DataFrame(0, index=pop_list, columns=cols, dtype=int)
    polymorphic_df = pd.DataFrame(0, index=pop_list, columns=cols, dtype=int)

    for te_cls, te_ord, pops_set in zip(classes_all, orders_all, popsets_all):
        if not pops_set:
            continue
        key = (te_cls, te_ord)
        if key not in fixed_df.columns:
            continue
        if len(pops_set) == 1:
            the_pop = next(iter(pops_set))
            if the_pop in fixed_df.index:
                fixed_df.at[the_pop, key] += 1
        else:
            for p in pops_set:
                if p in polymorphic_df.index:
                    polymorphic_df.at[p, key] += 1

    return fixed_df, polymorphic_df


def _sort_multiindex_cols(cols: pd.MultiIndex) -> pd.MultiIndex:
    class_rank = {"Class I": 0, "Class II": 1, "Unclassified": 2}
    tuples = sorted(cols.to_list(), key=lambda t: (class_rank.get(t[0], 99),) + _order_sort_key(t[1]))
    return pd.MultiIndex.from_tuples(tuples, names=cols.names)

# ===================== Single heatmap ============================

def make_me_classorder_heatmap(
    pop_list,
    newick_file,
    G_mei, var_ids_mei, G_mea, var_ids_mea,
    class_dict, order_dict, sfamily_dict,
    df_pop,
    FIG_OUTDIR,
    Sp
):

    """
    Heatmap du ratio Fixed/Polymorphic par population et par (TE class, TE order)
    - Couleur : ratio Fixed/Polymorphic (borné entre 0 et 1)
    - Texte   : "Fixed/Polymorphic" (ex: 12/10)
    - Organisation : Class I → Class II → Unclassified
    """
    fixed_df, poly_df = compute_me_fixed_polymorphic_counts_by_class_order(
        G_mei, var_ids_mei, G_mea, var_ids_mea,
        class_dict, order_dict, sfamily_dict,
        df_pop
    )
    #print(fixed_df)
    #print(poly_df)

    if fixed_df.empty and poly_df.empty:
        print("⚠️ No data to plot for single heatmap.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()


    # --- Harmonisation des lignes/colonnes ---
    cols_union = _sort_multiindex_cols(fixed_df.columns.union(poly_df.columns))
    fixed_df = fixed_df.reindex(index=pop_list, columns=cols_union, fill_value=0)
    poly_df  = poly_df.reindex(index=pop_list, columns=cols_union, fill_value=0)
    
    fixed_T = fixed_df.T
    poly_T  = poly_df.T


    # --- Calcul du ratio Fixed / Polymorphic ---
    ratio_df = pd.DataFrame(index=fixed_T.index)
    for pop in pop_list:
        f = fixed_T[pop].astype(float)
        p = poly_T[pop].astype(float)
        ratio = np.divide(p, (p+f), out=np.zeros_like(f), where=(p > 0))
        ratio_df[pop] = ratio.clip(0, 1)  # borné entre 0 et 1


    # --- Text "Polymorphic/(Polymorphic+Fixed)" ---
    annots = pd.DataFrame(index=fixed_T.index)
    for pop in pop_list:
        p = poly_T[pop].astype(int)
        f = fixed_T[pop].astype(int)
        #total = f + p
        annots[pop] = [f"{po}/{fo}" for po, fo in zip(p, f)]

    plot_df = ratio_df.fillna(0)


    # --- Load the Newick dendogram ---
    tree = Phylo.read(newick_file, "newick")

    # --- Tracer dendrogramme + heatmap dans une même figure ---
    fig = plt.figure(figsize=(10, 12)) #7, 12
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 5, 1.5], hspace=0.2)

    # Dendrogramme en haut
    ax_tree = fig.add_subplot(gs[0])
    Phylo.draw(tree, do_show=False, axes=ax_tree)
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])
    ax_tree.set_frame_on(False)
    ax_tree.set_xlabel("")
    ax_tree.set_ylabel("")

    # --- Color leaf labels in the dendrogram based on Crop_or_Wild ---
    color_map = {"Wild": "#158B22", "Cultivated": "#545EFA", "Rootstock": "orange"}
    color_map2 = {"Wild": "#158B22", "Cultivated": "#545EFA"}
    # Create a unique mapping Population → Crop_or_Wild (take the first occurrence)
    pop_to_crop = df_pop.drop_duplicates(subset="Population").set_index("Population")["Crop_or_Wild"]
    if Sp == "PPY":
        pop_to_crop = pd.Series(
            {
                "betu": "Rootstock",
                "pash": "Wild",
                "pyri_JP": "Cultivated",
                "Sand_CN": "Cultivated",
                "Sand_CN-SE": "Cultivated",
                "Sand_CN-SW": "Cultivated",
                "ussu": "Wild",
                "White": "Cultivated",
            },
            name="Crop_or_Wild"
        )
    for text_obj in ax_tree.findobj(match=lambda obj: isinstance(obj, plt.Text)):
        label = text_obj.get_text().strip()
        crop_status = pop_to_crop.get(label)
        # Apply color (default gray if unknown)
        text_obj.set_color(color_map.get(crop_status, "gray"))

    letter = ''
    if Sp == 'PCOM':
        letter = 'c'
    elif Sp == "PPY":
        letter = "d"
    fig.text(0, 0.9, letter, fontsize=18, fontweight="bold", va="top", ha="left")
 
    # --- Create legend handles manually ---
    if Sp == 'PCOM':
        color_map = color_map2
    legend_handles = [
        mlines.Line2D([], [], color=color, marker='o', linestyle='None', label=label)
        for label, color in color_map.items()
    ]

    fontsize = 10
    loc='best'
    if Sp == 'PPY':
        fontsize = 8
        
    # --- Add legend to your dendrogram axes ---
    ax_tree.legend(
        handles=legend_handles,
        loc=loc,
        frameon=False,
        fontsize=8
    )

    # --- Create heatmap ---
    #fig, ax = plt.subplots(figsize=(10, 6))
    ax_heat = fig.add_subplot(gs[1])
    #print(plot_df)
    
    sns.heatmap(
        plot_df,
        cmap="viridis", #"Reds_r"
        vmin=0, vmax=1,
        annot=annots,
        fmt="",
        annot_kws={"fontsize": fontsize},
        cbar_kws={"label": "Shared/Fixed ratio", "pad": 0.02},
        ax=ax_heat
    )

    # --- y axe
    classes = [idx[0] for idx in plot_df.index]
    orders  = [idx[1] for idx in plot_df.index]
    ax_heat.set_yticklabels(orders, fontsize=10)
    ax_heat.set_ylabel("", fontsize=11)

    # --- x axe
    # ---- Color x-axis labels based on df_pop["Crop_or_Wild"] ----
    # Get x tick labels from the heatmap (Population order)
    x_labels = ax_heat.get_xticklabels()
    pop_order = [lbl.get_text() for lbl in x_labels]
    # Map populations to crop type safely
    crop_status = pop_to_crop.reindex(pop_order)
    
    # Now colors will match correctly
    for lbl, status in zip(x_labels, crop_status):
        lbl.set_color(color_map.get(status, "black"))
    ax_heat.set_xticklabels(pop_list, rotation=0, fontsize=fontsize+1)
    ax_heat.set_xlabel("Population", fontsize=9)

    # --- Couleurs par classe ---
    row_colors = {"Class I": "darkgreen", "Class II": "#2ca02c", "Unclassified": "#7f7f7f"}
    for i, cls in enumerate(classes):
        ax_heat.get_yticklabels()[i].set_color(row_colors.get(cls, "black"))

    # --- Lignes séparatrices entre classes ---
    row_edges = [i for i in range(1, len(classes)) if classes[i] != classes[i-1]]
    x0, x1 = ax_heat.get_xlim()
    for y in row_edges:
        ax_heat.hlines(y, x0, x1, colors="white", linewidth=1.2, zorder=3)

    # --- En-têtes de classes sur la gauche ---
    block_centers, block_labels = [], []
    start = 0
    for i in range(1, len(classes) + 1):
        if i == len(classes) or classes[i] != classes[i - 1]:
            end = i
            block_centers.append((start + end - 1) / 2.0)
            block_labels.append(classes[start])
            start = i

    # --- Class headers on the far left ---
    ax_left = ax_heat.twinx()
    ax_left.set_ylim(ax_heat.get_ylim())
    ax_left.set_yticks(block_centers)
    ax_left.set_yticklabels(block_labels, fontsize=9, fontweight="bold")
    
    # Move axis to left side
    ax_left.yaxis.set_ticks_position("left")
    ax_left.yaxis.set_label_position("left")
    
    # Hide the right spine completely
    ax_left.spines["right"].set_visible(False)
    ax_left.spines["left"].set_visible(True)
    
    # Shift the left spine further left (outside the heatmap)
    ax_left.spines["left"].set_position(("outward", 70))  # décalage vers la gauche
    ax_left.tick_params(axis="y", length=0, pad=5)
    
    # Couleur cohérente par classe
    for lbl, cls in zip(ax_left.get_yticklabels(), block_labels):
        lbl.set_color(row_colors.get(cls, "black"))

    # Align all panels
    plt.tight_layout(rect=[0.18, 0.1, 0.95, 0.92])
    os.makedirs(FIG_OUTDIR, exist_ok=True)
    outpath = os.path.join(FIG_OUTDIR, "MEGAnE_ME_byClassOrder_heatmap.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"🗺️ Ratio heatmap saved: {outpath}")
    
    return fixed_df, poly_df, plot_df
 
    
    
def lineplots_check_quality(G_mei, G_mea, Sp, df_pop, pop_list, FIG_OUTDIR):
    '''
        Per-sample ME counts, missing gentoypes and depth
        For Supp Mat
    '''
    # Combine MEA and MEI counts (simply add them)
    if G_mei is not None and G_mea is not None:
        # Both are presence/absence genotype matrices (samples x variants)
        total_me = G_mei.sum(axis=1) + G_mea.sum(axis=1)
    elif G_mei is not None:
        total_me = G_mei.sum(axis=1)
    elif G_mea is not None:
        total_me = G_mea.sum(axis=1)
    else:
        total_me = pd.Series(dtype=float)
    print(total_me)
    
    # Use index directly since df_pop is already indexed by ID_vcf
    sample_to_pop = df_pop["Population"]
    sample_to_type = df_pop["Crop_or_Wild"]
    df_pop["Depth"] = df_pop["Depth"].astype(float)
    sample_to_depth = df_pop["Depth"]
    #print(df_pop)
    
    # Build dataframe: sample → population → crop status → ME count
    df_me_counts = pd.DataFrame({
        "Sample": total_me.index,
        "ME_count": total_me.values,
        "Population": sample_to_pop.reindex(total_me.index).values,
        "Crop_or_Wild": sample_to_type.reindex(total_me.index).values,
        "Depth": sample_to_depth.reindex(total_me.index).values
    }).dropna(subset=["Population", "Crop_or_Wild"])

    # Ensure populations appear in same order as heatmap
    df_me_counts["Population"] = pd.Categorical(df_me_counts["Population"], categories=pop_list, ordered=True)
    df_me_counts = df_me_counts.sort_values(["Population", "Sample", "Depth"])
    #print(df_me_counts)
    
    # Convert genotype matrices to presence/absence (0/1) format:
    # - 0 → no insertion (0/0)
    # - 1 → insertion present on at least one allele (0/1 or 1/1)
    # NaN values are treated as 0 (absence)
    G_mei_bin = (G_mei > 0).astype(int)
    G_mea_bin = (G_mea > 0).astype(int)
    sample = G_mei.index[0]
    print("Avant :", G_mei.loc[sample].sum())
    print("Après :", G_mei_bin.loc[sample].sum())
    # Count the total number of TE insertion sites per sample (presence/absence basis)
    total_me = G_mei_bin.sum(axis=1) + G_mea_bin.sum(axis=1)
    #ax_bar  = fig.add_subplot(gs[2])
    fig_suppMat, ax_bar = plt.subplots(figsize=(10, 4))
    # Plot barplot of ME counts per sample, grouped by population
    sns.lineplot(
        data=df_me_counts,
        x="Sample",
        y="ME_count",
        #order=pop_list,
        ax=ax_bar,
        #hue="Crop_or_Wild",
        #palette=color_map,
        color="blue",
        marker="o",
        markersize=3,
        linewidth=1,
        legend=False
    )
    ax_bar_twin = ax_bar.twinx()
    sns.lineplot(
        data=df_me_counts,
        x="Sample",
        y="Depth",
        color="red",
        linestyle="-",
        linewidth=1,
        marker="o",
        markersize=3,
        ax=ax_bar_twin,
        legend=False
    )
    
    # === Add vertical separators and population names ===
    pop_sizes = df_me_counts["Population"].value_counts().reindex(pop_list).fillna(0).astype(int)
    pos = 0
    group_centers = []
    for size in pop_sizes:
        # Calculer le centre du groupe (pour placer le nom)
        center = pos + size / 2 - 0.5
        group_centers.append(center)
        pos += size
    # Tracer les lignes verticales entre populations
    pos = 0
    for size in pop_sizes[:-1]:  # skip last
        pos += size
        ax_bar.axvline(x=pos - 0.5, color="#333333", linewidth=1.5, zorder=10)
    # Ajouter les noms des populations sous l'axe X (avec légère rotation)
    y_min, y_max = ax_bar.get_ylim()
    for pop_name, center in zip(pop_list, group_centers):
        ax_bar.text(center, y_min - 0.3 * (y_max - y_min),  # un peu sous la figure
                    pop_name, ha='center', va='top', fontsize=11, fontweight='medium',
                    rotation=30, rotation_mode='anchor')

    # # === Count missing genotypes (NaN or ./.) per sample ===
    # # Any NaN in genotype matrices means missing data.
    # # We'll sum across variants for each sample.
    # missing_mei = G_mei.isna().sum(axis=1)
    # missing_mea = G_mea.isna().sum(axis=1)
    # # Combine MEI and MEA missing counts per sample
    # total_missing = missing_mei + missing_mea
    # # === Plot as a line on the same axis (secondary color) ===
    # ax_bar.plot(range(len(total_missing)), total_missing,
    #             color="grey", linewidth=0.75, linestyle="-", markersize=3, label="Missing genotypes", marker="o")
    # Formatting
    ax_bar.set_xlabel("Samples", fontsize=10) #grouped by population as the heatmap
    ax_bar.set_ylabel("# TE insertions", color="blue", fontsize=10)
    ax_bar.tick_params(axis="y", colors="blue")
    ax_bar.tick_params(axis="x", labelsize=6, rotation=90)
    ax_bar.tick_params(axis="y", labelsize=7)
    ax_bar.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5)
    ax_bar_twin.set_ylabel("Depth", color="red")
    ax_bar_twin.tick_params(axis="y", colors="red")
    #if Sp == 'PPY':
    # ax_bar.legend(
    #     loc="upper left",
    #     bbox_to_anchor=(0.77, 1.09),
    #     borderaxespad=0,
    #     frameon=True, facecolor="white", framealpha=1.0, fontsize=10
    # )
    fig_suppMat.savefig(os.path.join(FIG_OUTDIR, "nbMEs_depth_per_sample.png"), dpi=300, bbox_inches="tight")
    

    #------------------------------------------------------------------------------------------------



def run_me_classorder_heatmap( VCF_MEI, VCF_MEA, pop_list, newick_file, FIG_OUTDIR, Sp, POP_FILE, CLASSIF_FILE ):
    """Generate the Class and Order heatmap"""
    class_dict, order_dict, df_pop, G_mei, var_ids_mei, G_mea, var_ids_mea, sfamily_dict = _load_all_inputs( VCF_MEI, VCF_MEA, POP_FILE, CLASSIF_FILE ) # G are presence/absence genotype matrix avc var_ids the variant IDs

    # version I kept
    #G_mei_filtered, G_mea_filtered = compute_maf_and_missing(G_mei, G_mea, maf_threshold=0, missing_threshold=1) #(MAF > 0.05, missing < 0.1)
    # version I do not keep
    #G_mei_filtered, G_mea_filtered = compute_maf_and_missing(G_mei, G_mea, maf_threshold=0, missing_threshold=0.1) #(MAF > 0.05, missing < 0.1)
    
    make_me_classorder_heatmap(
        pop_list,
        newick_file,
        G_mei, var_ids_mei, G_mea, var_ids_mea,
        class_dict, order_dict, sfamily_dict,
        df_pop,
        FIG_OUTDIR,
        Sp
    )
    lineplots_check_quality(G_mei, G_mea, Sp, df_pop, pop_list, FIG_OUTDIR)

# ===================== Fixed / Polymorphic by Class→Order =================
############################################################################







































def run_me_012_clustermap(
    outdir,
    class_dict, order_dict, df_pop,
    G_mei, var_ids_mei, G_mea, var_ids_mea,
    filename="MEGAnE_ME012_clustermap.png",
    min_present=1,                     # keep TE loci seen in >= min_present samples
    max_present_prop=1,             # drop TE loci present in ~all samples
    max_variants=None,                 # cap number of TE loci (speed/clarity)
    group_by_locus=True                # aggregate variants sharing the same ME locus
):
    """
    Clustermap of TE loci (rows) × samples (columns) on diploid genotype dosage.
      - Values: 0 (absent), 1 (heterozygote), 2 (homozygote)
      - Colors: 0 = white, 1 = gray (#BDBDBD), 2 = black
      - Both columns (samples) and rows (TEs) are clustered
      - Missing genotypes are imputed to 0 for clustering and display
      - Column color strip encodes Population (custom colors provided)
      - Row color strip encodes TE class (Class I / Class II / Unclassified)
      - Right-side legend panel shows separate legends for Population and TE class
      - Sample names are shown on x-ticks; TE labels are hidden on y-ticks
    """
    os.makedirs(outdir, exist_ok=True)

    # ---- Helpers: base class from IDs and locus ID from ID list ---------------------
    def infer_base_class_from_ids(ids):
        if not ids:
            return "Unclassified"
        cs = [str(class_dict.get(i, "Unclassified")).split("|")[0] for i in ids]
        if any(c == "I" for c in cs):   return "Class I"
        if any(c == "II" for c in cs):  return "Class II"
        return "Unclassified"

    # infer order from ME IDs (keeps 'LTR' priority; falls back to first non-NA)
    def infer_order_from_ids(ids):
        if not ids:
            return "NA"
        os_ = [str(order_dict.get(i, "NA")).split("|")[0] for i in ids]
        if "LTR" in os_:
            ord0 = "LTR"
        else:
            ord0 = next((o for o in os_ if o and o != "NA"), "NA")
        # If you already have canonicalize_order_label(...), use it; else keep ord0
        try:
            return canonicalize_order_label(ord0) if ord0 != "NA" else "NA"
        except NameError:
            return ord0

    def build_locus_id(ids, fallback):
        if ids:
            uniq = sorted({str(x) for x in ids if str(x) not in {"", "nan", "None", "NaN"}})
            if uniq:
                return "|".join(uniq)
        return fallback

    # ---- Build 0/1/2 dosage matrices (samples × variants) -------------------------
    X_mei = pd.DataFrame(
        _dosage_te(G_mei.values, "MEI"),
        index=G_mei.index,
        columns=[f"MEI_{c}" for c in G_mei.columns]
    ) if (G_mei is not None and not G_mei.empty) else pd.DataFrame()

    X_mea = pd.DataFrame(
        _dosage_te(G_mea.values, "MEA"),
        index=G_mea.index,
        columns=[f"MEA_{c}" for c in G_mea.columns]
    ) if (G_mea is not None and not G_mea.empty) else pd.DataFrame()

    if X_mei.empty and X_mea.empty:
        print("⚠️ No MEI/MEA genotypes for clustermap.")
        return

    # Align samples and concatenate
    all_samples = sorted(set(X_mei.index) | set(X_mea.index))
    X = pd.concat([X_mei.reindex(all_samples), X_mea.reindex(all_samples)], axis=1).astype(float)  # samples × variants

    # ---- Variant metadata: TE base class + locus ID --------------------------------
    varnames_mei = list(X_mei.columns)
    varnames_mea = list(X_mea.columns)
    te_class_mei = [infer_base_class_from_ids(var_ids_mei[j]) for j in range(len(var_ids_mei))] if len(var_ids_mei) else []
    te_class_mea = [infer_base_class_from_ids(var_ids_mea[j]) for j in range(len(var_ids_mea))] if len(var_ids_mea) else []
    locus_mei    = [build_locus_id(var_ids_mei[j], varnames_mei[j]) for j in range(len(var_ids_mei))] if len(var_ids_mei) else []
    locus_mea    = [build_locus_id(var_ids_mea[j], varnames_mea[j]) for j in range(len(var_ids_mea))] if len(var_ids_mea) else []
    order_mei    = [infer_order_from_ids(var_ids_mei[j]) for j in range(len(var_ids_mei))] if len(var_ids_mei) else []
    order_mea    = [infer_order_from_ids(var_ids_mea[j]) for j in range(len(var_ids_mea))] if len(var_ids_mea) else []

    meta = pd.DataFrame({
        "varname":  varnames_mei + varnames_mea,
        "TE_class": te_class_mei + te_class_mea,
        "TE_order": order_mei    + order_mea,
        "locus_id": locus_mei    + locus_mea
    }).set_index("varname").reindex(X.columns)

    # ---- Aggregate by locus if requested -------------------------------------------
    X_var = X.T  # variants × samples
    if group_by_locus and not meta["locus_id"].isna().all():
        X_var = X_var.groupby(meta["locus_id"]).max()             # genotype aggregation = max(0/1/2)
        row_classes = meta.groupby("locus_id")["TE_class"].agg(lambda s: s.value_counts().idxmax()).reindex(X_var.index).fillna("Unclassified")
        row_orders  = meta.groupby("locus_id")["TE_order"].agg(lambda s: s.value_counts().idxmax()).reindex(X_var.index).fillna("NA")
    else:
        X_var.index.name = "locus_or_variant"
        row_classes = meta["TE_class"].reindex(X_var.index).fillna("Unclassified")
        row_orders  = meta["TE_order"].reindex(X_var.index).fillna("NA")

    # ---- Filter loci to stabilize clustering ---------------------------------------
    present = (X_var >= 1).astype(float)
    n_samples = X_var.shape[1]
    upper_limit = n_samples if (max_present_prop is None) else (max_present_prop * n_samples)
    keep = (present.sum(axis=1) >= float(min_present)) & (present.sum(axis=1) <= upper_limit)
    X_var = X_var.loc[keep]
    row_classes = row_classes.loc[X_var.index]
    row_orders = row_orders.loc[X_var.index]

    if X_var.empty:
        print("⚠️ All TE loci filtered out for clustermap.")
        return

    if (max_variants is not None) and (X_var.shape[0] > max_variants):
        prev = present.loc[X_var.index].mean(axis=1)
        score = prev * (1.0 - prev)                               # favor mid-prevalence
        X_var = X_var.loc[score.sort_values(ascending=False).index[:max_variants]]
        row_classes = row_classes.loc[X_var.index]

    # ---- Impute missing to 0; coerce to {0,1,2} ------------------------------------
    X_plot = X_var.fillna(0.0).clip(0, 2).round().astype(int)     # rows = loci, cols = samples

    # ---- Column color strip: Population (use your exact colors) --------------------
    pops = df_pop["Population"].reindex(X_plot.columns).fillna("Unknown")
    for p in pops.unique():
        if p not in pop_color_map:
            pop_color_map[p] = "#9e9e9e"
    col_colors = pops.map(pop_color_map)

    # ---- Row color strip: TE class -------------------------------------------------
    class_colors = {"Class I": "#1f77b4", "Class II": "#2ca02c", "Unclassified": "#7f7f7f"}

    # ---- Row color strip: TE order --------------------------------------------------
    # Order colors (one color per distinct order)
    ord_levels      = sorted(row_orders.unique().tolist())
    order_palette   = sns.color_palette("tab20", n_colors=max(3, len(ord_levels)))
    order_color_map = {o: order_palette[i % len(order_palette)] for i, o in enumerate(ord_levels)}
    # Two-column row color strip (will render side-by-side)
    row_colors_df = pd.DataFrame({
        "TE class": row_classes.map(class_colors),
        "TE order": row_orders.map(order_color_map)
    }, index=X_plot.index)

    # ---- Discrete colormap for 0/1/2 ----------------------------------------------
    cmap_012 = ListedColormap(["#FFFFFF", "#BDBDBD", "#000000"])
    norm_012 = BoundaryNorm(boundaries=[-0.5, 0.5, 1.5, 2.5], ncolors=3, clip=True)

    # ---- Clustermap (default seaborn layout); show sample names on x-axis ----------
    g = sns.clustermap(
        X_plot,                           # rows = TE loci, cols = samples
        cmap=cmap_012, norm=norm_012,
        row_cluster=True, col_cluster=True,
        metric="euclidean", method="ward", #cityblock
        col_colors=col_colors,            # population strip
        row_colors=row_colors_df,            # TE class strip
        xticklabels=True, yticklabels=False,   # show sample names; hide TE labels
        cbar_kws={"label": "Genotype dosage (0/1/2)", "ticks": [0, 1, 2]}
    )
    g.cax.set_yticks([0, 1, 2])
    g.cax.set_yticklabels(["0", "1", "2"])

    # Format x-tick labels (samples): rotate for readability
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="center", va="top", fontsize=7)

    # ---- Right-side legend panel, with two separate legends (Population, TE class) --
    pop_handles   = [Patch(facecolor=pop_color_map[p], edgecolor="none", label=p) for p in sorted(set(pops.unique()))]
    class_handles = [Patch(facecolor=class_colors[c], edgecolor="none", label=c) for c in ["Class I", "Class II", "Unclassified"]]
    order_handles = [Patch(facecolor=order_color_map[o], edgecolor="none", label=o) for o in ord_levels]

    g.fig.canvas.draw()  # ensure positions are current
    heat_pos = g.ax_heatmap.get_position()   # [x0, y0, x1, y1] in figure coords
    panel_w  = 0.13                          # width of the right legend panel (figure fraction)
    gap      = 0.10                          # small horizontal gap

    panel_left = heat_pos.x1 + gap
    panel_bottom, panel_height = heat_pos.y0, (heat_pos.y1 - heat_pos.y0)
    panel_ax = g.fig.add_axes([panel_left, panel_bottom, panel_w, panel_height])
    panel_ax.axis("off")

    # Split panel into two sub-axes (Population on top, TE class below)
    top_ax  = g.fig.add_axes([panel_left, panel_bottom + panel_height*0.68, panel_w, panel_height*0.30]); top_ax.axis("off")
    mid_ax  = g.fig.add_axes([panel_left, panel_bottom + panel_height*0.36, panel_w, panel_height*0.30]); mid_ax.axis("off")
    bot_ax  = g.fig.add_axes([panel_left, panel_bottom + panel_height*0.04, panel_w, panel_height*0.30]); bot_ax.axis("off")
    
    top_ax.legend(handles=pop_handles,
                  labels=[h.get_label() for h in pop_handles],
                  title="Population", loc="center", frameon=False)
    
    mid_ax.legend(handles=class_handles,
                  labels=[h.get_label() for h in class_handles],
                  title="TE class", loc="center", frameon=False)
    
    bot_ax.legend(handles=order_handles,
                  labels=[h.get_label() for h in order_handles],
                  title="TE order", loc="center", frameon=False)

    # Save
    outpath = os.path.join(outdir, filename)
    g.fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"🧬 0/1/2 clustermap saved: {outpath}")

def run_me_012_clustermap_only():
    """Wrapper: load inputs and produce the 0/1/2 clustermap in one call."""
    class_dict, order_dict, df_pop, G_mei, var_ids_mei, G_mea, var_ids_mea, sfamily_dict = _load_all_inputs( VCF_MEI, VCF_MEA, POP_FILE, CLASSIF_FILE )
    run_me_012_clustermap(
        outdir=FIG_OUTDIR,
        class_dict=class_dict, order_dict=order_dict, df_pop=df_pop,
        G_mei=G_mei, var_ids_mei=var_ids_mei, G_mea=G_mea, var_ids_mea=var_ids_mea
    )





















# --------------------------- PLINK2 helpers ---------------------------

MODULE_LOADS = "module load all gencore/3 && module load plink2/2.00a5.10"

def _ensure_plink2():
    """Charge les modules et vérifie plink2."""
    cmd = f'{MODULE_LOADS} && plink2 --version'
    subprocess.run(["bash", "-lc", cmd], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def _run_plink_pca(vcf_path: str, outprefix: str,
                   maf: float = 0.05, geno: float = 0.10,
                   thin_frac: float | None = None, threads: int = 8,
                   memory_mb: int = 8192, n_pc: int = 10):
    """Exécute PLINK2 --pca et retourne chemins eigenvec/eigenval."""
    _ensure_plink2()
    os.makedirs(os.path.dirname(outprefix), exist_ok=True)

    parts = [
        MODULE_LOADS,
        "plink2",
        f"--vcf {vcf_path}",
        f"--maf {maf}",
        f"--geno {geno}",
        f"--pca {n_pc}",
        f"--threads {threads}",
        f"--memory {memory_mb}",
        "--allow-extra-chr",
        f"--out {outprefix}",
    ]
    if thin_frac is not None:
        parts.insert(-1, f"--thin {thin_frac}")

    eigvec_path = f"{outprefix}.eigenvec"
    eigval_path = f"{outprefix}.eigenval"
    
    if not (os.path.exists(eigvec_path) and os.path.exists(eigval_path)):
        cmd = " && ".join([parts[0], " ".join(parts[1:])])
        print("🚀 Running:", cmd.replace(MODULE_LOADS + " && ", ""))
        subprocess.run(["bash", "-lc", cmd], check=True)

    if not (os.path.exists(eigvec_path) and os.path.exists(eigval_path)):
        raise FileNotFoundError("Missing PLINK output (.eigenvec/.eigenval).")
        
    return eigvec_path, eigval_path

def _read_plink_pcs(eigvec_path: str, eigval_path: str):
    """
    Lit les PCs PLINK2 en détectant le format:
      - header présent: '#IID  PC1  PC2 ...' (ton cas) -> index='IID'
      - sinon: FID IID PC1 PC2 ...
    Retourne (pcs_df, evr%) avec PC1, PC2, ...
    """
    with open(eigvec_path, "r") as f:
        first_tokens = f.readline().strip().split()

    if first_tokens and first_tokens[0].lstrip("#") == "IID":
        eigvec = pd.read_csv(eigvec_path, sep="\t", header=0, dtype=str)
        eigvec = eigvec.rename(columns={eigvec.columns[0]: "IID"})
        eigvec = eigvec.set_index("IID")
        for c in eigvec.columns:
            eigvec[c] = pd.to_numeric(eigvec[c], errors="coerce")
    else:
        eigvec = pd.read_csv(eigvec_path, sep=r"\s+", header=None, dtype={0: str, 1: str})
        ncols = eigvec.shape[1]
        eigvec.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, ncols - 1)]
        eigvec = eigvec.drop(columns=["FID"]).set_index("IID")
        for c in eigvec.columns:
            eigvec[c] = pd.to_numeric(eigvec[c], errors="coerce")

    eigval = pd.read_csv(eigval_path, sep=r"\s+", header=None)[0].values
    #evr = (eigval / eigval.sum()) * 100.0 #over estimate warning to rm
    return eigvec, eigval

# --------------------------- TE helpers ---------------------------

def _build_me_matrix(G_mei, var_ids_mei, G_mea, var_ids_mea, FIG_OUTDIR):
    """Mat TE (samples × loci): agrégation par locus + filtres."""

    outFilename = os.path.join(FIG_OUTDIR, "MEGAnE_TE_PCA.tsv")
    #if not os.path.exists(outFilename):
    X_mei = (pd.DataFrame(_dosage_te(G_mei.values, "MEI"),
                          index=G_mei.index, columns=[f"MEI_{c}" for c in G_mei.columns])
             if (G_mei is not None and not G_mei.empty) else pd.DataFrame())
    X_mea = (pd.DataFrame(_dosage_te(G_mea.values, "MEA"),
                          index=G_mea.index, columns=[f"MEA_{c}" for c in G_mea.columns])
             if (G_mea is not None and not G_mea.empty) else pd.DataFrame())
    if X_mei.empty and X_mea.empty:
        return pd.DataFrame()
    samples = sorted(set(X_mei.index) | set(X_mea.index))
    X = pd.concat([X_mei.reindex(samples), X_mea.reindex(samples)], axis=1).astype(float)
    print(X)
    X.to_csv(outFilename, sep='\t')
    #X = pd.read_csv(outFilename, sep='\t', index_col=0)
    #print(X)

    def _locus_id(ids, fb):
        if not ids: return fb
        u = sorted({str(x) for x in ids if str(x) not in {"", "nan", "None", "NaN"}})
        return "|".join(u) if u else fb

    var_mei = list(X_mei.columns); var_mea = list(X_mea.columns)
    loc_mei = [_locus_id(var_ids_mei[j], var_mei[j]) for j in range(len(var_ids_mei))] if len(var_ids_mei) else []
    loc_mea = [_locus_id(var_ids_mea[j], var_mea[j]) for j in range(len(var_ids_mea))] if len(var_ids_mea) else []
    meta = pd.DataFrame({"varname": var_mei + var_mea,
                         "locus_id": loc_mei + loc_mea}).set_index("varname").reindex(X.columns)
    X_var = X.T.groupby(meta["locus_id"]).max()

    present = (X_var >= 1).astype(float)
    nS = X_var.shape[1]
    keep = (present.sum(axis=1) >= 2) & (present.sum(axis=1) <= (nS - 1))
    X_var = X_var.loc[keep]
    return X_var.T  # samples × loci

def _prep_and_pca(M: pd.DataFrame):
    """Impute moyenne par locus, drop var=0, PCA via SVD."""
    A = M.values.copy()
    m = np.nanmean(A, axis=0, keepdims=True)
    inds = np.where(np.isnan(A))
    if inds[0].size:
        A[inds] = np.take(m.ravel(), inds[1])
    var = A.var(axis=0, ddof=1)
    keep = var > 0
    A = A[:, keep]
    if A.shape[1] < 2:
        return None
    A = A - A.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    n = A.shape[0]
    ev = (S**2) / max(n - 1, 1)
    evr = ev / ev.sum()
    scores = U * S
    return scores[:, 0], scores[:, 1], float(evr[0]) * 100.0, float(evr[1]) * 100.0

# --------------------------- Plot PCA helpers ---------------------------

def _scatter_pca(ax, x, y, ids, pops, pop_color_map, title, evr1, evr2):
    def _color_of(p): return pop_color_map.get(p, "#9e9e9e")
    colors = pops.map(_color_of).values
    ax.scatter(x, y, c=colors, s=38, edgecolor="black", linewidth=0.6, alpha=0.7)

    # texts = []
    # for i, sample in enumerate(ids):
    #     t = ax.text(x[i], y[i], sample, fontsize=5)#, alpha=0.8)
    #     texts.append(t)
    # adjust_text(
    #     texts,
    #     x=x,
    #     y=y,
    #     ax=ax,
    #     arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.6),
    #     expand_points=(2.0, 2.0),
    #     expand_text=(1.2, 1.2),
    #     force_text=1
    # )
    
    ax.set_xlabel(f"PC1 ({evr1:.1f}% var.)")
    ax.set_ylabel(f"PC2 ({evr2:.1f}% var.)")
    ax.axhline(0, color="#dddddd", lw=0.8); ax.axvline(0, color="#dddddd", lw=0.8)
    ax.set_title(title)

# --------------------------- TE & SNP (independant, side by side) ---------------------------

def run_pca_TE_and_SNP( snp_vcf, pop_color_map, VCF_MEI_PlusAdmx, VCF_MEA_PlusAdmx, POP_FILE_PlusAdmx, FIG_OUTDIR, CLASSIF_FILE, Sp ):
    class_dict, order_dict, df_pop, G_mei, var_ids_mei, G_mea, var_ids_mea, sfamily_dict = _load_all_inputs( VCF_MEI_PlusAdmx, VCF_MEA_PlusAdmx, POP_FILE_PlusAdmx, CLASSIF_FILE )
    outdir = FIG_OUTDIR

    # ---- TE PCA ----
    M_te = _build_me_matrix(G_mei, var_ids_mei, G_mea, var_ids_mea, FIG_OUTDIR)
    if M_te.empty:
        print("⚠️ [TE PCA] No informative TE loci; aborting."); return
    p_te = _prep_and_pca(M_te)
    if p_te is None:
        print("⚠️ [TE PCA] Not enough informative TE features."); return
    pc1_te, pc2_te, evr1_te, evr2_te = p_te
    samples_te = [i for i in M_te.index if i in df_pop.index]
    pops_te = df_pop["Population"].reindex(samples_te).fillna("Unknown")

    # ---- SNP PCA ----
    outprefix = os.path.join(outdir, "MEGAnE_PLINK_SNP_PCA")
    eigvec_path, eigval_path = _run_plink_pca(
        vcf_path=snp_vcf, outprefix=outprefix,
        maf=0.05, geno=0.10, thin_frac=None,
        threads=8, memory_mb=8192, n_pc=20
    )
    #print(eigvec_path, eigval_path)
    pcs, evr_snp = _read_plink_pcs(eigvec_path, eigval_path)
    #print(pcs, evr_snp)
    samples_snp = [i for i in pcs.index if i in df_pop.index]
    pops_snp = df_pop["Population"].reindex(samples_snp).fillna("Unknown")
    pc1_snp = pcs.loc[samples_snp, "PC1"].values
    pc2_snp = pcs.loc[samples_snp, "PC2"].values
    #print(pc1_snp, pc2_snp)

    # ---- Couleurs / légende ----
    all_pops = pd.Index(pops_te.tolist() + pops_snp.tolist()).unique()
    handles = [Patch(facecolor=pop_color_map.get(p, "#9e9e9e"), edgecolor="none", label=p)
               for p in sorted(all_pops)]

    # ---- Figure ----
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    _scatter_pca(axs[0], pc1_te, pc2_te, samples_te, pops_te, pop_color_map,
                 "TE", evr1_te, evr2_te)
    _scatter_pca(axs[1], pc1_snp, pc2_snp, samples_snp, pops_snp, pop_color_map,
                 "SNP", evr_snp[0], evr_snp[1])

    fig.subplots_adjust(right=0.82, wspace=0.30)
    leg = fig.legend(handles=handles, title="Population", loc="center left",
                     bbox_to_anchor=(0.819, 0.5), frameon=True, fancybox=True,
                     fontsize=9, title_fontsize=10, borderaxespad=0.6)
    leg.get_frame().set_alpha(0.90)
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("#cccccc")

    letter = ''
    if Sp == 'PCOM':
        letter = "a"
    elif Sp == "PPY":
        letter = "b"
    fig.text(0, 0.9, letter, fontsize=18, fontweight="bold", va="top", ha="left")
 
    outpath = os.path.join(outdir, "MEGAnE_PCA_TE_vs_SNP.png")
    fig.savefig(outpath, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"📈 Side-by-side PCA figure saved: {outpath}")

















def concat_imgs(img1_path, img2_path, output_path, mode="horizontal"):
    img1 = Image.open(img1_path).convert("RGB")
    img2 = Image.open(img2_path).convert("RGB")

    # 🔹 Uniformiser la largeur si mode vertical
    if mode == "vertical":
        common_width = max(img1.width, img2.width)
        img1 = img1.resize((common_width, int(img1.height * common_width / img1.width)))
        img2 = img2.resize((common_width, int(img2.height * common_width / img2.width)))
        new_img = Image.new("RGB", (common_width, img1.height + img2.height), (255, 255, 255))
        new_img.paste(img1, (0, 0))
        new_img.paste(img2, (0, img1.height))

    elif mode == "horizontal":
        common_height = max(img1.height, img2.height)
        img1 = img1.resize((int(img1.width * common_height / img1.height), common_height))
        img2 = img2.resize((int(img2.width * common_height / img2.height), common_height))
        new_img = Image.new("RGB", (img1.width + img2.width, common_height), (255, 255, 255))
        new_img.paste(img1, (0, 0))
        new_img.paste(img2, (img1.width, 0))
    else:
        raise ValueError("mode must be either 'horizontal' or 'vertical'")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    new_img.save(output_path, dpi=(300, 300))
    print(f"✅ Combined figure saved to: {output_path}")





def run_te_samplewise_stacked_barplot(VCF_MEI, VCF_MEA, FIG_OUTDIR):
    """
    Sample-wise stacked barplot of genotype proportions.
    - Raw VCF GT strings preserved (0/0, 0|1, 0/., ./., etc.)
    - MEA: alleles 0<->1 swapped, separators and '.' unchanged
    - X-axis = samples
    - Y stacked bars = proportion of each genotype class for that sample
    """

    import gzip
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import defaultdict, Counter

    custom_colors = {
        "0/0": "lightgray",
        "0/.": "black",
        "1/1": "green",
        "0/1": "pink",
        "1/0": "purple",
        "./1": "blue",
        "1/.": "lightblue"
    }

    def parse_vcf_to_sample_gt_matrix(vcf_path, invert_mea):
        """
        Returns a dict sample -> list of raw GT strings (per locus), possibly inverted.
        """
        def invert_alleles(gt):
            out = []
            for ch in gt:
                if ch == "0":
                    out.append("1")
                elif ch == "1":
                    out.append("0")
                else:
                    out.append(ch)  # keep '.', '/', '|' as they are
            return "".join(out)

        sample_gts = defaultdict(list)
        samples = []

        with gzip.open(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    samples = cols[9:]
                    continue

                parts = line.strip().split("\t")
                fmt = parts[8].split(":")
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")

                fields = parts[9:]

                for samp, field in zip(samples, fields):
                    sub = field.split(":")
                    if gt_idx >= len(sub):
                        continue
                    gt = sub[gt_idx]

                    if invert_mea:
                        gt = invert_alleles(gt)

                    sample_gts[samp].append(gt)

        return sample_gts


    # Read MEI (unchanged) and MEA (inverted)
    gts_mei = parse_vcf_to_sample_gt_matrix(VCF_MEI, invert_mea=False)
    gts_mea = parse_vcf_to_sample_gt_matrix(VCF_MEA, invert_mea=True)

    # Merge locus lists per sample
    all_samples = sorted(set(gts_mei.keys()) | set(gts_mea.keys()))
    merged = {}

    for s in all_samples:
        merged[s] = (gts_mei.get(s, []) + gts_mea.get(s, []))

    # Compute per-sample proportions
    records = []
    all_genotypes = set()

    for s in all_samples:
        c = Counter(merged[s])
        total = sum(c.values()) if c else 1
        for gt, count in c.items():
            all_genotypes.add(gt)
            records.append({
                "sample": s,
                "gt": gt,
                "prop": count / total
            })

    df = pd.DataFrame(records)

    # Pivot to sample × genotype matrix
    df_pivot = df.pivot(index="sample", columns="gt", values="prop").fillna(0)

    # Sort genotypes for visual consistency
    df_pivot = df_pivot[sorted(df_pivot.columns)]
    gt_order = ["0/1", "1/.", "./1", "1/1", "0/0", "0/."]
    existing = [gt for gt in gt_order if gt in df_pivot.columns]
    df_pivot = df_pivot[existing]

    # Plot
    df_pivot.index = df_pivot.index.map(clean_sample_name)
    os.makedirs(FIG_OUTDIR, exist_ok=True)
    outpath = os.path.join(FIG_OUTDIR, "TE_samplewise_genotype_stacked_barplot.png")

    plt.figure(figsize=(max(10, len(all_samples) * 0.4), 4))

    # bottom = None
    # for gt in df_pivot.columns:
    #     if bottom is None:
    #         plt.bar(df_pivot.index, df_pivot[gt], label=gt)
    #         bottom = df_pivot[gt].copy()
    bottom = None
    for gt in df_pivot.columns:
        color = custom_colors.get(gt, None)  # default color if not defined
        plt.bar(df_pivot.index, df_pivot[gt], bottom=bottom, label=gt, color=color)
        bottom = df_pivot[gt] if bottom is None else bottom + df_pivot[gt]
        # else:
        #     plt.bar(df_pivot.index, df_pivot[gt], bottom=bottom, label=gt)
        #     bottom += df_pivot[gt]

    plt.xticks(rotation=90, fontsize=7)
    plt.ylabel("Proportion")
    plt.title("Genotype proportions per sample (MEI + MEA, raw VCF GT)")
    plt.legend(title="GT", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.margins(x=0)
    plt.tight_layout(rect=[0, 0, 0.80, 1])  # leave 20% space on the right for legend
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


    print(f"📊 Sample-wise genotype stacked barplot saved: {outpath}")


def clean_sample_name(name):
    """
    Extract substring between 'pear.' and '.sorted'
    If pattern not found, return original name.
    """
    if "pear." in name and ".sorted" in name:
        return name.split("pear.")[1].split(".sorted")[0]
    return name


# --------------------------- MAIN ---------------------------

def main():

    # ============================== PARAMETERS ===============================
    
    L_SPECIES = ['PCOM', 'PPY']
    
    for SPECIES in L_SPECIES:
        
        VCF_MEI       = f"../output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEI_jointcall.vcf.gz"
        VCF_MEA       = f"../output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEA_jointcall.vcf.gz"
        VCF_MEI_PlusAdmx       = f"../output/{SPECIES}/jointcall_out_PlusAdmx/{SPECIES}_MEI_jointcall.vcf.gz"
        VCF_MEA_PlusAdmx       = f"../output/{SPECIES}/jointcall_out_PlusAdmx/{SPECIES}_MEA_jointcall.vcf.gz"
        CLASSIF_FILE  = f"../data/{SPECIES}/TE_annotation_URGI/{SPECIES}_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
        POP_FILE      = f"../data/{SPECIES}/s01.individuals_table.txt"
        POP_FILE_PlusAdmx      = f"../data/{SPECIES}/s01.individuals_table_PlusAdmx.txt"
        FIG_OUTDIR    = f"../Figures/{SPECIES}"
        TSV_OUTDIR    = f"../output/{SPECIES}"

        if SPECIES == 'PCOM':
            #VCF_SNP = "../data/PCOM/SNPs/pear_Jul2024_ref_comm.Europe_noAdmix.Combine_chr.geno20.small.vcf.gz" # without admix
            VCF_SNP = "../data/PCOM/SNPs/pear_Jul2024_ref_comm.Europe.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf.gz"
            
            pop_color_map = {
                "comm_Dessert": "#442CF4",
                "comm_Perry":   "#648FFF",
                "pyra":         "#035b03",
                "cauc":         "#27bc40",
                "Unknown":      "#9e9e9e",
                "admixed": "#888888"
            }
            pop_list = ["cauc", "pyra", "comm_Perry", "comm_Dessert"]
            #(Outgroup,(cauc,(comm_Dessert,(pyra,comm_Perry))));
            #newick_file = "../data/PCOM/ASTRAL.root_ussu_Western.nwk" # wrong
            newick_file = "../data/PCOM/s02.root_ussu_Western.nwk" # wrong

        elif SPECIES == 'PPY':
            #VCF_SNP = "../data/PPY/SNPs/pear_Dec2023.Asia_noAdmix.Combine_chr.geno20.small.vcf.gz" # without admix
            VCF_SNP = "../data/PPY/SNPs/pear_Dec2023.Asia.Combine_chr.geno20_maf005.anno.syno.thin8k.vcf.gz"
            
            CLASSIF_FILE = f"../data/{SPECIES}/TE_annotation_URGI/{SPECIES}3_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
            POP_FILE = f"../data/{SPECIES}/s01.individuals_table.txt3"
            pop_color_map = {
                "betu": "#9A51FF",
                "pash": "#FFA07A",
                "pyri_JP": "#FF5800",
                "Sand_CN": "#960558",
                "Sand_CN-SW": "#FFB727",
                "Sand_CN-SE": "#E599F7",
                "ussu": "#5E4734",
                "White": "#FD79A8",
                "admixed": "#888888"
            }
            pop_list = ["ussu", "pash", "betu", "Sand_CN-SE", "pyri_JP", "Sand_CN", "Sand_CN-SW", "White"]
            #(Outgroup,(ussu,((betu,pash),(Sand_CN-SE,(pyri_JP,(Sand_CN,(Sand_CN-SW,White)))))));
            newick_file = "../data/PPY/root_cauc_Eastern_tree.nwk"

        # AF binning (used by both AFS display and χ²)
        NBINS_CHI2 = 4


        # Create out directories
        os.makedirs(FIG_OUTDIR, exist_ok=True); os.makedirs(TSV_OUTDIR, exist_ok=True)

        run_te_samplewise_stacked_barplot(VCF_MEI, VCF_MEA, FIG_OUTDIR)


        ##########################################################
        # PCAs
        #run_pca_TE_and_SNP( VCF_SNP, pop_color_map, VCF_MEI_PlusAdmx, VCF_MEA_PlusAdmx, POP_FILE_PlusAdmx, FIG_OUTDIR, CLASSIF_FILE, SPECIES )
        #plink2 --vcf ../output/PPY/jointcall_out/PPY_MEA_jointcall.vcf.gz --allow-extra-chr --vcf-half-call h --out here --threads 8 --pca 2
    
        ##########################################################
        # Class and Order heatmap, fixed and polymorphic loci
        #run_me_classorder_heatmap( VCF_MEI, VCF_MEA, pop_list, newick_file, FIG_OUTDIR, SPECIES, POP_FILE, CLASSIF_FILE )


        ##########################################################
        # Build AFS+χ² summary
        #run_afs_chi2_summary()
        #plot_tajima_vs_af(VCF_MEI, SPECIES, FIG_OUTDIR) # torm
        # Write χ² TSVs
        #run_chi2_tables()
        ##########################################################
        # Clustermap
        # run_me_012_clustermap_only()
        # concat_imgs(
        #     img1_path="../Figures/PCOM/MEGAnE_ME012_clustermap.png",
        #     img2_path="../Figures/PPY/MEGAnE_ME012_clustermap.png",
        #     output_path="../Figures/MEGAnE_ME012_clustermap_combined.png",
        #     mode="horizontal"
        # )
        ##########################################################


    ### Combine the two figures
    # # PCAs
    # concat_imgs(
    #     img1_path="../Figures/PCOM/MEGAnE_PCA_TE_vs_SNP.png",
    #     img2_path="../Figures/PPY/MEGAnE_PCA_TE_vs_SNP.png",
    #     output_path="../Figures/MEGAnE_PCA_TE_vs_SNP_combined.png",
    #     mode="horizontal"
    # )
    # # Heatmaps
    # concat_imgs(
    #     img1_path="../Figures/PCOM/MEGAnE_ME_byClassOrder_heatmap.png",
    #     img2_path="../Figures/PPY/MEGAnE_ME_byClassOrder_heatmap.png",
    #     output_path="../Figures/MEGAnE_ME_byClassOrder_combined.png",
    #     mode="horizontal"
    # )
    # #concat
    # concat_imgs(
    #     img1_path="../Figures/MEGAnE_PCA_TE_vs_SNP_combined.png",
    #     img2_path="../Figures/MEGAnE_ME_byClassOrder_combined.png",
    #     output_path="../Figures/main_figure_1.png",
    #     mode="vertical"
    # )
    # concat_imgs(
    #     img1_path="../Figures/main_figure_1.png",
    #     img2_path="../Figures/MEGAnE_vs_positive_genes.png",
    #     output_path="../Figures/main_figure.png",
    #     mode="vertical"
    # )
    # concat_imgs(
    #     img1_path="../Figures/PCOM/nbMEs_depth_per_sample.png",
    #     img2_path="../Figures/PPY/nbMEs_depth_per_sample.png",
    #     output_path="../Figures/nbMEs_depth_per_sample_combined.png",
    #     mode="horizontal"
    # )    

if __name__ == "__main__":
    main()
