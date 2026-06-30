#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MEGAnE_AFS_byPopulation.png — 4 panels (2x2), AF on X (linear), Count on Y (linear)

Each panel shows the allele-frequency spectrum for ONE Population.
Bars are side-by-side within each AF bin, colored by a simplified 4-class scheme:
 - LTR
 - Class I non-LTR
 - Class II
 - Unclassified  (rendered in grey)

MEI + MEA aggregated. 40 AF bins. NaN → "NA".
Top 4 populations by number of unique samples actually present in MEI/MEA (excluding 'Unknown').
"""

import os, gzip, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec

warnings.filterwarnings("ignore", category=UserWarning)
sns.set_style("whitegrid")

NBINS = 40

# ===================== Loaders =====================

def load_classification(classif_file):
    df = pd.read_csv(classif_file, sep="\t", dtype=str)
    for col in ["Seq_name", "class", "order"]:
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {classif_file}")
    # Normalize strings and fill
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
        df.loc[df[c].isin(["", "nan", "None", "NaN"]), c] = np.nan
    df["class"]  = df["class"].fillna("Unclassified")
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
        df = df.rename(columns={"#ID_vcf":"ID_vcf"})
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
        df.loc[df[c].isin(["", "nan", "None", "NaN"]), c] = np.nan
    if "Population" not in df.columns:
        raise ValueError("Missing 'Population' column")
    df["Population"] = df["Population"].fillna("Unknown")
    return df.set_index("ID_vcf")

# ===================== VCF parsing =====================

def parse_vcf_genotypes_with_ids(vcf_file, class_dict, order_dict, sfamily_dict):
    """Return: G (samples×variants), var_ids (list of MEI IDs per variant)."""
    samples = []
    variant_genos = []
    var_ids = []

    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                hdr = line.strip().split("\t")
                samples = hdr[9:]
                continue

            cols = line.strip().split("\t")
            if len(cols) < 10: continue
            info = cols[7]; fmt = cols[8]; genotypes = cols[9:]

            mei_ids = []
            for kv in info.split(";"):
                if kv.startswith("MEI="):
                    mei_ids = kv.split("=")[1].split("|")
                    break
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
                al = gt.replace("|","/").split("/")
                vals = [int(a) for a in al if a!="."]
                row.append(sum(vals) if len(vals)==2 else np.nan)
            variant_genos.append(row)

    idx = [s.split(".")[2] if len(s.split("."))>=3 else s for s in samples]
    if not variant_genos:
        return pd.DataFrame(index=idx), var_ids

    G = pd.DataFrame(variant_genos).T
    G.columns = [f"var{i}" for i in range(G.shape[1])]
    G.index = idx
    return G, var_ids

# ===================== AF & presence helpers =====================

def _dosage_te(G_values, type_name):
    """Return dosage of TE alleles per genotype matrix:
       - MEI: dosage = genotype (0/1/2)
       - MEA: dosage = 2 - genotype (since 0 = TE presence in MEA VCFs)
    """
    g = np.array(G_values, dtype=float)
    return g if type_name == "MEI" else (2.0 - g)

def af_per_variant_for_mask(G, type_name, sample_mask):
    """Allele frequency per variant (within masked samples)."""
    if G.empty or not np.any(sample_mask):
        return np.array([])
    d = _dosage_te(G.values, type_name)
    d = d[np.asarray(sample_mask), :]
    mask = ~np.isnan(d)
    n_called = mask.sum(axis=0)
    d[~mask] = 0.0
    te_sum = d.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        af = np.where(n_called > 0, te_sum / (2.0 * n_called), np.nan)
    af = af[~np.isnan(af)]
    return np.clip(af, 0.0, 1.0)

def variant_indices_present_in_mask(G, type_name, sample_mask):
    """Indices of variants present (dosage TE > 0 in >=1 masked sample)."""
    if G.empty or not np.any(sample_mask):
        return np.array([], dtype=int)
    d = _dosage_te(G.values, type_name)
    d = d[np.asarray(sample_mask), :]
    present = (d > 0).any(axis=0)
    return np.where(present)[0]

def unique_sample_count_in_pop(pop, G_mei, G_mea, df_pop):
    """Unique sample count in 'pop' present in MEI and/or MEA matrices."""
    pop_ids = set(df_pop.index[df_pop["Population"] == pop].astype(str))
    mei_ids = set(G_mei.index.astype(str)) if not G_mei.empty else set()
    mea_ids = set(G_mea.index.astype(str)) if not G_mea.empty else set()
    return len((mei_ids | mea_ids) & pop_ids)

# ===================== Simplified labeling (4 categories) =====================

def label_for_variant(ids, class_dict, order_dict, sfamily_dict):
    """Assign one simplified label per variant: LTR, Class I non-LTR, Class II, Unclassified."""
    if not ids:
        return "Unclassified"

    def get_cls(i):  return str(class_dict.get(i, "Unclassified")).split("|")[0]
    def get_odr(i):  return str(order_dict.get(i, "NA")).split("|")[0]

    orders  = [get_odr(i) for i in ids]
    classes = [get_cls(i) for i in ids]

    # LTR takes precedence
    if any(odr == "LTR" for odr in orders):
        return "LTR"

    # Class II
    if any(cls == "II" for cls in classes):
        return "Class II"

    # Class I non-LTR
    if any(cls == "I" for cls in classes):
        return "Class I non-LTR"

    # Unclassified
    return "Unclassified"

def build_labels_array(var_ids_list, class_dict, order_dict, sfamily_dict):
    return np.array([label_for_variant(ids, class_dict, order_dict, sfamily_dict)
                     for ids in var_ids_list], dtype=object)

def make_label_palette(labels):
    """Fixed palette for the four simplified categories; Unclassified forced to grey."""
    base_palette = {
        "LTR": (0.20, 0.60, 0.80),           # light blue
        "Class I non-LTR": (0.90, 0.60, 0.00), # orange
        "Class II": (0.30, 0.70, 0.30),      # green
        "Unclassified": (0.60, 0.60, 0.60),  # grey
    }
    # Keep only labels that occur; preserve canonical order
    order = ["LTR", "Class I non-LTR", "Class II", "Unclassified"]
    labels_present = [l for l in order if l in set(labels)]
    palette = {lab: base_palette[lab] for lab in labels_present}
    return palette, labels_present

# ===================== Plot helpers =====================

def counts_per_bin_label(af_vals, lab_vals, bins, label_order):
    """Return (nbins × nlabels) counts matrix."""
    if af_vals.size == 0 or (lab_vals is None) or (len(lab_vals) == 0):
        return np.zeros((len(bins) - 1, len(label_order)), dtype=int)
    bin_idx = np.digitize(af_vals, bins, right=False) - 1
    bin_idx = np.clip(bin_idx, 0, len(bins) - 2)
    dfc = pd.DataFrame({"bin": bin_idx, "label": lab_vals})
    counts = dfc.value_counts().reset_index(name="count")
    mat = counts.pivot(index="bin", columns="label", values="count") \
               .reindex(columns=label_order).fillna(0).astype(int)
    all_bins = pd.Index(range(len(bins) - 1), name="bin")
    mat = mat.reindex(index=all_bins, fill_value=0)
    return mat.values

def barplot_grouped_linear(ax, bins, counts_mat, labels, palette):
    """Side-by-side bars within each AF bin; X linear; Y linear."""
    nbins = len(bins) - 1
    nlab = len(labels)
    if nlab == 0:
        y = counts_mat.sum(axis=1)
        lefts = bins[:-1] + 0.05 * (bins[1:] - bins[:-1])
        widths = 0.90 * (bins[1:] - bins[:-1])
        ax.bar(lefts, y, width=widths, align='edge', color="grey",
               edgecolor="black", linewidth=0.3)
        return
    for b in range(nbins):
        left = bins[b]; bw = bins[b + 1] - bins[b]
        inner_left = left + 0.05 * bw
        group_w = 0.90 * bw
        bar_w = group_w / max(nlab, 1)
        for i, lab in enumerate(labels):
            x = inner_left + i * bar_w
            y = counts_mat[b, i]
            ax.bar(x, y, width=bar_w, align='edge',
                   color=palette.get(lab, (0.7, 0.7, 0.7)),
                   edgecolor="black", linewidth=0.3, label=lab if b == 0 else None)

def grouped_afs(ax, af_vals, labels=None, label_palette=None,
                show_xlabel=False, show_ylabel=False, title=None):
    bins = np.linspace(0.0, 1.0, NBINS + 1)
    ax.set_xlim(0.0, 1.0)

    if labels is None or (label_palette is None):
        af_np = np.asarray(af_vals, dtype=float)
        af_np = af_np[~np.isnan(af_np)]
        if af_np.size == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center")
            ax.set_axis_off()
            return
        af_np = np.clip(af_np, 0.0, 1.0)
        counts, _ = np.histogram(af_np, bins=bins)
        lefts = bins[:-1] + 0.05 * (bins[1:] - bins[:-1])
        widths = 0.90 * (bins[1:] - bins[:-1])
        ax.bar(lefts, counts, width=widths, align='edge', color="grey",
               edgecolor="black", linewidth=0.3)
    else:
        af_np = np.asarray(af_vals, dtype=float)
        lab_arr = np.array(labels, dtype=object)
        n = min(len(af_np), len(lab_arr))
        if n == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center")
            ax.set_axis_off()
            return
        af_np = af_np[:n]; lab_arr = lab_arr[:n]
        valid = ~np.isnan(af_np)
        af_np = np.clip(af_np[valid], 0.0, 1.0)
        lab_arr = lab_arr[valid]
        # enforce 4-class vocabulary only
        lab_arr = np.array([x if x in {"LTR","Class I non-LTR","Class II","Unclassified"} else "Unclassified"
                            for x in lab_arr], dtype=object)
        label_order = list(label_palette.keys())
        counts_mat = counts_per_bin_label(af_np, lab_arr, bins, label_order)
        barplot_grouped_linear(ax, bins, counts_mat, label_order, label_palette)

    if show_xlabel:
        ax.set_xlabel("Allele frequency (TE)", fontsize=10)
        ax.tick_params(axis='x', labelsize=9)
    else:
        ax.set_xlabel(""); ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("Count", fontsize=10)
        ax.tick_params(axis='y', labelsize=9)
    else:
        ax.set_ylabel(""); ax.set_yticklabels([])
    if title:
        ax.set_title(title, fontsize=12, pad=8)

# ===================== MAIN =====================

def main():
    # Paths
    vcf_mei = "../output/PCOM/jointcall_out/PCOM_MEI_jointcall.vcf.gz"
    vcf_mea = "../output/PCOM/jointcall_out/PCOM_MEA_jointcall.vcf.gz"
    classif_file = "../data/PCOM/TE_annotation_URGI/PCOM_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
    pop_file = "../data/PCOM/s01.individuals_table.txt"
    outdir = "/scratch/ss20440/MEGAnE_Yuqi/Figures"
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "MEGAnE_AFS_byPopulation.png")

    # Load metadata
    _, class_dict, order_dict, sfamily_dict = load_classification(classif_file)
    df_pop = load_population_table(pop_file)

    # Read VCFs
    G_mei, var_ids_mei = parse_vcf_genotypes_with_ids(vcf_mei, class_dict, order_dict, sfamily_dict)
    G_mea, var_ids_mea = parse_vcf_genotypes_with_ids(vcf_mea, class_dict, order_dict, sfamily_dict)

    # Build simplified labels per variant (independent of population)
    labels_mei = build_labels_array(var_ids_mei, class_dict, order_dict, sfamily_dict) if len(var_ids_mei) > 0 else np.array([], object)
    labels_mea = build_labels_array(var_ids_mea, class_dict, order_dict, sfamily_dict) if len(var_ids_mea) > 0 else np.array([], object)

    # Choose top-4 populations by number of unique samples actually present in MEI/MEA
    present_counts = {}
    for pop in df_pop["Population"].unique():
        if pop == "Unknown":
            continue
        present_counts[pop] = unique_sample_count_in_pop(pop, G_mei, G_mea, df_pop)
    top4 = [p for p, _ in sorted(present_counts.items(), key=lambda kv: kv[1], reverse=True)[:4]]
    if len(top4) < 4 and "Unknown" in df_pop["Population"].unique():
        top4 += ["Unknown"] * (4 - len(top4))

    # Prepare figure (2×2)
    fig = plt.figure(figsize=(22, 16))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.38, hspace=0.38)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]

    for ax, pop in zip(axes, top4):
        # Sample masks for this population (separately for MEI and MEA matrices)
        mask_mei = (df_pop["Population"].reindex(G_mei.index.astype(str)).fillna("Unknown").values == pop) if not G_mei.empty else np.array([], bool)
        mask_mea = (df_pop["Population"].reindex(G_mea.index.astype(str)).fillna("Unknown").values == pop) if not G_mea.empty else np.array([], bool)

        # Variant indices present in this population
        idx_mei_pop = variant_indices_present_in_mask(G_mei, "MEI", mask_mei) if not G_mei.empty else np.array([], dtype=int)
        idx_mea_pop = variant_indices_present_in_mask(G_mea, "MEA", mask_mea) if not G_mea.empty else np.array([], dtype=int)

        # AF arrays computed only on variants present in this population
        af_mei = af_per_variant_for_mask(G_mei.loc[:, G_mei.columns[idx_mei_pop]], "MEI", mask_mei) if idx_mei_pop.size else np.array([])
        af_mea = af_per_variant_for_mask(G_mea.loc[:, G_mea.columns[idx_mea_pop]], "MEA", mask_mea) if idx_mea_pop.size else np.array([])
        af_all = np.concatenate([af_mei, af_mea])

        # Labels aligned to selected variants
        labs_mei_sel = labels_mei[idx_mei_pop] if idx_mei_pop.size else np.array([], object)
        labs_mea_sel = labels_mea[idx_mea_pop] if idx_mea_pop.size else np.array([], object)
        labs_all = np.concatenate([labs_mei_sel, labs_mea_sel]) if (labs_mei_sel.size or labs_mea_sel.size) else None

        # Fixed palette for simplified classes (Unclassified = grey)
        label_palette, label_order = make_label_palette(labs_all if labs_all is not None else [])

        # Title with correct unique sample count + number of variants present
        n_mei_present = int(idx_mei_pop.size)
        n_mea_present = int(idx_mea_pop.size)
        n_samples = unique_sample_count_in_pop(pop, G_mei, G_mea, df_pop)
        title = f"{pop} (N={n_samples}; MEI={n_mei_present}, MEA={n_mea_present})"

        # Plot
        grouped_afs(ax, af_all, labels=labs_all, label_palette=label_palette,
                    show_xlabel=True, show_ylabel=True, title=title)

        # Legend (top-right, vertical)
        if label_palette:
            handles = [plt.Rectangle((0, 0), 1, 1, color=col, ec="black", lw=0.3)
                       for _, col in label_palette.items()]
            labels  = list(label_palette.keys())
            ax.legend(handles, labels, fontsize=9, frameon=True,
                      loc="upper left", bbox_to_anchor=(1.02, 1.00),
                      title="Category", ncol=1, borderaxespad=0.0)

    # Hide unused axes if fewer than 4 selected
    for ax in axes[len(top4):]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(outfile, dpi=350, bbox_inches="tight")
    plt.close()
    print(f"✅ Figure saved: {outfile}")

if __name__ == "__main__":
    main()
