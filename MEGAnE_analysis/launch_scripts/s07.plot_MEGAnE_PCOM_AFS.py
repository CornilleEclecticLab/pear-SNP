#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MEGAnE_AFS.png — 8-panel grid (4x2), AF on X (linear), Count on Y (log)

Panels (left→right, top→bottom):
a) LTR (bars by sFamily)
b) LINE + SINE (class I, orders LINE/SINE only)
c) DIRS                    -> class=='I' & order contains 'DIRS'
d) Class I (non-LTR, other)-> class=='I' & order ∉ {LTR, LINE, SINE, DIRS}
e) TIR                     -> class=='II' & order contains 'TIR'
f) MITE                    -> class=='II' & order contains 'MITE'
g) Helitron                -> class=='II' & order contains 'Helitron'
h) Other Class II + Unclassified -> class in {II, Unclassified} and (if II) order ∉ {TIR, Helitron, MITE}

Each panel: 3 rows (Wild-only, Cultivated-only, Both). Bars are side-by-side per label within each AF bin.
MEI + MEA aggregated. 40 AF bins (linear). NaN → "NA".
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
    if "Crop_or_Wild" not in df.columns:
        raise ValueError("Missing 'Crop_or_Wild' column")
    df["Crop_or_Wild"] = df["Crop_or_Wild"].fillna("Unknown").replace({"NA":"Unknown"})
    for col in ("Depth","ReadsAvgLength"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.set_index("ID_vcf")

# ===================== VCF parsing =====================

def parse_vcf_genotypes_with_taxa(vcf_file, class_dict, order_dict, sfamily_dict):
    samples = []
    variant_genos = []
    var_orders, var_classes, var_ids = [], [], []

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
            if not mei_ids: mei_ids = []

            if mei_ids:
                ords = set(str(order_dict.get(i, "NA")).split("|")[0] for i in mei_ids)
                clss = set(str(class_dict.get(i, "Unclassified")).split("|")[0] for i in mei_ids)
            else:
                ords = set(["NA"]); clss = set(["Unclassified"])
            var_orders.append(ords)
            var_classes.append(clss)
            var_ids.append(mei_ids)

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
        return pd.DataFrame(index=idx), var_orders, var_classes, var_ids

    G = pd.DataFrame(variant_genos).T
    G.columns = [f"var{i}" for i in range(G.shape[1])]
    G.index = idx
    return G, var_orders, var_classes, var_ids

# ===================== AF helpers =====================

def _dosage_te(G_values, type_name):
    g = np.array(G_values, dtype=float)
    return g if type_name == "MEI" else (2.0 - g)

def _af_all_variants_for_mask(G, type_name, sample_mask):
    """AF per variant for all variants, under a sample mask (X in [0,1], linear)."""
    if G.empty or not np.any(sample_mask): 
        return np.array([])
    d = _dosage_te(G.values, type_name)
    d = d[np.asarray(sample_mask), :]
    mask = ~np.isnan(d)
    n_called = mask.sum(axis=0)
    d[~mask] = 0.0
    te_sum = d.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        af = np.where(n_called>0, te_sum/(2.0*n_called), np.nan)
    af = af[~np.isnan(af)]
    af = np.clip(af, 0.0, 1.0)
    return af

# ===================== Crop groups =====================

def crop_masks_for_G(G, df_pop):
    """Return sample masks aligned to G.index."""
    if G.empty:
        return np.array([], bool), np.array([], bool), np.array([], bool)
    meta = df_pop["Crop_or_Wild"].reindex(G.index.astype(str)).fillna("Unknown")
    wild_mask = (meta.values == "Wild")
    cult_mask = (meta.values == "Cultivated")
    both_mask = wild_mask | cult_mask
    return wild_mask, cult_mask, both_mask

def crop_variant_index_sets(G, df_pop, type_name):
    """Indices of variants: Wild_only, Cultivated_only, Both (presence TE>0)."""
    if G.empty:
        return {'Wild_only': [], 'Cultivated_only': [], 'Both': []}
    wild_mask, cult_mask, _ = crop_masks_for_G(G, df_pop)
    d = _dosage_te(G.values, type_name)
    present_w = (d[wild_mask, :] > 0).any(axis=0) if wild_mask.any() else np.zeros(G.shape[1], bool)
    present_c = (d[cult_mask, :] > 0).any(axis=0) if cult_mask.any() else np.zeros(G.shape[1], bool)
    return {
        'Wild_only': np.where(present_w & ~present_c)[0].tolist(),
        'Cultivated_only': np.where(~present_w & present_c)[0].tolist(),
        'Both': np.where(present_w & present_c)[0].tolist(),
    }

# ===================== Category filtering by IDs =====================

def filter_ids_for_category(ids, category, class_dict, order_dict):
    """Keep only IDs that truly belong to the requested category."""
    if not ids:
        return []
    sel = []
    for i in ids:
        cls = str(class_dict.get(i, "Unclassified")).split("|")[0]
        odr = str(order_dict.get(i, "NA")).split("|")[0]
        if category == "LTR":
            cond = (odr == "LTR")
        elif category == "LINE_SINE":
            cond = (cls == "I") and (odr in {"LINE","SINE"})
        elif category == "DIRS":
            cond = (cls == "I") and ("DIRS" in odr)
        elif category == "ClassI_other_nonLTR":
            cond = (cls == "I") and (odr not in {"LTR","LINE","SINE","DIRS"})
        elif category == "TIR":
            cond = (cls == "II") and ("TIR" in odr)
        elif category == "MITE":
            cond = (cls == "II") and ("MITE" in odr)
        elif category == "Helitron":
            cond = (cls == "II") and ("Helitron" in odr)
        elif category == "ClassII_Unclassified_other":
            # Class II but not TIR/Helitron/MITE, or Unclassified
            cond = ((cls == "II") and (odr not in {"TIR","Helitron","MITE"})) or (cls == "Unclassified")
        else:
            cond = False
        if cond:
            sel.append(i)
    return sel

def mask_category_from_ids(var_ids_list, category, class_dict, order_dict):
    """True if the variant has at least one ID in the category."""
    return np.array([len(filter_ids_for_category(ids, category, class_dict, order_dict)) > 0
                     for ids in var_ids_list], dtype=bool)

# ===================== Labels =====================

def _mode_or_first(values):
    if not values: return "NA"
    s = pd.Series(list(values)).replace({np.nan:"NA"}).replace({"":"NA","nan":"NA","None":"NA","NaN":"NA"})
    return s.value_counts().index[0]

def build_variant_label(ids, class_dict, order_dict, sfamily_dict, category):
    """Category-consistent label (IDs filtered before aggregation)."""
    sel = filter_ids_for_category(ids, category, class_dict, order_dict)
    if not sel:
        return "NA"
    if category == "LTR":
        vals = [str(sfamily_dict.get(i, "NA")).split("|")[0] for i in sel]
        vals = ["NA" if v in ["","nan","None","NaN"] else v for v in vals]
        return _mode_or_first(vals)
    elif category == "ClassII_Unclassified_other":
        labeled = []
        for i in sel:
            cls = str(class_dict.get(i, "Unclassified")).split("|")[0]
            if cls == "II":
                labeled.append("Class II (other)")
            else:
                labeled.append("Unclassified")
        return _mode_or_first(labeled)
    else:
        vals = [str(order_dict.get(i, "NA")).split("|")[0] for i in sel]
        vals = ["NA" if v in ["","nan","None","NaN"] else v for v in vals]
        return _mode_or_first(vals)

def make_label_palette(labels):
    labels = [("NA" if (x in [None, "", "nan", "None", "NaN"]) else x) for x in labels]
    labels = list(dict.fromkeys(labels))
    if not labels: return {}, []
    n = len(labels)
    pal = sns.color_palette("tab20", n_colors=min(20, n))
    if n > 20:
        pal += sns.color_palette("husl", n_colors=n-20)
    return {lab: pal[i] for i,lab in enumerate(labels)}, labels

# ===================== Plot (grouped bars; X linear; Y log) =====================

def _counts_per_bin_label(af_vals, lab_vals, bins, label_order):
    if af_vals.size == 0 or (lab_vals is None) or (len(lab_vals)==0):
        return np.zeros((len(bins)-1, len(label_order)), dtype=int)
    bin_idx = np.digitize(af_vals, bins, right=False) - 1
    bin_idx = np.clip(bin_idx, 0, len(bins)-2)
    dfc = pd.DataFrame({"bin": bin_idx, "label": lab_vals})
    counts = dfc.value_counts().reset_index(name="count")
    mat = counts.pivot(index="bin", columns="label", values="count").reindex(columns=label_order).fillna(0).astype(int)
    all_bins = pd.Index(range(len(bins)-1), name="bin")
    mat = mat.reindex(index=all_bins, fill_value=0)
    return mat.values

def _barplot_grouped_linearX_logY(ax, bins, counts_mat, labels, palette):
    """Side-by-side bars within each bin; X linear; Y log."""
    nbins = len(bins)-1
    nlab = len(labels)
    if nlab == 0:
        y = counts_mat.sum(axis=1)
        lefts = bins[:-1] + 0.05*(bins[1:]-bins[:-1])
        widths = 0.9*(bins[1:]-bins[:-1])
        ax.bar(lefts, y, width=widths, align='edge', color="grey",
               edgecolor="black", linewidth=0.3, log=False)
        return
    for b in range(nbins):
        left = bins[b]; bw = bins[b+1]-bins[b]
        inner_left = left + 0.05*bw
        group_w = 0.90*bw
        bar_w = group_w / nlab
        for i, lab in enumerate(labels):
            x = inner_left + i*bar_w
            y = counts_mat[b, i]
            ax.bar(x, y, width=bar_w, align='edge',
                   color=palette.get(lab, (0.7,0.7,0.7)),
                   edgecolor="black", linewidth=0.3, label=lab if b==0 else None, log=False)

def grouped_afs(ax, af_vals, labels=None, label_palette=None, title=None,
                show_xlabel=False, show_ylabel=False, simple=False):
    bins = np.linspace(0.0, 1.0, NBINS+1)
    ax.set_xlim(0.0, 1.0)

    if simple or labels is None or (label_palette is None):
        af_np = np.asarray(af_vals, dtype=float)
        af_np = af_np[~np.isnan(af_np)]
        if af_np.size == 0:
            ax.text(0.5,0.5,"no data",ha="center",va="center"); ax.set_axis_off(); return
        af_np = np.clip(af_np, 0.0, 1.0)
        counts, _ = np.histogram(af_np, bins=bins)
        lefts = bins[:-1] + 0.05*(bins[1:]-bins[:-1])
        widths = 0.9*(bins[1:]-bins[:-1])
        ax.bar(lefts, counts, width=widths, align='edge', color="grey",
               edgecolor="black", linewidth=0.3, log=False)
    else:
        af_np = np.asarray(af_vals, dtype=float)
        lab_arr = np.array(labels, dtype=object)
        n = min(len(af_np), len(lab_arr))
        if n == 0:
            ax.text(0.5,0.5,"no data",ha="center",va="center"); ax.set_axis_off(); return
        af_np = af_np[:n]; lab_arr = lab_arr[:n]
        valid = ~np.isnan(af_np)
        af_np = np.clip(af_np[valid], 0.0, 1.0)
        lab_arr = lab_arr[valid]
        lab_arr = np.array(["NA" if x in [None,"","nan","None","NaN"] else x for x in lab_arr], dtype=object)
        label_order = list(label_palette.keys())
        counts_mat = _counts_per_bin_label(af_np, lab_arr, bins, label_order)
        _barplot_grouped_linearX_logY(ax, bins, counts_mat, label_order, label_palette)

    if show_xlabel:
        ax.set_xlabel("Allele frequency (TE)", fontsize=9)
        ax.tick_params(axis='x', labelsize=9)
    else:
        ax.set_xlabel(""); ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("Count", fontsize=9)
        ax.tick_params(axis='y', labelsize=9)
    else:
        ax.set_ylabel(""); ax.set_yticklabels([])
    if title: ax.set_title(title, fontsize=10, pad=8)

# ===================== MAIN =====================

def main():
    # Paths
    vcf_mei = "../output/PCOM/jointcall_out/PCOM_MEI_jointcall.vcf.gz"
    vcf_mea = "../output/PCOM/jointcall_out/PCOM_MEA_jointcall.vcf.gz"
    classif_file = "../data/PCOM/TE_annotation_URGI/PCOM_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
    pop_file = "../data/PCOM/s01.individuals_table.txt"
    outdir = "/scratch/ss20440/MEGAnE_Yuqi/Figures"
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "MEGAnE_AFS.png")

    # Load
    _, class_dict, order_dict, sfamily_dict = load_classification(classif_file)
    df_pop = load_population_table(pop_file)

    # Parse
    G_mei, var_orders_mei, var_classes_mei, var_ids_mei = parse_vcf_genotypes_with_taxa(
        vcf_mei, class_dict, order_dict, sfamily_dict
    )
    G_mea, var_orders_mea, var_classes_mea, var_ids_mea = parse_vcf_genotypes_with_taxa(
        vcf_mea, class_dict, order_dict, sfamily_dict
    )

    # Sample masks
    wmask_mei, cmask_mei, bmask_mei = crop_masks_for_G(G_mei, df_pop)
    wmask_mea, cmask_mea, bmask_mea = crop_masks_for_G(G_mea, df_pop)
    n_wild = int(wmask_mei.sum() or wmask_mea.sum())
    n_cult = int(cmask_mei.sum() or cmask_mea.sum())
    n_both = int((bmask_mei | bmask_mea).sum())

    # Variant indices per crop subset (presence-based)
    crop_idx_mei = crop_variant_index_sets(G_mei, df_pop, "MEI")
    crop_idx_mea = crop_variant_index_sets(G_mea, df_pop, "MEA")

    # Panel order (includes DIRS panel, and "Other Class II + Unclassified" at the end)
    cats = [
        ("a","LTR"),
        ("b","LINE_SINE"),
        ("c","DIRS"),
        ("d","ClassI_other_nonLTR"),
        ("e","TIR"),
        ("f","MITE"),
        ("g","Helitron"),
        ("h","ClassII_Unclassified_other"),
    ]

    # Variant masks per category (ID-filtered)
    masks = {}
    for _, cat in cats:
        masks[cat] = (
            mask_category_from_ids(var_ids_mei, cat, class_dict, order_dict),
            mask_category_from_ids(var_ids_mea, cat, class_dict, order_dict)
        )

    # Precompute labels per category
    def precompute_labels_for_category(cat, var_ids):
        return np.array([build_variant_label(ids, class_dict, order_dict, sfamily_dict, cat)
                         for ids in var_ids], dtype=object)
    labels_mei = {cat: precompute_labels_for_category(cat, var_ids_mei) for _, cat in cats}
    labels_mea = {cat: precompute_labels_for_category(cat, var_ids_mea) for _, cat in cats}

    # Figure 4x2 with larger spacing and vertical legends at top
    fig = plt.figure(figsize=(34, 18))
    gs_main = gridspec.GridSpec(2, 4, figure=fig, wspace=0.40, hspace=0.36)
    panel_positions = [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)]
    title_map = {
        "LTR":"LTR",
        "LINE_SINE":"LINE + SINE",
        "DIRS":"DIRS",
        "ClassI_other_nonLTR":"Class I (non-LTR, other)",
        "TIR":"TIR",
        "MITE":"MITE",
        "Helitron":"Helitron",
        "ClassII_Unclassified_other":"Other Class II + Unclassified",
    }

    for (letter, cat), pos in zip(cats, panel_positions):
        mask_mei, mask_mea = masks[cat]

        # 3×1 grid per panel
        gs_sub = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_main[pos], wspace=0.0, hspace=0.28)
        axes = [fig.add_subplot(gs_sub[i,0]) for i in range(3)]
        subset_defs = [("Wild", "Wild_only", n_wild),
                       ("Cultivated", "Cultivated_only", n_cult),
                       ("Both", "Both", n_both)]

        # Collect labels for legend across all three subsets
        labels_set = set()
        labs_mei_cat_all = labels_mei[cat]
        for key in ["Wild_only","Cultivated_only","Both"]:
            idx = np.intersect1d(np.where(mask_mei)[0], np.array(crop_idx_mei.get(key, []), dtype=int), assume_unique=False)
            if idx.size:
                labels_set.update(labs_mei_cat_all[idx])
        labs_mea_cat_all = labels_mea[cat]
        for key in ["Wild_only","Cultivated_only","Both"]:
            idx = np.intersect1d(np.where(mask_mea)[0], np.array(crop_idx_mea.get(key, []), dtype=int), assume_unique=False)
            if idx.size:
                labels_set.update(labs_mea_cat_all[idx])

        label_palette, label_order = make_label_palette(sorted(labels_set))

        # Draw the three rows
        for k, (ax, (subset_name, subset_key, n_samples)) in enumerate(zip(axes, subset_defs)):
            show_ylabel = True
            show_xlabel = (k == 2)

            idx_mei = np.intersect1d(np.where(mask_mei)[0], np.array(crop_idx_mei.get(subset_key, []), dtype=int), assume_unique=False)
            idx_mea = np.intersect1d(np.where(mask_mea)[0], np.array(crop_idx_mea.get(subset_key, []), dtype=int), assume_unique=False)

            # choose sample masks
            if subset_key == "Wild_only":
                sm_mei, sm_mea = wmask_mei, wmask_mea
            elif subset_key == "Cultivated_only":
                sm_mei, sm_mea = cmask_mei, cmask_mea
            else:
                sm_mei, sm_mea = (wmask_mei | cmask_mei), (wmask_mea | cmask_mea)

            # AF arrays (recompute safely on subset columns)
            def _af_select(G, typ, idx, smask):
                if G.empty or idx.size==0 or not np.any(smask): return np.array([])
                sub = G.loc[:, G.columns[idx]]
                return _af_all_variants_for_mask(sub, typ, smask)

            af_mei_sel = _af_select(G_mei, "MEI", idx_mei, sm_mei)
            af_mea_sel = _af_select(G_mea, "MEA", idx_mea, sm_mea)
            af_all = np.concatenate([af_mei_sel, af_mea_sel])

            # Subplot title with sample and variant counts
            title_this = f"{subset_name} (N={n_samples}; MEI={len(idx_mei)}, MEA={len(idx_mea)})"

            # Labels for this subset
            labs_mei_sel = labels_mei[cat][idx_mei] if idx_mei.size else np.array([], object)
            labs_mea_sel = labels_mea[cat][idx_mea] if idx_mea.size else np.array([], object)
            labs_all = np.concatenate([labs_mei_sel, labs_mea_sel]) if (labs_mei_sel.size or labs_mea_sel.size) else None

            if labs_all is None:
                grouped_afs(ax, af_all, labels=None, label_palette=None,
                            title=title_this, simple=True,
                            show_xlabel=show_xlabel, show_ylabel=show_ylabel)
            else:
                grouped_afs(ax, af_all, labels=labs_all, label_palette=label_palette,
                            title=title_this, simple=False,
                            show_xlabel=show_xlabel, show_ylabel=show_ylabel)

        # Legend: top-right and vertical (one column)
        if label_palette:
            handles = [plt.Rectangle((0,0),1,1, color=col, ec="black", lw=0.3)
                       for _, col in label_palette.items()]
            labels  = list(label_palette.keys())
            axes[0].legend(
                handles, labels, fontsize=8, frameon=True,
                loc="upper left", bbox_to_anchor=(1.02, 1.00),
                title="Labels", ncol=1, borderaxespad=0.0
            )

        # Panel header
        axes[0].text(0, 1.18, f"{letter}", transform=axes[0].transAxes,
                     fontsize=14, fontweight="bold")
        axes[0].text(0.10, 1.18, title_map.get(cat, cat),
                     transform=axes[0].transAxes, fontsize=12, fontweight="bold")

    plt.tight_layout()
    plt.savefig(outfile, dpi=350, bbox_inches="tight")
    plt.close()
    print(f"✅ Figure saved: {outfile}")

if __name__ == "__main__":
    main()
