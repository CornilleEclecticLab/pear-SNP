# -*- coding: utf-8 -*-
"""
REPET summary plots:
(a) TE composition, including "No TE coverage"
(b) TE length distribution (violin plot)
(c) TE sequence identity distribution (violin plot)
(d) Six-panel (2 × 3) plot showing TE density for each chromosome (100-kbp windows; solid black line)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from pathlib import Path

# -----------------------------
# Paths & constants
# -----------------------------
FIG_DIR = Path("../Figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

FAI_FILE = "../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat.fai"
GFF_FILE = "../data/PCOM/TE_annotation_URGI/PCOM_TE_annotation.gff3"
CLASSIF_FILE = "../data/PCOM/TE_annotation_URGI/PCOM_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"

WINDOW_SIZE = 100_000  # 100 kbp
OUTPUT_FILE = FIG_DIR / f"TEdensity_{WINDOW_SIZE//1000}kb.tsv"

# -----------------------------
# Helpers
# -----------------------------
def read_fai(fai_path: Path) -> pd.DataFrame:
    fai = pd.read_csv(
        fai_path, sep="\t", header=None,
        names=["chromosome", "length", "offset", "linebases", "linewidth"]
    )
    return fai[["chromosome", "length"]]

def extract_target(attr: str):
    m = re.search(r"Target=([^\s;]+)", attr or "")
    return m.group(1) if m else None

def extract_identity(attr: str):
    m = re.search(r"AlignIdentity=([\d\.]+)", attr or "")
    return float(m.group(1)) if m else np.nan

def keep_before_pipe(x):
    if isinstance(x, str):
        return x.split("|", 1)[0].strip()
    return x

def compute_local_density(df: pd.DataFrame, chrom, window_size=100_000, feature="TE",
                          chrom_length: int | None = None) -> pd.DataFrame:
    """Local density was calculated as the number of base pairs covered within each window."""
    sub = df[df["chromosome"] == chrom]
    if sub.empty and chrom_length is None:
        return pd.DataFrame()
    ch_len = chrom_length if chrom_length is not None else (sub["end"].max() if not sub.empty else 0)
    if ch_len == 0:
        return pd.DataFrame()

    bins = np.arange(0, ch_len + window_size, window_size)
    starts = sub["start"].values if not sub.empty else np.array([], dtype=int)
    ends = sub["end"].values if not sub.empty else np.array([], dtype=int)

    densities = []
    for start in bins[:-1]:
        end = start + window_size
        if starts.size:
            overlap_bp = np.maximum(0, np.minimum(ends, end) - np.maximum(starts, start) + 1).sum()
        else:
            overlap_bp = 0
        densities.append({
            "chromosome": chrom,
            "window_start_bp": start,
            f"{feature}_density_bp": int(overlap_bp),
        })
    return pd.DataFrame(densities)

def compute_and_store_te_density(gff: pd.DataFrame, chrom_order, chrom_sizes: pd.DataFrame,
                                 window_size: int, output_file: Path) -> pd.DataFrame:
    """Pre-compute (or reload) TE density only."""
    if output_file.exists():
        print(f"[INFO] Loading pre-computed TE density data from {output_file}")
        return pd.read_csv(output_file, sep="\t")

    print(f"[INFO] Computing TE densities (window_size {window_size/1000:.0f} kb)...")
    print(f"[DEBUG] Expected chromosome names : {chrom_order}")
    print(f"[DEBUG] Chromosomes present in the GFF : {sorted(gff['chromosome'].unique())}")

    records = []
    size_map = dict(zip(chrom_sizes["chromosome"], chrom_sizes["length"]))

    for ch in chrom_order:
        sub_te = compute_local_density(
            gff, ch, window_size, feature="TE", chrom_length=size_map.get(ch)
        )
        if not sub_te.empty:
            records.extend(sub_te.to_dict(orient="records"))

    df = pd.DataFrame(records)
    if df.empty:
        print("[WARN] No TE density records were found. Please check the chromosome names.")
        return pd.DataFrame(columns=["chromosome", "window_start_bp", "TE_density_bp"])

    df = df.groupby(["chromosome", "window_start_bp"], as_index=False).sum()
    df.to_csv(output_file, sep="\t", index=False)
    print(f"[INFO] TE density data saved to {output_file}")
    return df

# -----------------------------
# Main
# -----------------------------
def main():
    # --- Chr to keep ---
    chrom_order = [f"Chr{i}" for i in range(1, 18)]
    keep_chroms_raw = chrom_order[:]

    # --- Chr sizes ---
    chrom_sizes = read_fai(FAI_FILE)
    chrom_sizes = chrom_sizes[chrom_sizes["chromosome"].isin(keep_chroms_raw)].copy()
    chrom_sizes["Mbp_total"] = chrom_sizes["length"] / 1e6
    chrom_sizes["chromosome"] = pd.Categorical(chrom_sizes["chromosome"],
                                               categories=chrom_order, ordered=True)
    chrom_sizes = chrom_sizes.sort_values("chromosome")

    # --- TE GFF ---
    gff_cols = ["chromosome","source","type","start","end","score","strand","phase","attributes"]
    gff = pd.read_csv(GFF_FILE, sep="\t", comment="#", header=None,
                      names=gff_cols, dtype={"chromosome": str})
    gff = gff[gff["type"] == "match"].copy()
    gff["start"] = gff["start"].astype(int)
    gff["end"] = gff["end"].astype(int)
    gff = gff[gff["chromosome"].isin(keep_chroms_raw)].copy()
    gff["Seq_name"] = gff["attributes"].apply(extract_target)
    gff["identity"] = gff["attributes"].apply(extract_identity)

    # --- Classification ---
    classif = pd.read_csv(CLASSIF_FILE, sep="\t")
    for col in classif.columns:
        if classif[col].dtype == object:
            classif[col] = classif[col].apply(keep_before_pipe)
    gff = gff.merge(classif, on="Seq_name", how="left", suffixes=("", "_cls"))
    gff["order"] = gff["order"].fillna("Unclassified")

    # --- Pre-computing TE densities ---
    density_df = compute_and_store_te_density(gff, chrom_order, chrom_sizes,
                                              WINDOW_SIZE, OUTPUT_FILE)

    # --- TE coverage pour (a) ---
    gff["length"] = gff["end"] - gff["start"] + 1
    te_bp = (gff.groupby(["chromosome","order"])["length"].sum()
             .reset_index()
             .rename(columns={"length":"bp"}))
    te_bp["Mbp"] = te_bp["bp"] / 1e6
    te_bp["chromosome"] = pd.Categorical(te_bp["chromosome"], categories=chrom_order, ordered=True)

    # Remaining (no TE coverage)
    chrom_sizes = chrom_sizes.set_index("chromosome")
    te_total = te_bp.groupby("chromosome")["Mbp"].sum().reindex(chrom_sizes.index).fillna(0.0)
    remainder = (chrom_sizes["Mbp_total"] - te_total).clip(lower=0.0)

    # -----------------------------
    # Figure
    # -----------------------------
    sns.set(style="white")
    fig, axes = plt.subplots(2, 2, figsize=(22, 12))
    ax_a, ax_b, ax_c, ax_d = axes.flatten()

    # ---------- (a) TE composition ----------
    ax_a.grid(False)
    xcats = chrom_sizes.index.tolist()
    ax_a.bar(xcats, chrom_sizes.loc[xcats, "Mbp_total"].values,
             facecolor="none", edgecolor="black", linewidth=1.2, zorder=0)
    
    orders = te_bp["order"].unique().tolist()
    base_palette = sns.color_palette("tab10", n_colors=len(orders))
    color_map = {o: base_palette[i % len(base_palette)] for i, o in enumerate(orders)}
    
    # colors
    if "Unclassified" in color_map:
        color_map["Unclassified"] = (0.6, 0.6, 0.6)  # grey
    if "Maverick" in color_map:
        color_map["Maverick"] = (0, 0, 0)            # black
    if "SINE" in color_map:
        color_map["SINE"] = (1.0, 0.4, 0.7)          # pink (RGB)
    
    bottoms = {c: 0.0 for c in xcats}
    for o in orders:
        sub = te_bp[te_bp["order"] == o]
        heights = [float(sub.loc[sub["chromosome"] == c, "Mbp"].sum()) for c in xcats]
        ax_a.bar(xcats, heights,
                 bottom=[bottoms[c] for c in xcats],
                 color=color_map[o], edgecolor="black", linewidth=0.8, label=o, zorder=2)
        for c, h in zip(xcats, heights):
            bottoms[c] += h
    
    rem_heights = [float(remainder.loc[c]) for c in xcats]
    if any(rem_heights):
        ax_a.bar(xcats, rem_heights,
                 bottom=[bottoms[c] for c in xcats],
                 facecolor="white", edgecolor="black", linewidth=0.8,
                 hatch="//", label="No TE coverage", zorder=1)

    ax_a.set_ylabel("Length (Mbp)")
    ax_a.set_title(r"$\bf{(a)}$ TE composition within chromosome size")
    ax_a.tick_params(axis="x", rotation=45)

    handles, labels = ax_a.get_legend_handles_labels()
    if "No TE coverage" in labels:
        idx = labels.index("No TE coverage")
        handles.append(handles.pop(idx))
        labels.append(labels.pop(idx))
    ax_a.legend(handles, labels, title="Category",
                bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=True)

    # ---------- (b) Violin log10(TE length) ----------
    sns.violinplot(data=gff, x="chromosome", y=np.log10(gff["length"]), ax=ax_b,
                   inner=None, cut=0, scale="width", color="white")
    sns.boxplot(data=gff, x="chromosome", y=np.log10(gff["length"]), ax=ax_b,
                showcaps=True, boxprops={"facecolor":"white", "edgecolor":"black"},
                whiskerprops={"color":"black"}, medianprops={"color":"black"},
                showfliers=False, width=0.25)
    ax_b.set_ylabel("log10(TE length [bp])")
    ax_b.set_title(r"$\bf{(b)}$ TE length distribution")
    ax_b.tick_params(axis="x", rotation=45)

    # ---------- (c) Violin TE identity ----------
    sns.violinplot(data=gff, x="chromosome", y="identity", ax=ax_c,
                   inner=None, cut=0, scale="width", color="white")
    sns.boxplot(data=gff, x="chromosome", y="identity", ax=ax_c,
                showcaps=True, boxprops={"facecolor":"white", "edgecolor":"black"},
                whiskerprops={"color":"black"}, medianprops={"color":"black"},
                showfliers=False, width=0.25)
    ax_c.set_ylabel("TE identity (%)")
    ax_c.set_title(r"$\bf{(c)}$ TE identity distribution")
    ax_c.tick_params(axis="x", rotation=45)

    # ---------- (d) TE density : Automatically adjusted grid ----------
    nchroms = len(chrom_order)

    # Compute a grid large enough to include all chromosomes.
    ncols = 5
    nrows = int(np.ceil(nchroms / ncols))

    subfig = ax_d.figure.add_gridspec(nrows, ncols, left=0.55, right=0.95,
                                      bottom=0.05, top=0.45, wspace=0.25, hspace=0.3)
    sub_axes = [fig.add_subplot(subfig[i, j]) for i in range(nrows) for j in range(ncols)]

    max_te = (density_df["TE_density_bp"].max() / 1e6) if not density_df.empty else 1.0

    for i, ch in enumerate(chrom_order):
        ax_sub = sub_axes[i]
        lw = 1.5
        sub = density_df[density_df["chromosome"] == ch]

        if not sub.empty:
            ax_sub.plot(sub["window_start_bp"]/1e6,
                        sub["TE_density_bp"]/1e6,
                        color="black", linestyle="-", linewidth=lw)

        ax_sub.set_ylim(0, max_te*1.05)
        ax_sub.set_title(ch, fontsize=8, weight="bold")

        # position in the grid
        row_idx = i // ncols
        col_idx = i % ncols
        
        # ylabel :Apply the same rule to the first column and the entire second row.
        if col_idx != 0 or row_idx != 1:
            ax_sub.set_ylabel("")
        else:
            ax_sub.set_ylabel("TE density (Mbp/100 kb)")

        # Display x-axis labels only on the bottom row.
        if row_idx == nrows - 1:
            ax_sub.set_xlabel("Position (Mbp)")
        else:
            ax_sub.set_xlabel(""); ax_sub.set_xticklabels([])

    # Hide unused empty subplots.
    for j in range(nchroms, nrows*ncols):
        sub_axes[j].axis("off")

    ax_d.axis("off")
    ax_d.set_title(r"$\bf{(d)}$ Local TE density (100 kb windows)",
                   fontsize=12, weight="bold")

    # Legend
    custom_lines = [plt.Line2D([0], [0], color="black", linestyle="-", lw=2)]
    fig.legend(custom_lines, ["TE density"],
               loc="center", bbox_to_anchor=(0.5,0.5), frameon=True)

    fig.tight_layout()
    plt.subplots_adjust(right=0.95)
    fig.savefig(FIG_DIR / "REPET_summary.png", dpi=300)

if __name__ == "__main__":
    main()
