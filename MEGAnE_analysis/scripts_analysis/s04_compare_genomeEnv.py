#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ======================================================
# PARAMETERS
# ======================================================

SPECIES = ["PCOM", "PPY"]
BASE_DIR = "../output"
FIG_DIR = "../Figures"

OUT_PNG = f"{FIG_DIR}/rho_TE_gene_density_PCOM_PPY.png"
OUT_STATS = f"{FIG_DIR}/rho_TE_gene_density_PCOM_PPY.stats.txt"

os.makedirs(FIG_DIR, exist_ok=True)

CONDITIONS = [
    "Whole genome",
    "Upstream selected genes",
    "Upstream non-selected genes"
]

METRICS = [
    ("Rho", "Rho"),
    ("TE_density", "TE density"),
    ("gene_density", "Gene density")
]

# ======================================================
# STATISTICS UTILITIES
# ======================================================

def p_to_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def pairwise_tests(df, value_col):
    groups = {
        c: df[df["Condition"] == c][value_col]
        for c in CONDITIONS
    }

    comparisons = [
        (CONDITIONS[0], CONDITIONS[1]),
        (CONDITIONS[0], CONDITIONS[2]),
        (CONDITIONS[1], CONDITIONS[2]),
    ]

    raw_p = []
    pairs = []

    for g1, g2 in comparisons:
        _, p = mannwhitneyu(groups[g1], groups[g2], alternative="two-sided")
        raw_p.append(p)
        pairs.append((g1, g2))

    _, p_corr, _, _ = multipletests(raw_p, method="fdr_bh")

    return [(g1, g2, p) for (g1, g2), p in zip(pairs, p_corr)]


def add_stat_annotations(ax, df, value_col):
    results = pairwise_tests(df, value_col)

    y_max = df[value_col].max()
    y_step = y_max * 0.1
    y = y_max + y_step

    x_pos = {c: i for i, c in enumerate(CONDITIONS)}

    for g1, g2, p in results:
        x1, x2 = x_pos[g1], x_pos[g2]

        ax.plot(
            [x1, x1, x2, x2],
            [y, y + y_step / 3, y + y_step / 3, y],
            lw=1,
            c="black"
        )

        ax.text(
            (x1 + x2) / 2,
            y + y_step * 0.3, #/ 2,
            p_to_stars(p),
            ha="center",
            va="bottom",
            fontsize=9,
            #weight="bold"
        )

        y += y_step


def write_stats(fh, species, metric_label, df):
    fh.write(f"\n### {species} | {metric_label}\n")

    groups = {
        c: df[df["Condition"] == c]["Value"].dropna()
        for c in CONDITIONS
    }

    comparisons = [
        (CONDITIONS[0], CONDITIONS[1]),
        (CONDITIONS[0], CONDITIONS[2]),
        (CONDITIONS[1], CONDITIONS[2]),
    ]

    raw_p = []
    pairs = []

    for g1, g2 in comparisons:
        _, p = mannwhitneyu(groups[g1], groups[g2], alternative="two-sided")
        raw_p.append(p)
        pairs.append((g1, g2))

    _, p_corr, _, _ = multipletests(raw_p, method="fdr_bh")

    for (g1, g2), p in zip(pairs, p_corr):
        n1, n2 = len(groups[g1]), len(groups[g2])
        fh.write(
            f"{g1} vs {g2} : p_adj = {p:.3e} "
            f"({p_to_stars(p)}), n = {n1} vs {n2}\n"
        )


# ======================================================
# DATA PREPARATION
# ======================================================

def prepare_plot_df(df):
    df = df.copy()

    for col in ["Rho", "TE_density"]: #, "gene_density"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Rho", "TE_density"]) #, "gene_density"])
    
    df_all = df
    # # TIPs associated with AT LEAST one selected gene
    # df_sel = df[
    #     df["positive_selection"].str.contains("true", case=False, na=False)
    # ]
    # # TIPs associated with AT LEAST one non-selected gene
    # df_nonsel = df[
    #     df["positive_selection"].str.contains("false", case=False, na=False)
    # ]
    def has_true(x):
        return any(v.strip() == "true" for v in str(x).split(","))
    
    mask_true = df["positive_selection"].apply(has_true)
    
    df_sel = df[
        (df["gene_Gene_ID"].notna())
        & (df["gene_Gene_ID"] != "NA")
        & mask_true
    ]
    
    df_nonsel = df[
        (df["gene_Gene_ID"].notna())
        & (df["gene_Gene_ID"] != "NA")
        & ~mask_true
    ]

    out = []

    for col, label in METRICS:
        out.append(pd.DataFrame({
            "Value": df_all[col],
            "Metric": label,
            "Condition": CONDITIONS[0]
        }))
        out.append(pd.DataFrame({
            "Value": df_sel[col],
            "Metric": label,
            "Condition": CONDITIONS[1]
        }))
        out.append(pd.DataFrame({
            "Value": df_nonsel[col],
            "Metric": label,
            "Condition": CONDITIONS[2]
        }))

    return pd.concat(out, ignore_index=True)

def add_median_line(ax, df):
    """
    Draw a red line connecting the medians of each condition.
    """
    medians = (
        df.groupby("Condition")["Value"]
        .median()
        .reindex(CONDITIONS)
    )

    x = range(len(CONDITIONS))
    y = medians.values

    ax.plot(
        x, y,
        color="red",
        linewidth=1,
        marker="o",
        alpha=0.9,
        markersize=3,
        zorder=5
    )

def add_N_annotations(ax, df, y_offset=0.73):
    """
    Annotate each violin with N (number of TIPs) next to the plot.
    """
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    y_text = y_min + y_range * y_offset

    for i, cond in enumerate(CONDITIONS):
        n = df[df["Condition"] == cond].shape[0]
        ax.text(
            i,
            y_text,
            f"N={n}",
            ha="center",
            va="bottom",
            fontsize=9,
            color="dimgray"
        )


# ======================================================
# MAIN PLOTTING PIPELINE
# ======================================================

def main():

    stats_fh = open(OUT_STATS, "w")

    fig, axes = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=(18, 8),
        sharex="col"
    )

    for row, sp in enumerate(SPECIES):

        print(f"Processing {sp}...")

        tsv = f"{BASE_DIR}/{sp}/MEGAnE_summary_{sp}.tsv"
        df = pd.read_csv(tsv, sep="\t")

        plot_df = prepare_plot_df(df)

        for col, (_, metric_label) in enumerate(METRICS):

            ax = axes[row, col]
            sub = plot_df[plot_df["Metric"] == metric_label]

            write_stats(stats_fh, sp, metric_label, sub)

            sns.violinplot(
                data=sub,
                x="Condition",
                y="Value",
                inner="box",
                cut=0,
                linewidth=1,
                ax=ax
            )

            add_median_line(ax, sub)
            add_stat_annotations(ax, sub, "Value")
            add_N_annotations(ax, sub)

            if row == 0:
                ax.set_title(metric_label, fontsize=12) #, weight="bold")

            ax.set_xlabel("")
            ax.tick_params(axis="x", rotation=20)

        axes[row, 0].set_ylabel(
            sp,
            fontsize=15,
            weight="bold",
            rotation=90,
            labelpad=25
        )

    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300)
    plt.close()

    stats_fh.close()

    print(f"\n[OK] Figure saved → {OUT_PNG}")
    print(f"[OK] Stats written → {OUT_STATS}")


# ======================================================
# ENTRY POINT
# ======================================================

if __name__ == "__main__":
    main()
