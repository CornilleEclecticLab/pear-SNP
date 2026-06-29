#!/usr/bin/env python3
"""
s05_compare_Immunity.py

TIP-level enrichment analysis faithful to the UpSet plot.

Each row in the input table is one TIP (TE insertion polymorphism) record.
The boolean column "Immunity" flags TIPs that meet your criterion:
- TIP within 2 kb upstream of a positively selected gene involved in immunity (per your definition).

UpSet sets (4 flags):
- comm_Perry
- comm_Dessert
- pyra
- cauc

The UpSet plot displays only the top intersections (bars). To stay faithful to the figure,
we filter to the top N intersections by size (N=15, matching your plot).

Categorization (fine-grained as requested):
Wild-only split:
- Wild_cauc_only
- Wild_pyra_only
- Wild_pyra_cauc

Cultivated-only split:
- Cult_Perry_only
- Cult_Dessert_only
- Cult_Perry_Dessert

Mixed:
- Mix_Wild_Cultivated  (sum of remaining displayed intersections)

Shared-by-all:
- Shared_by_all (all four sets)
  NOTE: Shared_by_all is excluded from ALL statistical tests (as requested),
        but can still be reported in the descriptive summary.

Statistical analysis (excluding Shared_by_all):
- Global chi-square test across categories (multi-group)
- Pairwise tests across all category pairs (exploratory), with BH-FDR across all pairs:
    * Fisher's exact test (exact, good for small counts)
    * Chi-square test on 2x2 tables (approximate, shown for comparison)
- Targeted contrasts (in your chosen order), with BH-FDR across targeted tests only:
    * Fisher + Chi-square on 2x2 tables

Input:
- ../data/PCOM/Table_S14.tsv

Required columns:
- comm_Perry, comm_Dessert, pyra, cauc : TRUE/FALSE membership flags
- Immunity : TRUE/FALSE flag as described above

Run:
python3 s05_compare_Immunity.py
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests


# -------------------------------
# User parameters
# -------------------------------
INPUT = "../data/PCOM/Table_S14.tsv"
TOP_N_INTERSECTIONS = 15  # number of bars shown in your UpSet figure

# Targeted contrasts (YOUR order)
TARGETED_CONTRASTS = [
    ("Cult_Perry_only", "Cult_Dessert_only"),
    ("Cult_Perry_only", "Cult_Perry_Dessert"),
    ("Cult_Dessert_only", "Cult_Perry_Dessert"),
    ("Cult_Perry_Dessert", "Wild_pyra_cauc"),
]


# ===============================
# 1. Load data
# ===============================
df = pd.read_csv(INPUT, sep="\t")

# ===============================
# 2. Validate + convert TRUE/FALSE columns to booleans
# ===============================
required = ["comm_Perry", "comm_Dessert", "pyra", "cauc", "Immunity"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

for col in required:
    df[col] = df[col].astype(str).str.upper().eq("TRUE")

# Convenience variables
P = df["comm_Perry"]
D = df["comm_Dessert"]
Y = df["pyra"]
C = df["cauc"]

# ===============================
# 3. Build exact intersection encoding + filter to top intersections (UpSet bars)
# ===============================
df["Intersection_bits"] = (
    P.astype(int).astype(str) + "_" +
    D.astype(int).astype(str) + "_" +
    Y.astype(int).astype(str) + "_" +
    C.astype(int).astype(str)
)

def intersection_label(row) -> str:
    labels = []
    if row["comm_Perry"]:
        labels.append("comm_Perry")
    if row["comm_Dessert"]:
        labels.append("comm_Dessert")
    if row["pyra"]:
        labels.append("pyra")
    if row["cauc"]:
        labels.append("cauc")
    return "&".join(labels) if labels else "NONE"

df["Intersection_label"] = df.apply(intersection_label, axis=1)

top_bits = df["Intersection_bits"].value_counts().head(TOP_N_INTERSECTIONS).index
df = df[df["Intersection_bits"].isin(top_bits)].copy()

print("\nExact intersection sizes (TIP counts per exact combination):")
print(df["Intersection_label"].value_counts())

# Recompute convenience variables after filtering (safe)
P = df["comm_Perry"]
D = df["comm_Dessert"]
Y = df["pyra"]
C = df["cauc"]

cult = P | D
wild = Y | C
all_four = P & D & Y & C

# ===============================
# 4. Assign fine-grained categories
# ===============================
df["Category"] = None

# Shared by all (red)
df.loc[all_four, "Category"] = "Shared_by_all"

# Wild-only (green), split into subgroups
df.loc[(C & ~Y & ~cult & ~all_four), "Category"] = "Wild_cauc_only"
df.loc[(Y & ~C & ~cult & ~all_four), "Category"] = "Wild_pyra_only"
df.loc[(Y & C & ~cult & ~all_four), "Category"] = "Wild_pyra_cauc"

# Cultivated-only (blue), split into subgroups
df.loc[(P & ~D & ~wild & ~all_four), "Category"] = "Cult_Perry_only"
df.loc[(D & ~P & ~wild & ~all_four), "Category"] = "Cult_Dessert_only"
df.loc[(P & D & ~wild & ~all_four), "Category"] = "Cult_Perry_Dessert"

# Mixed (black): everything else among displayed intersections (excluding Shared_by_all)
df.loc[df["Category"].isna(), "Category"] = "Mix_Wild_Cultivated"

if df["Category"].isna().any():
    raise RuntimeError("Some rows were not assigned to any Category. Check logic.")

# ===============================
# 5. TIP-level descriptive summary (includes Shared_by_all)
# ===============================
summary = df.groupby("Category")["Immunity"].agg(
    total="count",
    immunity_TIPs="sum"
)
summary["proportion"] = summary["immunity_TIPs"] / summary["total"]

print("\nProportion of immunity-associated TIPs by UpSet category "
      "(TIP within 2 kb upstream; genes under positive selection):")
print(summary.sort_values("proportion", ascending=False))

# ===============================
# 6. Exclude Shared_by_all from statistical tests
# ===============================
df_test = df[df["Category"] != "Shared_by_all"].copy()

# ===============================
# 7. Global chi-square test across categories (excluding Shared_by_all)
# ===============================
contingency_global = pd.crosstab(df_test["Category"], df_test["Immunity"]).reindex(columns=[False, True], fill_value=0)
chi2_global, p_global, _, _ = chi2_contingency(contingency_global)
print("\nGlobal chi-square test p-value (excluding Shared_by_all):", p_global)

# ===============================
# 8. Helper functions for 2x2 tests (Fisher + Chi-square)
# ===============================
def build_2x2_table(data: pd.DataFrame, g1: str, g2: str) -> pd.DataFrame:
    """
    Build a 2x2 table for two categories:
    rows = [g1, g2]
    cols = [False, True] for Immunity
    """
    sub = data[data["Category"].isin([g1, g2])].copy()
    table = pd.crosstab(sub["Category"], sub["Immunity"]).reindex(columns=[False, True], fill_value=0)
    # Ensure row order
    table = table.reindex(index=[g1, g2], fill_value=0)
    return table

def run_fisher_and_chi2(table_2x2: pd.DataFrame) -> dict:
    """
    Run Fisher exact and chi-square on the same 2x2 table.
    Returns p-values and odds ratio (from Fisher).
    """
    oddsratio, p_fisher = fisher_exact(table_2x2.values)

    # Chi-square test on 2x2 table (approximate)
    chi2, p_chi2, _, _ = chi2_contingency(table_2x2.values)

    return {
        "p_fisher": p_fisher,
        "odds_ratio": oddsratio,
        "p_chi2": p_chi2,
    }

# ===============================
# 9. Exploratory pairwise tests across all category pairs (excluding Shared_by_all)
#    + BH-FDR correction across all pairs
# ===============================
categories = sorted(df_test["Category"].unique().tolist())

pairwise_rows = []
for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        g1, g2 = categories[i], categories[j]
        table = build_2x2_table(df_test, g1, g2)
        stats = run_fisher_and_chi2(table)

        pairwise_rows.append({
            "Category1": g1,
            "Category2": g2,
            "p_value_fisher": stats["p_fisher"],
            "Odds_ratio": stats["odds_ratio"],
            "p_value_chi2": stats["p_chi2"],
        })

pairwise_df = pd.DataFrame(pairwise_rows)

# BH-FDR across all exploratory pairwise Fisher tests
pairwise_df["p_adj_FDR_fisher"] = multipletests(pairwise_df["p_value_fisher"], method="fdr_bh")[1]

# BH-FDR across all exploratory pairwise Chi-square tests (separately)
pairwise_df["p_adj_FDR_chi2"] = multipletests(pairwise_df["p_value_chi2"], method="fdr_bh")[1]

print("\nPairwise Fisher tests among categories (excluding Shared_by_all) "
      "(exploratory, BH-FDR corrected across all pairs):")
print(pairwise_df.sort_values("p_adj_FDR_fisher")[["Category1", "Category2", "p_value_fisher", "Odds_ratio", "p_adj_FDR_fisher"]])

print("\nPairwise Chi-square tests among categories (excluding Shared_by_all) "
      "(exploratory, BH-FDR corrected across all pairs):")
print(pairwise_df.sort_values("p_adj_FDR_chi2")[["Category1", "Category2", "p_value_chi2", "p_adj_FDR_chi2"]])

# ===============================
# 10. Targeted contrasts (YOUR order): Fisher + Chi-square + BH-FDR over targets only
# ===============================
targeted_rows = []
available = set(categories)

for g1, g2 in TARGETED_CONTRASTS:
    if g1 == g2:
        print(f"\nSkipping invalid comparison: {g1} vs {g2} (same group)")
        continue
    if g1 not in available or g2 not in available:
        print(f"\nSkipping comparison {g1} vs {g2}: missing category in the filtered data.")
        continue

    table = build_2x2_table(df_test, g1, g2)
    stats = run_fisher_and_chi2(table)

    # TIP-level proportions for interpretation
    sub = df_test[df_test["Category"].isin([g1, g2])]
    props = sub.groupby("Category")["Immunity"].mean()

    print(f"\nTargeted test: {g1} vs {g2}")
    print("2x2 table (rows=category, cols=[False, True]):")
    print(table)
    print(f"Fisher p-value: {stats['p_fisher']:.6g} | Odds ratio: {stats['odds_ratio']:.4g}")
    print(f"Chi-square p-value: {stats['p_chi2']:.6g}")
    print("Proportions (TIP-level; mean of Immunity flag within category):")
    print(props)

    targeted_rows.append({
        "Category1": g1,
        "Category2": g2,
        "p_value_fisher": stats["p_fisher"],
        "Odds_ratio": stats["odds_ratio"],
        "p_value_chi2": stats["p_chi2"],
    })

if targeted_rows:
    targeted_df = pd.DataFrame(targeted_rows)

    # BH-FDR over targeted tests only (Fisher)
    targeted_df["p_adj_FDR_fisher"] = multipletests(targeted_df["p_value_fisher"], method="fdr_bh")[1]

    # BH-FDR over targeted tests only (Chi-square), separately
    targeted_df["p_adj_FDR_chi2"] = multipletests(targeted_df["p_value_chi2"], method="fdr_bh")[1]

    # Preserve YOUR requested order
    order_key = {(a, b): i for i, (a, b) in enumerate(TARGETED_CONTRASTS)}
    targeted_df["order"] = targeted_df.apply(lambda r: order_key.get((r["Category1"], r["Category2"]), 999), axis=1)
    targeted_df = targeted_df.sort_values("order").drop(columns=["order"])

    print("\nTargeted contrasts summary (BH-FDR corrected across targeted tests only):")
    print(targeted_df)
else:
    print("\nNo valid targeted contrasts were run.")
    
# ===============================
# EXTRA — Higher power test by pooling cultivated vs wild
# ===============================
# Rationale:
# - Some fine-grained categories have few Immunity=True counts (e.g., 3, 5, 7),
#   which reduces power in pairwise tests.
# - Pooling into two biologically motivated groups increases counts and power:
#     Cultivated_only = Cult_Perry_only + Cult_Dessert_only + Cult_Perry_Dessert
#     Wild_only       = Wild_cauc_only + Wild_pyra_only + Wild_pyra_cauc
# - Shared_by_all is excluded (as before). Mixed is excluded from this pooled test.

# Define pooled groups from existing Category labels
pool_map = {
    "Cult_Perry_only": "Cultivated_only",
    "Cult_Dessert_only": "Cultivated_only",
    "Cult_Perry_Dessert": "Cultivated_only",
    "Wild_cauc_only": "Wild_only",
    "Wild_pyra_only": "Wild_only",
    "Wild_pyra_cauc": "Wild_only",
}

df_pool = df[df["Category"].isin(pool_map.keys())].copy()
df_pool["Pooled"] = df_pool["Category"].map(pool_map)

# Build 2x2 table and run Fisher + Chi-square
pool_table = pd.crosstab(df_pool["Pooled"], df_pool["Immunity"]).reindex(columns=[False, True], fill_value=0)
pool_table = pool_table.reindex(index=["Cultivated_only", "Wild_only"], fill_value=0)

oddsratio_pool, p_fisher_pool = fisher_exact(pool_table.values)
chi2_pool, p_chi2_pool, _, _ = chi2_contingency(pool_table.values)

print("\nPooled comparison (higher power): Cultivated_only vs Wild_only")
print("2x2 table (rows=pooled group, cols=[False, True]):")
print(pool_table)

# TIP-level proportions for interpretation
props_pool = df_pool.groupby("Pooled")["Immunity"].mean()
print("Proportions (TIP-level; mean of Immunity flag within pooled group):")
print(props_pool)

print(f"Fisher p-value: {p_fisher_pool:.6g} | Odds ratio (Cultivated vs Wild): {oddsratio_pool:.4g}")
print(f"Chi-square p-value: {p_chi2_pool:.6g}")