#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate multiple UpSet plots for TE insertions polymorphism
across different gene-selection contexts:
1. Genome-wide
2. TEs near genes under positive selection
3. TEs near genes NOT under selection
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import from_indicators, UpSet
from itertools import combinations
import numpy as np





def split_data( df, POP_ORDER ):

    # Clean data and normalize the field to lowercase strings
    df = df.dropna(subset=["populations"])
    df = df[df["populations"].astype(str).str.lower() != "na"]
    df["populations"] = df["populations"].astype(str).str.strip()
    df = df[df["populations"] != ""]
    df["populations"] = df["populations"].str.split(",")
    df["positive_selection"] = df["positive_selection"].astype(str).str.lower()
   
    # Remove singleton TEs
    #df = df[df["populations"].str.len() > 1]
    
    df = build_boolean_indicators(df, POP_ORDER)
    
    ### ------- (1) Genome-wide (all TEs) ------- 
    df_all = df.copy()

    # Define a helper function to detect if *any* entry is true
    def has_true(x):
        # Split by commas, strip spaces, and check if any equals "true"
        return any(val.strip() == "true" for val in x.split(","))
    
    # Apply to get a boolean mask
    mask_true = df["positive_selection"].apply(has_true)
    mask_false = ~mask_true  # everything else
    
    ### ------- (2) TEs upstream of genes under positive selection ------- 
    df_selected = df[
        (df["gene_Gene_ID"].notna())
        & (df["gene_Gene_ID"] != "NA")
        & mask_true
    ]
    
    ### ------- (3) TEs upstream of genes NOT under positive selection ------- 
    df_nonselected = df[
        (df["gene_Gene_ID"].notna())
        & (df["gene_Gene_ID"] != "NA")
        & mask_false
    ]

    # --- Display subset info
    # print("\n--- Data Subsets ---")
    # print("Genome-wide (all):")
    # print(df_all.head(), "\n")
    # print(f"→ {len(df_all)} rows\n")
    
    # print("Genes under positive selection:")
    # print(df_selected.head(), "\n")
    # print(f"→ {len(df_selected)} rows\n")
    
    # print("Genes NOT under positive selection:")
    # print(df_nonselected.head(), "\n")
    # print(f"→ {len(df_nonselected)} rows\n")
    
    return df_all, df_selected, df_nonselected






def determine_color(combo, SPECIES):
    """Return a unique or mean color for a population combination."""
    combo = set(combo)

    if SPECIES == 'PCOM':
        wild = {"cauc", "pyra"}
        cultivar = {"comm_Dessert", "comm_Perry"}
        if combo == wild.union(cultivar):
            return "red"
        elif combo <= wild:
            return "#158B22"
        elif combo <= cultivar:
            return "#545EFA"
        else:
            return "black"
    elif SPECIES == 'PPY':
        wild = {"ussu", "pash"}
        cultivar = {"Sand_CN-SE", "pyri_JP", "Sand_CN", "Sand_CN-SW", "White"}
        rootstock = {"betu"}
        if combo == wild.union(cultivar, rootstock):
            return "red"
        if combo == rootstock.union(cultivar):
            return "pink"
        elif combo <= wild:
            return "#158B22"
        elif combo <= cultivar:
            return "#545EFA"
        elif combo <= rootstock:
            return "#545EFA"
        else:
            return "black"





def build_boolean_indicators(df, POP_ORDER):
    """Add boolean columns (True/False) for each population in POP_ORDER."""
    for pop in POP_ORDER:
        df[pop] = df["populations"].apply(
            lambda pops: pop in pops if isinstance(pops, list) else False
        )
    return df





		
def make_upset_plot(SPECIES, df_init, POP_ORDER, POP_COLORS, filename, legend, filename_df, Type):
    """
    Create an UpSet plot of TE polymorphisms found at 2kbp upstream of a set of genes
    """
    # Parameters
    d_upsetplot_color = {"#158B22" : "Wild",
                            "#545EFA": "Cultivated", 
                            "red": "Shared by all", 
                            "black": "mix Wild-Culti", 
                            "orange": "Rootstock", 
                            "pink": "mix Culti-RootS"
                        }

    # Add log scale for Rho
    # print(df_init)
    # df_init["Rho (log10 scale)"] = np.log10(df_init["Rho"].clip(lower=1e-10))
    # print(df_init["Rho (log10 scale)"])
    # print(df_init)
    
    # Prepare data
    ordered_pops = POP_ORDER[::-1]
    print(f"Ordered populations: {ordered_pops}")
    indicators = df_init[ordered_pops]
    print(f"Indicator matrix shape: {indicators.shape}")
    data = from_indicators(ordered_pops, data=indicators)
    print("Indicator matrix sample:")
    print(indicators.head())
    print("Unique MultiIndex?", data.index.is_unique)
    print("Duplicate count:", data.index.duplicated().sum())

    print(f"==== Before adding columns ===> Shape of data: {data.shape}")
    print(f"Index type: {type(data.index).__name__}")
    print(f"data: {data}")
    cols_to_copy = [col for col in df_init.columns if col not in data.columns]
    for col in cols_to_copy:
        data[col] = df_init[col].values
    print(f"==== After adding colums ===> Shape of data: {data.shape}")
    print(f"Index type: {type(data.index).__name__}")
    print(f"data: {data}")

    _key = data[ordered_pops].astype(int).astype(str).agg(''.join, axis=1)
    data["_key"] = _key
    print(f"data after adding keys: {data}")
   
    def percent_by_group(df, group_cols, value_name):
        counts = df.groupby(group_cols).size().reset_index(name="count")
        counts[value_name] = 100 * counts["count"] / counts.groupby(group_cols[0])["count"].transform("sum")
        return counts.drop(columns="count")

    # Remove duplicate genes per intersection (_key)
    data_gene_unique = data.drop_duplicates(subset=["_key", "gene_Gene_ID"])

    # df_class_pct = percent_by_group(data, ["_key", "class_curated"], "TE class (%)")
    # df_order_pct = percent_by_group(data, ["_key", "order_curated"], "TE order (%)")
    # df_go_pct    = percent_by_group(data, ["_key", "GO_Level2_Category"], "GO Level 2 (%)")
    df_class_pct = percent_by_group(data_gene_unique, ["_key", "class_curated"], "TE class (%)")
    df_order_pct = percent_by_group(data_gene_unique, ["_key", "order_curated"], "TE order (%)")
    df_go_pct    = percent_by_group(data_gene_unique, ["_key", "GO_Level2_Category"], "GO Level 2 (%)")
    df_immunity_pct    = percent_by_group(data_gene_unique, ["_key", "Immunity"], "Immunity (%)")

    class_dict = dict(zip(zip(df_class_pct["_key"], df_class_pct["class_curated"]), df_class_pct["TE class (%)"]))
    order_dict = dict(zip(zip(df_order_pct["_key"], df_order_pct["order_curated"]), df_order_pct["TE order (%)"]))
    go_dict    = dict(zip(zip(df_go_pct["_key"], df_go_pct["GO_Level2_Category"]), df_go_pct["GO Level 2 (%)"]))
    immunity_dict    = dict(zip(zip(df_immunity_pct["_key"], df_immunity_pct["Immunity"]), df_immunity_pct["Immunity (%)"]))
    
    data["TE class (%)"]   = [class_dict.get((k, c), 0.0) for k, c in zip(data["_key"], data["class_curated"])]
    data["TE order (%)"]   = [order_dict.get((k, o), 0.0) for k, o in zip(data["_key"], data["order_curated"])]
    data["GO Level 2 (%)"] = [go_dict.get((k, g), 0.0) for k, g in zip(data["_key"], data["GO_Level2_Category"])]
    data["Immunity (%)"] = [immunity_dict.get((k, g), 0.0) for k, g in zip(data["_key"], data["Immunity"])]

    print(f"data after adding proportions:\n{data}")
    
    print("\n==============================")
    print("🔍 VERIFYING GROUPBY COUNTS AND PROPORTIONS")
    print("==============================")
    
    # 1️⃣ Check total row counts per _key
    group_counts = data.groupby("_key").size().reset_index(name="total_count")
    print("\nTotal element counts per _key (first 5):")
    print(group_counts.head())
    
    # 2️⃣ Recompute counts directly from data (for validation)
    df_class_check = data.groupby(["_key", "class_curated"]).size().reset_index(name="count_direct")
    
    # Merge with df_class_pct to compare the computed percentages
    merged_check = pd.merge(
        df_class_pct,
        df_class_check,
        on=["_key", "class_curated"],
        how="left"
    )
    merged_check["recomputed_pct"] = 100 * merged_check["count_direct"] / merged_check.groupby("_key")["count_direct"].transform("sum")
    
    # 3️⃣ Check if original and recomputed percentages match closely
    merged_check["pct_diff"] = merged_check["TE class (%)"] - merged_check["recomputed_pct"]
    mean_diff = merged_check["pct_diff"].abs().mean()
    max_diff = merged_check["pct_diff"].abs().max()
    
    print(f"\nMean difference between stored and recomputed TE class (%): {mean_diff:.4f}")
    print(f"Max difference: {max_diff:.4f}")
    
    # 4️⃣ Show discrepancies if any
    discrepancies = merged_check[merged_check["pct_diff"].abs() > 1]
    if len(discrepancies) > 0:
        print(f"\n⚠️ Found {len(discrepancies)} groups with >1% mismatch:")
        print(discrepancies.head(10))
    else:
        print("\n✅ All class_curated group percentages match recomputed values!")
    
    # 5️⃣ Random sample: inspect one _key group in detail
    sample_key = merged_check["_key"].sample(1, random_state=42).iloc[0]
    print(f"\nDetailed verification for _key = {sample_key}")
    subset = merged_check[merged_check["_key"] == sample_key]
    print(subset[['_key', 'class_curated', 'TE class (%)', 'count_direct', 'recomputed_pct']])
    
    # 6️⃣ Check that per-_key sums are 100 ± tolerance
    sum_check = merged_check.groupby("_key")["recomputed_pct"].sum()
    bad_sums = sum_check[(sum_check < 99.9) | (sum_check > 100.1)]
    print(f"\nGroups with total percentage not close to 100%: {len(bad_sums)}")
    
    print("==============================\n")



    for col in data.columns:
        data[col] = data[col].apply(lambda x: tuple(x) if isinstance(x, (list, np.ndarray)) else x)
    
    data = data.drop_duplicates()
    data["positive_selection"] = data["positive_selection"].astype(str).str.lower()
    data["gene_Strand"] = data["gene_Strand"] = "(" + data["gene_Strand"].astype(str) + ")"
    cols = POP_ORDER[::-1] + [c for c in data.columns if c not in POP_ORDER[::-1]]
    data.to_csv(filename_df, sep='\t', index=False, columns=cols)

    upset = UpSet(
        data,
        sort_by= "cardinality", #"degree"
        show_counts=True,
        include_empty_subsets=False,
        sort_categories_by=None, #'cardinality'
        element_size=40,
        max_subset_rank=20 # for PPY to avoid long figure, cut by removing intersections with few intersections
        #min_subset_size=min_subset_size
    )

    # Choose color
    d_class_colors = {"I": "#1f77b4", "II": "#ff7f0e", "Unclassified": "lightgrey"}
    d_order_colors = {
        # Class I (Retrotransposons)
        "LTR": "blue",
        "DIRS": "#1f78b4", #derived from LTR
        "TRIM": "dodgerblue", #derived from LTR

        "LARD": "purple",

        "LINE": "green",
        "SINE": "limegreen",
        
        # Class II (DNA transposons)
        "TIR": "#ff7f0e",
        "MITE": "pink",
        "Helitron": "red",
    
        "Unclassified": "lightgrey"
    }
    
    #hue_order=["Developmental & Signaling", "Immunity", "Metabolism", "Response to Stimulus", "Cellular Process", "Ambiguous", "Molecular Function", "Cell Structure", "Other"],
    d_GO_level2_colors = {
        "Developmental & Signaling": "red",
        "Immunity": "orange",
        #"Metabolism": "#27AE60",
        "Response to Stimulus": "black",
        #"Cellular Process": "#3498DB",
        #"Molecular Function": "#9B59B6",
        #"Cell Structure": "#34495E",
        "Ambiguous": "#BDC3C7", #grey
        "Other": "lightgrey"
    }

    if legend == False:
        if SPECIES == 'PCOM':
            legend = True
        elif SPECIES == 'PPY':
            legend = False

    # Catplot "bar" = side-by-side bars by 'class', value = % (sum of weights)
    upset.add_catplot(
        kind="bar",
        value="TE class (%)",
        hue="class_curated",
        #estimator=sum,
        legend=legend,
        palette=d_class_colors,
        hue_order=["I", "II", "Unclassified"]
    )
    upset.add_catplot(
        kind="bar",
        value="TE order (%)",
        hue="order_curated",
        #estimator=sum,
        legend=legend,
        palette=d_order_colors,
        hue_order=["DIRS", "Helitron", "LINE", "MITE", "SINE", "LARD"], #"LTR",
    )
    upset.add_catplot(
        kind="bar",
        value="TE order (%)",
        hue="order_curated",
        #estimator=sum,
        legend=legend,
        palette=d_order_colors,
        hue_order=["LTR", "TIR", "TRIM", "Unclassified"]
    )
    
    if Type == 'selection':
        # catplot = upset.add_catplot(
        #     kind="bar",
        #     value="GO Level 2 (%)",
        #     hue="GO_Level2_Category",
        #     #estimator=sum,
        #     legend=legend,
        #     palette=d_GO_level2_colors,
        #     hue_order=["Developmental & Signaling", "Immunity", "Response to Stimulus", "Ambiguous"], #, "Other"
        #     #hue_order=["Developmental & Signaling", "Immunity", "Metabolism", "Response to Stimulus", "Cellular Process", "Ambiguous", "Molecular Function", "Cell Structure", "Other"],
        # )
        catplot = upset.add_catplot(
            kind="bar",
            value="Immunity (%)",
            hue="Immunity",
            #estimator=sum,
            legend=legend,
            palette={True: 'orange', False: 'grey'},
            hue_order=[True]
            )
    upset.add_catplot(
        kind="violin",
        value="Rho",
        inner="box",
        cut=0,
        linewidth=0.8
    )
    upset.add_catplot(
        kind="violin",
        value="TE_density",
        inner="box",
        cut=0,
        linewidth=0.8
    )
    # Apply dynamic colors to subsets
    for i in range(1, len(POP_ORDER) + 1):
        for comb in combinations(POP_ORDER, i):
            color = determine_color(comb, SPECIES)
            upset.style_subsets(
                present=list(comb), facecolor=color, edgecolor=color, linewidth=0.5, label=d_upsetplot_color[color]
            )

    if Type == 'selection':
        fig = plt.figure(figsize=(5, 5))
    else:
        fig = plt.figure(figsize=(10, 5))
    subplots = upset.plot(fig=fig)

    axes = upset.plot()
    
    print("\n--- Available keys in axes ---")
    print(axes.keys())

    if Type == 'selection':
        for key in ["extra4"]:
            ax_extra = axes.get(key)
            if ax_extra and ax_extra.get_legend():
                leg = ax_extra.get_legend()
                leg.set_bbox_to_anchor((-0.21, 0.7))
                leg.set_loc("upper right")
                leg.set_title("Immunity")
                break
        
    for key in ["extra3"]: #'matrix', 'shading', 'totals', 'intersections', , "extra3", "extra4", "extra5" "extra1", 
        ax_extra = axes.get(key)
        if ax_extra and ax_extra.get_legend():
            leg = ax_extra.get_legend()
            leg.set_bbox_to_anchor((-0.21, 1.02)) #((1.05, 1))
            leg.set_loc("upper right")
            leg._ncol = 5
            leg.set_title("TE order")
            break
        
    for key in ["extra2"]: #'matrix', 'shading', 'totals', 'intersections', , "extra3", "extra4", "extra5" "extra1", 
        ax_extra = axes.get(key)
        if ax_extra and ax_extra.get_legend():
            leg = ax_extra.get_legend()
            leg.set_bbox_to_anchor((-0.21, 1.2)) #((1.05, 1))
            leg.set_loc("upper right")
            leg._ncol = 5
            leg.set_title("TE order")
            break

    for key in ["extra1"]:
        ax_extra = axes.get(key)
        if ax_extra and ax_extra.get_legend():
            leg = ax_extra.get_legend()
            leg.set_bbox_to_anchor((-0.21, 1))
            leg.set_loc("upper right")
            leg._ncol = 5
            leg.set_title("TE class")
            break
    ax = plt.gca()
    for ax in fig.axes:
        bars = [p for p in ax.patches if isinstance(p, plt.Rectangle) and p.get_width() > p.get_height()]
        if len(bars) == len(POP_ORDER):
            for bar, pop in zip(bars, POP_ORDER[::-1]):
                bar.set_facecolor(POP_COLORS.get(pop, "gray"))
                bar.set_edgecolor(POP_COLORS.get(pop, "gray"))
                bar.set_linewidth(1)
            break

    #plt.suptitle(f"(# of MEs = {len(df_sub)})", fontsize=12)
    plt.savefig(filename, dpi=500, bbox_inches="tight", pad_inches=0.1)
    #plt.savefig(filename + ".svg", dpi=600, bbox_inches="tight", format="svg")
    plt.close(fig)

    print(f"✅ Saved figure: {filename}")






def combine_images(image_paths, titles, output_path, n_lines):
    
    """Combine multiple images (UpSet plots) into a single 2×2 figure."""
    size1 = 5
    size2 = 5
    if n_lines == 2:
        size1 = 7
        size2 = 8
    fig, axes = plt.subplots(n_lines, 2, figsize=(size1, size2))

    for ax, img_path, title in zip(axes.flatten(), image_paths, titles):
        if not os.path.exists(img_path):
            ax.text(0.5, 0.5, "Missing image", ha="center", va="center", color="red")
            ax.axis("off")
            continue
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title, fontsize=size1, fontweight="bold", loc="left")

    plt.tight_layout()
    plt.savefig(output_path, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(f"✅ Combined 2×2 figure saved as {output_path}")






def main():
    
    # ==========================
    # Main analysis
    # ==========================
    

    l_SPECIES = ['PCOM', 'PPY']
    
    for SPECIES in l_SPECIES:
    
        OUT_TSV_FILE = f"../output/{SPECIES}/MEGAnE_summary_{SPECIES}.tsv"
        OUTDIR = f"../Figures/{SPECIES}"
        os.makedirs(OUTDIR, exist_ok=True)
        OUT_FIG_FILE_suppMat = f"../Figures/{SPECIES}/MEGAnE_vs_all_others_genes.png"

        if SPECIES == 'PCOM':
            POP_ORDER = ["cauc", "pyra", "comm_Dessert", "comm_Perry"]
            POP_COLORS = {
                "cauc": "#27BC40",
                "pyra": "#035B03",
                "comm_Dessert": "#442CF4",
                "comm_Perry": "#648FFF",
            }
        elif SPECIES == 'PPY':
            POP_ORDER = ["ussu", "betu", "pash", "Sand_CN-SE", "pyri_JP", "Sand_CN", "Sand_CN-SW", "White"]  
            POP_COLORS = {
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
            
        # --- Load data
        df = pd.read_csv( OUT_TSV_FILE , sep="\t")
    
        # --- Split data into 3 categories
        df_all, df_selected, df_nonselected = split_data( df, POP_ORDER )
        
        # --- Generate individual UpSet plots
        #filename_all = f"{OUTDIR}/upset_genomewide.png" #Genome-wide TE polymorphism
        filename_selection = f"{OUTDIR}/upset_selected.png" #TE polymorphism 2kbp upstream of genes under positive selection
        #filename_noselection = f"{OUTDIR}/upset_nonselected.png" #TE polymorphism 2kbp upstream of other genes

        #filename_all_df = f"{OUTDIR}/upset_genomewide_{SPECIES}.tsv" #Genome-wide TE polymorphism
        filename_selection_df = f"{OUTDIR}/upset_selected_{SPECIES}.tsv" #TE polymorphism 2kbp upstream of genes under positive selection
        #filename_noselection_df = f"{OUTDIR}/upset_nonselected_{SPECIES}.tsv" #TE polymorphism 2kbp upstream of other genes
        
        #make_upset_plot(SPECIES, df_all, POP_ORDER, POP_COLORS, filename_all, True, filename_all_df, None)
        make_upset_plot(SPECIES, df_selected, POP_ORDER, POP_COLORS, filename_selection, False, filename_selection_df, 'selection')
        #make_upset_plot(SPECIES, df_nonselected, POP_ORDER, POP_COLORS, filename_noselection, True, filename_noselection_df, None)


    # For the main Figure 1x2
    filename_selection_PCOM = f"../Figures/PCOM/upset_selected.png"
    filename_selection_PPY = f"../Figures/PPY/upset_selected.png"
    OUT_FIG_FILE_main = f"../Figures/MEGAnE_vs_positive_genes.png"
    combine_images( [filename_selection_PCOM, filename_selection_PPY], ['e', "f"], OUT_FIG_FILE_main, 1)

    # For Supp Mat 2x2
    filename_all_PCOM = f"../Figures/PCOM/upset_genomewide.png"
    filename_noselection_PCOM = f"../Figures/PCOM/upset_nonselected.png"
    filename_all_PPY = f"../Figures/PPY/upset_genomewide.png"
    filename_noselection_PPY = f"../Figures/PPY/upset_nonselected.png"
    l_files = [filename_all_PCOM, filename_noselection_PCOM, filename_all_PPY, filename_noselection_PPY]
    OUT_FIG_FILE_suppMat_concat = f"../Figures/MEGAnE_vs_all_others_genes.png"
    combine_images( l_files, ['a', "b", 'c', "d"], OUT_FIG_FILE_suppMat_concat, 2)


if __name__ == "__main__":
    main()
