#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re, warnings, gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from matplotlib import gridspec
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

warnings.filterwarnings("ignore", category=UserWarning)
sns.set_style("whitegrid")

# ===================== Loaders =====================

def load_classification(classif_file):
    df = pd.read_csv(classif_file, sep="\t")
    for col in ["Seq_name", "class", "order"]:
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {classif_file}")
    df["class"] = df["class"].fillna("Unclassified").astype(str)
    df["order"] = df["order"].fillna("Unclassified").astype(str)

    length_col = None
    for cand in ["length","Length","consensus_length","Consensus_length","Seq_length","size"]:
        if cand in df.columns:
            length_col = cand
            break
    return df, length_col

def load_consensus_lengths_from_fasta(fasta_file):
    cons_len = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        cons_len[rec.id] = len(rec.seq)
    return cons_len

# ===================== VCF parsing =====================

def parse_vcf_summary(vcf_file, class_dict, order_dict, type_name):
    rows = []
    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10: continue

            chrom, _, _, _, _, _, _, info, fmt = cols[:9]

            # SVLEN
            svlen = None
            for kv in info.split(";"):
                if kv.startswith("SVLEN="):
                    try: svlen = abs(int(kv.split("=")[1]))
                    except: svlen = None

            # MEI= consensus list
            mei_ids = []
            for kv in info.split(";"):
                if kv.startswith("MEI="):
                    mei_ids = kv.split("=")[1].split("|")
                    break
            nTE = len(mei_ids) if mei_ids else 0

            # quick AF proxy (not used for AFS)
            fmt_fields = fmt.split(":")
            if "GT" not in fmt_fields: continue
            gt_idx = fmt_fields.index("GT")
            alt_count, called = 0, 0
            for g in cols[9:]:
                parts = g.split(":")
                if gt_idx >= len(parts): continue
                gt = parts[gt_idx]
                if gt in (".","./.",".|."): continue
                for a in gt.replace("|","/").split("/"):
                    if a == ".": continue
                    called += 1
                    a = int(a)
                    if type_name == "MEA":
                        if a == 0: alt_count += 1  # MEA: 0 = TE
                    else:
                        if a == 1: alt_count += 1  # MEI: 1 = TE
            if called == 0: continue
            freq = alt_count / called

            for mei in (mei_ids if mei_ids else ["NA"]):
                te_class = str(class_dict.get(mei,"Unclassified")).split("|")[0]
                te_order = str(order_dict.get(mei,"Unclassified")).split("|")[0]
                rows.append({
                    "chrom": chrom, "type": type_name,
                    "class": te_class, "order": te_order,
                    "svlen": svlen, "freq": freq, "nTE": nTE, "seq_name": mei
                })
    return pd.DataFrame(rows)

# ===================== VCF parsing =====================

def parse_vcf_genotypes_with_orders(vcf_file, order_dict, type_name):
    """
    Retourne:
      - G: DataFrame (samples × variants) dosage ALT (0/1/2) par individu
      - var_orders: list[set] des orders rencontrés pour chaque variante
    """
    samples = []
    variant_genos = []
    variant_orders = []
    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                hdr = line.strip().split("\t")
                samples = hdr[9:]
                continue
            cols = line.strip().split("\t")
            if len(cols) < 10: continue
            info = cols[7]; fmt = cols[8]
            genotypes = cols[9:]

            # orders from MEI list
            mei_ids = []
            for kv in info.split(";"):
                if kv.startswith("MEI="):
                    mei_ids = kv.split("=")[1].split("|")
                    break
            ords = set(str(order_dict.get(mei,"Unclassified")).split("|")[0] for mei in mei_ids) if mei_ids else set(["Unclassified"])
            variant_orders.append(ords)

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
                row.append(sum(vals) if len(vals)==2 else np.nan)  # ALT dosage
            variant_genos.append(row)

    if not variant_genos:
        G = pd.DataFrame(index=[s.split(".")[2] if len(s.split("."))>=3 else s for s in samples])
        return G, variant_orders

    G = pd.DataFrame(variant_genos).T  # samples × variants
    G.columns = [f"var{i}" for i in range(G.shape[1])]
    G.index = [s.split(".")[2] if len(s.split("."))>=3 else s for s in samples]
    return G, variant_orders

# ===================== Pop & palettes =====================

def load_population_table(pop_file):
    df = pd.read_csv(pop_file, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    if "#ID_vcf" in df.columns:
        df = df.rename(columns={"#ID_vcf":"ID_vcf"})
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    for col in ("Depth","ReadsAvgLength"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.set_index("ID_vcf")

def build_population_palette(df_pop):
    pops = sorted(pd.Series(df_pop["Population"].unique(), dtype=str))
    base = sns.color_palette("tab10", n_colors=max(3, len(pops)))
    return {p: base[i] for i,p in enumerate(pops)}

# ===================== PCA helpers =====================

def pca_coords_with_meta(G, df_pop):
    meta = df_pop.copy(); meta.index = meta.index.astype(str)
    X = G.copy(); X.index = X.index.astype(str)
    X = X.merge(meta[["Population","Depth","ReadsAvgLength"]],
                left_index=True, right_index=True, how="left")
    X["Population"] = X["Population"].fillna("Unknown")
    M = X.drop(columns=["Population","Depth","ReadsAvgLength"]).values
    if M.size==0 or M.shape[1] < 2: return None, None
    M = SimpleImputer(strategy="mean").fit_transform(M)
    if M.shape[1] < 2: return None, None
    p = PCA(n_components=2).fit(M)
    coords = p.transform(M)
    dfp = pd.DataFrame(coords, index=X.index, columns=["PC1","PC2"])
    dfp["Population"] = X["Population"].values
    dfp["Depth"] = X["Depth"].values
    dfp["ReadsAvgLength"] = X["ReadsAvgLength"].values
    exp = p.explained_variance_ratio_*100
    return dfp, exp

def scatter_pca(ax, dfp, exp, color_by, pop_palette=None, cbar_label=None, title_prefix=""):
    if dfp is None:
        ax.text(0.5,0.5,"PCA not available",ha="center",va="center"); ax.set_axis_off(); return

    if color_by == "Population":
        pal = dict(pop_palette or {}); pal.setdefault("Unknown",(0.8,0.8,0.8))
        sns.scatterplot(data=dfp, x="PC1", y="PC2", hue="Population",
                        palette=pal, s=42, edgecolor="k", linewidth=0.3, ax=ax, legend=True)
        ax.legend(loc="best", fontsize=7, frameon=True)
        ax.set_title(f"{title_prefix}Population\n(PC1 {exp[0]:.1f}%, PC2 {exp[1]:.1f}%)", fontsize=10)
    else:
        vals = pd.to_numeric(dfp[color_by], errors="coerce")
        invalid = vals.isna() | (vals == -1)
        if invalid.sum() > 0:
            ax.scatter(dfp.loc[invalid,"PC1"], dfp.loc[invalid,"PC2"],
                       c="lightgray", s=42, edgecolor="k", linewidth=0.3, label=f"{color_by} = -1")
        valid = ~invalid
        if valid.sum() > 0:
            vmin, vmax = float(vals[valid].min()), float(vals[valid].max())
            norm = Normalize(vmin=vmin, vmax=vmax)
            ax.scatter(dfp.loc[valid,"PC1"], dfp.loc[valid,"PC2"],
                       c=vals[valid], s=42, edgecolor="k", linewidth=0.3, cmap="viridis", norm=norm)
            cbar = plt.colorbar(ScalarMappable(norm=norm, cmap="viridis"), ax=ax, shrink=0.75, pad=0.01)
            cbar.set_label(cbar_label or color_by, fontsize=8)
        if invalid.sum() > 0:
            ax.legend(loc="best", fontsize=7, frameon=True)
        ax.set_title(f"{title_prefix}{color_by}\n(PC1 {exp[0]:.1f}%, PC2 {exp[1]:.1f}%)", fontsize=10)

    ax.set_xlabel("PC1", fontsize=9); ax.set_ylabel("PC2", fontsize=9)
    ax.tick_params(axis='both', labelsize=8)

# ===================== AFS helpers (AF = #TE alleles / 2N) =====================

def _dosage_te(G_values, type_name):
    g = np.array(G_values, dtype=float)
    return (2.0 - g) if type_name == "MEA" else g

def classify_variants(G, df_pop, type_name):
    meta = df_pop[["Population"]].copy(); meta.index = meta.index.astype(str)
    G = G.copy(); G.index = G.index.astype(str)
    pops_series = meta["Population"].reindex(G.index).fillna("Unknown")
    d = _dosage_te(G.values, type_name)
    present = (d > 0).astype(int)
    known = [p for p in pops_series.unique().tolist() if p != "Unknown"]
    pop_specific = {p: [] for p in known}; admixed = []
    if not known:
        any_present = present.sum(axis=0) > 0
        admixed = np.where(any_present)[0].tolist()
        return pop_specific, admixed, pops_series
    counts = np.vstack([(present[pops_series.values==p, :]).sum(axis=0) for p in known])
    support = (counts > 0).sum(axis=0)
    for j in range(G.shape[1]):
        if support[j] == 1:
            p_idx = np.where(counts[:, j] > 0)[0][0]
            pop_specific[known[p_idx]].append(j)
        elif support[j] >= 2:
            admixed.append(j)
    return pop_specific, admixed, pops_series

def build_afs_panels_with_stats(G, df_pop, type_name, top_k=4):
    G = G.copy(); G.index = G.index.astype(str)
    pop_specific, admixed_idx, pops_series = classify_variants(G, df_pop, type_name)
    total_loci = G.shape[1]; total_samples = G.shape[0]
    panels = {}

    def _af_from_matrix(d_te):
        d = np.array(d_te, dtype=float)
        mask = ~np.isnan(d)
        n_called = mask.sum(axis=0)
        d[~mask] = 0.0
        te_sum = d.sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            af = np.where(n_called > 0, te_sum / (2.0 * n_called), np.nan)
        return af[~np.isnan(af)]

    d_all = _dosage_te(G.values, type_name)
    pop_counts = {p: len(ix) for p, ix in pop_specific.items()}
    top_pops = [p for p,_ in sorted(pop_counts.items(), key=lambda kv: kv[1], reverse=True)[:top_k]]

    for p in top_pops:
        idxs = pop_specific[p]; nL = len(idxs)
        nS = int((pops_series.values == p).sum())
        if nL > 0 and nS > 0:
            cols = [G.columns[j] for j in idxs]
            sub = d_all[pops_series.values == p][:, [G.columns.get_loc(c) for c in cols]]
            af = _af_from_matrix(sub)
        else:
            af = np.array([])
        panels[f"{p}-specific"] = {"af": af, "nLocus": nL, "nSamples": nS}

    af_all = _af_from_matrix(d_all)
    panels["All populations"] = {"af": af_all, "nLocus": total_loci, "nSamples": total_samples}

    if len(admixed_idx) > 0:
        cols = [G.columns[j] for j in admixed_idx]
        sub = d_all[:, [G.columns.get_loc(c) for c in cols]]
        af_adm = _af_from_matrix(sub)
        panels["Admixed"] = {"af": af_adm, "nLocus": len(admixed_idx), "nSamples": total_samples}
    else:
        panels["Admixed"] = {"af": np.array([]), "nLocus": 0, "nSamples": total_samples}

    return panels, total_loci

def plot_afs_3x2_with_titles(panels, total_loci, axes, pop_palette, bins=None):
    if bins is None: bins = np.linspace(0,1,51)
    colors = dict(pop_palette)
    colors.setdefault("All populations",(0.25,0.25,0.25))
    colors.setdefault("Admixed",(0.6,0.6,0.6))
    specific_keys = [k for k in panels.keys() if k.endswith("-specific")]
    ordered = specific_keys[:4] + ["All populations","Admixed"]
    for idx, key in enumerate(ordered):
        ax = axes[idx]
        e = panels.get(key, {"af": np.array([]), "nLocus": 0, "nSamples": 0})
        af = e["af"]; nL = int(e["nLocus"]); nS = int(e["nSamples"])
        pctL = (100.0*nL/total_loci) if total_loci>0 else 0.0
        if af is None or len(af)==0 or nL==0:
            ax.text(0.5,0.5,f"{key}\n(nLocus=0, 0.0% ; nSamples={nS})",
                    ha="center", va="center", fontsize=9)
            ax.set_axis_off()
        else:
            col = colors.get(key.split("-")[0], colors.get(key, None))
            ax.hist(af, bins=bins, edgecolor="black", color=col)
            ax.set_title(f"{key}\n(nLocus={nL:,}, {pctL:.1f}% ; nSamples={nS:,})", fontsize=9)
        row, col_i = idx//2, idx%2
        if row < 2:
            ax.set_xlabel(""); ax.set_xticklabels([])
        else:
            ax.set_xlabel("Allele frequency (TE)", fontsize=8)
        if col_i == 0:
            ax.set_ylabel("Count", fontsize=8)
        else:
            ax.set_ylabel(""); ax.set_yticklabels([])
        ax.tick_params(axis='both', labelsize=8)
    for j in range(len(ordered), len(axes)): axes[j].axis("off")

# ===================== MAIN =====================

def main():
    ap = argparse.ArgumentParser(description="MEGAnE summary by TE order")
    ap.add_argument("--order", default="ALL", help="TE order filter (e.g., ALL, LTR, TIR, LINE, ...)")
    v = ap.parse_args()
    order_filter = v.order.upper()

    # Paths
    vcf_mea = "../output/PCOM/jointcall_out/PCOM_MEA_jointcall.vcf.gz"
    vcf_mei = "../output/PCOM/jointcall_out/PCOM_MEI_jointcall.vcf.gz"
    classif_file = "../data/PCOM/TE_annotation_URGI/PCOM_TEdenovoGr_sim_denovoLibTEs_PC_filtered.classif"
    fasta_file = "../data/PCOM/assembly/PyrusCommunis_BartlettDHv2.0.fasta.reformat"
    pop_file = "../data/PCOM/s01.individuals_table.txt"
    outdir = "/scratch/ss20440/MEGAnE_Yuqi/Figures"
    os.makedirs(outdir, exist_ok=True)

    # Output name
    if order_filter == "ALL":
        outfile = os.path.join(outdir, "MEGAnE_summary.png")
    else:
        outfile = os.path.join(outdir, f"MEGAnE_summary_{order_filter}.png")

    # Classif & dicts
    dfc, classif_len_col = load_classification(classif_file)
    cons_len_fa = load_consensus_lengths_from_fasta(fasta_file)
    class_dict = dict(zip(dfc["Seq_name"], dfc["class"]))
    order_dict = dict(zip(dfc["Seq_name"], dfc["order"]))

    # Résumés (a,b)
    df_mea = parse_vcf_summary(vcf_mea, class_dict, order_dict, "MEA")
    df_mei = parse_vcf_summary(vcf_mei, class_dict, order_dict, "MEI")

    if order_filter != "ALL":
        df_mea = df_mea[df_mea["order"] == order_filter]
        df_mei = df_mei[df_mei["order"] == order_filter]

    df = pd.concat([df_mea, df_mei], ignore_index=True)
    type_palette = {"MEA":"steelblue","MEI":"orange","Consensus":"grey"}

    # ===== Figure layout (espaces comme tu veux) =====
    fig = plt.figure(figsize=(26, 18))
    gs  = gridspec.GridSpec(2, 3, figure=fig, wspace=0.1, hspace=0.15)

    # (a)+(b) (colonne gauche)
    gs_ab = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0,0], hspace=0.18)
    ax_a  = fig.add_subplot(gs_ab[0,0])
    ax_b  = fig.add_subplot(gs_ab[1,0])

    # (c) AFS MEA 3×2 (haut centre)
    gs_c = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs[0,1], wspace=0.04, hspace=0.3)
    ax_c = [fig.add_subplot(gs_c[i,j]) for i in range(3) for j in range(2)]

    # (d) AFS MEI 3×2 (haut droite)
    gs_d = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs[0,2], wspace=0.04, hspace=0.3)
    ax_d = [fig.add_subplot(gs_d[i,j]) for i in range(3) for j in range(2)]

    # (e) PCA bloc centre bas -> **titres MEA**
    gs_e = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1,0], wspace=0.1, hspace=0.3)
    ax_e = [fig.add_subplot(gs_e[i,j]) for i in range(2) for j in range(2)]

    # (f) PCA bloc gauche bas  -> **titres MEI**
    gs_f = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1,1], wspace=0.1, hspace=0.3)
    ax_f = [fig.add_subplot(gs_f[i,j]) for i in range(2) for j in range(2)]

    # (g) Violin lengths par order (MEA/MEI/Consensus)
    ax_g = fig.add_subplot(gs[1,2])

    # ---------- (a) barplot counts ----------
    def natural_key(s): return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]
    if df.empty:
        ax_a.text(0.5,0.5,"No variants",ha="center",va="center"); ax_a.axis("off")
        ax_b.text(0.5,0.5,"No variants",ha="center",va="center"); ax_b.axis("off")
    else:
        chrom_counts = df.groupby(["chrom","type"]).size().reset_index(name="count")
        chrom_order = sorted(chrom_counts["chrom"].unique(), key=natural_key)
        sns.barplot(data=chrom_counts, x="chrom", y="count", hue="type",
                    order=chrom_order, palette=type_palette, ax=ax_a)
        ax_a.set_ylabel("Number of variants"); ax_a.set_xlabel("Chromosome")
        ax_a.set_xticklabels([re.sub(r"\D*","",t.get_text()) or t.get_text() for t in ax_a.get_xticklabels()], rotation=0)
        ax_a.text(0,1.02,"a", transform=ax_a.transAxes, fontsize=16, fontweight="bold")

        df["xlab"] = df["chrom"].apply(lambda s: re.sub(r"\D*","", s) or s)
        order_x = [re.sub(r"\D*","",c) or c for c in chrom_order]
        sns.violinplot(data=df, x="xlab", y="nTE", hue="type",
                       order=order_x, hue_order=["MEA","MEI"],
                       palette=type_palette, inner="quartile", cut=0, scale="width", ax=ax_b)
        ax_b.set_ylabel("Number of TE per variant"); ax_b.set_xlabel("Chromosome")
        ax_b.set_xticklabels(order_x, rotation=0)
        ax_b.text(0,1.02,"b", transform=ax_b.transAxes, fontsize=16, fontweight="bold")

    # ---------- PCA/AFS matrices par order ----------
    df_pop = load_population_table(pop_file)
    pop_palette = build_population_palette(df_pop)

    G_mea_all, var_orders_mea = parse_vcf_genotypes_with_orders(vcf_mea, order_dict, "MEA")
    G_mei_all, var_orders_mei = parse_vcf_genotypes_with_orders(vcf_mei, order_dict, "MEI")

    if order_filter == "ALL":
        is_selected_mea = pd.Series([True]*G_mea_all.shape[1], index=G_mea_all.columns) if G_mea_all.shape[1]>0 else pd.Series(dtype=bool)
        is_selected_mei = pd.Series([True]*G_mei_all.shape[1], index=G_mei_all.columns) if G_mei_all.shape[1]>0 else pd.Series(dtype=bool)
    else:
        is_selected_mea = pd.Series([(order_filter in s) for s in var_orders_mea], index=G_mea_all.columns) if G_mea_all.shape[1]>0 else pd.Series(dtype=bool)
        is_selected_mei = pd.Series([(order_filter in s) for s in var_orders_mei], index=G_mei_all.columns) if G_mei_all.shape[1]>0 else pd.Series(dtype=bool)

    G_mea = G_mea_all.loc[:, is_selected_mea] if not G_mea_all.empty else G_mea_all
    G_mei = G_mei_all.loc[:, is_selected_mei] if not G_mei_all.empty else G_mei_all

    # ---------- (c) AFS MEA ----------
    if (not G_mea.empty) and (G_mea.shape[1] > 0):
        panels_mea, total_mea = build_afs_panels_with_stats(G_mea, df_pop, "MEA", top_k=4)
        plot_afs_3x2_with_titles(panels_mea, total_mea, ax_c, pop_palette, bins=np.linspace(0,1,51))
    else:
        for a in ax_c: a.text(0.5,0.5,"MEA: no data",ha="center",va="center"); a.set_axis_off()
    ax_c[0].text(0,1.18,"c", transform=ax_c[0].transAxes, fontsize=16, fontweight="bold")

    # ---------- (d) AFS MEI ----------
    if (not G_mei.empty) and (G_mei.shape[1] > 0):
        panels_mei, total_mei = build_afs_panels_with_stats(G_mei, df_pop, "MEI", top_k=4)
        plot_afs_3x2_with_titles(panels_mei, total_mei, ax_d, pop_palette, bins=np.linspace(0,1,51))
    else:
        for a in ax_d: a.text(0.5,0.5,"MEI: no data",ha="center",va="center"); a.set_axis_off()
    ax_d[0].text(0,1.18,"d", transform=ax_d[0].transAxes, fontsize=16, fontweight="bold")

    # ---------- (e) PCA bloc centre (TITRES = **MEA** après inversion)
    dfp_mea, exp_mea = pca_coords_with_meta(G_mea, df_pop) if ((not G_mea.empty) and G_mea.shape[1]>1) else (None, None)
    ax_e[0].text(0,1.06,"e", transform=ax_e[0].transAxes, fontsize=16, fontweight="bold")
    scatter_pca(ax_e[0], dfp_mea, exp_mea, color_by="Population", pop_palette=pop_palette, title_prefix="MEA - ")
    scatter_pca(ax_e[1], dfp_mea, exp_mea, color_by="Depth", cbar_label="Depth", title_prefix="MEA - ")
    scatter_pca(ax_e[2], dfp_mea, exp_mea, color_by="ReadsAvgLength", cbar_label="ReadsAvgLength", title_prefix="MEA - ")
    ax_e[3].axis("off")
    
    # ---------- (f) PCA bloc gauche (TITRES = **MEI** après inversion)
    dfp_mei, exp_mei = pca_coords_with_meta(G_mei, df_pop) if ((not G_mei.empty) and G_mei.shape[1]>1) else (None, None)
    ax_f[0].text(0,1.06,"f", transform=ax_f[0].transAxes, fontsize=16, fontweight="bold")
    scatter_pca(ax_f[0], dfp_mei, exp_mei, color_by="Population", pop_palette=pop_palette, title_prefix="MEI - ")
    scatter_pca(ax_f[1], dfp_mei, exp_mei, color_by="Depth", cbar_label="Depth", title_prefix="MEI - ")
    scatter_pca(ax_f[2], dfp_mei, exp_mei, color_by="ReadsAvgLength", cbar_label="ReadsAvgLength", title_prefix="MEI - ")
    ax_f[3].axis("off")

    # ---------- (g) Violin par order (MEA/MEI + Consensus si ALL sinon filtré par order)
    if classif_len_col is not None:
        df_cons = dfc[["Seq_name","class","order",classif_len_col]].rename(columns={classif_len_col:"svlen"}).copy()
    else:
        df_cons = dfc[["Seq_name","class","order"]].copy()
        cons_len_map = load_consensus_lengths_from_fasta(fasta_file)
        df_cons["svlen"] = df_cons["Seq_name"].map(cons_len_map).astype(float)
    if order_filter != "ALL":
        df_cons = df_cons[df_cons["order"] == order_filter]
    df_cons["type"] = "Consensus"

    df_variants = pd.concat([df_mea.assign(type="MEA"),
                             df_mei.assign(type="MEI")], ignore_index=True, sort=False)
    df_g = pd.concat([df_variants[["order","svlen","type"]], df_cons[["order","svlen","type"]]],
                     ignore_index=True)

    # ordre par class I/II/Unclassified pour la lisibilité
    order_to_class = dict(zip(dfc["order"], dfc["class"]))
    grouped_orders = {}
    for o in df_g["order"].dropna().unique():
        grouped_orders.setdefault(order_to_class.get(o, "Unclassified"), []).append(o)
    classI_orders  = sorted(grouped_orders.get("I", []))
    classII_orders = sorted(grouped_orders.get("II", []))
    unclass_orders = sorted(grouped_orders.get("Unclassified", []))
    custom_order   = classI_orders + classII_orders + unclass_orders

    if df_g.dropna(subset=["svlen"]).empty:
        ax_g.text(0.5,0.5,"No length data",ha="center",va="center"); ax_g.axis("off")
    else:
        sns.violinplot(data=df_g.dropna(subset=["svlen"]), x="order", y="svlen", hue="type",
                       order=custom_order, hue_order=["MEA","MEI","Consensus"],
                       inner="quartile", scale="width", cut=0, palette=type_palette, ax=ax_g)
        ax_g.set_ylabel("Sequence length")
        ax_g.tick_params(axis='x', rotation=90, labelsize=9)
        ax_g.text(0,1.02,"g", transform=ax_g.transAxes, fontsize=16, fontweight="bold")
        # séparateurs de classes
        ymax = ax_g.get_ylim()[1]; ypos = ymax*1.01
        if classI_orders:
            ax_g.text((len(classI_orders)-1)/2, ypos, "Class I", ha="center", va="bottom", fontsize=11, fontweight="bold")
        if classII_orders:
            idx_sep1 = len(classI_orders)-0.5; ax_g.axvline(idx_sep1, color="black", linestyle="--", linewidth=1)
            ax_g.text(len(classI_orders)+(len(classII_orders)-1)/2, ypos, "Class II", ha="center", va="bottom", fontsize=11, fontweight="bold")
        if unclass_orders:
            idx_sep2 = len(classI_orders)+len(classII_orders)-0.5; ax_g.axvline(idx_sep2, color="black", linestyle="--", linewidth=1)
            ax_g.text(len(classI_orders)+len(classII_orders)+(len(unclass_orders)-1)/2, ypos, "Unclassified", ha="center", va="bottom", fontsize=11, fontweight="bold")

    plt.tight_layout()
    plt.savefig(outfile, dpi=350, bbox_inches="tight")
    plt.close()
    print(f"✅ Summary figure saved: {outfile}")

if __name__ == "__main__":
    main()
