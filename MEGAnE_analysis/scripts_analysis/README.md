# MEGAnE downstream analysis pipeline

This folder contains scripts used to analyse MEGAnE TE insertion polymorphisms in pear genomes and their relationships with gene proximity, positive selection, recombination rate, TE density, gene density, immunity-related genes, TFBS enrichment, and genomic enrichment by permutation tests.

---

## Workflow overview

```text
s00_merge_popFiles.py
        │
        ▼
s01_plot_MEGAnE.py
        │
        ▼
s02_create_summary.py
        │
        ▼
MEGAnE_summary_{SPECIES}.tsv
        │
        ├── s03_TEs_vs_genes_upsetplot.py
        ├── s04_compare_genomeEnv.py
        ├── s05_compare_Immunity.py
        ├── s06_launchAME.py
        ├── s08_reformat_files_forRegioneR.py
        │       └── s07_regionR.R
```

---

## Scripts

### s00_merge_popFiles.py

Merge and complete population metadata files.

This script adds missing `Crop_or_Wild` information to the PPY individual metadata table by mapping short population names to species names and then to crop/wild status.

Output:

```text
../data/PPY/s01.individuals_table.txt3
```

---

### s01_plot_MEGAnE.py

Generate main exploratory plots from MEGAnE joint-call VCF files.

This script processes PCOM and PPY MEGAnE outputs and can generate several downstream figures, including:

* sample-wise ME count plots
* TE/SNP PCA plots
* fixed vs polymorphic TE heatmaps
* allele-frequency spectrum plots
* chi-square residual summaries
* combined PCOM/PPY figures

Main inputs:

```text
./output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEI_jointcall.vcf.gz
./output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEA_jointcall.vcf.gz
./data/{SPECIES}/s01.individuals_table.txt
./data/{SPECIES}/TE_annotation_URGI/*classif
```

---

### s02_create_summary.py

Build the main MEGAnE summary table.

For each MEI and MEA variant, the script extracts:

* genomic coordinates
* TE ID and TE type
* carrier samples
* carrier populations
* TE class, order and family
* nearby genes within 2 kb upstream
* positive-selection status
* GO information
* immunity-related annotation
* recombination rate (`Rho`)
* local TE density
* local gene density

Output:

```text
../output/{SPECIES}/MEGAnE_summary_{SPECIES}.tsv
```

This is the central table used by most downstream scripts.

---

### s03_TEs_vs_genes_upsetplot.py

Generate UpSet plots showing TE insertion polymorphism sharing across populations.

The script creates UpSet plots for TE insertions:

* genome-wide
* upstream of genes under positive selection
* upstream of non-selected genes

It also adds side summaries for:

* TE class
* TE order
* recombination rate
* TE density
* immunity-associated genes for selected-gene analyses

Main outputs:

```text
../Figures/PCOM/upset_selected.png
../Figures/PPY/upset_selected.png
../Figures/MEGAnE_vs_positive_genes.png
../Figures/MEGAnE_vs_all_others_genes.png
```

---

### s04_compare_genomeEnv.py

Compare genomic environments of TE insertions.

This script compares three sets of TE insertions:

* whole genome
* upstream of selected genes
* upstream of non-selected genes

For each species, it compares:

* recombination rate (`Rho`)
* TE density
* gene density

Statistical comparisons are performed using Mann–Whitney U tests with Benjamini–Hochberg correction.

Outputs:

```text
../Figures/rho_TE_gene_density_PCOM_PPY.png
../Figures/rho_TE_gene_density_PCOM_PPY.stats.txt
```

---

### s05_compare_Immunity.py

Test whether immunity-associated TIPs are enriched in specific UpSet categories.

The script uses TIP-level information from the UpSet table and compares the proportion of immunity-associated TIPs across population-sharing categories.

It performs:

* descriptive summaries
* global chi-square test
* pairwise Fisher exact tests
* pairwise chi-square tests
* targeted contrasts
* BH-FDR correction

Main input:

```text
../data/PCOM/Table_S14.tsv
```

---

### s06_launchAME.py

Prepare upstream FASTA files and run AME TFBS enrichment analysis.

The script builds 2-kb upstream FASTA files for:

* selected genes with TIPs
* selected genes without TIPs
* non-selected genes

It then runs AME from the MEME suite using a plant JASPAR motif database.

Outputs:

```text
../output/{SPECIES}/TFBS/selected_with_TIP.bed
../output/{SPECIES}/TFBS/selected_with_TIP.fa
../output/{SPECIES}/TFBS/selected_no_TIP.bed
../output/{SPECIES}/TFBS/selected_no_TIP.fa
../output/{SPECIES}/TFBS/non_selected.bed
../output/{SPECIES}/TFBS/non_selected.fa
../output/{SPECIES}/TFBS/gene_counts.txt
../output/{SPECIES}/TFBS/AME/
```

---

### s07_regionR.R

Run permutation tests with regioneR.

This script tests whether selected genes overlap TE insertion polymorphisms more often than expected by chance.

It compares:

* selected-gene upstream regions vs population-specific TIPs
* selected-gene upstream regions vs all TIPs

The permutation test uses:

* `regioneR::permTest`
* `numOverlaps`
* `randomizeRegions`
* 1000 permutations
* one-sided alternative hypothesis: `greater`

Outputs:

```text
../Figures/PCOM/regioneR/pear_permtests_Dessert_vs_ALL.png
../Figures/PCOM/regioneR/pear_permtests_ALL_vs_Dessert.png
```

---

### s08_reformat_files_forRegioneR.py

Prepare input files for regioneR.

This helper script reformats several files into BED-compatible inputs.

It can:

* convert positive-selection TSV files to BED
* create chromosome size files from FASTA
* convert MEGAnE summary tables to population-specific BED files
* create a BED file containing all genes
* generate strand-aware 2-kb upstream BED files

Main outputs include:

```text
*.bed
*.sizes
MEGAnE_summary_{SPECIES}.tsv.{population}.bed
MEGAnE_summary_{SPECIES}.tsv.ALL.bed
```

---

## Main input files

```text
../output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEI_jointcall.vcf.gz
../output/{SPECIES}/jointcall_out_notAdmixed/{SPECIES}_MEA_jointcall.vcf.gz
../data/{SPECIES}/s01.individuals_table.txt
../data/{SPECIES}/genes/positive_selection_summary_table.*
../data/{SPECIES}/genes/*.gff
../data/{SPECIES}/TE_annotation_URGI/*classif
../data/{SPECIES}/TE_annotation_URGI/*.gff3
../data/{SPECIES}/Rho/*.bed
```

---

## Main output files

```text
../output/{SPECIES}/MEGAnE_summary_{SPECIES}.tsv
../Figures/{SPECIES}/upset_selected.png
../Figures/MEGAnE_vs_positive_genes.png
../Figures/MEGAnE_vs_all_others_genes.png
../Figures/rho_TE_gene_density_PCOM_PPY.png
../Figures/rho_TE_gene_density_PCOM_PPY.stats.txt
../output/{SPECIES}/TFBS/
../Figures/PCOM/regioneR/
```

---

## Software requirements

Python:

```text
pandas
numpy
matplotlib
seaborn
scipy
statsmodels
pyranges
goatools
upsetplot
Biopython
```

External tools:

```text
bedtools
samtools
MEME suite
```

R:

```text
regioneR
GenomicRanges
```

---

## Notes

The central file of the pipeline is:

```text
MEGAnE_summary_{SPECIES}.tsv
```

It integrates TE insertion calls, population information, gene proximity, positive-selection status, GO annotation, immunity flags, recombination rate, TE density and gene density.

Most downstream analyses are based on this summary table.
