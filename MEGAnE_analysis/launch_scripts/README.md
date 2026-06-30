# MEGAnE and REPET preprocessing pipeline

This repository contains helper scripts used to prepare REPET annotations for MEGAnE, perform mobile element discovery, generate joint VCFs, and summarize transposable element (TE) composition and allele-frequency spectra.

---

## Workflow

```
REPET consensus FASTA
          │
          ▼
s01.REPET2RM_fasta.py
          │
          ▼
RepeatMasker-formatted TE library (.fa)

REPET GFF3
          │
          ▼
s02.REPET2RM_gff.py
          │
          ▼
RepeatMasker .out annotation

Genome FASTA
          │
          ▼
s03.launchMEGAnE.sh
          │
          ▼
Genome indexing
(k-mer database + BLAST database)

CRAM files
          │
          ▼
s04.launchMEGAnE_part2_*.sh
          │
          ▼
MEGAnE call_genotype
(one directory per sample)

Sample outputs
          │
          ▼
s05.launchMEGAnE_part3_*.sh
          │
          ▼
Joint calling
(MEI + MEA VCFs)

Joint VCFs
          │
          ├── s06.plot_MEGAnE_*.py
          ├── s07.plot_MEGAnE_AFS.py
          └── s07.plot_MEGAnE_AFS_perPop.py
```

---

# Scripts

## s00.plot_REPET.py

Generate summary statistics for REPET annotations.

### Input

- Genome FASTA index (.fai)
- REPET GFF3 annotation
- REPET classification table

### Output

**REPET_summary.png**

Panels:

**(a)** TE composition, including the *No TE coverage* category

**(b)** TE length distribution

**(c)** TE sequence identity distribution

**(d)** Local TE density along each chromosome

TE density is calculated as the number of TE-covered base pairs within non-overlapping 100-kbp windows.

To speed up repeated executions, TE densities are pre-computed once and stored in:

```
TEdensity_100kb.tsv
```

The script automatically reloads this file if it already exists.

---

## s01.REPET2RM_fasta.py

Convert a REPET consensus library into a RepeatMasker-compatible FASTA library.

Header format:

```
>FamilyName    Class/Subclass    REPET
```

The script automatically:

- parses REPET classifications
- normalizes TE classes
- converts REPET nomenclature into RepeatMasker format

---

## s02.REPET2RM_gff.py

Convert REPET TE annotation (GFF3) into a RepeatMasker `.out` file.

The script reconstructs:

- TE class/family
- divergence
- genomic coordinates
- RepeatMasker-compatible formatting

This file is required by MEGAnE (`-repout` option).

---

## s03.launchMEGAnE.sh

Prepare reference files required by MEGAnE.

The script performs:

- genome k-mer database construction
- BLAST database generation

These resources are generated once per reference genome.

---

## s04.launchMEGAnE_part2_*.sh

Run **MEGAnE call_genotype** independently for each CRAM file.

Main steps:

- detect CRAM files automatically
- extract sample names
- skip already processed samples
- execute MEGAnE inside a Singularity container

Each sample receives its own output directory.

---

## s05.launchMEGAnE_part3_*.sh

Perform MEGAnE joint calling.

Steps:

1. collect all sample output directories
2. generate `dirlist.txt`
3. merge non-reference insertions (MEI)
4. merge reference polymorphisms (MEA)

Outputs:

- `*_MEI_jointcall.vcf.gz`
- `*_MEA_jointcall.vcf.gz`

---

## s06.plot_MEGAnE.py

Generate a global summary figure for MEGAnE results.

The figure includes:

- chromosome-wise variant counts
- number of TE consensuses per variant
- allele-frequency spectra
- PCA
- TE length distributions

Results can be produced either:

- genome-wide
- for a single TE order

---

## s07.plot_MEGAnE_AFS.py

Generate allele-frequency spectra grouped by TE categories.

Panels include:

- LTR
- LINE + SINE
- DIRS
- Class I (other)
- TIR
- MITE
- Helitron
- Other Class II + Unclassified

Variants are separated into:

- Wild-specific
- Cultivated-specific
- Shared (Both)

MEI and MEA variants are analysed together.

---

## s07.plot_MEGAnE_AFS_perPop.py

Generate allele-frequency spectra for individual populations.

The script:

- automatically selects the four largest populations
- computes allele frequencies independently for each population
- groups variants into four simplified TE categories:

- LTR
- Class I non-LTR
- Class II
- Unclassified

---

# Software

- MEGAnE
- REPET
- RepeatMasker
- BLAST+
- samtools
- Python
- pandas
- NumPy
- matplotlib
- seaborn

---

# HPC requirements

The pipeline is designed for SLURM clusters.

Execution relies on:

- Singularity
- BLAST+
- RepeatMasker
- samtools

Most computationally intensive steps are parallelized using multiple CPU cores.

---

# Output structure

```
output/
├── megane_kmer_set/
├── MEGAnE_result/
├── jointcall_out/
│   ├── *_MEI_jointcall.vcf.gz
│   └── *_MEA_jointcall.vcf.gz
└── Figures/
    ├── REPET_summary.png
    ├── MEGAnE_summary.png
    ├── MEGAnE_AFS.png
    └── MEGAnE_AFS_byPopulation.png
```