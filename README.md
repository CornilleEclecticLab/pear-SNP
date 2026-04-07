# Contrasting genomic routes to domestication in Occidental and Oriental pears
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2025.11.17.687327-blue)](https://doi.org/10.1101/2025.11.17.687327)


## The workflow for bioinformatic analysis used in this study (SNPs part)

Fig S1 The workflow of the study

```mermaid
%%{init: {'theme':'forest',
  "themeVariables": {
    "fontSize": "14px",
     "fontFamily": "Arial"
    }
}}%%

flowchart TD

ROW_DATA[["`**Raw Reads Data**
    (fastq)
    Published data from previous studies and newly generated data from this study`"]]

    --> QUALITY_CONTROL{{fastp, fastQC, MultiQC}}
    --  Remove redundant samples
        Remove low quality samples
        Remove&nbsp;samples&nbsp;with&nbsp;sequencing&nbsp;depth&nbsp;below&nbsp;9X
    
    --> CLEAN_DATA[["`**Clean Reads Data**
            (fastq)
            *N* = 674`"]]
    QUALITY_CONTROL
      ~~~
    gatk

    CLEAN_DATA
    --> gatk{{BWA2, Picard, GATK4, bcftools, Plink2 kinship}}

    --- L1[\Level 1 and 2 Filters/]

    --  Remove low quality sites
        Remove&nbsp;high&nbsp;SNP&nbsp;missing&nbsp;rate&nbsp;(>40%)&nbsp;samples
        Remove&nbsp;genetically&nbsp;closely&nbsp;related&nbsp;samples&nbsp;(KING&nbsp;kinship&nbsp;>&nbsp;0.354)
        Remove&nbsp;samples&nbsp;with&nbsp;unclear&nbsp;passport&nbsp;information
    
    --> FILTERED_GENOTYPE[["`**Filtered Genotype Data**
            (vcf)
            *N* = 396
            Ref. *Pyrus pyrifolia*
            and Ref. *Pyrus communis*`"]]
    
    -- -Only Invariant sites
    --> FILTERED_INVARIANT[["`**Filtered Invariants**
            *N* = 229`"]]

REF_GENOME_W[["`**Reference Occidental Pear Genome**
            (fasta)
            *Pyrus communis*`"]]
    --> gatk

REF_GENOME_E[["`**Reference Oriental Pear Genome**
            (fasta)
            *Pyrus pyrifolia*`"]]
    --> gatk

    
FILTERED_GENOTYPE
    --  -Only biallelic SNPs
    
    --> FILTERED_VARIANT[["`**Filtered Variants**
            *N* = 396`"]]

    --- L3[\Level 3 Filters/]
    -- Remove&nbsp;nonsynonymous&nbsp;substitution&nbsp;sites
       MAF&nbsp;<=&nbsp;0.05
       LD thinning

    --> UNLINKED_VARIANT[["`**Unlinked Variants**
            *N* = 396`"]]



UNLINKED_VARIANT
    --> splittree{{Splitstree}}
    --> tree([Neighbor-net tree])


UNLINKED_VARIANT
    --> plinkpca{{Plink PCA}}
    --> pca([PCA])


UNLINKED_VARIANT
    --> faststrcuture{{fastStructure}}
    --> structure([Population Structure])


REF_GENOME_BOTH[["`**Both Reference Genomes**
            (fasta)
            *Pyrus pyrifolia*
            and *Pyrus communis*`"]]
    --> genemap{{genemap, CentIER}}
    -.  Mask low mappability and centromeric regions
    .-> FINAL_VARIANT

REF_GENOME_BOTH
    --> sift_anno{{SIFT4G annotator}}
    -.-> FINAL_VARIANT

FILTERED_GENOTYPE
    --> FINAL_VARIANT
structure
    -.  Remove&nbsp;admixed&nbsp;samples&nbsp;and&nbsp;keep&nbsp;229&nbsp;pure&nbsp;samples&nbsp;belonging&nbsp;to&nbsp;12&nbsp;populations
        Identification&nbsp;of&nbsp;cultivated&nbsp;and&nbsp;wild&nbsp;populations
    .-> FINAL_UNLINKED_VARIANT[["`**Final Unlinked Variants**
            *N* = 229`"]]



structure
    -.  Remove&nbsp;admixed&nbsp;samples&nbsp;and&nbsp;keep&nbsp;229&nbsp;pure&nbsp;samples&nbsp;belonging&nbsp;to&nbsp;12&nbsp;populations
        Identification&nbsp;of&nbsp;cultivated&nbsp;and&nbsp;wild&nbsp;populations
    .-> FINAL_VARIANT[["`**Final Filtered Variants**
            *N* = 229
            Ref. *Pyrus pyrifolia*
            and Ref. *Pyrus communis*`"]]


REF_GENOME_BOTH
    --> eggNOG{{eggNOG-mapper, PlantTFDB, NLGenomeSweeper, TAIR, UniPort, FLOR-ID, GO, KEGG}}
    --> anno(["`**Genes Functional Annotation**`"])


FINAL_VARIANT
    -- Wild and cultivated populations
    --> sift4g{{SIFT4G}}
    --> sift(["`**Deleterious Variants Prediction**`"])


FINAL_VARIANT
    --> smc{{SMC++}}
    --> history(["` **Historical _Ne_** `"])


FINAL_VARIANT
    -- Cultivated and wild populations
       (focus on cultivars)
    --> omegaplus{{OmegaPlus}}
    --> omega(["`**_ω_**`"])

regionR{{regionR}}
    -. Permutation Test for Selection Regions
    .-> omega


regionR
    -. Permutation Test for Selection Regions
    .-> tajima

pixy2{{Pixy}} 
    -->tajima(["`**Tajima's *D***`"])

FINAL_VARIANT
    --> Dsutie{{Dsuite}}
    --> abba_baba(["`**D and f4-ratio statistic**`"])

FINAL_VARIANT
    --> FINAL_ALL_SITES[["`**Final All Sites**
            *N* = 229
            both invariant and variant sites`"]]


FILTERED_INVARIANT
    --> FINAL_ALL_SITES


FINAL_ALL_SITES
    --> pixy{{Pixy}}
    --> Pi(["`**_π_, _F_<sub>ST</sub>, _d<sub>XY</sub>_**`"])


UNLINKED_VARIANT
    --> FINAL_UNLINKED_VARIANT

FINAL_UNLINKED_VARIANT
    --> stacks{{Stacks Population}}
    --> statistics(["`**Population Genetics Statistics**`"])


FINAL_UNLINKED_VARIANT
    --  Merge&nbsp;18&nbsp;loquat&nbsp;samples&nbsp;as&nbsp;the&nbsp;outgroup
    --> svd{{SVDQuartets}}
    --> sptree(["`**Species Tree**`"])


svd
    -.-> modeling


FINAL_UNLINKED_VARIANT
    --> fastsimcoal{{fastsimcoal}}
    --> modeling(["` **Demographic Modeling** `"])


ABBA{{Dsuite}}
    -.-> modeling


FILTERED_VARIANT
    --> ld{{PopLDdecay}}
    --> LD(["`**LD decay**`"])


%% Subgraph section --------------------------------
subgraph Fig.1 Population Structure
    pca
    tree
    structure
end

subgraph Fig.2 Demographic History
    history
    Pi
    abba_baba
    sptree
    modeling
    statistics
end


subgraph Fig.4 Genetic Burden
    sift
end


subgraph "Fig.3&nbsp;Selection&nbsp;Footprint&nbsp;in&nbsp;Domestication"
    omega
    tajima
    anno
end



%% Legend section ------------------------------------

DATA[["`**Data**`"]]
    ~~~
software{{Software}}

filter[\Filters/]
    ~~~
analysis([Analysis])


subgraph Legend
    DATA
    software
    filter
    analysis
end


%% Configurations -----------------------------------
classDef empty fill:transparent,width:-1px,height:-1px,stroke:transparent;
classDef data fill:transparent,width:-1px,height:-1px,stroke:transparent;
classDef software fill:transparent,width:-1px,height:-1px,stroke:transparent;
classDef filter fill:transparent,width:-1px,height:-1px,stroke:transparent;
```


## The workflow for bioinformatic analysis used in this study (part integrating SNPs and TIPs)

```mermaid
%%{init: {'theme':'forest',
  "themeVariables": {
    "fontSize": "14px",
    "fontFamily": "Arial"}
}}%%

flowchart TD

%% Legend ----------------

DATA[["`**Data**`"]]
    ~~~
software{{Software}}

workflow[(Workflow)]
    ~~~
analysis([Analysis])


subgraph Legend
    DATA
    software
    workflow
    analysis
end


%% Workflow ----------------

ROW_DATA[["`**Raw Reads Data**
    (fastq)
    This study newly sequenced (major) and published data
    `"]]
    --> QUALITY_CONTROL[("`Quality Control Workflow`")]

    --> CLEAN_DATA[["`**Clean Reads Data**
            (fastq)`"]]
    --> bwa{{BWA2-MEM}}
    
    --> ALIGNMENT
    -- \- Exclude samples with depth < 18X or > 50X

    --> megane{{MEGAnE call_genotype}}


REFERENCE[["`**Each Reference Genomes**
    (fasta)
    *Pyrus pyrifolia*
    or *Pyrus communis*`"]]

  --> kemer{{MEGAnE build_kmerset}}
  --> megane

REFERENCE
    --> bastdb{{MEGAnE build_blastdb}}
    --> megane

REFERENCE
    --> megane

ALIGNMENT[["`**Aligned Reads**
    (cram)
    Ref. _Pyrus pyrifolia_
    and Ref. _Pyrus communis_`"]]

ALIGNMENT
    --> snp[(SNP Calling and Quantification Workflow)]
    --> population(["`Population Structure Analysis
     **SNPs**`"])
    --> population_con(["`Population Structure Comparison Analysis
    **SNPs + TIPs**`"])

snp --> selection

ASSEMBLY_W[["`**Occidental Reference Genome**
    (fasta)
    _Pyrus communis_`"]]
    -- TE library building
    --> repet{{TEdenovo REPET}}
    -- Built TE library
    --> teannt{{The second TEannot process}}
    -- Curated TE library
    --> pastec{{PASTEC}}
    -- Classified TE library
    --> teannt2{{TEannot REPET}}
    -- "TEs&nbsp;annotation&nbsp;of&nbsp;the&nbsp;reference&nbsp;genome
        and format conversion"
    -->TE_LIB_W[["`**TE Library for Occidental Genome**
            (fasta, gff3)`"]]
    --> megane

selection(["`Selection Footprint Analysis
    **SNPs**`"])
    --> selection_con(["`Selection Footprint Comparison Analysis
    **SNPs + TIPs**`"])

ASSEMBLY_E[["`**Oriental Reference Genome**
    (fasta)
    _Pyrus pyrifolia_`"]]
    -- TE library building
    --> repet0{{TEdenovo REPET}}
    -- Built TE library
    --> teannt0{{The second TEannot process}}
    -- Curated TE library
    --> pastec0{{PASTEC}}
    -- Classified TE library
    --> teannt1{{TEannot REPET}}
    -- "TEs&nbsp;annotation&nbsp;of&nbsp;the&nbsp;reference&nbsp;genome
        and format conversion"
    -->TE_LIB_E[["`**TE Library for Oriental Genome**
            (fasta, gff3)`"]]
    --> megane


megane
    -- Call and genotype polymorphic MEs
    --> TE_GENOTYPE[["`**TE Insertion / Absence Genotype**
            (vcf)
            MEI and MEA variants`"]]

    --> jointCall{{MEGAnE joint_calling_hs, reshape_vcf}}
    -- \- Exclude TEs with low confidence 
    --> jointTE[["`**TE Variants**
            (vcf)
            Ref. _Pyrus pyrifolia_
            and Ref. _Pyrus communis_`"]]
    --> population_con

jointTE
    --> selection_con




classDef empty fill:transparent,width:-1px,height:-1px,stroke:transparent;
classDef data fill:transparent,width:-1px,height:-1px,stroke:transparent;
classDef software fill:transparent,width:-1px,height:-1px,stroke:transparent;
classDef filter fill:transparent,width:-1px,height:-1px,stroke:transparent;
```  

