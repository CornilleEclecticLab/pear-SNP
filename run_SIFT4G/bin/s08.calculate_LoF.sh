#!/usr/bin/env bash

#SBATCH -J s08.calculate_LoF.sh
#SBATCH -o s08.log/s08.calculate_LoF.sh.%J.out
#SBATCH -e s08.log/s08.calculate_LoF.sh.%J.err
#SBATCH -c 18
#SBATCH --mem=64G

#Description:
    # Loss-of-function (LoF) annotation with snpEff + VEP
    
    # Requires: s08.a1.build_snpEff_db.sh
    # Usage: sbatch s08.calculate_LoF.sh PCOM
    #        sbatch s08.calculate_LoF.sh PPY


ORG_VERSION="${1:?Usage: sbatch s08.calculate_LoF.sh <PCOM|PPY>}"

source s00.config_${ORG_VERSION}.sh

module purge
module load java-jdk/11.0.9.1
module load bcftools/1.14
module load samtools/1.14
module load conda
. activate vep


############ Config ############
SNPEFF_JAR="${PROJECT_ROOT}/../run_Population_genetic_filter/bin/snpEff/snpEff.jar"
SNPEFF_CONFIG="${PROJECT_ROOT}/database/snpEff/snpEff.config"
SNPEFF_DB="$ORG_VERSION"

SNPEFF_MEM="3g"   # per chromosome, chromosomes run in parallel

# Sequence Ontology terms counted as LoF
LOF_TERMS="stop_gained|start_lost|stop_lost|splice_acceptor_variant|splice_donor_variant|frameshift_variant|exon_loss_variant|transcript_ablation"

SIFT_DIR="$PROJECT_ROOT/output/s02_annot/$ORG_VERSION"
LOF_OUTPUT_DIR="$PROJECT_ROOT/output/s08.calculate_LoF/$ORG_VERSION"
mkdir -p "$LOF_OUTPUT_DIR"


# Prepare GFF for VEP
VEP_GFF="$LOF_OUTPUT_DIR/$ORG_VERSION.vep.gff.gz"

grep -v "^#" "${GFF_FILE%.gff}.chr.gff" \
  | awk -F'\t' 'NF==9 && $3!="intron"' \
  | sed 's/;;/;/g' \
  | sort -t$'\t' -k1,1 -k4,4n -k5,5nr \
  | bgzip -c \
  > "$VEP_GFF"
tabix -f -p gff "$VEP_GFF"

# Index fasta
[ -s "${GENOME_FASTA_FULL}.fai" ] || samtools faidx "$GENOME_FASTA_FULL"


while read -r CHR; do
  (
  CHR_DIR="$LOF_OUTPUT_DIR/$CHR"
  mkdir -p "$CHR_DIR"
  PFX="$CHR_DIR/$ORG_VERSION.$CHR"

  # VCF
  bcftools view -G \
    "$INPUT_DIR/${ORG_VERSION}/"*."$CHR".*.vcf.gz \
    > "$PFX.sites.vcf"

  # snpEff
  java -Xmx${SNPEFF_MEM} -jar "$SNPEFF_JAR" \
    -c "$SNPEFF_CONFIG" \
    -csvStats "$PFX.snpEff_stats.csv" \
    -s "$PFX.snpEff_summary.html" \
    "$SNPEFF_DB" \
    "$PFX.sites.vcf" \
    > "$PFX.snpEff.vcf"

  # VEP
  vep \
    -i "$PFX.snpEff.vcf" \
    --gff "$VEP_GFF" \
    --fasta "$GENOME_FASTA_FULL" \
    --vcf \
    --symbol \
    --biotype \
    --no_stats \
    --force_overwrite \
    -o STDOUT \
    | bgzip -c > "$PFX.snpEff.vep.vcf.gz"
  tabix -f -p vcf "$PFX.snpEff.vep.vcf.gz"

  rm "$PFX.sites.vcf" "$PFX.snpEff.vcf"

  # VEP LoF: CHROM POS REF ALT GENE TRANSCRIPT CONSEQUENCE IMPACT PROTEIN_POS AA
  bcftools +split-vep "$PFX.snpEff.vep.vcf.gz" -a CSQ -d \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%Gene\t%Feature\t%BIOTYPE\t%Consequence\t%IMPACT\t%Protein_position\t%Amino_acids\n' \
    | awk -F'\t' -v OFS='\t' -v terms="$LOF_TERMS" \
        '$7=="protein_coding" && $8~("(^|&)(" terms ")(&|$)") {print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11}' \
    > "$PFX.vep_LoF.tsv"

  # snpEff LoF: CHROM POS REF ALT GENE TRANSCRIPT EFFECT IMPACT HGVS_P LOF_TAG
  bcftools query "$PFX.snpEff.vep.vcf.gz" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/LOF\n' \
    | awk -F'\t' -v OFS='\t' -v terms="$LOF_TERMS" '{
        n = split($5, anns, ",")
        for (i = 1; i <= n; i++) {
          split(anns[i], a, "|")
          # a[2] effect, a[3] impact, a[5] gene, a[7] transcript, a[8] biotype, a[11] HGVS.p, a[16] warnings
          # skip transcripts snpEff flags as broken (e.g. MULTIPLE_STOP_CODONS)
          if (a[8] == "protein_coding" && a[16] !~ /WARNING_TRANSCRIPT/ && a[2] ~ ("(^|&)(" terms ")(&|$)"))
            print $1, $2, $3, $4, a[5], a[7], a[2], a[3], a[11], $6
        }
      }' \
    > "$PFX.snpEff_LoF.tsv"

  # merge with SIFT4G LoF, one line per site
  awk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { k = $1 SUBSEP $2 SUBSEP $4; site[k] = $1 OFS $2 OFS $3 OFS $4
                          gene[k] = $5; vep[k] = (k in vep) ? vep[k] "," $7 : $7; next }
    FILENAME == ARGV[2] { k = $1 SUBSEP $2 SUBSEP $4; site[k] = $1 OFS $2 OFS $3 OFS $4
                          snp[k] = (k in snp) ? snp[k] "," $7 : $7; tag[k] = $10; next }
    FILENAME == ARGV[3] { k = $1 SUBSEP $2 SUBSEP $4; site[k] = $1 OFS $2 OFS $3 OFS $4
                          sift[k] = $9; if (!(k in gene)) gene[k] = $5; next }
    END {
      for (k in site) {
        v = (k in vep); s = (k in snp); f = (k in sift)
        print site[k], (k in gene) ? gene[k] : ".", v ? vep[k] : ".", s ? snp[k] : ".",
              (k in tag) ? tag[k] : ".", f ? sift[k] : ".", v, s, f, (v && s)
      }
    }' \
    "$PFX.vep_LoF.tsv" \
    "$PFX.snpEff_LoF.tsv" \
    "$SIFT_DIR/$CHR/$ORG_VERSION.$CHR.loss_of_funtion.tsv" \
    | sort -t$'\t' -k2,2n \
    > "$PFX.LoF_merged.tsv"

  # Extract SNP positions
  cut -f 1,2 "$PFX.vep_LoF.tsv"    | sort -u -k2,2n > "$PFX.vep_LoF.pos.txt"
  cut -f 1,2 "$PFX.snpEff_LoF.tsv" | sort -u -k2,2n > "$PFX.snpEff_LoF.pos.txt"
  awk '$13==1' "$PFX.LoF_merged.tsv"            | cut -f 1,2 > "$PFX.LoF_consensus.pos.txt"
  awk '$13==1 || $12==1' "$PFX.LoF_merged.tsv"  | cut -f 1,2 > "$PFX.LoF_combined.pos.txt"
  ) &
done < "$CHR_LIST"

wait


# Concatenate results across chromosomes
for SUFFIX in vep_LoF.tsv snpEff_LoF.tsv LoF_merged.tsv \
              vep_LoF.pos.txt snpEff_LoF.pos.txt LoF_consensus.pos.txt LoF_combined.pos.txt; do
  cat "$LOF_OUTPUT_DIR/"*/"$ORG_VERSION".*."$SUFFIX" > "$LOF_OUTPUT_DIR/merged.$SUFFIX"
done

sed -i '1i CHROM\tPOS\tREF\tALT\tGENE\tVEP_CONSEQUENCE\tSNPEFF_EFFECT\tSNPEFF_LOF_TAG\tSIFT_EFFECT\tIS_VEP\tIS_SNPEFF\tIS_SIFT\tIS_CONSENSUS' \
  "$LOF_OUTPUT_DIR/merged.LoF_merged.tsv"


# Summary
awk -F'\t' -v OFS='\t' -v org="$ORG_VERSION" 'NR > 1 {
    vep += $10; snp += $11; sift += $12; cons += $13
    both += ($12 && $13); sift_only += ($12 && !$13); cons_only += ($13 && !$12)
    comb += ($12 || $13)
  }
  END {
    print "REF", "VEP", "snpEff", "SIFT4G", "consensus(VEP&snpEff)",
          "SIFT4G&consensus", "SIFT4G_only", "consensus_only", "combined(SIFT4G|consensus)"
    print org, vep, snp, sift, cons, both, sift_only, cons_only, comb
  }' "$LOF_OUTPUT_DIR/merged.LoF_merged.tsv" \
  > "$LOF_OUTPUT_DIR/summary.LoF_counts.tsv"


awk -F'\t' 'NR > 1 && $13 == 1 {split($6, c, ","); split(c[1], t, "&"); print t[1]}' \
  "$LOF_OUTPUT_DIR/merged.LoF_merged.tsv" \
  | sort | uniq -c | sort -k1,1nr \
  > "$LOF_OUTPUT_DIR/summary.LoF_consensus_consequence.txt"

cat "$LOF_OUTPUT_DIR/summary.LoF_counts.tsv"
cat "$LOF_OUTPUT_DIR/summary.LoF_consensus_consequence.txt"
