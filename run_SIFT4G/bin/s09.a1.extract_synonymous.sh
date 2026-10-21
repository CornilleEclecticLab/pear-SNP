#!/usr/bin/env bash

#SBATCH -J s09.a1.extract_synonymous.sh
#SBATCH -o s09.log/s09.a1.extract_synonymous.sh.%J.out
#SBATCH -e s09.log/s09.a1.extract_synonymous.sh.%J.err
#SBATCH -c 18
#SBATCH --mem=16G

#Description:
    # Extract synonymous sites from the s08 snpEff + VEP annotated VCFs, as the
    # denominator of LoF burden (independent of SIFT4G).
    # A site is synonymous when both VEP and snpEff call synonymous_variant on a
    # protein_coding transcript.

    # Usage: sbatch s09.a1.extract_synonymous.sh PCOM
    #        sbatch s09.a1.extract_synonymous.sh PPY

set -euo pipefail

ORG_VERSION="${1:?Usage: sbatch s09.a1.extract_synonymous.sh <PCOM|PPY>}"

source s00.config_${ORG_VERSION}.sh

module purge
module load bcftools/1.14

LOF_OUTPUT_DIR="$PROJECT_ROOT/output/s08.calculate_LoF/$ORG_VERSION"


while read -r CHR; do
  (
  PFX="$LOF_OUTPUT_DIR/$CHR/$ORG_VERSION.$CHR"
  VCF="$PFX.snpEff.vep.vcf.gz"

  # VEP synonymous sites
  bcftools +split-vep "$VCF" -a CSQ -d \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%BIOTYPE\t%Consequence\n' \
    | awk -F'\t' -v OFS='\t' \
        '$5=="protein_coding" && $6~/(^|&)synonymous_variant(&|$)/ {print $1,$2,$3,$4}' \
    | sort -u \
    > "$PFX.vep_synonymous.tmp"

  # snpEff synonymous sites
  bcftools query "$VCF" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
    | awk -F'\t' -v OFS='\t' '{
        n = split($5, anns, ",")
        for (i = 1; i <= n; i++) {
          split(anns[i], a, "|")
          if (a[8] == "protein_coding" && a[16] !~ /WARNING_TRANSCRIPT/ && a[2] ~ /(^|&)synonymous_variant(&|$)/) {
            print $1, $2, $3, $4; break
          }
        }
      }' \
    | sort -u \
    > "$PFX.snpEff_synonymous.tmp"

  # Both tools agree
  comm -12 "$PFX.vep_synonymous.tmp" "$PFX.snpEff_synonymous.tmp" \
    | sort -t$'\t' -k2,2n \
    > "$PFX.synonymous_consensus.tsv"

  rm "$PFX.vep_synonymous.tmp" "$PFX.snpEff_synonymous.tmp"
  ) &
done < "$CHR_LIST"

wait

cat "$LOF_OUTPUT_DIR/"*/"$ORG_VERSION".*.synonymous_consensus.tsv \
  > "$LOF_OUTPUT_DIR/merged.synonymous_consensus.tsv"

wc -l "$LOF_OUTPUT_DIR/merged.synonymous_consensus.tsv"
