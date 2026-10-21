#!/usr/bin/env bash

#SBATCH -J s08.a1.build_snpEff_db.sh
#SBATCH -o s08.log/s08.a1.build_snpEff_db.sh.%J.out
#SBATCH -e s08.log/s08.a1.build_snpEff_db.sh.%J.err
#SBATCH -c 2
#SBATCH --mem=16G

#Description:
    # Build a snpEff database from the same GFF as VEP

set -euo pipefail

ORG_VERSION="${1:?Usage: sbatch s08.a1.build_snpEff_db.sh <PCOM|PPY>}"

source s00.config_${ORG_VERSION}.sh

module purge
module load java-jdk/11.0.9.1
module load samtools/1.14

SNPEFF_JAR="${PROJECT_ROOT}/../run_Population_genetic_filter/bin/snpEff/snpEff.jar"
SNPEFF_DB_DIR="${PROJECT_ROOT}/database/snpEff"
SNPEFF_CONFIG="${SNPEFF_DB_DIR}/snpEff.config"
GENOME_DIR="${SNPEFF_DB_DIR}/data/${ORG_VERSION}"
mkdir -p "$GENOME_DIR"


# GFF cleaning
grep -v "^#" "${GFF_FILE%.gff}.chr.gff" \
  | awk -F'\t' 'FNR==NR {chr[$1]; next} NF==9 && $3!="intron" && ($1 in chr)' "$CHR_LIST" - \
  | sed 's/;;/;/g' \
  | sort -t$'\t' -k1,1 -k4,4n -k5,5nr \
  > "$GENOME_DIR/genes.raw_phase.gff"

awk '
  BEGIN { FS = OFS = "\t" }
  function parent(attr) { match(attr, /Parent=[^;]+/); return substr(attr, RSTART + 7, RLENGTH - 7) }
  # Pass 1: collect CDS per transcript (input is sorted by start)
  FNR == NR {
    if ($3 == "CDS") { t = parent($9); n[t]++; s[t, n[t]] = $4; e[t, n[t]] = $5; st[t] = $7 }
    next
  }
  # Phase of each CDS = bases needed to reach the next codon start, counted 5 -> 3
  FNR == 1 {
    for (t in n) {
      L = 0
      if (st[t] == "-") { for (i = n[t]; i >= 1; i--) { ph[t, s[t, i]] = (3 - L % 3) % 3; L += e[t, i] - s[t, i] + 1 } }
      else              { for (i = 1; i <= n[t]; i++)  { ph[t, s[t, i]] = (3 - L % 3) % 3; L += e[t, i] - s[t, i] + 1 } }
    }
  }
  # Pass 2: rewrite phase
  $3 == "CDS" { $8 = ph[parent($9), $4] }
  { print }
' "$GENOME_DIR/genes.raw_phase.gff" "$GENOME_DIR/genes.raw_phase.gff" \
  > "$GENOME_DIR/genes.gff"

echo "CDS phases changed: $(paste <(cut -f8 "$GENOME_DIR/genes.raw_phase.gff") <(cut -f8 "$GENOME_DIR/genes.gff") | awk '$1!=$2' | wc -l)"
rm "$GENOME_DIR/genes.raw_phase.gff"


# Sequences
xargs samtools faidx "$GENOME_FASTA_FULL" < "$CHR_LIST" > "$GENOME_DIR/sequences.fa"


# Config
touch "$SNPEFF_CONFIG"
grep -q "^data.dir" "$SNPEFF_CONFIG" || echo "data.dir = ${SNPEFF_DB_DIR}/data/" >> "$SNPEFF_CONFIG"
grep -q "^${ORG_VERSION}.genome" "$SNPEFF_CONFIG" || echo "${ORG_VERSION}.genome : ${ORG}" >> "$SNPEFF_CONFIG"


# Build snpEff database
java -Xmx12g -jar "$SNPEFF_JAR" build \
  -gff3 -v \
  -noCheckCds -noCheckProtein \
  -c "$SNPEFF_CONFIG" \
  "$ORG_VERSION"
