#!/bin/bash

#SBATCH -J s01.build_sift4g_db.sh
#SBATCH -o s01.build_sift4g_db.sh.%J.out
#SBATCH -e s01.build_sift4g_db.sh.%J.err
#SBATCH -p long # this job need more than 24 hours
#SBATCH -c 8
#SBATCH --mem=36G

# @File     :   s01.build_sift4g_db.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/05/28 18:01:45
#Description:
    # 

set -euo pipefail

source s00.config.sh

# Check container and inputs 
if [ ! -f "$CONTAINER" ]; then
  echo "ERROR: File not found $CONTAINER"
  exit 1
fi
for f in "$GENOME_FASTA" "$PROTEIN_DB" "$GFF_FILE"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found $f"
    exit 1
  fi
done

# Create directory structure
echo "==> Creating directory under $OUTPUT_DIR"
mkdir -p \
  "$OUTPUT_DIR/chr-src" \
  "$OUTPUT_DIR/gene-annotation-src" 

# Move inputs
echo "==> Copying FASTA / GFF to corresponding directories"
mv "$GENOME_FASTA"      "$OUTPUT_DIR/chr-src/"
mv "$GTF_FILE"          "$OUTPUT_DIR/gene-annotation-src/"

# Generate configuration file
# echo "==> Generating SIFT4G configuration file: $CONFIG_FILE"
cat > "$CONFIG_FILE" <<EOF

GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=11
MITO_GENETIC_CODE_TABLENAME=Plant Plastid Code


PARENT_DIR=${OUTPUT_DIR}
ORG=${ORG}
ORG_VERSION=${ORG_VERSION}

# Running SIFT4G, this path works for the Dockerfile
SIFT4G_PATH=/sift4g/bin/sift4g

#PROTEIN_DB needs to be uncompressed
PROTEIN_DB=${PROTEIN_DB}


# Sub-directories, don't need to change
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=Log.txt
ZLOGFILE=Log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP

# Don't need to change
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord
EOF

echo "--- Configuration content below ---"
cat "$CONFIG_FILE"
echo "-------------------"

# Call container to build database
echo "==> Starting SIFT4G database build"
apptainer run \
  "$CONTAINER" \
    make-SIFT-db-all.pl \
      -config "$CONFIG_FILE"

echo "==> Completed! Database directory: $OUTPUT_DIR/${ORG_VERSION}"
