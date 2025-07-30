#!/usr/bin/env bash

#SBATCH -J s02.annotate_vcf_against_sift4g_db_PCOM.sh
#SBATCH -o s02.log/s02.annotate_vcf_against_sift4g_db_PCOM.sh.%J.out
#SBATCH -e s02.log/s02.annotate_vcf_against_sift4g_db_PCOM.sh.%J.err
#SBATCH -c 17
#SBATCH --mem=16G

# @File     :   s02.annotate_vcf_against_sift4g_db_PCOM.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi
# @Email    :   nieyuqi.cn@gmail.com
# @Time(CET):   2025/05/28 19:31:08
#Description:
    # 

source s00.config_PCOM.sh
source s00.load_bcftools.sh

ANNOT_OUPUT_DIR="$PROJECT_ROOT/output/s02_annot/$ORG_VERSION/"
mkdir -p $ANNOT_OUPUT_DIR


while read -r CHR; do
  (
  TEMP_DIR="$PROJECT_ROOT/temp/$ORG_VERSION/$CHR"
  TEMP_VCF="$TEMP_DIR/${ORG_VERSION}.${CHR}.vcf"
  
  mkdir -p "$TEMP_DIR"

  SIFT_DIR="${ANNOT_OUPUT_DIR}/${CHR}"
  
  # To save storage, only keep the first 8 columns of the VCF file
  # Use bcftools view to fresh the vcf AN, AC, AF fieelds
  bcftools view \
    "$INPUT_DIR/${ORG_VERSION}/"*."$CHR".*.vcf.gz \
    | grep -v "^##" \
    | cut -f 1-8 \
    > "$TEMP_VCF"

  apptainer run \
    "$CONTAINER" \
    java -jar ./SIFT4G_Annotator.jar \
    -c \
    -i "$TEMP_VCF" \
    -d  "$OUTPUT_DIR/$ORG_VERSION" \
    -r "$SIFT_DIR"

  # SIFT4g Option	Description
  # -c	To run on command line
  # -i	Path to your input variants file in VCF format
  # -d	Path to SIFT database directory
  # -r	Path to your output results folder
  # -t	To extract annotations for multiple transcripts (Optional)


  SIFT_XLS_RAW=$(ls "${SIFT_DIR}/"*_SIFTannotations.xls)
  SIFT_XLS="${SIFT_XLS_RAW}.cleaned"

  sed -i 's/\r$//'  $SIFT_XLS_RAW
  sed -i 's/\r$//'  $SIFT_XLS

  # Remove low-confidence results
  grep -v "Low confidence" $SIFT_XLS_RAW > $SIFT_XLS

  # Extract SIFT results
  grep "START\|STOP" $SIFT_XLS > $SIFT_DIR/$ORG_VERSION.$CHR.loss_of_funtion.tsv
  awk '$9=="NONSYNONYMOUS"' $SIFT_XLS | awk '$17=="DELETERIOUS"' > $SIFT_DIR/$ORG_VERSION.$CHR.deleterious.tsv
  awk '$9=="NONSYNONYMOUS"' $SIFT_XLS | awk '$17=="TOLERATED"'    > $SIFT_DIR/$ORG_VERSION.$CHR.tolerated.tsv
  awk '$9=="SYNONYMOUS"' $SIFT_XLS | grep -v "DELETERIOUS" | awk '$13!="NA"' > $SIFT_DIR/$ORG_VERSION.$CHR.synonymous.tsv
  awk '$8=="CDS"' $SIFT_XLS | awk '$13!="NA"' >$SIFT_DIR/$ORG_VERSION.$CHR.all.cds.tsv
  
  # Extract SNP positions
  cut -f 1,2 $SIFT_DIR/$ORG_VERSION.$CHR.loss_of_funtion.tsv > $SIFT_DIR/$ORG_VERSION.$CHR.loss_of_funtion.pos.txt
  cut -f 1,2 $SIFT_DIR/$ORG_VERSION.$CHR.deleterious.tsv > $SIFT_DIR/$ORG_VERSION.$CHR.deleterious.pos.txt
  cut -f 1,2 $SIFT_DIR/$ORG_VERSION.$CHR.tolerated.tsv > $SIFT_DIR/$ORG_VERSION.$CHR.tolerated.pos.txt
  cut -f 1,2 $SIFT_DIR/$ORG_VERSION.$CHR.synonymous.tsv > $SIFT_DIR/$ORG_VERSION.$CHR.synonymous.pos.txt
  cut -f 1,2 $SIFT_DIR/all.cds.tsv > $SIFT_DIR/$ORG_VERSION.$CHR.all.cds.pos.txt
  ) &
done < "$CHR_LIST"

wait

# Clean up temporary files
rm -r "$PROJECT_ROOT/temp/$ORG_VERSION"

# Concatenate annot_results across chromosomes
cat "$ANNOT_OUPUT_DIR/"*/*SIFTannotations.xls \
  | head -n 1 \
  > $ANNOT_OUPUT_DIR/merged.SIFTannotations.xls

for xls_file in "$ANNOT_OUPUT_DIR/"*/*SIFTannotations.xls; do
  if [ "$xls_file" != "$ANNOT_OUPUT_DIR/merged.SIFTannotations.xls" ]; then
    tail -n +2 "$xls_file" >> $ANNOT_OUPUT_DIR/merged.SIFTannotations.xls
  fi
done

cat "$ANNOT_OUPUT_DIR/"*/*.cleaned \
  | head -n 1 \
  > $ANNOT_OUPUT_DIR/merged.SIFTannotations.xls.cleaned

for cleaned_file in "$ANNOT_OUPUT_DIR/"*/*.cleaned; do
  if [ "$cleaned_file" != "$ANNOT_OUPUT_DIR/merged.SIFTannotations.xls.cleaned" ]; then
    tail -n +2 "$cleaned_file" >> $ANNOT_OUPUT_DIR/merged.SIFTannotations.xls.cleaned
  fi
done

cat "$ANNOT_OUPUT_DIR/"*/*loss_of_funtion.tsv \
> $ANNOT_OUPUT_DIR/merged.loss_of_funtion.tsv

cat "$ANNOT_OUPUT_DIR/"*/*deleterious.tsv \
> $ANNOT_OUPUT_DIR/merged.deleterious.tsv

cat "$ANNOT_OUPUT_DIR/"*/*tolerated.tsv \
> $ANNOT_OUPUT_DIR/merged.tolerated.tsv

cat "$ANNOT_OUPUT_DIR/"*/*synonymous.tsv \
> $ANNOT_OUPUT_DIR/merged.synonymous.tsv

cat "$ANNOT_OUPUT_DIR/"*/*all.cds.tsv \
> $ANNOT_OUPUT_DIR/merged.all.cds.tsv

# Merge SNP positions
cat "$ANNOT_OUPUT_DIR/"*/*loss_of_funtion.pos.txt \
> $ANNOT_OUPUT_DIR/merged.loss_of_funtion.pos.txt

cat "$ANNOT_OUPUT_DIR/"*/*deleterious.pos.txt \
> $ANNOT_OUPUT_DIR/merged.deleterious.pos.txt

cat "$ANNOT_OUPUT_DIR/"*/*tolerated.pos.txt \
> $ANNOT_OUPUT_DIR/merged.tolerated.pos.txt

cat "$ANNOT_OUPUT_DIR/"*/*synonymous.pos.txt \
> $ANNOT_OUPUT_DIR/merged.synonymous.pos.txt

cat "$ANNOT_OUPUT_DIR/"*/*all.cds.pos.txt \
> $ANNOT_OUPUT_DIR/merged.all.cds.pos.txt
