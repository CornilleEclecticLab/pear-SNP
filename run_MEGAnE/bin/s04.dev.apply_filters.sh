bcftools view \
    --apply-filters .,PASS \
    final.vcf.gz \
    --genotype ^miss \
    --include 'TYPE="snp"' \
    --output-type z \
    --output-file final.filtered.vcf.gz
