# Load R module (adjust as per your cluster)
module load r/4.3.1

# Run the R script
gff_file='../input/GWHBAOS00000000.no_blank.gff'
fasta_file='../input/GWHBAOS00000000.genome.fasta'
sample_file='../input/s01.6pops.3sample.outgroup_cauc.txt'
outgroup_pop='cauc'
chr_list='../input/s01.chr_list.txt'

vcf_prfix="../output/s01.get_vcf_from_variant_for_pops/*."
vcf_suffix=".subpops.vcf"