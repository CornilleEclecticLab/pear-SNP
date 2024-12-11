# Load the R module
module load r/4.3.1

# Run the R script
gff_file='../input/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.gff'
fasta_file='../input/assembly_Pyrus_communis_Bartlett_DH_Genome_v2.0/PyrusCommunis_BartlettDHv2.0.fasta'
sample_file='../input/s01.3pops.3sample.outgroup_ussu.txt'
outgroup_pop='ussu'

vcf_prefix="../output/s01.get_vcf_from_variant_for_pops/pear_Jul2024."
vcf_suffix=".subpops.vcf"
