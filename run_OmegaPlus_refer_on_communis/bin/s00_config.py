batch='pear_Jul2024'
input_merged_vcf_path = '/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter_refer_on_communis/output/pear_Jul2024_ref_comm.branch17.whole_pear/pear_Jul2024_ref_comm.Combine_chr.vcf.gz'

variant_vcf_list_filename = 's01.variant_vcf_list.txt'
chromosomes_list_filename = 's01.chr_list.txt'
individual_pop_map_filename = 's01.sample_tab_population.txt'

omega_plus_output_dir = 's03.generate_run_OmegaPlus'

cutoff_files_dir = 's03.find_FPR_cutoff_OmegaPlus'
cutoff_files_suffix = 'OmegaPlus.FPR0005.cutoff.txt'

load_bedtools = "module load bedtools"

gff_filename = "PyrusCommunis_BartlettDHv2.0.Chr.gene.slop1kb.gff"

cds_fasta_filename = "PyrusCommunis_BartlettDHv2.0.cds.fasta"

mask_pass_filename = "PyrusCommunis_BartlettDHv2.0.chr_list.fasta_centromere_range.txt.genmap_pass_subtract_centier.bed"

masked_omega_plus_output_dir = 's06.mask_results'
