batch='pear_Jul2024'
input_merged_vcf_path = '/shared/ifbstor1/projects/pear_snp3/pear/run_Population_genetic_filter_refer_on_communis/output/pear_Jul2024_ref_comm.branch17.whole_pear/pear_Jul2024_ref_comm.Combine_chr.vcf.gz'

variant_vcf_list_filename = 's01.variant_vcf_list.txt'
chromosomes_list_filename = 's01.chr_list.txt'
karyotype_filename = 'karyotype.txt'
individual_pop_map_filename = 's01.sample_tab_population.txt'

sweed_output_dir = 'run_SweeD_refer_on_communis/output/s03.generate_run_SweeD'
omegaplus_output_dir = 'run_OmegaPlus_refer_on_communis/output/s03.generate_run_OmegaPlus'
raisd_output_dir = 's02.generate_run_RAiSD'
raisd_files_name = 'RAiSD_Report.{pop}.{chr}.by_grid1000.{chr}'
raisd_files_by_win_name = 'RAiSD_Report.{pop}.{chr}.by_win50.{chr}'

cutoff_files_dir = 's03.find_FPR_cutoff_RAiSD'
cutoff_files_suffix = 'RAiSD.FPR005.cutoff.txt'

load_bedtools = "module load bedtools"

grid_window_size = 1000

gff_filename = "PyrusCommunis_BartlettDHv2.0.Chr.gene.gff"

cds_fasta_filename = "PyrusCommunis_BartlettDHv2.0.cds.fasta"
