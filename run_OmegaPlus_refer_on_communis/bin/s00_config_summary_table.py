eggNOG = '../../run_eggNOG-mapper/output/PyrusCommunis_BartlettDHv2.0.emapper.annotations'

# gff_ori = '../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.gff'
gff_filename = "../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop1kb.gff"
# gff_updown2k_filename = "../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop2kb_up_down.gff"
# gff_updown10k_filename = "../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop10kb_up_down.gff"
gff_upstream2k_filename = "../input/PyrusCommunis_BartlettDHv2.0.Chr.gene.slop2kb_upstream.gff"

positive_selection_output_list = '../input/s08.summary_list.txt'

chrID_map = '../input/PyrusCommunis_BartlettDHv2.0.chrID_map.txt'

collinearity = '../../run_MCScanX/output/s05.parse_MCScanX/intersection.collinearity.txt'

blast_annotation = '../../run_blast/output/s04.annotate_blast_result/PyrusCommunis_BartlettDHv2.0.pep.formatted.blastp.annotated.tsv'

interest_common_populations_combines = [['comm_Dessert', 'comm_Perry']]

go_terms_list = '../../run_eggNOG-mapper/output/s02.prepare_GO_KEGG/PyrusCommunis_BartlettDHv2.0.GO_GO_term.table.txt'
kegg_terms_list = '../../run_eggNOG-mapper/output/s02.prepare_GO_KEGG/PyrusCommunis_BartlettDHv2.0.KEGG_TERM2GENE.table.txt'
nl_genome_sweeper_annotation = '../../run_NLGenomeSweeper/output/PyrusCommunis_BartlettDHv2.0/NLGenomeSweeper/Final_candidates.bed.NLs.genes.txt'
