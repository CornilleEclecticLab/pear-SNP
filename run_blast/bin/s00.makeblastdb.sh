makeblastdb -in osa1_r7.gene_models.repr.pep.fa -dbtype prot -out Osa
makeblastdb -in TAIR10_pep_20110103_representative_gene_model -dbtype prot -out Ath

sed 's/^>tr|/>/g' UniProtKB_pear_Pyrus_taxid_3766.fa | sed 's/^>sp|/>/g' | sed 's/|/ | /g' > UniProtKB_pear_Pyrus_taxid_3766_format.fa
makeblastdb -in UniProtKB_pear_Pyrus_taxid_3766_format.fa -dbtype prot -out Pyr

