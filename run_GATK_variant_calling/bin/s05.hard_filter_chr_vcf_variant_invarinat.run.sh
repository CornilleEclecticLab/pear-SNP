cd s05.hard_filter_chr_vcf_variant_invarinat
sbatch -c 4 --mem=4G Pp01.s05.sh
sbatch -c 4 --mem=4G Pp02.s05.sh
sbatch -c 4 --mem=4G Pp03.s05.sh
sbatch -c 4 --mem=4G Pp04.s05.sh
sbatch -c 4 --mem=4G Pp05.s05.sh
sbatch -c 4 --mem=4G Pp06.s05.sh
sbatch -c 4 --mem=4G Pp07.s05.sh
sbatch -c 4 --mem=4G Pp08.s05.sh
sbatch -c 4 --mem=4G scaffolds.list.s05.sh
