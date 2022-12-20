cd s04.combine_chr_gvcf_to_vcf_allsite
sbatch -c 2 --mem=32G Pp01.s04.sh
sbatch -c 2 --mem=32G Pp02.s04.sh
sbatch -c 2 --mem=32G Pp03.s04.sh
sbatch -c 2 --mem=32G Pp04.s04.sh
sbatch -c 2 --mem=32G Pp05.s04.sh
sbatch -c 2 --mem=32G Pp06.s04.sh
sbatch -c 2 --mem=32G Pp07.s04.sh
sbatch -c 2 --mem=32G Pp08.s04.sh
sbatch -c 2 --mem=32G scaffolds.list.s04.sh
