#!/bin/bash

#SBATCH --job-name=s07.enrich
#SBATCH --output=s07_log/s07.enrich.%J.out
#SBATCH --error=s07_log/s07.enrich.%J.err
#SBATCH -c 4
#SBATCH --time=00:30:00



# @File     :   s07.run_enrichment.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi

#Description:
	# For multicore, set N_CORE in s07.x1.enrichment.R and raise -c to match.

module load r/4.3.1

mkdir -p s07_log

Rscript s07.x1.enrichment.R
