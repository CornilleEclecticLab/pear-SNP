#!/bin/bash

#SBATCH -J s08.b1.enrichment_pi_matched.sh
#SBATCH -o s08.b1.enrichment_pi_matched.sh.%J.out
#SBATCH -e s08.b1.enrichment_pi_matched.sh.%J.err

# @File     :   s08.b1.enrichment_pi_matched.sh
# @Version  :   1.0.0
# @Author   :   NIE Yuqi

set -euo pipefail
module purge
module load r/4.3.1

Rscript s08.x1.enrichment_pi_matched.R \
  --lineage Oriental \
  --dinv "../output/Dinvestigate/*/*_localFstats__50_25.txt" \
  --export_windows_only

echo "Then go to ../../run_Pixy_PiFstDxy/bin/s02.run_pixy_by_chromosome_for_local_D.py to calculate Pi."
