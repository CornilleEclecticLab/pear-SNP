#! /bin/bash
set -euo pipefail

. drive_DemoInfer.sh

# prepare bootstrap blocks
# myfunc_get_boot_blocks

## step5
export -f myfunc_bssbatach
bs=($(seq -w 1 50))

module load parallel
parallel -j 25 myfunc_bssbatach ::: "${bs[@]}"
