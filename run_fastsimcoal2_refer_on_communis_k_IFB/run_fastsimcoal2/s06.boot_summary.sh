#! /bin/bash
set -euo pipefail

. drive_DemoInfer.sh

set="west"

mkdir -p ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot_init/summary
perl ./scripts/fsc-selectbestrun.boot.pl \
    "./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot_init/*/*.run*/*/*.bestlhoods" \
    ./analysis/DemoInfer/demo_files/${set}_Bestrun_block_boot_init/summary

