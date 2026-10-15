#!/bin/bash

# v1.0.1
git clone https://github.com/grenaud/ROHan.git

ROHAN_BIN_DIR=$(dirname "$(realpath ./ROHan/bin/rohan)")

echo "export PATH=\"$ROHAN_BIN_DIR:\$PATH\"" > s00.load_rohan.sh
