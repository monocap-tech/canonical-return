#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ei/generated"
w="$PWD/ce-output/pristine"
mkdir -p "$out"
cp custody/ei/compare.py "$w/ei_compare.py"
python "$w/ei_compare.py" odd "$out" | tee "$out/compare.log"
