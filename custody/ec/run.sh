#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ec/generated"
w="$PWD/ce-output/pristine"
frozen="$PWD/custody/dz/generated/trial31"
python custody/ec/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/ec/ec_compare.py "$w/"
python "$w/ec_compare.py" odd "$out" "$frozen" | tee "$out/compare.log"
