#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/dm/generated"
w="$PWD/ce-output/pristine"
frozen="$PWD/custody/dk/generated/trial26"
python custody/dm/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/dm/dm_compare.py "$w/"
python "$w/dm_compare.py" odd "$out" "$frozen" | tee "$out/compare.log"
