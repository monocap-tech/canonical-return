#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ea/generated"
w="$PWD/ce-output/pristine"
frozen="$PWD/custody/dv/generated/trial31"
python custody/ea/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/ea/ea_compare.py "$w/"
python "$w/ea_compare.py" odd "$out" "$frozen" | tee "$out/compare.log"
