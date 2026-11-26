#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/dw/generated"
w="$PWD/ce-output/pristine"
frozen="$PWD/custody/dv/generated/trial31"
python custody/dw/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/dw/dw_compare.py "$w/"
python "$w/dw_compare.py" odd "$out" "$frozen" | tee "$out/compare.log"
