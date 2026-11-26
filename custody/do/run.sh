#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/do/generated"
w="$PWD/ce-output/pristine"
frozen="$PWD/custody/dk/generated/trial26"
python custody/do/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/do/do_compare.py "$w/"
python "$w/do_compare.py" odd "$out" "$frozen" | tee "$out/compare.log"
