#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/eb/generated"
w="$PWD/ce-output/pristine"
frozen="$PWD/custody/dz/generated/trial31"
python custody/eb/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/eb/eb_compare.py "$w/"
python "$w/eb_compare.py" odd "$out" "$frozen" | tee "$out/compare.log"
