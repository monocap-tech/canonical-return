#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ee/generated"
w="$PWD/ce-output/pristine"
python custody/ee/trial_prepare.py | tee "$out/trial-prepare.log"
cp custody/ee/census.py "$w/"
python "$w/census.py" odd "$out" | tee "$out/census.log"
