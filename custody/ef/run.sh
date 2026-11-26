#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ef/generated"
w="$PWD/ce-output/pristine"
p="$w/cr_int_t1d_ef_pilot"
python custody/ef/trial_prepare.py | tee "$out/trial-prepare.log"
python custody/ef/pilot_prepare.py | tee "$out/pilot-prepare.log"
for pair in 15873:16383 16385:18431 65025:65535; do
 lo="${pair%:*}"; hi="${pair#*:}"
 python "$p/bq_action.py" 1 "$hi" "$lo" | tee "$out/action-$lo-$hi.log"
 cp "$p/odd_U_${lo}_${hi}.npz" "$out/"
done
python custody/ef/check.py | tee "$out/check.log"
