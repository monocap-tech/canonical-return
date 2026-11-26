#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/eg/generated"
w="$PWD/ce-output/pristine"
p="$w/cr_int_t1d_eg_physical"
h="$w/cr_int_t1d_eg_head"
python custody/eg/trial_prepare.py | tee "$out/trial-prepare.log"
python custody/eg/prepare.py | tee "$out/prepare.log"
for pair in 16385:18431 65025:65535; do
 lo="${pair%:*}"; hi="${pair#*:}"
 cp "custody/ef/generated/odd_U_${lo}_${hi}.npz" "$out/"
done
for pair in 18433:26623 26625:34815 34817:43007 43009:51199 51201:59391 59393:65023; do
 lo="${pair%:*}"; hi="${pair#*:}"
 python "$p/bq_action.py" 1 "$hi" "$lo" | tee "$out/action-$lo-$hi.log"
 cp "$p/odd_U_${lo}_${hi}.npz" "$out/"
done
python "$h/bp_head_extend.py" 1 65536 | tee "$out/head.log"
cp "$h/odd_head_full.npz" "$out/"
python custody/eg/check.py | tee "$out/check.log"
