#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/en/generated"
w="$PWD/ce-output/pristine"
old="$w/cr_int_t1d_en_prior"
new="$w/cr_int_t1d_en_new"
mkdir -p "$out/new" "$out/merged"
(cd custody/em/generated && sha256sum -c sha256.txt)
python custody/en/prepare.py | tee "$out/prepare.log"
python custody/cm/cm_audit_corrector.py "$new/odd_corrector.npz" "$w/cr_int_t1d_bo/odd_joint.npz" "$old/odd_merged_corrector.npz" "$out/receiver-audit.json" | tee "$out/receiver-audit.log"
cp "$new/"{odd_corrector.npz,odd_U_tail.npz,odd_new.npy,by_action_check.py} "$out/new/"
for pair in 1:2047 2049:6143 6145:11263 11265:16383 16385:24575 24577:32767 32769:40959 40961:49151 49153:57343 57345:65535; do
 lo="${pair%:*}"; hi="${pair#*:}"
 python "$new/bq_action.py" 1 "$hi" "$lo" | tee "$out/new/action-$lo-$hi.log"
 cp "$new/odd_U_${lo}_${hi}.npz" "$out/new/"
done
python "$new/combine.py" | tee "$out/new/combine.log"
cp "$new/odd_U_full.npz" "$out/new/"
python "$new/by_action_check.py" | tee "$out/new/independent-action-check.log"
python custody/en/merge.py "$old" "$new" "$new" | tee "$out/merge.log"
cp "$new/"{odd_merged_corrector.npz,odd_combined_projection.npz,odd_merged_tail.npz} "$out/merged/"
python "$new/by_energy.py" 1 | tee "$out/energy.log"
cp "$new/odd_BQ_energy.npz" "$out/merged/"
python custody/en/check.py | tee "$out/check.log"
