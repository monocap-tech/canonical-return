#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/em/generated"
cb="$PWD/ce-output/pristine/cr_int_t1d_cb"
bo="$PWD/ce-output/pristine/cr_int_t1d_bo"
mkdir -p "$out/preparation"
python custody/em/trial_prepare.py | tee "$out/preparation.log"
python custody/em/targeted_selection.py | tee "$out/selection.log"
python "$cb/bq_corrector.py" 1 | tee "$out/corrector.log"
cp "$cb/odd_corrector.npz" "$out/"
python "$cb/bq_tail.py" 1 | tee "$out/tail.log"
cp "$cb/odd_U_tail.npz" "$out/"
python custody/cm/cm_audit_corrector.py "$cb/odd_corrector.npz" "$bo/odd_joint.npz" "$cb/odd_merged_corrector.npz" "$out/receiver-audit.json" | tee "$out/receiver-audit.log"
