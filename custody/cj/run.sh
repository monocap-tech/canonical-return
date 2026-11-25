#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/cj/generated"
w="$PWD/ce-output/pristine"
ca="$w/cr_int_t1d_ca"
frozen="$PWD/custody/ch/generated/expanded-8192-8192/cf-baseline.json"
cp custody/cj/cj_compare.py "$w/"
cp custody/cf/generated/{bu_trial.npz,joint_circle.npz} "$ca/"
cp custody/ch/generated/expanded-8192-8192/{even_joint_proposal.npz,even_joint_Uaction.npz} "$ca/"
sha256sum "$ca/"{bu_trial.npz,joint_circle.npz,even_joint_proposal.npz,even_joint_Uaction.npz} > "$out/frozen-inputs.sha256"
python "$w/cj_compare.py" even "$out" "$frozen" --control | tee "$out/comparison.log"
sha256sum -c "$out/frozen-inputs.sha256"
