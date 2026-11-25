#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ck/generated"
w="$PWD/ce-output/pristine"
cb="$w/cr_int_t1d_cb"
mkdir -p "$out/control12" "$out/trial13"
cp custody/ck/ck_compare.py "$w/"
sha256sum "$cb/"{bs_trial.npz,bs_circle.npz,odd_joint_proposal.npz,odd_joint_Uaction.npz,bw_full_proposal.npz} > "$out/control-inputs.sha256"
python "$w/ck_compare.py" odd "$out/control12" --control | tee "$out/control12/compare.log"
sha256sum -c "$out/control-inputs.sha256"
# Work exclusively in reconstructed scratch, after baseline authentication.
python custody/ck/prepare.py
(cd "$cb" && python bw_propose_full.py) | tee "$out/trial13/proposal.log"
python "$cb/bw_prepare.py" | tee "$out/trial13/prepare.log"
python "$cb/bs_circle.py" 32768 0 | tee "$out/trial13/circle0.log"
python "$cb/bs_circle.py" 32768 1 | tee "$out/trial13/circle1.log"
python "$cb/bs_circle.py" 32768 2 | tee "$out/trial13/circle2.log"
python "$cb/bs_circle.py" 32768 3 | tee "$out/trial13/circle3.log"
python "$cb/bs_combine.py" 32768 | tee "$out/trial13/combine.log"
python "$cb/bv_joint_action.py" odd 4096 | tee "$out/trial13/action.log"
python "$w/ck_compare.py" odd "$out/trial13" | tee "$out/trial13/compare.log"
cp "$cb/"{bw_full_proposal.npz,bs_trial.npz,bs_circle.npz,odd_joint_proposal.npz,odd_joint_Uaction.npz,circle_part0.npz,circle_part1.npz,circle_part2.npz,circle_part3.npz} "$out/trial13/"
