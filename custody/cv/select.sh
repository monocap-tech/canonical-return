#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/cv/generated"
cb="$PWD/ce-output/pristine/cr_int_t1d_cb"
python custody/cv/select.py | tee "$out/prepare.log"
(cd "$cb" && python cv_select.py) | tee "$out/selection.log"
cp "$cb/"{odd_new.npy,cv_selection.npz} "$out/"
python "$cb/bq_corrector.py" 1 | tee "$out/corrector.log"
cp "$cb/odd_corrector.npz" "$out/"
python "$cb/bq_tail.py" 1 | tee "$out/tail.log"
cp "$cb/odd_U_tail.npz" "$out/"
