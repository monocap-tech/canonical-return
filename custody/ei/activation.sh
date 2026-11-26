#!/usr/bin/env bash
set -euo pipefail
out="$PWD/custody/ei/generated"
w="$PWD/ce-output/pristine"
mkdir -p "$out/preparation" "$out/baseline" custody/eh/generated
python custody/ei/trial_prepare.py | tee "$out/preparation.log"
cp custody/ei/baseline.py "$w/ei_baseline.py"
python "$w/ei_baseline.py" odd "$out/baseline" "$PWD/custody/dz/generated/trial31" | tee "$out/baseline.log"
bash custody/eh/run.sh
bash custody/ei/run.sh
# Check that frozen evidence was not overwritten by receiver setup.
(cd custody/eg/generated && sha256sum -c sha256.txt)
