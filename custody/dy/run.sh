#!/usr/bin/env bash
set -euo pipefail
python custody/dy/trial_prepare.py | tee custody/dy/generated/trial-prepare.log
python custody/dy/screen.py | tee custody/dy/generated/screen.log
