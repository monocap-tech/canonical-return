#!/usr/bin/env bash
set -euo pipefail
mkdir -p custody/ek/generated/preparation
python custody/ek/trial_prepare.py | tee custody/ek/generated/preparation.log
python custody/ek/screen.py | tee custody/ek/generated/screen.log
