#!/usr/bin/env bash
set -euo pipefail
bash custody/ct/select.sh
bash custody/ct/action.sh
bash custody/ct/trial.sh
