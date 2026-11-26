#!/usr/bin/env bash
set -euo pipefail
bash custody/du/select.sh
bash custody/du/action.sh
bash custody/du/trial.sh
