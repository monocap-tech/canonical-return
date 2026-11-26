#!/usr/bin/env bash
set -euo pipefail
bash custody/dv/select.sh
bash custody/dv/action.sh
bash custody/dv/trial.sh
