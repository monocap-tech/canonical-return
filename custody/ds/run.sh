#!/usr/bin/env bash
set -euo pipefail
bash custody/ds/select.sh
bash custody/ds/action.sh
bash custody/ds/trial.sh
