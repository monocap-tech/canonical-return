#!/usr/bin/env bash
set -euo pipefail
bash custody/cy/select.sh
bash custody/cy/action.sh
bash custody/cy/trial.sh
