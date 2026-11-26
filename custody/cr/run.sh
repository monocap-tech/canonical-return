#!/usr/bin/env bash
set -euo pipefail
bash custody/cr/select.sh
bash custody/cr/action.sh
bash custody/cr/trial.sh
