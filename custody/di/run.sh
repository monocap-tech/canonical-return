#!/usr/bin/env bash
set -euo pipefail
bash custody/di/select.sh
bash custody/di/action.sh
bash custody/di/trial.sh
