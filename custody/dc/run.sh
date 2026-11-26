#!/usr/bin/env bash
set -euo pipefail
bash custody/dc/select.sh
bash custody/dc/action.sh
bash custody/dc/trial.sh
