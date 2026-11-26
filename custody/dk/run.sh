#!/usr/bin/env bash
set -euo pipefail
bash custody/dk/select.sh
bash custody/dk/action.sh
bash custody/dk/trial.sh
