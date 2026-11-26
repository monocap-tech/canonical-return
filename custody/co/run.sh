#!/usr/bin/env bash
set -euo pipefail
bash custody/co/select.sh
bash custody/co/action.sh
bash custody/co/trial.sh
