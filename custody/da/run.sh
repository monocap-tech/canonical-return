#!/usr/bin/env bash
set -euo pipefail
bash custody/da/select.sh
bash custody/da/action.sh
bash custody/da/trial.sh
