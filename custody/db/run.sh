#!/usr/bin/env bash
set -euo pipefail
bash custody/db/select.sh
bash custody/db/action.sh
bash custody/db/trial.sh
