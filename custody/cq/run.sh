#!/usr/bin/env bash
set -euo pipefail
bash custody/cq/select.sh
bash custody/cq/action.sh
bash custody/cq/trial.sh
