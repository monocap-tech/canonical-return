#!/usr/bin/env bash
set -euo pipefail
bash custody/dg/select.sh
bash custody/dg/action.sh
bash custody/dg/trial.sh
