#!/usr/bin/env bash
set -euo pipefail
bash custody/dn/select.sh
bash custody/dn/action.sh
bash custody/dn/trial.sh
