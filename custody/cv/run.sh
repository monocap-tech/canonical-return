#!/usr/bin/env bash
set -euo pipefail
bash custody/cv/select.sh
bash custody/cv/action.sh
bash custody/cv/trial.sh
