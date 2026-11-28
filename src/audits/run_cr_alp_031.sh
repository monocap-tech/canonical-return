#!/bin/sh
# Full reproducibility: original source intervals, 30k Gamma terms,
# 617 rigorous PSD Loewner endpoint groups, and independent rational LDL.
set -eu
cd "$(dirname "$0")"
python3 cr_alp_031_initial_rational_data.py
g++ -std=c++17 -O2 -o cr_alp_031_interval_core cr_alp_031_ten_mode_interval_core.cpp
./cr_alp_031_interval_core CR-ALP-031_initial_rational_intervals_1e25.txt 30000
python3 cr_alp_031_monotone_group_schedule_audit.py
python3 cr_alp_031_independent_ldl_verifier.py
