#!/bin/sh
# Reconstruct original signed eight-mode Weil interval certificate.
set -eu
cd "$(dirname "$0")"
python3 cr_alp_030_initial_rational_data.py
g++ -std=c++17 -O2 -o cr_alp_030_interval_core cr_alp_030_eight_mode_interval_core_fast.cpp
./cr_alp_030_interval_core CR-ALP-030_initial_rational_intervals_1e25.txt 30000
python3 cr_alp_030_independent_ldl_verifier.py
