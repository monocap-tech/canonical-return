#!/usr/bin/env bash
# Source-preserving CR-ALP-033 11-dimensional acceptance / 12th-mode lower-bound obstruction.
set -euo pipefail
cd "$(dirname "$0")"
python3 cr_alp_033_initial_rational_data.py
printf 'bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2  CR-ALP-033_initial_rational_intervals_1e25.txt\n' | sha256sum -c -
g++ -std=c++17 -O2 -o cr_alp_033_partial_core cr_alp_033_partial_11_mode_interval_core.cpp
./cr_alp_033_partial_core CR-ALP-033_initial_rational_intervals_1e25.txt 30000 | tee CR-ALP-033_partial_core_success.log
printf '7d4458f97326f831e5ad52bff4f9e8d8f43bae67bb0e564079abfdca985a3ab8  CR-ALP-033_refined_H_exact_matrix_intervals.tsv\n' | sha256sum -c -
python3 cr_alp_033_independent_refined_ldl_audit.py | tee CR-ALP-033_independent_exact_audit.log
