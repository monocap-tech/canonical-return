#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 src/audits/canonical_return_t1b_exact_audit.py
python3 src/audits/canonical_return_t1c_rational_enclosures.py
python3 src/audits/canonical_return_t1d_uniform_profile.py
python3 src/audits/canonical_return_t1d_b_25_collision.py
python3 src/audits/canonical_return_t1d_c_8_collision.py
python3 src/audits/canonical_return_t1d_d_31_schur.py
python3 src/audits/canonical_return_t1d_e_mixed_bound.py
python3 src/audits/canonical_return_t1d_f_parity_decoupling.py
python3 src/audits/canonical_return_t1d_g_sharp_rate.py
python3 src/audits/canonical_return_t1d_h_mellin_probe.py
python3 src/audits/canonical_return_t1d_i_residual_bracket.py
python3 src/audits/canonical_return_t1d_j_complement_probe.py
