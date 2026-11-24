# CR-ALP-038 reproducibility

Scope: original signed seven-prime and complete Gamma physical form at a=11/10, pole-null, DRAFT/NONPROMOTING. The real-analysis coefficient identities and previous fixed-collar tail bound are unformalized analytical premises; executable scripts verify their exact-rational consequences. **This archive does not prove RH or general fixed-collar positivity.**

## Minimal standard-library exact checks

```bash
python3 cr_alp_038_combined_infinite_tail_rational_audit.py
python3 cr_alp_038_full_interior_infinite_schur_audit.py
```

Expected: first command certifies two explicit rational near-null lines versus ALL high sine modes, second command verifies the true source-signed five-by-five parity blocks after a rank-one-plus-diagonal upper bound on the ENTIRE high-tail coupling. The second checks all ten interior modes jointly. Both use exact rational or outward-rounded 10^-25 integer arithmetic; no numeric eigenvalues enter an acceptance comparison.

`H_mid.tsv` is the verbatim CR-ALP-034 source interval matrix (SHA-256 `b6a31e623298272fdceac788fc56626ec704e280903ff6baff1dc8bfe0d45f92`). `initial.txt` is the verbatim CR-ALP-034 source interval input (SHA-256 `bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2`). Both scripts fail if hashes differ.

## Independent precursor check

`cr_alp_034_independent_midpoint_ldl_verifier.py` in `predecessor_034/` reproduces the predecessor 12-mode signed Hmid and physical mass positivity with the original SHA-256 inputs. Run from the `predecessor_034` directory, which includes the pinned original names.

## Optional numerical diagnostic

```bash
python3 combined_residual_pilot.py
python3 schur_band_pilot.py
```

Requires NumPy and SciPy. Only band-limited floating-point observations; never an upper bound on the entire infinite tail. The diagnostic's `B_band <= B_full`, so positivity of `H-4B_band` is only a necessary test for `H-4B_full`. The two source-specific exact arithmetic certificates do **not** depend on this pilot.

## Certified analytical consequence (conditional on CR-ALP-027/034)

For every nonzero `p+q` with `p` in the ten-dimensional interior span of `w_3,...,w_12` and `q` in the exact pole-null high-sine tail `Q_15000` within the original logarithmic form domain:

`W_original(p+q,p+q) > 10^-12 ||p||_2^2 + (1/8)||q||_2^2`.

**Open:** both boundary anchors, the intermediate 12-to-15000 physical signed shell, original complete fixed-collar positivity, global extension, Lean and RH.
