# CR-ALP-059 — actual global second prime moment (candidate, nonpromoting)

The prerequisite archive is the repository-controlled 90-part CR-ALP-058 ZIP, SHA-256 f5f69ea089b2a30a84d4d64be8f2c37a95a468e99e4897f24761836abae44bd2. No replacement witnesses, no frequency truncation, no mathematical promotion.

Given the retained 058 projected physical proxy r_p and exact rational shifts s_j, integrate v(t) = -sum_j w_j [r_p(t+s_j)+r_p(t-s_j)] with zero extension. Since v has the same reflection parity as r_p, integrate twice over [0,11/10]. For each panel the rational mesh splits all translated original knots, checks all active argument pieces and inherited complex disk conditions. The absolute complex v bound is <8*10^40, so each integrand (v^2, v e_i, v b) is bounded by 64*10^80. Exact degree-191 Gauss remainder <1e-12 follows by the same direct Cauchy argument as 058. Point enclosures use the unchanged 058 midpoint evaluator with the original forward-error bounds, outward binary128 arithmetic and MPFR basis sine intervals. This is conditional on those inherited arithmetic and continuum arguments; not a Lean proof.

The pole-null and ordinary low-sine projector gives N_p = ||v||^2 - sum_i <v,e_i>^2 - (<v,b>-sum_i c_i<v,e_i>)^2/st, with st=||b||^2-sum_i c_i^2>0.003. The second projection is scalar interval arithmetic, not a finite high-sine compression. The source reconstructs and rechecks the first Q85 projection and its proxy norm <2/5.

For the true N, let eps be the inherited 058 global proxy error, s=21/8, delta=s*eps+1e-13. The 1e-13 is the inherited rational-shift displacement bound, applied to the projected physical proxy, whose amplitude/derivative/jump estimates were documented in 058 VALIDATION.md. Then ||Q P r-Q P_rat r_p||<=delta and ||Q P_rat r_p||<=s*(2/5)+1e-13. The squared-norm transport is 2*(s*(2/5)+1e-13)*delta+delta^2; this does not require a separate spectral norm bound for the rational-shift operator.

## Fresh checkout replay

```sh
ROOT="$PWD"
WORK=$(mktemp -d)
python3 "$ROOT/src/audits/cr_alp_058_archive/reassemble.py" "$WORK/package.zip"
python3 -m zipfile -e "$WORK/package.zip" "$WORK"
(cd "$WORK/cr_alp_058" && sha256sum -c SHA256SUMS && python3 validate_geometry.py && python3 check_saved.py)
python3 "$ROOT/src/audits/cr_alp_059/mesh_second.py" "$WORK"
g++ -O2 -std=gnu++17 -I"$WORK/cr_alp_058" "$ROOT/src/audits/cr_alp_059/second.cpp" -o "$WORK/second" -lquadmath -l:libmpfr.so.6 -lgmp
for k in 1 2; do
  "$WORK/second" "$WORK/cr_alp_058/input_$k.txt" "$WORK/cr_alp_058/mesh.txt" "$WORK/cr_alp_058/gauss.txt" "$WORK/cr_alp_058/second_mesh.txt" > "$WORK/second_$k.log"
done
python3 "$ROOT/src/audits/cr_alp_059/check_second.py" "$WORK" "$WORK"
```

Acceptance is scoped only to the two frozen comparison witnesses; a failure of the comparison is not an original Weil negative direction. Neither the full 85-column moment bank nor the infinite Schur positivity theorem is established here. Until actual execution and saved output verification, N remains unevaluated.

Custody gate: commit source, immutable input hash and replay command before execution; archive the resulting logs and exact endpoint JSON on the same draft branch, not only in ephemeral storage. Historical standalone archive reconciliation remains open as recorded in the investigation-wide custody inventory.
