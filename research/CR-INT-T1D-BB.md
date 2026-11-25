# CR-INT-T1D-BB — Exact signed arithmetic quadratic metric

**Status:** analytical refinement and executed directed-integer metric pilot; no new full-sector spectral inertia. BA's conditional bound 1188 (608 even + 580 odd) and complement margin 3/2000 remain inherited. Fixed window (-2,2), original smooth double Mellin-pole-null source, all 24 signed prime powers. No Lean, remote CI or independent review.

## 1. Exact normalization

Reopen BA, AY, AX. Let lambda=399/50, J0=(lambda I-A)/4, B=(24/5)J0, c=8/35. P synthesizes restricted circular modes. In the joint inverse variational inequality take arithmetic trial v_B=B^(1/2)P Z/21. Its linear cross is c P*J0(T-RV), and its quadratic metric is K_P=P*(B+B^2)P/441. Writing D1=P*J0P and D2=(J0P)*(J0P), the normalized metric is K=K_P/c=D1/21+c D2. The old normalized metric D1 was an upper bound, not an equality. AM gives 0<=B<=(2394/125)I, so even without evaluating D2, K<=(2519/2625)D1< D1 in the nonzero directions. With BA's V'=V-cYZ, the effective normalized metric is K-c Q_Y, Q_Y=GY+Y*G*-Y*HY. In actual (unfactored) quadratic units this is K_P-c^2 Q_Y. The ridge is added only after this substitution. An arbitrary Y,Z still gives an upper comparator by the original variational inequality; no trusted solver is required.

## 2. Piecewise exact integration

A is the sum of 48 signed shifts with coefficients log(p)/sqrt(n) for the 24 prime powers n<=54. Restricted shifts change activity at the 48 points +/- (2-log n); with +/-2 these make 49 cells, whose strictly ordered endpoints are certified by rational log intervals. On each open cell J0 u_j=a_j cos(omega_j x)+b_j sin(omega_j x), omega_j=j*pi/8. Here u_j=cos(omega_j x)/2 or sin(omega_j x)/2, with j positive and fixed parity. Coefficients follow directly from addition formulae and exact activity tests. D2 is the sum of all four products of cellwise coefficients times exact elementary integrals of cos-cos, sin-sin, sin-cos and cos-sin. Use sum/difference frequencies and endpoint antiderivatives, with the zero frequency handled as cell length. Endpoint phase enclosure uses AX's rational Taylor and radix phase machinery; grid 10^9; split integer products, directed floors and ceilings. No quadrature or floating-point decisions in the metric certificate. D1 is computed independently by the same cell method and compared to the AX directed enclosure (using each entry's row/column error upper bound). Constant coordinate j=0 remains outside this evaluator. Pilot IDs are stored in bb_metric_{0,1}.npz, selected as 16 spaced low deficit indices and 48 spaced high deficit indices >10000 from inherited arrays. Selection is exploratory, not optimal or identical to AZ's 640 selected carriers.

## 3. Executed pilot

AQ/AS array hashes and AW manifest pass; inherited BA log checker passes all four recorded rational spectral checks (not a re-execution of the full spectral solve). The BB evaluator reran AM's arithmetic norm and high-tail check. Both 64-carrier parity runs passed; each used 49 cells and 4096 independent AX entry-overlap checks. Max full row interval width for K: even 5529/10^9; odd 39/7812500. Gershgorin lower bound on D1-K: even 116830093/10^9, odd 95586049/10^9.

Exact fraction-free Bareiss PD tests certify (13/20)D1-K positive definite on the even pilot set and (4/5)D1-K positive definite on the odd pilot set. For each rational factor a/b, form directed integer enclosure of a D1-b K, take symmetric integer midpoint, subtract maximal interval-error row sum on diagonal, and check all 64 leading principal minors by exact Bareiss elimination. Smaller trial factors that did not pass are not evidence of indefiniteness. These are relative metric certificates on the declared pilot sets ONLY, not full 640-carrier or original Weil positivity tests.

## 4. Reproduction and blocker

The inherited BA archive contains scripts, source arrays and logs, not az_carr.npy, az_H.npy, az_V.npy, ba_solve_checkpoint.npz or ba_metric_input.npz. Therefore no full BA checkpoint was available to reoptimize directly. The following executed locally: python bb_metric.py; python bb_relative.py; python ba_check.py; sha256sum -c array_sha256.txt; sha256sum -c aw_sha256.txt. The new arrays, logs and their manifest are bundled. BB full adapter bb_build.py and bb_projected.py have passed Python syntax compilation ONLY. They are not yet numerically exercised. Rebuild AZ and BA per their notes, preserve the actual accepted carrier order and input provenance, then run:

```
python bb_metric.py az_carr.npy
python bb_build.py
RIDGE_DENOM=64 python bb_projected.py
```

bb_build.py uses the same BA Y, E, F, G, H, C, V checks and computes the effective metric directly from K instead of adding a correction with duplicated D1 uncertainty. The BB metric enclosure error comes from directed rescaling to the 2^26 grid and maximal error row sum. The projected auditor retains BA's rational subspace acceptance check. This is a prepared path, not a new certificate. Need also check other ridges if the proposed matrix is ill-conditioned.

## 5. Next — BC

Regenerate the missing full AZ/BA checkpoint and execute the 640-carrier exact metric and spectral verification. If resources are limiting, store derived checkpoints rather than only their hashes to avoid repeated dense reconstruction. Compare against BA's odd 580 constraint count at beta<997/1000; even 608 inherited unless independently improved. Do not treat pilot metric reduction as a codimension reduction. Full fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.