# CR-INT-T1D-AU — Signed exterior Gram and iterated inverse refinement

**Status:** local directed-grid spectral audit PASS both parities, conditional on AP/AM original source comparisons, AQ/AS certified inherited arrays and AR ordinary binary64 exact-integer-dot arithmetic model. Not Lean, independent review, remote CI, whole positivity, or RH.

## Exact variational refinement

Keep AT's physical interval, pole constraints, low circular exclusions, baseline t=1/2, positive surplus, and complete original signed source lower comparator. Write B=T*(E*DE)^-1 T, M=T*T, C for the AS cross Gram, K=CC*, Q=D^-1, F=I-EE*, W=diag(w_j^2), and L=W^(1/2)K. The exterior synthesis Z=F Phi_J diag(sqrt(b_j)w_j) obeys Z*QZ=H=W^(1/2)(diag(b)-M-K)W^(1/2), 0<=H<=I, and Z*QET=L, by AT. Completing the square and testing the exterior dual gives for any trial matrix V

    B <= M-K - V*L - L*V + V*H V.

Set J0=I-H and P_r=sum_(k=0)^(r-1) J0^k, V_r=P_r L. Since HP_r=I-J0^r and P_r commutes with H, 2P_r-P_r H P_r=sum_(k=0)^(2r-1) J0^k. Thus the upper comparators M-K-L*(sum_(k=0)^(2r-1) J0^k)L decrease monotonically with r, without assuming H invertible. The actual audit uses r=2 as a floating proposal only; its rounded V is an arbitrary valid trial.

## Directed-grid implementation and arithmetic

Reuse the hash-checked AQ/AS arrays, not a new frequency scan. S=2^25, U=S/2. AR/AS interval routines bound M-M0/S by eta_M and C-C0/S by rho. Since max b<10, ||C||<sqrt(10)<7/2; e_K=rho(7+rho)+m/S bounds K-K0/S. Round w_j downward so b_j w_j^2<=1. Form H0 by rounding diag(b), then each of the two diagonal weight multiplications to grid S; symmetrize by entrywise minimum (each entry's rounding error <=2/S relative to the symmetric weighted input). Hence e_H=eta_M+e_K+(2m+1)/S. L0= floor(diag(w) K0)/grid has e_L=e_K+m/S. The trial proposal uses (2I-H0/S)L0/S, rounded on U grid. Let N bound its norm by the exact row/column norm inequality. HV is computed by exact integer dot product and rounded to S grid (error <=m/S). Form the trial gain by symmetric averaging of the two transposed quadratic-product approximations, then floor to S grid. Symmetric averaging does not enlarge the operator norm error. The resulting gain error is at most

    2N e_L + N^2 e_H + N m/S + m/S.

Add eta_M+e_K for the complete upper comparator enclosure. AR's exact-dot assertions ensure Cauchy bounds on absolute dot products below 2^53; assume ordinary binary64 multiply/add or FMA, not reduced precision. Cross-Gram accumulation checked for int64 bounds. Eigenvectors and eigenvalues are only proposals; the independent rational-grid Gram/residual check yields beta=(1+delta)(L_lambda+epsilon+eta)/(1-delta). The first exploratory even run used entrywise minimum rather than symmetric averaging for a nonsymmetric intermediate gain, so is NOT relied upon. Both final executions use corrected symmetric averaging.

## Executed results

| parity | deficit dimension | absorbed dimension | constraints incl. 512 | beta diagnostic |
|---|---:|---:|---:|---:|
| even | 4628 | 4022 | 1118 | .9948191842518328 |
| odd | 4691 | 4076 | 1127 | .9948661021245951 |

The exact beta fractions in au_check.log are strictly below 997/1000. Thus on the certified physical complement q0[f]>=3/2000 <E*DE f,f> >=3/2000 ||f||², conditional negative index <=2245 (320 fewer than AT). Initial 1024 low circular constraints retained. No earlier mixed dimensions subtracted. No original negative witness, fixed-window or global positivity proved.

Reproduce in this directory with Python 3, NumPy, SciPy: `sha256sum -c array_sha256.txt`; optionally regenerate via `python aq_scan.py` and `python as_surplus.py`; `OPENBLAS_NUM_THREADS=4 python au_exterior.py 0 2 > au_even.log; echo $? > even.exit`; corresponding parity 1 to au_odd.log and odd.exit; `python au_check.py`. The log checker confirms the particular executed results, not arbitrary new proposals. The analytical identities are not kernel certified.

## Next gate — AV: exact positive-surplus inverse on physical support

There is a stronger alternative to iterating only on AT's deficit-index exterior span. Use the surplus synthesis R=E*Phi_P diag(sqrt(alpha)), and Hs=diag(1-alpha)+R*R. Then E*DE=I+R diag(1-alpha)^-1 R*, and the finite Woodbury identity gives EXACTLY B=M-C Hs^-1 C*, with C=T*R. This retains the entire finite surplus comparator on physical support, not just its projection into the exterior deficit span. Parseval gives R*R<=diag(alpha), so Hs<=I. The inherited array check proves max alpha<23/24 and max b<10, hence Hs>=I/24. Therefore J_s=I-Hs is PSD with norm <=23/24 and

    B <= B_r=M-C (sum_(k=0)^(r-1) J_s^k) C*,
    0<=B_r-B<=240(23/24)^r I.

The exact-fraction checker verifies the remainder <1e-7 at r=512 and <2e-12 at r=768. These are truncation bounds only, not computed matrix inertia or rounding certificates. In practice propose a solve for Hs with residual or polynomial acceleration and certify its gain using the same arbitrary-trial variational principle. Requires the surplus-to-surplus physical Gram, not yet evaluated. It closes the auxiliary inverse approximation in principle, not original Weil positivity. Full source positivity, completed Suzuki–Connes transport, Phase IV and RH remain OPEN.