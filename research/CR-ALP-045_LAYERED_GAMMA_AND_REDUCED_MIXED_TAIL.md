# CR-ALP-045 — Layered Gamma lower bounds and a smaller mixed tail

DRAFT / NONPROMOTING. New exact arithmetic and classical analytic arguments conditional on CR-ALP-026/027/042 and CR-ALP-041/043/044 continuum and interval premises. No Lean, full-collar or RH acceptance.

## Source contract

Use the same original signed `Zeta23.EF.literatureRHS` (Weil-Lab blob `0034976fe2f1e5fa9d411108736160542eb94293`), interval I=(-11/10,11/10), physical L² metric, both pole moments zero, seven negative paired shifts at 2,3,4,5,7,8,9 and full Gamma symbol g(r)=Re psi(1/4+ir/2)-log pi. Q_N is the genuine pole-null ordinary high-sine subspace, the physical orthogonal complement of E_N. Return-012's exact proof and audit were not retrieved; none of its numerical weights or its prime bound are used. The original CR-ALP-042 support-aware ||P||<21/8 proof and arithmetic are reused and rerun. The attached NF662 congruence discussion is a different finite block, not an identification with this source family.

## 1. Layered Fourier leakage

Let J=N+1, p=314159/100000<pi by Machin alternating arctangent bounds, v=J p/(2a). In each ordinary sine parity separately, q=sum_{j>=J} c_j e_j, with the parity pole moment sum c_j nu_j/(nu_j²+1/4)=0. CR-ALP-027's exact subtraction formula bounds the squared Fourier transform at |r|<v by

    (4/a) (r²+1/4)² (1-r²/v²)^(-2) (2a/p)^6
       * sum_{j>=J, parity} j^(-6) * ||q||².

For the mesh-two sum, the integral test gives <= J^-6 + 1/(10 J^5). Opposite reflection parities have orthogonal Fourier transforms on every symmetric frequency interval, including complex coefficients, so the same bound on integrated leakage holds for arbitrary q in Q_N. The cumulative probability bound d_R is the sum over unit subintervals [k-1,k] of

    [4/(a p)] (2a/p)^6 [J^-6+1/(10J^5)]
    * [F(k)-F(k-1)] / [1-(k/v)^2]^2,
    F(r)=r^5/5+r^3/6+r/16.

Each term is rounded UP to a multiple of 10^-12, using exact rational arithmetic. The integral formula uses the numerator antiderivative and monotonic denominator bound on each interval. All evaluated endpoints are strictly below v. No numerical quadrature is asserted.

Take L_1=-6 and L_k=floor(10^9 log_lower(k))/10^9 -37/20-2/k for k>=2. The rational log lower uses 24 positive atanh terms after powers-of-two reduction. CR-ALP-026 gives g(r)>L_k whenever |r|>=k and g>-6 globally; the increasing staircase is a pointwise lower bound regardless of any monotonicity of g. The audit checks the increments are positive. Plancherel and the finite layer-cake expansion therefore give

    W(q,q) >= [-6 + sum_{k=2}^T (L_k-L_(k-1)) max(0,1-d_k) -21/8] ||q||².

The summation stops at the first d_k>=1, or last safe endpoint below v. New exact margins: Q66 >1/100000 (raw lower about 0.00001314628); Q67 >1/100; Q85 >1/4; Q100 >2/5; Q123 >3/5; Q172 >9/10. Q66 implies a nonpositive index at most 66 on the fixed pole-null collar. No finite E66 positivity or mixed Q66 result follows. Q85 has a much stronger tail margin and its own unevaluated Schur gate H85-4B85 PSD.

## 2. Independent simultaneous weighted comparison

For any integer 0<=s<5 apply weighted Cauchy-Schwarz to the exact cancellation, replacing ||q||² by sum (J/j)^s |c_j|² and the reciprocal-power sum by J^-s sum j^(s-6). The latter is bounded by J^-6 + 1/[2(5-s)J^5] on each parity. Let d_k^(s) be the resulting integrated upper bound. Without clipping the lower probability estimate (which can be negative), the layer-cake proof gives

    W(q,q) >= sum_{j>N} [L_T-21/8 - C (J/j)^s] |c_j|²,
    C=sum_{k=2}^T (L_k-L_(k-1)) d_k^(s).

For N=123,T=171,s=2 the exact baseline exceeds 13/20 and C<1/20. Thus the simultaneous diagonal lower weights are d_j=13/20-(124/j)^2/20, j>=124. These are not the Return-012 weights and do not reproduce its logarithmic growth; they are a modest increasing comparison with a finite limiting weight. For a residual r in physical high sine coordinates, the weighted inverse-form upper comparison is min_t sum_{j>N}|r_j-t c^pole_j|²/d_j in each parity. This follows by weighted Cauchy-Schwarz and minimizing over legal pole subtractions; no weighted residual numerical test is claimed here.

## 3. Downward fitted-band extension

Reuse the saved CR-ALP-041 rational energy basis (sigma=10^-14), CR-ALP-044 fitted Gram above j=172 through 30000, exact same rational fit and analytic remainder beyond 30000. The fitted coefficients are centered by a multiple of the true pole function; keeping the OLD centering and fit is legal for any smaller Q_N. For new N<172 simply add interval outer products from the newly included matching-parity j in (N,172]. No re-fit or finite-band-only assertion. The stored remainder bound is unchanged because its source and centering are unchanged. CR-ALP-043's digamma routine initially rejected new smaller frequencies: its atan argument-domain check was valid. Repair via the exact recurrence psi(z)=psi(z+n)-sum_{m=0}^{n-1}1/(z+m), n=ceil(3 Im z), when Im z<120. The inherited real-shift Euler-Maclaurin bound and atan interval domain then apply. No floating arithmetic participates in acceptance. Both predecessor basis and old saved interval-Gram check were rerun, as was CR-ALP-042's prime audit. The old matrix, basis and coefficient-generator SHA256 are pinned in the extension script.

Test each saved interval Gram against the first passing rational hundredth-grid scalar b I, using midpoint/radius LDL. Independent saved-matrix check uses integer Bareiss leading determinants of the row-radius-shifted comparison. Add the unchanged analytic remainder to b. The exact total rho bounds:

| N | odd-index rho | even-index rho | scalar lower mu (rounded for display) | outcome |
|---|---|---|---|---|
| 85 | 1.040537 | .427375 | .2574647 | insufficient |
| 100 | .670537 | .357375 | .4225946 | odd insufficient |
| 110 | .530537 | .297375 | .5192404 | odd insufficient |
| 111 | .530537 | .297375 | .5284470 | odd insufficient |
| 112 | .530537 | .287375 | .5373057 | PASS |
| 123 | .400537 | .257375 | .6323865 | PASS |
| 172 | .210537 | .187375 | .9709917 | PASS |

The energy-normalized residual bound ||Q_N Tp||² <= rho ||z||² and the low bound W(p)>=sigma||p||²+||z||² imply W(p+q)>=sigma||p||²+(1-rho/eta)||z||²+(mu-eta)||q||². Select strict rational mu and eta: N112 mu=537/1000, eta=107/200; N123 mu=3/5, eta=1/2; N172 mu=9/10, eta=23/100. Exact tests pass and yield strict mixed inequalities with tail margins 1/500, 1/10, 67/100 respectively, for nonzero p+q. Q112 is the smallest successful *tested symmetric mixed* cutoff in this run, not a proven optimum. It leaves a 98-dimensional intervening shell. The standalone Q66 index bound is stronger than the codimension 98 mixed bound, but describes a different positive sector. The prior CR-ALP-044 mixed result at Q172 is strengthened.

## 4. Reproduce and limits

From the package directory: `python cr_alp_041_energy_conditioned_schur_audit.py --verify`, `python cr_alp_042_support_aware_prime_tail_audit.py`, `python cr_alp_044_check_result.py`, `python cr_alp_045_layered_gamma_audit.py`, `python cr_alp_045_extend_fitted_band.py`, `python cr_alp_045_check_result.py`. All executed. The independent determinant verifier checks saved matrices and final inequalities; it does not independently reconstruct every analytic coefficient or old interval input. No Lean compilation, Return transfer, F41 family identification, full finite E66/E85/E112, full-collar or RH claim. Next: evaluate a complete low block with the Q85 or Q112 tail comparison; avoid mistaking smaller codimension for positivity of the missing directions.
