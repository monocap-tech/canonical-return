# CR-ALP-068 — all-resolvent pole centering and full-strength diagonal

DRAFT / NONPROMOTING. Source: CR-ALP-066 exact Gamma rank-one boundary identity on Q85, inherited CR-ALP-026/027/042/045/046 classical identities, same physical source and original signed prime norm <=21/8. No imported Return-012 constants. No Lean certification. This is a diagonal refinement in lieu of an explicit block mixed-matrix solve; the global estimate controls all mixed Gamma boundary terms by weighted Cauchy-Schwarz. No witness promotion.

## 1. Exact centering at every resolvent parameter

In each parity write v_j=j*pi/(2a), a=11/10, beta_m=2m+1/2, f_beta(v)=v/(v²+beta²), w_j=v_j^(-1/2). The true pole moment is sum_j f_(1/2)(v_j)c_j=0. Fix V=120<v_j for all j>=86. Subtract alpha_beta f_(1/2) from each row with alpha_beta=(V²+1/4)/(V²+beta²). Algebraically

  f_beta(v)-alpha_beta f_(1/2)(v)
   = f_beta(v) [(beta²-1/4)/(V²+beta²)] [(v²-V²)/(v²+1/4)].

Define b_beta=(beta²-1/4)/(V²+beta²) and t_j=(v_j²-V²)/(v_j²+1/4), both in [0,1]. In particular beta=1/2 is annihilated exactly. This transformation changes no quadratic form on the pole-null space. Apply Cauchy-Schwarz with the original positive factors f_beta(v_i) w_i:

  |sum_j f_beta(v_j)c_j|²
  = |sum_j f_beta(v_j)b_beta t_j c_j|²
  <= (sum_i f_beta(v_i)w_i) sum_j f_beta(v_j) b_beta² t_j² |c_j|²/w_j.

This uses the centering factors squared on the second factor, not a potentially invalid independently centered row sum. All parity cross terms remain included. The argument extends from finite sums by absolute convergence of the linear functionals in L2, and nonnegative summation.

## 2. Split the beta-tail bound

From 066, the total of the small-beta contributions beta<32 is <= S0*t_j², where S0 is the saved small_pole_centered rational (<0.001319). Indeed b_beta<=(32/120)² and t_j<=1, so the same small-beta budget applies with the extra t_j². For beta>=32+1/2, use the 066 unimodal lattice bound sum_i f_beta(v_i)w_i <= (a/sqrt(2)+1/32) beta^(-1/2). With the original factor (1-epsilon exp(-2a beta))<101/100, replace the resulting prefactor by

  C=(2/a)(101/100)(a/(7/5)+1/32).

The nonnegative sum in the large-beta diagonal correction is bounded by C t_j² sum_beta b_beta² beta^(-1/2) v_j^(3/2)/(beta²+v_j²). The unweighted-in-b sum is <= I= (22/7)/(2*(7/5)) by the same decreasing step-two integral argument as in 066. For a split B>=64, b_beta² <= A_B=(B²/(V²+B²))² on beta<=B and <=1 otherwise. The beta>B tail of the unweighted-in-b sum is <= min(I,v_j^(3/2)/(3*(B-2)^(3/2))), since each lattice term is bounded by half the integral over its preceding length-two interval, and beta^(-1/2)/(beta²+v²)<=beta^(-5/2). The first preceding interval lies above B-2. Consequently

  penalty_j <= t_j² [ S0 + C (A_B I + (1-A_B) min(I,v_j^(3/2)/(3*(B-2)^(3/2)))) ].

This is a valid upper bound for each j, simultaneously, even if B is selected separately for different frequency bands. No omitted beta or sine tail. For finite bands [j,J], use v_j >= j*(314159/100000)/(2a) to lower bound g, and v_J <= J*(22/7)/(2a) to upper bound t and v^(3/2). Use sqrt upper and lower via integer roots at scale 10^12. On the final infinite band use t<=1, b<=1 and S0+C I, i.e. the prior exact total <1.7. The code uses the minimum of the split candidates and the old universal budget before multiplying by the t upper squared.

## 3. Executed exact arithmetic

Audit src/audits/cr_alp_068_pole_centered_gamma.js reads the pinned 066 constants, verifies their exact algebraic recombination, evaluates 118 contiguous bands from j=86 to infinity, and stores exact rational endpoints in src/audits/CR-ALP-068_pole_centered_gamma_result.json. Log lower uses the same range-reduced positive atanh series and 1e-9 downward rounding as 066. The band Gamma lower D_j uses the inherited g(v)>log(k)-37/20-2/k for v>=k. The numerical split selection is not used: all split candidates are exact rational bounds.

The minimum certified D_j-21/8 is >0.265358630960 (near 226..229), hence minimum Gamma D_j>2.890358630960 > old gamma=2.88. Unlike 066, the new diagonal dominates the previous constant Gamma floor throughout Q85, not just at large j. The full-strength comparator C=Q85[D+P]Q85 is coercive without theta blending. It also dominates the earlier constant comparator in form order. For the explicit step bands:
- j=86..89: d>0.302590810741.
- j=851..960: d>0.990318735107.
- j=1251..1536: d>1.346374247468.
- j=2049..4096: d>1.822938456459.
- j=4097..8192: d>2.513943823640.
- j=8193..16384: d>3.206910683872.
- j>=16385: d>3.900012002176.

The early weights are not monotone; for a given residual cutoff use the minimum of all remaining band weights, not blindly the next band. The later listed weights do increase. The final constant cap is optional: retaining the original increasing symbol estimate with universal penalty gives a growing uncapped diagonal on the last band. The finite step comparator itself is bounded and coercive.

## 4. Relation to CR-ALP-067

The old fixed weighted response had an effective floor requirement >1.228745670493; the new post-2048 floor >1.822938456459 removes the old cap obstruction as a universal limitation. It does NOT establish that the old residual is unchanged or passes. Changing D changes the response energy and residual. For the exact old rational trial x, one may transport the existing response logs algebraically without quadrature: Delta=(D_new-D_old)x is a finite ordinary sine sum; E_new=E_old-<x,Delta>, y_new=y_old-Delta, ||y_new||²=||y_old||²-2 Re<y_old,Delta>+||Delta||², and <y_new,b>=<y_old,b>-<Delta,b>. Update each finite sine coordinate, then apply the same exact Q85 projection and Parseval tail checker with the new band weights. Use actual original trial bank, not a recomputed approximation. The new source proxy error transport uses the new global coercivity. If this test fails, select a new trial on the stronger comparator and certify it; no reason to repeat fifth-moment integration.

Executed BigInt arithmetic and exact source/output readback at fe026af6f57f6b5daa8a6cc9c4d1fb521b812cf4. This verifies numerical constants, not formal checking of the analytic identities. Weighted 067 native result remains UNDECIDED; constant response still running and original weighted trial bank recovery queued at start of this turn. No new computational jobs submitted. Odd witness, 1048 positive dimensions and conditional nonpositive-index cap 66 unchanged. Historical standalone archive reconciliation open.
