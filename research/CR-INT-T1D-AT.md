# CR-INT-T1D-AT — Exterior-constrained positive-energy inverse

**Disposition:** analytical refinement and executed local directed-grid certificate PASS in both parities, conditional on AP/AM original form comparisons, AQ/AS inherited certified arrays, and AR's ordinary binary64 exact-integer-dot arithmetic model. No independent review, Lean, remote CI, negative original witness, whole positivity, completed Suzuki–Connes transfer, Phase IV or RH.

## 1. Exact exterior penalty

Preserve I=(-2,2), the original physical L² norm, smooth two-Mellin-pole-null core, AP's parity-circle zero extension E, and AS's baseline t=1/2. Use the same deficit indices J and upper weights b_j, positive-surplus indices P disjoint from J, and rational alpha_k from AS. Set D diagonal on the circle, equal to 1/(1-alpha_k) on P and 1 elsewhere. Since alpha_k <= s_k/(t+s_k), this D is no larger than the actual AS surplus comparator, and

    q0[f] >= t(<E*DE f,f> - ||T*f||²),
    T e_j = sqrt(b_j) E*phi_j,

on the source with all j<1024 coefficient functionals removed. Define Q=D^-1=I-sum_(k in P) alpha_k phi_k tensor phi_k, and F=I-EE* (orthogonal exterior projection). For a physical g, completion of the square gives the exact formula

    <g,(E*DE)^-1 g> = <Eg,Q Eg> - inf_f ||D^(1/2)(Ef-QEg)||².

D>=I and F Ef=0 imply that the infimum is at least ||F Q Eg||². Hence

    (E*DE)^-1 <= E*Q E - E*Q F Q E.

This is a true exterior-constraint improvement over AS's inverse compression, not an assumed equality of the two inverses.

## 2. Reuse the AS cross Gram

Let M=T*T, C_jk=sqrt(b_j alpha_k)<u_j,u_k>, K=CC*, u_j=E*phi_j. Then T*E*Q E T=M-K. Because J and P are disjoint and E Tz is supported inside I, the exterior coefficient at deficit mode j satisfies

    <phi_j,F Q E Tz> = (Kz)_j / sqrt(b_j)

in a consistent linear-coordinate convention (all Gram entries real in each parity). Choose rational 0<=w_j<=min(1,b_j^(-1/2)), and W=diag(w_j²). Bessel's inequality on the circle yields

    ||F Q E Tz||² >= sum_(j in J) b_j w_j² |<phi_j,F Q E Tz>|² = ||W^(1/2) Kz||².

Thus the finite signed upper comparator is

    B=T*(E*DE)^-1 T <= M-K-K W K.

The upper comparator is PSD by this derivation; no additional approximate spectral premise is used. AS's absorption argument now applies with E*DE as the energy metric. On a certified subspace of the deficit coordinate space with upper form bound beta<1, impose T*f in that subspace. Codimension is at most |J|-dim(W_cert)+512 per parity, where W_cert is the certified subspace, distinct from the weight matrix W. Then q0[f]>=t(1-beta)<E*DE f,f>>=t(1-beta)||f||². No previous mixed directions are subtracted.

## 3. Executed fixed-point audit

`at_exterior_v2.py` imports AR's integer overlap, square-root, norm and exact-dot routines and AS's cross-interval construction. It reuses the supplied AQ/AS arrays, checked against the preceding bundle's SHA256 values; does not rerun the AQ frequency scan or AS surplus generation. The original AS 24-prime symbol, Gamma comparison and surplus constraints are unchanged.

S=2^24. Compute the integer matrices M0 and K0 on S grid using AS's directed cross enclosure; e_K=rho(8+rho)+m/S bounds ||K-K0/S||. For each b_j=B_j/10^8, set w_j=min(S,isqrt(floor(S²*10^8/B_j)))/S; verify b_j w_j²<=1 with Python integers. Let L=W^(1/2)K and L0 its approximate factor obtained by rounding diag(w)K0/S downward onto the S grid. Then e_L=e_K+m/S bounds ||L-L0/S||. A checked row/column norm bound N_L>=||L0/S|| supplies the much tighter error

    ||L*L-L0*L0/S²|| <= e_L(2 N_L+e_L).

The additional rounding of the exact integer factor product onto S grid adds m/S; total comparator error eta=eta_M+e_K+e_L(2 N_L+e_L)+m/S. The executable checks all integer dot products using AR's Cauchy absolute sum bound <2^53, conditional on ordinary binary64 multiply/add or FMA with no reduced precision. The approximate eigensolver only proposes a subspace; its rational-grid Gram and residual are independently checked with the same exact-dot bounds. The accepted beta is (1+delta)(L_lambda+epsilon+eta)/(1-delta), as in AR. The first exploratory even run (`at_exterior.py`) used a coarse factor norm and a too-large proposal threshold, and failed at beta>1; it is not evidence of a negative original direction. Revised runs use exact factor norm bounds and a proposal threshold adapted to eta. Both pass.

| Parity | Deficits | Surplus | Absorbed | Constraints incl. 512 | beta diagnostic |
|---|---:|---:|---:|---:|---:|
| even | 4628 | 12560 | 3862 | 1278 | .9959768876382481 |
| odd | 4691 | 12751 | 3916 | 1287 | .9960663693490266 |

The exact beta fractions are

- even 43377482670791775696287905481999830282157/43552700076858657734259170325391726346240
- odd 22211234390202563723990668249306161791270637/22298950224289358403613041609365394286444544

`at_check.py` checks both exact fractions <997/1000, inherited array hashes and counts from the executed logs. Conditional n_-(q0;D^0)<=2565, a reduction of 190 from AS, and q0[f]>=3/2000 <E*DE f,f> >=3/2000 ||f||² on the certified complement. The improvement is a sufficient codimension bound, not a census of negative original eigenvalues. Initial 1024 low-frequency constraints remain conservatively removed.

Reproduce from the bundle, Python 3 with NumPy/SciPy: `sha256sum -c array_sha256.txt`; optionally regenerate arrays by `python aq_scan.py` then `python as_surplus.py` (the latter imports the AQ prelude); `OPENBLAS_NUM_THREADS=4 python at_exterior_v2.py 0 > at_even.log`; `OPENBLAS_NUM_THREADS=4 python at_exterior_v2.py 1 > at_odd.log`; `python at_check.py`. Numerical proposals may vary and must be rechecked. The log checker records the executed result, not arbitrary reruns. No remote CI or Lean.

## 4. Next gate — AU: use the actual exterior Gram

The current Bessel subtraction uses the identity as an upper bound for a smaller exterior metric. Define V=F Phi_J diag(sqrt(b_j)w_j), where Phi_J synthesizes the deficit circle modes. Then V*Q V = W^(1/2)(diag(b)-M-K)W^(1/2)=H_ext, with 0<=H_ext<=I, and V*Q E T=W^(1/2)K=L. For any exterior dual trial vector v, the constrained inverse is bounded above by the unconstrained energy minus 2 Re<v,Lz>+<H_ext v,v> (equivalently subtract the gain 2 Re<v,Lz>-<H_ext v,v>). This follows from the variational distance to the physical range and Cauchy–Schwarz, or the exact exterior Schur formula. Taking v=Lz recovers at least L*L. Optimizing further, or using the positive polynomial lower inverse 2I-H_ext on the range, gives the stronger upper comparator M-K-L*(2I-H_ext)L. For singular H_ext the polynomial variational trial remains valid; no unsupported inverse is needed. This is an analytical follow-up, not yet a numerical inertia certificate. Retain other exact source deficits or evaluate the blanket low frequencies if the incremental exterior correction stalls. Full fixed-window positivity, all-support positivity, Suzuki–Connes transport, Phase IV and RH remain OPEN.