# CR-INT-T1D-AH — Scaled Gamma series, adaptive rational congruence and the weighted high inverse

**Disposition:** DERIVED / ANALYTIC with executed exact-fraction interval audit PASS, conditional on the AE/AF/AG original-source identities and Q's joint prime norm. Not Lean, remote CI, independent review, global fixed-window positivity, completed Suzuki–Connes transfer, Phase IV or RH.

## 1. Same source, larger low space

Keep the original physical L²((-2,2),dx) metric and pole-null form q₀=Γ−A, with the full digamma multiplier and all 24 genuine prime powers. The pole-free Fourier pairing equals the original Weil form on the two-Mellin-pole-null subspace only. In parity p=0,1, set w_n=nπ/4 and p_n(x)=cos(w_n(x+2)). Use n=p,p+2,...,p+34 (18 raw modes). Write f_i=p_{p+2i}−r_i p_p for i=1,...,17, where r_i=(1/4+w_p²)/(1/4+w_{p+2i}²). This is the exact parity moment correction of AF. The two resulting 17-dimensional spaces form V34, a subspace of the original pole-null logarithmic form domain. Raw-mode masses are 4 for n=0 and 2 otherwise. No substitution by periodic Weil tests is made.

## 2. Scaled evaluation of the complete Gamma tail

AF's absolute interval grid rounds very small unscaled Hurwitz-zeta tails before multiplication by large powers. On enlarging the basis this amplifies rounding drastically. The previous interval LDL failed to resolve the eighth pivot per parity; this is an enclosure failure, not a negative direction.

Take K=96, t=K+1/4, R=38. For integer s≥2 put Z_s=t^(s−1) ζ(s,t). Evaluate this normalized quantity with the Euler–Maclaurin remainder used by AF, through B_32:

    Z_s = 1/(s−1) + 1/(2t)
        + Σ_{k=1}^{16} B_{2k}(s)_{2k−1}/[(2k)! t^(2k)] + E_s,
    |E_s| ≤ |B_32|(s)_31/[(32)! t^32].

This follows by integration against the periodic Bernoulli polynomial with its standard sup bound. On the exact 10^−55 directed grid, use rho=(w/2)^2/t². The digamma tail after the first K summands is Σ_{r≥1} (−1)^(r−1) rho^r Z_{2r+1}; the absolute remainder after r=38 is ≤rho^39 Z_79. For b_sigma(w)/w² and the equal-frequency boundary entry, the nonexponential boundary tails are respectively (1/(2t))Σ_{r≥0} (−1)^r rho^r Z_{2r+2} and the same with coefficient r+1. The absolute remainder of either after r=38 is at most

    (Z_2/(2t)) rho^39 (40−39rho)/(1−rho)^2.

The omitted exponential part is at most exp(−8K−2) Z_2/(2t). The first K−1 positive boundary summands, all diagonal Gamma terms, and all original signed shift entries are retained as in AE/AF/AG. For distinct same-parity frequencies the raw off-diagonal entry (with j=0 omitted only for the subsequent exact pole-null congruence) is the divided difference of F_sigma(w)=b_sigma(w)+2w Σ κ_h sin(w log h). Opposite parity vanishes. This replaces repeated two-variable tail evaluations, not the form. All 17×17 physical pole-null energy matrices have entry interval widths <10^−47 in the executed audit.

## 3. Rational congruence rather than recursive interval elimination

Let Q be the actual energy Gram of one parity and G its physical Gram. A midpoint elimination proposes a unit lower triangular rational matrix U by exact rounded elementary row/column congruences. The proposal is not trusted as a certificate. Reevaluate C=U Q U^T directly from the original interval entries. For d_i=lower(C_ii)−Σ_{j≠i} upper(|C_ij|), the audit checks every d_i>0 and U unit triangular. Hence Q is positive and C≥diag(d) by the elementary symmetric diagonal dominance bound.

For any row functional v in the original basis, its energy-dual norm is bounded by sqrt(Σ_i |(Uv)_i|²/d_i). For the physical norm relative to energy, use the trace upper bound tr(diag(d)^−1 U G U^T), computed from the exact raw cosine masses and interval moment coefficients. The audit proves this trace <10^45 in each parity. Therefore q₀[f]>10^−45||f||² on V34. Printed upper-bound reciprocals are 3.9852×10^−45 even and 7.2566×10^−44 odd. These are conservative certified lower margins, not exact minimum eigenvalues. For the geometric tail, bound the raw coefficient ℓ¹ norm divided by sqrt(q₀[f]) by sqrt(Σ_i (Σ_m |c_im|)^2/d_i), with c_im the raw coefficients of the transformed basis. This does not require an absolute error tolerance smaller than the smallest Rayleigh quotient in the untransformed basis.

## 4. Mode-dependent high inverse and new cutoff

Take N0=27400. On the genuine pole-null high Neumann form sector Y, Q's ||A||<203/25 and AE's Γ≥D imply q₀[y]≥<D_h y,y>, with D_h(n)=mΓ(w_n)−203/25. Use the monotonicity of the digamma series. For an even integer n=a, put x=a π_lower/8 and t0=32+1/4. A certified lower bound on D_h(a) is

    mΓ(0)_lower + Σ_{j=0}^{31} x²/[(j+1/4)((j+1/4)^2+x²)]
                  + log(x/t0)_lower − 203/25.

The last term uses the decreasing integral test for the positive digamma tail. Round down to a 10^−6 rational grid. Choose even band starts 27400,27450,...,27950; 28000,28250,...,35750; then 36000,40000,48000,64000,96000,128000,256000,512000,1024000, with the last band extending to infinity. The 53 resulting lower values d(a) are positive and nondecreasing, and d(27400)=3268/10^6>3/1000. Thus the staircase lies below the actual diagonal on every band. It is not a replacement of the original signed form.

For either parity and p≥1, the exact decreasing integral test at even band endpoints gives

    Σ_{n in [a,b), matching parity} n^(−2p) ≤ [(a−2)^(1−2p)−(b−2)^(1−2p)]/[2(2p−1)],

where the second term is zero for the last infinite band. Divide by d(a), sum all bands and use w_n≥(785/1000)n to enclose the weighted inverse tails. The resulting weights 1/[n² d(a(n))] are nonincreasing, so Abel summation gives the original prime-phase cross bound

    |Σ_{n≥N0, matching parity} cos(nπ s/4)/[n² d(a(n))]|
      ≤ 1/[N0² d(N0) |sin(π s/4)|]

for nonzero phases s=log h±log k. The zero difference terms use the weighted integral sum directly. All nonzero denominators are certified by AF's interval trigonometry. The complete 24×24 ordered phase-pair sum and b_sigma(w)/w<8/5 give the leading weighted symbol norm upper diagnostic 0.08196542600470227. No independent/random phases assumption.

## 5. Relative mixed tail on all 34 low dimensions

AG's exact geometric expansion of the complete divided difference gives q₀(f,p_n)=Σ_{k≥0} [F_sigma(w_n) S_k(f)−T_k(f)]/w_n^(2k+2), with S_k=Σ_m c_m w_m^(2k), T_k=Σ_m c_m F_sigma(w_m) w_m^(2k). Compute the energy-dual bounds in the new congruence for k=0,...,7. For k≥8, with L a strict integer upper bound on the low frequencies and C the certified coefficient ℓ¹/energy bound, the raw mixed remainder is bounded by

    2 C L^16 sqrt(q₀[f]) (30/w_n^17 + 30 L/w_n^18).

Here |F_sigma(w)|≤30w from the full Gamma integral estimate and the genuine prime weight sum, and L/w_n<1/2. Use the weighted inverse tail sums in §4 and triangle inequality in weighted ℓ². Divide by sqrt(2) for physical normalized high cosines. The exact interval audit proves

    ||D_h^(−1/2) Xf||² < (16/25) q₀[f]

in each parity and hence on V34. The computed upper-bound diagnostics are 0.6279930448514953 and 0.6282890948477383. The high cutoff is lower and the low space larger than AG; the relative bound is deliberately weaker than AG's <11/50. The previous 28000 exploratory cutoff also passed, but the final executable certificate fixes 27400.

## 6. Conditional positivity on a larger sector, not whole positivity

For f∈V34, y∈Y, weighted Cauchy–Schwarz and q₀[y]≥<D_h y,y> imply

    q₀[f+y] > (1/5)(q₀[f]+<D_h y,y>)
             > (10^−45/5)||f||² + (3/5000)||y||²

for f+y nonzero (the inequalities follow also if one summand is zero). This is a positive form sector of the *same original* source. For smooth pole-null z whose projection onto the first 27400 Neumann modes lies in V34, set f=P_low z, y=z−f. Both are pole-null form-domain vectors, so the result applies without requiring the projections to be compact smooth. The codimension of these smooth-source constraints is at most 27400−34=27366, or 13683 in each parity. Therefore the original fixed-window pole-null smooth negative index is conditionally ≤27366. This is an upper bound, not an actual count. The intermediate modes and all their mixed interactions remain uncontrolled. No negative original witness was established.

## 7. Reproduction and next gate

Run `python3 canonical_return_t1d_ah_weighted_tail.py` with `canonical_return_t1d_ah_scaled_matrix.py` and the unchanged AF module beside it. The included audit_output.txt records the local execution. The strict checks use exact fractions, not the printed floating diagnostics. Source identities, analytic Gamma/tail formulas, Q's joint prime bound, interval arithmetic and spectral arguments remain subject to independent review. No Lean/remote CI claim.

**Next AI:** address the remaining intermediate sector by a source-preserving block elimination using the divided-difference structure, exact pole constraints, and a genuinely improved high inverse (not just another absolute-error rank estimate). In particular, the first 34 directions are not a model for the spectral sign of all intermediate directions; certify the Schur complement on that sector or find a genuine original smooth negative direction. Completed Suzuki–Connes transfer, Phase IV, all-support positivity and RH remain open.