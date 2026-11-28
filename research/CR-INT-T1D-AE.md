# CR-INT-T1D-AE — Neumann-resolvent bracketing and positive Gamma boundary restoration

**Disposition:** new DERIVED / ANALYTIC results, exact-fraction arithmetic audit locally PASS; conditional on Suzuki's original pole-null Fourier identity, the standard digamma expansion and (for the high-sector bound) T1D-Q's joint arithmetic norm. No Lean/CI, independent review, global spectral sign, negative original witness, completed Suzuki–Connes transfer, Phase IV or RH. This is an alternative to AD, not an identification of its prolate eigenfunctions with Neumann modes.

## 1. Original source and exact Gamma comparison

Retain I=(-2,2), H=L²(I,dx), D⁰=C_c^∞(I) with the two Mellin moments against e^{±x/2} equal to zero. Use the original unitary Fourier identity from Suzuki arXiv:2606.09096v4 §§1.1,2.4, as in R/AA/AD: q₀=Γ-A, with mΓ(ξ)=Re ψ(1/4+iξ/2)-log π and A=sum over all 24 true prime powers n≤54 of Λ(n)/sqrt(n) times the two compressed shifts ±log n. Here q₀ is the pole-free Fourier pairing on the full interval form space; it equals the original Weil form only on the pole-null subspace. Do not discard either original pole contribution on other vectors.

Use the ordinary Neumann Laplacian L_N on I, orthonormal basis e_n(x)=c_n cos(ω_n(x+2)), ω_n=nπ/4, c_0=1/2 and c_n=1/sqrt(2) for n≥1. Its parity under reflection is (-1)^n. Set a_j=2j+1/2. The real digamma series gives

mΓ(ξ)=mΓ(0)+Σ_{j≥0} (2/a_j) ξ²/(a_j²+ξ²).

Let R_j be the compressed whole-line resolvent, kernel exp(-a_j|x-y|)/(2a_j), and N_j=(a_j²+L_N)^(-1). R_j≤N_j: the compressed resolvent solves the Robin variational problem with the positive boundary form a_j(|u(-2)|²+|u(2)|²), whereas N_j has zero boundary penalty. Equivalently, restriction of the minimizing whole-line H¹ extension produces the Robin form. Thus the complete Gamma form satisfies Γ≥D=diag(mΓ(ω_n)), via monotone sums of positive summands after subtracting mΓ(0). This comparison preserves the physical metric and leaves -A unchanged. It is not an equality or a periodic replacement of the original Gamma multiplier.

## 2. Explicit positive boundary correction, including both parities

For σ=±1 put b_{j,σ}(x)=exp(-a_j(2-x))+σ exp(-a_j(2+x)). On a parity-σ f the difference (N_j-R_j)f solves (-∂²+a_j²)w=0, and w'(2)=a_j(R_j f)(2), since the free resolvent has derivative -a_j(R_j f)(2). Also b'_{j,σ}(2)=a_j(1-σ exp(-4a_j)), and (R_j f)(2)=<f,b_{j,σ}>/(4a_j) with first-argument linearity. Therefore

2a_j(N_j-R_j) on parity σ = b_{j,σ}⊗b_{j,σ}/[2(1-σ exp(-4a_j))].

The identity on the full interval form domain is Γ=D+B, where B is the positive form sum of these rank-one terms (two per j). The j=0 terms vanish *exactly* on D⁰, by the two Mellin moments; they are not omitted off that subspace. No boundedness of the infinite correction B on L² is asserted. Zero extensions of H¹(I) functions are in H^s(R), 0<s<1/2, so all their logarithmic form integrals converge; the series identity passes to these functions by monotone convergence. Compact smooth approximation and exact moment correction follow AC/AD, or use a bounded H^s extension/density argument and two fixed smooth functions whose moment matrix is invertible.

For equal parities σ=(-1)^m=(-1)^n, the exact j-th correction matrix entry is

B_j[m,n]=2a_j² c_m c_n (1-σ exp(-4a_j))/[(a_j²+ω_m²)(a_j²+ω_n²)].

Opposite parities vanish. Derivation: the moment against the right exponential is c_n a_j((−1)^n−exp(−4a_j))/(a_j²+ω_n²); the parity-matched b moment is twice this. Consequently the exact Gamma matrix is a known diagonal plus an explicit positive rank-two series. This uses all Gamma summands, not a finite resolvent replacement. On the first N Neumann modes the discarded correction beyond j≥K, K≥1, has operator norm ≤N/[2(K−3/4)]: each parity outer product has norm² bounded by 2N/a_j² using c_n²≤1/2 and 1−σ exp(−4a_j)≤2; sum by the integral test. The bound is conservative but effective and two-sided on the finite low block.

## 3. An 82-fold smaller conditional negative-index bound

The monotone digamma series gives mΓ(ω_n)≥mΓ(ω_32768) for n≥32768. For t_j=j+1/4, x=3*32768/8, π>3 and mΓ(0)>−43/8 (R), bound the first 32 terms from below at x. The positive decreasing tail sum is ≥ its integral from 32+1/4, equal to (1/2)log(1+(x/(32+1/4))²)>log(x/(32+1/4)). The exact-rational audit proves

−43/8 + Σ_{j=0}^{31} x²/[t_j(t_j²+x²)] + log(x/(32+1/4)) > 203/25+1/8.

The lower expression is approximately 8.253244659255. Import ||A||<203/25 from T1D-Q. Then q₀[f]> (1/8)||f||² on D⁰ orthogonal to e_0,...,e_32767, since Γ≥D. This is a true infinite-dimensional conditioned sector; it does not require the T1D-Y selected six-space or AA/AD penalties. Any strictly negative subspace injects into the first 32768 coefficient coordinates, proving n_-(q₀;D⁰)≤32768, and ≤16384 per parity. These are upper bounds, not measured indices. Conditional on Q, the reduction from AB's 2693201 is >82-fold (the hypotheses are different). No positivity of the low block or mixed Schur conclusion follows from this count alone.

## 4. Exact sign diagnostic on the same pole-null trial

In y=x+2∈(0,4), f(y)=cos(w y)−c, w=π/2, c=1/(1+π²). Reflection even. The e^{±x/2} moments vanish exactly by integrating exp(y/2)cos(wy); the moment ratio to the constant is (1/4)/(1/4+w²)=c. The physical mass is 2+4c². For s=log n<4,

∫_s^4 f(y)f(y−s)dy=(4−s)cos(ws)/2−sin(ws)/(2w)+2c sin(ws)/w+c²(4−s).

The signed arithmetic contribution is −2Σ κ_n times this correlation, all 24 genuine weights. The Neumann comparator energy is q_N[f]=2mΓ(w)+4c²mΓ(0)−<Af,f>. The standard-library audit certifies q_N[f]/||f||² ∈ (−0.13300,−0.13291), in particular <−13/100. This is a rigorous counterexample to *positivity of the lower surrogate*, NOT to Weil positivity.

The actual original Gamma boundary term for this f equals Σ_{j≥1}2a_j²(1−exp(−4a_j))[1/(a_j²+w²)−c/a_j²]²; j=0 is exactly zero. For j≥1, exp(−4a_j)≤exp(−10)<1/10000. Enclose j=1,...,63 by exact rational intervals; for j≥64 the positive difference in brackets is at most 1/a_j², so the positive discarded energy is ≤1/[2(64−3/4)]. The resulting genuine original Rayleigh quotient is enclosed by (0.02534068,0.02931959), hence strictly between 1/40 and 3/100. The trial belongs to the pole-null logarithmic form completion, not literally C_c^∞. Density with exact pole-moment corrections supplies genuine smooth tests with positive energy near it. No specific compact smooth negative trial is claimed. The result exposes a false-negative mechanism if the positive Gamma boundary terms are discarded. Floating-point exploratory eigenvalues of Neumann-only compressions were negative in both parities, but no other numerical spectral signs are certified here.

## 5. Boundary-restored finite Schur gate with exact pole constraints

Choose N≥32768, K finite ≥1. Let P project onto e_0,...,e_{N−1}, Q=1−P, B^(K)=Σ_{j<K}B_j, C_K=B^(K)−A, and F_K=D+C_K. Keep j=0 even off the constrained space; on pole-null vectors its quadratic contribution vanishes. F_K≤q₀ in the form sense. Write H=P F_K P, X=Q C_K P; no diagonal Gamma term mixes the blocks. On QH the high block satisfies F_K≥D_h:=diag(mΓ(ω_n)−203/25)_{n≥N}>0. Hence for x∈PH, y∈QH, q₀[x+y]≥<Hx,x>+2Re< Xx,y>+<D_h y,y>.

Let C:H→C² give the two original moments (one cosh, one sinh), C_l=CP, C_h=CQ. Enforce C_l x+C_h y=0 exactly. Completing the square in y with D_h gives the finite sufficient matrix test

H−X*D_h^(−1)X +(C_l−C_h D_h^(−1)X)* [C_h D_h^(−1)C_h*]^(−1) (C_l−C_h D_h^(−1)X) ≥ 0.

The 2×2 middle matrix is strictly positive: each parity has nonzero tail moment coefficients, with magnitude proportional to c_n/(1/4+ω_n²). Parity decouples it into two scalars. This test retains the physical constraints without a penalty parameter, and uses a mode-dependent high inverse rather than AD's scalar d_*^(-1). It is sufficient, not necessary; a negative eigenvalue of this lower matrix is not a negative original test. Positivity of F_K in the constrained Schur test would suffice for original pole-null positivity. The form F_K itself increases monotonically with K; the displayed sufficient matrix need not increase when K changes because its high-block lower bound is held fixed while the mixed block changes. A tighter high-block lower bound can also retain finite positive boundary terms. No numerical Loewner certification is given.

All matrix entries are explicit: A[m,n] is a sum of 24 cosine shift overlaps and their transposes, integrated over (log n,4) using elementary product-to-sum; no prime collision is ignored. The moment rows are integrals against cosh/sinh and the correction rows are the rational/exponential expression in §2. Infinite mixed weighted sums are effective: for a finite low vector x, ||x'||≤ω_{N−1}||x|| and ||x||∞≤sqrt((2N−1)/4)||x||. For M>N a single compressed shift g of x is piecewise H¹, with one interior jump bounded by ||x||∞. Integration by parts against e_m (m≥M) and sine Parseval give

||Q_M g|| ≤ ||x'||/ω_M + sqrt(8)/(π sqrt(M−1)) ||x||∞,

where Q_M omits the first M modes. Sum over the original 48 weighted shifts (total weight <25 by AA); for B^(K) use ||Q_M B^(K)x||≤||(B^(K)x)'||/ω_M with derivative operator bound Σ_{j<K,σ} ||b'_{j,σ}|| ||b_{j,σ}||/[2(1−σ exp(−4a_j))]. Thus the X tail beyond M has explicit O(M^(-1/2)) norm decay for fixed N,K. Multiply its squared bound by 1/(mΓ(ω_M)−203/25) to enclose the omitted X*D_h^(-1)X; C_h tails and cross tails follow from their explicit coefficients and Cauchy–Schwarz. This gives a finite certified evaluation procedure, not a claim that its necessary precision is affordable or that the actual matrix has been evaluated.

## 6. Audit and next gate

`python3 src/audits/canonical_return_t1d_ae_neumann.py` uses exact Fraction arithmetic: Machin/atanh bounds, Euler constant and digamma partial sums, exact prime-power census, sine/cosine Taylor enclosures, high cutoff and the same-vector sign comparison including positive boundary tail. It does not mechanically prove the resolvent order, correction kernel, source identity, Q operator bound, smooth density, or constrained Schur theorem. No floating diagnostics are used in strict comparisons.

**Next AE→AF:** evaluate the parity-separated boundary-restored constrained Schur criterion, or use the exact original Gamma correction to certify a negative smooth direction. Prioritize effective mixed-tail estimates and cancellation-resistant spectral arithmetic. The conditional index is bounded by 32768, not zero. Full fixed-window pole-null positivity, full pole carriers, support-independent Weil positivity, completed Suzuki–Connes transfer, Phase IV and RH remain open. Keep PR #2 draft.