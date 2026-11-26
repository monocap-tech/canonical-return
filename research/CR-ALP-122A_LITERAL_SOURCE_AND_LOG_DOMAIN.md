# CR-ALP-122A — Literal source and logarithmic-domain bridge

Research-only; classical analytic derivation, not native Lean certification. The exact CR-ALP-120 rank result is conditional. This note discharges the abstract domain/decomposition requirements of 122 for the historical source convention; numerical source attachment and tail-form-order replay remain separate.

## 1. Literal source readback

Reopened monocap-tech/weil-lab at research/iv068-original-source-kernel-gate. Git blobs match the historical pins:
- Def_Zeta23_ExplicitFormula.lean: 0034976fe2f1e5fa9d411108736160542eb94293.
- ActualZetaPhaseIVOriginalCrossSupportExtension.lean: 5ef51fd37337c466620eafccee2a27b77955746e.
- ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean: 14d8f01ce34f764900eccc68399d7323b4f4a25c.
- Def_Zeta23_Defs.lean: c5ebc771817c56114b552ec5ab9c300bb6ec02f1 (paperFT).

The actual IV011 theorem neutralIV011MixedCorrelation_eq_schwartzIntegral identifies the mixed correlation as k_{s,t}(u)=int t(x) conj(s(x-u)) dx, first-antilinear and second-linear. Therefore by Fubini paperFT(k)(r)=conj(F_s(r)) F_t(r) on the real line, with F_f(r)=int f(x) exp(irx) dx. At complex z it equals F_t(z) conj(F_s(conj z)). For the diagonal, the two poles sum to 2 Re(F_f(i/2) conj(F_f(-i/2))). The original literatureRHS uses the paired samples with the negative prime sign and (2 pi)^-1 times gammaBracket=Re digamma(1/4+ir/2)-log pi. IV011OriginalMixed_finitePrimeArithmetic gives the exact cutoff, not a positivity estimate. At a=11/10 the CR-ALP-024 audited log9<2a<log10 gives precisely n=2,3,4,5,7,8,9 with the von Mangoldt weights (zero for n=0,1 and 6). The prime part is a finite sum of compressed signed translation operators. No alternative physical H/V object is substituted. Source identity here concerns the pinned historical definitions, not every current IV receiver.

## 2. Domain of the extended form

Let I=(-a,a), H=L2(I,C) with zero extension J to R, and D={f in H: int log(1+|r|)|F_f(r)|^2 dr<infinity}. Give D the norm from weight w=1+log(1+|r|) with measure dr/(2pi). The Fourier transform is the paper convention. D is complete: the weighted space is complete, and the support condition is closed in L2. The CR-ALP-024 two-sided bound on g implies |g(r)|<=8+log(1+|r|); hence the archimedean sesquilinear form is continuous on D by weighted Cauchy-Schwarz. Each translation and each pole functional is bounded on H. Thus the exact smooth original form extends continuously to D, with the same expression. This is a continuous form extension, not an assertion that nonsmooth sources satisfy EF_lit on zeros.

Reflection Rf(x)=f(-x) is unitary on H and D; the multiplier g is real and even (classical digamma conjugation). The paired translations commute with reflection, as does the symmetric pole expression. Consequently opposite reflection parities have zero mixed pairing. Complexification is legitimate; the signed operator is real under conjugation, so the real symmetric sector matrix extends Hermitianly, and a positive real symmetric matrix gives positive complex energy by splitting real and imaginary parts.

## 3. Core density with exact pole constraints

For f in D, set f_s(x)=s^(-1/2) f(x/s), 1/2<=s<1, with zero extension. It is supported in [-sa,sa]. In frequency F_{f_s}(r)=s^(1/2)F_f(sr). These dilation operators are uniformly bounded on weighted L2, since w(u/s)<=w(u)+log2<=(1+log2) w(u). On continuous compactly supported frequency functions, strong convergence as s approaches 1 follows from dominated convergence; density and uniform boundedness extend it to all weighted L2. For fixed s, convolution with a smooth compactly supported approximate identity of radius epsilon<(1-s)a preserves interior support and converges in D by dominated convergence of the uniformly bounded Fourier multipliers. A diagonal sequence shows C_c^infty(I) dense in D, including the boundary-jump functions used below.

Set H0={f:F_+(f)=F_-(f)=0}, F_pm=int f exp(+-x/2), and D0=D intersect H0. The constraints are continuous on D. To correct the smooth approximants exactly, fix a nonnegative, nonzero even smooth bump chi supported in I and nonzero away from 0. The two functions chi cosh(x/2) and chi sinh(x/2) have diagonal moment matrix for the equivalent cosh/sinh constraints, with strictly positive diagonal entries. Subtract their linear combination with coefficients equal to the moment errors divided by these entries. The corrections tend to zero in D and yield exact constraints. Thus C_c^infty(I) intersect H0 is dense in D0. The same construction preserves a chosen reflection parity. This is a classical proof, not a reuse claim about a compiled density theorem.

## 4. E85/Q85 is an exact form-domain decomposition

Use the CR-ALP-046 projected sines v_j=e_j-c_j b/B, with b=cosh(x/2) for odd sine index and sinh(x/2) for even, B=||b||^2, c_j=<b,e_j> (real). Zero extension of each sine or complete pole function is piecewise C1 with finite endpoint jumps. Integration by parts on I gives |F_h(r)|<=C_h/|r| for |r|>=1 and |F_h(r)|<=||h||_1 otherwise. Consequently h belongs to D and even g F_h is in ordinary L2 because int_1^infinity log^2(1+r)/r^2 dr<infinity. This applies to every v_j, to their finite span E, and to all finite sine plus full pole frozen corrections of the recorded type.

The projected sines are independent: if sum_j alpha_j v_j=0 then a finite sine sum equals a scalar times b in I. Both sides are continuous on the closed interval; the sine sum vanishes at both endpoints, while b(a) is nonzero in each parity. Thus the scalar vanishes, and ordinary sine orthonormality forces all alpha_j=0. Dimensions are 43/42 for indices 1..85. The physical Gram is delta_ij-c_i c_j/B, positive definite. For f in H0, <v_j,f>=<e_j,f> because <b,f>=0. Thus the H-orthogonal complement of E inside H0 consists exactly of the functions with vanishing first 85 sine coefficients and both pole constraints, the Q85 of the historical receiver. Denote the physical orthogonal projections by P_E and Q=I-P_E on H0. Since E is finite dimensional and contained in D0, P_E:H0->D0 is bounded (finite sum of continuous H coefficients), and Q preserves D0. Hence D0=E direct-sum (D0 intersect QH0). No identification with the domain of a different spectral operator is needed.

The smooth core in section 3 also yields a dense tail core in form norm after applying Q, although these projected approximants need not themselves be smooth. No assumption that P_E preserves compactly supported smoothness is used.

## 5. Actual mixed coupling is physically L2-bounded

For p in E, define Gp=J^* paperFourierInverse(g F_p), using inverse (2pi)^-1 and the usual L2 Plancherel extension. Section 4 proves this is well-defined in H. Let P denote the actual bounded signed prime translation sum on I. On H0 the pole part vanishes. Define R p=Q (Gp+Pp). For q in D0 intersect QH0, Plancherel gives W(p,q)=<R p,q> (first-antilinear convention), with R:E->QH0 bounded because E is finite dimensional. Hence the true block decomposition W(p+q)=W(p)+2 Re<W(p,q)>+W(q) uses a bounded physical mixed functional. Here <W(p,q)> means the scalar W(p,q), not another inner product. This establishes the existence and exact analytic meaning of R required by 122, not the numerical identity of every archived coordinate of R.

## 6. Remaining gate

Under the historical Gamma growth premise, sections 2–5 prove domain inclusion, exact finite/tail decomposition, constrained smooth density and bounded mixed coupling. Literal source readback confirms the signs, normalization and input orientation. Still required: verify the archived energy-conditioned coordinate map and correction X against this R, all physical interval premises and the common K/B, and the full infinite-dimensional Gamma-tail comparison W_Q>=C. The source radius remains 11/10, and the pole-null restriction is unchanged. No fixed-collar positivity promotion, native Lean proof, unrestricted-support extension or RH claim is made.
