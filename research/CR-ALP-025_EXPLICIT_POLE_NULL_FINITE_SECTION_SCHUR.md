# CR-ALP-025 — Explicit physical pole-null sine section and positive infinite Weil tail

**2026-11-23; research-only; DRAFT/NONPROMOTING.** This is an analytical finite-section criterion, not an evaluated full-collar positivity certificate, a Lean kernel result, an F41 H/V transfer, or RH.

## 0. Original source contract

The unchanged original Zeta23.EF.literatureRHS is pinned in Weil-Lab, research/iv068-original-source-kernel-gate, WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean at Git blob 0034976fe2f1e5fa9d411108736160542eb94293; the original IV011 cutoff consumer is at blob 14d8f01ce34f764900eccc68399d7323b4f4a25c. The form has both poles, **negative paired** von Mangoldt samples at {2,3,4,5,7,8,9}, and genuine gammaBracket(r)=Re digamma(1/4+ir/2)-log(pi). In particular no positive surrogate is substituted.

CR-ALP-024 established analytically, not in Lean, that for the original two-pole-null sources on I=(-a,a), a=11/10,

\[
W(f,f)\ge E_{\log}(f)-16\|f\|_2^2,\qquad
E_{\log}(f)=\frac1{2\pi}\int_{\mathbb R}\log(1+|r|)|\widehat f(r)|^2\,dr.
\tag{25.1}
\]

The original signed form extends by closure to the logarithmic Fourier form domain. All claims below assume this previously stated analytic inequality.

## 1. An explicit physical, not F41, projected basis

Set b_+(x)=exp(x/2), b_-(x)=exp(-x/2) on I and F_±(f)=integral_I f(x)b_±(x) dx. Their genuine physical L² Gram is

\[
G=\begin{pmatrix}2\sinh a&2a\\2a&2\sinh a\end{pmatrix},
\qquad \det G=4(\sinh^2a-a^2)>0.
\]

Write \(\mathscr H_0=\ker F_+\cap\ker F_-\), and let \(\Pi_0=I-BG^{-1}B^*\) be its physical L² orthogonal projector, with B sending the two coordinate vectors to b_+,b_-.

Define, for j≥1,

\[
e_j(x)=a^{-1/2}\sin\frac{j\pi(x+a)}{2a},\quad
\nu_j=j\pi/(2a),\quad v_j=\Pi_0e_j.
\]

Every v_j is exactly two-pole-null. Every finite prefix of the v_j is linearly independent: a finite sine polynomial vanishes at the two endpoints, whereas a linear combination of b_+,b_- vanishing at both endpoints is zero. The **physical mass matrix** is

\[
M_N=(\langle v_i,v_j\rangle)_{i,j\le N}
=I_N-C_N^*G^{-1}C_N\succ0,\quad
(C_N)_{\sigma j}=F_\sigma(e_j).
\tag{25.2}
\]

This is an exact **rank-two physical Gram downdate**, not an orthonormality assumption. The moments are elementary:

\[
F_\sigma(e_j)=\frac{\nu_j}{\sqrt a}
\frac{e^{-\sigma a/2}-(-1)^j e^{\sigma a/2}}{\nu_j^2+1/4},\quad \sigma=\pm1.
\]

Let E_N=span(v_1,...,v_N), P_N its physical orthogonal projector in H_0, and Q_N=I-P_N. If q belongs to Q_N H_0, then \(\langle q,e_j\rangle=\langle q,\Pi_0e_j\rangle=0\) for all j≤N. Thus the projected tail has zero first N ordinary sine coefficients.

The zero-extended b_± have O(1/|r|) Fourier tails and belong to the logarithmic form closure of compact smooth interior functions (interior-cutoff approximations converge in H^s for 0<s<1/2); therefore Π_0 preserves the log form domain. Sine polynomials approximate compact interior smooth functions in H_0^1, hence in the weaker log form norm, and their Π_0 images are form dense in the pole-null domain.

## 2. Positive lower bound for the ENTIRE infinite-dimensional complement

Under the paper Fourier convention, direct integration gives

\[
\widehat e_j(r)=\frac{e^{-iar}}{\sqrt a}
\frac{\nu_j(1-(-1)^je^{2iar})}{\nu_j^2-r^2}.
\]

For |r|≤R and ν_j≥2R, its absolute value is ≤8/(3√a ν_j). Since every q∈Q_N H_0 has its first N sine coefficients zero, Cauchy–Schwarz and Σ_{j>N}j^{-2}≤1/N imply

\[
\frac1{2\pi}\int_{-R}^R|\widehat q(r)|^2\,dr
\le\delta_{N,R}\|q\|_2^2,\qquad
\delta_{N,R}=\frac{256aR}{9\pi^3N},
\quad \nu_{N+1}\ge2R.
\tag{25.3}
\]

Choose the **explicit** R_0=10^14 and N_0=250,000,000,000,000. With π>3 and a=11/10, ν_{N_0+1}>2R_0 and

\[
\delta_{N_0,R_0}<\frac{2816}{6075}<\frac12.
\]

CR-ALP-024's independently exact-audited log(1+R_0)>32 and its signed lower bound (25.1) then yield, for **every** q in Q_{N_0} H_0 intersected with the form domain,

\[
\boxed{W(q,q)\ge \mu\|q\|_2^2,\qquad
\mu=\frac{7088}{6075}>1.}
\tag{25.4}
\]

This is an explicit **strictly positive infinite-dimensional tail theorem**, not an entire collar PSD theorem. The N_0 prefix is enormous. In particular N_0 does not identify, or validate sufficiency of, the historical 242 F41 selected physical generators.

## 3. Original signed finite-matrix / Schur certificate

On H_0 the pole terms vanish. The compressed original form operator is A=Π_0(Γ−B_primes)Π_0, where Γ has multiplier Re digamma(1/4+ir/2)-log π and B_primes is precisely the sum of seven positive weights times both physical shifts. Each fixed v_j is in Dom(A): its Fourier transform is O(1/r), and |gammaBracket(r)|≤log(1+|r|)+8, so the multiplier image is L².

Form the finite **original signed Weil** and residual matrices

\[
H_N=(W(v_i,v_j))_{i,j\le N},\quad
z_j=\Pi_0Av_j,\quad
D_N=(\langle z_i,z_j\rangle)_{i,j\le N},\quad
B_N=D_N-H_NM_N^{-1}H_N.
\tag{25.5}
\]

B_N is PSD **identically**: \(x^*B_Nx=\|Q_N A\sum_jx_jv_j\|_2^2\). All norms, projections, and the Gram matrix are the actual physical L² objects, not frozen F41 H/V.

Writing f=p+q, p=Σx_jv_j, q∈Q_N H_0, the complete original signed form obeys

\[
W(f,f)\ge x^*H_Nx+
2\operatorname{Re}\langle Q_NAp,q\rangle+
\mu\|q\|_2^2
\ge x^*(H_N-\mu^{-1}B_N)x.
\]

Consequently a fully certified finite matrix inequality

\[
\boxed{H_{N_0}-\frac{6075}{7088}B_{N_0}\succeq0}
\tag{25.6}
\]

would prove nonnegativity of the **entire two-pole-null original fixed collar**. It has **NOT** been evaluated. The use of μ^{-1}I for the unknown tail inverse gives a *sufficient*, not a necessary, test. A certified negative vector of H_N itself would be an actual original signed-Weil negative source; a failure of the sufficient inequality is otherwise inconclusive.

## 4. Effective finite input enclosure and remaining difficulty

Let \(K_j=|v_j(-a^+)|+|v_j(a^-)|+\int_I|v_j'(x)|dx\). Integration by parts gives |\hat v_j(r)|≤K_j/|r|. If Γ_T cuts the genuine Gamma multiplier to |r|≤T and T≥1, then

\[
\boxed{\|\Gamma v_j-\Gamma_Tv_j\|_{L^2(I)}^2
\le\frac{K_j^2}{\pi T}
[(\log T+9)^2+2(\log T+9)+2].}
\tag{25.7}
\]

This provides an explicit and computable error modulus for every fixed matrix column after the compact-frequency digamma integrals are enclosed. The prime shifts have finite elementary physical overlap integrals, and Π_0 is norm one. Certified digamma interval arithmetic and quadrature would therefore enclose H_N and B_N, **in principle**. No real N_0-dimensional matrix, or its intervals, has been computed or certified.

Reflection splits the original full form analytically into even/odd blocks. The odd-index sine modes are even about zero, the even-index modes odd; the rank-two pole projection becomes rank one in each parity block. The receiving Lean parity gate has **not** been compiled.

## 5. Audit and governance

The exact stdlib Fraction audit, src/audits/cr_alp_025_finite_section_tail_audit.py, has been executed locally: **PASS** for the cutoff, rational leakage cap, positive μ, rank-two projection/Gram toy and both Schur success/failure examples. It does **not** mechanically prove the continuum Fourier, digamma, trace, or form-domain claims.

**CR-ALP-025 status:** explicit strictly positive infinite tail **HUMAN ANALYTIC RESULT**; full finite positive certificate **OPEN**; Lean kernel **OPEN**; F41 242 sufficiency **OPEN**; entire collar, arbitrary supports, and RH **OPEN**. Canonical-Return draft PR #3 stays **nonpromoting**. No Weil-Lab, Formalization, Hodge, F41/F42 or RH acceptance change.

**Next CR-ALP-026:** sharpen the Gamma/prime constant and source-specific sine-band concentration bounds; aim to reduce N dramatically, before attempting the actual H_N/B_N interval matrix gate.
