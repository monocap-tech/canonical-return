# CR-ALP-028 — Exact physical Gamma kernel, rank-two boundary basis, source-signed prime intervals, and a matrix-free Schur gate

**2026-11-23 PST — Canonical-Return research only; DRAFT / NONPROMOTING.** This stage advances the **original, complete signed Weil fixed-collar** calculation without evaluating or certifying the 15,000-dimensional low block. It proves mathematical identities from the classical digamma integral and explicit physical sine calculus; executes exact outward-rounded rational checks of selected ORIGINAL prime entries and conservative Fourier residual budgets; and runs a separate numerical (not certified) complete Weil 6×6 diagnostic. **No Lean kernel compilation, global fixed-collar positivity, F41 completeness, support extension or RH.**

## 0. Native source and retained tail

Pinned unchanged Weil-Lab source: `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean` blob `0034976fe2f1e5fa9d411108736160542eb94293` (`Zeta23.EF.literatureRHS`), and native paired prime cutoff `ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean` blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`.

On I=(-a,a), a=11/10, retain the original **two poles minus seven paired von Mangoldt translations plus full digamma/Gamma**. The seven n's are {2,3,4,5,7,8,9}, with w_n=Lambda(n)/sqrt(n)=log(p)/sqrt(p^k). On the two-pole-null subspace H0, the residues vanish *exactly*. CR-ALP-027 analytically gives W(q,q)>=1/4 ||q||² for the ENTIRE orthogonal complement of E_15000=span(Pi0 e_1,...,Pi0 e_15000), in the logarithmic Fourier form domain. That continuum theorem was not Lean-compiled.

## 1. A better coordinate system for the SAME exact E_N

Put e_j(x)=a^(-1/2) sin(nu_j(x+a)) and nu_j=j*pi/(2a). Let b_1=cosh(x/2), b_2=sinh(x/2). Pole-null is equivalent to being orthogonal in physical L² to both b_1,b_2. Define anchor k(j)=1 when j is odd, 2 when j is even; and for 3<=j<=N set

\[
\lambda_j=\frac{\nu_j(\nu_{k(j)}^2+1/4)}{\nu_{k(j)}(\nu_j^2+1/4)},\qquad
w_j=e_j-\lambda_j e_{k(j)}.
\tag{28.1}
\]

The two pole moments cancel separately and EXACTLY, and w_j vanishes at both endpoints. Set u_1=Pi0 e_1 and u_2=Pi0 e_2. Then v_j=Pi0 e_j=w_j+lambda_j u_{k(j)}, hence **without changing the low subspace**

\[
\boxed{E_N=\operatorname{span}\{w_3,\ldots,w_N,u_1,u_2\}.}\tag{28.2}
\]

This is N−2 zero-boundary H¹ functions plus just **two boundary-jump directions**. The exceptional anchors are explicit:

\[
\begin{gathered}
u_k=k\pi/(2a),\quad
c_1=\frac{2\nu_1\cosh(a/2)}{\sqrt a(\nu_1^2+1/4)},\quad
c_2=-\frac{2\nu_2\sinh(a/2)}{\sqrt a(\nu_2^2+1/4)},\\
B_1=a+\sinh a,\quad B_2=\sinh a-a,\quad
\gamma_k=c_k/B_k,\quad u_k=e_k-\gamma_k b_k,\quad\eta_k=c_k^2/B_k.
\end{gathered}\tag{28.3}
\]

For either parity, ordered as its w-columns followed by its anchor u_k, the EXACT physical mass Gram is

\[
\boxed{M_{N,k}=\begin{pmatrix}I+\lambda\lambda^\mathsf T&-\lambda\\-\lambda^\mathsf T&1-\eta_k\end{pmatrix},\quad
\delta_{N,k}=\frac1{1+\|\lambda\|_2^2}-\eta_k>0.}\tag{28.4}
\]

Positivity of delta is a finite Parseval/Bessel fact: the omitted sine coefficients of b_k are nonzero. The well-conditioned interior Gram is a **rank-one update** of identity in each parity; the only ill-conditioned boundary correction is rank two in total. For N=15,000, floating diagnostics give delta_odd≈2.52e−5 and delta_even≈5.02e−5. These are **not certified real enclosures** and should NOT be used as an asserted numerical lower bound.

## 2. Exact original paired prime translations: no quadrature

For 0<=t<=2a, define S_t f(x)=f(x−t) with zero extension and U_ij(t)=<e_i,S_t e_j>+<e_i,S_{−t} e_j>. Direct integration of the two original sine functions gives

\[
\boxed{U_{ij}(t)=\begin{cases}
\frac{2a-t}{a}\cos(\nu_jt)+\frac{\sin(\nu_jt)}{a\nu_j},&i=j,\\
\frac{2[\nu_i\sin(\nu_jt)-\nu_j\sin(\nu_it)]}{a(\nu_i^2-\nu_j^2)},&i\ne j\ \text{and }i\equiv j\pmod2,\\
0,&i\not\equiv j\pmod2.
\end{cases}}\tag{28.5}
\]

Every original SIGNED prime matrix entry is now the exact finite expression

\[
\boxed{P_{ij}=-\sum_{n\in\{2,3,4,5,7,8,9\}}\frac{\Lambda(n)}{\sqrt n}\,U_{ij}(\log n).}\tag{28.6}
\]

For interior w_i=e_i−lambda_i e_anchor, expand (28.6) by four terms. The separate exact arithmetic audit encloses the **original signed PRIME-ONLY diagonal entries** as follows (strict certified rational bounds):

| Original source direction | Signed prime-only term, rigorous interval |
|---|---|
| w_3 | (3/5, 7/10) |
| w_4 | (3/10, 2/5) |
| w_5 | (1/100, 1/40) |
| w_6 | (−1/4, −1/10) |
| w_7 | (−2/5, −1/4) |
| w_8 | (−3/5, −2/5) |

All opposite-parity prime cross entries are EXACTLY zero. The auditor uses integer/Fraction outward rounding at 10^-38, Machin arctangent π enclosures, atanh logarithms, integer sqrt bounds and Taylor trig remainders. These values are **not the complete Weil form**; the Gamma contribution can cancel their signs. No negative RH witness is claimed.

## 3. NEW exact finite-interval physical-space identity for the COMPLETE original Gamma term

The classical digamma partial-fraction identity (α=1/4, t=r/2) implies

\[
\operatorname{Re}\psi(1/4+ir/2)-\psi(1/4)
=2\int_0^\infty\rho(s)(1-\cos(rs))\,ds,
\qquad \rho(s)=\frac{e^{-s/2}}{1-e^{-2s}}.
\tag{28.7}
\]

The complete original Gamma symbol is g(r)=Re psi(1/4+ir/2)−log(pi), so g(0)=−EulerGamma−pi/2−3log(2)−log(pi). By Plancherel and positivity/Tonelli of the difference form,

\[
\boxed{W_\Gamma(f,f)=g(0)\|f\|_2^2+
\int_0^\infty\rho(s)\|f-S_sf\|_2^2\,ds.}\tag{28.8}
\]

This identity is exact and does **not** make W_Gamma positive, because g(0)<0. If supp f ⊂ (−a,a), the shifts for s>=2a are disjoint. Their scalar tail has the EXACT integral

\[
\int_{2a}^\infty\rho(s)\,ds
=\operatorname{artanh}(e^{-a})+\arctan(e^{-a}).
\tag{28.9}
\]

Define c_a=g(0)+2[artanh(e^(-a))+arctan(e^(-a))]. Then **all Gamma integration reduces to one compact interval**:

\[
\boxed{W_\Gamma(f,f)=c_a\|f\|_2^2+
\int_0^{2a}\rho(s)\|f-S_sf\|_2^2\,ds.}\tag{28.10}
\]

Polarization yields, for any pole-null finite physical basis vectors φ_i,φ_j (interior or anchor), the **complete original signed Weil matrix**

\[
\boxed{H_{ij}=c_aM_{ij}+
\int_0^{2a}\rho(s)[2M_{ij}-U_{\varphi_i\varphi_j}(s)]ds
-\sum_{n\in\mathscr P}w_n U_{\varphi_i\varphi_j}(\log n),}\tag{28.11}
\]

where U_fg(s)=<f,S_sg>+<f,S_−sg>. For w_i,w_j all U values are the explicit sine expressions (28.5) with rank-two anchor substitutions. Cross terms with u_1,u_2 are likewise elementary finite-interval sine/exponential overlaps. This is a **source-preserving one-dimensional finite-interval physical quadrature problem**; numerical evaluation of digamma at a Fourier grid is no longer necessary.

At s=0, rho(s)~1/(2s), but the bracket has a removable/integrable singularity. For interior w∈H¹_0, ||w−S_sw||≤s||w'||, so its Gamma integrand is O(s). For boundary anchor u, translation increments have squared L² norm O(s), giving a bounded Gamma integrand. The calculation is therefore well-defined on the CR-ALP-025 logarithmic form domain and the finite section.

## 4. Quantified near-zero and residual gates

For 0<s<=1, elementary 1−e^(−2s)>=2s e^(−2s) and e^(3/2)<5 imply rho(s)<=5/(2s). For any δ<=1 and interior w∈H¹_0,

\[
\boxed{\int_0^\delta\rho(s)\|w-S_sw\|_2^2ds
\le\frac54\delta^2\|w'\|_2^2.}\tag{28.12}
\]

For a bounded-variation zero-extended exceptional anchor u, with TV(u) its total variation and B=||u||∞, the inequalities ||u−S_su||_1<=s TV(u) and ||u−S_su||∞<=2B give

\[
\boxed{\int_0^\delta\rho(s)\|u-S_su\|_2^2ds
\le5\delta\,\mathrm{TV}(u)B.}\tag{28.13}
\]

Cross integrals follow from Cauchy–Schwarz in the positive translation-difference energy. Thus near-zero interval bounds are explicit; one does not have to numerically divide two nearly equal overlap values at s≈0.

Write the original full physical L² operator, acting on supported finite basis functions in x∈I, as

\[
\boxed{\mathcal T f(x)=c_af(x)+\int_0^{2a}\rho(s)
[2f(x)-f(x+s)-f(x-s)]ds
-\sum_{n\in\mathscr P}w_n[f(x+\log n)+f(x-\log n)].}\tag{28.14}
\]

The integral converges in L²(I) on all specified finite piecewise-smooth columns (zero-extended endpoint jumps yield an O(s^1/2) L² increment). Project this *original* operator with Π0 and choose **any** finite matrix C_N producing a proposed within-E_N approximant g=V C_N x to p=Vx. Define the physical residual form

\[
\mathcal R_{C_N}(x)=\|\Pi_0\mathcal T (Vx)-V C_Nx\|_{L^2(I)}^2.
\tag{28.15}
\]

For Q_N the true physical orthogonal complement, ||Q_N A p||² <= R_C(x), because Q_N VC_Nx=0 and Π0 is orthogonal. By CR-ALP-027's positive entire tail W(q,q)>=1/4||q||², the mathematically exact **sufficient original signed finite test** is

\[
\boxed{H_N-4R_{C_N}\succeq0\quad\Longrightarrow\quad
W(f,f)\ge0\text{ for every two-pole-null original fixed-collar }f.}\tag{28.16}
\]

This avoids constructing the operator-column Gram D_N. For optimal C_N=M_N^(−1)H_N, R_C equals the true cross-tail B_N. Approximate sparse or matrix-free C_N can also be used, at the cost of a larger residual. **NOT evaluated or certified for N=15,000.** The physical rank-two Gram and the entire infinite-tail bound remain essential; finite H>=0 alone is not sufficient.

## 5. Additional rigorous Fourier envelope for residual computation

As an independent fallback, the true Fourier multiplier m_a(r)=g(r)−2 sum_n w_n cos(r log n) satisfies

\[
|m_a(r)|\le\log(1+|r|)+16,\quad |m_a'(r)|<26.
\tag{28.17}
\]

The first bound follows from CR-ALP-024's two-sided digamma estimate and S_a<177/50; the second uses |ψ'(1/4+ir/2)|<=sum_(n>=0)(n+1/4)^(-2)<18, the factor 1/2 in the r derivative and 2(2a)S_a<15.576.

For T>=2nu_N the N−2 interior columns satisfy |hat w_j(r)|<=8 beta_j/(3sqrt(a)r²), beta_j=nu_j+lambda_j nu_anchor<=2nu_j. The two boundary anchors instead satisfy |hat u_k(r)|<=K_k/|r|, where their explicit BV constants obey **rigorous arithmetic bounds K_1<5, K_2<9**.

Define K_0(T),K_2(T) as Fourier tail quadratic-form upper envelopes of ||hat p||² and |m_a|²||hat p||² respectively, in the raw {w_j,u_k} coefficient coordinates. At **T=10^6, N=15,000**, a separate exact Fraction audit gives diagonal block ceilings:

| Fourier envelope | 14,998 interior columns | Two boundary columns |
|---|---:|---:|
| K_0 (unweighted tail) | <1/70000 | <1/14000 |
| K_2 (squared-symbol tail) | <1/70 | <3/40 |

These are **coefficient-coordinate Loewner upper bounds, not lower bounds in normalized physical mass**. Boundary Gram conditioning must be included before interpreting spectral accuracy. The direct compact translation kernel (28.10) is preferable for the low Weil form, while either physical or Fourier residual enclosures may be used for (28.16).

## 6. Actual complete signed Weil 6×6 pilot — NUMERICAL DIAGNOSTIC ONLY

The separate `cr_alp_028_full_signed_6x6_diagnostic.py` evaluates (28.11) by 45-digit mpmath quadrature for the genuine six **interior pole-null** directions w_3,...,w_8, retaining all seven prime coefficients and the complete Gamma bracket. It produces `CR-ALP-028_original_signed_6x6_diagnostic.csv` with the complete original H and true physical M. These are **NOT outward-rounded interval certificates**, and should not be promoted.

Approximate full original signed diagonals:

| Direction | Gamma numerical | Original prime numerical | Complete original W numerical |
|---|---:|---:|---:|
| w_3 | −0.64436255125 | +0.64787555550 | +0.00351300425 |
| w_4 | −0.32426100907 | +0.32566201737 | +0.00140100830 |
| w_5 | −0.00005284316 | +0.01675778936 | +0.01670494621 |
| w_6 | +0.18941092823 | −0.18395701107 | +0.00545391716 |
| w_7 | +0.38038011095 | −0.30745509607 | +0.07292501488 |
| w_8 | +0.51737294338 | −0.50050965945 | +0.01686328393 |

The *numerically computed* 6×6 generalized eigenvalues of (H,M) are approximately **2.2579e−8, 2.2779e−7, 1.1648e−4, 1.9202e−4, 1.9389e−2, 8.5761e−2**. The near-cancellation between original Gamma and prime terms makes **uncertified floating PSD especially unreliable**. A strict interval lower bound on the smallest generalized eigenvalue and the entire cross-tail residual is indispensable; this diagnostic establishes neither.

An independent high-precision spot check comparing g(r) to the infinite translation-difference integral at r=0,2,7,20 matched to about 1e−17 after truncating the integral at s=80 with 1/4-unit partitions. Those checks are *diagnostic*, while (28.7) is a separate exact analytical derivation.

## 7. Audit / acceptance ledger

- SOURCE PIN: original Zeta23.EF.literatureRHS two poles, negative paired seven Λ terms, complete Gamma bracket; receiving source blobs re-fetched and unchanged.
- HUMAN ANALYTIC: exact physical 2-boundary rank decomposition; true physical mass Gram; original paired-prime sine formula; original Gamma translation-difference identity with exact closed-form tail; finite-support signed operator formula; physical residual Schur sufficiency.
- EXECUTED exact Fraction source prime audit: `src/audits/cr_alp_028_prime_interval_and_gram_audit.py` (six outward-rational original prime diagonals; parity; pole-anchor regression) — **PASS**.
- EXECUTED exact Fraction residual envelope audit: `src/audits/cr_alp_028_residual_envelope_audit.py` (derivative, endpoint BV constants, 15k/10^6 tail envelopes, abstract residual Schur) — **PASS**.
- EXECUTED complete source-signed **6×6 DIAGNOSTIC**, not certified: `src/audits/cr_alp_028_full_signed_6x6_diagnostic.py` and accompanying CSV.
- OPEN: rigorous Gamma quadrature/enclosure of even the small full signed block; true H_15000; certified residual R_C; full fixed-collar two-pole-null positivity; arbitrary support extension; Lean compiler and `#print axioms`; F41 completeness; Connes/Suzuki transfer; RH.

**Next CR-ALP-029:** build a singularity-subtracted, validated one-dimensional quadrature for (28.11), first certify a **complete original signed 6×6 physical H/M** as an interval object, then accelerate source-attached matrix-free action toward the 15k Schur residual test. The verified 6×6 test does not itself establish entire-collar positivity. No acceptance-ledger promotions or receiving Weil-Lab writes.