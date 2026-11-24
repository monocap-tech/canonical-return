# CR-ALP-026 — A sharper original Gamma estimate reduces the positive-tail finite section to 500,000 modes

**2026-11-23 PST — Canonical-Return, research only, DRAFT/NONPROMOTING.** This is a human analytical argument coupled to a standalone exact-rational audit. The continuum digamma, Fourier, and closed-form steps have **not** been Lean-compiled. No evaluation of the resulting signed-Weil finite matrices is claimed; fixed-collar positivity, F41 completeness, all-support positivity and RH remain open.

## 0. Source contract retained

The original pinned `Zeta23.EF.literatureRHS` in `monocap-tech/weil-lab`, `research/iv068-original-source-kernel-gate`, file `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean` (Git blob `0034976fe2f1e5fa9d411108736160542eb94293`) has **two original pole residues, the negative paired original von Mangoldt terms, and the complete original** `gammaBracket(r)=Re digamma(1/4+ir/2)-log(pi)`. The genuine fixed physical support radius is `a=11/10`, so precisely `n=2,3,4,5,7,8,9` may occur in the prime part. Work on the original physical Hilbert space `L²(-a,a;C)` and its two-pole-null subspace H0.

From CR-ALP-025 retain the physical Dirichlet sine basis

\[
 e_j(x)=a^{-1/2}\sin(j\pi(x+a)/(2a)),\quad \nu_j=j\pi/(2a),\quad v_j=\Pi_0e_j,
\]

where Π0 is the **exact physical L²** projection onto the two-pole-null space. Let P_N be the physical projection onto E_N=span(v_1,...,v_N) and Q_N=I-P_N on H0. For q in Q_N H0, the first N ordinary sine coefficients vanish. When ν_{N+1}≥2R, CR-ALP-025's elementary sine Fourier estimate gives

\[
 p_R(q):=\frac{1}{2\pi\|q\|_2^2}\int_{-R}^R|\widehat q(r)|^2\,dr
 \le \delta_{N,R}:=\frac{256aR}{9\pi^3 N}.
\tag{26.1}
\]

For q=0, the statements below are trivial. This is a **uniform bound on every vector of the infinite-dimensional complement**, not a selected 242-mode statistic.

## 1. New sharp high-frequency digamma lower bound: an elementary total-variation lemma

Take α=1/4, t≥1, and the positive real rational function

\[
 h_t(x)=\frac{x+\alpha}{(x+\alpha)^2+t^2},\qquad x\ge0.
\]

The classical digamma partial-fraction identity (the same input used in CR-ALP-024) implies

\[
 \operatorname{Re}\psi(\alpha+it)
 =\lim_{M\to\infty}\left(\log M-\sum_{n=0}^{M-1}h_t(n)\right).
\]

For each integer n, the elementary fundamental-theorem-of-calculus estimate bounds the error between the left endpoint and cell integral by total variation on that cell:

\[
 \left|h_t(n)-\int_n^{n+1}h_t(x)\,dx\right|
 \le\int_n^{n+1}|h_t'(x)|\,dx.
\]

Since t≥α, h_t first increases up to x=t−α, attaining 1/(2t), and then decreases to 0. Thus

\[
 \int_0^\infty |h_t'(x)|\,dx
 =\frac1t-\frac{\alpha}{\alpha^2+t^2}\le\frac1t.
\]

Also

\[
 \int_0^M h_t(x)\,dx
 =\tfrac12\log\frac{(M+\alpha)^2+t^2}{\alpha^2+t^2}.
\]

Passing to the limit therefore proves the **continuum analytical inequality**

\[
 \boxed{\operatorname{Re}\psi(\tfrac14+it)
 \ge\tfrac12\log(\tfrac1{16}+t^2)-\tfrac1t,
 \quad t\ge1.}\tag{26.2}
\]

For the **unmodified original Gamma bracket** `g(r)=Re digamma(1/4+ir/2)-log(pi)` this yields, whenever |r|≥R≥2,

\[
 g(r)\ge\log R-\log(2\pi)-2/R
 >\boxed{\log R-2-2/R}.\tag{26.3}
\]

The last inequality uses π<22/7 and e²>7>44/7, hence log(2π)<2. At **all** real frequencies, the same series shows Re ψ(1/4+it) is nondecreasing in |t|; the audited special value at zero gives

\[
 \boxed{g(r)\ge g(0)>-6\qquad(r\in\mathbb R).}\tag{26.4}
\]

No claim that the Gamma bracket is everywhere positive is made.

## 2. Sharper exact original prime budget

The genuine seven weights obey the **strictly sharpened** rational bound

\[
 \boxed{S_a=\sum_{n\in\{2,3,4,5,7,8,9\}}\Lambda(n)/\sqrt n<71/20.}\tag{26.5}
\]

This is certified by 40-term rational atanh logarithm enclosures and integer square-root floors at denominator 10^6 in the standalone exact audit. The full original paired-prime contribution consequently satisfies

\[
 W_{\rm prime}(q,q)\ge-2S_a\|q\|_2^2>-\frac{71}{10}\|q\|_2^2.
\tag{26.6}
\]

The two original residues vanish exactly for q in H0. No sign has been changed.

## 3. **Explicit 500,000-mode tail strictly positive**

Set the much smaller exact integers

\[
 \boxed{R=30{,}000,\qquad N=500{,}000.}\tag{26.7}
\]

The frequency separation ν_{N+1}≥2R follows already from π>3 and a=11/10. The exact atanh enclosures verify `log 3>1`, `log 10>16/7`; hence

\[
 \log R=\log3+4\log10>71/7,
 \quad g(r)>L:=\frac{57}{7}-\frac1{15000}=\frac{854993}{105000}
 \quad(|r|\ge R).
\]

The Fourier leakage bound (26.1) gives, using π>3,

\[
 p_R(q)\le\delta_{N,R}<d:=\frac{256(11/10)\,30000}{243\,500000}
 =\frac{704}{10125}.
\]

Apply (26.3) to the high-frequency mass (at least 1−p_R), (26.4) to the remaining low-frequency mass, and (26.6) to the unchanged prime sum. For every nonzero q in the original log-form domain and Q_N H0,

\[
\begin{aligned}
 \frac{W(q,q)}{\|q\|_2^2}
 &\ge L(1-p_R(q))-6p_R(q)-71/10\\
 &>L(1-d)-6d-71/10\\
 &=\boxed{\frac{63181553}{1063125000}>\frac1{20}}.
\end{aligned}\tag{26.8}
\]

Therefore the **entire infinite-dimensional physical orthogonal complement** of the explicit 500,000-mode pole-null prefix has a certified-in-arithmetic, **analytically proved** lower margin:

\[
 \boxed{W(q,q)\ge\tfrac1{20}\|q\|_2^2
 \quad\forall q\in Q_{500000}\mathscr H_0\cap\mathcal D_{\log}.}\tag{26.9}
\]

This replaces CR-ALP-025's 250-trillion-mode positive tail with a **500-million-fold smaller specified prefix**. It still contains 500,000 modes, and a dense 500,000×500,000 numerical interval matrix would be very costly.

**Margin/size tradeoff.** At the same R=30,000, choosing N=1,000,000 instead gives the exact rational lower margin 585899089/1063125000>1/2, so the corresponding original signed tail obeys `W(q,q)≥(1/2)||q||²`. This doubles the prefix but improves the later conservative Schur penalty tenfold.

## 4. Fully original signed finite Schur obligations — still unevaluated

Retain precisely the physical metric and source-signed compressed operator of CR-ALP-025:

\[
 M_N=(\langle v_i,v_j\rangle),\quad
 H_N=(W(v_i,v_j)),\quad z_j=\Pi_0 A v_j,\quad
 D_N=(\langle z_i,z_j\rangle),\quad
 B_N=D_N-H_NM_N^{-1}H_N\succeq0.
\]

Here A is the closed-form operator of the actual original signed Gamma **minus** paired seven-prime shifts, compressed to the two-pole-null Hilbert space. The mass matrix is `M_N=I-C_N*G^{-1}C_N` with the genuine pole Gram G, not an F41 substitute.

Because Q_N A Q_N≥µ I by (26.9), the infinite-dimensional Schur estimate gives the following **sufficient**, genuinely finite certificate:

\[
 \boxed{H_{500000}-20B_{500000}\succeq0
 \quad\Longrightarrow\quad
 W(f,f)\ge0\ \forall f\in\mathscr H_0\cap\mathcal D_{\log}.}\tag{26.10}
\]

The alternative 1,000,000-mode version is `H_1000000−2B_1000000⪰0`. **Neither finite matrix inequality has been evaluated.** Failure of either sufficient condition alone would not imply a negative Weil vector. A genuinely certified negative vector of H_N, after pole-preserving form approximation, would furnish a negative smooth source.

Reflection reduces each physical matrix into even/odd blocks as in CR-ALP-025, but does not reduce the mathematical obstacle to 242 F41 modes. Effective matrix-free numerical methods with rigorous residual enclosures may be preferable to storing the dense 500k×500k table.

## 5. Reproduction, custody and next gate

Standalone script: `src/audits/cr_alp_026_sharp_gamma_tail_audit.py`. Run with Python stdlib only. It **checks** the original seven weights and tightened rational bound, independent atanh log comparisons, π/exponential comparison, frequency separation, leakage bound, exact µ fractions, both N options and finite Schur algebra. It does **not** prove the digamma infinite series, Fourier concentration, Hilbert-form closure, or calculate the genuine signed H_N/B_N.

**CR-ALP-026 result:** 500,000 explicit projected sine directions suffice to isolate every possible nonpositive direction of the full original **two-pole-null** fixed collar, via a uniform strictly positive **infinite tail**. This is a major quantitative improvement but **not** original-collar positivity. Unrestricted two-pole residues, growing/unbounded supports, global RH, and F41 242-mode sufficiency remain OPEN. No receiving Weil-Lab, Formalization, F41/F42, provenance, Hodge or RH acceptance is changed.

**CR-ALP-027:** pursue the best rigorous numerical Schur gate using source-specific symbolic prime matrices, spectral-interval bounds and certified matrix-free residual norms, or reduce the prefix further by exact Gamma-symbol and finite-time concentration spectral estimates. Do not promote without the complete original signed tail cross-term control.