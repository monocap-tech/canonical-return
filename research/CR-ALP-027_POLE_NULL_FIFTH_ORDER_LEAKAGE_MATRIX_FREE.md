# CR-ALP-027 — Exact two-pole cancellation sharpens Fourier leakage to fifth order; explicit 15,000-mode positive tail

**2026-11-23 PST — Canonical-Return research only / DRAFT / NONPROMOTING.** This note proves a new *human analytic* finite-section reduction for the **original** signed Weil quadratic form at physical support radius \(a=11/10\). The accompanying Python script checks all finite rational cutoffs and algebraic cancellation identities, but the functional-analytic proof below has **not been Lean-compiled**. No finite positive matrix, genuine residual norm, fixed-collar positivity, or RH claim is made.

## 0. Exact source contract and previous stages

Pinned native interface: `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, Git blob `0034976fe2f1e5fa9d411108736160542eb94293`, `Zeta23.EF.literatureRHS`; source-native paired prime cutoff `ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`. The original form has **both** pole residues, **negative** paired \(\Lambda(n)/\sqrt n\) samples, and the full \(\operatorname{Re}\psi(1/4+ir/2)-\log\pi\) Gamma bracket.

On \(I=(-a,a)\), \(a=11/10\), only the seven prime powers \(\mathscr P=\{2,3,4,5,7,8,9\}\) appear in the original von Mangoldt term, with \(\Lambda(p^k)=\log p\). Let \(\mathscr H_0\) be the genuine two-pole-null subspace \(F_+(f)=F_-(f)=0\) of \(L^2(I;\mathbb C)\). Keep CR-ALP-025's physical orthogonal pole projector \(\Pi_0=I-BG^{-1}B^*\), where \(G=\begin{psmallmatrix}2\sinh a&2a\\2a&2\sinh a\end{psmallmatrix}\), and its projected Dirichlet basis \(v_j=\Pi_0e_j\),
\[
 e_j(x)=a^{-1/2}\sin(\nu_j(x+a)),\quad \nu_j=j\pi/(2a).
\]
The physical projected mass matrix remains \(M_N=I-C_N^*G^{-1}C_N\), not the historical F41 physical H/V.

For the actual orthogonal complement \(q\in Q_N\mathscr H_0\), \(Q_N=I-P_{\operatorname{span}(v_1,\dots,v_N)}\), one has \(\langle q,e_j\rangle=0\) for every \(j\le N\). Thus, in the ordinary sine expansion, \(q=\sum_{j>N}c_je_j\) with \(\sum|c_j|^2=\|q\|_2^2\).

## 1. NEW: exact pole-null cancellation of the first low-frequency sine term

Set \(\varepsilon_j=(-1)^j\). The two genuine pole moments of the *unprojected* sine modes are, **exactly**,
\[
 F_\sigma(e_j)=\frac{\nu_j}{\sqrt a}\frac{e^{-\sigma a/2}-\varepsilon_j e^{\sigma a/2}}{\nu_j^2+1/4},\qquad \sigma=\pm1.
\tag{27.1}
\]
For the tail \(q=\sum_{j>N}c_je_j\), the two equations \(F_+(q)=F_-(q)=0\) (valid by \(L^2\) continuity) therefore imply the **separate** parity identities
\[
 \boxed{\sum_{\substack{j>N\\j\ \mathrm{even}}}\frac{c_j\nu_j}{\nu_j^2+1/4}=0,\qquad
 \sum_{\substack{j>N\\j\ \mathrm{odd}}}\frac{c_j\nu_j}{\nu_j^2+1/4}=0.}\tag{27.2}
\]
The Fourier transform of a single sine mode, in the paper convention \(\widehat f(r)=\int f(x)e^{irx}dx\), is
\[
 \widehat e_j(r)=\frac{e^{-iar}}{\sqrt a}\bigl(1-\varepsilon_je^{2iar}\bigr)\frac{\nu_j}{\nu_j^2-r^2},
\tag{27.3}
\]
with the continuous removable limits at \(r=\pm\nu_j\). Since the factor in parentheses depends *only on parity*, we may subtract (27.2) separately to obtain, whenever \(|r|\le R\) and \(\nu_{N+1}\ge2R\), the **exact** identity
\[
 \boxed{\widehat q(r)=\frac{e^{-iar}}{\sqrt a}\sum_{j>N}c_j\bigl(1-\varepsilon_je^{2iar}\bigr)
 \frac{\nu_j(r^2+1/4)}{(\nu_j^2-r^2)(\nu_j^2+1/4)}.}\tag{27.4}
\]
The \(j^{-1}\) leading terms present in CR-ALP-025/026 **cancel**. The residual terms are \(O((r^2+1/4)j^{-3})\). Convergence and termwise subtraction are valid: the pole-moment terms are absolutely summable by Cauchy–Schwarz and \(\sum j^{-2}<\infty\), and the residual series is even more rapidly convergent.

For \(|r|\le R\) with the separation condition, \(|1-\varepsilon_je^{2iar}|\le2\), \(\nu_j^2-r^2\ge3\nu_j^2/4\), and \(\nu_j^2+1/4\ge\nu_j^2\). Hence
\[
 |\widehat q(r)|\le\frac{8}{3\sqrt a}(r^2+1/4)\|q\|_2
 \left(\sum_{j>N}\nu_j^{-6}\right)^{1/2}.
\]
Using \(\sum_{j>N}j^{-6}\le1/(5N^5)\), Plancherel, and direct integration gives the **NEW fifth-order tail concentration theorem**
\[
 \boxed{p_R(q):=\frac1{2\pi\|q\|_2^2}\int_{-R}^{R}|\widehat q(r)|^2dr
 \le d^{(5)}_{N,R}:=\frac{4096a^5}{45\pi^7N^5}
 \left(\frac{R^5}{5}+\frac{R^3}{6}+\frac R{16}\right).}\tag{27.5}
\]
Unlike the previous \(O(R/N)\) estimate, this is \(O((R/N)^5)\) at fixed radius. This is an *exact two-pole-null consequence* and is **not valid** for arbitrary nonnull sine tails. No assumptions about the actual zeta zeros enter it.

## 2. Tightened original Gamma and seven-prime arithmetic; concrete 15,000-mode theorem

CR-ALP-026's human digamma inequality for \(|r|\ge R\ge2\) yields \(g(r)>\log R-\log(2\pi)-2/R\), and its global bound \(g(r)>-6\). The **same genuine seven source prime powers** can be enclosed more sharply without changing a coefficient:
\[
 \boxed{S_a=\sum_{n\in\mathscr P}\Lambda(n)/\sqrt n<177/50,\quad
 \log(2\pi)<37/20,\quad\log(10000)>921/100.}\tag{27.6}
\]
The exact companion audit proves the first inequality via 40-term rational atanh/log bounds and integer sqrt floors, the second via \(\pi<22/7\) and a finite **lower** Taylor sum for \(e^{37/20}>44/7\), and the third by the same rational atanh log construction. Thus the **full original** paired negative prime contribution is bounded below by \(-177/25\) physical mass, and for \(|r|\ge R=10000\), the original full Gamma bracket satisfies
\[
 g(r)> L:=\frac{921}{100}-\frac{37}{20}-\frac{2}{10000}
 =\frac{36799}{5000}.
\]
Take the explicit **N=15,000**, for which \(\nu_{N+1}>2R\) follows from \(\pi>3\) and \(3(N+1)>4aR\). The exact rational fifth-order bound evaluated with \(\pi>3\) gives \(d^{(5)}_{N,R}<1/500\). Therefore, for every \(q\in Q_{15000}\mathscr H_0\) in the original log-form domain,
\[
 \begin{aligned}
 \frac{W(q,q)}{\|q\|_2^2}
 &\ge L(1-p_R(q))-6p_R(q)-\frac{177}{25}\\
 &>L\left(1-\frac1{500}\right)-\frac6{500}-\frac{177}{25}\\
 &=\frac{632701}{2500000}>\frac14.
 \end{aligned}
\]
Thus the principal new finite-section theorem is
\[
 \boxed{\forall q\in Q_{15000}\mathscr H_0\cap\mathcal D_{\log},\qquad W(q,q)\ge\frac14\|q\|_2^2.}\tag{27.7}
\]
This improves CR-ALP-026's **500,000-mode, margin 1/20** tail to a **15,000-mode, margin 1/4** positive tail: **33⅓-fold fewer modes and fivefold stronger coercivity margin**. Both margins are lower bounds for an infinite-dimensional complement; neither certifies the low block.

Optional exact-audited size/margin options: \((R,N)=(8000,12000)\) yields \(W(q,q)>\frac1{30}\|q\|_2^2\); \((9000,14000)\) yields \(W(q,q)>\frac3{20}\|q\|_2^2\). These use the exact rational log lower bounds rather than the single coarse rational \(L\) displayed above. The 15,000-mode option is preferred for its much smaller Schur penalty.

## 3. NEW: direct Fourier-multiplier form of the *entire* source-signed prime + Gamma operator

For every pole-null compact smooth physical source \(f\), the original exact correlation has \(\widehat{k_f}(r)=|\widehat f(r)|^2\). Since all higher prime powers vanish by support, Fourier inversion of the seven retained paired correlations gives
\[
 \boxed{W(f,f)=\frac1{2\pi}\int_{\mathbb R}m_a(r)|\widehat f(r)|^2\,dr,\qquad
 m_a(r)=\operatorname{Re}\psi\!\left(\frac14+\frac{ir}{2}\right)-\log\pi
 -2\sum_{n\in\mathscr P}\frac{\Lambda(n)}{\sqrt n}\cos(r\log n).}\tag{27.8}
\]
This is precisely the **original signed** operator restricted to the fixed collar and pole-null space. All seven von Mangoldt multiplicities and both \(\pm\log n\) evaluations are retained; no positivity is inserted. The symbol is real and even, with \(m_a(r)=O(\log(2+|r|))\).

For the native physical projected finite basis, this gives a fully explicit one-dimensional integral for every signed entry
\[
 H_{ij}=\frac1{2\pi}\int_{\mathbb R}m_a(r)\widehat v_i(r)\overline{\widehat v_j(r)}\,dr.
\tag{27.9}
\]
The elementary Fourier functions \(\widehat e_j\), \(\widehat{e^{\pm x/2}1_I}\), and the rank-two projector completely specify \(\widehat v_j\). In physical space, the original prime operator is exactly the bounded sum of the seven zero-extended translation pairs, while the Gamma component is the full scalar digamma Fourier multiplier. This furnishes **matrix-free operator actions** (seven physical shifts + one Fourier multiplier + a rank-two physical projection); it does not by itself yield rigorous global numerical spectral enclosures.

## 4. Remaining ORIGINAL full-collar signed Schur gate, now 15k

With genuine physical mass \(M_N\), signed matrix \(H_N=(W(v_i,v_j))\), operator-column Gram \(D_N=(\langle\Pi_0Av_i,\Pi_0Av_j\rangle)\), and residual
\[
 B_N=D_N-H_NM_N^{-1}H_N\succeq0,
\]
we obtain from (27.7) the concrete sufficient original signed certificate
\[
 \boxed{H_{15000}-4B_{15000}\succeq0
 \quad\Longrightarrow\quad W(f,f)\ge0\quad
 \forall f\in\mathscr H_0\cap\mathcal D_{\log}.}\tag{27.10}
\]
**NOT evaluated.** A computationally useful but potentially more conservative *matrix-free scalar alternative* is to certify
\[
 \lambda_{\min}(M_N^{-1/2}H_NM_N^{-1/2})\ge 4t_N^2,
 \quad t_N\ge\|Q_NA P_N\|_{L^2\to L^2},
\tag{27.11}
\]
where \(t_N\) is a rigorously *upper-enclosed* residual operator norm, not a numerical estimate without a tail bound. The finite original signed matrix and its cross-tail residual are still unknown. Positivity of \(H_N\) alone is insufficient. Failure of (27.10) or (27.11) is inconclusive. A certified negative finite vector of \(H_N\) would yield a negative original signed form test after pole-preserving smooth approximation.

## 5. Custody, audit scope and next steps

- Audit: `src/audits/cr_alp_027_pole_null_fifth_order_audit.py` (Python stdlib `Fraction`/integer square roots). Executed locally, checks exact source seven weights, tighter logarithmic/Taylor constants, pole-parity rational resolvent cancellation regressions, Fourier frequency separation, fifth-order \(d^{(5)}\), tail margin and the 12k/14k/15k tradeoff.
- **Human analytical, not Lean compiled:** physical sine/pole relation, continuum resolvent cancellation, Plancherel integration, digamma lower bound, domain stability, and full original Fourier formula. No new receiving `#print axioms`/CI acceptance.
- **Open:** the 15,000-dimensional signed low-frequency matrix, its actual physical Gram and cross-tail residual intervals, the entire fixed collar PSD, arbitrary support extension, completeness of frozen 242 F41 generators, global RH.
- **Governance:** Canonical-Return PR #3 stays DRAFT/NONPROMOTING; no receiving Weil-Lab, Formalization, Hodge, provenance, F41/F42, or RH acceptance mutation.

**Next CR-ALP-028:** implement a *verified* matrix-free \(H_N\)-action and a rigorous residual-norm enclosure that respects the rank-two physical metric, full digamma symbol and seven prime shifts. Use (27.10)/(27.11) as the acceptance gate; do not confuse an approximate lowest Ritz value with a universal spectral lower bound.