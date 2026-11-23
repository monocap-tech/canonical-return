# CR-ALP-024 — The ORIGINAL signed Weil form has finite negative index on every fixed support collar

**Research-only / Canonical-Return, 2026-11-23 PST.** **Status:** directly reopened original IV010/IV011/Zeta23 definitions; new analytic lower and upper logarithmic Fourier-form estimates, an **effective negative-index bound** for the original radius 11/10, and a separately executed exact-rational audit. **NOT** a new Lean compilation, a positive-definiteness theorem on the complete collar, an F41 numerical H/V acceptance, arbitrary support-extension positivity, a Connes/Suzuki domain identification, or RH.

## 0. Exact native-source and audit custody

- Source of the actual explicit formula: `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, blob **`0034976fe2f1e5fa9d411108736160542eb94293`**. Its `Zeta23.EF.literatureRHS k` is the **two pole evaluations minus both-sign von Mangoldt terms plus the full digamma/Gamma integral**. The genuine `gammaBracket r` is \(a(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi\).
- Literal receiving mixed pairing: `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportExtension.lean`, blob **`5ef51fd37337c466620eafccee2a27b77955746e`**, defines `neutralIV010Mixed` as that `literatureRHS` evaluated on the true physical correlation.
- Literal IV011 finite support prime formula: `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob **`14d8f01ce34f764900eccc68399d7323b4f4a25c`**, establishes paired \(\pm\log n\) sampling and original pole − prime + Gamma split.
- Native digamma-series input: `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_GammaFacts_Series.lean` / theorem sources on `GammaFacts/Series` and `Zeta23.WeilEF.digamma_conj`; the classical partial-fraction series of \(\psi\) and the classical special value \(\psi(1/4)=-\gamma-\pi/2-3\log2\) are the stated external analytic identities. They are **not** fresh Lean theorems of CR-ALP-024.
- Original reviewed Phase IV support collar: \(I=(-a,a)\), \(a=11/10\). The frozen historical F41 selected **physical** H/V is not used as if it were the signed Weil operator.
- Reopened and re-executed the actual uploaded CR-ALP-023 exact audit: PASS, 74 negative / 167 positive original oscillator physical moments. That historical finite-mode result is independent of this new **arbitrary-source** theorem.

## 1. An exact fixed-radius operator decomposition without changing the Weil sign

Let \(f\in C_c^\infty(I;\mathbb C)\), extend it by zero, put \(k=f*\widetilde f\) with \(\widetilde f(x)=\overline{f(-x)}\), and adopt the *source's* paper Fourier convention
\[
\widehat f(r)=\int_{\mathbb R}f(x)e^{irx}dx,\qquad
\frac1{2\pi}\int|\widehat f(r)|^2dr=\|f\|_2^2.
\]
Writing the actual original signed EF directly gives
\[
\boxed{\begin{aligned}
W(f,f)={}&\frac1{2\pi}\int_{\mathbb R}a(r)|\widehat f(r)|^2dr
\\&+2\operatorname{Re}\bigl(F_-(f)\overline{F_+(f)}\bigr)
-\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\bigl(k(\log n)+k(-\log n)\bigr),
\end{aligned}}\tag{24.1}
\]
where \(F_\pm(f)=\int f(x)e^{\pm x/2}dx\) and \(a(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi\). This preserves **both poles, the original negative prime sign, the exact von Mangoldt multiplicities, and the complete Gamma term**. It is not a new PSD kernel.

Since \(\operatorname{supp}k\subset[-2a,2a]\) and the rational interval check proves
\[\boxed{\log9<2a=11/5<\log10,}\]
only the **seven** prime powers \(n\in\{2,3,4,5,7,8,9\}\) contribute. With the **correct** \(\Lambda(4)=\Lambda(8)=\log2\), \(\Lambda(9)=\log3\), put
\[
S_a=\sum_{n\in\{2,3,4,5,7,8,9\}}\frac{\Lambda(n)}{\sqrt n}.
\]
An exact rational atanh-log / integer-square-root enclosure proves \(\boxed{S_a<4}\). Cauchy–Schwarz for translations gives \(|k(\pm\log n)|\le\|f\|_2^2\), hence
\[
\boxed{|W_{\rm prime}(f,f)|\le2S_a\|f\|_2^2<8\|f\|_2^2.}\tag{24.2}
\]
The two original poles obey \(|F_\pm|^2\le2\sinh(a)\|f\|_2^2\), so
\[
|W_{\rm pole}(f,f)|\le4\sinh(a)\|f\|_2^2<8\|f\|_2^2;\tag{24.3}
\]
\(e^{11/10}<4\) is separately audited using exact Taylor-series rational tails. For the two-pole-null subspace \(F_+=F_-=0\), the residue term vanishes exactly.

## 2. Elementary logarithmic growth bound on the ACTUAL Gamma bracket

**New analytic lemma (with classical digamma identities as input):** for every \(r\in\mathbb R\),
\[
\boxed{\log(1+|r|)-8\ \le\ a(r)\ \le\ \log(1+|r|)+3.}\tag{24.4}
\]

**Detailed derivation.** Put \(\alpha=1/4,t=|r|/2\). The classical digamma partial-fraction series yields
\[
\operatorname{Re}\psi(\alpha+it)
=-\gamma+\sum_{n=0}^\infty\left(\frac1{n+1}-\frac{n+\alpha}{(n+\alpha)^2+t^2}\right).
\]
For \(t\ge1\), set \(N=\lfloor t\rfloor\ge1\). On \(n\le N\), replace \((n+\alpha)/((n+\alpha)^2+t^2)\) by its upper bound \((n+\alpha)/t^2\); the sum of these upper bounds is at most
\[
\frac{(N+1)(N/2+\alpha)}{t^2}\le\frac{(N+1)(2N+1)}{4N^2}\le\frac32.
\]
On \(n>N\), each summand is at least \(-3/[4(n+1)(n+\alpha)]\ge-3/(4n^2)\), and the tail loss is at most \(3/(4N)\le3/4\). Since \(H_{N+1}\ge\log(N+2)>\log t\) and \(\gamma<1\),
\[
\operatorname{Re}\psi(\alpha+it)>\log t-\frac{13}{4}.
\]
The inequalities \(\log((1+|r|)/t)\le\log3<2\) and \(\log\pi<3/2\) imply the lower bound in (24.4) for \(|r|\ge2\). On \(|r|<2\), each partial-fraction summand increases with \(t\); the exact quarter-value \(\psi(1/4)=-\gamma-\pi/2-3\log2\) gives \(a(r)\ge a(0)>-6\), while \(\log(1+|r|)<2\). This completes the lower bound.

For the upper bound, split the same series at \(N\). Terms through \(N\) are bounded above by \(1/(n+1)\), and the **positive part** of each tail summand is at most \(t^2/n^3\). As \(\sum_{n>N}n^{-3}\le1/(2N^2)\) and \(t<N+1\), the upper tail is \(\le2\); \(H_{N+1}\le1+\log(N+1)\). Thus \(\operatorname{Re}\psi(\alpha+it)\le3+\log(1+t)\), and subtracting \(\log\pi\ge0\) proves the upper bound when \(t\ge1\). For \(0\le t<1\), split at \(N=1\), obtaining \(\operatorname{Re}\psi(\alpha+it)\le H_2+1/2=2\). This completes (24.4). No unverified positivity assumption on zeta zeros is used.

Set the **genuine logarithmic Fourier energy**
\[
\boxed{\mathcal E_{\log}(f)=\frac1{2\pi}\int_{\mathbb R}\log(1+|r|)|\widehat f(r)|^2dr.}\tag{24.5}
\]
Plancherel and (24.4) give
\[
\boxed{\mathcal E_{\log}(f)-8\|f\|_2^2\le W_\Gamma(f,f)\le\mathcal E_{\log}(f)+3\|f\|_2^2.}\tag{24.6}
\]
Therefore the **complete original signed** Weil quadratic form on the fixed radius satisfies
\[
\boxed{\begin{aligned}
W(f,f)&\ge\mathcal E_{\log}(f)-16\|f\|_2^2&&\text{if }F_+(f)=F_-(f)=0,\\
W(f,f)&\ge\mathcal E_{\log}(f)-24\|f\|_2^2&&\text{for arbitrary }f\in C_c^\infty(I).
\end{aligned}}\tag{24.7}
\]
The bounds are one-sided coercivity **modulo physical \(L^2\)**, not positive definiteness.

## 3. New effective FINITE negative-index theorem (all old-core sources)

Let \(V\) be **any linear subspace** of the smooth, original-physical pole-null sources supported in \(I\), such that \(W(f,f)\le0\) for **every** \(f\in V\). The definition does not assume that \(V\) is spanned by the F41 selected generators. For unit \(f\in V\), (24.7) implies
\[
\mathcal E_{\log}(f)\le16.
\]
For \(R>0\) put \(\mathrm{high}_R(f)=(2\pi)^{-1}\int_{|r|>R}|\widehat f(r)|^2dr\). Then
\[
\boxed{\mathrm{high}_R(f)\le\frac{16}{\log(1+R)}.}\tag{24.8}
\]
Choose the **exact integer** cutoff \(\boxed{R_0=10^{14}}\). The exact atanh-series enclosure \(\log10>16/7\) gives
\[
\log(1+R_0)>14\log10>32,
\]
so **every** unit vector of \(V\) has **more than one-half** of its physical Fourier energy inside \([-R_0,R_0]\).

The positive time-frequency concentration operator \(T_{a,R}=1_I\mathcal F^{-1}1_{[-R,R]}\mathcal F1_I\) is trace class, with **exact** trace
\[
\boxed{\operatorname{tr}T_{a,R}=\frac{2aR}{\pi}.}\tag{24.9}
\]
For any \(m\) orthonormal vectors in \(V\), summing their strictly >1/2 low-frequency masses and using \(T\ge0\) yields
\[
\frac m2<\operatorname{tr}T_{a,R_0}=\frac{2aR_0}{\pi},\qquad
m<\frac{4aR_0}{\pi}<\frac{22}{15}10^{14}
\]
(using \(a=11/10,\pi>3\)). Hence:
\[
\boxed{\dim_{\mathbb C}V\le146{,}666{,}666{,}666{,}666.}\tag{24.10}
\]
This is a **universal, completely explicit finite upper bound** on the dimension of any nonpositive subspace in the original two-pole-null Schwartz collar \(I=(-11/10,11/10)\). It is deliberately conservative; it is **not** a count of actual negative eigenvalues and **not** a positivity certificate. The same proof with the unrestricted \(-24\) lower bound and \(R_1=10^{21}\) proves
\[
\boxed{\dim_{\mathbb C}V_{\rm unrestricted}\le
1{,}466{,}666{,}666{,}666{,}666{,}666{,}666}\tag{24.11}
\]
for any nonpositive subspace of the **full**, pole-nonnull fixed collar.

**Functional-analytic corollary.** On the closure of \(C_c^\infty(I)\) in the Hilbert norm \(\|f\|_2^2+\mathcal E_{\log}(f)\), the original signed quadratic form is closed and lower semibounded after a constant physical-mass shift: (24.4) supplies two-sided form-norm equivalence, while (24.2)–(24.3) are bounded physical \(L^2\) perturbations. The weighted log-frequency norm embeds **compactly** into \(L^2(I)\) (bounded physical support + uniformly decaying Fourier tails, Riesz–Kolmogorov). Its associated self-adjoint Friedrichs form operator therefore has **compact resolvent** and only finitely many nonpositive eigenvalues counting multiplicity. These are standard functional-analysis implications of the proved inequalities, **not** a newly formalized Connes/Suzuki operator-domain identification. One may also rely solely on the explicit dimension argument (24.8)–(24.11), which does not assume completion identification.

## 4. Why this is a REAL advancement but NOT a full-support theorem

This is a positive result about the **actual fixed-collar source-wide infinite-dimensional form**, not just a finite hand-selected packet span: all sufficiently high Fourier-frequency physical energy is costly in the original Gamma term, and the original finite-prime and two-pole terms are bounded perturbations. A negative/nonpositive obstruction cannot have arbitrarily large linear dimension at fixed support.

**However:** the bound (24.10) is much larger than 242, gives **no eigenvectors**, does not show the negative index is zero, and gives no identification with F41's historical finite selected source span. In particular, the exact F41 physical 241 nonzero packet moments or the positive 2D/4D source blocks of CR-ALP-018–023 do not establish positivity of the remaining low-frequency subspace. Enlarging the support grows the prime cutoff and changes the bounded-perturbation constants; no uniform support-independent conclusion follows from (24.10).

The even/odd Schwartz source blocks are orthogonal under the full original signed Weil form by the **analytical** parity theorem of CR-ALP-020–023. Both parity blocks have the same finite-index property separately. **Lean parity certification remains OPEN**: no \`lean\` or \`lake\` executable was available locally, no receiving file was compiled, and no \`#print axioms\` was obtained. No new axiom is assumed in the analytic proof.

## 5. Gate ledger / reproduction

| Obligation | CR-ALP-024 status |
|---|---|
| Actual source EF poles − paired-primes + full Gamma | **PINNED** |
| Correct cutoff \(\{2,3,4,5,7,8,9\}\) / \(S_a<4\) | **RATIONAL PASS** |
| Genuine Gamma \(\log(1+|r|)\) two-sided estimate | **ANALYTIC PROOF + RATIONAL GATES** |
| Source-wide complete-form lower bounds \(-16\), \(-24\) | **ANALYTIC PASS** |
| Effective fixed-collar finite negative-index bound | **ANALYTIC + RATIONAL TRACE PASS** |
| Compact-resolvent extension in logarithmic form domain | **STANDARD FUNCTIONAL ANALYSIS COROLLARY** |
| Exact uploaded CR-ALP-023 audit rerun | **LOCAL PASS** |
| All-old collar signed Weil positivity / zero negative index | **NOT ESTABLISHED** |
| Actual frozen F41 242 physical H/242×4 V promotion | **NOT PERFORMED** |
| Receiving Lean compile / axiom audit | **NOT PERFORMED** |
| RH / arbitrary support-extension theorem | **NOT CLAIMED** |

**Standalone exact audit:** \`src/audits/cr_alp_024_fixed_collar_log_spectral_audit.py\`, Python stdlib Fraction and integer root, independently executed **PASS**. This checks all numerical cutoffs, seven actual von Mangoldt terms, Gamma special-value estimates and bounds, pole Taylor comparison, and exact rational finite-index trace threshold. It does **not** mechanically prove the infinite digamma series, Plancherel or trace formula. Those are the detailed separate human analytic derivations above.

**Next CR-ALP-025:** construct a **computably convergent finite-section Schur/min-max test** for the actual \(\log\)-form and compact bounded-prime operator, using explicit physical test basis and certified cross-term estimates. Demonstrate, without invoking F41 physical sign counts as a substitute, what finite matrix and **tail operator norm** would be sufficient to certify an entire fixed collar. The original RH/unbounded-support issue remains untouched until that operator-tail gate is genuinely discharged.

**Governance:** Remain in Canonical-Return draft PR #3. No receiving Weil-Lab, F41/F42, RH/provenance or Hodge ledger changes.