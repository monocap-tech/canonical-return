# CR-ALP-004 — Quantitative (M) → smooth signed-Weil lower bound, with upstream Lean dependency

**Lineage:** Canonical-Return / CR-ALP; continuation of 001–003.
**Date:** 2026-11-23 PST.
**Status:** ANALYTIC / CONDITIONAL TRANSFER; direct-source Lean inspection completed; NO NEW RH PROOF, NO NEW LEAN ACCEPTANCE, NO THETA-LIFT ↔ WEIL OPERATOR INTERTWINER.

## 0. Fixed source contracts

- Alpöge, *A remark* (20 November 2026), p. 3 item 5 and p. 4 item 8: https://alpo.ge/remark.pdf.
- Rubinstein et al., *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, **v2 working draft 22 November 2026**:
  https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf.
  Used **exactly**: (I.2)–(I.4); (I.29)–(I.30); (II.9), (II.13)–(II.19), (II.24), (II.28)–(II.36); Part III §III.11.
  **Critical:** (M) is quoted/reduced from the Lean development in v2; not proved outright in that human-readable draft.
- Public source mirror at frozen commit **7cbf7a0c32fca1badf05931879a37c3f8a2b03b4**, repository https://github.com/CocoSgt/riemann-hypothesis-lean. The precise source paths and declarations appear below. Inspected source, **not compiled here**.
- CR-ALP-002: analytic identity for the full original signed zeta Weil form on the ramp; CR-ALP-003: fixed-\(X\) smooth mollifier passage.

The central question now is **quantitative**: does the source implication (M) become a meaningful lower bound for an **admissible smooth compact** original Weil test, with smoothing error smaller than square-root scale? The answer is **YES, conditional on (M) and Estimates B/C, with explicit pole-preserving smooth tests**. It is an analytic reexpression, not an independent proof of (M).

## 1. Exact conditional chain with the actual quantifiers

Fix \(0<\varepsilon<1/52\). For the designated scales \(x\in[Y,Y^{3/2}]\), Rubinstein chooses **auxiliary primes \(q,d\) depending on \(Y\)** satisfying Estimates B/C, and for each \(x\), a finite abelian group \(H\) of cardinality \(n\), sets \(T,T_2\), \(m\) maps, cell weights \(w_{u,j,i}\), and cutoffs \(a_t\). Put
\[
\delta_H=|T_2|/n,\qquad B=|T|mb\delta_H/2>0.
\]
We intentionally denote the **column density** by \(\delta_H\) and the **mollifier width** by \(\eta_x\) below, to avoid the collision in prior notes.

Define the actual mean excess and its scale-dependent error:
\[
A(x)=\frac1n\sum_{u\in H}\sum_{t\in T}\sum_{j=1}^m(P_{u,t,j}-M_{u,t,j}).
\]
The principal input (M) is
\[
A(x)\le \mathcal E,\qquad \frac{\mathcal E}{B}\le
x^{C\varepsilon}\left(\frac xm+\frac{xm}{k}+m+\frac{k}{n^2}\right)
=:\mathfrak b(x).
\tag{4.1}
\]
Here \(\mathcal E\ge0\), and \(C\) is the source's absolute exponent constant; \(m,k,n\) are **not free choices** but produced by the arithmetic construction with the size constraints.

Rubinstein (II.36), **using Estimate C and the weight/argument replacement bounds**, proves
\[
A(x)=B R_{q,d}(x)+\Delta(x),\quad
|\Delta(x)|\le C_\varepsilon B(1+\log x).
\tag{4.2}
\]
Thus
\[
R_{q,d}(x)\le\mathfrak b(x)+C_\varepsilon(1+\log x).
\tag{4.3}
\]
The exact character identity (II.9), retaining the three twists, is
\[
R_{q,d}(x)=G(x)+G_{\chi_q}(x)-G_{\chi_d}(x)-G_{\chi_q\chi_d}(x).
\]
**Estimate B**, for the **same** auxiliary \(q,d\), supplies
\(|G_{\chi_q}|,|G_{\chi_d}|,|G_{\chi_q\chi_d}|\le x^{1/2+\varepsilon}\);
therefore
\[
\boxed{G(x)\le
\mathfrak b(x)+C_\varepsilon(1+\log x)+3x^{1/2+\varepsilon}.}
\tag{4.4}
\]
There is **only an upper bound**, not a two-sided estimate. No pointwise positivity of the weighted table, its cell Gram kernels, or the twisted zeta Weil forms was used to derive (4.4).

Let
\[
C_\zeta=2+\gamma-\log(4\pi)>0,\quad
D(x)=\sum_{h\ge1}\frac{x^{-2h}}{(2h)(2h+1)}.
\]
The complete original pole–prime–Gamma identity from CR-ALP-002 is
\[
\mathcal W(k_x)=C_\zeta(x+1)-2G(x)+2\log(2\pi)-2D(x),
\quad
k_x(y)=\big(xe^{-|y|/2}-e^{|y|/2}\big)_+.
\tag{4.5}
\]
Combining (4.4) and (4.5) gives the **exact one-sided conditional projection**:
\[
\boxed{
\begin{aligned}
\mathcal W(k_x)\ge{}&
C_\zeta(x+1)+2\log(2\pi)-2D(x)\\
&-2\mathfrak b(x)-2C_\varepsilon(1+\log x)
-6x^{1/2+\varepsilon}.
\end{aligned}}
\tag{4.6}
\]
This is an inequality **for a specific \(x\)-dependent test family**, conditional on arithmetic inputs (M), B and C, not global positivity for all original Weil sources.

### 1a. Critical balance is preserved without reinterpretation

From the source's stated scales
\[
m\asymp x^{1/2-\varepsilon/16},\
k\asymp x^{9/2+10\varepsilon},\
n\asymp x^{2+8\varepsilon},
\]
the four cost exponents, before the common \(x^{C\varepsilon}\) factor, are
\[
\frac xm:\quad \tfrac12+\tfrac\varepsilon{16},\qquad
\frac{xm}{k}:\quad -3-\tfrac{161}{16}\varepsilon,\qquad
m:\quad \tfrac12-\tfrac\varepsilon{16},\qquad
\frac{k}{n^2}:\quad \tfrac12-6\varepsilon.
\]
The \(\sqrt{x}\) scale is exactly the balance of \(x/m\) and \(m\), with the \(k/n^2\) contribution at similar order. The numerical powers here are the **orders printed in the source**, not a free optimization of the Weil test.

## 2. New uniform smooth-source bridge: the cusp is NOT the decisive loss

The fixed-\(x\) continuity statement of CR-ALP-003 left an apparent obstacle: could smoothing the nonsmooth ramp destroy the \(x^{1/2+\varepsilon}\)-quality bound?

We can now eliminate **that particular obstacle** for the *original scalar signed Weil functional*, by choosing a source-dependent shrinking mollifier.

Let \(\rho\ge0\) be any real even \(C_c^\infty(-1,1)\) function of unit integral. For \(x\ge4\), set
\[
\eta_x=x^{-2},\qquad
f_x(u)=e^{u/2}\mathbf1_{[0,\log x]}(u),\quad
f_x^{\mathrm{sm}}=\rho_{\eta_x}*f_x,\quad
k_x^{\mathrm{sm}}=f_x^{\mathrm{sm}}*\widetilde f_x^{\mathrm{sm}}.
\]
Then \(f_x^{\mathrm{sm}}\in C_c^\infty\), and
\[
k_x^{\mathrm{sm}}=\sigma_{\eta_x}*k_x,\qquad
\sigma_\eta=\rho_\eta*\widetilde\rho_\eta.
\]
Translate \(f_x^{\mathrm{sm}}\) by \(-\frac12\log x\) if support must lie symmetrically in \((-a,a)\); the autocorrelation does not change. The support radius required is
\[
a_x>\tfrac12\log x+x^{-2},
\]
so **support still escapes to infinity**. This is not uniform-support positivity.

### 2a. Uniform elementary error bounds

Put \(K_x=(x+1)/2\). The exact ramp \(k_x\) has Lipschitz constant at most \(K_x\). Since \(\sigma_\eta\) is a probability measure supported in \([-2\eta,2\eta]\),
\[
\|k_x^{\mathrm{sm}}-k_x\|_\infty\le 2K_x\eta=(x+1)\eta.
\]
The difference has Lipschitz constant at most \(2K_x\).

**Prime part.** Only \(n\le N_x:=xe^{2\eta}\) contribute, and trivially \(\Lambda(n)\le\log n\). Thus
\[
\begin{aligned}
|\Delta\mathcal W_{\rm prime}|
&\le 2(x+1)\eta\sum_{n\le N_x}\frac{\Lambda(n)}{\sqrt n}\\
&\le 4(x+1)\eta\sqrt{N_x}\log N_x
=O(x^{3/2}\eta\log(ex)).
\end{aligned}
\tag{4.7}
\]
**Both poles.** For \(m_\eta=\int\rho_\eta(v)e^{v/2}dv=\int\rho_\eta(v)e^{-v/2}dv\), evenness gives
\(1\le m_\eta\le\cosh(\eta/2)\). Therefore
\[
\mathcal P(k_x^{\mathrm{sm}})=2m_\eta^2(x-1)\log x,\quad
|\Delta\mathcal P|\le 2(x-1)\log x\,\sinh^2(\eta/2)
=O(x\eta^2\log x).
\tag{4.8}
\]
Here the first inequality is exact; for \(0<\eta\le1/16\) the final bound uses \(\sinh(\eta/2)=O(\eta)\).

**Complete Gamma part.** Write \(d_\eta=k_x^{\mathrm{sm}}-k_x\) and \(w(t)=e^{-t/2}/(1-e^{-2t})\). Its change is
\[
-(\gamma+\log4\pi)d_\eta(0)+
2\int_0^\infty[d_\eta(0)e^{-t/2}-d_\eta(t)]w(t)\,dt.
\]
For \(0<t<2\eta\), the numerator is \(O(K_xt)\) by the common Lipschitz bound and \(1-e^{-t/2}\le t/2\). For \(2\eta\le t\le1\), it is \(O(K_x\eta)\) and \(w(t)=O(1/t)\). For \(t\ge1\), use \(w(t)=O(e^{-t/2})\). Therefore
\[
|\Delta\Gamma|=O(x\eta[1+\log(1/\eta)]).
\tag{4.9}
\]
All implicit constants in (4.7)–(4.9) are **absolute**, depending at most on the unit-mass/support conditions for the mollifier; they do not depend on \(x\).

Combining the three independent components,
\[
\boxed{
|\mathcal W(k_x^{\mathrm{sm}})-\mathcal W(k_x)|
=O\left(x^{3/2}\eta\log(ex)+x\eta^2\log x+
x\eta[1+\log(1/\eta)]\right).
}
\tag{4.10}
\]
In particular for \(\eta_x=x^{-2}\):
\[
\boxed{
|\mathcal W(k_x^{\mathrm{sm}})-\mathcal W(k_x)|
=O\left(\frac{\log(ex)}{\sqrt x}\right).
}
\tag{4.11}
\]
**This is uniform in \(x\ge4\).** Thus it is negligible compared with the arithmetic \(x^{1/2+\varepsilon}\) scale. It does not imply that the completed Connes trace operator admits these tests; their two pole moments remain nonzero.

### 2b. Smooth signed-Weil conditional lower bound

For the same \((x,q,d,n,m,k,T,T_2)\) satisfying (M), B, C, the source
\(f_x^{\mathrm{sm}}\in C_c^\infty\) obeys
\[
\boxed{
\begin{aligned}
\mathcal W(k_x^{\mathrm{sm}})\ge {}&
C_\zeta(x+1)+2\log(2\pi)-2D(x)\\
&-2\mathfrak b(x)-2C_\varepsilon(1+\log x)
-6x^{1/2+\varepsilon}
-O\left(\frac{\log(ex)}{\sqrt x}\right).
\end{aligned}}
\tag{4.12}
\]
The sign of the last \(O\) is shorthand for subtracting a nonnegative error bounded by an absolute constant times \(\log(ex)/\sqrt x\).

**What is genuinely new here:** the analytic map from Alpöge's assumed one-sided arithmetic inequality to a lower bound for a family of **actual smooth compact original signed Weil sources** can be made quantitatively lossless at square-root scale. This is a *conditional analytical bridge*, not an independent estimate for those sources.

## 3. Lean provenance: upstream source result is now identified

Frozen public Lean mirror, commit
\`7cbf7a0c32fca1badf05931879a37c3f8a2b03b4\`.

| Source path | Exact declaration | What it does; what it does not do |
|---|---|---|
| \`src/Formalism/StatementOC.lean\` | \`def Sz\`, \`def Gz\` | Confirms \(G_z=\sum\Lambda(n)(x/n-1)-x\log x+(1+\gamma)x\) at the **literal Lean definition**. |
| \`src/Formalism/SizeWindowBrackets/GrowthBoundGivesRiemannHypothesis.lean\` | \`def topExcess\`, \`GrowthInequalityDataExists\`, \`GrowthBoundAtPoint\` | Conditional mean-excess and size/character data; no direct Weil operator. |
| \`src/Formalism/SizeWindowBrackets/ChainSlopeWindows/LipschitzWeightAndMinimum.lean\` | \`theorem proposition_at_point\` | Converts packaged growth-data hypotheses to \(G_z\le x^{1/2+(C+3)\varepsilon}\). |
| \`src/Formalism/SizeWindowBrackets/SChain.lean\` | \`def ChainCore\`, \`ChainInputs\`, \`ChainTop\` | Records the **upstream quantitative arithmetic data**: slope, divisibility, size windows, losses, identification and cutoffs. |
| \`src/Formalism/SizeWindowBrackets/ChainSlopeWindows/LevelOfCongruence.lean\` | \`theorem uniform_bound_statement\` | Derives weighted mean bound under its 11 hypotheses. |
| \`src/Formalism/SizeWindowBrackets/ChainSlopeWindows/LipschitzWeightAndMinimum.lean\` | \`theorem top_from_weighted_mean\` | Identifies weighted-mean main terms with the prime table's \`topExcess\`. |
| \`src/Formalism/SizeWindowBrackets/ChainSlopeWindows/SubThreeMulLogTwoLe.lean\` | \`theorem joint_error_bound\` | Supplies the one-sided upper error under the source's premises. |
| \`src/Formalism/SizeWindowBrackets/ChainSlopeWindows/PChain.lean\` | \`theorem chain_top_bound\` | **Actual source proof** constructs \`topExcess ≤ x^(CErr*ε) * errorBracketSizeWindows\` by calling \`joint_error_bound top_from_weighted_mean ... (uniform_bound_statement ...)\` and packages the growth-data existence. |
| same \`PChain.lean\` | \`theorem chain_at_alt\`, \`chain_weighted_mean_to_unit_weight\` | Converts that package into a bound on \(G_z\) and the conditional RH-level implication. |

This pin separates the **one-sided table bound producer** (\`chain_top_bound\` through \`joint_error_bound\`) from the later **size/character analytic consumer** (\`proposition_at_point\`). It resolves the previous *file locator* gap, not the mathematical correspondence gap. No claim that the full dependency graph has been independently compiled or all imported published theorems reaudited.

In particular **the candidate "projection" cannot identify**
\(\texttt{WeightedMeanData.Data}\)'s arithmetic-lattice slope/inner-product datum with the original smooth Weil source metric \(\|f_x^{\mathrm{sm}}\|_2^2\), which is close to \(x-1\), **without a constructed map preserving both numerator and positive native denominator**. Existing results prove arithmetic inequalities feeding the Weil test, not equality of the arithmetic Hilbert spaces and the Weil Hilbert completion.

## 4. Noncircularity and falsifiers

1. **PASS analytic conditional:** (M)+B+C \(\Rightarrow\) (4.12), with no use of RH as a premise and no omitted pole/Gamma terms.
2. **PASS uniformity:** explicit \(\eta_x=x^{-2}\) controls the smooth-source replacement to \(O(x^{-1/2}\log x)\) in the *original scalar signed functional*. Does not establish trace-class operator intertwinement.
3. **OPEN arithmetic mechanism:** a new theorem deriving (M) from a simpler Weil/spectral principle, or a genuine map from theta-lift weighted slopes to Weil Gram pencils.
4. **FALSIFIER:** any proposed operator transfer must match the nonzero two-pole moments, all character and exceptional-prime corrections, the support radius \(\sim\frac12\log x\), and the quantitative baseline \(C_\zeta x\), not merely the sign \(\mathcal W\ge0\).

**Next — CR-ALP-005:** inspect the *input producer* \`uniform_bound_statement\` / \`joint_error_bound\` and extract the precise rank-trace inequalities and their positive metrics; then compare, entry by entry, with the log-\(L^2\) source Gram pencil on the smooth ramp family. If no invariant-preserving carrier map is constructible, close that projection route as **ANALYTIC REPRESENTATION ONLY** rather than an operator equivalence.

No modifications of Weil-Lab, original RH/provenance, or Hodge RH acceptance ledgers.
