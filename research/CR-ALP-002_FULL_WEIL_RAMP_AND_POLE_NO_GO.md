# CR-ALP-002 — Exact full signed Weil ramp, pole moments, and spectral baseline

**Jurisdiction:** independent Canonical-Return / Alpöge projection.
**Opened:** 2026-11-22 PST / 2026-11-23 UTC.
**Standing:** DERIVED / ANALYTIC, with executable high-precision regression PASS. **No Lean certification, completed Suzuki/Connes operator-domain theorem, new RH theorem, or equivalence with Alpöge's theta lifts.**

## 1. Pinned sources and comparison

- Alpöge, *A remark* (20 November 2026), https://alpo.ge/remark.pdf, printed p. 3, items 4–5: the older Weil/determinant experiment and the October 20 switch to the weighted von Mangoldt statistic. Historical motivation only; not a transfer theorem.
- Rubinstein and collaborators, *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, **version 2, working draft 22 November 2026**, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, §I.1 equations (I.2)–(I.4), one-sided interval lemma, and §I.5 (M). The paper derives the analytic implication but still quotes (M) from Lean. It does not itself prove the main arithmetic estimate.
- Canonical-Return original full signed Weil Gamma convention: archive/space/10-t1c.md, §3 (complete archimedean diagonal), with the prime/pole conventions in archive/space/09-t1b.md. These are local research records, not authority for an extension to operator completions.

All computations in this note use **the complete original signed form** with two poles, negative von Mangoldt prime term and full archimedean Gamma term. No zero-location assumption is required for the identities below.

## 2. Exact **single-source** factorization found

Write \(X>1\), \(L=\log X\), and set
\[
f_X(u)=e^{u/2}\mathbf 1_{[0,L]}(u),\qquad
\widetilde f_X(u)=\overline{f_X(-u)}.
\]

Then \(f_X\) is a compactly supported \(L^2\) function, but **not** smooth. Its autocorrelation is precisely the CR-ALP-001 ramp:
\[
\boxed{
(f_X*\widetilde f_X)(y)
=\left(Xe^{-|y|/2}-e^{|y|/2}\right)_+
=k_X(y).}
\]

**Proof.** For \(0\le y\le L\),
\[
\int_0^{L-y}e^{(u+y)/2}e^{u/2}\,du
=e^{y/2}(e^{L-y}-1)=Xe^{-y/2}-e^{y/2};
\]
the correlation is even and vanishes outside \([-L,L]\). Consequently,
\[
\|f_X\|_2^2=k_X(0)=X-1,\quad
\widehat k_X(t)=|\widehat f_X(t)|^2
=\frac{|\sqrt X-e^{itL}|^2}{t^2+\frac14}\ge0.
\]

This is a **strict improvement over the earlier positive tent-mixture argument**: the logarithmic prime ramp is an explicit single autocorrelation, not just a positive definite mixture.

### 2a. Exact obstruction to the existing CR pole-null trace lane

The two Mellin pole moments of the source are
\[
F_+(f_X)=\int f_X(u)e^{u/2}du=X-1,
\qquad
F_-(f_X)=\int f_X(u)e^{-u/2}du=L.
\]
Neither vanishes. The even correlation satisfies
\[
\boxed{
\int_{\mathbb R}k_X(y)e^{y/2}dy
=\int_{\mathbb R}k_X(y)e^{-y/2}dy
=(X-1)\log X>0.}
\]

For a source \(g\) with both Mellin pole moments zero, its autocorrelation has zero corresponding exponential moment. Therefore **no sequence** of such pole-null autocorrelations can converge to \(k_X\) in the weighted \(L^1(e^{|y|/2}dy)\) topology: the continuous moment functional would give both 0 and \((X-1)L>0\).

This is an actual **NO-TRANSFER** result for the naive identification with CR-T1B's two-pole-null Connes trace lane, not merely an unverified domain issue. Other, explicitly corrected or non-pole-null carrier mappings remain open.

## 3. Full signed Weil distribution on the ramp

Use the even logarithmic test normalization
\[
\mathcal W(k)=\mathcal P(k)-2\sum_{n\ge2}
\frac{\Lambda(n)}{\sqrt n}k(\log n)+\Gamma(k),
\]
where
\[
\mathcal P(k)=\int_{\mathbb R}
k(y)(e^{y/2}+e^{-y/2})\,dy
\]
includes **both zeta poles**, and the complete Gamma contribution is
\[
\Gamma(k)=-(\gamma+\log 4\pi)\,k(0)
+2\int_0^\infty[k(0)e^{-t/2}-k(t)]
\frac{e^{-t/2}}{1-e^{-2t}}\,dt .
\]
This is the convention used in T1C. All these expressions converge for the explicit continuous, piecewise-smooth ramp \(k_X\); this is **a distributional quadratic-form value**, not yet admission of nonsmooth \(f_X\) to Suzuki's differentiated operator or Connes's semilocal trace class.

The prime term is exactly \(-2S(X)\), where
\[
S(X)=\sum_{n\le X}\Lambda(n)(X/n-1).
\]
The pole term is
\[
\boxed{\mathcal P(k_X)=2(X-1)\log X.}
\]

### 3a. Direct Gamma integration (no zero sum and no RH)

Since for \(0<t<L\),
\([k_X(0)e^{-t/2}-k_X(t)]=e^{t/2}-e^{-t/2}\),
and for \(t>L\) it is \((X-1)e^{-t/2}\), the complete integral splits into
\[
2\int_0^L \frac{dt}{1+e^{-t}}
+2(X-1)\int_L^\infty\frac{e^{-t}}{1-e^{-2t}}\,dt.
\]
Both pieces are elementary:
\[
\boxed{
\Gamma(k_X)=-(\gamma+\log4\pi)(X-1)
+2\log\frac{X+1}{2}
+2(X-1)\operatorname{atanh}(X^{-1}).}
\]
No archimedean cross term is dropped.

Thus the **full signed Weil distribution** at the explicit ramp is the exact prime-side identity
\[
\boxed{\mathcal W(k_X)=
2(X-1)\log X-2S(X)+\Gamma(k_X).}
\]
This is a complete pole–prime–Gamma attachment at the **distribution level**. It neither promotes the source to the pole-null trace domain nor proves a positivity bound.

## 4. Exact arithmetic–spectral renormalization

Following Rubinstein (I.2), let
\[
G(X)=S(X)-X\log X+(1+\gamma)X.
\]
Define constants and a convergent correction
\[
C_\zeta=2+\gamma-\log(4\pi)
=\sum_{\rho}\frac1{\rho(1-\rho)}
\approx0.04619141793>0,
\]
\[
D(X)=\sum_{m\ge1}\frac{X^{-2m}}{(2m)(2m+1)}
=1-X\operatorname{atanh}(X^{-1})
-\frac12\log(1-X^{-2})>0.
\]
The sum over \(\rho\) counts nontrivial zeros with multiplicities; its evaluation as \(C_\zeta\) is a standard consequence of the completed zeta product and **does not require RH**.

Substituting the above Gamma formula and simplifying gives the central exact identity
\[
\boxed{
\mathcal W(k_X)
=C_\zeta(X+1)-2G(X)+2\log(2\pi)-2D(X).
}
\tag{CR-ALP-002}
\]
It is derived directly from the original pole–prime–Gamma side; no assumption about zeros is used.

**Independent zero-side check (classical explicit formula).** At a nontrivial zero \(\rho\), evaluate the entire Fourier transform at \(t=-i(\rho-\frac12)\):
\[
\widehat k_X\left(-i(\rho-\frac12)\right)
=\frac{X+1-X^\rho-X^{1-\rho}}{\rho(1-\rho)}.
\]
The multiset symmetry \(\rho\leftrightarrow1-\rho\) and Rubinstein (I.4), including the exact trivial-zero tail \(-D(X)\), reproduce the same equation after summing. This is a **consistency identity**, not a replacement for the prime-side proof. The zeros need not lie on the critical line.

## 5. The *quantitative* projection of Alpöge's one-sided criterion

From the exact identity, for any \(0<\theta<1\),
\[
\boxed{
G(X)\le X^\theta
\quad\Longleftrightarrow\quad
\mathcal W(k_X)\ge
C_\zeta(X+1)+2\log(2\pi)-2D(X)-2X^\theta.
}
\]
Rubinstein's §I.1 proves that this one-sided estimate, for \(\theta=\tfrac12+\varepsilon\) over arbitrarily large intervals \([Y,Y^{3/2}]\), excludes off-critical-line zeros. Here the right side is the **required sharp renormalized lower bound**, not mere nonnegativity.

**Important quantitative obstruction:** \(\mathcal W(k_X)\ge0\) would imply only
\[
G(X)\le \tfrac12 C_\zeta(X+1)+\log(2\pi)-D(X),
\]
an \(O(X)\) upper bound, far too weak. By contrast, the one-sided criterion requires control of \(\mathcal W(k_X)\) at its \(C_\zeta X\) baseline with error of square-root quality. This separates **positivity** from **near-baseline rigidity** on the particular ramp family. The family norm itself grows as \(\|f_X\|_2^2=X-1\).

This gives a precise, falsifiable arithmetic-to-Weil **projection of the analytical criterion**, **not** a derivation of its required lower bound from the theta-lift construction. Repackaging \(G(X)\) with an exact equality does not strengthen the underlying estimate.

## 6. What is actually blocked and what remains

**PASS / DERIVED:** explicit nonsmooth single source; pole moments; exact two-pole, prime and Gamma terms; algebraic normalized identity; finite-test quadrature/audit regressions.

**NO TRANSFER:** CR-T1B pole-null trace lane cannot admit this ramp directly, nor through weighted \(L^1\)-convergent pole-null autocorrelations. The obstruction is the nonzero exponential pole moment, not an incidental lack of smoothness.

**OPEN:** rigorous form-domain approximation as required by each specific Suzuki/Connes operator; a corrected pole-preserving/non-pole-null trace representation; a noncircular derivation of the square-root-quality baseline rigidity from an *intermediate* Alpöge theorem; the exact correspondence of Rubinstein's **averaged, character-weighted** (M) with any candidate Weil form. Version 2 of Rubinstein leaves (M) dependent on unexpanded Lean inputs (§III.11).

**No claim:** RH, independent whole-Weil positivity, spectral operator equivalence, improvement on Alpöge's verified proof, or Lean theorem promotion.

## 7. Next — CR-ALP-003

A different transfer mechanism is required, because the direct pole-null path is blocked.

1. Pin the **intermediate Lean** inequality behind Rubinstein's (M) with every character average, auxiliary prime, scale and existential quantifier intact.
2. Study *pole-preserving* rather than pole-null operator attachments, using smooth approximants to \(f_X=e^{u/2}1_{[0,\log X]}\) and explicitly retaining their moment corrections.
3. Test whether the Hilbert-space maps in Rubinstein §III.8 (moment/rank/trace data) can be compared to CR's Gram-pencil invariants without conflating their distinct physical metrics.
4. Search for a new independent proof of the **renormalized baseline inequality**. If none is found, classify the arithmetic–Weil identity as a useful exact reexpression, not a simpler proof.

Source-level research only; no action on Weil-Lab, original RH/provenance, or Hodge ledgers.
