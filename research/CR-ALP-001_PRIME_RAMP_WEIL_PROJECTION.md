# CR-ALP-001 — Alpöge prime-ramp to signed-Weil kernel: first projection audit

**Jurisdiction:** Canonical-Return, independent arithmetic-to-Weil research branch.
**Date:** 2026-11-22 PST / 2026-11-23 UTC.
**Standing:** SOURCE-ANCHORED; DERIVED / ANALYTIC; OPEN TRANSFER. **No new RH proof, global Weil positivity, Lean certification, or identification with theta lifts is claimed.**

## 1. Sources and precise question

**Question.** Can a specific intermediate arithmetic estimate from Alpöge's proof, rather than its already-certified RH conclusion, be transported noncircularly into a statement about the original signed Weil quadratic form?

- Alpöge, *A remark*, 20 November 2026, https://alpo.ge/remark.pdf: PDF p. 3, timeline items 4–5 (earlier determinant route to Weil-type positivity, then weighted von Mangoldt ramp); p. 4, items 6, 8–10, 14–16 (another positivity trial, theta-lift chain, divisibility, averaging). This is an author's retrospective timeline, **not** a proof that the two routes are equivalent.
- Lean/Mathlib joint statement, 21 November 2026, https://leanprover-community.github.io/blog/posts/riemann-hypothesis/: external certification of the submitted final RH proof, **not** of this new projection.
- Michael O. Rubinstein, *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, working draft v2, 22 November 2026, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper.pdf: §I.1 (I.2)–(I.4), the one-sided criterion; §I.5 (I.21), (I.25)–(I.26), and the main inequality (M). Rubinstein explicitly reports that the inequality (M) is **quoted** from the Lean development rather than proved in his draft.
- Canonical-Return T1B/T1C: exact selected pole-null Hermitian pencil and positive two-dimensional restrictions, stored at archive/space/09-t1b.md and archive/space/10-t1c.md. Neither gives global complement coercivity.

## 2. Original prime statistic and its Mellin kernel

For \(X>1\), define

\[
S(X)=\sum_{n\ge2}\Lambda(n)(X/n-1)_+
=\sum_{n\le X}\Lambda(n)(X/n-1),\qquad
G(X)=S(X)-X\log X+(1+\gamma)X.
\]

Here \(\gamma\) is Euler's constant. These are Rubinstein's (I.2).

For \(\operatorname{Re}s>1\),

\[
\int_n^\infty (X/n-1)X^{-s-1}\,dX
=\frac{n^{-s}}{s(s-1)};
\quad
\int_1^\infty S(X)X^{-s-1}\,dX
=-\frac{\zeta'(s)}{\zeta(s)s(s-1)}.
\]

**CANONICAL / RECOVERY:** This is the elementary Mellin identity already used in Rubinstein (I.3); not a new result.

## 3. Exact even prime-ramp test kernel

Put \(L=\log X\). Define the continuous, even, compactly supported function

\[
\boxed{k_X(y)=\big(Xe^{-|y|/2}-e^{|y|/2}\big)_+,\qquad y\in\mathbb R.}
\]

Its support is \([-L,L]\). For each integer \(n\ge1\),

\[
\boxed{\frac{k_X(\log n)}{\sqrt n}=(X/n-1)_+},
\quad\text{hence}\quad
\boxed{S(X)=\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}k_X(\log n).}
\]

Under the **even-test normalization** in which the prime part of the original signed Weil functional equals
\(-2\sum_{n\ge2}\Lambda(n)n^{-1/2}k(\log n)\),
the prime contribution at \(k_X\) is **exactly \(-2S(X)\)**. This statement attaches only the prime term: both poles and the entire Gamma contribution remain to be evaluated.

### 3a. Positive definite as a test kernel

With \(\widehat k(\tau)=\int_\mathbb R k(y)e^{-i\tau y}\,dy\), direct integration gives for real \(\tau\),

\[
\boxed{
\widehat{k_X}(\tau)
=\frac{X+1-2\sqrt X\cos(\tau\log X)}{\tau^2+1/4}
=\frac{|\sqrt X-e^{i\tau\log X}|^2}{\tau^2+1/4}\ge0.
}
\]

Thus \(k_X\) is positive definite in the *classical Fourier sense*. This is **not** positivity of the signed Weil quadratic form.

There is an independent real-space proof. Let
\(h(t)=Xe^{-t/2}-e^{t/2}\) on \(0\le t\le L\).
Then \(h(t)\ge0\), \(h''=h/4\), \(h(L)=0\), and \(-h'(L)=\sqrt X\). Twice integrating yields

\[
\boxed{
k_X(y)=\sqrt X(L-|y|)_+
 +\frac14\int_0^L h(t)(t-|y|)_+\,dt.
}
\]

Each tent \((t-|y|)_+\) is the autocorrelation of an indicator function of an interval of length \(t\). This positive mixture gives a second elementary proof of positive definiteness.

**DERIVED / ANALYTIC:** The normalization equality and the Fourier/tent factorizations are exact elementary identities; no novelty beyond this restricted observation is claimed. No Lean or external mathematical review has occurred.

### 3b. Smooth-source admissibility is NOT inherited

The kernel has cusps at \(0\) and at \(\pm L\), so it is **not** the smooth compact autocorrelation of the kind used for T1B's Connes trace attachment. Nor has it been shown to arise from a *single* smooth source satisfying the two Mellin pole-null constraints. A nonnegative mixture of nonsmooth autocorrelations does not automatically enter the original smooth trace domain. Precise smoothing, pole moments, trace class, and Gamma continuity are **OPEN**.

The support must also grow with \(X\). If a source \(f\) is supported in \((-a,a)\), its autocorrelation is supported in \((-2a,2a)\). Therefore any approximation preserving the entire ramp support requires \(a\) approximately at least \(\frac12\log X\); the fixed-radius T1C result cannot be reused uniformly as \(X\to\infty\).

## 4. Why this does not reproduce Alpöge's key inequality

Rubinstein's §I.1 reduces RH to a one-sided upper bound
\(G(X)\le X^{1/2+\varepsilon}\) on appropriately long, arbitrarily large intervals. The main arithmetic estimate (M) is an *averaged, character-weighted, place-dependent* bound, not a universal bound for arbitrary Weil tests.

The Fourier numerator of \(k_X\) contains an **\(X+1\) background** and a **\(\sqrt X\) oscillatory cross-term**. Positivity of the test kernel alone does not cancel the background or supply a square-root-quality estimate.

For diagnostic purposes **only**, if one assumes all zeros lie on the critical line, a zero with spectral parameter \(\tau\) is weighted by

\[
\widehat{k_X}(\tau)
=\frac{X+1-2\sqrt X\cos(\tau\log X)}{\tau^2+1/4}.
\]

This conditional spectral display is not an ingredient of a noncircular proof. The research obligation is to identify a legitimate **renormalization or quantitative bound for the oscillatory part** from an intermediate arithmetic theorem, retaining all pole/Gamma terms and support quantifiers.

## 5. Falsifiable tests / next gate

1. **Full-form completion:** calculate both original pole contributions and full Gamma term on justified smooth approximants to \(k_X\), and compare their leading terms to the normalization of \(G(X)\). An uncontrolled \(O(X)\) background blocks the desired square-root projection.
2. **Trace-domain test:** determine whether there is a legitimate smooth, two-Mellin-pole-null source for the required Connes trace application, or report **NO TRANSFER**.
3. **Lean source pin:** isolate a genuinely intermediate inequality corresponding to Rubinstein's (M), with precise auxiliary characters, averaging, scales, existential quantifiers, and dependencies. Do not use the final RH theorem as a premise.
4. **Uniform-support test:** track the expanding \(a\gtrsim \tfrac12\log X\) and quantify losses from complementary coupling.
5. **Counterexample test:** seek sources with identical ramp prime data but different pole/Gamma contributions; failure would refute any prime-only full-form dictionary.

## 6. Claim ledger

| Statement | Status |
|---|---|
| Alpöge announcement and joint Lean verification | CANONICAL external evidence |
| Rubinstein \(S,G\), Mellin kernel, criterion | CANONICAL / RECOVERY |
| \(k_X(\log n)/\sqrt n=(X/n-1)_+\) | DERIVED / ANALYTIC |
| Fourier factorization / positive tent mixture | DERIVED / ANALYTIC; symbolic check only |
| Full signed Weil form at \(k_X\) | OPEN, nonsmooth-source qualification |
| Two-pole-null semilocal trace transfer | OPEN |
| Square-root upper bound from signed Weil form | OPEN |
| Alpöge theta-lift mechanism equivalent to CR construction | NOT CLAIMED |
| RH/global Weil positivity/Lean certification from CR | NOT CLAIMED |

**Next:** CR-ALP-002, a full pole-and-Gamma normalization and regularization audit, coupled with an exact intermediate-Lean-inequality inventory. The current work is a **new research branch**, not a modification of Weil-Lab or any earlier accepted proof ledger.
