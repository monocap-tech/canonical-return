# Canonical Return — Phase T1D research index

**Scope:** original signed Weil form on selected, two-Mellin-pole-null, four-cell source families at support radius \(a=2\). These notes are **analytical research** with independent finite symbolic/rational checks, **not Lean-certified, not complete operator-domain comparisons, not global Weil positivity, and not RH**. The public repository's `main` branch remains unchanged until independent review of [draft PR #2](https://github.com/monocap-tech/canonical-return/pull/2).

| Research note | Claim with exact support restriction | Next obstruction |
|---|---|---|
| [T1D-A](CR-INT-T1D-A.md) | Every real nonzero smooth profile of width ≤1/16 has positive **selected 2D** Weil form, lower margin 9/20; first extra prime 29 | Further prime collisions |
| [T1D-B](CR-INT-T1D-B.md) | Every real profile of width <log(9/8) retains **selected 2D** lower margin 2/25 through prime 25 | First spacing-two collision (prime power 8) |
| [T1D-C](CR-INT-T1D-C.md) | Every real profile of width <log(31/27) has **selected 2D** lower margin 1/40, including prime power 8 | Prime 31; infinite complement |
| [T1D-D](CR-INT-T1D-D.md) | Exact prime-31 mixed-only Gram perturbation; **nonnegative** profiles of width <log(11/9) retain generalized phase splitting >1/16 (no wider minimum-sign claim). Exact *conditional* source-core Schur criterion for the infinite complement. | Actual bound on mixed (q(u_\pm,y)) and (q(y,y)), pole carrier, completed-domain dictionary |

The first five arithmetic activations beyond the base 3,9,27 are 29, 25, 8, 31, 23, 32; the next spacing-two event is 11. All thresholds are *possible*, profile-dependent interactions; the full signed form includes Gamma and both Mellin poles, even though the two poles vanish on this particular selected source.

## Audits

Run `bash tools/check_all.sh` from the repository root on this PR branch (Python 3 + SymPy). This executes the T1B/T1C regressions and T1D-A/B/C/D finite algebra and interval tests. **No GitHub Actions workflow is installed.**

## What an actual cross-canonical theorem must add

The finite-source subspace (E_\phi) has (mu_\pm>0) only in the certified support ranges. On the pole-null smooth core decompose (mathcal D_a^0=E_\phi\oplus F_\phi) orthogonally in the **native (L^2) metric**. Conditional on (mu_\pm>0), positivity of the full pole-null core requires and is equivalent to:

\[
q(y,y)\ge \frac{|q(u_+,y)|^2}{\mu_+}+\frac{|q(u_-,y)|^2}{\mu_-}\qquad\forall y\in F_\phi.
\]

**This inequality is not proved.** Moreover (E_\phi) has both pole moments zero, so it does not contain the two pole directions in the finite conditioning space of Andrews–Connes–La Gatta §8.6; an operator-domain and conditioning dictionary is independently required. Neither a positive finite compression nor its phase gap controls the infinite-dimensional spectral bottom.

## T1D-E/F/G — mixed-source and sharp-rate frontier

| Stage | Research note | Proven source-level statement | Unproved boundary |
|---|---|---|---|
| **E** | [T1D-E](CR-INT-T1D-E.md) | Explicit full signed Weil *one-sided* mixed L² representer and a positive rank-at-most-two Schur penalty | Actual infinite-dimensional complement lower bound |
| **F** | [T1D-F](CR-INT-T1D-F.md) | For the fixed even shape, the two projected residuals have opposite parity, with zero mixed Gram entry; the *penalty* norm decays as O(1/log(1/ell)) | Certified finite-width residual Gram values and coercivity |
| **G** | [T1D-G](CR-INT-T1D-G.md) | Exact **nonzero limits of the residual-square norms**, giving the *sharp* rate \(\log(1/\ell)\|\mathcal R_\ell\|\to V_{\psi}+V_{3,+}+C_{+}>18/5\). A fixed low-frequency-energy smooth test persists inside all selected narrow-source complements. | Certified \(V_\psi\), actual width-10⁻⁴ residual norms and positive complement coercivity |

The T1D-G exact rational/SymPy script is `src/audits/canonical_return_t1d_g_sharp_rate.py` and is included in `tools/check_all.sh` along with earlier checks. The Fourier rescaling/concentration argument is documented in the research note; **running a finite algebra audit does not certify that analytic proof in Lean**.

**Next — T1D-H:** perform certified integral enclosures for the *actual frozen-width* projected \(\beta_\pm^2\) and \(V_\psi\), then independently test a coercive estimate for the *actual* pole-null complementary Weil form (or record the exact mismatch with Connes's semilocal conditioned Sonin carrier).
