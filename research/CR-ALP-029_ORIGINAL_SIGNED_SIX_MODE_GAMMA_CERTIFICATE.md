# CR-ALP-029 — Exact positive-kernel Gamma expansion certifies an original signed Weil 6×6 block

**2026-11-23 PST — Canonical-Return. Research-only; DRAFT/NONPROMOTING.** Original fixed-collar two-pole-null signed Weil form; no full 15,000-mode Schur certificate, receiving Lean kernel theorem, support extension, or RH. This stage goes beyond the uncertified CR-ALP-028 six-mode diagnostic by producing a finite, independently reproducible **outward-integer/rational positive-definiteness certificate** on the six *genuine original signed Weil* directions `w_3,...,w_8`.

## 0. Original source custody

The original signed explicit formula is `Zeta23.EF.literatureRHS`, `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, Git blob `0034976fe2f1e5fa9d411108736160542eb94293`. Its original paired prime cutoff is `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`. The original form consists of **both original poles minus seven paired Λ(n)/sqrt(n) prime-power translations plus full Gamma bracket**. The pole terms vanish exactly on the two-pole-null physical sources used here. The seven prime powers are exactly `{2,3,4,5,7,8,9}`; no positive surrogate is inserted.

Retain the original physical interval `I=(-a,a)`, `a=11/10`, the zero-extended normalized sine functions `e_j(x)=a^{-1/2} sin(ν_j(x+a))`, `ν_j=jπ/(2a)`, and the **same exact pole-null interior combinations** from CR-ALP-028:

\[
 w_j=e_j-\lambda_j e_{k(j)},\qquad
 k(j)=\begin{cases}1&j\text{ odd},\\2&j\text{ even},\end{cases}\quad
 \lambda_j=\frac{\nu_j(\nu_{k(j)}^2+1/4)}{\nu_{k(j)}(\nu_j^2+1/4)}
 \quad(j=3,\ldots,8).
\]

They are smooth in the open collar and vanish at both boundaries; zero extension belongs to `H¹(ℝ)`, hence to the logarithmic Fourier form domain. The physical metric is **not** replaced: `M_{ij}=δ_{ij}+1_{i≡j(2)}λ_iλ_j`. The odd/even parity blocks decouple in the complete original signed form.

## 1. A convergent positive-matrix expansion of the entire original Gamma symbol

The classical digamma integral used already in CR-ALP-028 gives the exact energy identity (under the paper Fourier convention):

\[
 g(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi,\qquad
 W_\Gamma(f,f)=g(0)\|f\|_2^2+\int_0^\infty \frac{e^{-s/2}}{1-e^{-2s}}\|f-\tau_sf\|_2^2\,ds.
\]

Here `g(0)=-γ-π/2-3log(2)-log π`. Expand the strictly positive scalar kernel:

\[
 \frac{e^{-s/2}}{1-e^{-2s}}=\sum_{k=0}^\infty e^{-\beta_k s},\qquad \beta_k=2k+\frac12.
\]

Tonelli for the diagonal positive energies and polarization for mixed entries give **without any numerical quadrature**

\[
 \boxed{W_\Gamma(f,f)=g(0)\|f\|_2^2+\sum_{k=0}^\infty J_k(f,f),\quad
 J_k(f,f)=\int_0^\infty e^{-\beta_k s}\|f-\tau_sf\|_2^2ds\ge0.}\tag{29.1}
\]

Every `J_k` is positive semidefinite **as a Hermitian form on any physical source subspace**. The integral series is finite for these H¹ zero-extended sources (translation Lipschitz estimate at small s; exponential decay at large s). This is the analytic reason a partial Gamma sum is a **Loewner lower bound** on the *complete* Gamma term, even though `g(0)<0` and the *original prime part has its original negative sign*.

## 2. A closed sine-Laplace formula for each positive kernel

Write `S_s f(x)=f(x-s)` and `U_{ij}(s)=⟨e_i,S_s e_j⟩+⟨e_i,S_{-s}e_j⟩`; `U` is given in CR-ALP-028 (28.5), and is zero for `s≥2a`. Elementary Laplace integration uses `ν_j·2a=jπ` and gives, for every β>0,

\[
 \int_0^\infty e^{-\beta s}U_{ij}(s)ds
 =\mathbf1_{i=j}\frac{2\beta}{\beta^2+\nu_i^2}
 +\mathbf1_{i\equiv j\ (2)}
 \frac{2\nu_i\nu_j(1-(-1)^ie^{-2a\beta})}{a(\beta^2+\nu_i^2)(\beta^2+\nu_j^2)}.
\]

Since `⟨e_i-S_s e_i,e_j-S_s e_j⟩=2δ_{ij}-U_{ij}(s)` (appropriate Hermitian polarization), the exact matrix for `J_β` on the sine basis is

\[
 \boxed{
 J_{\beta,ij}=
 \mathbf1_{i=j}\frac{2\nu_i^2}{\beta(\beta^2+\nu_i^2)}
 -\mathbf1_{i\equiv j\ (2)}
 \frac{2\nu_i\nu_j(1-(-1)^ie^{-2a\beta})}{a(\beta^2+\nu_i^2)(\beta^2+\nu_j^2)}.
 }\tag{29.2}
\]

The matrix is the exact Gram of a translation difference integrated against a positive weight. It **need not appear entrywise positive**; its quadratic form is PSD. The numerical regression separately compared (29.2) to direct correlation Laplace quadrature for `(i,j,β)=(3,3,1/2),(3,5,1/2),(4,6,5/2),(3,4,1/2),(7,7,41/2)`; all differences were below `10^-11`. Those spot checks are diagnostics; (29.2) is a closed analytic identity.

Map `J_{β}` to the six w-directions by the exact sine coefficients `e_j-λ_j e_{k(j)}`. The original signed prime matrix `P` comes from the separately reexecuted CR-ALP-028 **outward-rational** seven-prime audit (not from floating point).

## 3. An EXACT finite positivity criterion with no Gamma-tail estimate

Let `H_K` be the genuine **lower matrix**

\[
 \boxed{H_K=g(0)M+P_{\rm original}+\sum_{k=0}^{K-1}J_{\beta_k}^{(w)},\quad
 H_{\rm full}=H_K+\sum_{k\ge K}J_{\beta_k}^{(w)}\succeq H_K.}\tag{29.3}
\]

This preserves the full source-sign prime contributions and the actual physical metric. **If** the real-symmetric (complex Hermitian) finite matrix `H_K` is positive definite, so is the *complete original signed Weil six-mode matrix*, without any numerical upper bound on the uncomputed Gamma tail.

The executed audit takes `K=30,000`. It constructs a directed, exact-integer interval implementation with grid denominator `10^42`, starting with 38-digit rigorous enclosures from the original prime auditor. No normal floating-point arithmetic determines a sign. It rigorously encloses `π` via Machin arctan, `log` via atanh with explicit remainder, `sqrt(n)` by integer `isqrt`, `exp(a)` via rational Taylor + geometric tail, and the Euler constant `γ` via the classical alternating Euler–Maclaurin harmonic inequality

\[
 D_6 < H_n-\log n-\gamma < D_4,\quad
 D_4=\frac1{2n}-\frac1{12n^2}+\frac1{120n^4},\quad
 D_6=D_4-\frac1{252n^6},\quad n=256.
\]

This classical Euler–Maclaurin error-sign assertion, like the digamma identity, is **human mathematical input**, not a theorem newly compiled by Lean. All arithmetic consequences are checked with directed integer/rational enclosures.

The exact interval `LDLᵀ` pivot lower bounds (listed below as **rounded-down displays**, while the executable works with exact rational endpoints) are:

| Original six-mode parity block | Pivot 1 | Pivot 2 | Pivot 3 |
|---|---:|---:|---:|
| `w_3,w_5,w_7` | `>0.0035130` | `>0.000034706` | `>0.0000018372` |
| `w_4,w_6,w_8` | `>0.0014010` | `>0.000027311` | `>0.000031161` |

The interval LDL pivots are all rigorously positive. **Thus `H_30000 ≻ 0`, and (29.3) gives `H_full ≻ 0`** for the whole *six-mode original signed Weil* carrier.

In fact the audit certifies the stronger quantitative *physical* generalized eigenvalue inequality:

\[
 \boxed{\forall f\in\operatorname{span}_\mathbb{C}\{w_3,w_4,w_5,w_6,w_7,w_8\},\quad
 W_{\rm original}(f,f)>5\cdot10^{-9}\|f\|_{L^2(I)}^2\quad(f\ne0).}\tag{29.4}
\]

For each 3×3 parity block with `H_30000=L D Lᵀ`, the exact proof bounds `||L^{-1}||_F²` and `tr(M)` by directed rational intervals and uses

\[
 H_{30000}\succeq\frac{\min_iD_{ii}}{\|L^{-1}\|_F^2}I,
 \qquad 0\prec M\preceq\operatorname{tr}(M)I.
\]

The exact rational lower quotient is `>5.5822×10^-9`, so the displayed `5×10^-9` is a deliberately conservative certified constant. This is not merely reproducing the previous 45-digit unvalidated 6×6 eigenvalue diagnostic `≈2.2579×10^-8`.

## 4. Reproduction, scope and gate ledger

**Executable:** `src/audits/cr_alp_029_exact_gamma_series_6x6_audit.py` (requires the same-directory `cr_alp_028_prime_interval_and_gram_audit.py` and the Python standard library). Command: `python src/audits/cr_alp_029_exact_gamma_series_6x6_audit.py 30000`. The complete exact lower-pivot fractions and rational eigenvalue lower quotient are emitted as `CR-ALP-029_certified_six_mode_result.json` in the script directory. This audit was executed successfully in the present Python runtime. No `lean`/`lake` kernel output is claimed.

**PASS**: literal native prime signs, exact seven prime powers, independently source-certified prime 6×6 component, original physical pole-null basis and metric, exact Gamma kernel finite lower matrix, certified 30,000-term directed-integer arithmetic, all six positive pivots, and physical generalized eigenvalue `>5×10^-9`.

**ANALYTIC, NOT LEAN**: source-to-Gamma identity; positivity and convergence of the infinite Gamma translation-difference series; sine Laplace formula; Euler–Maclaurin error bound; full-form Loewner comparison. These are not discharged by the executable itself.

**OPEN**: the additional `u_1,u_2` physical boundary anchors; any other modes `w_9,...,w_15000`; all mixed low-low and low-high original signed Weil couplings; the true `15,000` finite signed matrix and cross-tail Schur residual; full fixed-collar positivity; arbitrary-support extension; Connes/Suzuki transfer; F41 physical H/V completeness; formal Lean acceptance; RH. Positivity of a six-mode principal subblock **does not** prove positivity of any of these larger domains. No receiving Weil-Lab, F41/F42, Formalization, Hodge, provenance, or RH ledger modification.

**Next CR-ALP-030**: exploit the Gamma PSD-series monotone lower matrix and the exact physical parity/Gram rank-one structure to build *block-certified* matrix-free lower bounds on progressively larger original low sections; detect near-zero modes and verify the two exceptional anchors rather than assuming six-mode positivity transfers. In parallel, produce certified **upper** residual norms for the complete 15k Schur gate. Draft remains nonpromoting.
