# CR-ALP-030 — A complete original signed Weil eight-mode certificate, including both boundary anchors

**Canonical-Return; November 23, 2026 (America/Los_Angeles). RESEARCH / DRAFT / NONPROMOTING.** Stage CR-ALP-030 extends CR-ALP-029's certified six interior modes by the **two genuine pole-null physical boundary-anchor directions**. It preserves the entire original signed Weil explicit-formula expression. The new executable verifies an eight-dimensional finite **form-domain** inequality, *not* all fixed-collar Schwartz sources or the 15,000-mode finite Schur gate. The analytic identities below are human mathematical inputs, **not** newly kernel-certified Lean results.

## 0. Original source and mathematical scope

- Receiving original EF: `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, Git blob `0034976fe2f1e5fa9d411108736160542eb94293`, literal `Zeta23.EF.literatureRHS`. Its terms are precisely two poles, **minus** the two-sign prime-power contributions, and the **full** `gammaBracket(r)=Re digamma(1/4+ir/2)-log(pi)`.
- Exact original IV011 finite cutoff source: `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`.
- Original fixed interval: `I=(-a,a)`, `a=11/10`. The seven prime powers are **exactly** `2,3,4,5,7,8,9`, with weights `Lambda(n)/sqrt(n)`, retaining all repeated von Mangoldt weights.
- On the two-pole-null logarithmic Fourier form-domain closure, original pole residues vanish. Boundary anchors need not be compactly supported smooth after zero extension; the result concerns the **closed original signed quadratic-form extension**, not a false claim that `literatureRHS`'s original compactly smooth input condition directly holds for a discontinuous endpoint representative.

## 1. Exact physical eight-dimensional section

Let the zero-extended physical normalized Dirichlet sine functions be

\[
e_j(x)=a^{-1/2}\sin(\nu_j(x+a)),\qquad \nu_j=\frac{j\pi}{2a}.
\]

Write `b_1(x)=cosh(x/2)`, `b_2(x)=sinh(x/2)` on `I`, and let `Pi_0` be the physical `L²(I)` orthogonal projection onto the common kernel of the two original pole functionals. Define

\[
\lambda_j=\frac{\nu_j(\nu_{k(j)}^2+1/4)}{\nu_{k(j)}(\nu_j^2+1/4)},\quad
k(j)=\begin{cases}1&j\text{ odd},\\2&j\text{ even},\end{cases}
\quad w_j=e_j-\lambda_je_{k(j)}\quad (3\le j\le8),
\]

\[
B_1=a+\sinh a,\quad B_2=\sinh a-a,\quad
c_k=\langle e_k,b_k\rangle,\quad \gamma_k=c_k/B_k,\quad
u\!\operatorname{anchor}_k:=u_k=e_k-\gamma_k b_k=\Pi_0e_k.
\]

The physical source basis, ordered by parity, is

\[
\mathscr E_8
=\underbrace{\operatorname{span}_{\mathbb C}(w_3,w_5,w_7,u_1)}_{\text{even physical functions}}
\oplus
\underbrace{\operatorname{span}_{\mathbb C}(w_4,w_6,w_8,u_2)}_{\text{odd physical functions}}.
\tag{30.1}
\]

All eight directions are **exactly two-pole-null**. The first six vanish at both endpoints. The two `u_k` are the genuine boundary-discontinuous (on zero extension) physical anchors, included without substituting a different space. The exact physical mass blocks are

\[
\boxed{M_k=
\begin{pmatrix}
I_3+\lambda\lambda^T&-\lambda\\
-\lambda^T&1-c_k^2/B_k
\end{pmatrix}\succ0.}\tag{30.2}
\]

The auditor independently proves the physical Gram intervals are positive definite. This is a **true physical L² metric**, not the frozen F41 `H/V` matrices and not an invented orthonormal basis.

## 2. Source-signed seven-prime entries, including anchor interactions

Write

\[
U_{fg}(t)=\langle f,\tau_tg\rangle+\langle f,\tau_{-t}g\rangle,\qquad
P_{fg}=-\sum_{n\in\{2,3,4,5,7,8,9\}}
\frac{\Lambda(n)}{\sqrt n}U_{fg}(\log n).
\tag{30.3}
\]

The interior sine-sine pairings come from CR-ALP-028's exact sine translation formulas, checked by its independent original prime interval audit. All new anchor-to-interior and anchor-to-anchor terms use the *same original* negative signed seven-prime pairing, with `u_k=e_k-γ_k b_k` and exact elementary overlap integrals.

For example, with `b_p(x)=exp(px)` on `I`, `p=±1/2`, and `t=log n`, the exact paired translation is calculated by

\[
\begin{aligned}
C_t(e_j,b_p)&=
\frac{\nu_j e^{-pa}\cos(\nu_jt)-p e^{-pa}\sin(\nu_jt)
-(-1)^j\nu_j e^{p(a-t)}}{\sqrt a(\nu_j^2+p^2)},\\
C_t(b_p,e_j)&=
\frac{(-1)^je^{pa}\bigl[-p\sin(\nu_jt)-\nu_j\cos(\nu_jt)\bigr]
+\nu_j e^{p(t-a)}}{\sqrt a(\nu_j^2+p^2)},
\end{aligned}\tag{30.4}
\]

where `C_t(f,g)=<f,τ_t g>`, so `U_{e_j,b_p}=C_t(e_j,b_p)+C_t(b_p,e_j)`. Likewise, for `p+r ≠ 0`,

\[
U_{b_p,b_r}(t)=
(e^{-pt}+e^{-rt})\frac{e^{(p+r)a}-e^{(p+r)(-a+t)}}{p+r},
\tag{30.5}
\]

with the continuous `p+r=0` limit `(e^{-pt}+e^{-rt})(2a-t)`. Even/odd hyperbolic combinations then yield all required pairings. The inputs `e^{±a/2}`, `sqrt(n)`, `log n`, `π`, and trigonometric evaluations are enclosed by the original outward-rational machinery; *no floating sign check* enters the prime block.

## 3. Complete original Gamma energy: certified finite Loewner lower matrix

The classical digamma series gives, for supported functions in the logarithmic form domain,

\[
W_\Gamma(f,f)=g(0)\|f\|_2^2+\sum_{m=0}^{\infty}J_m(f,f),
\quad g(0)=-\gamma-\pi/2-3\log 2-\log\pi,
\]

\[
J_m(f,f)=\int_0^\infty e^{-\beta_m s}\|f-\tau_sf\|_2^2\,ds\ge0,
\qquad \beta_m=2m+\frac12.
\tag{30.6}
\]

Every `J_m` is a **positive-semidefinite Hermitian form**, including on boundary anchors. Consequently the *full* authentic signed matrix is **Loewner above**

\[
\boxed{H_K=g(0)M_8+P_{\mathrm{original}}+
\sum_{m=0}^{K-1}J_m\,,\quad H_{\rm full}\succeq H_K.}\tag{30.7}
\]

The interval evaluator uses `K=30000`. The `J_m` matrices are elementary because the Laplace-convolution kernel is `exp(-β|x-y|)`, and the primitive carrier consists only of sines and `exp(±x/2)`. The sine-sine entry is the CR-ALP-029 rational-exponential formula. The new sine-exponential term follows by integrating

\[
(\mathcal K_\beta e_j)(x)=
\frac{2\beta e_j(x)}{\beta^2+\nu_j^2}+
\frac{\nu_j}{\sqrt a(\beta^2+\nu_j^2)}
\left(e^{-\beta(x+a)}-(-1)^je^{-\beta(a-x)}\right).
\tag{30.8}
\]

For the exponential-exponential terms, the closed double integral is

\[
\langle b_p,\mathcal K_\beta b_r\rangle
=\int_{-a}^a\int_{-a}^a e^{px+ry-\beta|x-y|}\,dy\,dx,
\tag{30.9}
\]

and for `β>|p|,|r|` it simplifies to

\[
\begin{aligned}
\frac{2\beta I(p+r)}{\beta^2-p^2}
&-\frac{e^{-(p+r)a}(1-e^{-2a(\beta-r)})}{(\beta+p)(\beta-r)}\\
&-\frac{e^{(p+r)a}(1-e^{-2a(\beta+r)})}{(\beta-p)(\beta+r)},
\quad I(z)=\int_{-a}^a e^{zx}dx.
\end{aligned}\tag{30.10}
\]

At `β=1/2` the apparent denominator zeros are handled by **exact integral limits**, not numerical regularization. Equations (30.8)–(30.10), together with the exact physical projection, generate all 30,000 original Gamma matrices. No term from the infinite digamma series is dropped on the **upper** side: every omitted term is PSD, proving the lower-matrix comparison (30.7) analytically.

## 4. New exact eight-mode positive-definiteness and physical margin

The locally **executed and independently checked** certificate uses outward fixed-point integers with denominator `10^25`, checked signed 256-bit arithmetic (overflow aborts), and exact integer floor/ceiling for every operation. The 89 initial source-attached constants come from strict Fraction/Taylor/atanh/Machin/`isqrt`/Euler–Maclaurin bounds. The compiled C++ runner creates all 32 within-parity signed interval entries at `K=30000`. A **separate Python `Fraction` verifier** reads the emitted interval table and source physical Gram data; no floating-point eigenvalue determines the result.

It establishes the stronger inequality

\[
\boxed{H_{30000}-10^{-8}M_8\succ0}\tag{30.11}
\]

by exact interval `LDLᵀ` with **all four pivots strictly positive in each parity block**. Strict verified lower pivot values, shown *rounded down for presentation*:

| Physical block, shifted by `10^-8 M` | Pivot 1 | Pivot 2 | Pivot 3 | Pivot 4 |
|---|---:|---:|---:|---:|
| `w3,w5,w7,u1` | `>0.0035129902` | `>0.00003464588` | `>0.000000845895` | `>0.0341430003` |
| `w4,w6,w8,u2` | `>0.0014009908` | `>0.00002725858` | `>0.00002981279` | `>0.0474679041` |

All strict signs are decided by exact integer comparisons on full fraction endpoints, not on the rounded displays. Thus (30.7)–(30.11) prove the new *analytical-plus-arithmetic* original signed Weil form-domain statement:

\[
\boxed{
W_{\rm original}(f,f)\ge10^{-8}\|f\|_{L^2(I)}^2
\quad\forall f\in \mathscr E_8.
}\tag{30.12}
\]

This expands CR-ALP-029's six interior directions to both genuine physical boundary anchors, and quantifies the generalized lower margin in the actual physical mass metric. It does not state that the actual full operator on the entire collar is positive.

**Independent consistency test:** all 32 enclosed finite-matrix entries were compared diagnostically against a separate floating physical-overlap/translation quadrature and finite 30,000-term Gamma implementation; the largest absolute discrepancy was below `4×10^-15`. That comparison is **diagnostic only**, not the source of exact sign certification.

## 5. Reproduction and custody

Retained executable artifacts:

1. `cr_alp_028_prime_interval_and_gram_audit.py`: prior exact original seven-prime and physical source expressions;
2. `cr_alp_030_initial_rational_data.py`: rigorous input enclosures for complete physical eight-mode Gram, original signed seven-prime matrix, Gamma constant and anchors;
3. `cr_alp_030_eight_mode_interval_core_fast.cpp`: checked signed 256-bit, outward-integer summation of **30,000** original Gamma PSD kernels and direct shifted LDL positivity test;
4. `cr_alp_030_independent_ldl_verifier.py`: separately verifies all emitted interval entries and **both** shifted physical LDL blocks using exact Python `Fraction` arithmetic;
5. `CR-ALP-030_initial_rational_intervals_1e25.txt`, `CR-ALP-030_H30000_exact_matrix_intervals.tsv`, and `CR-ALP-030_certified_eight_mode_result.json`: exact input, finite matrix enclosures, and machine-readable acceptance snapshot;
6. `run_cr_alp_030.sh`: one-command complete source regeneration, compilation and verification, using Python 3 standard library, `g++` and installed Boost Multiprecision headers.

From the artifact directory run `./run_cr_alp_030.sh`. This was executed **PASS** in the local container. It is *not* a Lean `#print axioms` certificate, and no new receiving `lean`/`lake` compile occurred.

## 6. Acceptance ledger and next gate

| Gate | CR-ALP-030 status |
|---|---|
| Original EF poles − seven paired genuine prime weights + full Gamma | Source-pinned, unchanged |
| Exact two-pole-null original eight-mode physical basis | Analytic and source-audited |
| Two true boundary anchors and physical Gram | Exact interval PASS |
| All original signed anchor-to-sine and anchor-to-anchor prime interactions | Exact interval PASS |
| 30,000-term positive-kernel Gamma lower matrix | Exact outward-integer PASS, analytic positive-kernel proof |
| Original **eight-mode** physical lower bound `10^-8` | Exact shifted `LDLᵀ` PASS; independent `Fraction` check PASS |
| Complete **15,000-mode** original finite Schur criterion | **OPEN** |
| Cross-tail coupling / residual and all remaining modes | **OPEN** |
| Original entire fixed-collar Schwartz positivity | **OPEN** |
| Original RH, arbitrary-support, Connes/Suzuki transfer | **OPEN** |
| Lean kernel compile, receiving Weil-Lab/F41/F42 acceptance | **OPEN/UNTOUCHED** |

**Next CR-ALP-031:** use the positive-semidefinite Gamma partial series together with exact physical parity blocks to certify adaptive finite-block enlargements `E_8→E_10→E_12…` and the *actual block Schur complements* (including how nearly-null physical directions couple to new modes). Separately derive computable upper bounds on the cross-tail residual and a scaling strategy for the complete `E_15000` sufficient criterion. No conclusion for whole-collar positivity or RH may be promoted without these gates.

**Governance:** retain draft/nonpromoting Canonical-Return PR #3. No writes to receiving Weil-Lab, Formalization, Hodge, F41/F42 or RH acceptance ledgers.
