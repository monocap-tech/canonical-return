# CR-ALP-034 — Loewner-convex Gamma midpoint: complete original signed twelve-mode certificate

**2026-11-23, America/Los_Angeles. Canonical-Return research only / DRAFT / NONPROMOTING.** The exact finite arithmetic was executed by a checked-int256 C++ lower-matrix builder and a separate Python `Fraction` verifier. The continuum Fourier/Gamma identities, Loewner convexity, and form domain are human analytical inputs, **not newly compiled Lean theorems**. Neither whole-collar positivity, the 15,000-mode Schur gate, F41 historical sufficiency, arbitrary-support extension, nor RH is asserted.

## 0. Genuine native source and retained custody

- Original signed Weil explicit formula: `monocap-tech/weil-lab`, `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, blob `0034976fe2f1e5fa9d411108736160542eb94293`, `Zeta23.EF.literatureRHS`.
- Native finite-prime support/correlation: `WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`.
- Physical collar `I=(-a,a)`, `a=11/10`; both original pole evaluations vanish **exactly** on the two-pole-null sources. The genuine prime contribution remains **minus** paired translations at `{2,3,4,5,7,8,9}` weighted by `Λ(n)/sqrt(n)`, with `Λ(4)=Λ(8)=log 2`, `Λ(9)=log 3`. The full Gamma symbol is `Re ψ(1/4+ir/2)-log π`.
- Source-indexed rational starting intervals, original signed prime cross entries, full Gamma constant, physical mass and two boundary anchors are the **same authenticated CR-ALP-033 inputs**, SHA256 `bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2`. This stage changes only the Gamma grouped lower-bound *rule*, not the source, physical subspace or metric.

Define the original Dirichlet sines `e_j(x)=a^-1/2 sin(jπ(x+a)/(2a))`, `ν_j=jπ/(2a)`. With `k(j)=1` for odd `j` and `k(j)=2` for even `j`, set

\[
\lambda_j=\frac{\nu_j(\nu_{k(j)}^2+1/4)}{\nu_{k(j)}(\nu_j^2+1/4)},\quad
w_j=e_j-\lambda_j e_{k(j)},\quad u_k=\Pi_0e_k.
\]

All `w_j` and both `u_k` are genuinely two-pole-null. The exact physical metric `M`, **not** an identity surrogate, is retained. The complete original twelve-mode target is

\[
E_{12}=\operatorname{span}\{w_3,w_5,w_7,w_9,w_{11},u_1\}
\oplus\operatorname{span}\{w_4,w_6,w_8,w_{10},w_{12},u_2\}.
\]

CR-ALP-033 had certified only the eleven-mode subspace without `w11`: a grouped *endpoint* Gamma lower bound produced a **negative fifth pivot** after adding it. That was not a negative full Weil witness, since the omitted Gamma terms were PSD.

## 1. NEW source-preserving Loewner convexity inequality

The original complete Gamma term has, analytically,

\[
W_\Gamma(f,f)=g(0)\|f\|_2^2+\sum_{m=0}^\infty J_{\beta_m}(f,f),\quad
\beta_m=2m+\tfrac12,\quad
J_\beta(f,f)=\int_0^\infty e^{-\beta s}\|f-\tau_sf\|_2^2\,ds.
\]

For all form-domain `f` and `β>0`, the Laplace kernels are PSD and *Loewner-convex* in `β`:

\[
\frac{d^2}{d\beta^2}J_\beta(f,f)
=\int_0^\infty s^2e^{-\beta s}\|f-\tau_sf\|_2^2ds\ge0.
\]

Consequently, for **each** group `[L,U)` of consecutive integer `m`, `ℓ=U−L`, the midpoint `\bar\beta=\frac1\ell\sum_{m=L}^{U-1}\beta_m=L+U-\tfrac12` yields **Jensen's Loewner lower bound**

\[
\boxed{\sum_{m=L}^{U-1}J_{2m+1/2}\ \succeq\ \ell\,J_{L+U-1/2}.}\tag{34.1}
\]

This is stronger than the earlier decreasing-Loewner endpoint bound `ℓ J_{2(U−1)+1/2}` because `\bar\beta\le\beta_{U-1}`. Importantly, it needs **one** exactly enclosed midpoint kernel per group, **not** a separate costly summation of the thousands of entries in that group.

The CR-ALP-033 schedule is reused literally: the first **30,000** original Gamma terms individually, and **4,219** groups covering `[30,000,2,000,000)`. The midpoint parameter is an **exact positive rational** with numerator `(2L+2U−1)` and denominator `2`, so every sine/exponential Laplace formula remains within checked directed fixed-point integer arithmetic. For grouped β≥60000, the factor `exp(-2aβ)` is rigorously enclosed by `[0,10^-25]` as in the predecessor. Every unprocessed `m≥2,000,000` Gamma kernel is nonnegative and safely omitted from the *lower* form.

Define the signed lower matrix on `E12`

\[
\boxed{\widetilde H_{34}=g(0)M_{12}+P_{\rm original}
+\sum_{m=0}^{29999}J_{\beta_m}|_{E12}
+\sum_{[L,U)\in\mathcal G}(U-L)J_{L+U-1/2}|_{E12}.}\tag{34.2}
\]

The original form satisfies `W|E12 ⪰ \widetilde H_34` by positive-semidefinite Jensen, **not** by assuming positivity of the uncomputed signed operator.

## 2. NEW complete twelve-mode positivity with actual physical metric

The checked 256-bit outward-integer builder constructs **all 72 ordered within-parity matrix entries**, including both anchors and the newly added original `w11` and `w12`. The independent `Fraction` verifier reads the exact matrix-bank file and the untouched starting source intervals, enforces symmetric intersection of outward intervals, and checks the physical mass Grams by rational `LDLᵀ`.

At the stronger conservative shift `σ=3/10^12`, every shifted exact interval `LDLᵀ` pivot of `\widetilde H_34−σM_12` is **strictly positive** in both parity blocks (6/6 even and 6/6 odd). The important previously unresolved fifth even pivot has an exact lower endpoint of `3583183386415032639/10^25 > 0`; this is a pivot of the *shifted* lower matrix, not a generalized eigenvalue itself. The last even pivot is also strictly positive; the complete **six-dimensional** even test has been checked, not merely `w11` in isolation.

Thus, from the analytical Gamma-series identity and the exact arithmetic,

\[
\boxed{\forall f\in E_{12}\setminus\{0\},\quad
W_{\rm original}(f,f)>\frac3{10^{12}}\|f\|_{L^2(I)}^2.}\tag{34.3}
\]

The prior negative fifth pivot belonged only to the **weaker endpoint lower approximation** and is superseded by a positive pivot for the Jensen-midpoint lower approximation. No actual signed negative Weil witness was found, and no RH hypothesis is needed for (34.3).

## 3. NEW sharper **two-sided** Gamma remainder on finite H¹ sections

The midpoint lower approximation also admits a better upper error bound, without assuming any finite-block positivity. As in CR-ALP-032, for physical `f∈H¹(I)`, let

\[
D(f)=\int_I|f'(x)|^2dx,\qquad T(f)=|f(-a)|^2+|f(a)|^2.
\]

The zero-extension translation estimate implies `\|f−\tau_sf\|_2²≤2s²D(f)+4sT(f)`. Differentiating twice under the convergent Laplace integral gives

\[
0\preceq J_\beta''\preceq\frac{48}{\beta^5}D+\frac{24}{\beta^4}T.
\]

Taylor's inequality about each midpoint, combined with `\sum_{m=L}^{U-1}(\beta_m−\bar\beta)^2=\ell(\ell^2−1)/3`, yields

\[
0\preceq\sum_{m=L}^{U-1}J_{\beta_m}-\ell J_{\bar\beta}
\preceq\frac{8\ell(\ell^2-1)}{\beta_L^5}D
+\frac{4\ell(\ell^2-1)}{\beta_L^4}T.
\tag{34.4}
\]

The infinite tail from index `K=2,000,000` is bounded, exactly as in CR-ALP-032, by `c_D(K)=4/β_K³+1/β_K²` and `c_T(K)=4/β_K²+2/β_K`. The independent verifier adds **each positive rational coefficient rounded upward separately at `10^-30` precision**; including this finite tail, the resulting rational bounds are

\[
 c_D\le\frac{62636009953159599}{10^{30}}<\frac1{15\times10^{12}},
\quad
 c_T\le\frac{500008263973209505225640}{10^{30}}<\frac1{1,999,000}.
\]

Therefore, **on every finite physical H¹ section** built using this same source (not just on `E12`), the original signed form is squeezed as

\[
\boxed{
0\preceq W_{\rm original}-\widetilde H_{34}
\preceq\frac{D}{15\times10^{12}}+\frac{T}{1,999,000}.
}\tag{34.5}
\]

The `T` matrix has rank at most two on the finite projected physical basis; it vanishes on all interior `w_j`. This remainder is much sharper in the derivative form than CR-ALP-033's `D/10^12` bound. **It still does not control the mixed Schur coupling against arbitrary log-form vectors in the infinite-dimensional complement of the 15,000-mode section.**

## 4. Exact audit custody, reproduction, and acceptance ledger

`src/audits/cr_alp_034_midpoint_12_mode_interval_core.cpp` is a self-contained checked-int256 implementation of (34.2) using the pinned source interval input. It executes `30,000` individual kernels, the same `4,219` groups with exact midpoint beta, and all genuine signed original prime entries. Its acceptance check is the shifted physical test `\widetilde H_34−(3/10^12)M_12≻0` in **both** parity blocks.

`src/audits/cr_alp_034_independent_midpoint_ldl_verifier.py` independently uses only `Fraction` and integer arithmetic, verifies every physical and shifted pivot, and certifies the (34.5) upper-error arithmetic. It writes the machine-readable result JSON. Local execution status: **C++ checked-integer PASS; independent Python Fraction PASS**. The uploaded input SHA256 is `bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2`; the generated midpoint matrix SHA256 is `b6a31e623298272fdceac788fc56626ec704e280903ff6baff1dc8bfe0d45f92`.

**Status:** CR-ALP-034 source-wide Gamma Jensen/Loewner convexity analytically proved; full genuine original signed twelve-dimensional physical section **EXACT ARITHMETIC PASS**; full 15,000-mode Schur PSD and its infinite tail interaction **OPEN**; receiving Lean kernel proof **OPEN**; arbitrary collar support-extension, F41 242 sufficiency and RH **OPEN**. Retain PR #3 in DRAFT/NONPROMOTING; no receiving Weil-Lab, Formalization, Hodge, F41/F42, provenance or RH acceptance ledger mutations.

**Next CR-ALP-035:** decide between scalable hierarchical finite-section Schur computation and a direct rigorous mixed-residual bound on the 15,000-mode complement. The Jensen midpoint technique removes a major Gamma grouping cost, but a growing list of positive finite sections alone does **not** establish positivity of the complete collar.
