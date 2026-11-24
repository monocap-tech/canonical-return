# CR-ALP-036 — The pole-null first Gamma resolvent vanishes across the high-sine cut; source-preserving signed-prime tail bound

**Research date 2026-11-23 PST; Canonical-Return, DRAFT / NONPROMOTING.** This stage verifies exact rational consequences of the original fixed-collar Weil form and derives an exact source-specific high-sine coefficient formula. It does **not** evaluate the whole original signed cross-tail residual, certify the 15,000-mode Schur matrix, recompile Lean, or prove RH. The continuum Fourier and resolvent identities are analytical proofs, not new machine-checked theorems.

## 0. Original source, sign, and accepted predecessor

The literal original `Zeta23.EF.literatureRHS` is pinned in `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, path `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, Git blob `0034976fe2f1e5fa9d411108736160542eb94293`. It has the two original poles, **minus** both-sign von Mangoldt evaluations, and the **complete** digamma bracket `Re digamma(1/4+ir/2)-log π`. The genuine finite cutoff source is `ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`. At physical radius `a=11/10`, the seven genuine prime powers are `2,3,4,5,7,8,9`, with full coefficient sum `S_a<177/50` from CR-ALP-027. No invented PSD prime term is used.

CR-ALP-034 proved analytically with an exact-rational numerical certificate `W(f,f)>(3/10^12)||f||²` on its genuine pole-null `E12` of two anchors plus `w_3,...,w_12`; CR-ALP-035 proved the **exact** high-sine residual projector

\[
Q_N=R_N-R_NB(B^*R_NB)^{-1}B^*R_N,
\qquad R_N=I-S_N,
\]

where `B` is the physical two-pole span `{cosh(x/2),sinh(x/2)}`. CR-ALP-027 analytically supplies `W(q,q)≥(1/4)||q||²` on `Q_15000 H0`. This stage concerns the true **operator cross-term** for the first twelve directions, not positivity of all intermediate physical modes.

## 1. Exact source-indexed sine coefficients of the full signed operator

Use `e_j(x)=a^{-1/2}sin(ν_j(x+a))`, `ν_j=jπ/(2a)`, and for every `k≥3` put `κ=1` if `k` odd, else `κ=2`,

\[
\lambda_k=\frac{\nu_k(\nu_\kappa^2+1/4)}{\nu_\kappa(\nu_k^2+1/4)},\qquad
w_k=e_k-\lambda_ke_\kappa.
\]

The original complete Gamma term decomposes `g(0)I+Σ_{m≥0}J_{β_m}`, with `β_m=2m+1/2` and `J_β(f,g)=∫_0^∞e^{-βs}⟨f-τ_sf,g-τ_sg⟩ds`. The exact offdiagonal ordinary sine matrix, for `j≠i`, is zero between opposite parity and, at the same parity,

\[
\langle e_j,J_\beta e_i\rangle=
-\frac{2\nu_j\nu_i}{a}
\frac{1-(-1)^j e^{-2a\beta}}
{(\beta^2+\nu_j^2)(\beta^2+\nu_i^2)}.
\tag{36.1}
\]

For **every** `j>N≥k`, the physical mass term `g(0)⟨e_j,w_k⟩` is exactly zero, and the original signed prime-plus-Gamma coefficient becomes

\[
\boxed{a_{jk}:=\langle e_j,Tw_k\rangle
=-\sum_{n\in\{2,3,4,5,7,8,9\}}
\frac{\Lambda(n)}{\sqrt n}
[U_{jk}(\log n)-\lambda_kU_{j\kappa}(\log n)]
-\frac{2\nu_j}{a}\sum_{m=0}^\infty
\frac{1-(-1)^j e^{-2a\beta_m}}{\beta_m^2+\nu_j^2}
\,\Delta_k(\beta_m),}
\tag{36.2}
\]

when `j,k` have matching parity; it is zero otherwise. Here

\[
U_{ji}(t)=\frac{2(\nu_j\sin(\nu_i t)-\nu_i\sin(\nu_jt))}{a(\nu_j^2-\nu_i^2)},
\qquad
\Delta_k(\beta)=\frac{\nu_k}{\beta^2+\nu_k^2}
-\lambda_k\frac{\nu_\kappa}{\beta^2+\nu_\kappa^2}.
\]

The term `−Σ` is the **original negative paired prime contribution**. The Gamma series is not a substitution for the original symbol: its exact equality is the previously established digamma integral representation. For fixed `j,k` its series is absolutely convergent. In the original pole-null complement, (35.3) evaluates the resulting sequence with its exact rank-two correction.

## 2. NEW exact cancellation of the first original Gamma resolvent

A rational simplification of the **genuine pole-null** coefficient yields

\[
\boxed{
\Delta_k(\beta)=
\frac{\nu_k(\nu_k^2-\nu_\kappa^2)(\beta^2-1/4)}
{(\nu_k^2+1/4)(\beta^2+\nu_k^2)(\beta^2+\nu_\kappa^2)}.}
\tag{36.3}
\]

Thus `Δ_k(1/2)=0` **exactly**, without any inequality or numerical tolerance. The first genuine digamma resolvent (`m=0`, `β_0=1/2`) has **zero** high-sine coefficients against every original pole-null interior `w_k`:

\[
\boxed{\forall j>N\ge k,\quad
\langle e_j,J_{1/2}w_k\rangle=0
\quad\Longrightarrow\quad Q_NJ_{1/2}w_k=0.}
\tag{36.4}
\]

This is a source-preserving cancellation. It does not assert that *all* remaining Gamma resolvents vanish, nor that the full residual is zero. It also does not cover the two exceptional boundary anchors automatically.

## 3. NEW strict original signed prime-only tail bound: no numerical quadrature

The original signed prime operator `P` has exactly the **minus** seven paired shift terms. For same parity and `j> N`, its high-sine coefficient on `w_k` is

\[
p_{jk}=-\sum_n\omega_n\bigl[U_{jk}(t_n)-\lambda_kU_{j\kappa}(t_n)\bigr],
\qquad (\omega_n,t_n)=(\Lambda(n)/\sqrt n,\log n).
\]

Define `Θ_k=Σ_n ω_n[sin(ν_k t_n)−λ_k sin(ν_κ t_n)]`. For the corresponding original pole `b`, its sine coefficient is

\[
c_j=\varepsilon_k\frac{2b(a)}{\sqrt a}
\frac{\nu_j}{\nu_j^2+1/4},
\qquad \varepsilon_k=+1\ (k\text{ odd, }b=\cosh(x/2)),\quad
\varepsilon_k=-1\ (k\text{ even, }b=\sinh(x/2)).
\]

Choose `t_k^P=−Θ_k/(ε_k√a b(a))`. Since the **actual** `Q_N` annihilates both original pole functions, 

\[
\|Q_NPw_k\|_2\le\|R_N(Pw_k-t_k^P b)\|_2.
\]

The signed shifts' `1/ν_j` coefficient is cancelled exactly by this scalar pole correction. Using `|sin|≤1`, `ν_j≥2ν_k`, and `S_a<177/50`, one obtains the concrete pointwise inequality

\[
|p_{jk}-t_k^Pc_j|
\le \frac{8S_a}{3a}
 \left[\frac{\nu_k+\lambda_k\nu_\kappa}{\nu_j^2}
 +\frac{\nu_k^2+\lambda_k\nu_\kappa^2}{\nu_j^3}\right]
+\frac{S_a(1+\lambda_k)}{2a\nu_j^3}.
\tag{36.5}
\]

Here `λ_k>0`, and `λ_k<κ/k+a²/(9kκ)<1` by `3<π<22/7`. Substituting rational upper bounds `ν_i<11i/(7a)` and `1/ν_j<2a/(3j)` yields entirely rational constants `A_k,B_k>0` with

\[
|p_{jk}-t_k^Pc_j|\le A_k/j^2+B_k/j^3.
\]

By Minkowski and the elementary integral test, the **entire infinite original prime residual** satisfies

\[
\|Q_{15000}Pw_k\|_2
\le \frac{A_k}{\sqrt{3(15000)^3}}
+\frac{B_k}{\sqrt{5(15000)^5}}.
\tag{36.6}
\]

The executable `Fraction` audit computes strict outward-rational upper bounds for every `k=3,...,12`, proving **each bound is `<1/30000`** and the `k=3` bound is `<1/100000`. Because the physical Gram on the ten interior columns is **exactly** `I+λλᵀ` in each parity, so `||Σx_kw_k||²≥Σ|x_k|²`, summing the ten column-error squares produces the **source-wide ten-dimensional operator estimate**

\[
\boxed{\|Q_{15000}P|_{\operatorname{span}\{w_3,\ldots,w_{12}\}}\|_{L^2\to L^2}<\frac1{18000}.}
\tag{36.7}
\]

This is a **genuine**, rigorous bound on the original signed seven-prime part of the infinite Schur residual, not a finite-dimensional toy. It does **not** bound the Gamma part or the two boundary-anchor columns. The crude prime-only bound exceeds the `≈8.66×10^-7` sufficient coupling scale associated with CR-ALP-034's tiny twelve-mode mass-relative positivity margin, so it is **insufficient to certify the combined Schur bound** even before adding Gamma and anchors.

## 4. Explicit certified Gamma-index truncation moduli, and why anchors differ

From (36.1), with `|1−(−1)^j e^{-2aβ}|≤2`, `β_m≥2m`, and `M≥2`, one gets the exact **fixed high-index pair** remainder

\[
\boxed{\left|\sum_{m=M}^{\infty}
\langle e_j,J_{\beta_m} w_k\rangle\right|
\le\frac{\nu_j(\nu_k+\lambda_k\nu_\kappa)}
{12a(M-1)^3}.}
\tag{36.8}
\]

The audit certifies the RHS `<10^-8` for the original source-indexed sample pairs `(j,k)=(15001,3),(15001,7),(15001,11),(15002,12)` at `M=30000`. This is a **per-coefficient** Gamma-series modulus. It grows with `j` and therefore **cannot be summed naively over every high sine index**; to certify a full residual one must combine a finite high-index band with a separate regularity/energy tail bound.

For the two authentic physical anchors `u_κ=Π_0e_κ`, their zero extensions have endpoint jumps. They belong to the logarithmic form and original operator `L²` domain but not automatically the zero-extended `H¹` domain. They do, however, have the elementary physical budgets

\[
D(u_\kappa)=\int_I|u'_\kappa|^2<25,
\qquad T(u_\kappa)=|u_\kappa(a)|^2+|u_\kappa(-a)|^2<4.
\tag{36.9}
\]

For example, the original pole norms are `||b_e||²=a+sinh(a)`, `||b_o||²=sinh(a)−a`; `sinh(a)−a>a³/6`, `sinh(a)<2`, `cosh(a/2)<6/5`, `sinh(a/2)<3/5`, and `3<π<22/7` suffice for the strict rational comparisons in the audit. These bounds apply specifically to `u_1,u_2` as physical pole projections, not arbitrary anchors.

Using `Jβ(e_j,e_j)≤2ν_j²/β³`, `Jβ(u,u)≤4D(u)/β³+4T(u)/β²` and Cauchy–Schwarz, for each **fixed** omitted sine index `j`,

\[
\left|\sum_{m\ge M}J_{\beta_m}(e_j,u_\kappa)\right|
\le 2\sqrt2\nu_j\bigl(\sqrt{D(u_\kappa)}\sum_{m\ge M}\beta_m^{-3}
 +\sqrt{T(u_\kappa)}\sum_{m\ge M}\beta_m^{-5/2}\bigr).
\tag{36.10}
\]

At `M=30000`, `j=15001,15002`, the source-budget rational auditor gives the conservative **per-index** bound `<1/100`. This is too coarse for the low-to-tail Schur test; it is a transparent quantification of why the boundary anchors require special treatment. It is **not** a numerical measurement of their actual couplings.

## 5. Verification, precise acceptance ledger, and next gate

Standalone audit: `src/audits/cr_alp_036_source_prime_cross_tail_audit.py` (Python standard-library `Fraction`, `isqrt`) was locally **EXECUTED PASS** for each of the ten original interior signed-prime norm estimates, the **physical** 10-column norm bound (36.7), rational resolvent cancellation regressions, Gamma fixed-index truncation budgets and the anchor derivative/trace budgets. All displayed inequalities are exact rational comparisons; the sine shift formula, pole projector, Parseval/infinite-series justification, and original Gamma-source attachment are analytical inputs, **not** new Lean kernel theorems. The algebraic Gamma cancellation is verified with finite rational parameter samples and separately derived as the identity (36.3), not formally quantified over real parameters by the executable.

**New certified mathematical boundary:** We now have an explicit bound on the **original prime-only low-to-infinite-tail operator** on every one of the ten certified interior modes. The original Gamma coupling and both anchor columns are still open, as is the complete `H_15000−4\mathcal B_15000 ⪰0` criterion and the intervening 12→15000 signed shell. A small ten-column prime-only bound does not imply the full fixed collar is nonnegative. RH, global support extension, Weil-Lab/F41 accepting states and Lean compilation remain **OPEN**. Keep Canonical-Return PR #3 **draft/nonpromoting**, with no changes to Weil-Lab, Formalization, Hodge, provenance or RH ledgers.

**Next CR-ALP-037:** Combine exact first-resolvent cancellation with a *uniform* multi-index Gamma high-sine tail enclosure. Treat the two physical boundary anchors by an endpoint-renormalized fractional regularity estimate or a certified high-sine coefficient summation with its own tail bound. Only after both are available, compare the combined original signed 12-to-15000 coupling with the tiny physical positivity margin; separately certify the intervening finite shell.