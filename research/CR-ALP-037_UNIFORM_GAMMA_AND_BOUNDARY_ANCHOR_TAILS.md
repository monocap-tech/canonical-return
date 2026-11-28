# CR-ALP-037 — Uniform original Gamma infinite-tail bounds and endpoint-renormalized boundary anchors

**Date:** 2026-11-23, America/Los_Angeles. **Status:** analytical continuum reduction plus locally executed exact `Fraction` inequality audit; **DRAFT / NONPROMOTING**. No Lean compilation, entire 15,000-mode Schur certificate, original two-pole-null fixed-collar positivity, arbitrary-radius transfer, or RH proof.

## 0. Original source, fixed space and predecessor custody

We retain the original `Zeta23.EF.literatureRHS` at `monocap-tech/weil-lab` branch `research/iv068-original-source-kernel-gate`, path `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, source blob `0034976fe2f1e5fa9d411108736160542eb94293`. Its two genuine poles, negative paired von Mangoldt interactions at the original seven powers `{2,3,4,5,7,8,9}`, and complete digamma symbol `Re ψ(1/4+ir/2)−log π` are retained. Write `P=−Σ_{n} (Λ(n)/√n)(τ_{log n}+τ_{−log n})` and `T=T_Γ+P` restricted to `I=(−a,a), a=11/10`. Pole-null sources annihilate both original pole evaluations.

The original pole functions are `b_e(x)=cosh(x/2)` and `b_o(x)=sinh(x/2)`; `Π_0` projects to their joint orthogonal complement. The physical Dirichlet sines are `e_j(x)=a^{-1/2}sin(ν_j(x+a))`, `ν_j=jπ/(2a)`. The original source modes from CR-ALP-028 are `w_k=e_k−λ_k e_{κ(k)}` for `3≤k≤12`, with `κ(k)=1` for odd k and `2` for even k and

\[
λ_k=\frac{ν_k(ν_κ^2+1/4)}{ν_κ(ν_k^2+1/4)},\qquad
u_κ=Π_0e_κ=e_κ−γ_κ b_κ,\quad γ_κ=\frac{\langle e_κ,b_κ\rangle}{\|b_κ\|_2^2}.
\]

The actual physical twelve-dimensional space is `E_12=span{w_3,...,w_12,u_1,u_2}`. Its ordinary physical metric must be retained: the `w` Gram is `I+λλᵀ` in each parity but the full `E_12` Gram includes nontrivial anchor cross-entries.

CR-ALP-035 gives the exact infinite-dimensional pole-null tail projector at `N=15000`:

\[
Q_N=R_N-R_NB(B^*R_NB)^{-1}B^*R_N,\quad R_N=I-S_N,
\quad Q_Nb_e=Q_Nb_o=0.
\]

CR-ALP-036 already established, by a separately executed exact rational audit, `||Q_15000 P|_{span w_3,...,w_12}||_{L²→L²}<1/18000`. CR-ALP-027 gives an *analytical*, not Lean-certified, positive tail form lower bound `W(q,q)≥(1/4)||q||²` on `Q_15000 H_0`. CR-ALP-034 provides a separate exact-arithmetic finite twelve-mode `W|E_12>(3/10^12)M_12` certificate. None of these previous statements evaluates the full 15,000-dimensional signed Schur matrix.

## 1. NEW complete infinite-frequency Gamma estimate for all ten interior directions

The complete genuine Gamma term has the analytically established expansion

\[
T_Γ=g(0)I+\sum_{m=0}^\infty J_{β_m},\qquad β_m=2m+\tfrac12,
\qquad J_β(f,f)=\int_0^\infty e^{-βs}\|f-τ_sf\|_2^2\,ds.
\]

The equality is in the original Fourier multiplier/operator domain, with zero extension of the compact physical sources. For `j>N≥k≥3` in the same parity, CR-ALP-036 gives

\[
\langle e_j,T_Γw_k\rangle
=-\frac{2ν_j}{a}\sum_{m=1}^\infty
\frac{(1-σ_j e^{-2aβ_m})\Delta_k(β_m)}{ν_j^2+β_m^2},\quad σ_j=(-1)^j,
\]

where the exact genuine pole-null cancellation is

\[
\Delta_k(β)=\frac{ν_k(ν_k^2-ν_κ^2)(β^2-1/4)}{(ν_k^2+1/4)(β^2+ν_k^2)(β^2+ν_κ^2)}.
\]

In particular `Δ_k(1/2)=0`, so the first genuine Gamma resolvent has **zero** high-sine coupling. For each fixed source/parity set

\[
S_k=\sum_{m=1}^\infty (1-σ_j e^{-2aβ_m})\Delta_k(β_m),
\quad \text{which converges absolutely},
\]

and subtract the corresponding scalar multiple `t_k^Γ b_κ` whose coefficient is `−(2ν_j/a) S_k/(ν_j^2+1/4)`. This is permitted because the *exact original* `Q_N` annihilates `b_κ`. The difference is

\[
\langle e_j,T_Γw_k-t_k^Γb_κ\rangle
=\frac{2ν_j}{a(ν_j^2+1/4)}\sum_{m=1}^\infty
\frac{(1-σ_j e^{-2aβ_m})\Delta_k(β_m)(β_m^2-1/4)}{ν_j^2+β_m^2}.
\tag{37.1}
\]

This term is **nonnegative**. Put `C_k=ν_k(ν_k²−ν_κ²)/(ν_k²+1/4)<ν_k`; then `0≤Δ_k(β)(β²−1/4)≤C_k` for `β≥1/2`. At `m≥1`, `β_m≥5/2`, so `e^{-2aβ_m}≤e^{-11/2}<1/200` (finite exact Taylor inequality), hence `1−σ_j e^{-2aβ_m}<201/200`. The decreasing-series integral comparison gives

\[
\sum_{m=1}^\infty\frac1{ν_j^2+β_m^2}
\le\sum_{m=1}^\infty\frac1{ν_j^2+4m^2}
\le\int_0^\infty\frac{dx}{ν_j^2+4x^2}
=\frac{π}{4ν_j}.
\]

Therefore, **simultaneously for every omitted j**,

\[
\boxed{0\le \langle e_j,T_Γw_k-t_k^Γb_κ\rangle
<\frac{201}{200}\frac{k}{j^2}.}\tag{37.2}
\]

Opposite parity coefficients vanish exactly. Parseval and `Σ_{j>N}j^{-4}≤1/(3N³)` yield

\[
\|Q_NT_Γw_k\|_2^2
<\frac{(201/200)^2k^2}{3N^3}.
\]

The ten-column physical Gram satisfies `||Σ_{k=3}^{12}x_kw_k||²≥Σ|x_k|²` (rank-one updates `I+λλᵀ` by parity). Since `Σ_{k=3}^{12}k²=645`, the **entire infinite-dimensional** original Gamma tail, without frequency or resolvent truncation, satisfies at `N=15000`

\[
\boxed{\|Q_{15000}T_Γ|_{\mathrm{span}\{w_3,\ldots,w_{12}\}}\|_{L²\to L²}
<\frac1{120000}.}\tag{37.3}
\]

Combining (37.3) with the CR-ALP-036 source-signed seven-prime bound **by the triangle inequality**, retaining all negative original weights,

\[
\boxed{\|Q_{15000}T|_{\mathrm{span}\{w_3,\ldots,w_{12}\}}\|_{L²\to L²}
<\frac1{120000}+\frac1{18000}=\frac{23}{360000}.}\tag{37.4}
\]

The rational script independently checks `(201/200)²·645/(3·15000³)<1/120000²`, with no floating-point comparison.

## 2. NEW exact endpoint-renormalized anchor Gamma coefficient identity

Both genuine anchor sources `u_κ=e_κ−γ_κ b_κ` have zero-extension endpoint jumps. The interior `H¹` argument of CR-ALP-035 is **not** applied to them. Instead use the exact physical Laplace kernel `K_βf(x)=∫_{-a}^a e^{-β|x-y|}f(y)dy` and `J_β=2I/β−K_β`. For any pole function `b∈{b_e,b_o}` and sine index j of matching parity, let

\[
c_j=\langle e_j,b\rangle=\frac{2b(-a)ν_j}{\sqrt a(ν_j^2+1/4)},\quad
I_b(β)=\int_0^{2a}e^{-βt}b(-a+t)dt.
\]

An exact integration of `K_β e_j` against `b`, including both physical endpoint exponentials, gives

\[
\boxed{\langle e_j,J_βb\rangle
=\frac{2ν_j^2}{β(β^2+ν_j^2)}c_j
-\frac{2ν_j}{\sqrt a(β^2+ν_j^2)}I_b(β).}\tag{37.5}
\]

Set `E_β=I_b(β)−b(-a)/β`. Then **algebraically**

\[
\boxed{\langle e_j,J_βb\rangle-\frac{c_j}{β}
=-\frac{2ν_j}{\sqrt a(ν_j^2+β^2)}
\left[E_β+\frac{b(-a)(β^2+1/4)}{β(ν_j^2+1/4)}\right].}\tag{37.6}
\]

The term `c_j/β` is the dangerous harmonic contribution associated with the physical endpoint jump. Crucially, it is **exactly aligned with the original pole function**, so it can be removed before applying `Q_N`.

Choose `M=10000`, and subtract the **finite** scalar pole correction `H_Mb`, `H_M=Σ_{m=0}^{M}1/β_m`, from `Σ_m J_{β_m}b`. For every `j>N=15000`, combine (37.6) on `m≤M` with the original (37.5) for `m>M`. The coefficients of the new remainder are globally square-summable, and no infinite divergent pole coefficient is introduced.

For `B0≥||b||∞`, `B1≥||b'||∞`, the source-specific error satisfies

\[
|E_β|\le B_1/β^2+B_0 e^{-2aβ}/β.
\]

Using `e^a>3`, `e^{4a}>80`, and the elementary reciprocal-square sum,

\[
\sum_{m≥0}|E_{β_m}|
<\frac92B_1+\frac{160}{237}B_0=:S_E.
\]

Furthermore, for `ν_j≥ν_{15000}>2M`, the high-β remainder obeys

\[
\sum_{m>M}\frac{ν_j}{β_m(β_m^2+ν_j^2)}
\le \frac1{ν_j}\left[\frac12\log\frac{ν_j}{2M}+\frac9{16}\right].
\]

The first term is a harmonic comparison through `m≈ν_j/2`; the second comes from the reciprocal-cubic tail and `floor(ν_j/2)≥ν_j/3`. The fixed-source ratio satisfies

\[
\log\frac{ν_j}{2M}\le \frac1{14}+\log\frac{j}{15000},
\]

using `π<22/7` and `log(1+t)≤t`. Put

\[
S_β=(M+1)(M+\tfrac12)+\tfrac14(2+M/2),\quad
ν_* =\frac{3(15000)}{2a},
\]

\[
C_0=2S_E+\frac{2B_0 S_β}{ν_*^2}+\frac{27B_0}{8}+\frac{3B_0}{14},\qquad C_1=3B_0.
\]

Then **uniformly for every `j>15000`** after the finite pole subtraction,

\[
\boxed{\left|\sum_m\langle e_j,J_{β_m}b\rangle-H_Mc_j\right|
\le \frac{11}{15j}\left(C_0+C_1\log\frac{j}{15000}\right).}\tag{37.7}
\]

Since the right-hand side divided by j is decreasing on `j≥15000`, a single exact integral proves the full square-tail bound

\[
\boxed{\left\|R_N\left(\sum_{m≥0}J_{β_m}b-H_Mb\right)\right\|_2^2
\le \left(\frac{11}{15}\right)^2\frac{C_0^2+2C_0C_1+2C_1^2}{N}.}\tag{37.8}
\]

The sine component `e_κ` of the actual anchor can likewise be corrected by a **convergent** pole-aligned resolvent sum; the resulting full Gamma high-sine residual has column norm below `(201/200)κ/√(3N³)<10^{-6}`. The scalar `g(0)u_κ` contributes *zero* after the actual projection because `u_κ∈E_N`.

At the original `a=11/10`, exact Taylor/π/square-root rational enclosures of the genuine anchor coefficients prove `|γ_1|<3/5`, `|γ_2|<5/3`. For the cosh pole use `B0=6/5, B1=3/10`; for sinh use `B0=3/5, B1=3/5`. The standalone exact `Fraction` audit proves, from (37.8), the genuine **full infinite-Gamma** anchor-column estimates

\[
\boxed{\|Q_{15000}T_Γu_1\|_2<\frac1{20},\qquad
\|Q_{15000}T_Γu_2\|_2<\frac19.}\tag{37.9}
\]

These are **individual original physical source-column bounds**, not the operator norm on the entire physically conditioned span `E_12`.

## 3. NEW complete signed-prime anchor estimate via true zero-extension BV

The original seven-prime operator remains `P=−Σ_n ω_n(τ_{t_n}+τ_{-t_n})`, with `S_a=Σ_n ω_n<177/50`. For a zero-extended physical source `u` of bounded variation, `Var_R(Pu)≤2S_a Var_R(u)`; restricting to `I` and extending by zero cannot increase the total variation of a compactly supported BV function. Integration by parts in the original physical sine basis gives

\[
|\langle e_j,Pu\rangle|\le\frac{2S_a\operatorname{Var}_{\mathbb R}(u)}{\sqrt a\,ν_j}
<\frac{2a}{3j}\,2S_a\operatorname{Var}_{\mathbb R}(u).
\]

With the exact two authentic anchors `u_κ=e_κ−γ_κb_κ`, `Var_R(e_κ)=2κ/√a`, and the endpoint jumps of `b_κ`, the previously proved rational envelopes yield

\[
\operatorname{Var}_{\mathbb R}(u_1)<4,\qquad
\operatorname{Var}_{\mathbb R}(u_2)<8.
\]

Consequently the infinite j-sums are bounded directly, with the authentic signed seven prime powers unchanged:

\[
\boxed{\|Q_{15000}Pu_1\|_2<\frac15,\qquad
\|Q_{15000}Pu_2\|_2<\frac25.}\tag{37.10}
\]

Combining these with (37.9) by the triangle inequality gives **complete original signed Gamma-plus-prime column estimates**

\[
\boxed{\|Q_{15000}Tu_1\|_2<\frac14,\qquad
\|Q_{15000}Tu_2\|_2<\frac{23}{45}.}\tag{37.11}
\]

The negative prime signs have **not** been replaced by positive interactions; their absolute values enter only a conservative upper norm bound.

## 4. Exact scope: all 12 columns now bounded, but the physical Schur gate is OPEN

| Source sector | Full Gamma infinite-tail residual | Original seven-prime residual | Combined signed residual |
|---|---:|---:|---:|
| Ten interior directions `span(w_3,...,w_12)` | `<1/120000` (operator norm on interior physical span) | `<1/18000` (prior genuine operator norm) | `<23/360000` (operator norm on interior physical span) |
| `u_1` boundary anchor | `<1/20` (individual column) | `<1/5` (individual column) | `<1/4` (individual column) |
| `u_2` boundary anchor | `<1/9` (individual column) | `<2/5` (individual column) | `<23/45` (individual column) |

**The three rows are not interchangeable:** the two anchor bounds are *not* normalized by the true physical mass Gram and must not be combined as a twelve-dimensional `L²→L²` operator certificate without accounting for the anchor/interior Gram off-diagonals. The conservative bounds are far larger than the ~`8.66×10^-7` sufficient coupling threshold inferred from the extremely small `3×10^-12` CR-ALP-034 twelve-mode margin and the `1/4` CR-ALP-027 entire-tail coercivity margin.

Thus, while all twelve *source columns* now have explicit, **uniform over every high sine index** residual budgets, the original genuine signed `H_15000−4\mathcal B_15000⪰0` remains **UNEVALUATED**. Nor is the signed physical intermediate shell from `E_12` to `E_15000` certified. There is no original negative Weil witness and no claim that the crude residual bounds are close to optimal. The endpoint renormalization is exact and source-preserving; it is not a substitute for a verified Schur lower bound.

## 5. Exact arithmetic audit and next investigative gate

Standalone audit `src/audits/cr_alp_037_uniform_gamma_tail_audit.py`, Python standard-library `Fraction`, locally **EXECUTED PASS**: finite Taylor exponential inequalities, ten-interior uniform complete Gamma squared norm, carry-forward arithmetic for the genuine signed seven-prime bound, exact rational enclosures of both *genuine* anchor coefficients, `M=10000` endpoint-renormalized Gamma whole-tail squared norm inequalities, true zero-extended BV prime-anchor high-sine norms and combined signed column ceilings.

The analytic identities (37.1), (37.5)–(37.8), positivity/Parseval, displacement of zero extensions and the original digamma source attachment are mathematical arguments rather than new Lean-compiled declarations. The script checks the rational ceilings after applying them and does not magically formalize any continuum lemma. No receiving Weil-Lab, F41/F42, Formalization, Hodge/provenance, or RH acceptance record is modified; keep canonical-return PR #3 **DRAFT / NONPROMOTING**.

**Next CR-ALP-038:** Use the exact physical `E_12` Gram and source-specific signed *combined* entries to tighten the two boundary-anchor residuals; then evaluate an actual blockwise signed Schur comparison against an intermediate shell. The issue is now **strength and conditioning** of the quantitative cross-tail bounds—not the former absence of a globally summable Gamma anchor estimate.
