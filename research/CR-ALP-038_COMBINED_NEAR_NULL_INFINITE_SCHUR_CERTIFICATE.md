# CR-ALP-038 — Original signed Weil: two explicit near-null lines couple positively to the full infinite tail

**Research date:** 2026-11-23 PST. **Status:** new source-preserving classical analytic Schur theorem, backed by locally executed exact-`Fraction` arithmetic and an independent numerical diagnostic. **DRAFT / NONPROMOTING.** This is not a native Lean theorem, not a positivity result on all `E_12`, not a certificate of `E_15000`, not an arbitrary-support result and not RH.

## 0. Native source and certified inputs

The genuine original fixed-collar signed Weil form is the pole-null restriction of `Zeta23.EF.literatureRHS` at the pinned Weil-Lab branch `research/iv068-original-source-kernel-gate`, path `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, blob `0034976fe2f1e5fa9d411108736160542eb94293`. It contains the **original two poles**, both-sign **negative** von Mangoldt prime interactions for exactly `n∈{2,3,4,5,7,8,9}` at `a=11/10`, and the complete Gamma multiplier `Re ψ(1/4+ir/2)−log π`. No sign change or PSD surrogate is introduced.

Let `I=(−a,a)`, `a=11/10`, and `b_e=cosh(x/2)`, `b_o=sinh(x/2)`. The pole-null subspace is `H_0=ker B*`, with physical Dirichlet sines `e_j(x)=a^−1/2 sin(ν_j(x+a))`, `ν_j=jπ/(2a)`. Define exact pole-null interior sources

\[
w_k=e_k-\lambda_ke_{\kappa(k)},\qquad
\kappa(k)=\begin{cases}1,&k\text{ odd},\\2,&k\text{ even},\end{cases}\qquad
\lambda_k=\frac{\nu_k(\nu_{\kappa(k)}^2+1/4)}{\nu_{\kappa(k)}(\nu_k^2+1/4)}.
\]

The `w_k` have zero physical endpoint traces. They belong to the original logarithmic form domain and operator domain; their zero extensions are `H¹` and have Fourier decay `O(r^−2)`, although they are not themselves necessarily `C_c^∞`. `T` is the complete source-signed Gamma-minus-seven-prime operator and `Q_N` is CR-ALP-035's exact **rank-two-corrected** pole-null high-sine projector.

The following predecessor statements are ANALYTICAL with locally executed exact arithmetic; neither is Lean certified:

1. **CR-ALP-027**: `W(q,q)≥(1/4)||q||²` for all original pole-null `q∈Q_15000 H_0` in the logarithmic form domain.
2. **CR-ALP-034**: `W(p,p)>(3/10^12)||p||²` for all nonzero `p∈E_12=span{w_3,...,w_12,u_1,u_2}`. Its original interval input and midpoint matrix hashes were rechecked locally in CR-ALP-038; the unchanged independent `Fraction` verifier passed.

The original seven-weight upper bound `S_a=Σ_{n∈{2,3,4,5,7,8,9}} Λ(n)/√n <177/50` and `e^(11/2)>200` are retained from the exact predecessor audits. The subsequent bounds are **on the full, not truncated, sine tail**.

## 1. New near-null rational sources

Fix precisely the following **rational** coefficient vectors; the irrational physical `λ_k` remain exact as defined above:

\[
\begin{aligned}
p_o&=\frac1{10^7}(6872105w_3-6735272w_5+2338391w_7-188293w_9-20112w_{11}),\\
p_e&=\frac1{10^7}(-7021868w_4+6604165w_6-1934708w_8-26547w_{10}+43236w_{12}).
\end{aligned}
\]

These were selected *numerically* from the near-null generalized eigenvectors of CR-ALP-034's original physical `H_mid/M` and rounded to seven decimal places; the exact rational certificate below **does not depend on numerical eigenvalues**. They lie in opposite physical reflection parities and are mutually orthogonal. They are also orthogonal to `Q_15000 H_0` by construction.

Write either source `p` as a finite ordinary sine expansion `p=Σ_{i∈J}d_i e_i`, including the anchor index `κ=1` or `2`, where

\[d_k=x_k\ (k\ge3),\qquad d_\kappa=-\sum_k x_k\lambda_k.
\]

The exact original pole constraint is `Σ_i d_i ν_i/(ν_i²+1/4)=0`. Define the **signed** leading coefficient and the three unsigned high-frequency moments:

\[D(p)=\sum_i d_i\nu_i,\qquad C_s(p)=\sum_i |d_i|\nu_i^s,\quad s=0,2,3.\]

**The essential new observation:** after subtracting the exact original pole-aligned `j^-1` coefficient, the leading **Gamma and prime** `j^-2` terms are both governed by the same small signed scalar `D(p)`. This cancellation is destroyed if each `w_k` is bounded separately before combining them.

## 2. Original signed prime high-sine coefficients with common leading scalar

For `j>N≥12` in the same parity, the original paired shift has the exact source-indexed formula

\[U_{ji}(t)=\frac{2}{a}\frac{\nu_j\sin(\nu_it)-\nu_i\sin(\nu_jt)}{\nu_j^2-\nu_i^2}.\]

The authentic prime coefficient is `p_j=−Σ_n ω_n Σ_i d_i U_ji(log n)`, `ω_n=Λ(n)/√n>0`. Put `Θ=Σ_n ω_n Σ_i d_i sin(ν_i log n)`, and let `c_j=〈e_j,b〉` for the matching original pole `b∈{cosh(x/2),sinh(x/2)}`. There is a unique scalar `t^P_p` making `t^P_pc_j=−(2Θ/a)ν_j/(ν_j²+1/4)`; thus `Q_N(t^P_pb)=0` **exactly**.

Expand `1/(ν_j²−ν_i²)=1/ν_j²+ν_i²/[ν_j²(ν_j²−ν_i²)]`. Since `ν_i/ν_j≤12/15001<1/2` and hence `1/[1−(ν_i/ν_j)²]<2`, the original signed-prime remainder obeys **for every** high `j`:

\[
\boxed{\left|p_j-t^P_pc_j\right|\le\frac{2S_a|D(p)|}{a\nu_j^2}+\frac{S_aC_0(p)/2+8S_aC_2(p)}{a\nu_j^3}.}\tag{38.1}
\]

This is a bound for all seven genuine negative prime-power interactions; `S_a` was used only for the upper error envelope, never to replace the original signed form.

## 3. Entire Gamma-resolvent tail with the same signed scalar

For the genuine `β_m=2m+1/2`, define `Δ_p(β)=Σ_i d_iν_i/(β²+ν_i²)`. Pole-nullity gives `Δ_p(1/2)=0`, so the first Gamma resolvent couples to no high sine mode. For `m≥1`, the exact source-indexed expansion from CR-ALP-036/037 allows a pole-aligned scalar `t_p^Γ` to be subtracted:

\[
\langle e_j,T_Γp-t_p^Γb\rangle=\frac{2\nu_j}{a(\nu_j^2+1/4)}\sum_{m=1}^\infty
\frac{(1-(-1)^je^{-2a\beta_m})(\beta_m^2-1/4)\Delta_p(\beta_m)}{\nu_j^2+\beta_m^2}.
\]

The key identity, preserving all signed cancellations, is

\[\Delta_p(\beta)=\frac{D(p)}{\beta^2}-\sum_i\frac{d_i\nu_i^3}{\beta^2(\beta^2+\nu_i^2)}.\]

Consequently, for `β≥5/2`,

\[\left|(\beta^2-1/4)\Delta_p(\beta)\right|\le |D(p)|+\frac{C_3(p)}{\beta^2}.\]

Because `e^{-2aβ_m}<1/200` for `m≥1`, `Σ_{m≥1}1/(ν_j²+β_m²)≤π/(4ν_j)` and `Σ_{m≥1}β_m^{-2}<1/2`, the **complete infinite Gamma series**, without finite-resolvent truncation, satisfies

\[
\boxed{\left|\langle e_j,T_Γp-t_p^Γb\rangle\right|
\le \frac{(201/200)\pi |D(p)|}{2a\nu_j^2}
+\frac{(201/200)C_3(p)}{a\nu_j^3}.}\tag{38.2}
\]

The constant Gamma multiplier `g(0)I` contributes no high-sine coefficient for interior sources with index ≤12; both `t_p^Γb` and `t_p^Pb` vanish under the **actual** rank-two-corrected `Q_N`.

## 4. Full infinite-dimensional original signed Schur residual bound

With `3<π<22/7`, `ν_j>3j/(2a)`, define the following **exact positive quantities** (the rational audit uses an outward bound on `π`):

\[
\begin{aligned}
A(p)&=\left(\frac{2S_a}{a}+\frac{201\pi}{400a}\right)|D(p)|\left(\frac{2a}{3}\right)^2,\\
B(p)&=\left(\frac{S_aC_0(p)}{2a}+\frac{8S_aC_2(p)}a+\frac{201C_3(p)}{200a}\right)\left(\frac{2a}{3}\right)^3.
\end{aligned}
\]

Combining (38.1) and (38.2), for the original signed `T=T_Γ+P`,

\[|\langle e_j,Tp-(t^P_p+t^Γ_p)b\rangle|\le A(p)/j^2+B(p)/j^3,\quad j>15000.\]

By the *actual* original pole-null projector, `Q_Nb=0`, Parseval and the integer-integral tests `Σ_{j>N}j^-4≤1/(3N³)`, `Σ_{j>N}j^-6≤1/(5N^5)`, we obtain the genuinely uniform bound

\[
\boxed{\|Q_NTp\|_2\le\frac{A(p)}{\sqrt{3N^3}}+\frac{B(p)}{\sqrt{5N^5}}.}\tag{38.3}
\]

The squared source norm has the exact physical expression `||p||²=Σ x_k²+(Σ_k x_kλ_k)²≥Σ x_k²`. The executable source-specific `Fraction` audit uses a 32-term alternating Machin arctangent enclosure

\[\pi=16\arctan(1/5)-4\arctan(1/239)\]

and explicitly rounds rational bounds outwards. **No floating-point arithmetic enters any acceptance comparison.**

| Original source | Outward-rational full `Q_15000 T` norm upper | Certified `W(p,p)` lower from CR-ALP-034 | Upper for `4||QTp||²` |
|---|---:|---:|---:|
| `p_o` | `<1.793851×10^-8` | `>2.942809×10^-12` | `<1.287161×10^-15` |
| `p_e` | `<2.837677×10^-8` | `>2.900018×10^-12` | `<3.220964×10^-15` |

Every displayed decimal is rounded in the conservative direction; the **actual comparisons are exact rational inequalities**. The raw numerical `D` values are `≈+0.00172671112954` (odd) and `≈−0.00512119830804` (even), versus unsigned `Σ_i|d_i|ν_i` of order `10`. Both survive rational approximation to the near-null eigenvectors.

## 5. NEW genuinely mixed finite-plus-infinite positivity theorem

For `q∈Q_15000 H_0∩𝒟_log`, let `q_o,q_e` be its two physical reflection parities. The original signed operator, the true pole-null space and `Q_N` respect this parity split. If `f=αp_o+βp_e+q`, the exact Weil form is bounded, using the certified CR-ALP-027/034 predecessors and the new **entire** cross-tail bounds, by

\[
W(f,f)\ge m_o|\alpha|^2+m_e|\beta|^2+	frac14\|q\|_2^2
-2c_o|\alpha|\|q_o\|_2-2c_e|\beta|\|q_e\|_2,
\]

where `m_o>(3/10^12)Σ_odd x_k²`, `m_e>(3/10^12)Σ_even x_k²` and `c_{o/e}` are the rational upper bounds from (38.3). The exact audit checks the **stronger** inequalities `m_o−8c_o²>28/10^13` and `m_e−8c_e²>28/10^13`.

Young's inequality `2c|α|r≤(1/8)r²+8c²|α|²` therefore proves the explicit **strict quantitative mixed-form theorem**

\[
\boxed{W(\alpha p_o+\beta p_e+q,\alpha p_o+\beta p_e+q)
>\frac{28}{10^{13}}(|\alpha|^2+|\beta|^2)+\frac18\|q\|_2^2,}\tag{38.4}
\]

for every nonzero permitted combination (the strict sign at zero is interpreted in the usual way). Because `p_o,p_e` have opposite reflection parity and both are orthogonal to the whole high tail, this is a genuine **infinite-dimensional positive mixed subspace**, not merely separate tests on the finite vectors.

**Precisely NOT claimed:** this does *not* establish full `E_12⊕Q_15000` positivity, positivity of the other 10 directions in `E_12` against the high tail, the enormous signed physical shell `E_15000⊖E_12`, the complete fixed-collar original signed form, support extension, Suzuki–Connes transfer, Lean compilation, or RH. The two rational p-vectors are form-domain test sources, not an asserted complete basis or an F41 selected-mode certificate.

## 6. STRONGER RESULT — the ENTIRE ten-dimensional interior span couples positively to the whole infinite tail

The two near-null-line inequalities are not the endpoint. The same source-preserving coefficient envelope can be elevated to a **matrix inequality for every simultaneous combination of all ten interior directions**.

For either five-element parity index set `K`, set `p=Σ_{k∈K}x_kw_k`, define the exact signed vector `D_k=ν_k−λ_kν_κ`, and define the conservative positive source-column bounds

\[
\begin{aligned}
c_{0,k}&=1+\lambda_k,\qquad
c_{2,k}=\nu_k^2+\lambda_k\nu_\kappa^2,\qquad
c_{3,k}=\nu_k^3+\lambda_k\nu_\kappa^3,\\
\alpha&=\left(\frac{2S_a}{a}+\frac{201\pi}{400a}\right)\left(\frac{2a}{3}\right)^2\frac1{\sqrt{3N^3}},\\
\beta_k&=\left(\frac{S_a c_{0,k}}{2a}+\frac{8S_a c_{2,k}}a+\frac{201c_{3,k}}{200a}\right)
\left(\frac{2a}{3}\right)^3\frac1{\sqrt{5N^5}},\quad N=15000.
\end{aligned}
\]

The signed cancellation is retained for **every** coefficient vector `x` in the first term, while the conservative second term uses a coordinatewise triangle inequality:

\[
\boxed{\|Q_N T(Vx)\|_2\le\alpha|D^\mathsf T x|+\sum_{k\in K}\beta_k|x_k|
\le\alpha|D^\mathsf T x|+\Bigl(\sum_k\beta_k^2\Bigr)^{1/2}\|x\|_2.}\tag{38.5}
\]

Therefore, by `(r+s)^2≤2r²+2s²`, the complete source-specific residual has a **Loewner upper bound**

\[
\boxed{\mathcal B_N^{\mathrm{interior}}
\preceq 2\alpha^2DD^\mathsf T+2\Bigl(\sum_k\beta_k^2\Bigr)I.}\tag{38.6}
\]

This matrix is **rank one plus diagonal** in each reflection-parity block. Its rank-one direction is precisely the signed `D` responsible for the near-null cancellation; no separate worst-case column norms are multiplied into a fictitious full-dimensional bound.

Now combine CR-ALP-034's original **full signed** 12-mode Jensen midpoint lower matrix `H_mid` with CR-ALP-027's positive entire tail `μ=1/4`. Young's inequality with `1/8` reserved as a positive tail margin requires only the **explicit ten-interior-mode finite test**

\[
\boxed{H_{\mathrm{mid},K}-16\alpha^2DD^\mathsf T
-16\left(\sum_k\beta_k^2\right)I
-10^{-12}M_K\succ0.}\tag{38.7}
\]

**This test now PASSES in both parity blocks by rigorous outward-rounded integer interval `LDLᵀ` at scale `10^-25`.** The independent original CR-ALP-034 midpoint matrix and physical mass inputs were byte-verified against their fixed SHA-256 values before use. No floating-point comparison enters any positivity pivot.

| True original signed parity block | Interval `LDLᵀ` positive pivots | Certified smallest pivot lower |
|---|---:|---:|
| Five odd-index interior `w_3,w_5,w_7,w_9,w_11` | 5/5 | `>8×10^-7` |
| Five even-index interior `w_4,w_6,w_8,w_10,w_12` | 5/5 | `>5×10^-6` |

The positive pivots refer to the **shifted lower matrix (38.7)** in the displayed coefficient basis. They are *not* asserted to be generalized eigenvalues, nor are the numerically obtained near-null spectral values promoted as exact enclosures.

Because the original form, physical metric, two pole constraints and `Q_N` respect the parity decomposition, the two blockwise `LDLᵀ` certificates combine. We obtain the strictly stronger **mixed finite-plus-infinite coercivity theorem**

\[
\boxed{W(p+q,p+q)>10^{-12}\|p\|_{L^2(I)}^2+\frac18\|q\|_{L^2(I)}^2}
\tag{38.8}
\]

for every **nonzero** `p+q` with `p∈span_C{w_3,...,w_12}` and `q∈Q_15000 H_0∩𝒟_log`. This includes *all ten* interior physical directions and *every* high-frequency pole-null tail vector. It remains **conditional on the unformalized classical continuum identities and predecessor tail coercivity**, while the finite arithmetic comparison is fully exact and reproducible.

**What is still missing:** both boundary-anchor directions `u_1,u_2`, the entire intervening signed physical sector `E_15000⊖E_12`, and their cross-couplings. Thus (38.8) does NOT prove `E_12⊕Q_15000`, let alone the whole fixed collar or RH. The two-line theorem of §5 is retained as a more transparent independent subcertificate.

## 7. Numerical pilot: diagnostic, separately marked

A separate SciPy diagnostic evaluated the original complete signed Gamma resolvent through the **exact digamma partial-fraction identity** plus exponentially small endpoint correction; it used the authentic paired prime shift coefficients from CR-ALP-028 and removed the pole-aligned `j^-1` tail exactly. For the rational near-null interior vectors, the numerical combined norm in the first 15,000 omitted same-parity sine indices was `~4.24×10^-10` odd and `~1.28×10^-9` even. These are band-limited **noncertified** measurements, not upper bounds for the full infinite tail. A separate 12-dimensional numerical band-Schur diagnostic found the generalised minimum eigenvalue of `H_mid−4 B_band` still positive. Because `B_band⪯B_full`, that positivity is **only a necessary, NOT sufficient** check for the true Schur inequality. The exact-rational two-line theorem above does not depend on the pilot.

## 8. Audit, provenance, and next gate

The new source-attached `src/audits/cr_alp_038_full_interior_infinite_schur_audit.py` executes an independent *outward-integer interval LDL* on both complete five-mode interior blocks, with all operations rounded outward at `10^-25` and a verified `10^-12` physical-metric shift. It passed after the exact scalar `Fraction` cancellation audit was completed. The executable standard-library audit `src/audits/cr_alp_038_combined_infinite_tail_rational_audit.py` was **locally executed PASS**. It verifies the Machin rational enclosure, each exact physical λ and signed D interval, rational bounds for `C0,C2,C3`, pointwise envelope constants, **the entire infinite j-sum**, and both strict Young-Schur margins. The existing CR-ALP-034 independent `Fraction` verifier was **rerun PASS** with original `Hmid` matrix SHA-256 `b6a31e623298272fdceac788fc56626ec704e280903ff6baff1dc8bfe0d45f92` and source-input SHA `bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2`; CR-ALP-036/037 exact audit scripts were rerun PASS as independent source custody checks. The Gamma/prime real-analysis coefficient formulas, infinite Plancherel and operator-domain handling remain *human analytical proofs*, **not Lean kernel-checked**. No receiving Weil-Lab changes.

**Next — CR-ALP-039:** address both genuine boundary anchors `u_1,u_2` with a source-preserving pole-renormalized interval residual matrix, and separately treat the still-uncertified signed physical shell `E_15000⊖E_12`. The ten interior directions plus the complete infinite tail are now certified as a joint positive sector; this is not yet the entire fixed collar.
