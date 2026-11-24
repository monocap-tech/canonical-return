# CR-ALP-032 — Original signed Weil: two-sided monotone Gamma enclosure and structured hierarchical Schur

**Canonical-Return, 2026-11-23 (America/Los_Angeles). RESEARCH / DRAFT / NONPROMOTING.** This stage proves an analytic upper as well as lower approximation to the complete original signed Weil form on every finite source section whose members are H¹ in the *physical interval*; it attaches this to CR-ALP-031's genuine ten-mode physical certificate and proves a new **strict upper Rayleigh witness** using entirely outward-rational data. It also identifies the source's **single rank-two prime displacement per parity** and diagonal-minus-rank-one Gamma matrix structure. It **does not** certify the full 15,000-mode signed Schur matrix or its coupling to the infinite complementary subspace, nor does it prove fixed-collar positivity, support extension, Lean kernel acceptance, or RH.

## 0. Fixed original source and inherited acceptance state

Original source: `monocap-tech/weil-lab`, `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, blob `0034976fe2f1e5fa9d411108736160542eb94293`, native `Zeta23.EF.literatureRHS`. The exact original IV011 arithmetic cutoff is blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`. The physical interval is `I=(-a,a)`, `a=11/10`. The form consists of both original pole residues, **minus** the seven genuinely weighted paired prime-power correlations at `n=2,3,4,5,7,8,9`, and the full Gamma bracket `Re ψ(1/4+ir/2)−log π`.

Restrict to the exact common kernel of the two pole moments. All formulas below refer to the **closed logarithmic Fourier form**, so the two CR-ALP-030–031 endpoint-discontinuous (on zero extension) anchors are legitimate *form-domain* sources, not incorrectly classified as literal compact smooth EF tests. The previous CR-ALP-031 analytically plus exact-arithmetically certified the entire true pole-null ten-mode physical block `E10=span(w3,w5,w7,w9,u1,w4,w6,w8,w10,u2)` and `W(f,f)≥10^-11 ||f||²` there. It did **not** certify a full source-wide complement.

## 1. Exact matrix-free structures: prime displacement rank TWO, Gamma diagonal-minus-rank-one

Let the genuine Dirichlet sine basis on `I` be

\[
e_j(x)=a^{-1/2}\sin(\nu_j(x+a)),\qquad \nu_j=j\pi/(2a).
\]

For each native prime-power shift `t_n=log n`, its paired translation matrix `U_{ij}(t_n)=⟨e_i,τ_{t_n}e_j⟩+⟨e_i,τ_{-t_n}e_j⟩` vanishes for opposite parity. The exact original-IV formula, for same parity and distinct indices, is

\[
U_{ij}(t)=\frac{2}{a}\frac{\nu_i\sin(\nu_jt)-\nu_j\sin(\nu_it)}{\nu_i^2-\nu_j^2}.
\]

Define the **original signed** prime matrix `P=−Σ_{n∈{2,3,4,5,7,8,9}} Λ(n)n^-1/2 U(t_n)`; all repeated Mangoldt prime-power weights are preserved. With `X=diag(ν_j²)` and `b_j=Σ_n Λ(n)n^-1/2 sin(ν_j log n)`, the **entire seven-prime block** on each parity has the exact displacement identity

\[
\boxed{XP-PX=-\frac2a(\nu b^T-b\nu^T),\qquad \operatorname{rank}(XP-PX)\le2.}\tag{32.1}
\]

This is **rank-two displacement**, not rank-two of the prime matrix itself. The diagonal entries remain the literal original paired translation values. The 40 off-diagonal identities for source indices 1–10 have been independently checked against CR-ALP-028's outward-rational prime interval bank.

For `β_m=2m+1/2` and the *complete original Gamma* positive-square kernels `J_m(f,g)`, the sine restriction similarly has, for equal parity,

\[
\boxed{J_{m,ij}=\delta_{ij}\frac{2\nu_i^2}{\beta_m(\beta_m^2+\nu_i^2)}
-\frac{2(1-(-1)^i e^{-2a\beta_m})}{a}\frac{\nu_i}{\beta_m^2+\nu_i^2}\frac{\nu_j}{\beta_m^2+\nu_j^2}.}\tag{32.2}
\]

Opposite-parity entries vanish. Thus **each Gamma kernel is a diagonal-minus-one-rank-one matrix per parity**. The entire projected physical finite section uses `w_j=e_j−λ_j e_{k(j)}` with exact pole-null rank-one corrections, plus the two genuine projected exponential anchor columns. These formulas enable matrix-free operations while preserving the physical Gram; they do **not** furnish a certified spectrum for 15,000 modes on their own.

## 2. NEW universal TWO-SIDED tail bound for the actual Gamma series

For `f∈H¹(I)` (possibly with nonzero boundary traces, zero-extended outside I), put

\[
D(f)=\int_I|f'(x)|^2dx,\qquad T(f)=|f(-a)|^2+|f(a)|^2.
\]

The zero extension has distributional derivative

\[
\partial_x(1_I f)=1_I f'+f(-a)\,\delta_{-a}-f(a)\,\delta_a.
\]

For all `s>0`, the fundamental theorem of calculus for the absolutely continuous interior and its two jumps, followed by Minkowski, yields

\[
\|f-\tau_s f\|_2\le s\|f'\|_{L²(I)}+\sqrt{s}\bigl(|f(-a)|+|f(a)|\bigr).
\]

It follows, for every `β>0`, that

\[
\boxed{0\le J_\beta(f,f)=\int_0^\infty e^{-\beta s}\|f-\tau_s f\|_2^2ds
\le\frac4{\beta^3}D(f)+\frac4{\beta^2}T(f).}\tag{32.3}
\]

If both endpoint traces vanish, one has the sharper `J_β(f,f)≤2D(f)/β³`. Since `β_m=2m+1/2` and the positive functions `(2x+1/2)^-p` decrease, for any integer `K≥1`

\[
0\preceq\sum_{m=K}^\infty J_m
\preceq c_D(K)D+c_T(K)T,
\]

where the explicit rigorous rational coefficients are

\[
\boxed{c_D(K)=\frac4{\beta_K^3}+\frac1{\beta_K^2},\qquad
c_T(K)=\frac4{\beta_K^2}+\frac2{\beta_K}.}\tag{32.4}
\]

For the boundary-zero subspace, use `c_0(K)=2/β_K³+1/(2β_K²)` instead. At **K=30,000**, the exact Fraction audit verifies

\[
c_D<\frac1{3\,000\,000\,000},\qquad c_T<\frac1{29\,000},\qquad
c_0<\frac1{7\,000\,000\,000}.
\]

These estimates are two-sided *on any finite H¹ section*: the full signed Weil finite matrix is above its genuine signed finite Gamma partial sum, and below that partial sum plus the displayed **positive-semidefinite derivative and endpoint Grams**. They do **not** directly upper-bound the cross action against arbitrary logarithmic-form vectors in `Q_15000`, since these vectors need not be H¹.

## 3. NEW two-sided error for the existing 617 Loewner groups

CR-ALP-031 already computed a rigorously enclosed **lower** matrix

\[
\widetilde H= g(0)M+P+\sum_{m=0}^{29999}J_m+
\sum_{[L,U)\in\mathcal G}(U-L)J_{U-1},
\]

where `\mathcal G` consists of 617 consecutive integer groups covering `[30000,350000)` and each group was validated by a separate exact group-schedule audit. The Loewner-decreasing family `J_m` yields **both** endpoint inequalities

\[
(U-L)J_{U-1}\preceq\sum_{m=L}^{U-1}J_m\preceq(U-L)J_L.
\]

To control the gap without calculating another matrix, for `β≤β'` note

\[
e^{-\beta s}-e^{-\beta's}\le(\beta'-\beta)s e^{-\beta s}.
\]

Inserting (32.3)'s derivation into the difference integral gives the precise positive-form bound

\[
0\preceq J_\beta-J_{\beta'}\preceq(\beta'-\beta)
\left(\frac{12}{\beta^4}D+\frac8{\beta^3}T\right).
\]

Thus a group with `ℓ=U−L` contributes an *upper error relative to the already certified lower endpoint* bounded by

\[
\ell\bigl(J_L-J_{U-1}\bigr)
\preceq\frac{24\ell(\ell-1)}{\beta_L^4}D+
\frac{16\ell(\ell-1)}{\beta_L^3}T.
\]

Sum these **exact rational** coefficients over all 617 groups, then add the tail (32.4) from `m≥350000`. The independent Fraction audit proves the source-preserving **two-sided matrix/form inequality**

\[
\boxed{\widetilde H\ \preceq\ W_{\mathrm{original}}|_{E_N}\ \preceq\
\widetilde H+\frac1{180\,000\,000\,000}D_N+\frac1{320\,000}T_N}\tag{32.5}
\]

for any finite `E_N` of the same pole-null H¹ physical basis, **provided** the lower matrix `\widetilde H` is built by exactly that 30,000+617 grouping construction for `E_N`. At present only the source-exact `\widetilde H` for `E_10` has been computed and certified; the general formula is a construction, not a numerical certificate for larger N.

The actual rational coefficient sums are smaller: `cD,total≈5.340831690×10^-12`, `cT,total≈3.100214659×10^-6`. These decimals are explanatory only; the inequalities above are exact. The derivative Gram `D_N` and endpoint trace Gram `T_N` are **not invented identity matrices**.

## 4. Exact sparsity of the physical derivative and boundary matrices

Within a parity block, let `k=1` for the even physical anchor or `k=2` for the odd anchor, `B_1=a+sinh a`, `B_2=sinh a−a`, `c_k=<e_k,b_k>`, `γ_k=c_k/B_k`, `b_1=cosh(x/2)`, `b_2=sinh(x/2)`. Since `w_j` has zero endpoint traces **and is exactly pole-null**, integration by parts gives

\[
\langle w_j',b_k'\rangle=-\tfrac14\langle w_j,b_k\rangle=0.
\]

The physical derivative Gram follows exactly:

\[
\begin{aligned}
D(w_i,w_j)&=\nu_i^2\delta_{ij}+\lambda_i\lambda_j\nu_k^2,\\
D(w_j,u_k)&=-\lambda_j\nu_k^2,\\
D(u_k,u_k)&=\nu_k^2+\frac{c_k^2}{2B_k}+\frac{c_k^2 B_{3-k}}{4B_k^2}.
\end{aligned}\tag{32.6}
\]

The boundary trace matrix is supported **only on the two anchor coordinates**:

\[
T(u_1,u_1)=2\gamma_1^2\cosh^2(a/2),\qquad
T(u_2,u_2)=2\gamma_2^2\sinh^2(a/2),
\]

and all `T(w_i,·)=0` and `T(u_1,u_2)=0`. The source interval bank certifies `D` strictly positive on each five-dimensional parity block; the endpoint Gram is exactly rank two on the combined ten-mode section. Its two nonzero diagonal magnitudes are about `0.84926945` and `1.68889315`, respectively (diagnostic display). This is the **rank-two boundary correction** that will remain small-memory as N grows, though its magnitude can still be large relative to an extremely small spectral margin.

## 5. NEW certified near-null original full-form Rayleigh witness

Use the **existing CR-ALP-031 independently certified grouped lower matrix and its original physical Gram**. Choose the following **explicit rational coefficient vectors**, each with common denominator `10^12`:

- Even block in order `(w3,w5,w7,w9,u1)`:
  `(-727226995401, 643702820298, -182669753899, 9336596462, 37487631) / 10^12`.
- Odd block in order `(w4,w6,w8,w10,u2)`:
  `(-743143741936, 620073220342, -145278032240, -1238582207, -82669993) / 10^12`.

The auditor **re-evaluates both Rayleigh numerators and physical denominators using only exact Fraction arithmetic on the prior outward intervals**, and uses (32.5) for a mathematically safe *upper* bound on the omitted **complete** original Gamma energy. Every strict comparison is made before decimal display. The result is

\[
\boxed{\frac{W(f_{\rm even},f_{\rm even})}{\|f_{\rm even}\|_2^2}
<\frac1{2\,500\,000\,000}=4\cdot10^{-10},\qquad
\frac{W(f_{\rm odd},f_{\rm odd})}{\|f_{\rm odd}\|_2^2}
<\frac1{250\,000\,000}=4\cdot10^{-9}.}\tag{32.7}
\]

The exact-interval test's computed *upper* displays with the rounded-up simple coefficients are about `3.433×10^-10` and `2.894×10^-9`. These are certificates of unusually **small positive directions**, not negative witnesses. Combined with the already-certified CR-ALP-031 full-form lower bound, (32.7) yields

\[
\boxed{10^{-11}\le\lambda_{\min}(W_{\rm original}|_{E_{10}},M_{10})
<4\times10^{-10}.}\tag{32.8}
\]

This is a new **two-sided** spectral enclosure for the complete original signed ten-mode physical form, conditional on the original continuum identities. It does not claim that the original collar's global lowest eigenvalue obeys the same bounds; enlarging the section can lower the Rayleigh infimum.

## 6. Genuine hierarchical Schur architecture and its open gates

For `E_10 ⊂ E_M ⊂ E_15000` (all in the genuine physical pole-null form domain), use *physical* orthogonal shells. To prove the entire fixed collar, a finite block certificate must be combined with the already-derived `Q_15000` positivity `W(q,q)≥(1/4)||q||²` and a certified upper bound for its coupling to the entire `E_15000`.

Equation (32.5) can rigorously enclose **finite shell** signed matrices and finite-shell cross terms: for arbitrary finite `f,g`, positivity of the Gamma tail and Cauchy–Schwarz in the positive tail form imply

\[
|W(f,g)-\widetilde H(f,g)|^2\le U(f,f)U(g,g),\quad
U=\frac{D}{180\,000\,000\,000}+\frac{T}{320\,000}.
\]

This permits an *interval-certified hierarchical finite Schur recursion* without blindly trusting rounded eigenvalues. **But the bound is not an operator-norm estimate against arbitrary `Q_15000` logarithmic-form directions**. It cannot currently discharge the 15,000-mode-to-infinite-tail residual; any claim that it does would be a false promotion.

Furthermore the explicit Rayleigh upper bound `4×10^-10` shows why crude off-diagonal errors can swamp actual small positive spectral margins. The raw 30,000-term boundary upper coefficient is about `3.33×10^-5` and cannot by itself close a `10^-11` spectral gate. The monotone grouped tail and source-specific trace structure are necessary, not cosmetic.

**Next CR-ALP-033:** use the exact rank-two prime displacement / diagonal-plus-rank-one Gamma shell structure to develop a certified interval **block Schur transform**, first checking a 10→12-mode shell with full signed cross entries. Then obtain genuinely effective `E_15000↔Q_15000` residual bounds, perhaps through source-specific Fourier regularity/weighted form estimates; the finite H¹ tail enclosure here is not sufficient for that final step.

## 7. Verification and custody

Executed `src/audits/cr_alp_032_structured_schur_tail_audit.py`: **PASS**, using Python stdlib `Fraction`, the pinned prior source interval bank, and *no* numerical quadrature for acceptance. Specifically: 40 original sine-prime displacement checks, symbolic-rational Gamma rank-one structure, exact 617-group/320,000-index coverage, rational sum-to-integral tail majorants, source-attached H¹ derivative and rank-two trace Grams with ten positive interval pivots, both exact rational physical Rayleigh witness upper comparisons, and source SHA256 consistency.

Retained artifacts: this note; `cr_alp_032_structured_schur_tail_audit.py`; `CR-ALP-032_structured_spectral_result.json`; execution log; and a reproduction package copying the unchanged predecessor `cr_alp_028_prime_interval_and_gram_audit.py`, `CR-ALP-031_initial_rational_intervals_1e25.txt`, and `CR-ALP-031_grouped_H_exact_matrix_intervals.tsv`. The audit's `float(...)` outputs are **display only**; all strict decisions precede them in exact rational arithmetic.

**Not established:** source-wide original Weil PSD; sufficiency of F41's selected 242 physical sources; all 15,000-mode signed positivity; whole-space Schur residual; arbitrary support extension; Suzuki/Connes transfer; RH; or Lean kernel certification. **Governance:** retain Canonical-Return PR #3 as DRAFT/NONPROMOTING; no receiving Weil-Lab, F41/F42, Formalization, Hodge or RH ledger modification.