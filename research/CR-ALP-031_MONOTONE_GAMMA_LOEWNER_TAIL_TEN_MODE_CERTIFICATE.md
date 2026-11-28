# CR-ALP-031 — Monotone Loewner Gamma-tail grouping certifies the genuine original signed ten-mode Weil block

**Canonical-Return, November 23, 2026 (America/Los_Angeles). RESEARCH / DRAFT / NONPROMOTING.** Stage CR-ALP-031 expands CR-ALP-030's true physical eight-generator two-pole-null block by the exact two new original projected sine combinations `w_9` and `w_10`. It preserves the source-pinned original two-pole-minus-seven-paired-prime-plus-full-Gamma Weil form. It supplies a human-analytic/independently exact-arithmetic ten-mode result, **not** an all-collar PSD theorem, an F41 242-mode sufficiency theorem, an original RH proof, a Lean compilation, or acceptance-ledger promotion.

## 0. Source provenance and inherited identities

Native source `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, Git blob `0034976fe2f1e5fa9d411108736160542eb94293`, `Zeta23.EF.literatureRHS`. The original negative prime pairing and seven exact von Mangoldt weights are pinned by `ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean`, blob `14d8f01ce34f764900eccc68399d7323b4f4a25c`. The full digamma symbol remains `Re ψ(1/4+ir/2) − log π`. The exact classical Gamma-series identity and physical source/Gram equations are derived in CR-ALP-028–030; the analytic bridge has *not* been freshly Lean-compiled.

Let `I=(-a,a)`, `a=11/10`, and set `ν_j=jπ/(2a)`, `e_j(x)=a^{-1/2}sin(ν_j(x+a))` in the physical interval. For each `j≥3` let `k(j)=1` for odd j and 2 for even j, and

\[
λ_j=\frac{ν_j(ν_{k(j)}^2+1/4)}{ν_{k(j)}(ν_j^2+1/4)},\quad
w_j=e_j-λ_je_{k(j)},\quad u_k=\Pi_0 e_k.
\]

The complete genuine ten-mode physical pole-null space is

\[
\mathscr E_{10}=\operatorname{span}_{\mathbb C}\{w_3,w_5,w_7,w_9,u_1\}
\oplus\operatorname{span}_{\mathbb C}\{w_4,w_6,w_8,w_{10},u_2\}.
\tag{31.1}
\]

The first block is physically even, the second odd. In either parity block the actual **physical** Gram remains the exact low-rank formula

\[
M=\begin{pmatrix}I_4+λλ^T&-λ\\-λ^T&1-c_k^2/B_k\end{pmatrix},\quad
B_1=a+\sinh a,\quad B_2=\sinh a-a,
\tag{31.2}
\]
with `c_k=⟨e_k,b_k⟩`, `b_1=cosh(x/2)`, `b_2=sinh(x/2)`. All ten sources are exactly pole-null. In contrast to the interior `w_j`, anchors `u_1,u_2` have nonzero endpoint jumps after zero extension and are initially elements of the **closed logarithmic form domain**, not falsely asserted compact smooth inputs to the literal EF.

## 1. The crucial failed straightforward enlargement (retained)

Using the exact source-enclosed original signed prime block and **only** the first `K=30,000` positive-semidefinite Gamma kernels inherited from CR-ALP-030, the ten-mode lower matrix `H_30000` is not positive. Independently exact interval LDL for the first four *interior* modes in each parity block gives:

| Finite lower-matrix Schur pivot | Rigorous sign | Approximate enclosure (display only) |
|---|---|---|
| even, adding `w_9` after `w_3,w_5,w_7` | negative | `(-6.304386e-5, -6.304379e-5)` |
| odd, adding `w_10` after `w_4,w_6,w_8` | negative | `(-2.597876e-3, -2.597873e-3)` |

These are rigorous signs **of the truncated Loewner lower matrix**, NOT negative vectors for the complete original Weil form. The omitted Gamma kernels are PSD and can change the sign. A numerical generalized-eigenvalue diagnostic for `H_30000` gives minimum approximately `−4.54×10^-9` even and `−4.14×10^-9` odd, close enough to zero that omission matters. No RH contradiction may be inferred from the failed partial-matrix test.

An exploratory individual-term `100,000`-kernel exact run was too expensive and did not finish in the allotted execution interval. A separate floating-point continuation is retained solely as a diagnostic and predicts recovery of positivity with more positive Gamma kernels. The certificate below **does not depend on** those floating-point estimates.

## 2. NEW Loewner-monotone grouped Gamma tail lemma

The complete original Gamma energy on any supported log-form source `f` has the classical, exact positive-square decomposition

\[
W_\Gamma(f,f)=g(0)\|f\|_2^2+\sum_{m=0}^\infty J_m(f,f),\qquad
J_m(f,f)=\int_0^\infty e^{-\beta_m s}\|f-\tau_s f\|_2^2\,ds,\quad
\beta_m=2m+\tfrac12.
\tag{31.3}
\]

For all `m≤n`, because the squared difference is nonnegative pointwise and `β_m≤β_n`, one has the **operator-form Loewner inequality** `J_m ⪰ J_n ⪰ 0`, without a Fourier approximation or RH assumption. Thus, for any integers `0≤L<U`,

\[
\boxed{\sum_{m=L}^{U-1}J_m\ \succeq\ (U-L)J_{U-1}.}
\tag{31.4}
\]

This is the missing device: replace a long run of individually evaluated positive-semidefinite Gamma kernels by a **single endpoint kernel with a certified lower weight**. It is mathematically a lower sum for a Loewner-decreasing operator-valued sequence, not merely a quadrature guess.

## 3. Source-preserving 30,000 + 617-block effective lower matrix

Construct a deterministic integer partition of `[30000,350000)` by `start=30000`, `length=max(1,floor(start/250))`, `end=min(350000,start+length)`, advancing `start=end`. The exact independent schedule audit proves there are **617 blocks**, covering every integer index from `30000` through `349999` exactly once.

Define the **finite** lower matrix

\[
\boxed{\widetilde H= g(0)M + P_{\rm seven\ original\ signed\ primes}
+\sum_{m=0}^{29999}J_m
+\sum_{[L,U)\text{ in 617 groups}}(U-L)J_{U-1}.}
\tag{31.5}
\]

The full original two-pole-null signed Weil form satisfies the rigorous form-order inequality

\[
\boxed{H_{\mathrm{original}}|_{\mathscr E_{10}}\succeq\widetilde H.}
\tag{31.6}
\]

The infinite remainder `m≥350000` is dropped **only from the lower side**, as it is PSD. Thus this test never assumes omitted Gamma energy vanishes. The exact physical formulas for every sine/sine, sine/exponential and exponential/exponential `J_m`, including both anchors, are unchanged from CR-ALP-030, as are original negative paired prime weights.

For grouped endpoints `m≥30000`, their exact source factor `q_m=e^{-2aβ_m}=e^{-a(4m+1)}` lies inside `[0,10^-25]`: the original enclosure has `e^-a<2/5<1/2`, and `(1/2)^100<10^-25` by the integer inequality `2^100>10^25`. Every grouped endpoint therefore uses the outward fixed-point range `[0,10^-25]`; the neglected exponential is NOT simply replaced by a falsely exact zero. The first 30,000 terms use the original exact recurrence `q_m=e^{-a}(e^{-4a})^m` with outward rounding. The complete matrix is accumulated in **checked signed 256-bit integers** at scale `10^25`; overflow aborts the certification.

## 4. Exact positive-definite ten-mode result

The executed C++ 30,000-term plus 617-block interval summation produced a complete 50-entry within-parity matrix interval bank (25 entries per block). Source constants, the original prime interactions and the physical Gram are outward enclosed by an independent Python rational generator.

In both physical parity blocks, *all five* exact interval LDL pivots of `\widetilde H - (1/10^11) M` have strictly positive lower endpoints. Independently, Python `Fraction` recomputed the LDL decomposition from the emitted matrix and physical Gram intervals and verified all strict signs. Rounded-down displays of the rational lower pivot values:

| `\widetilde H - 10^-11 M` | Pivot 1 | Pivot 2 | Pivot 3 | Pivot 4 | Pivot 5 |
|---|---:|---:|---:|---:|---:|
| even `w3,w5,w7,w9,u1` | `>0.00351300420` | `>0.00003472567` | `>0.000002234070` | `>0.000002240396` | `>0.02625193` |
| odd `w4,w6,w8,w10,u2` | `>0.00140100822` | `>0.00002733954` | `>0.00003197728` | `>0.001939440` | `>0.04149313` |

These are *pivots of an unnormalized physical-basis congruence*, not directly eigenvalues. The interval verifier checks the real strict inequalities using the full exact fractions. The physical Gram is separately positive definite, and both blocks include all genuine original signed cross-entries.

By (31.6), the result is the **human-analytic-plus-executed exact-arithmetic** lower bound

\[
\boxed{\forall f\in\mathscr E_{10},\quad W_{\rm original}(f,f)\ge 10^{-11}\|f\|_{L^2(I)}^2.}
\tag{31.7}
\]

It genuinely enlarges the previous eight-mode source space (not merely replaces the basis). Its certified lower constant is smaller than CR-ALP-030's `10^-8`, reflecting the newly exposed near-null directions. As a **diagnostic only**, floating-point eigenvalues of the finite Loewner lower matrix `\widetilde H` are approximately `1.53×10^-10` (even) and `2.62×10^-9` (odd), comfortably above the rigorously selected `10^-11` threshold. The exact interval signs, not these floating values, establish (31.7).

## 5. Independent custody and reproduction

The retained reproducible files are:

- `cr_alp_028_prime_interval_and_gram_audit.py`: authentic seven von Mangoldt weights and elementary original signed prime overlap intervals; also controls `π`, trig, log and physical projection moments.
- `cr_alp_031_initial_rational_data.py`: exact rational 129-constant input generator for both 5×5 physical Gram and original signed prime/Gamma offset blocks.
- `cr_alp_031_ten_mode_interval_core.cpp`: 30,000 exact Gamma kernels, 617 monotone endpoint-weight groups, checked signed 256-bit outward interval summation; writes 50 genuine signed lower-matrix enclosures.
- `cr_alp_031_monotone_group_schedule_audit.py`: separate exact integer/rational test for 617 group endpoints, coverage, strict upper bound on exponential factors, and full symmetric matrix interval bank.
- `cr_alp_031_independent_ldl_verifier.py`: independent Python `Fraction` physical-gram/shifted LDL proof and results JSON.
- `CR-ALP-031_initial_rational_intervals_1e25.txt`, `CR-ALP-031_grouped_H_exact_matrix_intervals.tsv`, `CR-ALP-031_certified_ten_mode_result.json` and concise executed logs; `run_cr_alp_031.sh` reproduces the certificate with Python 3, `g++` and Boost Multiprecision headers.
- `cr_alp_031_numerical_gamma_tail.py`: **diagnostic only** floating-point continuation; never used for rational certification.
- `CR-ALP-031_H30000_exact_matrix_intervals.tsv`: genuine truncated-only source-original matrix, used to document how a negative lower-matrix pivot can be repaired by a PSD Gamma tail. It is NOT a negative signed Weil witness.

**Executed arithmetic status:** 129 outward source constants PASS; all original seven-prime interval checks PASS; `30,000+617` checked-int Loewner lower matrix PASS; independent group/coverage/formula-enclosure audit PASS; ten exact rational physical-metric shifted LDL pivots PASS. The classical digamma sum identity, paper Fourier attachment, operator-valued Loewner monotonicity and zero-extension form domain statements remain written **human proofs**, not Lean kernel theorems. In particular no machine-checked proof of original RH or source-wide signed positivity is asserted.

## 6. Acceptance gates and next investigation

| Obligation | CR-ALP-031 status |
|---|---|
| Original EF complete poles − seven paired primes + full Gamma | Source-pinned, signs unchanged |
| True physical 10-mode two-pole-null source with both anchors | Exact source arithmetic + analytic definition |
| Original 30,000-kernel 10-mode lower matrix | **INDEFINITE**, exact interval observed; not a Weil counterexample |
| NEW monotone Gamma Loewner interval grouping | Analytically proved + exact coverage PASS |
| Genuine original signed 10-mode physical margin `10^-11` | **EXACT INTERVAL + INDEPENDENT RATIONAL PASS** |
| Entire original two-pole-null fixed collar and 15,000-mode Schur | **OPEN** |
| Actual 10-mode-to-infinite-tail coupling or finite Schur residual | **OPEN** |
| F41 original 242-mode physical sufficiency; global RH | **OPEN** |
| Lean kernel and receiving acceptance ledgers | **NOT PERFORMED** |

**Next CR-ALP-032:** use the successful monotone Loewner grouping to extend source-aware parity blocks further, but prioritize certified low-to-high cross-coupling norms and a realistic *hierarchical* tail/Schur strategy. Requiring a single 15,000×15,000 dense LDL is not the only possible acceptance architecture; any hierarchical alternative must still control all omitted directions in the original signed form. PR #3 remains DRAFT/NONPROMOTING, with no receiving Weil-Lab, F41/F42, Formalization, Hodge or RH acceptance modifications.
