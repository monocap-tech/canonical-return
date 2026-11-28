# CR-ALP-015 — Source-defined cyclic-margin audit and a certificate-level stop

**Research jurisdiction:** Canonical-Return / independent Alpöge projection, 2026-11-23 PST.  
**Disposition:** SOURCE-INSPECTED; new finite quantitative implication DERIVED; standalone exact symbolic/rational audit **PASS**. **Full complex arithmetic-to-Weil operator transfer: NOT ESTABLISHED; current phase-certificate strategy: STOP.** No new Lean compilation, no evaluated arithmetic `D7` triple, no RH theorem or global original Weil positivity theorem.

## 0. Why this is a completion gate rather than another analogy

CR-ALP-004 gave a **conditional scalar analytic map** from the arithmetic prime ramp to a family of full, signed, smooth Weil test values. CR-ALP-010–014 attempted the stronger claim that arithmetic period operators transport **native metrics, mixed complex phases, and eventually the complete signed Weil numerator**. The present stage tests the only defensible outstanding phase criterion: whether an actual arithmetic **labelled** cyclic component is independently bounded away from zero by more than the proven truncation/complement error.

Evidence is limited to the following directly inspected public Lean source at frozen **`CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4`**:

| Exact source path/declaration | What the contract actually provides | What it does *not* state |
|---|---|---|
| `src/Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean` / `D7.BuRational` | `Bu_u(jN_u φ, jEE x) = embLL(Bk_u φ x)` on the indicated number-field vector spaces | A nonzero labelled three-word cyclic trace; a denominator, height, or phase-margin bound on the **full normalized `KD`** |
| `src/Formalism/SharpRestriction/GlobalTheoremInputs.lean` / `BesselId` | `BuRational` plus **one common target linear isometry** `V`, wordwise linear isometries `U_u` and nonzero scalars `a_u` with `Bu_u = a_u besS_u(U_u -, V -)` | The nonvanishing of `Bu_u` itself or of any cyclic product; the exact value of the full-period cyclic trace |
| Same file / `LinkDIV_at` | Arithmetic lattice, heights, and a **per-word conditional** `Bk_u ≠ 0 → htMap(Bk_u) ≤ h_u`, as well as `BesselId` | An identity turning *three* Bk maps into a labelled complex cross cycle; a nonzero real/imaginary component |
| Same file / `LinkS_at`, `LinkS` | Upper bounds on distinct-word **squared HS** overlaps, full diagonal defect, and full-target dimension, **subject to the displayed dimension-comparability antecedent** | A third-order labelled cyclic trace or its lower bound |
| `src/Formalism/SharpRestriction/GlobalSchurBounds/Chosen.lean` / `GamL`, `ChProviso_at` | Native Gamma and source-dimension comparability | An oriented phase or a relative phase-error bound |
| `src/Formalism/SharpRestriction/GlobalSchurBounds/TopDensityNGeSixteen.lean` / `RamShareF`, `MatS_at`, `MatSF` | Complementary-rank bound and frame/rank/slope reductions based on HS second moments | A specified triple's complex phase |

These are **source-contract limitations**, not a proof that no other source file or the full arithmetic construction could imply a margin. The source is not independently recompiled in CR. The relevant original research context is Alpöge's *A remark* (20 Nov 2026) and Rubinstein et al.'s working reconstruction v2 (22 Nov 2026); neither is substituted here for the frozen Lean declarations.

## 1. New mathematical result: the sigma-only phase *certificate* is asymptotically larger than the largest possible cycle

Use three distinct words `i,j,k`. Denote their full normalized mixed operator blocks by `K_uv`, their ramified restrictions by `Khat_uv=R_u K_uv R_v`, and the native source ranks by `d_u≥1`. Let

\[
\sigma=e^{-x^{\varepsilon/8}}\in(0,1],\qquad
\Gamma=\exp(c\,m^5 f\,x^\varepsilon),\qquad c>0,\quad m\ge2,\ f\ge1.
\]

With `Dstar>0`, `C0≥1` and the *simultaneous* `LinkS`, `ChosenProviso`, `RamShareF` hypotheses, set

\[
\beta=\min\left(1,\frac{C_0\Gamma^2}{D_*}\right),\qquad
\mathscr C(K)=\frac{D_*}{d_id_jd_k}\operatorname{tr}(K_{ij}K_{jk}K_{ki}).
\]

CR-ALP-013's dimension-sensitive but **sigma-only** completion-error *upper budget* is

\[
\boxed{E_{13}=6(1+\sigma)\sigma\sqrt{\beta+\sigma\sqrt\beta},\qquad
|\mathscr C(K)-\mathscr C(\widehat K)|\le E_{13}.}\tag{15.1}
\]

The **same source** upper HS bounds force a maximum envelope on the ramified reference:

\[
\boxed{|\mathscr C(\widehat K)|\le
M:=\frac{(1+\sigma)^{3/2}}{\sqrt{D_*}}.}\tag{15.2}
\]

By `ChProviso_at`, every source dimension obeys \(d_u\Gamma^2/D_*\ge1/C_0\). Because \(d_u\ge1\),

\[
D_*\ge\frac{\Gamma^2}{C_0},\qquad
D_*\beta=\min(D_*,C_0\Gamma^2)\ge\frac{\Gamma^2}{C_0}.
\]

Therefore, without any unproved arithmetic phase hypothesis,

\[
\boxed{
\frac{E_{13}}{M}
\ge\frac{6\Gamma\sigma}{\sqrt{C_0(1+\sigma)}}
\longrightarrow+\infty.}\tag{15.3}
\]

Indeed \(\Gamma\sigma\ge\exp(32c\,x^\varepsilon-x^{\varepsilon/8})\to\infty\), by the *actual native Gamma definition* and positive family constant. For all sufficiently large `x`, the **known error upper budget** `E13` therefore exceeds **even the largest mathematically permitted** \(|\mathscr C(\widehat K)|\).

**Consequent certificate-level STOP:** No separate nonzero real, imaginary, or modulus margin \(\mu\le|\mathscr C(\widehat K)|\) can satisfy the *sufficient test* \(E_{13}<\mu\) at large scales. This is stronger than merely saying a margin was not found: **a margin alone cannot rescue this particular sigma-only error certificate**. It does **not** say the actual unknown error exceeds the true cycle; a sharper *parameter-sensitive* or phase-specific estimate could still succeed.

### 1a. The rank-integrality improvement does not change this verdict

Under the **substantive proper-complement rank law** \(1\le q_u\le\sigma^2d_u\), CR-ALP-014 sharpened the **absolute** budget to

\[
E_{14}=12\sqrt{C_0^2+C_0}\,\sigma^2.
\]

Now \(d_u\ge\sigma^{-2}\), whence \(D_*\ge\Gamma^2/(C_0\sigma^2)\) and

\[
\boxed{
\frac{E_{14}}M
\ge\frac{12\sqrt{C_0+1}\,\Gamma\sigma}{(1+\sigma)^{3/2}}
\longrightarrow+\infty.}\tag{15.4}
\]

Again this concerns the **upper-bound expression**. The best independent spectral-envelope estimate gives only \(|\mathscr B(K)-\mathscr B(\widehat K)|\le2(1+\sigma)^{3/2}\), which also cannot distinguish a nonzero cyclic sign without more structure.

**Important methodological restriction:** (15.3)–(15.4) rule out phase certification using the *coarse eventual `sigma` bounds* as assembled in CR-ALP-013/014. They do **not** rule out an argument using the actual per-word `e0`, `de`, dimensions, cancellation between mixed blocks, arithmetic integrality, or an independently specified exact period identity. Those must be proved and tested separately.

## 2. Rationality + common Bessel gauge + exact diagonal isometry still do not force a cyclic margin

The absence of a nonzero cyclic margin in the inspected source statements has an explicit **finite, rational-coefficient** manifestation. This example is a *model of the pairing/finite-frame schemas*, **not** an instantiated `D7` or a model of every arithmetic level/divisibility assumption.

Take three labelled source lines \(C_u=\mathbb C\), target \(\mathcal H=\mathbb C^3\), and **one fixed field** \(L=\mathbb Q(i)\) with standard full coordinate lattices. For every integer \(n\ge4\), let

\[
\begin{aligned}
v_0&=(1,0,0)^T,\\
v_1&=(5/13,12/13,0)^T,\\
v_2^{\pm}(n)&=\left(7/25,\ \pm i\frac{24}{25}\frac{2n}{n^2+1},\ \frac{24}{25}\frac{n^2-1}{n^2+1}\right)^T.
\end{aligned}\tag{15.5}
\]

Every vector is **unit** and has all coefficients in the *same* fixed \(\mathbb Q(i)\); the three are independent. Set \(B_u(z)=zv_u\), \(K_{uv}=B_u^*B_v\), and the bilinear pairing
\(\mathrm{Bu}_u(\phi,F)=F^Tv_u\phi\). If `besS` is represented by the same bilinear coefficients, the **pairing-equality component** of `BesselId` is realized with `a_u=1`, `U_u=V=I`, while the algebraic map `Bk_u` is given by the identical \(\mathbb Q(i)\)-bilinear coefficients. Thus the **form** of `BuRational` and the source's Bessel attachment are satisfied on the chosen algebraic vector spaces; **each `Bk_u` is nonzero**.

All native source diagonals are **exactly** 1; all pairwise squared overlaps obey \(|K_{uv}|^2<1/3\). Consequently the **numerical** generic `InputS`/`LinkS_at` frame inequalities hold with `d_u=1`, `Dstar=Dall=3`, `Gamma²=3`, `C0=1`, `e0=eta0=eta1=0`. This Gamma is illustrative; it is **not the exponential native `GamOf` for a `D7` family**.

But their gauge-invariant labelled three-cycle is

\[
\boxed{
\operatorname{tr}(K_{01}K_{12}K_{20})
=\frac{49}{4225}\ \pm\ i\frac{4032n}{21125(n^2+1)}.}\tag{15.6}
\]

The imaginary part is strictly nonzero for each finite `n`, and tends to zero as \(n\to\infty\), **despite one fixed number field, three nonzero algebraic pairings, common Bessel target gauge, exact diagonal isometry, and uniform overlap upper bounds**. Complex conjugation exchanges the signs without changing the entire frame spectrum. No uniform lower bound on the labelled oriented imaginary part can be inferred from this *finite source-schema information alone*.

There is a stronger elementary degeneracy: the entirely rational triples with
\(v_2=(7/25,0,24/25)^T\) and \(v_2=(7/25,-72/125,96/125)^T\) have **nonzero positive** and **nonzero negative** real cyclic traces, respectively, while preserving exact diagonal isometry and the same numerical HS upper bound. Three orthogonal unit vectors instead yield exactly zero cyclic trace, although all three individual `Bk` coefficients are nonzero. The recorded exact audit checks all cases.

**Rationality is not a quantitative phase margin.** If *all actual full* \(K_{uv}\) entries were known to lie in \(N^{-1}\mathbb Z[i]\), then a **nonzero** imaginary triple trace would obey \(|\operatorname{Im}\mathcal T|\ge N^{-3}\). But `BuRational` does not assert such a denominator bound for the *normalized full period kernels*: its equality applies to `Bu/Bk` on specified number-field input spaces, while `KD` includes Hilbert-adjoint, period and positive metric normalizations. Even a certified bound on `N` would still need an independent proof of **nonvanishing** of the relevant real/imaginary part. `LinkDIV_at` provides some **conditional per-word height information** for nonzero `Bk`, not a labelled triple nonvanishing/denominator transfer.

## 3. Decision, precise reopening conditions, and next productive research lane

**Source-specific positive statement retained:** `BesselId` preserves the normalized cyclic phase between the **restricted** `Bu` and `besS` carriers, under its pairing-equality witnesses. The full period-vs-ramified phase error has explicit upper bounds. The source uses native isometries onto `ramCu` (CR-ALP-011 correction).

**Certified limitation:** The joint `LinkS` + `ChosenProviso` + `RamShareF` `sigma`-only bounds **cannot certify** a nonzero full cyclic sign/phase at the native Gamma scale, even with a hypothetical maximal reference margin. They do **not** refute an actual arithmetic phase identity or a stronger arithmetic estimate.

**No new source extraction:** No explicit **evaluated, nondegenerate three-word** arithmetic `D7` instance and no independently sourced nonzero labelled cyclic phase component were extracted. This is a bounded source inspection, not a theorem of absence throughout the repository. The actual `BuRational` and lattice-divisibility contracts do not currently attach the *complete original signed Weil numerator* with its prime-power, pole and Gamma terms to a native-metric arithmetic source map.

**Disposition of the attempted correspondence:**

| Lane | State as of CR-ALP-015 |
|---|---|
| Original signed Weil prime-ramp scalar identity | **RETAINED** (CR-ALP-001/002) |
| Conditional quantitative smooth signed Weil family | **RETAINED, CONDITIONAL** on arithmetic estimate and documented smooth analysis (CR-ALP-004) |
| Internal ramified `Bu ↔ besS` cyclic-phase covariance | **DERIVED** under explicit BesselId pairing-equality component (CR-ALP-010) |
| Dimension-normalized full/ramified cyclic absolute difference | **DERIVED, SOURCE-CONDITIONAL** (CR-ALP-013/014) |
| Current coarse-error budget certifies nonzero phase | **STOP: IMPOSSIBLE FOR THIS BUDGET AT LARGE SCALE** |
| Actual full theta-lift-to-original-Weil `H/M` intertwiner | **NOT ESTABLISHED** |
| Global Weil positivity, RH proof, newly compiled Lean result | **NOT CLAIMED** |

**Reopening conditions for the full operator lane:** (i) an actual evaluated arithmetic **three-word** `D7` instance with native `Pu`, `Bu`, `KD` and lattice/embedding normalizations; (ii) a **source-derived nonzero labelled real or imaginary cyclic margin**, or an exact structural phase equation; (iii) an **improved error bound below that margin at the correct normalization**, not the sigma-only upper budget; and (iv) separately, a map preserving **both** the original signed Weil numerator `H` and physical logarithmic `L²` metric `M` with its admissibility and support growth. Without these, phase analogy is not operator equivalence.

**Suggested next independent work: CR-ALP-016, analytic-return gate.** Focus on the legitimate scalar ramp-to-original-Weil route, including source-bound admissibility, two-pole moments and support-extension obstructions, rather than generating another abstract frame similarity. Continue F41/F42 Formalization and Phase IV under their own ledgers; **no governance promotion** from this CR draft.

## 4. Reproducibility

Companion file: `src/audits/cr_alp_015_phase_margin_exact_audit.py` (Python standard library + SymPy). **Executed locally: PASS.** Tests exact `Q(i)` bilinear coefficients, nonzero per-word Bk models, positive Gram determinants, all offdiagonal HS-square bounds, opposite complex-conjugate cycles, cycles of each real sign and zero, exact diminishing phase margins at `n=4,5,10,101,1000,10000`, denominator separation, and squared inequalities comparing source-derived error **budgets** to their maximal reference-cycle **envelopes**. This is a mathematical finite-algebra regression, **not** Lean/CI acceptance, evaluation of the actual arithmetic theta-lifts, or an RH certificate.