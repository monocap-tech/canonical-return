# CR-ALP-013 — Dimension-normalized cyclic traces converge; phase does not follow without a margin

**Canonical-Return / Alpöge branch · 2026-11-23 PST**  
**Status:** SOURCE-CONTRACT INSPECTED; DERIVED quantitative finite-Hilbert lemma; independent exact rational/SymPy audit **EXECUTED LOCALLY: PASS**. **No new Lean theorem, actual numeric arithmetic three-word instance, global signed Weil positivity, arithmetic-to-Weil operator equivalence, or independent RH result.**

## 0. Objective and immutable source pins

CR-ALP-012 established a mixed-complement HS bound from two **separate** arithmetic source contracts: the full-period diagonal near-isometry under \`LinkS_at\`, and the exponentially small *relative ramified rank loss* under \`RamShareF\`. CR-ALP-013 asks exactly what remains after division by the rapidly growing source dimensions. It separates **convergence of a normalized complex number** from **preservation of its argument**.

Directly inspected source: the frozen public Lean mirror \`CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4\`:

- \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`, around lines 145–160: \`hsDiagFull\`, \`LinkS_at\`, \`LinkS\`. The full, unprojected diagonal HS-square defect satisfies \(\|KD_{uu}-I\|_{\rm HS}^2\le e_0d_u\), while distinct-word \(\|KD_{uv}\|_{\rm HS}^2\le(1+\eta_0)d_ud_v/D_*\) **only if the displayed dimension-comparability antecedent holds**. Asymptotically \(\sqrt{e_0}\le e^{-x^{\varepsilon/8}}\) and \(\eta_0\le e^{-x^{\varepsilon/8}}\).
- \`src/Formalism/SharpRestriction/GlobalSchurBounds/Chosen.lean\`, around lines 55–60: \`ChProviso_at\` supplies that comparability, \(C_0^{-1}\le d_u\Gamma^2/D_*\le C_0\), with \(\Gamma=\mathrm{GamL}(L)>0\).
- \`src/Formalism/SharpRestriction/GlobalSchurBounds/TopDensityNGeSixteen.lean\`, around lines 858–1000: \`RamShareF\`, \`RamRank_at\`, \`ramShareF_of\`. The first gives \(q_u\le d_ed_u\) and \(\sqrt{d_e}\le e^{-x^{\varepsilon/8}}\), with **nonnegative** \(d_e\). The second gives the rank ratio \(r_u/d_u=((T_u-2)/(T_u-1))^f\). \`ramShareF_of\` uses \`RamRankF\` and \`TuSize\`, not \`LinkS\` alone.
- \`src/Formalism/SharpRestriction/GlobalSchurBounds/LinkMin.lean\`: \`ramRankF_min\` and \`linkGLOBAL7_items\` show the rank-loss chain is consumed by the conditional global frame reduction, under the named upstream arithmetic hypotheses.
- \`src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean\`, around lines 339–348, 509–525: \(K_{uv}=KD_{uv}\) is the **complete** normalized complex period cross block; \(R_u=\Pi_u\) is the orthogonal projection onto the ramified source.
- \`src/Formalism/LocalRingCardinality/DataParamsS.lean\`: \`SizeData₂.hsSq\` is exactly \`hsSqD\` of the full \(KD\).

The source contracts are statements/theorem uses read from Lean text, **not newly recompiled or audited for all imported assumptions here**. No actual three-word period-matrix instance has been extracted.

## 1. Dimension-normalized mixed-complement theorem

Use three distinct labelled words \(i,j,k\). Put
\[
d_u=\dim_{\mathbb C}C_u\ge1,\quad q_u=\operatorname{rank}(I-R_u),\quad
K_{uv}=A_u^\dagger A_v,\quad \widehat K_{uv}=R_uK_{uv}R_v
\]
with the canonical positive Kmat scaling folded into \(A_u\). Let
\[
\delta_{uv}=\|K_{uv}-\widehat K_{uv}\|_{\rm HS},\quad
e_u=\|K_{uu}-I\|_{\rm HS}.
\]
Assume \(D_*>0,C_0\ge1,\Gamma>0\), the *actual antecedent* of \`LinkS_at\` and \`ChProviso_at\`, and the simultaneous frame/ramified contracts
\[
\begin{aligned}
&0\le e_0,\ d_e,\quad e_u^2\le e_0d_u,\quad
q_u\le d_ed_u,\\
&\|K_{uv}\|_{\rm HS}^2
\le \frac{(1+\eta_0)d_ud_v}{D_*}
\quad (u\ne v),\\
&\frac1{C_0}\le \frac{d_u\Gamma^2}{D_*}\le C_0,\qquad 1+\eta_0\ge0.
\end{aligned}
\tag{13.1}
\]
The nonnegativity of \(1+\eta_0\) is automatic with three distinct positive-dimension words satisfying the offdiagonal inequalities, but is explicit here for clarity. Define
\[
\boxed{\beta=\min\!\left(1,\frac{C_0\Gamma^2}{D_*}\right).}
\]
Because \(d_u\) is a positive **integer**, \(1/d_u\le1\); comparability also gives \(1/d_u\le C_0\Gamma^2/D_*\). Thus \(1/d_u\le\beta\).

The full-complement bound from CR-ALP-012 gives after dividing by \(d_ud_v\):
\[
\begin{aligned}
\frac{\delta_{uv}^2}{d_ud_v}
&\le(d_e+\sqrt{d_ee_0})
 \left(\frac1{d_u}+\frac1{d_v}
 +\sqrt{\frac{e_0}{d_u}}+\sqrt{\frac{e_0}{d_v}}\right)\\
&\le\boxed{2(d_e+\sqrt{d_ee_0})(\beta+\sqrt{e_0\beta})}.
\end{aligned}
\tag{13.2}
\]
The estimate uses **both** full diagonal control and the relative-rank theorem; an individual one is insufficient.

Introduce the *labelled dimension-normalized cyclic trace*
\[
\boxed{
\mathscr C_{ijk}(K)=
\frac{D_*}{d_id_jd_k}
\operatorname{tr}(K_{ij}K_{jk}K_{ki}),}
\tag{13.3}
\]
and analogously \(\mathscr C_{ijk}(\widehat K)\). It is invariant under independent wordwise unit-modulus gauges. The three-term telescoping identity, trace Hölder inequalities and \(\|\widehat K_{uv}\|_{\rm HS}\le\|K_{uv}\|_{\rm HS}\) imply
\[
\boxed{
|\mathscr C_{ijk}(K)-\mathscr C_{ijk}(\widehat K)|
\le (1+\eta_0)
\sum_{\rm cyclic}\frac{\delta_{uv}}{\sqrt{d_ud_v}}.}
\tag{13.4}
\]
Indeed, in each term the three \(\sqrt{d_ud_v}\) dimension factors multiply to \(d_id_jd_k\), cancelling the normalization; the other two HS upper bounds supply \(1/D_*\).

Set \(\sigma=\exp[-x^{\varepsilon/8}]\), with \(0<\sigma\le1\). Under **both** eventual family bounds \(\sqrt{e_0}\le\sigma\), \(\sqrt{d_e}\le\sigma\), \(\eta_0\le\sigma\), (13.2) yields
\[
\frac{\delta_{uv}}{\sqrt{d_ud_v}}
\le 2\sigma\sqrt{\beta+\sigma\sqrt\beta}.
\]
Therefore
\[
\boxed{
|\mathscr C_{ijk}(K)-\mathscr C_{ijk}(\widehat K)|
\le 6(1+\sigma)\sigma\sqrt{\beta+\sigma\sqrt\beta}
\le12\sqrt2\,\sigma.}
\tag{13.5}
\]
**Positive result:** The two *dimension-normalized complex cyclic traces* converge together at an explicitly exponentially small, dimension-uniform **absolute** error, when \`LinkS\`, \`ChProviso\` and \`RamShareF\` hold on the same family and scale. This is a new algebraic consequence of those source contracts; it is **not** a Lean-certified theorem.

### 1a. Beware the different spectral normalization

The same offdiagonal HS estimates alone imply
\[
\boxed{|\mathscr C_{ijk}(K)|,\ |\mathscr C_{ijk}(\widehat K)|
\le\frac{(1+\eta_0)^{3/2}}{\sqrt{D_*}}.}
\tag{13.6}
\]
Thus if \(D_*\to\infty\), the quantities whose absolute difference is small may **both** tend to zero. An estimate of the form \(|z-\widehat z|\to0\) says nothing about their limiting arguments when \(|\widehat z|\to0\).

An alternative scale with order-one maximum,
\[
\mathscr B_{ijk}=\frac{D_*^{3/2}}{d_id_jd_k}\operatorname{tr}(K_{ij}K_{jk}K_{ki})
=\sqrt{D_*}\,\mathscr C_{ijk},
\]
satisfies only
\[
\boxed{
|\mathscr B(K)-\mathscr B(\widehat K)|
\le6(1+\sigma)\sigma\sqrt{D_*}\,
\sqrt{\beta+\sigma\sqrt\beta}.}
\tag{13.7}
\]
A limit of zero on this stronger scale requires an **extra growth inequality** connecting \(D_*,\Gamma,\sigma\); it has **not** been established from the contracts audited in this phase. This is why the distinction between *absolute normalized convergence* and a *relative phase-preserving estimate* matters.

## 2. Rigorous phase-margin gate, not phase promotion

If \(\mathscr C_{ijk}(\widehat K)\ne0\), then a sufficient phase criterion is
\[
\boxed{
6(1+\sigma)\sigma\sqrt{\beta+\sigma\sqrt\beta}
<|\operatorname{Im}\mathscr C_{ijk}(\widehat K)|}
\tag{13.8}
\]
to preserve the sign of the imaginary cyclic component; replace \(\operatorname{Im}\) by \(\operatorname{Re}\) for a real sign test. More generally an error strictly smaller than \(|\mathscr C(\widehat K)|\) keeps the full cycle nonzero and constrains its angular discrepancy.

**No lower bound on the reference oriented cyclic component is present in the inspected \`LinkS_at\`, \`RamShareF\` or \`ChProviso\` source contracts.** They bound magnitudes from above. The comparison to CR-ALP-007's *negative real signed Weil triangle* cannot supply this arithmetic margin without first constructing the very signed numerator-and-metric map at issue; using the Weil answer as a premise would be circular.

## 3. Exact family falsifier: relative rank tends to zero while full orientation flips

The absence of a margin is substantive, not a technicality. The companion regression checks the following **abstract finite positive-frame family**.

Let \(p>2\) range through odd primes. Set
\[
f=2,\qquad
d=(p-1)^2,\quad r=(p-2)^2,\quad q=2p-3,\quad
\frac qd=1-\left(\frac{p-2}{p-1}\right)^2\longrightarrow0.
\tag{13.9}
\]
This exactly reproduces the **numeric shape** of \`RamRank_at\` at \(T=p\), with f even and a **nonzero proper** complement. Take three source spaces \(C_i=\mathbb C^d\), common target \(\mathcal H=\mathbb C^{3d}\), \(D_*=D_{\rm all}=3d\), \(\Gamma=\sqrt3\), \(C_0=1\), and \(\eta_0=\eta_1=e_0=0\).

Construct a positive global Gram by using one ramified source coordinate per word with three-vector Gram
\[
F=\begin{pmatrix}1&t&t\\t&1&t\\t&t&1\end{pmatrix},
\quad t=1/10,
\]
and one **complement** coordinate per word with Gram
\[
G_s=
\begin{pmatrix}
1&a&b\\a&1&s\,ic\\b&-s\,ic&1
\end{pmatrix},
\quad a=1/5,b=1/6,c=1/7,\quad s=\pm1.
\]
All remaining \(r-1\) ramified and \(q-1\) complementary coordinates are orthogonal across words and to these two active sectors. Positivity of \(F,G_s\) gives a genuine Gram factorization into \(B_i^{(s)}:C_i\to\mathcal H\). The full global Gram is unitarily permuted to \(F\oplus G_s\oplus I_{3(d-2)}\). Consequently **both signs have identical complete frame spectra**.

Every full diagonal \(K_{ii}=I_d\), and **every restricted ramified mixed Gram is independent of \(s\)**, with its nonzero cyclic trace \(t^3=1/1000\). The entire offdiagonal HS-square data are equal for the two signs:
\[
\|K_{ij}\|_{\rm HS}^2=t^2+|(G_s)_{ij}|^2\le1/20
\le d^2/D_*=d/3.
\]
The comparability antecedent is exact:
\(d\,\Gamma^2/D_*=1\). The relative rank deficit \(d_e=q/d\to0\); for **any** \(\sigma>0\), choosing sufficiently large odd prime \(p\) gives \(\sqrt{d_e}\le\sigma\). No statement about an actual \`GamOf\` or the global \`ChosenCode\` count is intended.

Nonetheless, for all \(p\),
\[
\boxed{
\operatorname{tr}(K_{01}^{(s)}K_{12}^{(s)}K_{20}^{(s)})
=t^3+s\,\frac i{210},
\qquad
\operatorname{tr}(\widehat K_{01}\widehat K_{12}\widehat K_{20})=t^3.}
\tag{13.10}
\]
The full complex phase has **opposite nonzero imaginary orientation** for \(s=+1\) and \(-1\), even though the relative complement shrinks arbitrarily and the restricted data, full diagonal isometries, offdiagonal HS², rank law, target/source dimensions and *entire spectral moments* agree.

After the normalization (13.3),
\[
\mathscr C^{(s)}=\frac3{d^2}\left(t^3+s\,\frac i{210}\right),
\quad
\widehat{\mathscr C}=\frac{3t^3}{d^2}.
\]
Their difference tends to zero, **exactly as (13.5) predicts**, while their phases do not converge to one another. The ratio \(\operatorname{Im}\mathscr C^{(s)}/\operatorname{Re}\mathscr C^{(s)}=s\,100/21\) is independent of \(p\). This is a rigorous countermodel to a phase inference from the indicated **numerical** source contracts alone, **not a counterexample to the complete actual arithmetic D7/lattice/divisibility construction**.

## 4. Reproducibility, disposition, next gate

**Executed local audit:** \`src/audits/cr_alp_013_normalized_cycle_audit.py\`, using Python \`Fraction\` and SymPy exact algebra. It checks positive \(F,G_s\), complex cycles, rank law at \(p=3,5,11,101,1009\), offdiagonal LinkS-shaped bounds, comparability, complementary HS inequalities, normalized cycle orientation, and the (13.5) universal bound for exact rational \(\sigma,\beta\) samples. **PASS**. No new GitHub Actions or Lean compiler result is asserted.

| Obligation | Status |
|---|---|
| Quantitative dimension-normalized \(\mathscr C(K)-\mathscr C(\widehat K)\) | **DERIVED / ANALYTIC PASS** |
| Uniform explicit bound \(12\sqrt2\,e^{-x^{\varepsilon/8}}\) | **DERIVED from joint source hypotheses** |
| Source rank law / near-isometry / ChProviso dictionary | **SOURCE-INSPECTED** |
| Phase without an independent nonzero arithmetic margin | **NOT IMPLIED** (exact countermodel) |
| Natural order-one spectral normalization | **GROWTH-DEPENDENT; OPEN** |
| Concrete substantive D7 with evaluated three-word \(KD\) | **NOT EXTRACTED** |
| Original signed Weil metric-and-numerator intertwiner | **OPEN** |
| New Lean theorem / RH proof | **NOT CLAIMED** |

**Next — CR-ALP-014:** inspect the specific \`GamOf\`, \(D_*\), source-dimension and window relations against the stronger-scale budget (13.7), and search the full period/Bessel/lattice source for a **nonzero labelled cyclic phase margin**, not only rank or HS-square bounds. If no arithmetic phase margin is source-defined, stop promoting any full operator correspondence: the scalar arithmetic-to-smooth-original-Weil conditional test-family projection survives independently of this phase gate.

**Governance:** research-only Canonical-Return draft; no changes to Weil-Lab, provenance or Hodge acceptance ledgers.
