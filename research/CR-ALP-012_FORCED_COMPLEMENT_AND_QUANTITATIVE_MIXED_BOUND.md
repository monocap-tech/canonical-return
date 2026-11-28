# CR-ALP-012 — Source rank law forces a proper ramified complement; quantitative mixed-block control

**Project:** Canonical-Return, independent Alpöge projection, draft PR #3.  
**Date:** 2026-11-23 PST.  
**Standing:** directly SOURCE-INSPECTED formal contracts; new ANALYTIC conditional estimates; EXACT FINITE REGRESSION PASS executed locally. **Not a new Lean compile, evaluated arithmetic period triple, global Weil result, or independent RH proof.**

## 0. A correction to the CR-ALP-011 fork

CR-ALP-011 left two possibilities for the substantive period family: perhaps the Bessel-visible ramified subspace was the entire full period source, or perhaps it had a complement. The actual **rank formula used in the global arithmetic argument** now settles that fork **under its explicit large-family hypotheses**: for nonzero source rank, exponent \(f\ge1\) and auxiliary level \(T>2\), the complement is **necessarily proper and nonzero**. The later framework also bounds its **relative** rank exponentially, which is a positive quantitative input not exploited in the earlier no-go examples.

This does *not* identify a concrete numerical multiword D7 family. The theorem schemas and proof uses have been directly read at the public frozen mirror.

## 1. Pinned source declarations

Frozen public Lean mirror: CocoSgt/riemann-hypothesis-lean @ **7cbf7a0c32fca1badf05931879a37c3f8a2b03b4**.

- src/Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean: D7.duS = complex finrank of Mu0; D7.eMu_isom and D7.s_Cu_isom are the canonical isometry conditions. src/Formalism/LocalRingCardinality/DataForms.lean: A.du is the complex finrank of Mram.
- src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean, around lines 509–525: ramCu, its orthogonal projection PiU, full normalized Kmat/KD, projected diagonal hsDiagD.
- src/Formalism/SharpRestriction/GlobalTheoremInputs.lean, around lines 143–160: LinkS_at includes the **full** diagonal HS-square error and off-diagonal HS² estimates, conditional on dimension comparability; LinkS wraps these into asymptotically small e0, eta0, eta1.
- src/Formalism/SharpRestriction/GlobalSchurBounds/TopDensityNGeSixteen.lean, around lines 891–1000: **RamRank_at**, **RamRankF**, **ramShare_at_of**, **ramShareF_of**, and FrameRealF_of''. The exact family rank law and relative-defect bound appear there.
- src/Formalism/SharpRestriction/GlobalSchurBounds/LinkMin.lean, around lines 659–664 and the LinkGLOBAL7 chain: **ramRankF_min** consumes J.LvlF and supplies RamRankF. The substantive frame reduction **uses this source rank formula**, not only a hypothetical toy example.
- src/Formalism/SharpRestriction/GlobalSchurBounds/JointsL.lean: two_lt_level_prime under the appropriate arithmetic setting. src/Formalism/SharpRestriction/GlobalSchurBounds/Chosen.lean: ChosenProviso contains positive/comparable per-word source dimensions.
- src/Formalism/SharpRestriction/GlobalSchurBounds/TopDensityNGeSixteen.lean: MatKF, RamShareF, MatSF, and matSF_of combine source-frame data, relative rank loss and the generic InputSG formulation.

No full end-to-end theorem compilation or construction of explicit numerical period matrices was performed for CR.

## 2. Exact rank formula and forced proper complement

For a D7 with RamRank_at, denote
\[
d_u=\dim_\mathbb C(C_u)=\mathrm{duS}(u),\quad
r_u=\dim_\mathbb C(\mathrm{ramCu}_u)=A_u.\mathrm{du},\quad
q_u=d_u-r_u.
\]
Source compatibility (eMu_isom, s_Cu_isom) identifies the ramified rank exactly with A.du; this is a **native metric**, not an arbitrary identification.

RamRank_at states
\[
\boxed{
r_u=d_u\left(1-\frac1{T_u-1}\right)^f
=d_u\left(\frac{T_u-2}{T_u-1}\right)^f.}
\tag{12.1}
\]
Consequently, **assuming \(d_u>0,\ f\ge1,\ T_u>2\)**, one has
\[
\boxed{
0<q_u=d_u\left[1-
\left(\frac{T_u-2}{T_u-1}\right)^f\right]<d_u.}
\tag{12.2}
\]
Since \(q_u\) is an integer, \(q_u\ge1\). Thus the full period source cannot equal the Bessel-visible ramified subspace **in the nonzero-rank, nonexceptional-level portion of this rank-law family**. The all-ramified hostLevelRecord toy from CR-ALP-011 has its own construction and should not be substituted for this family.

Bernoulli's inequality gives a simple relative defect bound:
\[
\boxed{
0<\frac{q_u}{d_u}
=1-(1-1/(T_u-1))^f
\le\frac f{T_u-1}.}
\tag{12.3}
\]
The source theorem ramShare_at_of obtains the coarser \(q_u/d_u\le 2f/b\) with its size window, and ramShareF_of packages this at large \(x\) as
\[
\boxed{
q_u\le d_e\,d_u,\qquad
d_e\ge0,\qquad
\sqrt{d_e}\le\exp[-x^{\varepsilon/8}].}
\tag{12.4}
\]
Here \(d_e\) is a **relative dimension-loss parameter**. In particular \(d_e\le e^{-2x^{\varepsilon/8}}\) when the source hypotheses apply. Since \(q_u\ge1\), (12.4) entails \(d_e d_u\ge1\): **exponentially small relative deficit need not mean a vanishing absolute complement**. No inference that every source is eventually all-ramified is justified.

## 3. New quantitative mixed-complement estimate from the full diagonal input

Let \(C_u\) be finite-dimensional Hilbert spaces, \(P_u:C_u\to\mathcal H\), \(c_u>0\) be the source Kmat scaling, and set
\[
A_u=c_u^{-1/2}P_u,\qquad K_{uv}=A_u^*A_v.
\]
Let \(R_u\) be the genuine orthogonal projection onto ramCu and \(Q_u=I-R_u\), with \(q_u=\operatorname{rank}Q_u\). Suppose the **full** diagonal defect obeys
\[
\|K_{uu}-I\|_{\mathrm{HS}}\le e_u.
\]
Write \(\widehat K_{uv}=R_uK_{uv}R_v\),
\(\Delta_{uv}=K_{uv}-\widehat K_{uv}\), and
\(\delta_{uv}=\|\Delta_{uv}\|_{\mathrm{HS}}\).

A direct Hilbert-space argument yields the **uniform dimension-explicit estimate**
\[
\boxed{
\delta_{uv}^2\le
(1+e_v)(q_u+\sqrt{q_u}\,e_u)
+(1+e_u)(q_v+\sqrt{q_v}\,e_v).}
\tag{12.5}
\]

**Derivation:** \(\|A_u\|_{\mathrm{op}}^2
=\|K_{uu}\|_{\mathrm{op}}\le1+e_u\). Also
\[
\|A_uQ_u\|_{\mathrm{HS}}^2
=\operatorname{tr}(Q_uK_{uu}Q_u)
\le q_u+\sqrt{q_u}\,e_u
\]
by Hilbert–Schmidt Cauchy–Schwarz. The two omitted blocks
\(Q_uK_{uv}\) and \(R_uK_{uv}Q_v\) have orthogonal left source ranges. Their HS² add. Applying \(\|X^*Y\|_{\rm HS}\le\|X\|_{\rm HS}\|Y\|_{\rm op}\) to the first, and the reversed HS/op estimate to the second, gives (12.5). The estimate holds even when the ramified subspace is proper.

Now include the **actual shape** of LinkS_at and RamShareF:
\[
q_u\le d_e d_u,\qquad
e_u^2\le e_0d_u,\qquad e_0,d_e\ge0.
\]
Then
\[
q_u+\sqrt{q_u}e_u
\le(d_e+\sqrt{d_e e_0})d_u
\]
and \(e_u\le\sqrt{e_0d_u}\). Thus
\[
\boxed{
\delta_{uv}^2\le
(d_e+\sqrt{d_e e_0})
\left[d_u(1+\sqrt{e_0d_v})+
d_v(1+\sqrt{e_0d_u})\right].}
\tag{12.6}
\]
This is a **genuine quantitative conditional complement estimate**, not merely a missing hypothesis. It is deduced from the source-defined *full* diagonal near-isometry plus rank deficit. It does **not** automatically provide the scale needed to determine a cyclic phase, because \(d_u\) can be extremely large and the source provides no nonzero oriented cyclic margin.

### 3a. A precise phase-margin interface

For a labelled triple \(i,j,k\),
\[
\mathcal T(K)=\operatorname{tr}(K_{ij}K_{jk}K_{ki}),
\quad
\mathcal T(\widehat K)=\operatorname{tr}(\widehat K_{ij}\widehat K_{jk}\widehat K_{ki}).
\]
The telescoping/Schatten argument gives
\[
\boxed{
\begin{aligned}
|\mathcal T(K)-\mathcal T(\widehat K)|\le{}&
\delta_{ij}\|K_{jk}\|_{\rm HS}\|K_{ki}\|_{\rm HS}\\
&+\|\widehat K_{ij}\|_{\rm HS}\delta_{jk}\|K_{ki}\|_{\rm HS}\\
&+\|\widehat K_{ij}\|_{\rm HS}\|\widehat K_{jk}\|_{\rm HS}\delta_{ki}.
\end{aligned}}
\tag{12.7}
\]
Substituting (12.6) and LinkS_at's off-diagonal bounds produces an explicit computable error budget in \(d_i,d_j,d_k,D_*,e_0,d_e,\eta_0\). A source-certified sign/phase follows only if that error is **strictly smaller than** the appropriate nonzero reference real or imaginary cyclic margin. No such margin was found in the Lean rank/HS inputs.

**Interpretive correction:** CR-ALP-011's claim that full diagonal near-isometry does not *by itself* make the complement disappear is still correct. However, the **additional source rank law and RamShareF** give a rigorous *coarse quantitative mixed-complement bound*. The distinction between those statements must be preserved.

## 4. The actual Lean arithmetic branch and the remaining gap

We now know the *formal* rank/construction chain:
\[
\mathrm{J.LvlF}\xrightarrow{\mathrm{ramRankF\_min}}\mathrm{RamRankF}
\xrightarrow[\mathrm{TuSize}]{\mathrm{ramShareF\_of}}\mathrm{RamShareF},
\]
which feeds the frame reduction together with LinkS, MatKF and other hypotheses. This is **not** only a free-standing algebraic model; these exact theorem names occur in the source's global reduction (LinkMin.lean, linkGLOBAL7_items).

This settles the substantive source **fork**: when the arithmetic rank and setting hypotheses hold, the ramified complement is not zero. It is *relatively thin* but still required in the full period Gram; the appropriate next test is quantitative, not an all-ramified shortcut.

**Still missing:** (i) a fully instantiated three-word arithmetic operator family with complex mixed K values; (ii) a source-defined nonzero cyclic phase margin, or a structural phase identity; (iii) a verified estimate sufficiently sharp **after the actual growing frame dimensions and normalizations**; (iv) a separate map preserving the *signed original Weil numerator* \(H\), *physical metric* \(M\), prime/pole/Gamma corrections and expanding support. The conditional arithmetic-to-smooth-Weil test-family **reexpression** in CR-ALP-004 remains valid but is not an independent proof of the hard arithmetic estimate.

## 5. Reproducibility and gate

The attached standalone Python regression uses exact Fraction/SymPy calculations for three integer rank-law examples (including a case with small relative but large absolute complement) and the exact 6D frame of CR-ALP-011. It directly checks the HS² compression identity, zero full diagonal defect with nonzero complementary mass, the bound (12.5) on all blocks, and the opposite full cyclic orientations. **Executed locally: PASS**. This is not a Lean compile or a proof of the source-wide analytic estimate.

**Next — CR-ALP-013:** insert the *actual per-word source dimension/scale windows* into (12.6) and the off-diagonal LinkS bound to produce a **normalized phase-error estimate**. Determine whether the error can be smaller than a source-attached cyclic margin; if no such margin is available, classify full complex operator transport as **UNESTABLISHED** rather than extrapolating from moment/rank bounds. Keep the Weil scalar conditional analytical projection separate.

**Governance:** research-only CR draft; no writes to Weil-Lab, provenance or Hodge acceptance ledgers, no global Weil positivity or RH certificate.
