# CR-ALP-008 — Exact phase-blindness witness for the arithmetic frame input

**Canonical-Return / Alpöge projection** · 2026-11-23 PST.
**Standing:** SOURCE-INSPECTED / EXACT FINITE COUNTERMODEL / LOCAL RATIONAL AUDIT PASS. **No actual instantiated three-word theta-lift operator triple; no new Lean theorem; no RH, global Weil positivity, or operator-equivalence claim.**

## 0. Precise question and pinned primary evidence

Can the arithmetic **rank and Hilbert–Schmidt second-moment hypotheses** used in Alpöge's theta-lift proof determine a **gauge-invariant three-cell cyclic phase**, of the kind isolated in CR-ALP-007 for the original signed Weil form?

Primary paper: Levent Alpöge, *A remark*, 20 November 2026, https://alpo.ge/remark.pdf. Independent working exposition: Rubinstein et al., *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, v2, 22 November 2026, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, especially §III.8 (S), Lemma T1 / trace and rank bounds, §III.11.

Source-checked public Lean mirror: \`CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4\`. Exact file/declaration anchors:

- \`src/Formalism/SizeWindowBrackets/FrameTraceRatioFinEquiv.lean\`: \`MatrixLatticeEstimates.KpG\` = \(\sqrt{\gamma_i\gamma_j}B_i^*B_j\), \`KdG\`, \`QmG\`, and \`InputSG\` (off-diagonal HS-square bounds, diagonal defects, dimension comparison).
- \`src/Formalism/SizeWindowBrackets/TraceSquareAndDiagBounds.lean\`: \`HilbertSchmidtNormSquare.Kp\`, \`Qm\`, \`InputS\`, \`TraceIdentities\`, \`F1\`, \`F2\`; **different** \`InputS\` and \`InputSG\` have the same finite-frame shape.
- \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`: \`LinkS_at\`, \`LinkIRR\`, \`BesselId\`. The quantitative \`LinkS_at\` is **conditional** on its dimension-comparability antecedent and is not identical to the generic \`InputS\` instance below.
- \`src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean\`: \`Kmat\`, \`KD\`, \`hsSqD\`, with a complex cross-operator \((P_u)^*P_v\) and Hilbert–Schmidt-square evaluation.
- \`src/Formalism/LocalRingCardinality/DataRows.lean\`: word-indexed Hilbert source/target types and \`Rows.Data.Bu\`.
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/ExampleClassGroupCardPos.lean\`: concrete example \`exampleRecordOfData\` uses \`wordFin := inferInstanceAs (Fintype Unit)\`, so this **specific example contains only one word** and cannot furnish the required arithmetic triple.
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/RankRecordRowDimensionZero.lean\`: \`rankRecord\` also has \`wordFin := Fintype Unit\`.
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/NZRowE.lean\`: \`hostLevelRecord\` can take a *parameterized* finite word type \`W\`, with arbitrary prescribed pairing input; the definition does **not** supply three explicit evaluated period matrices.

These are direct source observations. Search/extraction has **not** established that no three-word arithmetic instance exists elsewhere; it has established that the inspected concrete one-word records cannot serve as such an instance.

## 1. The exact two-moment phase-blindness countermodel

For \(0<r<\frac12\) and \(s\in\{+1,-1\}\), define the **real Hermitian** matrices
\[
\boxed{
G_s(r)=
\begin{pmatrix}
1&r&r\\
r&1&s r\\
r&s r&1
\end{pmatrix}.}
\tag{8.1}
\]
All are strictly positive definite, e.g. by strict row diagonal dominance \(2r<1\). Hence each is a Gram matrix of **three genuine vectors** in a common complex Hilbert space.

For the explicit rational value \(r=\frac14\), one can choose target \(\mathcal H=\mathbb C^3\), source lines \(C_i=\mathbb C\), and columns
\[
B_0=(1,0,0)^\top,\quad
B_1=(1/4,\sqrt{15}/4,0)^\top,
\]
\[
B_2^{(+)}=(1/4,3/(4\sqrt{15}),3/\sqrt{10})^\top,\quad
B_2^{(-)}=(1/4,-5/(4\sqrt{15}),\sqrt{5/6})^\top.
\]
A direct dot-product check yields \(B_i^*B_j=(G_s(1/4))_{ij}\) exactly.

Choose all frame weights \(\gamma_i=1\). With the **typed** Lean \`Kp\` convention, the mixed cross block is the \(1\times1\) scalar \(K_{ij}=B_i^*B_j\), and the target frame is \(Q_s=\sum_i B_iB_i^*\).

These two frames have **identical**:

- diagonal native source Grams \(K_{ii}=1\), diagonal defects \(\|K_{ii}-1\|_{\rm HS}^2=0\);
- every off-diagonal squared HS overlap \(\|K_{ij}\|_{\rm HS}^2=r^2\);
- full target/source dimensions, positive ranks (=3), and first two trace moments;
- frame isotropy excess \(\Phi_3=2r^2\).

The exact identities are
\[
\boxed{\operatorname{tr}Q_s=3,\quad
\operatorname{tr}(Q_s^2)=3+6r^2,\quad
\Phi_3(Q_s)=2r^2,\quad
\det G_s=1-3r^2+2sr^3>0.}
\tag{8.2}
\]

**Explicit source-condition check:** for \(r=1/4\), with \(D=D_*=3\), \(\eta_0=\eta_1=e_0=0\), \(d_i=1\), \(\gamma_i=1\), both **satisfy the literal generic finite** \`InputS\`/\`InputSG\` inequalities: \(1/16\le 1/3\), diagonal defect \(0\le0\), and \(D\le D_*\). This is a **model of the finite frame input schema**, not an instantiated \`D7\` that establishes \`LinkS_at\`, \`LinkIRR\`, \`BesselId\`, or the deeper rank/slope theorems.

Yet the three-cell cyclic trace is
\[
\boxed{
\mathcal T_{012}=\operatorname{tr}(K_{01}K_{12}K_{20})=sr^3.
}
\tag{8.3}
\]
At \(r=1/4\),
\[
\boxed{\mathcal T_{012}^{(+)}=+1/64,\qquad
\mathcal T_{012}^{(-)}=-1/64.}
\tag{8.4}
\]

Independent line phase rotations \(B_j\mapsto e^{i\theta_j}B_j\) change \(K_{ij}\) by \(e^{-i\theta_i}e^{i\theta_j}\), but these cancel in the cyclic product. Thus (8.4) represents **genuinely different labelled gauge-invariant phase data**, despite all first/two-moment agreement. For every fixed target dimension one can choose arbitrarily small \(r>0\), making \(\Phi_3=2r^2\) arbitrarily small **without changing the opposite cyclic phases**. Their *magnitudes* then tend to zero; no fixed nonzero phase gap is asserted.

**Finite NO-GO lemma:** No map that takes **only** diagonal native source masses, pairwise HS overlap squares, rank, and \(\operatorname{tr}Q,\operatorname{tr}Q^2\) as its inputs can return the correct \(\mathcal T_{012}\) for every positive finite frame of the kind described above. The two inputs are identical while the outputs differ.

## 2. A precise third-order observable distinguishes them

The third trace moment is
\[
\boxed{
\operatorname{tr}(Q_s^3)
=\operatorname{tr}(G_s^3)
=3+18r^2+6sr^3.}
\tag{8.5}
\]
For \(r=1/4\):
\[
\boxed{\operatorname{tr}(Q_+^3)=135/32,\qquad
\operatorname{tr}(Q_-^3)=129/32.}
\]
Their difference is \(3/16\). In this very specific three-line model, **the first genuinely phase-sensitive trace information appears at third order**.

For general operator-valued multiword frames, the identity
\[
\operatorname{tr}(Q^3)=
\sum_{i,j,k}\operatorname{tr}(K_{ij}K_{jk}K_{ki})
\]
still holds (with the appropriate weight conventions). But this total real moment ordinarily captures **only the sum of cyclic real contributions**, not each individual complex three-word phase. A labelled complex cyclic trace or stronger pairing certificate would be needed for a full phase dictionary.

The CR-ALP-007 original signed Weil source triple has a negative gauge-invariant labelled product \(\Xi<0\). The abstract arithmetic frame models prove that **both phase signs are compatible with all the inspected two-moment input data**. They do **not** identify the Weil \(\Xi\) with \(\mathcal T\), which live in different mathematical constructions and may have different units/normalizations.

## 3. Source extraction verdict and candidate transfer hypothesis

**Extracted and typed:** \`Kmat\`/\`KD\` from arithmetic period maps, the source Hilbert-space types \`Cu(u)\`, \`Mu(u)\`, \`EE\`, \`Kp\`/\`InputS\`, \`LinkS_at\`, \`BesselId\`, and a concrete example \`D7\` whose word type is \`Unit\`.

**Not extracted:** a published, explicitly evaluable arithmetic **three-word** \`D7\` instance with its \`Pu\` or \`Bu\` matrices, numerical \(\operatorname{tr}(K_{01}K_{12}K_{20})\), lattice/character normalization, or a Weil test-source map. Parameterized \`hostLevelRecord\` is not an evaluated example, and \`InputS\` alone cannot supply phase.

**Noncircular next hypothesis:** A true arithmetic-to-Weil projection would require a specified map \(J\) from the relevant arithmetic carrier to admissible original smooth Weil sources and formulas preserving **(i)** the native metric, **(ii)** the signed numerator including all poles/Gamma/prime-power and character terms, **(iii)** cyclic mixed phases, and **(iv)** quantitative support/scale errors. A frame-moment equality alone is **insufficient**.

### CR-ALP-009 — hard source gate

1. Search actual instantiated multiword arithmetic records, not merely a general \`D7\` record constructor or \`Unit\)-word example. Pin source-defined \`Pu\`/\`Bu\` matrices or report that their evaluation is unsupported.
2. If available, compute a **labelled third-order cyclic trace** \(\operatorname{tr}(K_{uv}K_{vw}K_{wu})\), with the native \(\gamma\) normalizations and phase gauge recorded.
3. Inspect the actual arithmetic \`BesselId\`/lattice input for an independent **phase-coherence** constraint. Do not infer it from \`LinkS_at\` or \(\operatorname{tr}(Q^2)\).
4. Attempt a *source-defined* numerator-and-metric-preserving Weil carrier for that triple. If unavailable, record **NO PHASE TRANSFER FROM TWO-MOMENT DATA**, not **NO POSSIBLE PROJECTION**.

## 4. Reproducibility and custody

The companion standalone Python audit \`src/audits/cr_alp_008_frame_phase_exact_audit.py\` uses **only \`fractions.Fraction\`** for all tested identities, the \(3\times3\) positivity minors, first three trace moments, two opposed cyclic products, literal generic \`InputS\` inequalities, and all eight real sign gauges. **Executed locally: PASS** for \(r=1/4,1/8,1/100,1/10000\). This is a finite algebra regression, not Lean kernel acceptance or a mathematical statement about Alpöge's actual arithmetic family.

No writes to Weil-Lab, RH/provenance, or Hodge-RH ledgers; draft research branch only.
