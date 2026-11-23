# CR-ALP-009 — All spectral moments can miss an oriented arithmetic cyclic phase

**Canonical-Return / Alpöge projection — 2026-11-23 PST**  
**Standing:** source-inspected Lean comparator; finite exact analytic lemma; local Gaussian-rational audit **PASS**. **NO newly compiled Lean theorem; NO computed three-word theta-lift period triple; NO full Weil/operator equivalence or new RH result.**

## 1. Research question, source pins, and the extraction gate

CR-ALP-008 showed that the *first two* trace moments of a positive arithmetic-style frame cannot determine a real signed three-cycle. CR-ALP-009 asks whether computing **all spectral moments** would fix the complex *oriented* three-cycle, and whether the actual theta-lift source provides a stronger, explicitly evaluable phase.

The relevant external mathematical sources remain Alpöge, *A remark* (20 Nov 2026), https://alpo.ge/remark.pdf, and Rubinstein et al., *The Alpöge–Claude proof ... explained*, v2 (22 Nov 2026), https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, §III.8 and §III.11. This note is primarily an independent finite-algebra investigation informed by the **actual Lean text** at the immutable public mirror commit
\`CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4\`:

- \`src/Formalism/SizeWindowBrackets/FrameTraceRatioFinEquiv.lean\`: \`KpG\`, \`QmG\`, \`InputSG\`. Type-correct \(K_{ij}=\sqrt{\gamma_i\gamma_j}B_i^*B_j\); \`InputSG\` bounds **squared Hilbert–Schmidt magnitudes**, diagonal defects, and target dimensions.
- \`src/Formalism/SizeWindowBrackets/TraceSquareAndDiagBounds.lean\`: analogous \`Kp\`, \`Qm\`, \`InputS\`.
- \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`: \`LinkS_at\` **conditional on dimension comparability**; \`BesselId D\` asserts \`D.BuRational\` and a **common target isometry** \(V\), plus **word-dependent nonzero complex scalars** \(a_u\) and **word-dependent source isometries** \(U_u\), such that \(Bu_u(\phi,F)=a_u\,besS_u(U_u\phi,VF)\). This is additional source structure, **not** numerical evaluation of a three-word cyclic trace.
- \`src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean\`: \`Kmat\`/\`KD\` are actual normalized **complex** cross-operator schemas from period maps \(P_u^*P_v\).
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/ExampleClassGroupCardPos.lean\`: \`exampleRecordOfData\` has \`wordFin := Fintype Unit\`.
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/RankRecordRowDimensionZero.lean\`: \`rankRecord\` has \`wordFin := Fintype Unit\`.
- \`src/Development/AbelianVarietyTheory/Definitions/SemilocalFrames/GlobalRowsOfSizeParameters.lean\`, \`src/Development/AbelianVarietyTheory/Definitions/LatticeFrames/HypOver.lean\`: inspected concrete constructors \`recd\` also pass \`(W := Unit)\`.
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/NZRowE.lean\`: \`hostLevelRecord\` has a **parameterized** finite word type \(W\), with a caller-supplied \`per\` pairing. This is **not an evaluated concrete arithmetic triple**.

**Source-extraction result:** No specific **instantiated, numerically evaluable three-word** \(D7\) \`Pu\`/\`Bu\` operator triple was established by these searches. The concrete inspected example families are single-word, and the general constructor is parametric. This is **not proof that no such triple exists elsewhere** in the very large Lean repository. The prior plan to compute an actual arithmetic cyclic trace remains unsatisfied.

## 2. Stronger theorem: the complete frame spectrum is phase-orientation blind

Let
\[
a=\tfrac15,\quad b=\tfrac16,\quad c=\tfrac17,\quad S=a^2+b^2+c^2=\tfrac{3889}{44100}.
\]
For \(s=\pm1\), take two labelled Hermitian matrices
\[
\boxed{
G_s=\begin{pmatrix}
1&a&b\\
a&1&s\,ic\\
b&-s\,ic&1
\end{pmatrix}.}
\tag{9.1}
\]
They are strictly positive: their upper-left \(2\times2\) minor is \(1-a^2>0\), and
\[
\boxed{\det G_s=1-S=\tfrac{40211}{44100}>0.}
\tag{9.2}
\]
Hence, by Gram factorization, each arises from three actual unit vectors \(B_i^{(s)}:\mathbb C\to\mathbb C^3\), \((B_i^{(s)})^*B_j^{(s)}=(G_s)_{ij}\). With frame weights \(\gamma_i=1\), the corresponding \(Q_s=\sum_i B_i^{(s)}(B_i^{(s)})^*\) have the **same nonzero spectra as \(G_s\)**.

The characteristic polynomial is **identical** for both signs:
\[
\boxed{\det(\lambda I-G_s)=(\lambda-1)^3-S(\lambda-1),}
\tag{9.3}
\]
so both have the **same complete spectrum**
\[
\boxed{\operatorname{spec}(Q_s)=\{1-\sqrt S,\;1,\;1+\sqrt S\}.}
\tag{9.4}
\]
In particular, **every spectral trace moment, determinant, rank, condition number, spectral measure and eigenvalue-based frame statistic agrees**. Examples:
\[
\operatorname{tr}(Q_s)=3,\quad
\operatorname{tr}(Q_s^2)=3+2S,\quad
\operatorname{tr}(Q_s^3)=3+6S.
\tag{9.5}
\]
All pairwise absolute overlaps are also identical: \(\{|G_{01}|,|G_{02}|,|G_{12}|\}=\{a,b,c\}\), and the literal generic finite-frame \`InputS\`/\`InputSG\` inequalities hold with \(D=D_*=3\), \(\eta_0=\eta_1=e_0=0\), because every \(|G_{ij}|^2\le 1/3\) and every diagonal defect vanishes.

**But the oriented three-word cyclic phase is opposite:**
\[
\boxed{
\mathcal T_{012}(G_s)
=(G_s)_{01}(G_s)_{12}(G_s)_{20}
=s\,i\,abc=s\,\frac{i}{210}.}
\tag{9.6}
\]
Thus \(\mathcal T_{012}(G_+)=+i/210\) whereas \(\mathcal T_{012}(G_-)=-i/210\). This cyclic product is unchanged under \(G_{ij}\mapsto\bar\omega_i\omega_jG_{ij}\) for independent unit complex phases \(\omega_i\). Unlike the weaker \(r,r,r\) examples, the three **different** magnitudes \(a,b,c\) prevent even a nontrivial permutation of the labels from preserving the entire edge-magnitude pattern.

**Exact mechanism:** \(G_-=\overline{G_+}\) entrywise. Complex conjugation is an **antiunitary**, rather than complex-linear, change of the underlying labelled frame data. It preserves all real spectral information while reversing the imaginary oriented cycle. The two matrices are unitarily similar as **unlabelled Hermitian operators** because they have equal spectra, but **not equivalent as labelled Gram data under independent line phases and edge-magnitude-preserving relabellings**.

### 2a. Precise NO-GO and scope

> No universally valid rule that uses only the complete spectrum of \(Q\), any/all of its trace moments, the labelled pairwise HS overlap **magnitudes**, diagonal source masses, and ranks can recover the **oriented imaginary part** of a labelled cyclic trace \(\mathcal T_{ijk}=\operatorname{tr}(K_{ij}K_{jk}K_{ki})\) for every finite positive frame.

The two positive examples agree on every listed input and disagree on the requested output. The conclusion is **finite and unconditional**, but deliberately **restricted to those inputs**.

It does **not** refute:
- the existence of a richer arithmetic phase recovered from the actual \`Pu\`/\`Bu\` operators, \`BesselId\`, Galois or complex-embedding data;
- a **real-valued** quantitative arithmetic-to-Weil inequality without a full complex operator map;
- the previously certified **negative real** cyclic product of CR-ALP-007. In a three-line model, the **real part** of a cyclic trace can contribute to \(\operatorname{tr}Q^3\), so the present imaginary-orientation no-go is **strictly different** from the CR-ALP-008 second-moment no-go;
- any theorem about the actual Alpöge theta-lift family. \`InputS\` is generic finite frame geometry, **not** \`LinkS_at\` plus \`LinkIRR\` plus \`BesselId\` plus divisibility/arithmetical hypotheses.

In particular, an antiunitary conjugation might be prohibited or fixed by the arithmetic embedding and pairing conventions. That must be **checked on an actual source instance**, not guessed.

## 3. Consequences for the CR-ALP transfer design

This phase settles a hard negative possibility: escalating from second moments to **all unlabelled spectral moments** still cannot recover a **full complex orientation-sensitive correspondence**. The proposed operator-intertwiner must therefore include **labelled mixed cross operators and an orientation/phase convention**, in addition to the complete original signed Weil numerator and native physical metric.

The practical research program should now separate:

1. **Analytic one-sided projection (established, conditional):** CR-ALP-004 transfers the arithmetic upper bound on \(G(x)\) into a sharp lower bound for a specified family of smooth original signed Weil tests, including all prime, pole and Gamma terms. No theta-lift isometry is needed for this statement.
2. **Full structural/operator projection (open):** requires real source-defined maps matching complete native metric \(M\), signed numerator \(H\), labelled cross interactions and normalized cyclic phases, with nonzero Mellin pole moments and support growth retained. Generic \`InputS\`, even enlarged to the full spectrum, is insufficient.

**Next CR-ALP-010:** pin the exact *additional arithmetic source input* that could control oriented cyclic phases: audit \`BesselId\` witnesses together with \`Pu\`/\`Bu\` and Galois/complex-embedding normalizations, seeking a genuine, explicit multiword family. If the source only proves moment/size inequalities and provides no evaluated cyclic traces, classify **FULL OPERATOR PROJECTION: NOT ESTABLISHED; ANALYTIC TEST-FAMILY PROJECTION: RETAINED**. Do not construct dummy numeric operators and call them arithmetic lifts.

## 4. Executable audit / governance

\`src/audits/cr_alp_009_conjugate_phase_exact_audit.py\` uses only Python's \`Fraction\` and exact Gaussian-rational pairs. It checks Hermitian structure and leading minors, exact determinant and characteristic polynomial, first three and 14 higher trace moments, the literal finite \`InputS\` inequalities, cyclic traces \(\pm i/210\), all eight real sign gauges, and the fact that the unequal labelled edges forbid nontrivial magnitude-preserving permutations. **Executed locally: PASS**. This is an algebraic regression, not Lean acceptance or an evaluation of the actual arithmetic theta lifts.

No writes to Weil-Lab or accepted RH/provenance/Hodge records; CR draft investigation only.
