# CR-ALP-007 — Three-source gauge-invariant Weil holonomy and the missing arithmetic phase

**Jurisdiction:** Canonical-Return / Alpöge-projection branch, not Weil-Lab certification.  
**Date:** 2026-11-23 (America/Los_Angeles).  
**Standing:** DERIVED / ANALYTIC with **executed pure-rational interval audit PASS**; external Lean source inspection only. **No new Lean compilation, actual numerical theta-lift two-word instance, global Weil positivity, or independent RH proof.**

## 0. What changes at three cells?

CR-ALP-006 showed a strictly negative *mixed* signed Weil pairing between physically orthogonal ramp sources. The sign of a **single** cross entry is not invariant under independent phase changes of the sources: the two-cell sign can be flipped by replacing one source by its negative. Consequently it would be incorrect to present that sign, by itself, as an intrinsic obstruction to a unitary correspondence.

The first finite diagnostic that retains phase information **modulo independently rephasing each labelled source line** is a three-cell cyclic product. This note constructs one from the **full original signed Weil form**, establishes its negative sign using rational inequalities, and compares its transformation law with the actual arithmetic frame source definitions.

## 1. Source custody: exactly what Alpöge's Lean development records

Primary: Alpöge, *A remark* (20 Nov 2026), https://alpo.ge/remark.pdf. Human-readable working reconstruction: Rubinstein et al., **v2 dated 22 Nov 2026**, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, §III.8 (frame, rank, trace) and §III.11 (arithmetic inputs). Mirror source frozen at **CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4**. Source inspected, **not newly compiled**:

- \`src/Formalism/LocalRingCardinality/DataRows.lean\`: \`Rows.Data\` contains actual source types \`Word\`, \`Mu u\` and common Hilbert target \`EE\`; its \`Bu\` is a word-dependent bilinear pairing \`Mu u → EE → ℂ\`. This is **schematic record data**, not a concrete finite list of numerical period matrices.
- \`src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean\`: \`Kmat s u v = (c0(u)c0(v))^(-1/2) (Pu u)^* (Pu v)\` with native source spaces \`Cu u\`, \`Cu v\`, and common target; \`KD\` is the instantiated version and \`hsSqD\` takes its Hilbert–Schmidt norm square. The source **retains a complex mixed operator**, although later estimates retain only the norm square.
- \`src/Formalism/SizeWindowBrackets/TraceSquareAndDiagBounds.lean\`: for finite matrices \`Kp B ga i j = sqrt(ga_i ga_j) B_i^* B_j\`, \`Qm=sum ga_i B_i B_i^*\`, and \`InputS\` bounds the off-diagonal \`hsSq(Kp)\` and diagonal defects. The order \`B_i^*B_j\` is type-correct.
- \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`: \`LinkS_at\` has a **dimension-comparability antecedent**; \`BesselId\` only asserts existence of nonzero complex scalars and isometries matching \`Bu\` to an arithmetic Bessel pairing. It does not fix a canonical complex phase for each word.
- \`src/Formalism/SizeWindowBrackets/JointDPublishedSizes.lean\`: definitions \`KpGTransport\` and \`InputSTransport\` explicitly *state* frame transport under word-dependent nonzero complex rescaling, with \`ga_i\` compensated by \`|c_i|^{-2}\`. A \`def ... : Prop\` is a **statement**; no proof of these precise declarations was independently identified or compiled for CR.

**Extraction status:** The actual operator **schema**, its normalized cross block, and the native source/target Hilbert types are now pinned. A specific instantiated \`D7\` with two or three concrete numeric \`Word\` operators, their matrices, and an independent Weil-source map **has not been extracted**. CR will not substitute the ramp prototypes for those arithmetic operators.

## 2. Frame gauges: what the second moment cannot see

For any finite family of matrices \(B_i:C_i\to\mathcal H\) with weights \(\gamma_i>0\), define
\[
K_{ij}=\sqrt{\gamma_i\gamma_j}\,B_i^*B_j,\qquad
Q=\sum_i\gamma_i B_iB_i^*.
\]
Choose a common unitary \(V\) on \(\mathcal H\), and nonzero wordwise scalars \(c_i\in\mathbb C\). Put
\[
B'_i=c_i VB_i,\quad \gamma'_i=\gamma_i/|c_i|^2,\quad\omega_i=c_i/|c_i|.
\]
Direct matrix calculation gives
\[
\boxed{K'_{ij}=\overline{\omega_i}\omega_j K_{ij}},\quad
Q'=VQV^*,\quad \|K'_{ij}\|_{\rm HS}^2=\|K_{ij}\|_{\rm HS}^2.
\]
Thus **all** trace/rank and off-diagonal HS-square data can be unchanged while the phase of a *particular labelled mixed entry* changes. This is a covariance of the **source schema**, not an assertion that every rephasing corresponds to a distinct permitted arithmetic lattice realization.

For three distinct labelled cells \(i,j,k\), the cyclic operator product is well typed:
\[
\boxed{
\mathcal T_{ijk}
 =\operatorname{tr}_{C_i}(K_{ij}K_{jk}K_{ki})
\in\mathbb C.
}
\tag{7.1}
\]
Its gauge factors cancel, so \(\mathcal T'_{ijk}=\mathcal T_{ijk}\). Moreover,
\[
|\mathcal T_{ijk}|
\le \|K_{ij}\|_{\rm HS}\|K_{jk}\|_{\rm HS}\|K_{ki}\|_{\rm HS}.
\tag{7.2}
\]
For a \`LinkS\`-style estimate \(\|K_{ij}\|_{\rm HS}^2\le(1+\eta_0)d_i d_j/D_*\), the magnitude obeys
\[
|\mathcal T_{ijk}|\le(1+\eta_0)^{3/2}d_i d_j d_k/D_*^{3/2},
\]
provided the size-comparability antecedent and \(D_*>0,\eta_0\ge-1\) hold. **No phase or real-part sign follows from the bound.** Actual arithmetic phase information would require computing \(\mathcal T\), not only the \(HS^2\) values.

## 3. Three disjoint original physical Weil sources

Let \(L=\log9\), \(f(u)=e^{u/2}\mathbf1_{[0,L]}(u)\), and define \(f_j(u)=f(u-jL)\), \(j=0,1,2\). These are nonsmooth but compact \(L^2\) sources, with pairwise disjoint interiors:
\[
\boxed{M_{ij}=\langle f_i,f_j\rangle=8\,\delta_{ij}.}
\tag{7.3}
\]
Each **individual** source has the same autocorrelation \(k_9\), the same complete diagonal signed Weil value \(d\), the same pole/prime/Gamma diagonal, and the identical ramp samples \((9/n-1)_+\).

Let \(H_{ij}\) be the original signed Weil Hermitian pairing, defined for these real sources by the full even symmetrized autocorrelation
\[
k_{ij}(t)=\tfrac12[k_9(t-(j-i)L)+k_9(t+(j-i)L)].
\]
CR-ALP-002's complete two-pole, negative prime-power, full-Gamma normalization is used—no surrogate PSD kernel. The matrix is real symmetric Toeplitz:
\[
\boxed{
H=\begin{pmatrix}
d&c_1&c_2\\
c_1&d&c_1\\
c_2&c_1&d
\end{pmatrix}.
}
\tag{7.4}
\]
The previous phase derived, and the new independent audit rechecks,
\[
\boxed{-0.132426<c_1<-0.132425,\qquad 0.450<d<0.456.}
\tag{7.5}
\]

### 3a. New exact \(2L\)-separation entry, including *all* prime powers

For \(i=0,j=2\), \(k_{02}(t)=\tfrac12 k_9(|t-2L|)\) for \(t\ge0\); its support is exactly \([L,3L]\) and \(k_{02}(0)=0\).

**Two poles**:
\[
\mathcal P_{02}=8L(9+1/9)=\frac{656}{9}L.
\]

**Prime powers**: for \(n=p^h\), the coefficient of \(-\Lambda(n)\) is
\[
v_2(n)=\frac{2k_{02}(\log n)}{\sqrt n}
=
\begin{cases}
0,&n\le9\text{ or }n>729,\\
1-9/n,&9<n\le81,\\
81/n-1/9,&81<n\le729.
\end{cases}
\]
Writing \(c_p^{(2)}=\sum_{h\ge1}v_2(p^h)\), this includes **129 primes with nonzero coefficients**, and gives the *entire* prime part
\(-\sum_{p\le727}c_p^{(2)}\log p\).
Exact spot checks: \(c_2^{(2)}=4285/1536\), \(c_3^{(2)}=16/9\).

**Full Gamma**: the diagonal constant is zero since \(k_{02}(0)=0\), but the integral is not. For \(L<t<2L\) and \(2L<t<3L\), respectively, the integrand \(-2 k_{02}(t)e^{-t/2}/(1-e^{-2t})\) is
\[
-\frac{1-9e^{-t}}{1-e^{-2t}},
\qquad
-\frac{81e^{-t}-1/9}{1-e^{-2t}}.
\]
Elementary antiderivatives give
\[
\boxed{
\Gamma_{02}=
-\frac89L-\frac12\log\frac{82}{81}
+\frac1{18}\log\frac{6643}{6642}
+\frac92\log\frac{50}{41}
+\frac{81}{2}\log\frac{3650}{3731}.
}
\]
Therefore the **exact full signed mixed entry** is
\[
\boxed{
\begin{aligned}
c_2={}&72\log9-\sum_{p\le727}c_p^{(2)}\log p\\
&-\frac12\log\frac{82}{81}
+\frac1{18}\log\frac{6643}{6642}
+\frac92\log\frac{50}{41}
+\frac{81}{2}\log\frac{3650}{3731}.
\end{aligned}}
\tag{7.6}
\]
A standard-library **pure rational interval audit** (atanh series with an explicit geometric tail for every logarithm) proves
\[
\boxed{-0.001978<c_2<-0.001977<0.}
\tag{7.7}
\]
This negative sign is **much smaller** than \(c_1\) and cannot be inferred from nearby-case continuity alone. The full prime-power calculation is essential.

### 3b. A gauge-invariant, nontrivial signed cycle

Because \(H_{01}=H_{12}=c_1<0\) and \(H_{20}=c_2<0\),
\[
\boxed{H_{01}H_{12}H_{20}=c_1^2c_2<0.}
\tag{7.8}
\]
Under independent source phases \(f_j\mapsto e^{i\theta_j}f_j\), the entries transform as \(H_{ij}\mapsto e^{-i\theta_i}e^{i\theta_j}H_{ij}\); the three factors cancel in the product. Its phase \(\pi\) is the **labelled three-line Bargmann-type cyclic invariant**. The normalized version
\[
\boxed{\Xi=\frac{H_{01}H_{12}H_{20}}{H_{00}H_{11}H_{22}}<0}
\]
is also invariant under positive independent source scalings. The coarse certified bound is
\[
-0.00039<\Xi<-0.00036.
\]
This is invariant under **diagonal gauge transformations of the labelled source lines**, *not* under arbitrary changes of the three-dimensional basis. It is not a new global operator invariant independent of the specified source geometry.

This strengthens CR-ALP-006: an individual negative \(H_{01}\) can be rephased away, but one **cannot make all three** mixed entries simultaneously positive real by individual source rephasings.

### 3c. Positive three-dimensional restriction despite nontrivial holonomy

The exact rational audit also confirms \(d>0.450\) and \(|c_1|<0.132426\), \(|c_2|<0.001978\). Gershgorin's strict diagonal-dominance lower bound yields
\[
\lambda_{\min}(H)\ge d-\max(2|c_1|,|c_1|+|c_2|)
>0.450-0.264852=0.185148.
\]
Since \(M=8I_3\),
\[
\boxed{\lambda_{\min}(H,M)>0.02314>0.}
\tag{7.9}
\]
The **original** signed Weil restriction is positive on this selected three-dimensional span. The negative cyclic invariant measures its labelled off-diagonal geometry, **not** a negative eigenvalue or violation of Weil positivity.

### 3d. Smooth original Weil tests

Choose any nonnegative even \(C_c^\infty\) mollifier \(\rho_\eta\), and translate its three smoothed ramp sources by offsets \(j(L+2\eta)\), \(j=0,1,2\). Their supports are still disjoint, so their **physical mixed metric remains zero exactly**, and each source has the same individual autocorrelation. As \(\eta\downarrow0\), all full signed mixed pole/prime/Gamma functionals converge to (7.4) by the same fixed-support Lipschitz/compactness argument as CR-ALP-003/006. For all sufficiently small \(\eta\), the **three mixed signs, negative cyclic product and positive three-dimensional restriction persist**. This is a strict analytic **existence** result; no explicit \(\eta\)-threshold or completed semilocal/Suzuki domain has been certified.

## 4. A more discriminating compatibility conjecture

The correct next test is *not* to equate a signed Weil cross entry with an arithmetic **Hilbert–Schmidt norm square**, which is real and nonnegative, nor to identify \(H_{ij}\) simply with the positive physical Gram \(\langle f_i,f_j\rangle\).

Instead, freeze a **concrete** arithmetic triple of word operators \(B_i\), their native Hilbert structures and weights \(\gamma_i\), and compute both:
\[
\mathcal T_{012}^{\rm arith}=
\mathrm{tr}(K_{01}K_{12}K_{20})
\quad\text{and}\quad
\Xi_{012}^{\rm Weil}=
\frac{H_{01}H_{12}H_{20}}{H_{00}H_{11}H_{22}}.
\]
A proposed transfer must identify an *explicit, normalization-preserving rule* relating them, including any character, pole/Gamma, sign and scale corrections. Equality of these two numbers is **not** presently conjectured as a theorem: they may have different dimensions and the first vanishes for orthogonal positive frame blocks even where the second is negative. The first objective is a falsifiable **gauge-covariant dictionary**, not assumed equality.

**Key asymmetry:** Alpöge's published \(LinkS\) input controls only \(\|K_{ij}\|_{\rm HS}^2\) and diagonal defects; this cannot determine the phase of the arithmetic cyclic trace. The complete complex mixed operators, or a separate phase-coherent pairing input from \(BesselId\), would be necessary. Source \`BesselId\` permits word-dependent nonzero complex scalars, so a phase-normalization condition has not been automatically supplied.

## 5. Gate result / next

| Obligation | Result |
|---|---|
| Actual \`Rows.Data.Bu\`, \`KD\`, \`Kp\` **type and normalization** source anchors | SOURCE INSPECTED |
| Concrete numeric arithmetic 3-word triple and full \`Kp\` operators | **NOT EXTRACTED** |
| Original three-cell physical metric \(M=8I_3\) | EXACT |
| Full original three-cell signed \(H\), including *129-prime* distant interaction | EXACT ANALYTIC / RATIONAL AUDIT PASS |
| Negative **gauge-invariant** three-line cyclic product | EXACT / RATIONAL AUDIT PASS |
| Positive generalized three-dimensional restriction | EXACT BOUND / RATIONAL AUDIT PASS |
| Smooth-source persistence | ANALYTIC EXISTENCE, not explicit threshold |
| Theta-lift–Weil operator intertwiner / RH theorem | OPEN / NOT CLAIMED |

**Next — CR-ALP-008:** source-assemble a *genuinely instantiated* arithmetic \`D7\` / frame triple (not a schematic record), certify either a computable complex cyclic trace or a source-level obstruction to extracting it, and state a genuine transfer that preserves both \(M\) and \(H\) including phase. If the arithmetic sources only provide \(LinkS\) norm-square bounds, report **NO PHASE TRANSFER** under those hypotheses.

**Custody:** research-only CR branch; no promotion to Weil-Lab, Hodge, or any formal RH claim.
