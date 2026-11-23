# CR-ALP-010 — Bessel phase coherence versus period-operator completion

**Canonical-Return / independent Alpöge-projection branch.** Opened 2026-11-23 PST.  
**Standing:** SOURCE-INSPECTED + DERIVED FINITE-ALGEBRA LEMMA + EXECUTED EXACT SYMBOLIC AUDIT PASS. **No new Lean compilation, no numerical multiword theta-lift period instance, no arithmetic–Weil operator intertwiner, no RH or global Weil positivity claim.**

## 0. Question and primary source custody

CR-ALP-008–009 demonstrated that positive-frame rank, Hilbert–Schmidt squared magnitudes, and even every eigenvalue of an unlabelled frame fail to determine the *oriented complex* three-word cyclic trace. Does the actual \`BesselId\` statement supply more, and precisely where would its information cease to determine the **full** arithmetic period-kernel \`KD\`?

Original sources: Alpöge, *A remark* (20 November 2026), https://alpo.ge/remark.pdf; Rubinstein et al., *The Alpöge–Claude proof ... explained* (working draft v2, 22 November 2026), https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, §III.8 and §III.11. Exact Lean source inspected at the pinned **public mirror** \`CocoSgt/riemann-hypothesis-lean@7cbf7a0c32fca1badf05931879a37c3f8a2b03b4\`:

1. \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`, declaration \`BesselId\` (near line 65):
   \[
     D.\mathrm{BuRational}\;\wedge\;
     \exists V:\mathrm{EE}\simeq_{\mathrm{iso},\mathbb C}\mathrm{EE},
     \quad \forall u,\ \exists a_u\in\mathbb C^\times,\
     \exists U_u:\mathrm{Mu}_u\simeq_{\mathrm{iso},\mathbb C}\mathrm{Mu}_u,\
     \mathrm{Bu}_u(\phi,F)=a_u\,\mathrm{besS}_u(U_u\phi,VF).
   \]
   **One target isometry \(V\) is common to every word**, while the nonzero scalar \(a_u\) and source isometry \(U_u\) are word-dependent. This is a real coherence condition, not merely a list of independent norm-square inequalities.
2. \`src/Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean\`: source field \`sEE : what₂.EE →ₗᵢ⋆[ℂ] s.HH\` is **conjugate-linear isometric** (near line 307); \`eMu : ∀ u, what₂.Mu u ≃ₗ[ℂ] (mem u).A.Mram\` (line 165) and \`s_Cu : ∀ u, s.Cu u ≃ₗ[ℂ] (mem u).A.Mu0\` (line 233) are **linear equivalences, not declared isometries**. The same file defines
   \[
   \mathrm{besS}_u(\phi,F)
   =\langle D.sEE(F),\,D.s.Pu_u(E_u\phi)\rangle
   \]
   (near line 353), where \(E_u\) is the canonical composite via the ramified subspace inclusion and the two linear equivalences. It also requires \`Bu_linL\` and \`Bu_linR\` (lines 173–174); \`BuRational\` is an additional arithmetic lattice rationality statement (near line 367).
3. \`src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean\`: the full \`Kmat\`/\`KD\` uses the **entire** \((D.s.Pu_u)^* D.s.Pu_v\), with a positive normalization involving \`c0\`; \`hsSqD\` takes its Hilbert–Schmidt norm square. The map \(E_u:\mathrm{Mu}_u\to D.s.Cu_u\) in \`besS\` need not cover the whole \`Cu_u\`.
4. \`src/Formalism/SharpRestriction/GlobalSchurBounds/TopDensityNGeSixteen.lean\`: \`blockReal_at_of\` unpacks \`BesselId\` and records an embedding \(\varepsilon\) restricted to \`IsIdOrConj\`; the chosen infinite place must be tracked. This is **not** a free license to swap conjugations after fixing a source embedding.
5. \`src/Formalism/LocalToGlobal/SiegelLatticeTheory/FlagStabilizerTensors/KerJoinS.lean\`: \`besselId_of_joins\` derives the pairing-compatibility component via a common target isometry and wordwise maps, under explicit hypotheses. **This research note does not claim an independent Lean proof of the complete \`BesselId\` plus its \`BuRational\` component.**

No concrete evaluated three-word arithmetic \`D7\` operator instance has been recovered in this phase.

## 1. A *positive* normalized cyclic phase-preservation lemma inside arithmetic

In finite-dimensional complex orthonormal coordinates, represent the **bilinear** (linear in each argument) forms \(\mathrm{besS}_u\) and \(\mathrm{Bu}_u\) by matrices \(T_u,A_u:C_u\to E\):
\[
\mathrm{besS}_u(\phi,F)=F^{\mathsf T}T_u\phi,\qquad
\mathrm{Bu}_u(\phi,F)=F^{\mathsf T}A_u\phi.
\]
The bilinear representation is appropriate here: the source's \`sEE\` is **conjugate-linear**, and Lean's inner product in the displayed \`besS\` formula has the compensating conjugation; \`Bu_linL\`/\`Bu_linR\` are separately required.

Under the *pairing-equality component* of \`BesselId\`, with a single complex-linear unitary target \(V\), wordwise complex-linear unitary sources \(U_u\), and nonzero \(a_u\),
\[
\boxed{A_u=a_uV^{\mathsf T}T_uU_u.}
\tag{10.1}
\]
Because \(V^{\mathsf T}\) is unitary when \(V\) is unitary, the type-correct mixed source Gram is
\[
\boxed{A_u^\dagger A_v
=\overline{a_u}a_v\,U_u^\dagger (T_u^\dagger T_v)U_v.}
\tag{10.2}
\]
For three labelled words \(u,v,w\), define
\[
\mathfrak T(A)=
\operatorname{tr}\bigl[(A_u^\dagger A_v)
(A_v^\dagger A_w)(A_w^\dagger A_u)\bigr].
\]
Then unitary cancellations and cyclicity of the trace give
\[
\boxed{\mathfrak T(A)
=|a_ua_va_w|^2\mathfrak T(T).}
\tag{10.3}
\]
The multiplier is **positive real**. Therefore, for \(\mathfrak T(T)\ne0\), its **oriented cyclic phase is preserved** by the entire \`BesselId\` pairing gauge, despite independent wordwise scalars.

More strongly, if each \(\|T_j\|_{\mathrm{HS}}>0\),
\[
\boxed{
\frac{\mathfrak T(A)}
{\|A_u\|_{\mathrm{HS}}^2\|A_v\|_{\mathrm{HS}}^2\|A_w\|_{\mathrm{HS}}^2}
=
\frac{\mathfrak T(T)}
{\|T_u\|_{\mathrm{HS}}^2\|T_v\|_{\mathrm{HS}}^2\|T_w\|_{\mathrm{HS}}^2}.}
\tag{10.4}
\]
This is an **exact finite-algebraic preservation law** on the **same native \(\mathrm{Mu}\) and \(\mathrm{EE}\) spaces**. It corrects a possible overreading of CR-ALP-009: \`BesselId\` does supply a coherent relative phase-preserving dictionary between \`Bu\` and \`besS\`; it does **not** assign the numerical phase, nor identify either with the original signed zeta Weil form.

### 1a. Executable finite regression

The accompanying exact SymPy audit chooses three genuinely complex source columns \(T_u\), a nontrivial unitary \(V\), three different nonzero complex \(a_u\), and three source phases \(U_u\). It checks (10.1)–(10.4) symbolically, including a nonreal cyclic trace. **PASS**. This validates the finite identity, not its Lean theorem status.

## 2. The separate, still-open attachment: restricted Bessel data versus full \`KD\`

The actual Lean \`besS\` sees the **restriction**
\[
\boxed{T_u\ \sim\ (D.s.Pu_u)\,E_u,}
\tag{10.5}
\]
up to the **common target antiunitary identification** induced by \`sEE\`. By contrast, the normalized full arithmetic \`KD\` sees
\[
\boxed{K^{\mathrm{full}}_{uv}\ \sim\
(D.s.Pu_u)^\dagger D.s.Pu_v}
\tag{10.6}
\]
on all of \(D.s.Cu_u,D.s.Cu_v\).

The source composite \(E_u\) lands in the ramified subspace and uses linear equivalences not declared isometric. Even if an idealized \(E_u\) **were an isometric embedding**, knowing \(P_uE_u\) would **not** determine \(P_u^\dagger P_v\) on the orthogonal complement. Neither \`BesselId\` alone nor \`BuRational\` supplies that missing completion information.

### 2a. Exact orientation-reversal *compression* countermodel

This is a **universal finite Hilbert-model falsifier**, not an instantiated arithmetic \(D7\).

Take target \(\mathcal H=\mathbb C^2\), three source spaces \(C_u=\mathbb C^2\), and a **fixed isometric inclusion**
\(E:\mathbb C\hookrightarrow\mathbb C^2,\ E(z)=(z,0)\). Define
\[
c_0=\binom10,\qquad c_1=\binom11,\qquad c_2=\binom1i,
\quad
\boxed{P_u(t)=\begin{bmatrix}c_u&t\overline{c_u}\end{bmatrix}},
\quad t\ge0.
\tag{10.7}
\]
Then \(P_u(t)E=c_u\) for **every \(t\)**. All bilinear restricted Bessel/row pairings can be held **exactly the same** for every \(t\), with \(a_u=1\), \(U_u=V=I\). The *pairing-equality component* of \`BesselId\` is trivially satisfied for the chosen models. Its \`BuRational\` and arithmetic divisibility/lattice hypotheses are **not** being modelled.

But the full period cross operators \(K^{\rm full}_{uv}(t)=P_u(t)^\dagger P_v(t)\) have the labelled cyclic trace
\[
\boxed{\mathfrak T_{\rm full}(t)
=(1+t^2)^3+i(1-t^2)(1+t^2)^2.}
\tag{10.8}
\]
The normalized **full** trace equals
\[
\boxed{
\Theta_{\rm full}(t)=
\frac{\mathfrak T_{\rm full}(t)}
{\prod_{u=0}^2\operatorname{tr}(P_u(t)^\dagger P_u(t))}
=\frac14+\frac{i(1-t^2)}{4(1+t^2)}.}
\tag{10.9}
\]
In particular,
\[
\boxed{\Theta_{\rm full}(0)=\frac14+\frac i4,\qquad
\Theta_{\rm full}(2)=\frac14-\frac{3i}{20}.}
\tag{10.10}
\]
The **oriented full cyclic imaginary phase changes sign**, even though the entire family of **restricted bilinear pairings** is unchanged. All maps are finite and the inclusion \(E\) is exactly isometric. The obstruction is **unseen complementary source directions**, not merely norm loss or missing smoothness.

The companion audit derives the *identity for symbolic real \(t\)*, not just spot values. **PASS**.

**Boundary:** The toy \(P_u(t)\) need not satisfy the actual \(D7\) \`LinkS_at\` near-isometry, arithmetic lattice conditions, \`LinkIRR\`, \`LinkDIV\`, or the native period-map construction. Thus (10.7) is a no-go for transfer from **restricted Bessel pairing data alone**, not a counterexample to all Alpöge hypotheses in combination.

## 3. The new preservation architecture and error gate

There are now **three distinct arrows**, not two:

\[
\boxed{
\mathrm{Bu}\ \overset{\mathrm{BesselId}}{\longleftrightarrow}\
\mathrm{besS}
\quad\Big|\quad
\mathrm{besS}\longrightarrow \mathrm{full\ KD}
\quad\Big|\quad
\mathrm{full\ arithmetic\ data}\longrightarrow(H_{\rm Weil},M_{\rm phys}).
}
\]

- **First arrow: finite algebraic phase-preservation PASS**, assuming the stated pairing-equality component of \`BesselId\`, common target \(V\), and nonzero HS source masses. This does *not* mean we evaluated the actual cyclic trace.
- **Second arrow: NOT JUSTIFIED BY \`BesselId\` ALONE**, due to the incomplete/non-isometric \(E_u\). Explicit compression countermodel shows normalized phase can reverse while the restricted pairings stay constant.
- **Third arrow: OPEN**, because the signed original Weil numerator, pole/Gamma/prime phases, native log-\(L^2\) metric and growing support are not carried by the arithmetic pairing identity.

**Quantitative test for any future proposed completion map.** For comparable cross blocks \(\widehat K_{ij}\) and \(K_{ij}\), with \(\Delta_{ij}=K_{ij}-\widehat K_{ij}\), the typed telescoping identity and Schatten Hölder inequalities give
\[
\begin{aligned}
|\operatorname{tr}(K_{01}K_{12}K_{20})-
\operatorname{tr}(\widehat K_{01}\widehat K_{12}\widehat K_{20})|
\le{}&
\|\Delta_{01}\|_{\rm HS}\|K_{12}\|_{\rm HS}\|K_{20}\|_{\rm HS}\\
&+\|\widehat K_{01}\|_{\rm HS}\|\Delta_{12}\|_{\rm HS}\|K_{20}\|_{\rm HS}\\
&+\|\widehat K_{01}\|_{\rm HS}\|\widehat K_{12}\|_{\rm HS}\|\Delta_{20}\|_{\rm HS}.
\end{aligned}
\tag{10.11}
\]
If this error is smaller than the **absolute imaginary part** of the reference cycle, the sign of that imaginary part is preserved. When comparing to CR-ALP-007's *real signed* Weil triangle, one would instead require a correctly normalized **signed** carrier map and an error smaller than the certified **real** cyclic margin. Merely comparing Hilbert–Schmidt positive kernels does not provide that map.

## 4. Status and next gate

| Obligation | Standing |
|---|---|
| \`BesselId\` common-target/unitary/scalar source contract | SOURCE INSPECTED |
| Bilinear normalized three-cycle under \`BesselId\` | DERIVED / ANALYTIC PASS |
| Exact symbolic finite regression | EXECUTED PASS |
| Actual \`besS\` uses conjugate-linear \`sEE\` | SOURCE INSPECTED |
| Source-to-period map \(E_u\) proven metric-isometric and onto | **NOT ESTABLISHED** |
| Restricted-pairing determination of full \`KD\` cycle | **DISPROVED IN FINITE ABSTRACT MODEL** |
| Concrete multiword arithmetic \(D7\) period matrices | NOT RECOVERED |
| Original Weil \((H,M)\) numerator-and-metric transfer | OPEN |
| New RH theorem / formal Lean certification | NOT CLAIMED |

**CR-ALP-011 objective:** inspect the actual arithmetic **\`ramCu\`/inclusion and complementary projection** together with the \`LinkS_at\` near-isometry hypothesis, to see whether a *quantitative* mixed-complement estimate can close the gap in (10.5)–(10.6). Separately demand a genuine signed-form attachment to Weil's source Gram pencil; passing the Bessel phase-coherence lemma alone is explicitly insufficient.

This work changes no Weil-Lab, provenance, Hodge or RH acceptance state.
