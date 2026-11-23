# CR-DE-001 — Deninger: weight-one Hodge pairing, corrected time reversal, and finite-metric transport tests

**Independent Canonical Return research note · 2026-11-22 PST / 2026-11-23 UTC.** Outcome: **CANONICAL RECOVERY** of Deninger’s conditional rigidity, **ANALYTIC NECESSARY COMPATIBILITY LAW** for anti-equivariant transport, and **FINITE-DIMENSIONAL SIGNATURE COUNTERMODEL**. **No new arithmetic cohomology constructed; no Lean compilation; no unrestricted Weil positivity or RH.** The original RH/PRED/HODGE STOP statuses remain in force.

## 1. Exact primary-source custody

- Christopher **Deninger**, [*There is no “Weil-”cohomology theory with real coefficients for arithmetic curves*, arXiv:2204.02714v2 (27 Jan 2023)](https://arxiv.org/html/2204.02714v2), §2 properties **2.1–2.7**, especially 2.3 (zeros ↔ θ eigenvalues), 2.4–2.6 (trace-weight-one cup pairing), **2.7** (conjectural Hodge star, positivity), **Theorem 2.13** (specific no-real-cohomology theorem), and subsequent discussion of quaternionic target category \(K_{\mathbb R}\).

- Masanori **Morishita**, [arXiv:2508.15971v5](https://arxiv.org/html/2508.15971v5), Theorem **3.6**, anti-equivariant topological flow map; [CR-MO-001](https://chatgpt.com/space/page_b6e89919038c8191bfbdd6d202544491) showed only **orbitwise normalized-Haar unitarity**, not cohomological metric transport.

- Connes–Consani/CCM sources and [CR-CC-001](https://chatgpt.com/space/page_9bd218795c188191b89f97a8f592b814): full original Weil smooth-source window recovery; semilocal Sonin isomorphisms with generally **nonunitary** native metric; the 2026 spectral Schur strategy remains conditional.

- Original RH [PRED-SEAM-16 (2026-09-14)](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_16_DENINGER_EQUALITY_20260914.json) already identified Deninger’s conjectural Hodge metric and absence of a transport theorem. `monocap-tech/hodge-rh/state/RH_HODGE_HEAD.json` checked read-only: **HODGE-R133 STOP**, PAIR-CT/HVAR/FINV/LOC/CPH gates open. This note creates **no** re-entry.

## 2. Deninger’s conditional Hodge mechanism — distinguish statements

Deninger’s proposed complex cohomology has (a) an operator \(\theta\) acting as a **derivation for the cup product** (2.1); (b) \(H^0\simeq\mathbb C\) with \(\theta=0\); (c) \(H^2\simeq\mathbb C(-1)\) with \(\theta=1\) after trace (2.4); and (d) **conjecturally** \(H^1\) with \(\theta\)-eigenvalues equal to all nontrivial Dedekind-zeta zeros, counted with algebraic multiplicity (2.3).

Deninger’s additional **hope** (property **2.7**, not a proved existence theorem) is an **antilinear** automorphism \(J=*\) of \(H^1\) such that \(J^2=-I\), \(J\theta=\theta J\), and \(h(x,y)=\operatorname{tr}(x\smile Jy)\) is **positive-definite Hermitian**. These properties and the derivation rule imply, on an appropriate common invariant domain,

\[\boxed{h(\theta x,y)+h(x,\theta y)=h(x,y).}\]

Equivalently, \(A=\theta-\tfrac12 I\) is **skew-symmetric** with respect to \(h\): \(h(Ax,y)+h(x,Ay)=0\). For any nonzero \(\theta\)-eigenvector \(v\) with eigenvalue \(\rho\), \((2\operatorname{Re}\rho-1)h(v,v)=0\); positive definiteness implies \(\operatorname{Re}\rho=1/2\). Combined with the independently assumed zeros/eigenvalues identification (2.3), this would give RH. This is **Deninger’s existing conditional argument** and was already understood in PRED-SEAM-16.

## 3. Exact finite-dimensional Hodge-metric test

On a finite-dimensional \(\theta\)-invariant subspace with matrix \(\Theta\), write the Hermitian metric as a positive-definite matrix \(M\). The cohomological weight-one identity becomes the **Lyapunov equation**

\[\boxed{\Theta^*M+M\Theta=M,\qquad M>0.}\]

With \(A=\Theta-\tfrac12 I\), \(A^*M+MA=0\). Then \(M^{1/2}AM^{-1/2}\) is skew-Hermitian; in particular the finite spectral generator is diagonalizable with spectrum on \(\operatorname{Re}\rho=1/2\). Merely obtaining an ordinary positive Gram matrix without constructing **the same** \(\Theta\), a common invariant domain, and the Lyapunov equation establishes no Deninger relation.

**Sharp countermodel when positivity is absent:** for any \(\varepsilon\in\mathbb R\setminus\{0\}\), put

\[\Theta=\begin{pmatrix}\frac12+\varepsilon&0\\0&\frac12-\varepsilon\end{pmatrix},\quad M=\begin{pmatrix}0&1\\1&0\end{pmatrix}.\]

Then \(\Theta^*M+M\Theta=M\) **exactly**, but \(M\) has signature \((1,1)\), and \(\Theta\) has off-line eigenvalues. Thus the **derivation/weight-one relation alone is not enough** to force the critical line; positive definiteness, on the actual zero-bearing carrier, is essential. This is an elementary signature obstruction, not an arithmetic counterexample or a counterexample to RH.

## 4. Newly explicit canonical-return compatibility law: time reversal requires an affine half-weight shift

Morishita’s topological map reverses flow. On individual normalized orbit circles, the corresponding generator \(D=-i\partial_t\) is **weight zero** (skew after insertion of \(i\)) and obeys ordinary sign reversal. Deninger’s proposed cohomological generator \(\Theta\), by contrast, satisfies a **weight-one** Hermitian relation. These are differently typed objects.

**No-go theorem (conditional, algebraic):** Suppose nonzero positive Hermitian spaces \((H_D,h_D)\) and \((H_A,h_A)\) have common-domain generators \(\Theta_D,\Theta_A\) each satisfying \(h_i(\Theta_i x,y)+h_i(x,\Theta_i y)=h_i(x,y)\). Suppose \(T:H_D\to H_A\) is an isometric linear injection preserving these domains. Then **a naive, unshifted anti-flow intertwiner**

\[\Theta_A T=-T\Theta_D\]

**cannot exist** on any nonzero invariant subspace. Indeed, from the target relation and the intertwiner, \(h_A(Tx,Tx)=-h_D(x,x)\), whereas isometry gives \(h_A(Tx,Tx)=h_D(x,x)>0\). Hence a contradiction for every \(x\ne0\).

**Correct compatibility law (not an existence claim):** Define the centered/skew generators \(A_i=\Theta_i-\tfrac12 I\). A positive metric-preserving *time reversal* would instead have to satisfy

\[\boxed{A_AT=-TA_D\quad\Longleftrightarrow\quad\Theta_AT=T(I-\Theta_D).}\]

Indeed \(\Theta_AT=T(I-\Theta_D)\) is compatible with the weight-one identity by direct substitution; at the group level it reads \(e^{t\Theta_A}T=e^tTe^{-t\Theta_D}\), i.e. an **\(e^t\) weight twist** rather than raw sign reversal. This is a necessary compatibility condition for *such a cohomological metric-preserving transport*. It is **not** asserted that Morishita, Connes–Consani, or Deninger has constructed \(T\).

**Falsifiable first test for any proposed cross-projection:** a purported isometric lift from Morishita orbit dynamics to two weight-one Deninger Hodge carriers that retains \(\Theta_A T=-T\Theta_D\) without the affine shift is already algebraically impossible. Conversely, merely writing the correct shifted equation does not provide a cohomology, a Hodge star, an eigenvalue identification, or a Weil-form map.

## 5. Quaternionic structure is part of the actual target, not optional decoration

Deninger’s property **2.7** requires an *antilinear* \(J\) with \(J^2=-I\). Therefore \(J\) sends a \(\theta\)-eigenvector of eigenvalue \(\rho\) to one with eigenvalue \(\overline\rho\), and any finite-dimensional **\(J\)-invariant** complex subspace has even complex dimension; in particular the exact real spectral eigenspace \(\rho=1/2\) is quaternionic. A naive complex conjugation \(J^2=+I\), or a real Hilbert model lacking the necessary twisting, does not satisfy the stipulated structure.

Deninger **Theorem 2.13** excludes a **specific functorial real-valued arithmetic Weil-cohomology theory** satisfying the required Artin-multiplicity identification (2.12), using symplectic representations with root number \(-1\). This is **not** a theorem forbidding all real Hilbert spaces or all positive real quadratic forms. A twisted complex/quaternionic target in his \(K_{\mathbb R}\) category remains *possible*, not constructed.

## 6. Effect of the previously recovered Weil-Lab/CCM data

- From [CR-SZ-003](https://chatgpt.com/space/page_f22d2c1da1f0819181937662f715a667), the full original Weil pairing on a selected two-bump subspace has matrix \(\begin{psmallmatrix}q&c\\c&q\end{psmallmatrix}\), with \(c<0\) and **unknown sign of \(q\pm c\)**. A negative off-diagonal entry does **not** establish an indefinite full metric. Without a cohomological \(\Theta\), this subspace cannot be used in the Lyapunov test.

- From [CR-CC-001](https://chatgpt.com/space/page_9bd218795c188191b89f97a8f592b814), finite-place Sonin transport \(\theta_S\) is a bounded Hilbert isomorphism but generally not an isometry of the native metric. To compare a positive target metric to a source one must use the **exact pullback Gram operator** \(M_S=\theta_S^*\theta_S\) (on the correct Hilbert carrier), and independently check \(A^*M_S+M_SA=0\) if claiming Deninger-type centered flow rigidity. These laws do not follow from nonunitarity bounds alone.

- From [CR-MO-001](https://chatgpt.com/space/page_b6e89919038c8191bfbdd6d202544491), matched prime-circle time reversal is compatible with **centered unitary orbit generators**, but does not itself supply an \(H^1\) class or global cup product. Identical orbit lengths and Fourier powers are not a cross-prime sesquilinear form.

## 7. Transfer-admission contract and relation to original RH corpus

A true cross-projection proof toward Deninger would have to construct: **(i)** an independent arithmetic/trace cohomological \(H^1\) and generator \(\theta\) with the correct zeta-zero spectrum (without RH as input); **(ii)** cup product, trace and antilinear quaternionic \(J\) with strictly positive metric; **(iii)** a source-defined map from the actual signed Weil/Schwartz carrier with both poles, Gamma and prime terms that preserves the full form, plus adequate domain and support-quantifier control; **(iv)** if importing Morishita anti-flow, the **affine weight-one twist**, not naïve generator negation; **(v)** if importing CCM Sonin spaces, the precise nonunitary metric change and a valid completed-domain intertwiner.

Even the resulting Deninger realization would generally **bypass** original STRICT/BAL via an independent RH implication, rather than supply the selected Bombieri first-Picard/endpoint theorems. That distinction is explicit in [PRED-SEAM-16](https://github.com/monocap-tech/provenance/blob/main/rh/checkpoints/RH_PRED_SEAM_16_DENINGER_EQUALITY_20260914.json). The checked HODGE-R133 state is STOP and contains different PAIR-CT/HVAR/FINV obligations. Nothing in this page modifies either status.

## 8. Disposition and cross-canonical next step

**RECOVERY:** Deninger’s conditional \(\theta-1/2\) skew/Hodge rigidity and the 2023 no-real-cohomology obstruction.

**EXACT ANALYTIC NECESSITY:** finite Lyapunov test \(\Theta^*M+M\Theta=M\); a naive sign-reversing isometric intertwiner between two positive **weight-one** generators is impossible; correct affine compatibility is \(\Theta_AT=T(I-\Theta_D)\).

**EXACT NEGATIVE TOY:** indefinite weight-one metrics allow symmetric off-critical eigenvalues, so one cannot discard strict positivity or replace it with finite support graphs.

**OPEN:** real arithmetic cohomology candidate in twisted complex/quaternionic category, original Weil-form/Hodge pairing equality, source carrier and completion, uniform positivity, spectral realization, RH.

**Next milestone: CR-INTEGRATION-001** — cross-compare the four author return findings (Suzuki, Morishita, Connes–Consani, Deninger) against exactly five invariants: *source-domain, full signed sesquilinear form, physical metric/weight, arithmetic phase data, and support/flow parameter*. Classify diagrams as PROVED, CONDITIONAL, BLOCKED or UNKNOWN, and identify **one noncircular possible new theorem** to pursue. Do not accidentally promote a mathematical ‘necessary interface’ to a proof of an existing relationship.
