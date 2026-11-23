# CR-ALP-011 — Ramified isometry corrected; full-diagonal near-isometry does not control mixed complements

**Jurisdiction:** Canonical-Return / independent Alpöge-projection research branch.  
**Date:** 2026-11-23, America/Los_Angeles.  
**Standing:** Lean **source contracts directly inspected**, finite analytic statements **derived**, exact SymPy regression **executed locally: PASS**. **No new Lean compilation, no actual arithmetic three-word period computation, no completed operator equivalence, no RH/whole-Weil positivity theorem.**

## 0. Correction to CR-ALP-010: the restricted attachment **is isometric**

The prior CR-ALP-010 source review correctly identified a potential **complement** gap, but incorrectly left the metric status of the map \(E_u\) into the ramified source subspace uncertain. The frozen Lean source actually contains *additional D7 structure fields*:

- \`src/Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean\`, within \`D7\`, around lines 165, 233, **288–305**: \`eMu : Mu u ≃ₗ[ℂ] Mram\`, \`s_Cu : s.Cu u ≃ₗ[ℂ] Mu0\`, **\`eMu_isom\`** and **\`s_Cu_isom\`**. The latter two assert equality of the relevant native complex inner products with the same form \`(fm u).ipC 0\`; \`sEE\` (near line 307) is a conjugate-linear isometry.
- \`src/Formalism/SharpRestriction/CuspSpanNormBounds/RelevantZeroGramTerms.lean\`, near lines 509–525: \`ramCu D u\` is exactly the inverse image of \`Mram\` under \`s_Cu\`; \`PiU D u\` is its orthogonal projection; \`KD D u v\` is the **whole** period-operator Gram; \`hsDiagD\` controls a \`PiU\`-compressed diagonal.
- \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`, near lines 143–155: in contrast, \`hsDiagFull\` is the **full** diagonal error \(\|KD_{uu}-I\|_{\rm HS}^2\), and \`LinkS_at\` bounds full mixed HS squares under a dimension-comparability antecedent. It does **not** contain a mixed-complement smallness or phase-margin clause.

Define the genuine source map, by the source's own constructions,
\[
E_u:\mathrm{Mu}(u)\longrightarrow s.Cu(u),\qquad
E_u\phi=(s_Cu(u))^{-1}(\mathrm{incl}(eMu(u)\phi)).
\]
The two inner-product compatibility fields give
\[
\boxed{\langle E_u\phi,E_u\psi\rangle_{s.Cu}
=\langle\phi,\psi\rangle_{\mathrm{Mu}(u)},\qquad
\mathrm{ran}(E_u)=\mathrm{ramCu}(u).}
\tag{11.1}
\]
Thus \(E_u\) is an **isometric linear equivalence onto the ramified subspace** for every well-formed \(D7\). Its surjectivity onto the **whole** \(s.Cu(u)\) is separate and equivalent to \(\mathrm{ramCu}(u)=s.Cu(u)\). CR-ALP-010 was amended explicitly with this erratum.

This is a source-contract deduction, not a new Lean kernel theorem established in CR.

## 1. Full period Gram versus restricted Bessel Gram

Write \(P_u=s.Pu(u):C_u\to\mathcal H\), \(R_u=\Pi_u\), \(Q_u=I-\Pi_u\), and
\[
K_{uv}=(c_{0,u}c_{0,v})^{-1/2}P_u^\dagger P_v
\]
with the positive scalar \(c_{0,u}\) from the Lean \`Kmat\` normalization. The bilinear \`besS\` pairing, through its common antiunitary \`sEE\`, encodes \(P_uE_u\), so the **restricted** normalized source Gram is unitarily represented by \(E_u^\dagger K_{uv}E_v\).

Because \(E_u\) is an isometry onto \(R_uC_u\), this restricted Gram preserves the native physical metric **on the ramified subspace**. However, it sees only the \(R_uK_{uv}R_v\) block.

The exact block decomposition is
\[
\boxed{K_{uv}=
R_uK_{uv}R_v+R_uK_{uv}Q_v
+Q_uK_{uv}R_v+Q_uK_{uv}Q_v.}
\tag{11.2}
\]
The four blocks are orthogonal in the Hilbert–Schmidt inner product, and hence
\[
\boxed{
\delta_{uv}^2:=
\|K_{uv}-R_uK_{uv}R_v\|_{\rm HS}^2
=\|K_{uv}\|_{\rm HS}^2
-\|R_uK_{uv}R_v\|_{\rm HS}^2.
}
\tag{11.3}
\]
The **actual new transfer obligation** is a bound on these **mixed-complement** norms \(\delta_{uv}\), relative to a nonzero cyclic phase/sign margin, or a proof that every \(R_u=I\). Neither follows solely from the shape of \`BesselId\` or \`LinkS_at\`.

### 1a. Full near-isometry can actually preserve a large complement

If \(\|K_{uu}-I\|_{\rm HS}\le\epsilon_u\), then by the triangle inequality on \(Q_uC_u\),
\[
\boxed{
\|Q_uK_{uu}Q_u\|_{\rm HS}\ge
\sqrt{\dim(Q_uC_u)}-\epsilon_u.
}
\tag{11.4}
\]
So a very small **full** diagonal defect does **not** imply the complement has small period mass; for positive complement dimension it keeps a nearly isometric block. This is the precise direction of the information loss: the full diagonal near-isometry is compatible with substantial vectors beyond the Bessel-visible subspace.

## 2. Stronger exact countermodel: all \`LinkS_at\`-shaped inequalities and fixed *nonzero* restricted Bessel cycle

This is an **abstract finite Hilbert-space model of the numerical frame hypotheses**, **not an instantiated arithmetic \(D7\)**. In particular the actual \`GamOf\` function, period-map geometry, \`BuRational\`, lattice/divisibility inputs and global family hypotheses are not represented.

Take three source spaces \(C_i=\mathbb C^2\), target \(\mathcal H=\mathbb C^6\), and the **same** rank-one ramified inclusion \(E_i(z)=(z,0)\) for each word.

Set \(t=1/10\), \(a=1/5\), \(b=1/6\), \(c=1/7\). Let the restricted three-vector Gram be
\[
F=\begin{pmatrix}1&t&t\\ t&1&t\\t&t&1\end{pmatrix}\succ0,
\qquad \mathrm{spec}(F)=\{1+2t,1-t,1-t\}.
\]
For \(s\in\{\pm1\}\), let the **orthogonal complementary** vectors have positive Gram
\[
G_s=\begin{pmatrix}1&a&b\\a&1&s\,ic\\b&-s\,ic&1\end{pmatrix}\succ0,
\quad
\det G_s=40211/44100.
\]
Both source systems can be realized in \(\mathbb C^3_F\oplus\mathbb C^3_{G_s}\), hence in a standard \(\mathbb C^6\), by Gram factorization. Each \(P_i^{(s)}:\mathbb C^2\to\mathbb C^6\) has one restricted and one complementary column, chosen orthogonal across those two target sectors. The resulting **type-correct full cross blocks** are
\[
\boxed{
K^{(s)}_{ij}=(P_i^{(s)})^\dagger P_j^{(s)}
=\begin{pmatrix} F_{ij}&0\\0&(G_s)_{ij}\end{pmatrix}.
}
\tag{11.5}
\]

**Exact contract checks.**

1. Every \(K_{ii}^{(s)}=I_2\). Therefore the **entire full diagonal** near-isometry defect is **zero**, not merely small; the compressed one is zero too.
2. The \(E_i\) are exact isometries and the **entire restricted Bessel Gram** is the **same nontrivial** \(F\) for both signs, including its positive cyclic trace \(t^3=1/1000\). Take identical restricted pairing coefficients with \(a_i=1\) and \(U_i=V=I\) to model the pairing-equality part of \`BesselId\`.
3. Each full off-diagonal HS square is \(t^2+|(G_s)_{ij}|^2\), **independent of \(s\)**, with maximum \(1/100+1/25=1/20\).
4. The **literal numerical inequality pattern** of \`LinkS_at\` holds with \(d_i=\mathrm{duS}=2\), \(D_{\rm all}=D_*=6\), \(\Gamma=\sqrt3\), \(C_0=1\), and \(\eta_0=\eta_1=e_0=0\): comparability ratio \(d_i\Gamma^2/D_*=1\); \(1/20<d_i d_j/D_*=2/3\); target dimension \(6=D_*\).
5. The two full frames even have the **same entire spectrum**. Their six-by-six source Gram matrices, after column reordering, are \(F\oplus G_s\), with spectra
\[
\boxed{\{\,1+2t,\ 1-t,\ 1-t,\ 1,\ 1-\sqrt S,\ 1+\sqrt S\,\}},
\quad S=3889/44100,
\]
independent of \(s\). Thus all frame trace moments, determinants and rank data agree.

**Nevertheless** the full labelled three-word cyclic traces differ:
\[
\boxed{\mathcal T_{\rm full}^{(s)}
=\mathrm{tr}(K_{01}^{(s)}K_{12}^{(s)}K_{20}^{(s)})
=t^3+s\,\frac{i}{210}.}
\tag{11.6}
\]
The full diagonal masses are \(2\) per word, so their normalized cycles are
\[
\boxed{\Theta_{\rm full}^{(s)}=
\frac{1}{8000}+s\,\frac{i}{1680}.}
\tag{11.7}
\]
Their oriented imaginary parts have opposite signs although the **restricted Bessel data, full diagonal near-isometries, all off-diagonal HS squares, entire spectrum, source/target dimensions and arithmetic-style comparability numbers are identical**.

The omitted mixed block is \(\mathrm{diag}(0,(G_s)_{ij})\), with
\[
\delta_{ij}^2=|(G_s)_{ij}|^2>0.
\tag{11.8}
\]
This shows exactly why the diagonal near-isometry gives no mixed-complement estimate. A stronger arithmetic input can still rule out such alternatives; we have **not** modelled the full D7 family.

## 3. A special all-ramified source family — and why it is not a nontrivial frame certificate

Source inspection identified a genuine **source-dependent route around the complement**:

- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/NZColC.lean\`: theorem \`atTN_Mram\` states that for \(T\ne2\), the example \`atTN(...)\` has \(\mathrm{Mram}=\top\).
- \`src/Formalism/LocalRingCardinality/SpanPowerMembership/NZRowE.lean\`: \`memC\` uses precisely such an \`atTN\`, and \`hostLevelRecord\` has an explicit hypothesis \`hT2:\forall u,\ T(u)\ne2\`. Hence for **this constructed record**, \(\mathrm{ramCu}(u)=C_u\) and \(R_u=I\), so **there is no complement**.

However, the very same \`NZRowE.lean\` constructor chooses
\[
\boxed{\texttt{sDataN.Pu := fun \_ => 0}.}
\]
Thus its **unmodified** period operators and \(KD\) are zero; its diagonal defect is \(\|0-I_{C_u}\|_{\rm HS}^2=\dim C_u\). For nonzero \(C_u\), it cannot itself demonstrate a small-error \`LinkS_at\` at \(e_0<1\). This is a **toy/source-model limitation**, **not evidence that the later substantive arithmetic chain fails**: other records or later transports can supply different period operators.

This establishes the **exact fork**:
- If \(\mathrm{ramCu}(u)=C_u\) for the **actual nontrivial arithmetic family**, the restricted Bessel Gram determines the full \(KD\) up to the already-pinned positive normalizations and common antiunitary convention. The internal arithmetic normalized cyclic phase preservation then holds.
- Otherwise, one needs an explicit **mixed-complement bound** or source-specific structural relation; \`LinkS_at\` alone does not provide it.

Neither fork supplies the **original signed Weil numerator \(H_{\rm Weil}\)**, its native \(L^2\) metric, or its pole–prime–Gamma support-extension conditions.

## 4. Conditional gate and next

For three words \(i,j,k\) let \(\widehat K_{ij}=R_iK_{ij}R_j\). Put \(\delta_{ij}=\|K_{ij}-\widehat K_{ij}\|_{\rm HS}\). The telescoping/Schatten bound already recorded in CR-ALP-010 gives
\[
\begin{aligned}
|\mathcal T(K)-\mathcal T(\widehat K)|\le{}&
\delta_{ij}\|K_{jk}\|_{\rm HS}\|K_{ki}\|_{\rm HS}\\
&+\|\widehat K_{ij}\|_{\rm HS}\delta_{jk}\|K_{ki}\|_{\rm HS}\\
&+\|\widehat K_{ij}\|_{\rm HS}\|\widehat K_{jk}\|_{\rm HS}\delta_{ki}.
\end{aligned}
\tag{11.9}
\]
Therefore a genuine arithmetic phase-preservation certificate could be obtained if the right side is **strictly smaller than the known oriented imaginary (or relevant real) margin** of the restricted cycle after all positive normalizations. Note that if the restricted cycle is real, even arbitrarily small uncontrolled imaginary corrections can change the **sign of its imaginary part**; a nonzero margin hypothesis is indispensable.

**Certification table**

| Obligation | CR-ALP-011 result |
|---|---|
| Canonical \(E_u:\mathrm{Mu}_u\to\mathrm{ramCu}_u\) metric preservation | **SOURCE CONTRACT PASS** (corrects ALP-010) |
| Exact HS mixed-complement block identity | **ANALYTIC PASS** |
| Full diagonal near-isometry \(\Rightarrow\) negligible complement | **FALSE in general** |
| \`LinkS_at\`-shaped bounds + fixed \`BesselId\` restricted pairing \(\Rightarrow\) full cyclic orientation | **FALSE in exact abstract model** |
| All-ramified \`hostLevelRecord\` example | **SOURCE PASS**, but its specified \(Pu=0\) makes it degenerate for nonzero near-isometric frame tests |
| Exact 6D finite symbolic regression | **LOCAL PASS**, not Lean/CI |
| Actual multiword substantive \(D7\) frame with \(\mathrm{ramCu}\) census | **NOT EXTRACTED** |
| Original signed Weil \((H,M)\) transfer | **OPEN** |

**Next — CR-ALP-012 (actual-source fork):** identify the precise **nontrivial** \(D7\) family used by the published \`LinkS\`/global argument; determine \(\dim\mathrm{ramCu}(u)\), \(\dim C_u\), and whether \(\mathrm{ramCu}=\top\). If proper, seek an explicit mixed-block estimate \(\delta_{ij}\) for genuine \`Pu\`, not generic spectral moments. Independently require a signed numerator-and-metric-preserving map to the original Weil source class. A full arithmetic phase transport alone is not a Weil equivalence.

No changes to Weil-Lab, original RH/provenance, Hodge proof ledgers or formal acceptance status.
