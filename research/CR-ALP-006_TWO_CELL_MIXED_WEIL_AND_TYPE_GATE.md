# CR-ALP-006 — Two-cell mixed Weil coupling and the typed arithmetic frame boundary

**Independent Canonical-Return / Alpöge projection** · opened 2026-11-23 PST.  
**Disposition:** **DERIVED / ANALYTIC** finite two-source lemma; **AUDIT-PASS** exact rational interval script executed locally; **SOURCE-INSPECTED** actual Lean operator schema; **NO ARITHMETIC TWO-WORD TRANSFER CLAIM**. No Lean compilation, completed-operator equivalence, or new RH claim.

## 0. Exact primary-source scope

- Levent Alpöge, *A remark*, 20 Nov 2026, https://alpo.ge/remark.pdf.
- Rubinstein et al., *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, **v2, 22 Nov 2026**, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, **§III.8 printed pp. 138–142**, Statement (S), Lemma T1, Corollary T2, especially (III.219)–(III.228). The manuscript quotes its arithmetic operator inputs from Lean and is not itself a proof of all of them.
- Frozen public Lean mirror **7cbf7a0c32fca1badf05931879a37c3f8a2b03b4** in \`CocoSgt/riemann-hypothesis-lean\`, source checked:
  1. \`src/Formalism/SizeWindowBrackets/TraceSquareAndDiagBounds.lean\`: typed \`Kp\`, \`Kd\`, \`Qm\`, \`TraceIdentities\`, \`InputS\`.
  2. \`src/Formalism/SharpRestriction/GlobalTheoremInputs.lean\`: \`LinkS_at\`, \`LinkIRR\`, \`BesselId\`.
  3. \`src/Formalism/LocalRingCardinality/DataRows.lean\`: source \`Rows.Data.Bu : ∀ u, Mu u → EE → ℂ\`, with separate native complex Hilbert structures and actual word-dependent sources.
  4. \`src/Formalism/LocalRingCardinality/DataParamsS.lean\`: \`SizeData₂.hsSq\` delegates to source \`S.hsSqD\`.

**Important notation hygiene:** The draft's Lemma T1 prints \(K_{ww'}=B_wB_{w'}^\dagger\) although it declares \(B_w:C_w\to\mathcal H\) with word-dependent \(C_w\); that literal product would not generally compose. The *type-correct cross block* is \(K_{ww'}=B_w^\dagger B_{w'}:C_{w'}\to C_w\). This is precisely the orientation in the Lean definition
\[
\texttt{Kp B ga u u'}=\sqrt{\gamma_u\gamma_{u'}}\, (B_u)^\dagger B_{u'},
\]
with dimensions \(\operatorname{Matrix}(\operatorname{Fin}(d_u))(\operatorname{Fin}(d_{u'}))\).
The source defines \(Q=\sum_u\gamma_u B_u B_u^\dagger\) on the common target. The trace-square identity using cross blocks is valid with this orientation. **This is an exposition notation/type correction, not an allegation that the checked Lean proof contains the misprint.** The arithmetic \`LinkS_at\` input is moreover an implication whose dimension-comparability antecedent must be established.

## 1. A concrete two-source *original signed Weil* test

Let \(L=\log9\), and use the *actual Alpöge-ramp source* from CR-ALP-002 with cutoff \(a=9\), weight \(w=1\):
\[
f(u)=e^{u/2}\mathbf1_{[0,L]}(u),\qquad
g(u)=f(u-L).
\]
The individual autocorrelations are **identical** and reproduce every sampled prime table value \((9/n-1)_+\). Their native physical \(L^2\) Gram matrix is
\[
\boxed{M=\begin{pmatrix}8&0\\0&8\end{pmatrix}.}
\tag{6.1}
\]
The supports meet only at one endpoint. Thus the sources are exactly physically orthogonal (as \(L^2\) functions).

To evaluate the **genuine mixed** signed Weil form, use polarization:
\[
k(t)=(9e^{-|t|/2}-e^{|t|/2})_+,\quad
k_{\mathrm{mix}}(t)=\tfrac12[k(t-L)+k(t+L)].
\]
The sesquilinear cross value is \(H_{12}=\mathcal W(k_{\mathrm{mix}})\) with the **complete two poles, negative von Mangoldt prime sum and archimedean Gamma** from CR-ALP-002. This is not the form of a different positive surrogate. The symmetrized cross test is even, \(k_{\mathrm{mix}}(0)=0\), and
\[
k_{\mathrm{mix}}(t)=\tfrac12 k(|L-t|),\quad 0\le t\le 2L.
\]

### 1a. All contributions are explicit

The two poles contribute
\[
\boxed{\mathcal P_{12}=\frac{80}{3}\log9.}
\]

For \(n=p^j\le81\), the weighted even cross-prime sample satisfies
\[
\frac{2k_{\mathrm{mix}}(\log n)}{\sqrt n}=
\begin{cases}
3(n-1)/n,&2\le n\le9,\\
(81-n)/(3n),&9<n\le81,\\
0,&n>81.
\end{cases}
\]
Therefore the **full prime contribution** is the finite, exact expression
\[
\boxed{\mathcal P^{\rm prime}_{12}=-\sum_{p\le79}c_p\log p,}
\quad c_p=
\sum_{\substack{j\ge1\\p^j\le81}}
\begin{cases}
3(p^j-1)/p^j,&p^j\le9,\\
(81-p^j)/(3p^j),&9<p^j\le81.
\end{cases}
\tag{6.2}
\]
This includes every prime power, not only primes or the \(n=9\) interaction. Example exact coefficients: \(c_2=533/64,\ c_3=16/3\).

Because \(k_{\mathrm{mix}}(0)=0\), the Gamma diagonal constant proportional to Euler's constant vanishes, but the **Gamma integral remains nonzero**:
\[
\begin{aligned}
\Gamma_{12}
&=-\int_0^{2L} k(|L-t|)\frac{e^{-t/2}}{1-e^{-2t}}\,dt\\
&=\boxed{-3\log5-\frac{27}{2}\log\frac{50}{41}
+\frac13\log9+\frac16\log\frac{82}{81}.}
\end{aligned}
\tag{6.3}
\]
For \(0<t<L\), the integrand simplifies to \(3/(1+e^{-t})\); for \(L<t<2L\) it reduces to \((27e^{-t}-1/3)/(1-e^{-2t})\). This is an **elementary exact closed form**, not a numerical approximation to the Gamma term.

Combining all three pieces:
\[
\boxed{
H_{12}=27\log9-\sum_{p\le79}c_p\log p
-3\log5-\frac{27}{2}\log\frac{50}{41}
+\frac16\log\frac{82}{81}.
}
\tag{6.4}
\]

### 1b. Rigorous sign and positive finite restriction

A pure rational interval auditor, using the positive atanh series and explicit geometric remainder for every logarithm, proves
\[
\boxed{-0.132426<H_{12}<-0.132425<0.}
\tag{6.5}
\]
This is a **rigorous negative original mixed Weil pairing between physically orthogonal** nonnegative \(L^2\) ramp sources, including *all* prime powers, poles, and Gamma.

For either individual diagonal, the full original expression is
\[
H_{11}=H_{22}=16\log9-\frac{39}{4}\log2-4\log3-\frac85\log5-\frac47\log7
-8\gamma-8\log\pi+10\log5-16\log4.
\tag{6.6}
\]
Using the elementary Euler-constant bounds at \(N=100\) and the classical rational bounds \(223/71<\pi<22/7\), the same **exact rational** auditor certifies
\[
\boxed{0.450<H_{11}=H_{22}<0.456.}
\tag{6.7}
\]
Therefore the *actual two-dimensional* Hermitian pencil has
\[
H=\begin{pmatrix}d&c\\c&d\end{pmatrix},\quad
M=8I_2,\quad d>0.450,\quad -0.132426<c<-0.132425,
\]
and both generalized eigenvalues are positive. In particular,
\[
\boxed{\lambda_{\min}(H,M)=\frac{d+c}{8}>0.0396.}
\tag{6.8}
\]
The **negative mixed interaction does not negate positivity on this selected two-dimensional span**. This is a local pencil, not any conclusion about the complement or whole Weil operator.

### 1c. The result survives passage to smooth compact original sources

Choose a real even nonnegative \(C_c^\infty(-1,1)\) unit-mass mollifier \(\rho_\eta\). Let \(f_\eta=\rho_\eta*f\) and \(g_\eta(u)=f_\eta(u-b_\eta)\) with \(b_\eta=L+2\eta\). Both lie in \(C_c^\infty\); their supports are disjoint up to boundary, so \(\langle f_\eta,g_\eta\rangle=0\) **exactly**. Their individual autocorrelations remain identical. As \(\eta\downarrow0\), all diagonal and symmetrized mixed autocorrelations converge to the above ramps; support is bounded for small \(\eta\), prime evaluations and two-pole integrals converge, and the full Gamma integral converges using the uniform Lipschitz bound near zero. Thus
\[
M_\eta\to8I_2,\qquad H_\eta\to H,
\]
so **for all sufficiently small \(\eta>0\)**, the smoothed sources are physically orthogonal, have **strictly negative original mixed signed Weil pairing**, and possess a **strictly positive two-dimensional signed restriction**. This is an *existence consequence*, not an explicit kernel-certified \(\eta\) or operator-completion theorem.

## 2. A counterintuitive positive result: rank-one physical metrics can always be fitted

For two abstract normalized rank-one arithmetic vectors with inner product \(\alpha\in\mathbb C\), \(r=|\alpha|\le1\), put
\[
z=\sqrt{16r^2+9}-4r\in[1,3],\qquad b=2\log z.
\]
Then \(9e^{-b/2}-e^{b/2}=8r\). By choosing the unit-modulus phase of the second translated source, the *normalized* physical Gram matrix of the two ramps can be made exactly \(\bigl[\begin{smallmatrix}1&\alpha\\\bar\alpha&1\end{smallmatrix}\bigr]\) (subject to the chosen inner-product convention), without changing **either diagonal** arithmetic cell observable.

This is a *universal two-rank-one metric fit*; it is **not evidence that the original arithmetic spaces coincide with the Weil space**, because a corresponding signed Hermitian numerator, integral lattice, character factors and support constraints are entirely unspecified. Indeed, §1 shows that physically orthogonal (\(\alpha=0\)) ramps can carry a nonzero negative original signed cross numerator.

This dual result sharpens the transfer gate: **matching a two-cell positive metric is relatively cheap; matching the complete signed mixed form is the genuinely discriminating obligation.**

## 3. Status of the requested *actual arithmetic* two-word test

The source-provided arithmetic \`Rows.Data.Bu\` is indexed by the record's \`Word\` type and maps \(\mathrm{Mu}(u)\times\mathrm{EE}\to\mathbb C\), with a global target \(\mathrm{EE}\). The typed matrix definition \`Kp\` shows what to calculate once two actual \`Word\` values and their normalized operators are specified:
\[
K_{uv}=\sqrt{\gamma_u\gamma_v}\,B_u^\dagger B_v,\qquad
\|K_{uv}\|_{\rm HS}^2,\quad
\|K_{uu}-I\|_{\rm HS}^2.
\]
**No specific instantiated two-word arithmetic operator pair with numeric cross HS data was extracted or calculated in this phase.** The paper's (S) and Lean's \`LinkS_at\` provide conditional *bounds* for those overlaps, not a concrete two-word matrix or a source-to-Weil intertwiner.

No evidence justifies taking the very explicit **prototype ramp cells** of §1 as the *actual theta-lift period operators*. The direct arithmetic-to-Weil **mixed-carrier** question remains **OPEN**, rather than falsely marked PASS. The rank-one fitting theorem in §2 is a feasibility test only.

## 4. Audit and next objective

**Executed local exact-rational audit:** \`src/audits/cr_alp_006_exact_cross_gram_audit.py\` (Python standard library, \`Fraction\`; all comparisons exact). It derives all 22 prime coefficients from prime powers \(n\le81\), encloses log terms by a rational atanh/geometric series, encloses \(\gamma\) by the harmonic-number inequalities and \(\pi\) by Archimedes' rational bounds, then certifies the mixed negative sign, diagonal positive sign, and the two-dimensional pencil lower bound. **PASS**. Decimal display is cosmetic only.

**Next CR-ALP-007 — source-defined mixed-operator extraction:** select an actually instantiated \`D7\` record and two \`Word\` indices with their natural \`Mu\` and \`EE\` metrics; compute their \`Kp\` / \`hsSqD\` overlap and pairing \`Bu\`, track the \`BesselId\` source maps, and ask whether the resulting signed target can be transported to **both** entries of a physical Weil Gram pencil. If only norm/trace bounds can be obtained, close the operator-conjugacy conjecture at this stage as **NOT ESTABLISHED**, retaining the valid restricted arithmetic-to-Weil reexpression of CR-ALP-004.

This source-only CR continuation changes no Weil-Lab, provenance or Hodge acceptance ledger.
