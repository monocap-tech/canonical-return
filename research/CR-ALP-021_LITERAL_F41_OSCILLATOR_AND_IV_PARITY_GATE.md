# CR-ALP-021 — Literal F41 oscillator attachment and original IV010 parity-block gate

**Canonical-Return, 2026-11-23 PST. Research-only draft.**  
**Standing:** EXACT frozen Lean source inspected; analytic parity/correlation and source-row derivations; standalone exact rational interval audit **LOCALLY EXECUTED: PASS**. **No new Lean compilation, new theorem in receiving Weil-Lab, F41 242×242 / 242×4 entry audit, original signed Weil whole-core positivity, support extension, or RH.**

## 0. Frozen source custody

This phase explicitly reopens and uses these **immutable original sources**:

- \`monocap-tech/weil-lab@a431f903254000bb5ad8edbc276ccc8599da98f5\`, \`WeilDefect/Arithmetic/ActualZetaOriginalJ242SourceFunctions.lean\` (blob **53935773feb3bba3c10e699cd8f862f9dd590b76**), defines all 242 literal \`nf545OriginalSelectedRealGenerator\` values: 241 rational-frequency/order sine-or-cosine monomials, followed by \`Real.sinh (x/2)\`.
- At the same commit \`ActualZetaOriginalNF407FrozenData.lean\` (blob **86ca70c8a2a8ba24ce4204081fef74d4ac543646**) gives first oscillator node \((37/64,0)\); \`ActualZetaOriginalNF407SourceBinding.lean\` (blob **342c6cffee6cee2719dfe35572a9e05aa8979726**) proves the exact frozen-node equality and defines four original critical generators (32 sine terms plus one \(\sinh(x/2)\)).
- \`ActualZetaOriginalPhysicalSolveResidual.lean\` (blob **2cf06688dc7e63bbb0c797102260559884325855**) defines **physical**, not signed Weil, \`neutralNF606PhysicalH\` (242×242 original selected \(L²\) Gram) and \`neutralNF606PhysicalV\` (242×4 selected/critical \(L²\) mixed integrals) on \((-11/10,11/10)\). Their historical numerical certification is separate and is not updated.
- Receiving branch \`research/iv068-original-source-kernel-gate\`: \`ActualZetaPhaseIVOriginalCrossSupportExtension.lean\` (blob **5ef51fd37337c466620eafccee2a27b77955746e**) defines authentic \`neutralIV010Mixed b f g = Zeta23.EF.literatureRHS (neutralWindowCorrelation b ...)\`; \`ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean\` (blob **14d8f01ce34f764900eccc68399d7323b4f4a25c**) proves the mixed-correlation raw Schwartz integral and decomposes the **original signed** form into **poles − paired prime powers + complete Gamma**.

No Lean source has been modified or recompiled in the receiving repository.

## 1. Actual IV010/IV011 parity-block identity: source-level proof

Let \(\mathsf Rf(x)=f(-x)\). Let \(f,o\in C_c^\infty((-B,B);\mathbb C)\), with \(f\) **even**, \(o\) **odd**, and embed these sources into the usual Schwartz class. The **literal IV011 physical mixed correlation** (the actual conjugation convention) is
\[
C_{f,o}(t)=\int_\mathbb R o(x)\,\overline{f(x-t)}\,dx.
\]
By substituting \(x\mapsto-x\),
\[
\boxed{C_{f,o}(-t)=-C_{f,o}(t),\quad C_{f,o}(0)=0.}\tag{21.1}
\]
There is no approximation or sign inference here.

Each of the **original** zeta explicit-formula components is an even distribution on this correlation input: the pole terms at \(\pm i/2\) occur in a pair, the von Mangoldt samples are explicitly \(C(\log n)+C(-\log n)\), and the genuine Gamma/digamma bracket is even in the Fourier variable. Consequently the literal signed original pair vanishes:
\[
\boxed{\operatorname{neutralIV010Mixed}(B,f,o)=0.}\tag{21.2}
\]
The same is true reversing the inputs by Hermitian symmetry. This applies to smooth sources for which the authentic IV011 mixed formula is admissible; **not** automatically to the sharply restricted analytic NF545 functions as signed-Weil domain vectors.

**Lean formalization gate (NOT compiled here):** Prove reflection identities for \`neutralWindowCorrelation\` under \`toLp\` from the existing \`neutralIV011MixedCorrelation_eq_schwartzIntegral\`, then show the pole, paired-prime and \`gammaBracket\` contributions each vanish for odd input correlation. The receiving library must supply or prove \`gammaBracket (-r)=gammaBracket r\`, the Fourier parity/integrability lemmas, and window support compatibility. A theorem claiming \`neutralIV010Mixed\` parity before these obligations are discharged would be unverified.

## 2. A STRICTLY NONZERO attachment to an actual frozen F41 oscillator

CR-ALP-020 established that the **evenized** packet sources \(E_p=h_p+\mathsf Rh_p\), \(p=2,3\), have exactly zero physical inner products with all frozen 242 selected and four critical odd generators. That does not imply that the **raw** \(h_p\) have zero moments. Here we certify an exact nonzero **actual first** selected row, without inventing a new F41 column.

Fix
\[
b=10^{-8},\quad \eta=10^{-14},\quad\ell=\log2,\quad
\omega=2\pi(37/64)=37\pi/32.
\]
Use CR-ALP-020's genuine centered raw packet
\[
\begin{gathered}
g_b(t)=e^{t/2}\mathbf1_{[0,b]}(t)-e^{t/2-b}\mathbf1_{[b,2b]}(t),\\
\widetilde h_2(x)=g_b(x+\ell+b)-2g_b(x-\ell+b),\qquad
h_2=\rho_\eta*\widetilde h_2,
\end{gathered}
\]
for **any** real even nonnegative unit-mass \(\rho_\eta\in C_c^\infty(-\eta,\eta)\).

The **literal frozen** first source row is
\[
\boxed{\varphi_0(x)=\sin(\omega x),\quad\omega=37\pi/32.}\tag{21.3}
\]
It is precisely \`nf545OriginalSelectedRealGenerator (0:Fin 242)\` by the pinned frozen node \((37/64,0)\). The original physical integral is over \((-11/10,11/10)\), but \(h_2\) is supported strictly inside.

Set \(\theta=\omega\ell\), \(u=t-b\). From exact rational atanh series bounds \(2/3<\log2<3/4\),
\[
\boxed{2\pi/3<\theta<\pi,\quad \cos\theta<-1/2,\quad\sin\theta>0.}\tag{21.4}
\]
Use the *actual two subpacket values*, for \(0<u<b\),
\[
g_b(b-u)=e^{(b-u)/2},\qquad
g_b(b+u)=-e^{-(b-u)/2}.
\]
Define
\[
\begin{aligned}
S_b&=\int_{-b}^{b}g_b(b+u)\sin(\omega u)\,du
=-2\int_0^b\cosh((b-u)/2)\sin(\omega u)\,du,\\
C_b&=\int_{-b}^{b}g_b(b+u)\cos(\omega u)\,du
=2\int_0^b\sinh((b-u)/2)\cos(\omega u)\,du.
\end{aligned}
\]
The exact rational \(\pi\)-enclosure \(3<\pi<22/7\) yields \(3<\omega<4\) and \(\omega b<1\). By \(\sin y\ge y/2\) on \(0\le y\le1\), \(\cosh\ge1\), and \(\cos(\omega u)>0\) here,
\[
\boxed{S_b\le-\omega b^2/2<0,\qquad C_b>0.}\tag{21.5}
\]
The **exact unsmoothed selected moment** is
\[
\begin{aligned}
v_{0,2}^{\mathrm{raw}}
&=\int\widetilde h_2(x)\sin(\omega x)dx\\
&=-\cos\theta\,S_b-3\sin\theta\,C_b
<-\omega b^2/4<-3b^2/4.
\end{aligned}
\tag{21.6}
\]

**No mollifier specialization is needed.** Because \(\rho_\eta\) is real, even and nonnegative,
\[
v_{0,2}=\int h_2(x)\varphi_0(x)dx=
\mu_\rho(\omega)\,v_{0,2}^{\rm raw},
\quad
\mu_\rho(\omega)=\int \rho_\eta(y)\cos(\omega y)dy.
\]
Since \(\omega\eta<1\), \(\mu_\rho(\omega)>1/2\) by \(\cos y\ge1-y^2/2\). Thus
\[
\boxed{
v_{0,2}<-\frac38\,b^2
=-\frac{3}{8\cdot10^{16}}<0.
}\tag{21.7}
\]
This is an **actual frozen NF545 first physical row** with a certified strict sign for **every** permitted mollifier. It is not a computed new entry of the historical F41 242×4 \`V\`, because that frozen matrix uses NF587 critical generators, not \(h_2\).

### 2a. The source has both an inside AND an outside component

Let \(S_{242}\) be the physical \(L²(-11/10,11/10)\) span of the frozen selected 242 modes. Since \(\varphi_0\in S_{242}\) and
\(\|\varphi_0\|_2^2\le11/5\), Cauchy and (21.7) yield the **nonzero selected projection lower bound**
\[
\boxed{\|P_{S_{242}}h_2\|_2^2>
\frac{45}{704}b^4>0.}\tag{21.8}
\]
Independently CR-ALP-020 established
\[
\boxed{\operatorname{dist}(h_2,S_{242})^2\ge\tfrac1{10}\|h_2\|_2^2,\quad
\operatorname{dist}(h_3,S_{242})^2\ge\tfrac1{5}\|h_3\|_2^2.}\tag{21.9}
\]
These statements are **compatible**: our raw \(h_2\) lies partly within and partly outside the actual finite odd F41 selected span. Nothing implies the frozen 242 basis spans arbitrary even functions or that it is complete in the full odd infinite-dimensional physical space.

The selected pole row \(j=241\) has zero moment against \(h_2\) by its exact two original Mellin pole zeros and the identity \(\sinh(x/2)=(e^{x/2}-e^{-x/2})/2\); the first oscillator row, by contrast, is strictly negative. This distinguishes a *true* nonzero attachment from a trivial nonvanishing of the packet norm.

## 3. Reproducibility and disposition

Standalone source-specific auditor \`src/audits/cr_alp_021_literal_f41_oscillator_sign_audit.py\` was **executed locally: PASS**. It uses only \`fractions.Fraction\` for a rigorous rational logarithm enclosure via atanh, Machin alternating rational \(\pi\) bounds, the strict angle inequalities, the rational coefficient/mollifier inequalities and the selected projection floor. It does **not** claim to computer-prove the continuum sign identities or numerically sample a mollifier—the human proof above provides those.

| Obligation | CR-ALP-021 |
|---|---|
| Frozen first F41 oscillator source and true frequency | **SOURCE-PINNED** |
| Original IV010/IV011 literal correlation / signed pole-prime-Gamma split | **SOURCE-PINNED** |
| Smooth even/odd literal mixed Weil cross parity | **ANALYTIC EXACT; Lean formalization OPEN** |
| New even packet physical mixed rows against all 246 original odd generators | **EXACT ZERO (CR-ALP-020 retained)** |
| At least one actual raw packet *nonzero* frozen selected row | **\(v_{0,2}<-(3/8)b^2\), analytic/rational PASS** |
| Projection of \(h_2\) onto actual frozen selected span is nonzero | **PASS with absolute floor \(45b^4/704\)** |
| Original raw packet has nonzero physical complement | **10%/20% squared-mass floors retained** |
| Full 242-row numeric interval enclosure / historical V rewrite | **NOT PERFORMED** |
| New receiving Lean theorem compiled / axiom audit | **NOT PERFORMED** |
| Full even+odd old-domain Schur, original RH / global positivity | **OPEN** |

**Next — CR-ALP-022:** formalize the exact reflection transform of \`neutralWindowCorrelation\` and the evenness of literal \`Zeta23.EF.gammaBracket\` on the receiving Lean IV010/IV011 Schwartz core. Separately build a *fixed explicit mollifier* with certified interval evaluator for several nonzero original F41 oscillator rows, while respecting frozen source custody. No historical F41 numerical H/V or RH ledgers may be overwritten without a separate source/evidence gate.

**Governance:** Canonical-Return draft research only. No writes to Weil-Lab, F41/F42, original RH/provenance/Hodge acceptance states.
