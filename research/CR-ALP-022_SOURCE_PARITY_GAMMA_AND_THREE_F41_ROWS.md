# CR-ALP-022 — Complete original IV Weil even/odd parity interface and THREE source-exact NF545 oscillator signs

**Weil / Canonical-Return independent research branch, 2026-11-23 PST.**  
**Standing:** direct immutable Lean **source inspection**, new **human analytic proofs**, independently executed **exact Fraction/SymPy audit PASS**. **NOT** a Lean compilation, new receiving IV010/IV011 theorem, frozen 242-row numerical F41 acceptance, all-old-core support extension, or RH proof.

## 0. Primary source custody

Reopened exact sources and blob anchors (do not replace the historical physical data):

**Frozen F41 revision:** \`monocap-tech/weil-lab@a431f903254000bb5ad8edbc276ccc8599da98f5\`:
- \`WeilDefect/Arithmetic/ActualZetaOriginalJ242SourceFunctions.lean\` blob \`53935773feb3bba3c10e699cd8f862f9dd590b76\`: exact \`nf545OriginalSelectedRealGenerator\`; the first frequency group has nodes \((37/64,k)\), \(k=0,\ldots,13\).
- \`WeilDefect/Arithmetic/ActualZetaOriginalNF407FrozenData.lean\` blob \`86ca70c8a2a8ba24ce4204081fef74d4ac543646\`: explicitly stores the frozen node sequence.
- \`WeilDefect/Arithmetic/ActualZetaOriginalNF407SourceBinding.lean\` blob \`342c6cffee6cee2719dfe35572a9e05aa8979726\`: the four frozen critical generators retain their exact 32 sine terms and \`sinh(x/2)\).
- \`WeilDefect/Arithmetic/ActualZetaOriginalPhysicalSolveResidual.lean\` blob \`2cf06688dc7e63bbb0c797102260559884325855\`: frozen \`neutralNF606PhysicalH\` and \`V\` are **physical \(L^2\)** moment matrices, NOT signed Weil matrices.

**Receiving actual original IV source:** \`monocap-tech/weil-lab\`, branch \`research/iv068-original-source-kernel-gate\`:
- \`WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportExtension.lean\` blob \`5ef51fd37337c466620eafccee2a27b77955746e\`: \`neutralIV010Mixed B f g\` is \`Zeta23.EF.literatureRHS (neutralWindowCorrelation B ...)\`.
- \`WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean\` blob \`14d8f01ce34f764900eccc68399d7323b4f4a25c\`: \`neutralIV011MixedCorrelation_eq_schwartzIntegral\` has actual first-slot conjugation \(\int t(x)\overline{s(x-\tau)}\,dx\); \`neutralIV011OriginalMixedPrime\` samples at **both** \(\pm\log n\); the full decomposition is **pole − prime + arch**.
- \`WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean\` blob \`0034976fe2f1e5fa9d411108736160542eb94293\`: the ACTUAL \`Zeta23.EF.gammaBracket(r)\` is \(\operatorname{Re}\psi(1/4+ir/2)-\log\pi\), and \`literatureRHS\` carries both pole evaluations, both prime signs, and the full archimedean integral.
- \`WeilDefect/External/Zeta23/Theorems/Thm_Zeta23_WeilEF_digamma_conj.lean\` blob \`27f17b3b6e522f59d663486e851921b5f4c7dedd\`: pre-existing theorem \`Zeta23.WeilEF.digamma_conj\` on \`Complex.integerComplement\`.
- \`WeilDefect/External/Zeta23/Theorems/Thm_Zeta23_WeilEF_gammaR_bracket.lean\` blob \`9b43c6390b38b3377da7a8b01c65b86af4d6fba4\`: pre-existing \`Zeta23.WeilEF.gammaR_bracket\`, which also proves Gamma bracket evenness by swapping its \(\pm it\) terms.

These are **source declarations**, not a newly accepted compiled parity theorem. The local container has no \`lean\` or \`lake\` binary.

## 1. Complete parity gate: exact mathematical proof and actual Lean attachment

Let \(B>0\); let \(f,o\in C_c^\infty((-B,B);\mathbb C)\) be **even** and **odd**, respectively, and regard them as Schwartz sources. The literal IV011 correlation (the actual original convention) is
\[
C_{f,o}(\tau)=\int_\mathbb R o(x)\,\overline{f(x-\tau)}\,dx.
\]
Substitute \(x\mapsto-x\), use \(o(-x)=-o(x)\) and \(f(-x)=f(x)\):
\[
\boxed{C_{f,o}(-\tau)=-C_{f,o}(\tau).}\tag{22.1}
\]
The original Weil distribution is **even** in its correlation argument. Verify all **three** source-defined terms separately:

1. **Prime (already a direct source-level vanishing):** Every summand of \`neutralIV011OriginalMixedPrime\` is proportional to \(C(\log n)+C(-\log n)=0\). Therefore \(\boxed{\mathcal P_{\rm mixed}(f,o)=0}\). No prime cutoff or PNT estimate is required.
2. **Two poles:** By direct variable reflection, the genuine \`Zeta23.paperFT\` of an odd \(C\) satisfies \(h(-z)=-h(z)\) for every complex \(z\). Thus \(\boxed{h(i/2)+h(-i/2)=0}\).
3. **Full Gamma:** The ACTUAL \`EF.gammaBracket(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi\) is an even **real** function. This follows from the existing \`WeilEF.digamma_conj\` because \(1/4+ir/2\) has positive real part, away from integer poles. Alternatively, the already proved \`WeilEF.gammaR_bracket\` expresses it as the symmetric \(\Gamma_\mathbb R\) log-derivative sum at \(1/2\pm ir\). For odd compact \(C\), \(h(-r)=-h(r)\) on the real line, so the full convergent Gamma integral \(\int_\mathbb R h(r)\,\gammaBracket(r)\,dr\) equals its negative and is zero.

Therefore, **without replacing the original signed arithmetic distribution**:
\[
\boxed{\mathrm{neutralIV010Mixed}(B,f,o)=0,\quad
\mathrm{neutralIV010Mixed}(B,o,f)=0.}\tag{22.2}
\]
The second follows also from the receiving \`neutralIV010Mixed_hermitian\` theorem when both supports fit \(B\).

### Lean kernel formalization gate, explicitly NOT passed

A *type-exact target* for the receiving branch is:
\[
\forall B\ge0,\ f,o:\mathcal S(\mathbb R,\mathbb C),
\ \bigl[\mathrm{supp}(f),\mathrm{supp}(o)\subset[-B,B],\
f(-x)=f(x),\ o(-x)=-o(x)\bigr]
\Longrightarrow \mathrm{neutralIV010Mixed}(B,f,o)=0.
\]
The proof dependency order is:
- **L1:** instantiate \`neutralIV011MixedCorrelation_eq_schwartzIntegral B B (le_refl B) ...\`, prove the reflection formula in the genuine windowed Schwartz/L² carrier using \`integral_neg_eq_self\` or the measure-preserving negation substitution;
- **L2:** prove oddness of \`paperFT\` of the odd correlation for **complex** Fourier arguments and its integrability at \(\pm i/2\);
- **L3:** make \`gammaBracket (-r)=gammaBracket r\` explicit from already compiled \`WeilEF.digamma_conj\` (or \`gammaR_bracket\`), including the \`integerComplement\` side condition \(\operatorname{Re}(1/4+ir/2)=1/4\);
- **L4:** show \`neutralIV011OriginalMixedArch=0\` by real-line \(\tau\mapsto-\tau\) integral parity; ensure the original Fourier integral and Gamma bracket product is integrable;
- **L5:** rewrite with \`neutralIV011OriginalMixed_arithmeticSplit\`, and use receiving Hermitian symmetry for reversed inputs.

No \`sorry\`, invented axiom, or additional assumption is inserted into the real receiving Lean repository in this research phase. **No kernel compilation is claimed.** The above is a proof roadmap with a full human proof.

## 2. Canonical explicit mollifier and THREE certified nonzero frozen oscillator rows

CR-ALP-021 certified only the frozen first selected row against the unsymmetrized \(p=2\) packet. Here we obtain **THREE** independently signed source-exact rows.

Choose the *specific reproducible* smooth bump
\[
\boxed{
\rho_\eta(t)=
\frac{\mathbf 1_{|t|<\eta}\exp\bigl[-1/(1-(t/\eta)^2)\bigr]}
{\displaystyle\int_{-\eta}^{\eta}\exp\bigl[-1/(1-(v/\eta)^2)\bigr]\,dv},
\quad \eta=10^{-14}.}
\tag{22.3}
\]
It is real, even, positive on \((-\eta,\eta)\), compactly supported and \(C^\infty\). **No numerical quadrature of its normalizing integral is needed:** the strict moment certificates below hold for **every** even nonnegative mass-one mollifier of at most this width, hence in particular (22.3).

Let \(b=10^{-8}\), \(\ell=\log2\), \(\omega=37\pi/32\), and
\[
g_b(t)=e^{t/2}\mathbf1_{[0,b]}(t)-e^{t/2-b}\mathbf1_{[b,2b]}(t),\quad
\widetilde h_2(x)=g_b(x+\ell+b)-2g_b(x-\ell+b),\quad
h_2=\rho_\eta*\widetilde h_2.
\]
The first *three exact frozen NF545 generators* (rows \(j=0,1,2\)) are
\[
\boxed{\varphi_0(x)=\sin(\omega x),\quad
\varphi_1(x)=x\cos(\omega x),\quad
\varphi_2(x)=x^2\sin(\omega x).}\tag{22.4}
\]
This is literal frozen NF407 node \((37/64,k)\), \(k=0,1,2\); **not** a newly chosen basis.

### 2a. Exact first-order packet expansion with an explicit rigorous error

For any real odd \(C^2\) function \(F\) on the packet neighborhoods,
\[
v_F^{\mathrm{raw}}
=\int_{-b}^{b}g_b(b+u)\,\left[-F(\ell-u)-2F(\ell+u)\right]\,du.
\]
Write \(m_0=\int_{-b}^b g_b(b+u)\,du\) and \(m_1=\int_{-b}^b u\,g_b(b+u)\,du\). From the **exact two packet values** \(g_b(b-u)=e^{(b-u)/2}\), \(g_b(b+u)=-e^{-(b-u)/2}\),
\[
\begin{aligned}
m_0&=2\int_0^b\sinh((b-u)/2)\,du,\qquad
\frac{b^2}{2}\le m_0\le\frac{b^2}{2}+\frac{b^4}{8},\\
m_1&=-2\int_0^b u\cosh((b-u)/2)\,du,\qquad
-b^2-\frac{b^4}{4}\le m_1\le-b^2.
\end{aligned}
\tag{22.5}
\]
These elementary bounds use \(1\le\cosh(b/2)\le1+b^2/4\) for \(b\le1\).

Taylor expanding the bracket at \(u=0\) and using the uniform bounds \(|F|\le2\), \(|F'|\le8\), \(|F''|\le50\) for all three *frozen functions* on \(|x|\le1\), together with \(\int u^2|g_b(b+u)|du\le(4/3)b^3\), gives
\[
\boxed{\left|v_F^{\mathrm{raw}}
-b^2\left(F'(\ell)-\frac32F(\ell)\right)\right|<103b^3.}
\tag{22.6}
\]
Convolution with any permitted \(\rho_\eta\) changes the moment by no more than \(72b\eta\), because \(\|\widetilde h_2\|_1<9b\) and \(\sup|F'|\le8\). Thus
\[
\boxed{
\left|v_F-
b^2\left(F'(\ell)-\frac32F(\ell)\right)\right|
<103b^3+72b\eta<10^{-4}b^2.}
\tag{22.7}
\]

### 2b. Three strict source-indexed negative signs

The executed exact rational auditor derives \(2/3<\ell<3/4\), \(3<\omega<4\), and
\[
\boxed{\frac{3\pi}{4}<\theta:=\omega\ell<\frac{7\pi}{8},\quad
\cos\theta<-\frac7{10},\quad
\sin\theta>\frac38.}
\tag{22.8}
\]
The strict sine bound uses \(\sin^2(\pi/8)=(2-\sqrt2)/4>9/64\).

Set \(A_j=\varphi'_j(\ell)-(3/2)\varphi_j(\ell)\). Exact differentiation and elementary monotonicity give
\[
\begin{aligned}
A_0&=\omega\cos\theta-\tfrac32\sin\theta<-\tfrac{213}{80},\\
A_1&=(1-\tfrac32\ell)\cos\theta-\theta\sin\theta<-\tfrac58,\\
A_2&=(2\ell-\tfrac32\ell^2)\sin\theta+\omega\ell^2\cos\theta<-\tfrac4{15}.
\end{aligned}
\tag{22.9}
\]
Combining (22.7) with these strict margins proves the **new three-row certificate**:
\[
\boxed{
\begin{aligned}
v_{0,2}:=\int_{-11/10}^{11/10}h_2(x)\varphi_0(x)\,dx&<-2b^2,\\
v_{1,2}:=\int_{-11/10}^{11/10}h_2(x)\varphi_1(x)\,dx&<-\tfrac12b^2,\\
v_{2,2}:=\int_{-11/10}^{11/10}h_2(x)\varphi_2(x)\,dx&<-\tfrac15b^2.
\end{aligned}}
\tag{22.10}
\]
All are **strictly nonzero** for the specific standard mollifier (22.3) and, in fact, uniformly over all even nonnegative unit-mass mollifiers at this support width. This substantially strengthens CR-ALP-021's \(j=0\) result \(\int h_2\varphi_0<-3b^2/8\).

**Type note:** These are **three NEW physical mixed moments between our packet and three real frozen F41 selected oscillator functions**. They are not entries of the historical 242×4 \`neutralNF606PhysicalV\`, which is defined against the **four different NF587 critical generators**. The historical physical H/V numerical acceptance remains untouched.

## 3. Certification status

| Obligation | CR-ALP-022 |
|---|---|
| Literal IV010/IV011 original signed correlation and prime/pole/Gamma conventions | **SOURCE-INSPECTED** |
| Source-authored digamma conjugation / Gamma bracket symmetry ingredients | **EXIST in pinned Lean source** |
| Full even/odd mixed signed Weil pairing | **HUMAN ANALYTIC PASS; not Lean compiled** |
| Literal receiving Lean theorem / \`#print axioms\` | **NOT ESTABLISHED** |
| Canonical explicit smooth mollifier | **SPECIFIED**, positive, even and Cc∞ |
| Exact source-indexed F41 physical rows j=0,1,2 | **THREE strict negative signs, human analytic PASS** |
| Self-contained rational algebra auditor | **EXECUTED LOCAL PASS** |
| Whole 242-entry mixed attachment / frozen historical F41 H/V | **NOT PERFORMED / NOT MODIFIED** |
| Entire old Schwartz collar Schur / arbitrary support extension | **OPEN** |
| RH / global signed Weil positivity / full operator domain | **NOT CLAIMED** |

**Reproducibility:** \`src/audits/cr_alp_022_three_frozen_rows_exact_audit.py\`, using Python \`Fraction\` and SymPy, checks exact atanh \(\log2\), Machin \(\pi\), the frozen frequency, three derivative formulae, rigorous uniform Taylor/mollifier budgets, and all three final inequalities. PASS. The parity and analytic Taylor premises are displayed as human proofs here; no computer-verified integral theorem is claimed.

**Next — CR-ALP-023:** attempt to discharge the receiving \`gammaBracket\` evenness and \`paperFT\` oddness as *actual* Lean lemmas against the pinned IV source, with an explicit kernel/axiom report if a working Lean environment becomes available. Extend the per-row moment auditor to several additional frozen frequency groups using exact symbolic derivatives and certified interval arithmetic; do not conflate new source moments with the frozen F41 numerical H/V matrices. The essential infinite-dimensional signed Weil Schur/support-extension bound remains open.

**Governance:** canonical-return draft/nonpromoting. No write to Weil-Lab, F41/F42, historical RH/provenance or Hodge acceptance ledgers.
