# CR-ALP-005 — Frame-potential comparison and the missing mixed-metric carrier

**Jurisdiction:** Canonical-Return / Alpöge branch.  
**Date:** 2026-11-23 (America/Los_Angeles).  
**Claim state:** SOURCE-PINNED + DERIVED / EXACT FINITE MODEL + AUDIT-PASS (SymPy). **No new Lean compilation, no theorem about Alpöge's actual theta lifts, no noncircular proof of the arithmetic inequality (M), no operator equivalence, no RH advance.**

## 0. Source custody and the question tested

This stage asks whether the **rank/trace mechanism** in the Alpöge formal proof has a mathematical quantity compatible with the original physical Weil Gram pencil, and whether the CR-ALP-003 *cell-diagonal source dictionary* contains enough information to reconstruct that quantity.

The answer is **mixed**: there is an exact common **finite frame-potential identity**, but **the existing arithmetic table cannot determine the mixed Gram data needed to use it**. We provide an explicit one-parameter family showing that information loss. This is a no-go for an *inference from cellwise data alone*, **not** a no-go for a richer joint theta-lift/Weil correspondence.

Sources reopened:

- Levent Alpöge, *A remark*, 20 November 2026, https://alpo.ge/remark.pdf, printed pp. 3–4: determinant/Weil-positivity experiments (timeline items 4, 6), switch to the prime ramp (item 5), iterated theta lifts, norms, divisibility, Bost slopes and relative trace formula (items 8–18).
- Rubinstein et al., *The Alpöge–Claude proof of the Riemann hypothesis explained for the Riemann zeta function*, **v2 working draft, 22 November 2026**, https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf, **§III.8, pp. 136–144**, especially (III.217)–(III.227), (III.233), (III.241); §III.11, (H8). This exposition explicitly distinguishes proved linear-algebra reductions from three **quoted** trace/rank/pairing statements and other imported arithmetic construction inputs.
- Source-checked public Lean mirror, immutable ref **7cbf7a0c32fca1badf05931879a37c3f8a2b03b4**, https://github.com/CocoSgt/riemann-hypothesis-lean, with exact paths:
  - src/Formalism/SharpRestriction/GlobalTheoremInputs.lean — definitions \`LinkS_at\`, \`LinkS\`, \`LinkIRR\`, \`BesselId\`;
  - src/Formalism/RankGapBounds/MeanBoundSecondForm.lean — \`wmean\`, \`Data.s\`, \`Data.SlopeIneq\`, \`Data.MeanBoundSecondForm\`;
  - src/Formalism/SizeWindowBrackets/ChainSlopeWindows/LevelOfCongruence.lean — \`each_unit_bound_statement\`, \`weighted_mean_bound_statement\`, \`uniform_bound_statement\`;
  - src/Formalism/SizeWindowBrackets/ChainSlopeWindows/PChain.lean — \`chain_top_bound\`.

The paper and Lean source have different proof presentations. In particular, Rubinstein's §III.8 **Corollary T2** is an independent human-readable reconstruction, not a declaration with the same name in the Lean repo. The source uses other statements and the \`FrameLA1\` / \`TheoremF\` pathway. Do not silently identify them.

## 1. The actual rank/trace geometric mechanism

Let \(\mathcal H\) be a finite-dimensional complex Hilbert space of dimension \(D>0\), and let \(B_w:C_w\to\mathcal H\) be linear operators. Set
\[
B=\bigoplus_w B_w,\quad Q=BB^*=\sum_w B_wB_w^*,\quad M=B^*B.
\]
Then \(Q\ge0\), \(M\ge0\), and the exact mixed-block identities are
\[
\boxed{
\operatorname{tr}Q=\sum_w\|B_w\|_{\rm HS}^2=\operatorname{tr}M,
\qquad
\operatorname{tr}(Q^2)=\sum_{w,w'}\|B_w^*B_{w'}\|_{\rm HS}^2=\operatorname{tr}(M^2).
}
\tag{5.1}
\]
(The cross-Gram block is \(B_w^*B_{w'}:C_{w'}\to C_w\); use this type-correct orientation when checking source conventions.)

For \(\operatorname{tr}Q>0\), define the **frame isotropy excess relative to the specified target dimension**:
\[
\boxed{\Phi_D(Q)=\frac{D\,\operatorname{tr}(Q^2)}{(\operatorname{tr}Q)^2}-1
=\frac{D\left\|Q-\frac{\operatorname{tr}Q}{D}I\right\|_{\rm HS}^2}{(\operatorname{tr}Q)^2}\ge0.}
\tag{5.2}
\]
This is an exact finite linear-algebra identity. It is invariant under **unitary** changes of the target Hilbert coordinates, and under unitary changes of source coordinates that preserve the operators' domains. It is **not** invariant under arbitrary nonunitary source rescalings or arbitrary Hermitian-pencil congruences unless the source metric is transported and the frame operators are renormalized accordingly.

The Cauchy–Schwarz inequality for the nonzero eigenvalues of \(Q\) yields
\[
\boxed{\operatorname{rank}Q\ge
\frac{(\operatorname{tr}Q)^2}{\operatorname{tr}(Q^2)}
=\frac{D}{1+\Phi_D(Q)}.}
\tag{5.3}
\]
This is a **rank/capacity** conclusion, not a sign bound on a different Hermitian Weil numerator.

Rubinstein §III.8, Lemma T1(c), under the stated trace-diagonal and dimension hypotheses, obtains \(\operatorname{tr}(Q^2)\le(1+5\varkappa)(\operatorname{tr}Q)^2/D\), i.e. \(\Phi_D(Q)\le5\varkappa\). Alongside deletion bounds and the pairing/integrality hypotheses, this supports the high-rank image-space control entering the mean-slope inequality. It does **not** say the different images are orthogonal. The source explicitly says they need not be. For the actual arithmetic \(B_w\), the key inputs \`LinkS\` (off-diagonal HS overlap and diagonal defect), \`LinkIRR\` (rank loss), and \`BesselId\` (pairing attachment) operate in a **different** Hilbert carrier from the logarithmic \(L^2\) tests of CR.

**Important level distinction:** the finite-frame identity (5.1) is a universal linear-algebra theorem. The actual quantitative bounds \`LinkS\`, \`LinkIRR\` and \`BesselId\` belong to the constructed arithmetic family; they are not consequences of (5.1) or of the Weil prime-ramp sampling identity.

## 2. Exact no-go: cellwise arithmetic observables do not determine mixed source Gram data

CR-ALP-003 represents each positive-part arithmetic cell weight \((a/n-w)_+\) as an autocorrelation sample of a logarithmic \(L^2\) source. The representation is **not unique as a joint Hilbert-space realization of many cells**.

Take the already-used elementary cell example \(a=9,\ w=1,\ L=\log9\), and
\[
f(u)=e^{u/2}\mathbf1_{[0,L]}(u),\qquad
f_b(u)=f(u-b),\qquad b\ge0.
\]
For **every** \(b\), both individual cells have exactly the same autocorrelation \(k=f*\widetilde f=f_b*\widetilde f_b\); therefore their entire diagonal original Weil tests \(\mathcal W(k)\), prime samples, diagonal Gamma and pole values, and diagonal physical masses are identical. In particular,
\[
\|f\|_2^2=\|f_b\|_2^2=8,\qquad
\frac{k(\log n)}{\sqrt n}=(9/n-1)_+
\]
for all \(n\ge1\), independent of \(b\).

But their **cross inner product** varies:
\[
\boxed{
\langle f,f_b\rangle=
\begin{cases}
9e^{-b/2}-e^{b/2},&0\le b\le\log9,\\
0,&b\ge\log9.
\end{cases}}
\tag{5.4}
\]
The two-source *native physical* Gram matrix is
\[
M(b)=\begin{pmatrix}8&c_b\\c_b&8\end{pmatrix},
\qquad c_b=\langle f,f_b\rangle.
\tag{5.5}
\]
The Gram diagonals, and every separate arithmetic/Weil cell observable, are **independent of \(b\)**; the mixed metric \(M(b)\) is not.

Exact cases:
- \(b=0:\ c_b=8\), rank \(M=1\), \(\Phi_{D=2}=1\) for a fixed two-dimensional ambient target.
- \(b=\log3:\ c_b=2\sqrt3\), rank \(M=2\), \(\det M=52\), and \(\boxed{\Phi_2=3/16}\).
- \(b=\log9:\ c_b=0\), rank \(M=2\), \(\boxed{\Phi_2=0}\).

The value \(\Phi_2=3/16\) is computed exactly from \(2\operatorname{tr}(M^2)/(\operatorname{tr}M)^2-1\). This is a **strict information-loss witness**: the cellwise samples and even the complete **diagonal** signed Weil values fail to fix the pairwise Gram overlaps, ranks, or frame potential. All arithmetic cell rows remain unchanged.

### Smooth-source survival

Let \(\rho_\delta\) be a compact real mollifier as in CR-ALP-003 and put \(g_\delta=\rho_\delta*f\). The translated functions \(g_{\delta,b}(u)=g_\delta(u-b)\) have **the same individual autocorrelation for every \(b\)**, including identical complete signed Weil diagonal values. At \(b=0\) the pair is equal (rank one), while for \(b>L+2\delta\) supports are disjoint (rank two, zero cross mass). For \(b=\log3\) and \(\delta\downarrow0\), the cross mass converges to \(2\sqrt3\). Thus the ambiguity is **not a nonsmooth artifact**. No statement about completed Suzuki/Connes operator domains is deduced.

**Scope of the no-go.** Individual *source* pole moments \(F_\pm(f_b)\) do scale by \(e^{\pm b/2}\); the **diagonal autocorrelation pole term** and full diagonal Weil test remain invariant because the product is unchanged. A richer dictionary retaining individual pole moments and all **mixed** data may detect the relative translations. Our conclusion concerns deductions from **cellwise scalar autocorrelation/table data alone**.

## 3. Even perfectly isotropic frames do not determine the signed numerator

Consider a purely abstract 2D finite Hilbert model with synthesis map \(B=I_2\) and metric \(M=I_2\). Then \(Q=I_2\) and \(\Phi_2(Q)=0\), the strongest possible frame-isotropy statement. But the same \(B\) and same \(M\) admit distinct Hermitian form numerators, e.g.
\[
H_+=\operatorname{diag}(1,1),\qquad
H_-=\operatorname{diag}(-1,1).
\]
The first quadratic form is positive and the second is indefinite. Hence the frame-moment bound by itself says **nothing about the sign of an independently specified Weil numerator**. This is an abstract falsifier for a logical implication; it is **not** a counterexample to positivity of the actual zeta Weil form.

The actual arithmetic mechanism obtains a mean-slope/divisibility inequality using its frame/rank data **plus additional arithmetic and lattice hypotheses**. The independent Weil positivity criterion requires a separately attached numerator, native metric, full poles/prime/Gamma and support control. No cross-theory equality of the required forms is known.

## 4. Lean/working-draft correspondence: what must travel together

| Mathematical requirement | Actual source anchor | Required CR transfer |
|---|---|---|
| Diagonal and off-diagonal HS moments | \`LinkS_at\` in \`GlobalTheoremInputs.lean\`; manuscript §III.8 (S), (III.219)–(III.221) | Full **mixed** inner-product matrix, not only each \(k_{a,w}\) |
| Rank-loss / retained image dimension | \`LinkIRR\`; manuscript (I), (III.222) | Correct common target dimension and projective losses |
| Actual pairing attachment | \`BesselId\`; manuscript (P), (III.223) | A real source-defined intertwiner; not redefinition of \(W\) |
| Weighted slopes, height, divisibility | \`Data.s\`, \`Data.SlopeIneq\`; manuscript (III.241), (H8) | Arithmetic lattice plus a genuine pairing to the same map |
| Original signed Weil form and native metric | CR-ALP-002–004; CR-T1B/T1C | \(H_{\rm Weil}\) **and** \(M_{\rm physical}\), both poles, Gamma and prime phases |
| \(x\)-uniform analytic control | CR-ALP-004 smoothing (source-level) | Support \(a_x\sim\tfrac12\log x\) plus mixed/complement control |

The manuscript §III.8 and §III.11 expressly flag that the arithmetic pairings used in the lattice, lift/size, and trace accounts must be shown to correspond to **one and the same map with consistent norms**; its discussion of O25(e) says that identification is not checked in the **human-readable reconstruction**. This is a **review obligation**, not a declaration that the independently checked Lean proof is invalid.

## 5. Disposition and next independent theorem candidate

**PASS — exact finite algebra:** typed frame traces, nonnegative isotropy excess and rank lower bound; explicit translated-ramp family and exact mixed-Gram variation; qualitative smooth persistence; abstract independent-numerator countermodel.

**PASS — source anchored:** manuscript's three named arithmetic rank/trace/pairing statements and the corresponding Lean definition locations at a frozen mirror. Exact source snippets inspected; CR did **not** run Lean.

**NO TRANSFER from the existing cell-table dictionary alone:** it fails to determine the required mixed Gram metric or frame potential; positive source kernels and diagonal entries do not imply the theta-lift Hilbert mean/slope estimate or signed Weil positivity.

**OPEN:** a *joint* source-defined map from the actual arithmetic \(B_w\) family to a specified Weil test carrier, preserving (i) target/source native metrics, (ii) the mixed blocks \(B_w^*B_{w'}\), (iii) the correct **signed** Weil numerator including poles, Gamma and characters, (iv) the arithmetic integral lattice/slope data, and (v) uniform control as support grows. A completely abstract congruence is insufficient.

**Next CR-ALP-006 — mixed-pairing transport protocol.** Freeze a **two-cell** candidate from the actual (M) input (including character and prime-place labels). Require the arithmetic target \(\mathcal H\), actual operators \(B_w\), and their cross HS moment. Independently construct two smooth Weil source functions and compute the full **mixed** mass and signed Gram. Attempt an explicit intertwiner that preserves both. If even a finite-cell candidate cannot preserve the normalized cross moment under the correct metrics, record NO TRANSFER for that projection before any infinite-dimensional speculation.

No promotion of Weil-Lab, provenance, Hodge or other RH custody states. This is a strict finite falsification gate.
