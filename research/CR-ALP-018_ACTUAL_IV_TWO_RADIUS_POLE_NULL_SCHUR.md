# CR-ALP-018 — Literal original IV010–IV011 two-support Schur certificate using two prime-indexed pole-null packets

**Canonical-Return, independent analytical source work — 2026-11-23 (America/Los_Angeles).**

**Standing:** **NEW exact-source two-function theorem, analytic proof with executed rational inequality audit PASS**. The constructed \(2\times2\) block is the **complete ORIGINAL signed zeta Weil form** on an explicitly specified nested-support smooth physical source pair. **NOT** the actual frozen 242×242 F41 numerical H / 242×4 V selected bank, not a Schur inequality on the entire old support core, not an old-to-new support theorem for arbitrary tests, not a new Lean kernel proof or RH.

## 0. Crucial live receiving-source recovery

Unlike CR-ALP-017, which had only a Weil-Lab source-review snapshot, this stage independently located the **actual live receiving repository**:

- **\`monocap-tech/weil-lab\`, branch \`research/iv068-original-source-kernel-gate\`** (not \`monocap-tech/weil\`). Source file \`WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportExtension.lean\`, contents blob SHA **\`5ef51fd37337c466620eafccee2a27b77955746e\`**: \`neutralIV010Mixed\` is literal \`Zeta23.EF.literatureRHS(neutralWindowCorrelation...)\`; \`neutralIV010OldPlusNew_realBlock\` and \`neutralIV010OriginalOneNewSource_nonnegative_of_crossBound\` encode the **complex Hermitian** one-old/one-new support extension if \(A,D\ge0\), \(|B|^2\le AD\).
- Same branch, \`WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean\`, blob SHA **\`14d8f01ce34f764900eccc68399d7323b4f4a25c\`**: \`neutralIV011MixedCorrelation_eq_schwartzIntegral\`, \`neutralIV011MixedCorrelation_vanish_far\`, \`neutralIV011OriginalMixedPrime_eq_finite\`, \`neutralIV011OriginalMixed_arithmeticSplit\` and \`neutralIV011OriginalOneNew_nonnegative_iff_crossBound\` retain the true **first-slot conjugation**, source support cutoff \(a_{\rm old}+a_{\rm new}\), and **pole − prime + arch**.
- The opened source-review file \`IV_CR_TRANSFER_II_ORIGINAL_RADIUS_OPERATOR.md\` explicitly fixes the original physical core \(C_c^\infty((-11/10,11/10);\mathbb C)\), source metric \(\int\overline u v\), and the full Gamma/contact symbol \(m_\Gamma(0)=-\gamma-3\log2-\pi/2-\log\pi\). It is an **analytic review document**, not independent certification of an operator closure.
- \`monocap-tech/weil-lab\` F41-001 source contract (PR #920) defines a **different** frozen 242-dimensional original physical selected H and four-column V with thousands of required entry checks. The two explicitly constructed functions below have **not** been identified with those selected columns; no frozen 242×242 H/V values were imported or fabricated.

These exact IV010/IV011 Lean declarations were **inspected**, not newly compiled or accepted by this CR phase. Their availability as source code does not certify the separate new inequalities below.

## 1. Two native physical sources, with genuinely nested support

Fix a single small step and mollifier width
\[
\boxed{b=10^{-6},\qquad\eta=10^{-12}.}
\tag{18.1}
\]
Let \(\rho_\eta\) be **any** real even nonnegative unit-mass \(C_c^\infty(-\eta,\eta)\) mollifier, and define
\[
g_b(u)=e^{u/2}\mathbf1_{[0,b]}(u)
-e^{u/2-b}\mathbf1_{[b,2b]}(u),\quad
g_\eta=\rho_\eta*g_b.
\]
For the two different primes \(p=2,3\), let
\[
L_p=2\log p,\qquad
h_p(u)=g_\eta(u+\log p+b)
-p\,g_\eta(u+\log p+b-L_p).
\tag{18.2}
\]
These are exactly the centered **compressed two-pole-null sources** of ALP-017 with \(X=p^2\). Each is a nonzero, real, compactly supported smooth Schwartz function; both satisfy
\[
\boxed{F_+(h_p)=F_-(h_p)=0.}
\tag{18.3}
\]
The two packets of \(h_p\) lie near \(\pm\log p\). Their smoothed packet width is \(w=2(b+\eta)\). Exact rational atanh/log enclosures certify
\[
\boxed{
\operatorname{supp}h_2\subset(-4/5,4/5),\quad
\operatorname{supp}h_3\subset(-11/10,11/10),\quad
\operatorname{supp}h_2\cap\operatorname{supp}h_3=\varnothing.
}
\tag{18.4}
\]
Moreover \(h_3\) is **not** supported inside \((-4/5,4/5)\). This is a genuine old-radius \(a=4/5\) to new-radius \(b_{\rm IV}=11/10\) single-new-physical-mode test.

Writing \(G=\|g_\eta\|_2^2\), the unchanged native physical metric is exactly
\[
\boxed{
M=\begin{pmatrix}5G&0\\0&10G\end{pmatrix},\qquad
G\ge b-2\eta>0.
}
\tag{18.5}
\]
The lower bound comes from the explicit unsmoothed positive interior \([\eta,b-\eta]\), on which \(g_\eta(u)=m_\eta e^{u/2}\), \(m_\eta=\int\rho_\eta(v)e^{-v/2}dv\ge1\).

## 2. Entire mixed prime term vanishes **on the authentic original IV011 source**

Because the source packets are centered at \(\pm\log2\) and \(\pm\log3\), the mixed correlation of \(h_2,h_3\) is supported only within \(w\) of
\[
\boxed{\pm\log(3/2),\qquad \pm\log6.}
\tag{18.6}
\]
The original IV011 **mixed** support cutoff is \(a+b_{\rm IV}=19/10\), so \(n\ge7\) makes **zero** contribution since \(\log7>19/10\). For \(n=2,3,4,5\), exact rational logarithm gaps give
\[
\log(3/2)+w<\log2,\quad \log5<\log6-w.
\]
The only integer at the second mixed peak is \(n=6\), but \(\Lambda(6)=0\). Therefore **with both signs and correct prime-power multiplicities**
\[
\boxed{
\mathcal W^{\rm mixed}_{\rm prime}(h_2,h_3)=0
\quad\text{EXACTLY}.
}
\tag{18.7}
\]
Both crossed two-pole terms vanish individually because each source is two-pole-null. The mixed full original signed Weil entry thus is **purely the complete archimedean Gamma interaction**; there is no omitted prime/correction term.

The two *diagonal* sources have a different arithmetic profile: the only nonzero original prime samples are \(n=4\) for \(h_2\) and \(n=9\) for \(h_3\), with \(k_{h_p}(\log p^2)=-pG\), hence
\[
\boxed{\mathcal W_{\rm prime}(h_p,h_p)=2\log p\,G>0.}
\tag{18.8}
\]
No pointwise claim about positivity on arbitrary old or new support tests follows.

## 3. Complete original Gamma **mixed** bound and genuine Schur determinant

Let \(k_{23}(\tau)=\int h_3(x)h_2(x-\tau)\,dx\), the actual IV011 first-slot-conjugated physical mixed correlation (the sources are real). Physical disjointness gives \(k_{23}(0)=0\). Under the original full Gamma/contact formula of ALP-002, ALP-017 and the opened IV source-review,
\[
\mathcal W_{\rm Gamma}(h_2,h_3)=
-\int_0^\infty[k_{23}(t)+k_{23}(-t)]
\frac{e^{-t/2}}{1-e^{-2t}}\,dt.
\tag{18.9}
\]
All points of its support have \(|t|>\log(3/2)-w>2/5\), and
\[
w_\Gamma(t)=e^{-t/2}/(1-e^{-2t})<2\qquad(t>2/5),
\]
using \(e^{4/5}>2\). Convolution \(L^1\) gives
\[
\big|\mathcal W_{\rm Gamma}(h_2,h_3)\big|
\le2\,\|h_2\|_1\|h_3\|_1.
\]
Since \(\|\rho_\eta*g_b\|_1\le\|g_b\|_1
=4\sinh(b/2)<3b\) and \(h_2,h_3\) consist of packets with weights \((1,2)\), \((1,3)\), respectively,
\[
\boxed{
|B|:=|\mathcal W_{\rm IV}(h_2,h_3)|<216b^2
=\frac{27}{125000000000}.
}
\tag{18.10}
\]
This is the **actual original crossed pole−prime+Gamma entry**, not a generic overlap proxy. In particular the mixed Weil entry is *small* even though the support transfer is genuine.

Each diagonal \(\mathcal W_{\rm IV}(h_p,h_p)\) has zero pole, strictly positive sole prime power, and the complete Gamma lower bound from ALP-017:
\[
\Gamma(k_{h_p})\ge
\left[m_\Gamma(0)+e^{-1/2}\log(1/w)\right]\|h_p\|_2^2.
\]
The exact rational inequalities \(m_\Gamma(0)>-6\), \(e^{-1/2}>3/5\), \(\log(1/w)>11\) give
\[
\boxed{A:=\mathcal W_{\rm IV}(h_2,h_2)>\frac35(5G),\qquad
D:=\mathcal W_{\rm IV}(h_3,h_3)>\frac35(10G).}
\tag{18.11}
\]
Combining (18.5), (18.10) and (18.11), the **actual original two-source Hermitian matrix**
\[
\boxed{H_{\rm IV}=\begin{pmatrix}A&B\\\overline B&D\end{pmatrix}}
\]
satisfies
\[
\boxed{
A>0,\ D>0,\quad AD>|B|^2,\qquad
H_{\rm IV}\succ0.
}
\tag{18.12}
\]
The exact audit certifies a stronger native-metric inequality:
\[
\boxed{
\mathcal W_{\rm IV}(c_2h_2+c_3h_3,\ c_2h_2+c_3h_3)
>\frac{5999}{10000}\,
\|c_2h_2+c_3h_3\|_2^2
\qquad
\forall(c_2,c_3)\ne(0,0)\in\mathbb C^2.
}
\tag{18.13}
\]
**Proof of constant:** \(G\ge b-2\eta\) and \(\sqrt{50}>7\) imply
\[
\frac{|B|}{\sqrt{(5G)(10G)}}<
\frac{216b^2}{7(b-2\eta)}<\frac1{10000}.
\]
For two complex coefficients, \(2|B||c_2c_3|\le
(|B|/\sqrt{50}G)(5G|c_2|^2+10G|c_3|^2)\). Add the two diagonal margins. This is the correct IV010 complex Hermitian Schur condition, not a comparison of the two diagonal signs alone.

The \(a=4/5\) old source \(h_2\) has positive full original signed energy; \(h_3\) is genuinely new at \(b_{\rm IV}=11/10\). The **live Lean IV010/IV011 conditional criterion's mathematical premises are explicitly met for this constructed finite old span and one new mode** at the original signed-form level. No Lean theorem with these particular \(h_p\) as instantiated \`SchwartzMap\` arguments was added or compiled.

## 4. Exact relationship to F41 H/V and the remaining problem

The original F41 H/V physical numerical inventory uses a **242×242** selected matrix and **242×4** mixed matrix with independent entry and residual obligations. The two functions above are **not** an identified pair of F41 selected basis vectors, and the \(H_{\rm IV}\) of (18.12) is **not** that numerical \(H\). An actual transfer would need a source-defined embedding into the frozen selected bank, preservation of the physical \(L^2\) metric and all original signed entries, a quantitative cross estimate against **every** old-domain direction, and the genuine infinite-dimensional Schur complement. None is provided by this two-source certificate.

The receiving fixed-radius source core also cannot become a universal support-extension theorem merely by repeating this positive old/new example. Outside the special prime-separated boundary packets, additional original prime terms and uncontrolled mixed Gamma contributions reappear.

**Gate result:**
| Task | Outcome |
|---|---|
| Live receiving IV010/IV011 source and exact original sign convention | **SOURCE-PINNED** |
| One strictly new physical source at radius 11/10 over radius 4/5 | **EXACT** |
| Both original Mellin pole moments | **EXACTLY ZERO** |
| Complete old/new cross-prime sum | **EXACTLY ZERO** |
| Complete crossed Gamma bound | **ANALYTIC with exact rational constants** |
| Old/new **full original signed Weil** Schur determinant | **STRICT PASS** |
| Native physical \(M\) and generalized positive margin \(5999/10000\) | **PASS on selected 2D span** |
| 242×242 original F41 selected H/V matrix identification | **NOT ESTABLISHED** |
| All-old-core Schur estimate / general collar support-extension | **OPEN** |
| New Lean build / full domain equivalence / RH | **NOT CLAIMED** |

**Next CR-ALP-019:** source-attach this *explicit genuine IV010 old/new two-source block* to a receiving formal Schwartz test with a standardized mollifier, and attempt an actual source-defined mixed representer bound against the pre-existing F41 selected physical modes. Do not conflate the 2D newly constructed span with F41's accepted 242-row bank. If the latter is inaccessible or incompatible, pursue a finite **prime-indexed packet-family** positive span as a separate local mathematical theorem, preserving the nonpromotion gate.

**Reproducibility:** \`src/audits/cr_alp_018_two_source_iv_schur_audit.py\` executed locally: **PASS** using only \`fractions.Fraction\`; no CAS floating-point inequalities. The continuum convolution/Gamma identities are the accompanying human analytic proof, not a Lean formalization.

**No writes to Weil-Lab / F41 / F42 / historical RH or provenance acceptance ledgers.**
