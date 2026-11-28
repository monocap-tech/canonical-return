# CR-ALP-023 — All 241 frozen F41 oscillator moments have certified strict signs; remaining original-IV parity Lean gate

**Date:** 2026-11-23, America/Los_Angeles.  
**Jurisdiction:** Canonical-Return analytical research. **Status:** SOURCE-INSPECTED + DERIVED universal packet-moment enclosure + independently EXECUTED **exact outward-rounded integer/Fraction audit PASS**. **Not a Lean-compiled theorem, not frozen F41 physical H/V reacceptance, not an infinite-dimensional signed Weil Schur bound, not RH.**

## 0. Exact source custody

This stage extends CR-ALP-021–022 from three nonzero physical rows to **all 241 original oscillator rows** without substituting a basis:

- Frozen original receiver **\`monocap-tech/weil-lab@a431f903254000bb5ad8edbc276ccc8599da98f5\`**:
  - \`WeilDefect/Arithmetic/ActualZetaOriginalJ242ExponentialSignature.lean\`, blob **\`b56b2dc6c482ca79e8f023d556bd52edf5e24f48\`**: EXACT 67 rational frequency groups and their contiguous monomial multiplicities \(\sum m_f=241\).
  - \`WeilDefect/Arithmetic/ActualZetaOriginalJ242SourceFunctions.lean\`, blob **\`53935773feb3bba3c10e699cd8f862f9dd590b76\`**: the literal generator is \(x^{k_j}\sin(2\pi f_j x)\) if \(k_j\) is even and \(x^{k_j}\cos(2\pi f_j x)\) if odd, with row 241 equal to \(\sinh(x/2)\).
  - \`WeilDefect/Arithmetic/ActualZetaOriginalNF407FrozenData.lean\`, blob **\`86ca70c8a2a8ba24ce4204081fef74d4ac543646\`**: independent flattened 241-node order matches those 67 groups.
  - \`WeilDefect/Arithmetic/ActualZetaOriginalPhysicalSolveResidual.lean\`, blob **\`2cf06688dc7e63bbb0c797102260559884325855\`**: the frozen historical \(\mathrm H_{242×242}\) and \(\mathrm V_{242×4}\) are **native physical \(L^2\)** moments of selected/critical generators; our packet moments are NEW row values, not entries of that old \(\mathrm V\).
- Real original IV receiver \`research/iv068-original-source-kernel-gate\`:
  - \`ActualZetaPhaseIVOriginalCrossSupportExtension.lean\` blob \`5ef51fd37337c466620eafccee2a27b77955746e\`, with literal \`neutralIV010Mixed\`;
  - \`ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean\` blob \`14d8f01ce34f764900eccc68399d7323b4f4a25c\`, genuine \`neutralIV011MixedCorrelation_eq_schwartzIntegral\` and pole − prime + arch split;
  - \`WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean\` blob \`0034976fe2f1e5fa9d411108736160542eb94293\`, literal \`Zeta23.EF.gammaBracket\` and \`literatureRHS\`;
  - \`WeilDefect/External/Zeta23/Theorems/Thm_Zeta23_WeilEF_digamma_conj.lean\` blob \`27f17b3b6e522f59d663486e851921b5f4c7dedd\`, already-proved digamma conjugation on \`integerComplement\`;
  - \`WeilDefect/External/Zeta23/Theorems/Thm_Zeta23_WeilEF_gammaR_bracket.lean\` blob \`9b43c6390b38b3377da7a8b01c65b86af4d6fba4\`, already-proved symmetric \(\Gamma_\mathbb R\) bracket.

The canonical 67 frequency groups were compared **independently against BOTH frozen Lean signatures** with the same FNV-1a64 transcription checksum **\`29ec35e3e29c3ee9\`**, and both yielded **241** oscillator nodes. The checksum is a transcription consistency test, not a cryptographic source attestation; immutable Git blob SHAs are above. **No live receiving Lean compilation was executed** because this container has neither \`lean\` nor \`lake\`.

## 1. Uniform exact physical moment theorem for the genuine frozen 241 rows

Keep the **same genuine two-pole-null raw p=2 source** used in CR-ALP-020–022:
\[
\begin{aligned}
&b=10^{-8},\quad \eta=10^{-14},\quad \ell=\log2,\\
&g_b(t)=e^{t/2}\mathbf1_{[0,b]}(t)-e^{t/2-b}\mathbf1_{[b,2b]}(t),\\
&\widetilde h_2(x)=g_b(x+\ell+b)-2g_b(x-\ell+b),\qquad h_2=\rho_\eta*\widetilde h_2,
\end{aligned}\tag{23.1}
\]
where **any** real, even, nonnegative, mass-one \(\rho_\eta\in C_c^\infty(-\eta,\eta)\) is allowed. It fits the true physical radius \((-11/10,11/10)\). Both original Mellin pole moments vanish exactly.

For original frozen oscillator index \(0\le j<241\), let \((f_j,k_j)\) be its actual integer-ordered NF545 signature, \(\omega_j=2\pi f_j\), and
\[
\varphi_j(x)=
\begin{cases}
x^{k_j}\sin(\omega_jx),&k_j\text{ even},\\
x^{k_j}\cos(\omega_jx),&k_j\text{ odd}.
\end{cases}
\]
Set the NEW physical mixed moment \(v_j=\int_{-11/10}^{11/10}h_2(x)\varphi_j(x)\,dx\). The exact universal leading coefficient is
\[
\boxed{A_j=\varphi_j'(\ell)-\tfrac32\varphi_j(\ell).}\tag{23.2}
\]
For \(r=3/4\), the genuine source support obeys \(|x|<r\). Put
\[
\begin{aligned}
M_{0,j}&=r^{k_j},\\
M_{1,j}&=k_jr^{k_j-1}+\omega_jr^{k_j},\\
M_{2,j}&=k_j(k_j-1)r^{k_j-2}+2k_j\omega_jr^{k_j-1}+\omega_j^2r^{k_j},
\end{aligned}
\]
where terms with a zero \(k_j\) prefactor are omitted. Then the **rigorous uniform enclosure** is
\[
\boxed{
\left|\frac{v_j}{b^2}-A_j\right|
\le\mathcal E_j:=
\frac{3b^2}{8}M_{0,j}
+\frac{b^2}{4}M_{1,j}
+2bM_{2,j}
+9\frac{\eta}{b}M_{1,j}.}\tag{23.3}
\]
Proof: in centered packet coordinates \(u=t-b\), the source's exact two half-packets satisfy \(g_b(b-u)=e^{(b-u)/2}\) and \(g_b(b+u)=-e^{-(b-u)/2}\) for \(0<u<b\), so
\[
v_j^{\rm raw}=\int_{-b}^{b}g_b(b+u)\bigl[-\varphi_j(\ell-u)-2\varphi_j(\ell+u)\bigr]du.
\]
Its two moments satisfy \(b^2/2\le m_0\le b^2/2+b^4/8\) and \(-b^2-b^4/4\le m_1\le-b^2\). Taylor's theorem bounds the quadratic remainder by \(2M_{2,j}b^3\) using \(\int u^2|g_b(b+u)|du\le4b^3/3\). The mollifier changes the pairing by at most \(\eta\|\widetilde h_2\|_1M_{1,j}<9b\eta M_{1,j}\). Combining gives (23.3), with **no quadrature** and no specialization of \(\rho_\eta\).

### 1a. FULL frozen signature: every oscillator entry strictly NONZERO

The new executable auditor evaluates (23.2) and (23.3) for **every** historical row \(0\le j<241\), retaining the actual order. It uses:
- \(\log2\) enclosed by a rational atanh series with rigorous tail;
- \(\pi\) enclosed by rational Machin alternating series;
- \(\theta_j=2\pi f_j\log2\) reduced modulo \(2\pi\) using an exact integer interval;
- \(24\)-term sine/cosine Taylor interval arithmetic with explicit next-term remainder and Lipschitz phase-radius error;
- **outward-rounded fixed-point integer arithmetic at \(10^{-52}\) resolution**, so no arbitrary float/round-to-nearest decisions enter the proof;
- rational outward bounds \(\mathcal E_j\) for **every** \((f_j,k_j)\).

**Independent local exact execution results:**
\[
\boxed{\begin{aligned}
\#\{j<241:v_j<0\}&=74,\\
\#\{j<241:v_j>0\}&=167,\\
\#\{j<241:\text{sign unresolved}\}&=0.
\end{aligned}}\tag{23.4}
\]
In fact the narrowest *certified* \(b^{-2}\)-normalized margin (row \(j=12\)) exceeds \(0.0773\), and the auditor verifies **strictly in exact integers**
\[
\boxed{\forall j<241,\qquad |v_j|>\frac7{100}b^2>0.}\tag{23.5}
\]
The entire conclusion is **uniform over every permitted mollifier**, including the canonical smooth bump selected in CR-ALP-022. The first three rows are again negative and satisfy all previously proved stronger thresholds.

**Independent cross-check:** 27 source-indexed rows, covering multiple frequency blocks through \(j=240\), were separately evaluated by \(65\)-digit high-precision direct quadrature of the **unsmoothed raw packet**. All 27 signs matched the exact interval certificates. The mpmath checks are diagnostic only; the formal sign proof is the outward-rational audit plus the analytic Taylor enclosure.

The selected \(j=241\) source is the **actual \(\sinh(x/2)\)** pole generator, and its new physical mixed moment is **exactly zero**, independently of the 241 signs:
\[
\boxed{v_{241}=\int h_2(x)\sinh(x/2)\,dx
=\tfrac12[F_+(h_2)-F_-(h_2)]=0.}\tag{23.6}
\]
Thus the new *242-row* physical mixed-vector pattern is **241 rigorously nonzero entries with assigned signs, and one exactly zero pole entry**. It is **NOT** the historical 242×4 physical V array, whose four columns are original NF587 critical sources.

The presence of these 241 nonzero pairings is consistent with CR-ALP-020's exact lower bound that at least **10% of the physical squared norm** of the raw packet lies OUTSIDE the frozen odd selected span, via its even component. Neither a complete span decomposition nor a full signed Weil operator interchange is asserted.

## 2. Actual original Weil even/odd formalization: precisely what is and is NOT closed

The **ANALYTIC** source-core result remains true:
\[
\boxed{\mathrm{neutralIV010Mixed}(B,f_{\rm even},o_{\rm odd})=0}
\]
for compact smooth supported Schwartz sources. The literal \`neutralWindowCorrelation\` is odd under \(\tau\mapsto-\tau\); the full original \`Zeta23.EF.literatureRHS\` annihilates odd test kernels. It is useful to factor the proposed receiving Lean proof into **two separate reusable statements** rather than attempting a monolithic Schur rewrite:

**EF lemma:** under sufficient compactness / smoothness / integrability,
\[
(\forall t,\ k(-t)=-k(t))\Longrightarrow
\mathrm{Zeta23.EF.literatureRHS}(k)=0.
\]
The two pole Fourier evaluations cancel, each von Mangoldt summand \(k(\log n)+k(-\log n)\) is zero, and the entire Gamma integral vanishes since \(\mathrm{gammaBracket}(-r)=\mathrm{gammaBracket}(r)\), derivable from the **already-source-proved** \`Zeta23.WeilEF.digamma_conj\` or symmetric \`gammaR_bracket\).

**Correlation lemma:** under IV011's **actual support hypotheses**, \(\mathrm{neutralWindowCorrelation}(B,f_{\rm even}.toLp,g_{\rm odd}.toLp)\) is odd, using \`neutralIV011MixedCorrelation_eq_schwartzIntegral\`, reflection substitution, and the literal conjugation order.

Compose these two to conclude \(\mathrm{neutralIV010Mixed}=0\); the reversed entry then follows by \`neutralIV010Mixed_hermitian\` already in the receiving source.

**Not completed in CR-ALP-023:** no Lean binary exists in the local container, no Lean file was compiled, no \`#print axioms\` report was produced, and the receiving branch was **not modified**. A human proof and source-level imported lemma do not equal a kernel-checked theorem. The outstanding Lean obligations are reflection of the true L² window representative, odd \`Zeta23.paperFT\` for complex arguments and its Gamma-integrability, and type-specific \`gammaBracket_neg\`.

## 3. Disposition / exact reproducibility

| Gate | CR-ALP-023 disposition |
|---|---|
| 67 native frequencies / 241 exact original oscillator rows | **DOUBLE-SOURCE PINNED; FNV checksum MATCH** |
| Uniform derived moment error formula for every frozen row | **ANALYTIC PASS** |
| Source-indexed strict nonzero signs for all 241 rows | **EXACT INTERVAL PASS (74 negative / 167 positive)** |
| Global minimum absolute physical moment >\(7b^2/100\) | **EXACT INTERVAL PASS** |
| Last actual selected sinh-pole physical row | **EXACT ZERO** |
| 27 independent high-precision direct raw integral checks | **DIAGNOSTIC PASS** |
| Compact smooth source and exact two original Mellin pole cancellations | **RETAINED from CR-ALP-020–022** |
| Literal EF odd-kernel and IV010 mixed parity | **HUMAN ANALYTIC PASS / LEAN OPEN** |
| Historical frozen F41 242×242 H and 242×4 V entry acceptance | **UNTOUCHED** |
| Entire old-collar Schur / all-even / RH | **OPEN / NOT CLAIMED** |

**Executable audit:** \`src/audits/cr_alp_023_all_f41_rows_exact_audit.py\`, local executable, no third-party dependencies, **PASS**. The locally run version SHA256 is \`098b686a34398cc983f79f8a5cbc4365bc7a07fbdd17043d131dccbfd1285106\`; the GitHub copy is from a text-staging transfer and contains the same source checks (line-ending normalization may change the byte SHA). The GitHub source is committed on the Canonical-Return CR-ALP branch. No unseen CI runs claimed.

**Next — CR-ALP-024:** obtain an actual Lean runtime and discharge the EF odd-kernel theorem, gammaBracket parity and native IV011 correlation parity with an axiom/kernel check; then study mixed Weil cross-term bounds against genuinely **arbitrary odd source-core directions**, not the frozen physical moment vector alone. This is the remaining phase-specific Schur hurdle. Keep PR draft/nonpromoting, and do not promote any claim to Weil-Lab/F41/F42/RH/provenance acceptance.
