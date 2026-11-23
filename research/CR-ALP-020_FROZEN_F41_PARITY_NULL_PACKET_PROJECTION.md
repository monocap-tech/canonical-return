# CR-ALP-020 — Frozen F41 parity-null physical columns and even original-Weil packet Schur certificate

**Canonical-Return research-only, 2026-11-23.** SOURCE-PINNED + DERIVED EXACT PARITY + EXACT LOCAL RATIONAL AUDIT PASS. **Not** a new Lean compile, F41 numerical H/V acceptance, all-old-core Schur, completed Connes/Suzuki domain theorem, global Weil positivity or RH.

## 1. Frozen source custody: actual F41 functions are all odd

Directly inspected **monocap-tech/weil-lab@a431f903254000bb5ad8edbc276ccc8599da98f5**:
- \`WeilDefect/Arithmetic/ActualZetaOriginalJ242SourceFunctions.lean\`, blob **53935773feb3bba3c10e699cd8f862f9dd590b76**: frozen \`nf545OriginalSelectedRealGenerator j\` for j<241 is \(x^{k_j}\sin(2\pi\omega_jx)\) if \(k_j\) even, and \(x^{k_j}\cos(2\pi\omega_jx)\) if \(k_j\) odd, with original NF407 rational frequency/order; j=241 is \(\sinh(x/2)\).
- \`WeilDefect/Arithmetic/ActualZetaOriginalNF407SourceBinding.lean\`, blob **342c6cffee6cee2719dfe35572a9e05aa8979726**: each of the four \`neutralNF587CriticalGenerator\` functions is exactly **32 original sine terms plus the original \(\sinh(x/2)\)**, with frozen rational coefficients.
- \`WeilDefect/Arithmetic/ActualZetaOriginalNF407NormalizedComplexSource.lean\`, blob **301acf2f7ef626f094a0f9a48a135048f27a3402**: the original NF588 complexified source explicitly reduces to these same real generators.
- \`WeilDefect/Arithmetic/ActualZetaOriginalPhysicalSolveResidual.lean\`, blob **2cf06688dc7e63bbb0c797102260559884325855**: \`neutralNF606PhysicalH\` is the physical \(L²(-11/10,11/10)\) **242×242 Gram**, while \`neutralNF606PhysicalV\` is the physical **242×4 mixed matrix** of the above source functions. **Neither is the signed Weil matrix.**

Every one of these 242 selected plus four critical functions is **ODD**:
\[
\boxed{\varphi_j(-x)=-\varphi_j(x)\ (j<242),\qquad
\chi_k(-x)=-\chi_k(x)\ (k<4).}\tag{20.1}
\]
This is a literal parity fact for the frozen native source, independent of any H/V numerical enclosures.

For the **actual signed Weil** cross-form convention, also directly reopened \`monocap-tech/weil-lab\` branch \`research/iv068-original-source-kernel-gate\`: \`ActualZetaPhaseIVOriginalCrossSupportExtension.lean\` (blob **5ef51fd37337c466620eafccee2a27b77955746e**) and \`ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean\` (blob **14d8f01ce34f764900eccc68399d7323b4f4a25c**). Its \`neutralIV010Mixed\` is the literal original zeta correlation form, with original crossed poles, negative prime sum, complete Gamma, first-slot conjugation, and complex Schur geometry. This note DOES NOT claim a new Lean compile of those declarations.

## 2. A type-correct exact source-moment formula (not a numerical F41 substitution)

Retain CR-ALP-019's parameters \(b=10^{-8}\), \(\eta=10^{-14}\), and ANY real even nonnegative unit-mass \(C_c^\infty(-\eta,\eta)\) mollifier \(\rho_\eta\). Let
\[
g_b(u)=e^{u/2}1_{[0,b]}(u)-e^{u/2-b}1_{[b,2b]}(u),\quad
g=\rho_\eta*g_b,\quad G=\|g\|_2^2\ge b-2\eta>0.
\]
The centered prime packet is
\[
\boxed{h_p(x)=g(x+\ell_p+b)-p\,g(x-\ell_p+b),\quad\ell_p=\log p.}\tag{20.2}
\]
For p=2,3 its support lies strictly within the F41 interval \((-11/10,11/10)\), and both Mellin pole moments \(F_\pm(h_p)=0\) exactly.

For ANY frozen selected generator \(\varphi_j\), the new native PHYSICAL mixed moment (not part of the historical V) has the exact identity
\[
\boxed{v_{j,p}=\int_{-11/10}^{11/10}\varphi_j(x)h_p(x)dx
=\int g(t)[\varphi_j(t-\ell_p-b)-p\varphi_j(t+\ell_p-b)]dt.}\tag{20.3}
\]
A closed-form frequency evaluator follows. With \(\widehat f(\xi)=\int f(x)e^{i\xi x}dx\), \(\alpha=1/2+i\xi\),
\[
\boxed{\begin{aligned}
\widehat g_b(\xi)&=\frac{(e^{\alpha b}-1)(1-e^{-b}e^{\alpha b})}{\alpha},\\
\widehat h_p(\xi)&=\widehat\rho_\eta(\xi)\widehat g_b(\xi)
e^{-i\xi(\ell_p+b)}(1-p e^{2i\xi\ell_p}).
\end{aligned}}\tag{20.4}
\]
The \(\alpha=0\) singularity is removable. For an oscillator source \(x^{k_j}\sin(\omega_jx)\), \(k_j\) even, the moment is the imaginary part of \(i^{-k_j}\partial_\xi^{k_j}\widehat h_p(\omega_j)\); for odd \(k_j\), \(x^{k_j}\cos(\omega_jx)\), it is the real part. The selected j=241 pole row vanishes identically by \(F_\pm(h_p)=0\). Each critical generator's moment is the original 32-term weighted sine-frequency sum; its \(\sinh\) part vanishes.

**This is a source-exact recipe, not 242 numeric interval certificates.** A fixed explicit \(\rho_\eta\) and original frequency/critical arrays must be attached before any such entry audit. No original F41 H or V has been overwritten.

## 3. Exact parity-null attachment and source-span nonmembership

Let reflection be \((\mathsf Rh)(x)=h(-x)\), and define
\[
\boxed{E_p=h_p+\mathsf Rh_p,\quad p=2,3.}\tag{20.5}
\]
These are nonzero **EVEN** real \(C_c^\infty\) sources still strictly within the frozen F41 \(11/10\) radius. The earlier \(4/5\) radius contains \(E_2\) but not \(E_3\). Both original Mellin moments vanish exactly because reflection interchanges \(F_+\) and \(F_-\).

By the exact frozen odd-source formulas (20.1), the **entire original physical generator inventory has zero mixed moments**:
\[
\boxed{\int_{-11/10}^{11/10}E_p(x)\varphi_j(x)\,dx=0\ \forall j<242,
\qquad
\int_{-11/10}^{11/10}E_p(x)\chi_k(x)\,dx=0\ \forall k<4.}\tag{20.6}
\]
This is an exact all-246 source-level null **new-column** certificate, NOT an assertion that the historical odd/odd F41 H/V entries vanish.

The original \(h_p\) themselves cannot lie in the frozen selected physical span \(S_{242}\). On the two disjoint packet sectors, \(E_p\) is identified with \(g-p\mathsf R_bg\) and its reflected copy, where \(\mathsf R_bg(u)=g(2b-u)\). By the reverse triangle inequality,
\[
\|E_p\|_2^2=2\|g-p\mathsf R_bg\|_2^2\ge2(p-1)^2G.
\]
Since \(P_{\rm even}h_p=E_p/2\) and \(S_{242}\subset L^2_{\rm odd}\),
\[
\boxed{\frac{\operatorname{dist}(h_p,S_{242})^2}{\|h_p\|_2^2}
\ge\frac{(p-1)^2}{2(1+p^2)}
=\begin{cases}1/10,&p=2,\\1/5,&p=3.\end{cases}}\tag{20.7}
\]
At least 10% or 20% of the squared physical mass of these raw packet sources is outside the entire frozen odd F41 selected span, regardless of numerical basis conditioning. This does NOT allege F41 attempted to span all even physical sources.

For **smooth compactly supported** even \(E\) and odd \(o\), the authentic signed Weil cross correlation \(C(E,o;\tau)=\int o(x)\overline{E(x-\tau)}dx\) is odd in \(\tau\); the complete original pole/prime/Gamma distribution is even. Therefore \(W_{\rm original}(E,o)=0\) **on the common smooth source core**. This is NOT a domain or Weil-mixed assertion for sharply truncated analytic F41 selected functions, which have only been attached in native L².

## 4. A positive even two-prime ORIGINAL signed Weil block, in the fixed F41 radius

The even sources \(E_2,E_3\) have disjoint supports, so their native physical metric is
\[
M_{\rm even}=\operatorname{diag}(m_2,m_3),\qquad
m_p=\|E_p\|_2^2\ge2(p-1)^2G.
\]
In particular \(m_2\ge2G\), \(m_3\ge8G\), and \(\sqrt{m_2m_3}\ge4G\).

Each autocorrelation has support only within \(w=2(b+\eta)\) of \(0,\pm2\log p\). Hence the **only** possible genuine original diagonal prime-power sample is \(n=p^2\). Its sign after symmetrization is not presumed. The absolute original prime coefficient is \(2\log p/p<1\) for p=2,3, and Cauchy gives \(|k_p(\log p^2)|\le m_p\). Therefore the **complete signed prime part is bounded below by \(-m_p\)**.

The full original Gamma/contact distribution remains
\[
\Gamma(k)=m_\Gamma(0)k(0)+
2\int_0^\infty[k(0)-k(t)]\,
\frac{e^{-t/2}}{1-e^{-2t}}dt.
\]
The integrand is nonnegative for an autocorrelation, and \(k_p(t)=0\) throughout \([w,1]\). Exact rational bounds \(m_\Gamma(0)>-6\), \(e^{-1/2}>3/5\), \(\log(1/w)>16\) imply
\[
W_\Gamma(E_p,E_p)>\frac{18}{5}m_p,\qquad
\boxed{W_{\rm original}(E_p,E_p)>\frac{13}{5}m_p.}\tag{20.8}
\]

Mixed \(E_2,E_3\) correlations live only near \(\pm\log(3/2)\) and \(\pm\log6\). Both crossed pole terms vanish; all original prime-power samples vanish (6 is not a prime power). Thus the mixed entry is the **complete original Gamma** alone, and with \(\inf|t|>1/4\), \(w_\Gamma(t)<3\), \(\|E_p\|_1<6(1+p)b\):
\[
\boxed{|W_{\rm original}(E_2,E_3)|<1296b^2.}\tag{20.9}
\]
Its M-normalized magnitude is \(<1296b^2/[4(b-2\eta)]<10^{-4}\), as certified in the exact rational auditor. Applying the genuine complex-Hermitian Schur bound gives
\[
\boxed{
H_{\rm even}:=(W_{\rm original}(E_p,E_q))_{p,q=2,3}
>\frac{25999}{10000}M_{\rm even}.}\tag{20.10}
\]
Equivalently, the full original signed Weil energy of **every nonzero complex combination** \(z_2E_2+z_3E_3\) exceeds \((25999/10000)\) times its actual physical \(L²\) mass. No surrogate PSD kernel has been introduced.

## 5. Disposition and next

| Gate | Standing |
|---|---|
| Actual frozen 242+4 source generator formulas | **SOURCE-PINNED** |
| All frozen functions odd | **EXACT / universal over 246 indices** |
| Exact all-row formula for unsymmetrized packet mixed moments | **ANALYTIC, no numerical enclosures** |
| New even packet physical mixed columns against ALL frozen modes | **EXACT ZERO** |
| Quantitative nonmembership of raw h_2,h_3 in frozen F41 span | **10%, 20% squared-mass floors** |
| Even/odd original signed cross pairing on *smooth* core | **ANALYTIC EXACT** |
| Even two-prime genuine signed Weil Schur margin | **\(25999/10000\), exact-rational numerical gates** |
| Original 242×242 H/242×4 V matrix numerical acceptance | **UNTOUCHED** |
| Entire old-core/infinite-even-complement Weil positivity | **OPEN** |
| New Lean kernel proof / RH | **NOT CLAIMED** |

**Audit:** \`src/audits/cr_alp_020_f41_parity_packet_audit.py\`, local Fraction regression **EXECUTED PASS**, proves exact rational radius/log gaps, mass fractions, Gamma/prime/Schur numerical inequalities; the analytic parity, convolution and Gamma identities are the separate human proof, **not** Lean compiled.

**CR-ALP-021 candidate:** formalize literal Schwartz reflection and even/odd source cross-block vanishing for \`neutralIV010Mixed\` in the receiving original Lean source, preserving original correlation and pole/prime/Gamma conventions. Separately fix a canonical \(\rho_\eta\) and certify **at least one nonzero** F41 oscillator physical mixed row of the raw packet \(h_p\) by (20.4) with original frequency/order, instead of inventing a 242-vector. No writes to Weil-Lab, F41/F42, historical RH/provenance or Hodge acceptance ledgers. Keep CR draft/nonpromoting.
