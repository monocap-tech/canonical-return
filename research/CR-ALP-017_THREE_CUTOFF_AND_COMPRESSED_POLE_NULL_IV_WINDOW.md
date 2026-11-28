# CR-ALP-017 — Exact three-cutoff arithmetic identity and a compressed pole-null test inside Phase IV's original radius

**Jurisdiction:** Canonical-Return / independent analytical return after the CR-ALP-015 operator-phase STOP.  
**Date:** 2026-11-23 (PST).  
**Standing:** new elementary ANALYTIC source-core lemmas; independent self-contained algebra/rational audit **EXECUTED LOCALLY: PASS**. **No Lean kernel certification, original Phase IV Schur/complement inequality, support extension to arbitrary radii, Connes/Suzuki completed-domain admission, or independent RH proof.**

## 0. Source custody and new distinction

Reopened, not guessed:
- CR repo, draft PR #3: \`research/CR-ALP-002_FULL_WEIL_RAMP_AND_POLE_NO_GO.md\`, \`CR-ALP-003_SMOOTH_CELL_GRAM_AND_LEAN_INTERMEDIATE.md\`, \`CR-ALP-004_QUANTITATIVE_SMOOTH_WEIL_TRANSFER.md\`, \`CR-ALP-007_TRIANGLE_GAUGE_HOLONOMY_AND_ARITHMETIC_PHASE.md\`, and \`CR-ALP-016_TWO_POLE_NULL_SMOOTH_WEIL_AND_PRIME_WINDOW.md\` plus its exact auditor.
- Weil-Lab **source-review snapshot** \`IV_CR_TRANSFER_II_ORIGINAL_RADIUS_OPERATOR.md\`, dated 2026-11-23, describes the actual proposed physical source core \(D=C_c^\infty((-11/10,11/10);\mathbb C)\), antilinear-first-slot metric, \`WeilDefect/Arithmetic/ActualZetaPhaseIVOriginalCrossSupportArithmeticCutoff.lean\` with seven prime powers \(2,3,4,5,7,8,9\), and the complete crossed two-pole/full-prime/full-Gamma source dictionary. That snapshot is a human analytical **candidate** and not a finished Lean/IV certificate. The cited private repo ref \`165981950e193a8b6ea219f7b31e379db495befa\` could not be independently fetched via the currently accessible GitHub connector; statements about this specific IV core are based on the reopened exact review, not represented as a newly verified live branch.
- Classical additional external input for §2's optional eventual positivity: the *known* prime number theorem plus \(\sum_{n\le y}\Lambda(n)/n=\log y-\gamma+o(1)\). This is an **external classical theorem**, not derived in CR or from Alpöge's arithmetic hypothesis (M).

The central result is a **two-track repair**: (A) the long-shift source of ALP-016 has an exact three-cutoff arithmetic identity, avoiding a need to invent a fourth prime statistic; (B) the same two-moment annihilator with a **short translation step** gives a *different* pole-null smooth source that actually fits the fixed original IV physical radius. Neither track alone establishes support-independent positivity.

## 1. Exact long-shift three-cutoff identity (all \(X>1\))

Let \(f_X(u)=e^{u/2}\mathbf1_{[0,L]}(u)\), \(L=\log X\), \(q=\sqrt X\). For the disjoint three-shift source from ALP-016,
\[
h_X=f_X-(q+q^{-1})\,f_X(\cdot-L)+f_X(\cdot-2L),
\quad k_{h_X}=h_X*\widetilde h_X.
\]
Put
\[
\boxed{\alpha_X=3X+9+\frac3X,\quad
\beta_X=3+\frac3X,\quad \chi_X=\frac1X.}
\tag{17.1}
\]
The ramp \(k_Y(t)=(Y e^{-|t|/2}-e^{|t|/2})_+\) satisfies, **as an even function on all real \(t\)**,
\[
\boxed{k_{h_X}=\alpha_X k_X-\beta_X k_{X^2}+\chi_X k_{X^3}.}
\tag{17.2}
\]
**Proof:** for \(n=e^t\ge1\), write \(P_Y(n)=(Y/n-1)_+\), so \(k_Y(t)/\sqrt n=P_Y(n)\). The shift-one term of the correlation obeys
\[
\frac{k_X(t-L)}{\sqrt n}
=\frac1qP_{X^2}(n)-(q+q^{-1})P_X(n),
\]
and the shift-two term obeys
\[
\frac{k_X(t-2L)}{\sqrt n}
=P_X(n)-\Big(1+\frac1X\Big)P_{X^2}(n)+\frac1X P_{X^3}(n).
\]
For \(t\ge0\), the positive shifts \(k_X(t+L),k_X(t+2L)\) vanish. Substituting the full correlation
\((2+(q+q^{-1})^2)k_X(t)-2(q+q^{-1})k_X(t-L)+k_X(t-2L)\)
and collecting coefficients gives (17.2), extended by evenness. The algebra is checked piecewise on \([1,X]\), \([X,X^2]\), \([X^2,X^3]\), beyond \(X^3\), and in the exact rational sample loop \(X=9,n=1,\dots,800\).

The source's two-pole cancellation is reflected **in the kernel identity**:
\[
\boxed{\alpha_X(X-1)-2\beta_X(X^2-1)+3\chi_X(X^3-1)=0.}
\tag{17.3}
\]
Consequently **the original signed Weil distribution, including all pole, prime-power and full Gamma terms**, satisfies by linearity
\[
\boxed{\mathcal W(k_{h_X})=
\alpha_X\mathcal W(k_X)-\beta_X\mathcal W(k_{X^2})
+\chi_X\mathcal W(k_{X^3}).}
\tag{17.4}
\]
This is an exact signed distributional identity, **not** an assertion that the three RHS kernels share a positive target Hilbert metric or a completed operator.

### 1a. The full prime-power window collapses to three ordinary \(S\)-values

Define \(S(Y)=\sum_{n\le Y}\Lambda(n)(Y/n-1)\), and \(G(Y)=S(Y)-Y\log Y+(1+\gamma)Y\), exactly as CR-ALP-002–004. Equation (17.2) gives
\[
\boxed{
\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,
k_{h_X}(\log n)
=\alpha_X S(X)-\beta_X S(X^2)+\chi_X S(X^3).
}
\tag{17.5}
\]
This **corrects the impression** that wholly new *types* of prime sums are required: the mixed windows are exactly a signed linear combination of the **same** ramp statistic at **three different scales**. ALP-016's statement that a sole bound at \(X\) cannot control the new window remains correct.

Substitute the exact original ramp identity
\[
\mathcal W(k_Y)=C_\zeta(Y+1)-2G(Y)+2\log(2\pi)-2D(Y),
\quad
D(Y)=\sum_{m\ge1}\frac{Y^{-2m}}{(2m)(2m+1)}.
\]
Writing \(\mathcal J(X)=\alpha_XG(X)-\beta_XG(X^2)+\chi_XG(X^3)\), the **complete** result is
\[
\boxed{
\begin{aligned}
\mathcal W(k_{h_X})={}&
C_\zeta\left(X^2+9X+9+X^{-1}\right)
+2\log(2\pi)\left(3X+6+X^{-1}\right)\\
&-2\mathcal J(X)
-2\left[\alpha_XD(X)-\beta_XD(X^2)+\chi_XD(X^3)\right].
\end{aligned}}
\tag{17.6}
\]

**Noncircular sufficient signed arithmetic gate:** an upper bound on \(\mathcal J(X)\) with coefficient below \(C_\zeta X^2/2\) suffices for eventual positivity. Separate *upper* bounds on \(G(X),G(X^2),G(X^3)\) alone do **not** upper-bound \(\mathcal J\), because the \(G(X^2)\) coefficient is **negative**: one additionally needs a **lower** bound at \(X^2\), or a direct bound on the signed combination. No term may be silently discarded.

**Classical unconditional comparison (not Alpöge input):** the known PNT/Mertens asymptotics give \(G(Y)=o(Y)\), hence
\[
\boxed{\mathcal J(X)=o(X^2),\qquad
\mathcal W(k_{h_X})=C_\zeta X^2+o(X^2)>0
\text{ for all sufficiently large }X.}
\tag{17.7}
\]
This is **only a single explicitly parameterized positive source family**, not Weil's universal criterion. In particular it does **not** prove a new prime-number theorem or RH; it uses the already known PNT. (Classical reference: https://mathoverflow.net/questions/437114/proving-mertens-theorem-using-the-prime-number-theorem and standard PNT \(\psi(Y)\sim Y\).)

For an even compact nonnegative unit-mass mollifier of width \(\eta_X=X^{-4}\), the convolution preserves both Mellin pole zeros *exactly*, and (17.2) holds also after autocorrelation smoothing. Applying CR-ALP-004's **uniform single-ramp** complete-form error estimate separately at \(X,X^2,X^3\), with coefficients \(\alpha_X\asymp X,\beta_X\asymp1,\chi_X=X^{-1}\), gives
\[
\boxed{
|\mathcal W(k_{h_X,\eta_X})-\mathcal W(k_{h_X})|
=O(X^{-1/2}\log(eX)),\quad X\ge4.
}
\tag{17.8}
\]
The largest contribution is \(X^{-1}\cdot (X^3)^{3/2}\eta_X\log(X^3)\). Thus the *long-shift* smooth family also has (17.7). Its centered support radius still grows like \(\frac32\log X+\eta_X\), and **cannot** inhabit any fixed support collar as \(X\to\infty\).

## 2. A distinct *short-shift* pole-null source that fits the original Phase IV radius

The disjoint-translate choice \(b=L\) in ALP-016 was convenient for an exact \(3\times3\) Gram, **but is not necessary for two-pole cancellation**. Let \(0<b<L/2\) be any **short** step, define
\[
s_b=e^{b/2}+e^{-b/2},\qquad
\boxed{h_{X,b}=f_X-s_b f_X(\cdot-b)+f_X(\cdot-2b).}
\tag{17.9}
\]
The translation polynomial \((1-e^{b/2}T_b)(1-e^{-b/2}T_b)\) kills both \(F_\pm\), since \(T_b\) multiplies the moments by \(e^{\pm b/2}\). Thus
\[
\boxed{F_+(h_{X,b})=F_-(h_{X,b})=0.}
\tag{17.10}
\]
The source has support **hull** \([0,L+2b]\), rather than \([0,3L]\), so its minimal centered support radius is
\[
\boxed{a_{\rm needed}=\frac L2+b.}
\tag{17.11}
\]

There is an even stronger exact boundary-packet factorization. Set
\[
g_b(u)=e^{u/2}\mathbf1_{[0,b]}(u)
-e^{u/2-b}\mathbf1_{[b,2b]}(u).
\]
Because \(\mathrm e^{u/2}\) is annihilated on every overlap interior, the three overlapping ramps cancel there, leaving
\[
\boxed{h_{X,b}=g_b-\sqrt X\,g_b(\cdot-L).}
\tag{17.12}
\]
One checks \(F_+(g_b)=0\) and \(F_-(g_b)=b(1-e^{-b})\). Thus (17.12) separately kills the remaining minus moment. The two narrow packets have disjoint supports for \(2b<L\). Also \(\|g_b\|_2^2=2\sinh b\), so \(\|h_{X,b}\|_2^2=2(X+1)\sinh b>0\).

For **any** real even nonnegative unit-mass mollifier \(\rho_\eta\in C_c^\infty(-\eta,\eta)\), put \(g_\eta=\rho_\eta*g_b\) and \(h_\eta=g_\eta-\sqrt X\,g_\eta(\cdot-L)\). It remains nonzero \(C_c^\infty\), with \(F_\pm(h_\eta)=0\) **exactly** and
\[
\|h_\eta\|_2^2=(X+1)\|g_\eta\|_2^2,
\qquad
\mathrm{supp}(h_\eta)\subset[-\eta,L+2b+\eta].
\tag{17.13}
\]
No delicate numerical smoothing-continuity argument is needed for the following sign certificate.

### 2a. Explicit radius \(a_{\rm IV}=11/10\), \(X=9\)

Using the *reopened* original-IV-radius source-review snapshot, set
\[
\boxed{
X=9,\quad L=\log9,\quad
b=10^{-6},\quad \eta=10^{-12},
\quad w=2(b+\eta).
}
\]
Translation by \(-(L/2+b)\) centers the smooth test, with
\[
\boxed{
a_{\rm needed}=\tfrac12\log9+b+\eta<11/10=a_{\rm IV}.
}
\tag{17.14}
\]
The strict radius inequality is independently checked by rational atanh-series log enclosures. In contrast, the original **disjoint** three-shift test at \(X=9\) requires \(\frac32\log9>11/10\) and lies **outside** this fixed IV core. This is a **different test**, not a shortening of the same function by translation.

The smoothed packet \(g_\eta\) has support width at most \(w\), and the two packets are separated by \(L\). Their even autocorrelation \(k=h_\eta*\widetilde h_\eta\) is supported in
\[
[-w,w]\cup[L-w,L+w]\cup[-L-w,-L+w].
\]
All original prime-power positions \(\log n\) satisfy
\[
w<\log2,\quad \log8<L-w,\quad L+w<\log10.
\]
Therefore the **only** possible original prime contribution is at \(n=9=3^2\), with the exact sign
\[
\boxed{
k(\log9)=-3\|g_\eta\|_2^2<0,\quad
\mathcal W_{\rm prime}(k)=
-2\frac{\log3}{3}k(\log9)
=2\log3\|g_\eta\|_2^2>0.
}
\tag{17.15}
\]
Every other prime-power term vanishes **exactly** (including 2,3,4,5,7,8 and all \(n\ge10\)). The crossed original two-pole contribution vanishes exactly because \(F_\pm(h_\eta)=0\).

### 2b. Direct analytic **full-Gamma positive lower bound**

The original signed full Gamma expression of CR-ALP-002 can be rewritten (no change of conventions) as
\[
\boxed{
\Gamma(k)=m_\Gamma(0)k(0)
+2\int_0^\infty[k(0)-k(t)]w_\Gamma(t)\,dt,
}
\quad
w_\Gamma(t)=\frac{e^{-t/2}}{1-e^{-2t}},
\tag{17.16}
\]
where
\[
m_\Gamma(0)=-\gamma-3\log2-\log\pi-\frac\pi2>-6.
\]
The identity follows from
\[
2\int_0^\infty(1-e^{-t/2})w_\Gamma(t)\,dt
=4\int_0^1\frac{du}{(1+u)(1+u^2)}
=\log2+\frac\pi2.
\]
It agrees with the original Gamma symbol in the **source-review candidate** \`IV_CR_TRANSFER_II_ORIGINAL_RADIUS_OPERATOR.md\`, without replacing the pole term or making a spectral isometry claim.

Because \(k\) is the real autocorrelation of a nonzero \(L^2\) source,
\(k(t)\le k(0)=\|h_\eta\|_2^2\) for every \(t\), by Cauchy–Schwarz. Therefore the integrand in (17.16) is **nonnegative** for every \(t>0\). Moreover \(k(t)=0\) on the entire interval \([w,1]\), since \(w\ll1<L-w\). For \(0<t\le1\),
\[
w_\Gamma(t)\ge \frac{e^{-1/2}}{2t}.
\]
All following inequalities have **independent exact rational certificates** in the attached audit:
\[
e^{-1/2}>\frac35,\qquad
\log(1/w)>11,\qquad
m_\Gamma(0)>-6.
\]
Hence
\[
\boxed{
\Gamma(k)\ge
\left[m_\Gamma(0)+e^{-1/2}\log(1/w)\right]k(0)
>\frac35\,k(0).
}
\tag{17.17}
\]
Together with zero poles and the **strictly positive** sole prime-power term (17.15), this proves the genuinely smooth original signed-form inequality
\[
\boxed{
\mathcal W(k_{h_\eta})>
\frac35\,\|h_\eta\|_2^2>0,\qquad
h_\eta\in C_c^\infty((-11/10,11/10))
}
\tag{17.18}
\]
after centering, **for every** such mollifier \(\rho_\eta\) of the stated width. The bound is *direct*—it does not depend on the possibly crude smoothing-error budget of ALP-016. The exact auditor certifies the rational inequalities and finite prime cutoffs; the continuum Gamma identity, Cauchy and support argument are given here as a human analytic proof, **not independently Lean-compiled**.

**Scope:** The constructed test belongs to the specified *smooth original IV source core* as described in the reopened Phase IV review, and satisfies both algebraic pole-null conditions. This **does not** certify any additional closed/semilocal Connes/Suzuki domain, or an actual selected old/new matrix \`H/V\` datum, Schur complement bound, neutral-shell persistence or global RH.

## 3. What the Phase IV support gate now actually says

At fixed radius \(a_{\rm IV}=11/10\), the family \(f_X\) or \(h_{X,b}\) can fit only if
\[
\boxed{\frac12\log X+b+\eta<a_{\rm IV}.}
\]
For the original **long** three-shift source, the condition is
\(\frac32\log X+\eta<a_{\rm IV}\). Neither allows \(X\to\infty\) at **fixed \(a_{\rm IV}\)**. In particular, the long-shift \(X=9\) test fails fixed-radius inclusion, while the new short-shift \(X=9\) test passes with an explicit buffer.

At a variable radius \(a(X)\), controlling the three-cutoff \(\mathcal J(X)\), both pole moments, complete Gamma/smoothing errors and **cross-window old/new mixed pairings** is necessary for any Phase IV support-extension claim. The new fixed-radius positivity result is on **one carefully constructed smooth direction** and supplies no infinite-dimensional complement or uniform old-collar inclusion inequality.

### Evidence classification
| Obligation | Status |
|---|---|
| Exact three-cutoff kernel identity and prime-sum reduction | **ANALYTIC / algebra audit PASS** |
| Exact long-filter two-pole cancellation and complete signed \(\mathcal W\) identity | **PASS** |
| Classical PNT corollary for the long-filter *individual family* | **CONDITIONAL ON KNOWN CLASSICAL PNT; not a CR arithmetic proof** |
| Uniform smooth long-filter bound with \(\eta_X=X^{-4}\) | **ANALYTIC, using CR-ALP-004** |
| Fixed \(11/10\)-radius compressed smooth pole-null source \(X=9\) | **ANALYTIC / radius certificate PASS** |
| Full original signed value \(\mathcal W(k)>(3/5)\|h_\eta\|_2^2\) on that one test | **DIRECT ANALYTIC PASS / rational finite audit** |
| Actual Phase IV Schur / all-complement support extension | **NOT DISCHARGED** |
| Live receiving Lean branch / completed Connes-Suzuki operator domain | **NOT INDEPENDENTLY VERIFIED** |
| New RH theorem or global Weil positivity | **NOT CLAIMED** |

**Next CR-ALP-018:** verify this explicit *compressed* test against the receiving Phase IV native physical inner-product and selected-cell \(H/V\) grid, including the complete mixed representer, and measure its interaction with the actual complementary subspace. A positive single vector is not a positive old-collar Schur complement. Reopen live receiving repository access or use a pinned current file export before any promotion. Alternatively formulate a variable-radius support extension with a sourced uniform old/new cross-operator estimate; absent this, retain the local test only.

**No changes to Weil-Lab/F41/F42/Phase IV acceptance or external RH/provenance ledgers.**
