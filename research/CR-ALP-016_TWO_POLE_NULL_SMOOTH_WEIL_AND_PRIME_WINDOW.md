# CR-ALP-016 — Exact two-pole-null smooth signed Weil test and prime-window expansion

**Jurisdiction:** Canonical-Return, independent analytical return; handoff candidate for Phase IV **only** under the established evidence-transfer gate.  
**Date:** 2026-11-23 PST.  
**Standing:** DERIVED / EXACT RATIONAL AUDIT PASS (locally executed) / source-preserving finite test. **No Lean kernel run, global Weil positivity, Connes/Suzuki core admission, arithmetic estimate (M), independent RH proof, or support-extension theorem.**

## 0. Motivation and original custody

CR-ALP-015 stopped the *full arithmetic-to-Weil complex operator correspondence* on insufficient evidence. The remaining valid source channel is the *conditional scalar* reconstruction from CR-ALP-002–004. One exact gap has a natural mathematical test: the source ramp \(f_X\) has nonzero Mellin moments at both original zeta poles. Can one **explicitly produce a genuinely smooth, two-pole-null source** while preserving the complete **original signed Weil form**, and determine the *arithmetic price* of doing so?

The source convention is exactly the one in:
- CR-ALP-002, \`research/CR-ALP-002_FULL_WEIL_RAMP_AND_POLE_NO_GO.md\`: Mellin moments \(F_{\pm}\), full original two-pole / von Mangoldt / complete Gamma Weil distribution.
- CR-ALP-003, \`research/CR-ALP-003_SMOOTH_CELL_GRAM_AND_LEAN_INTERMEDIATE.md\`: even mollifier and fixed-support continuity for the complete form.
- CR-ALP-004, \`research/CR-ALP-004_QUANTITATIVE_SMOOTH_WEIL_TRANSFER.md\`: **conditional** arithmetic \(G(X)\) upper bound transferred to an individual pole-nonnull smooth ramp; **no** global Weil-positive conclusion.
- CR-ALP-007, \`research/CR-ALP-007_TRIANGLE_GAUGE_HOLONOMY_AND_ARITHMETIC_PHASE.md\`: complete three-translated-ramp original signed \(3\times3\) Weil Gram, independently rechecked here.

No substitute spectral positivity functional is introduced.

## 1. Universal three-shift two-pole annihilator

For arbitrary \(X>1\), put \(L=\log X>0\), \(f_X(u)=e^{u/2}\mathbf1_{[0,L]}(u)\), and define its translates \(f_j(u)=f_X(u-jL)\), \(j=0,1,2\). Let \(q=\sqrt X\) and \(s=q+q^{-1}\). Define
\[
\boxed{h_X=f_0-sf_1+f_2.}\tag{16.1}
\]

The two original Mellin pole moments are
\[
F_+(g)=\int_{\mathbb R}g(u)e^{u/2}\,du,\qquad
F_-(g)=\int_{\mathbb R}g(u)e^{-u/2}\,du.
\]
Translation by \(b\) multiplies \(F_\pm\) by \(e^{\pm b/2}\). Since \(F_+(f_X)=X-1\) and \(F_-(f_X)=L\), we obtain
\[
\begin{aligned}
F_+(h_X)&=(X-1)(1-sq+q^2)=0,\\
F_-(h_X)&=L(1-s/q+q^{-2})=0.
\end{aligned}
\tag{16.2}
\]
The polynomial \(1-sz+z^2=(z-q)(z-q^{-1})\) is the **unique degree-two symmetric translation filter, up to a scalar, that kills both exponential moments**. Two noncoincident translates alone cannot do this, because their \(2\times2\) moment matrix has nonzero determinant.

The translates have disjoint support interiors, so their native \(L^2\) Gram is \((X-1)I_3\) and
\[
\boxed{\|h_X\|_2^2=(X-1)(s^2+2)>0.}\tag{16.3}
\]
The real autocorrelation \(k_{h_X}=h_X*\widetilde h_X\) is positive definite in the Fourier/Bochner sense:
\(\widehat{k_{h_X}}(t)=|\widehat f_X(t)|^2|1-se^{-itL}+e^{-2itL}|^2\ge0\).
**This is not a theorem about the sign of the original signed Weil form.**

The original zeta pole functional of the correlation is identically zero:
\[
\boxed{\mathcal P(k_{h_X})=
\int_{\mathbb R}k_{h_X}(y)(e^{y/2}+e^{-y/2})\,dy
=2F_+(h_X)F_-(h_X)=0.}\tag{16.4}
\]
Therefore the exact original signed value here consists of the **full signed prime and Gamma parts**, without silently discarding pole terms.

## 2. A complete finite signed Weil certificate at \(X=9\)

Set \(X=9\), \(q=3\), \(L=\log9\), \(s=10/3\). Then
\[
\boxed{h=f_0-\tfrac{10}{3}f_1+f_2,\qquad
\|h\|_2^2=\tfrac{944}{9}.}\tag{16.5}
\]

Let \(d,c_1,c_2\) be the **complete original signed Weil** diagonal, adjacent and distance-\(2L\) mixed entries computed in CR-ALP-007 (including both poles, every prime power and the full Gamma term):
\[
H=
\begin{pmatrix}d&c_1&c_2\\c_1&d&c_1\\c_2&c_1&d\end{pmatrix}.
\]
Exact polarization gives
\[
\boxed{\mathcal W(k_h)=
(2+s^2)d-4s\,c_1+2c_2
=\tfrac{118}{9}d-\tfrac{40}{3}c_1+2c_2.}\tag{16.6}
\]

The standalone accompanying audit **independently reconstructs** the three entries by enumerating all prime powers through 729 and using rational atanh/log enclosures, a Machin-series rational bracket for \(\pi\), and an Euler–harmonic rational enclosure. It certifies
\[
\boxed{7.7143<\mathcal W(k_h)<7.7145.}\tag{16.7}
\]
The numerical approximation is \(\mathcal W(k_h)\simeq7.71435689678134\), given only for orientation; the inequalities are rigorously rational.

The two-pole cancellation also checks *directly at the original-form level*: the diagonal, adjacent, and distant original pole coefficients are respectively \(16L,\frac{80}{3}L,\frac{656}{9}L\), so
\[
\frac{118}{9}\cdot16L-\frac{40}{3}\cdot\frac{80}{3}L+2\cdot\frac{656}{9}L=0.
\]
**Result:** an explicit nonzero two-pole-null \(L^2\) source with a strictly positive value of the **original signed Weil distribution**. This is a *one-source finite test*, not a support-independent positivity result or an RH proof.

## 3. The arithmetic price: the original cutoff \(X\) becomes \(X^3\)

For real \(t\ge0\) the original ramp correlation \(k_X(t)=(Xe^{-|t|/2}-e^{|t|/2})_+\) yields the exact identity
\[
\boxed{
k_{h_X}(t)=(2+s^2)k_X(t)-2s\,k_X(t-L)+k_X(t-2L).
}\tag{16.8}
\]
The terms \(k_X(t+L)\) and \(k_X(t+2L)\) vanish for \(t\ge0\), so (16.8) retains every mixed cross interaction, not only diagonal cells.

For \(X=9\), set \(r_j(n)=k_9(\log n-jL)/\sqrt n\), \(j=0,1,2\), with \(k_9\) even and zero outside \([-L,L]\). The exact prime sample is
\[
\frac{k_h(\log n)}{\sqrt n}
=\tfrac{118}{9}r_0(n)-\tfrac{20}{3}r_1(n)+r_2(n).
\tag{16.9}
\]
For \(n\ge2\), these are the completely explicit **rational** pieces:
\[
r_0(n)=
\begin{cases}9/n-1,&n\le9,\\0,&n>9,\end{cases}
\quad
r_1(n)=
\begin{cases}
3(1-1/n),&n\le9,\\
27/n-1/3,&9<n\le81,\\
0,&n>81,
\end{cases}
\]
\[
r_2(n)=
\begin{cases}
0,&n\le9,\\
1-9/n,&9<n\le81,\\
81/n-1/9,&81<n\le729,\\
0,&n>729.
\end{cases}
\tag{16.10}
\]
The signs of the *kernel samples* genuinely vary:
\[
\boxed{\frac{k_h(\log9)}3=-\frac{160}{9},\qquad
\frac{k_h(\log81)}9=\frac89.}
\]
The exact full prime part is
\[
-2\sum_{n\ge2}\Lambda(n)\frac{k_h(\log n)}{\sqrt n}.
\]
The finite prime-power inventory has **129 distinct prime bases** with nonzero aggregated coefficients (through primes \(\le727\)); the unsmoothed kernel vanishes at \(n=729=X^3\).

**Transfer obstruction:** the original CR-ALP-004 conditional estimate on the **single** weighted sum \(S(X)\), or on \(G(X)\), does **not** by itself bound the new signed mixed sums on \(X<n<X^3\). Pole-nullification changes the arithmetic object. It is a nontrivial, exact new source—not a convergence of pole-null tests to the old nonnull \(k_X\). Hence it does **not** contradict CR-ALP-002's weighted-\(L^1\) impossibility for such an approximation.

A support-extension theorem must control all those new mixed signed prime windows and the Gamma term; simply importing the \(S(X)\) estimate is invalid.

## 4. Effective smooth and two-pole-null witness

Choose **any** real, even, nonnegative unit-mass mollifier \(\rho_\eta\in C_c^\infty(-\eta,\eta)\), with the explicit width
\[
\boxed{\eta=10^{-9}.}
\]
Put
\[
h_\eta=\rho_\eta*h,\quad
k_{h_\eta}=h_\eta*\widetilde{h_\eta}
=(\rho_\eta*\widetilde\rho_\eta)*k_h.
\]
Since convolution multiplies each Mellin pole moment by \(\int\rho_\eta(v)e^{\pm v/2}dv\), both moments remain **exactly zero**:
\[
\boxed{h_\eta\in C_c^\infty(\mathbb R),\qquad
F_+(h_\eta)=F_-(h_\eta)=0.}\tag{16.11}
\]
The source support lies in \([-\eta,3L+\eta]\). Translation to center it produces an equally pole-null smooth test with support
\[
\boxed{[-\tfrac32\log9-\eta,\ \tfrac32\log9+\eta].}\tag{16.12}
\]
The autocorrelation support is inside \([-3L-2\eta,3L+2\eta]\).

### 4a. An explicit uniform smooth-sign certificate (not just existence)

For \(k_9\), the global Lipschitz constant is at most \(5\). The full even three-shift formula gives
\[
\operatorname{Lip}(k_h)\le5\left(\frac{118}{9}+\frac{40}{3}+2\right)
=\frac{1280}{9}<143.
\]
Let \(\Delta=k_{h_\eta}-k_h\). With \(\operatorname{supp}(\rho_\eta*\widetilde\rho_\eta)\subset[-2\eta,2\eta]\), we have
\[
D:=\|\Delta\|_\infty<286\eta,\qquad
\operatorname{Lip}(\Delta)<286.
\]
For \(\eta=10^{-9}\), the full original signed Weil difference can be bounded **without choosing the shape of \(\rho_\eta\)**:

- Both poles: \(|\mathcal P(\Delta)|<216D\), using \(e^{(3L)/2}=27\).
- Entire prime-power sum: at most \(n\le729\) contributes, \(\Lambda(n)\le\log n<8\), hence absolute error \(<11648D\).
- Full Gamma diagonal and integral: \(<5D+2300\eta+200D+12D\), splitting the integral into \(t\in(0,2\eta)\), \(t\in(2\eta,1)\), and \(t\ge1\); use \(e^{-t/2}/(1-e^{-2t})\le2/t\) for \(0<t\le1\), \(\le2e^{-t/2}\) for \(t\ge1\), and \(\log(1/(2\eta))<25\).

Therefore
\[
\boxed{
|\mathcal W(k_{h_\eta})-\mathcal W(k_h)|
<(216+11648+5+200+12)D+2300\eta
<3{,}500{,}000\eta=0.0035.
}\tag{16.13}
\]
Combining (16.7) with this *explicit universal smooth* error budget gives
\[
\boxed{\mathcal W(k_{h_\eta})>7.7>0.}\tag{16.14}
\]
This is a genuine smooth compact-support test satisfying **both algebraic pole-null conditions**, with the full original signed form and an effective mollifier-width threshold (at the stated \(\eta\)). The exact finite audit certifies the numerical constants; the Gamma Lipschitz/integral split is the separately stated elementary analytic proof.

It does **not** independently certify admission to a particular *completed Connes or Suzuki operator core*, including any additional semilocal or boundary hypotheses beyond smoothness, compact support, and the two pole moments.

## 5. Phase IV interface: precise benefit and explicit limit

| Transfer gate | CR-ALP-016 evidence | Status |
|---|---|---|
| Precise mathematical theorem | \(C_c^\infty\) pole-null source at \(\eta=10^{-9}\) with original signed \(\mathcal W(k)>7.7\) | **ANALYTIC / EXACT AUDIT PASS** |
| Preservation of original signed form | Complete zeta pole, von Mangoldt **all prime powers**, complete Gamma, and source \(L^2\) mass | **PASS for this one test** |
| Actual outstanding IV obligation | Control **arbitrary expanding supports** and all mixed negative prime windows, not merely \(X=9\) | **NOT DISCHARGED** |
| Independent verification | Standalone rational finite-prime and log-bound audit run locally; no Lean or original F41 physical H/V grid compilation | **LOCAL AUDIT PASS ONLY** |
| Connes/Suzuki operator-core membership | Necessary algebraic pole conditions and \(C_c^\infty\) hold; further analytic hypotheses untested | **OPEN** |
| Conditional Alpöge scalar inequality transfer | The new kernel requires prime windows through \(X^3\), not bounded by the sole \(G(X)\) input | **NO DIRECT TRANSFER** |

**Handoff disposition:** **ELIGIBLE AS A FALSIFIABLE LOCAL POLE-NULL TEST / NOT ELIGIBLE AS A PHASE IV SUPPORT-EXTENSION THEOREM.** This should remain a Canonical-Return research artifact until independently checked against IV's exact admissible class and support collar constraints. Formalization F41/F42 certification and original physical H/V matrices are unaffected.

**Next — CR-ALP-017:** formulate a *noncircular* family-level support-extension criterion for the three-shift pole-null filter, stating the exact additional mixed arithmetic estimate and uniform full Gamma/smoothing bounds required for \(X\to\infty\). Identify whether any existing Phase IV collars can carry the centered support of radius \(\frac32\log X+\eta\), or record a failed geometric preservation test. No equivalence or RH claim.

**Executable:** \`src/audits/cr_alp_016_pole_null_return_exact_audit.py\`. All sign/log comparisons are exact \`fractions.Fraction\` with rigorous atanh series tails, Euler-harmonic and Machin-\(\pi\) bounds; floating prints are display only.
