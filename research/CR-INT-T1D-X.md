# CR-INT-T1D-X — Genuine central–satellite Weil coupling, the five-prime compensation, and the second Schur gate

**Canonical Return · 2026-11-23 (America/Los_Angeles) · independent smooth-source investigation.**  
**Disposition:** **DERIVED / ANALYTIC** source-preserving, exact rational/SymPy audit of central-to-satellite prime and Gamma cross estimates; explicit **original-Weil positive two-mode Gram in each parity sector**, plus an **undischarged** second Schur complement. The quantitative statement for the *Schur-corrected* cross pairing imports the separately unreviewed T1D-I global residual bound. **Not** a certified upper bound for the entire compact Birman–Schwinger spectrum, the full physical complement, the completed semilocal Sonin operator, Lean, Phase IV, or RH. Follows [T1D-W](CR-INT-T1D-W.md).

## 1. Freeze the unmodified source and the physical complement

Use Suzuki's **original complete signed Weil functional**, [*Weil's quadratic form via the screw function*, arXiv:2606.09096v4, §1.1, §2.3–2.4](https://arxiv.org/html/2606.09096v4), and the physical Hilbert space \(H=L^2((-2,2),dx)\). Write \(q(y,z)=W(y*\widetilde z)\), linear in its first argument, and retain all zeta-pole, signed prime-power and Gamma terms. Let

\[
\mathcal D^0=\{f\in C_c^\infty(-2,2):\int f(x)e^{\pm x/2}\,dx=0\},\quad
\mathscr H_F=(E\oplus P)^\perp,\quad F_s=\mathcal D^0\cap\mathscr H_F,
\]
with unchanged selected T1C pair \(E=\operatorname{span}\{u_+,u_-\}\) and physical two-pole carriers \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\). The original selected-space cross correction is \(\mathcal R\ge0\) (rank at most two); define \(\mathfrak s(y,z)=q(y,z)-\langle\mathcal R y,z\rangle\). T1D-R's full-complement problem is exactly \(\mathfrak s\ge0\) on all of \(F_s\), equivalent at the form-completion level to \(\sup\sigma(\mathscr B_F)\le1\).

Keep the actual **fixed** T1D-W compact smooth mollifier \(\rho_*\), even and nonnegative with unit integral, supported in \([-10^{-7},10^{-7}]\). For \(\sigma=+\) (even) or \(-\) (odd), put

\[
\begin{aligned}
 b_+(x)&=(1-16x^2)^3\mathbf1_{|x|<1/4},&
 b_-(x)&=x(1-16x^2)^3\mathbf1_{|x|<1/4},\\
 f_\sigma&=(\partial_x^2-\tfrac14)b_\sigma,&
 y_\sigma&=f_\sigma*\rho_*.
\end{aligned}\tag{X.1}
\]

These are the \(k=3\) vectors of the two parity spaces \(Y_\sigma=\operatorname{span}\{y_{\sigma,3},y_{\sigma,4}\}\) already certified by T1D-W, not substituted examples. They satisfy both Mellin moments **exactly**, and their support \(J=[-1/4-\epsilon,1/4+\epsilon]\) is physically disjoint from \(E\).

## 2. New true complementary satellite, still inside the original interval

With \(d=\log3\), define

\[
\boxed{z_\sigma=\frac{T_dy_\sigma+T_{-d}y_\sigma}{\sqrt2},\qquad T_tf(x)=f(x-t).}\tag{X.2}
\]

Its two satellite supports lie in \(\pm[d-1/4-\epsilon,d+1/4+\epsilon]\subset(-2,2)\), are mutually disjoint, and are disjoint from every original narrow selected bump and the complete central space \(Y_+\oplus Y_-\). Translations preserve both Mellin pole-null equations; reflection preserves them while exchanging their signs. Therefore \(z_+\) is even, \(z_-\) is odd and

\[
\boxed{z_\sigma\in F_s\cap (Y_+\oplus Y_-)^\perp,\quad
\|z_\sigma\|_2=\|y_\sigma\|_2,\quad \langle y_\sigma,z_\sigma\rangle=0.}\tag{X.3}
\]

The existing *physical* decomposition is consequently genuine, not an informal Fourier projection.

## 3. Exact prime census: the crucial **five-prime** satellite compensation

The original arithmetic term is \(-A\), with
\[
A=\sum_{n=p^k\le54}\kappa_nP_I(T_{\log n}+T_{-\log n})P_I,\quad \kappa_n=\frac{\Lambda(n)}{\sqrt n}>0.
\]
The checker enumerates **all 24 prime powers**, verifies every strict support comparison by exact rational bounds for logarithms, and proves:

* **Central to one satellite:** only \(n=2,3,4\) contribute, because \(|\log(n/3)|<1/2+2\epsilon\).
* **Between the two satellites:** only \(n=7,8,9,11,13\) contribute, because \(|\log(n/9)|<1/2+2\epsilon\).
* **Within the central or an individual satellite copy:** no prime self-pairing survives, since support diameter \(1/2+2\epsilon<\log2\).

To retain the nontrivial overlap amplitudes, define the *true mollified autocorrelation* \(R_\sigma(t)=\langle y_\sigma,T_ty_\sigma\rangle\), \(m_\sigma=\|y_\sigma\|_2^2\). The exact original prime contributions are

\[
\boxed{\begin{aligned}
 p_\sigma &:=-\frac{q_{\rm prime}(y_\sigma,T_d y_\sigma)}{m_\sigma}
  =\kappa_3+\kappa_2\frac{R_\sigma(\log(3/2))}{m_\sigma}
   +\kappa_4\frac{R_\sigma(\log(4/3))}{m_\sigma},\\
 w_\sigma &:=\frac{q_{\rm prime}(T_dy_\sigma,T_{-d}y_\sigma)}{m_\sigma}
  =-\sum_{n\in\{7,8,9,11,13\}}\kappa_n
 \frac{R_\sigma(|\log(9/n)|)}{m_\sigma}.
\end{aligned}}\tag{X.4}
\]

All the \(R_\sigma\) on the right are **actual smooth-source** pairings, not raw replacements. For verification, \(f_\sigma\) is a rational polynomial inside \([-1/4,1/4]\), so its raw autocorrelation at \(0\le t\le1/2\) is an **exact rational polynomial**:
\[
R^0_\sigma(t)=\int_{t-1/4}^{1/4}f_\sigma(x)f_\sigma(x-t)\,dx.
\]
The audit constructs that polynomial symbolically. Since \(\|y_\sigma-f_\sigma\|_2\le\epsilon\|f_\sigma'\|_2\), all shifted autocorrelations differ from their raw counterparts by at most \(2\epsilon\sqrt{M^0_\sigma K^0_\sigma}+\epsilon^2 K^0_\sigma\). Direct rational log-series, integer-root, and interval evaluation then certify

\[
\boxed{\begin{array}{c|cc}
 &p_\sigma&w_\sigma\\\hline
+&(0.674180724764,\ 0.674188109421)&(0.278607148935,\ 0.278632756572)\\
-&(0.705659054647,\ 0.705669062832)&(-0.097918626680,\ -0.097888838619)
\end{array}}\tag{X.5}
\]

The decimals displayed are outward diagnostics from rational intervals. In particular the clean rational assertions are \(0.67<p_+<0.68\), \(0.70<p_-<0.71\), \(0.27<w_+<0.29\), and \(-0.11<w_-<-0.09\).

**Crucial pitfall:** retaining only the visually obvious prime-9 interaction between the two satellite copies would give \(-\kappa_9=-\log3/3\approx-0.3662\), but the **complete** prime sum is actually \(+0.2786\) in the even sector and \(-0.0979\) in the odd sector. The contributions from \(7,8,11,13\) are necessary. Omitting them produces a *spurious negative finite determinant* and cannot support any original Weil claim.

## 4. A new rigorous bound for nonlocal Gamma coupling

On disjoint supports, the original Suzuki Gamma kernel is \(-w(|x-t|)\) with \(w(s)=e^{-s/2}/(1-e^{-2s})\). Write \(y_\sigma=(\partial^2-\tfrac14)B_\sigma\) where \(B_\sigma=b_\sigma*\rho_*\). Moving the differential operator onto the smooth kernel twice (both source functions are compact smooth) yields

\[
|\gamma(y_\sigma,T_t y_\sigma)|
\le H(t-1/2-2\epsilon)\|B_\sigma\|_1^2,
\quad H(s):=(\partial_s^2-\tfrac14)^2w(s)
=\sum_{j\ge1}[2j(2j+1)]^2e^{-(2j+1/2)s}.\tag{X.6}
\]

The \(j=0\) term vanishes **exactly** after the two pole-annihilating differential operators. This is why the naive bound \(w(\text{separation})\|y\|_1^2\) wastes almost all the cancellation. For \(q=e^{-2s}\) the remaining series has the closed positive generating-function majorant
\[
H(s)\le16\frac{q(1+11q+11q^2+q^3)}{(1-q)^5}
+16\frac{q(1+4q+q^2)}{(1-q)^4}
+4\frac{q(1+q)}{(1-q)^3}.
\]

An exact rational Taylor lower enclosure for \(e^{1.18}\) and \(e^{3.38}\) proves
\(H(d-1/2-2\epsilon)<231\), \(H(2d-1/2-2\epsilon)<3\). The original polynomial \(L^1\) norms are **exact**:
\[
\|b_+\|_1=8/35,\quad \|b_-\|_1=1/64,
\]
while the original polynomial masses from T1D-V are \(M_+^0=28194976/15015\) and \(M_-^0=184586/3465\). Since convolution by \(\rho_*\) cannot increase \(L^1\) mass, and the mollified physical masses have a strict rational lower bound, (X.6) proves **uniformly in both parities**

\[
\boxed{\frac{|\gamma(y_\sigma,T_dy_\sigma)|}{m_\sigma}<\frac1{100},\qquad
\frac{|\gamma(T_dy_\sigma,T_{-d}y_\sigma)|}{m_\sigma}<\frac1{1000}.}\tag{X.7}
\]

The audit computes the much tighter absolute upper diagnostics \(0.006428\) (even) and \(0.001059\) (odd) for the *near* coupling. These are rigorous bounds, not FFT samples. The proof uses the **full original Gamma distribution** and the two-pole differential annihilation, with no new signed term.

## 5. A positive two-mode Weil Gram, but a large off-block interaction

Normalize \(\widehat y_\sigma=y_\sigma/\sqrt{m_\sigma}\), \(\widehat z_\sigma=z_\sigma/\sqrt{m_\sigma}\). From T1D-W's **independently audited smooth Gamma diagonal** for the same \(k=3\) fixed mollifier, we safely import
\[
0.89<\gamma[\widehat y_+]<0.91,
\qquad1.18<\gamma[\widehat y_-]<1.20.
\tag{X.8}
\]

Equations (X.4), (X.7), and translation invariance of the complete Gamma form give the **exact full Weil** physical-orthonormal \(2\times2\) Gram formulas
\[
\boxed{\begin{aligned}
q(\widehat y_\sigma,\widehat y_\sigma)&=g_\sigma,\\
q(\widehat z_\sigma,\widehat z_\sigma)&=g_\sigma+w_\sigma+\gamma(T_dy_\sigma,T_{-d}y_\sigma)/m_\sigma,\\
q(\widehat y_\sigma,\widehat z_\sigma)&=\sqrt2\left[-p_\sigma+\gamma(y_\sigma,T_dy_\sigma)/m_\sigma\right].
\end{aligned}}\tag{X.9}
\]

These are **complete original Weil** entries; the Mellin pole terms vanish by exact source constraints. Rational interval propagation gives
\[
\begin{array}{c|ccc}
 &q(\widehat y,\widehat y)&q(\widehat z,\widehat z)&|q(\widehat y,\widehat z)|\\\hline
+&(0.89,0.91)&(1.159,1.201)&(\sqrt2\cdot0.66,\sqrt2\cdot0.69)\\
-&(1.18,1.20)&(1.069,1.111)&(\sqrt2\cdot0.69,\sqrt2\cdot0.72)
\end{array}
\]
so the true Gram determinant is **strictly positive in both parities**:
\[
\boxed{\det Q_{\{\widehat y_+,\widehat z_+\}}>\frac7{100},\qquad
\det Q_{\{\widehat y_-,\widehat z_-\}}>\frac15.}\tag{X.10}
\]

This establishes positive original Weil form on the **four-dimensional** authentic span of these two central/satellite pairs. By contrast, the off-block mixed Weil entries themselves have large magnitudes: they are bounded below by \(\sqrt2\cdot0.66\) and \(\sqrt2\cdot0.69\). **A positive central Gram is not decoupled from the true remainder.** (X.10) relies on the missing wing-channel arithmetic being included; it is *not* a theorem on the entire physical complement.

### Actual selected-space Schur remainder: a conditional nonzero coupling certificate

T1D-W directly established on the entire central two-spaces
\(\mathcal R[\widehat y_+]<4/1575\), \(\mathcal R[\widehat y_-]<9/4375\) without invoking T1D-I. **If in addition** T1D-I's separately derived, unreviewed global estimate \(\|\mathcal R\|<155/182\) is correct, physical Cauchy–Schwarz gives \(|\langle\mathcal R\widehat y,\widehat z\rangle|<1/20\). Then from (X.9):
\[
\boxed{\begin{aligned}
|\mathfrak s(\widehat y_+,\widehat z_+)|&>\frac{437}{500}>\frac{17}{20},\\
|\mathfrak s(\widehat y_-,\widehat z_-)|&>\frac{229}{250}>\frac9{10}.
\end{aligned}}\tag{X.11}
\]

The **original Weil q cross** inequalities and (X.10) are independent of I; only the numbers in (X.11) import I. These genuine off-block witnesses prove the original central subspace \(Y=Y_+\oplus Y_-\) is *not* a decoupled block for \(\mathfrak s\), conditional on that residual estimate.

### Necessary satellite threshold if the entire complement were positive

The central diagonal satisfies \(0<\mathfrak s[\widehat y_+]<0.91\), \(0<\mathfrak s[\widehat y_-]<1.20\) from W's positive Schur result. If \(\mathfrak s\ge0\) held on the entire complement, the \(2\times2\) Schur determinant with these actual \(\widehat z\) vectors would **necessarily** require
\[
\boxed{\mathfrak s[\widehat z_+]>
\frac{(437/500)^2}{91/100}=\frac{190969}{227500}>0.83,
\qquad
\mathfrak s[\widehat z_-]>
\frac{(229/250)^2}{6/5}=\frac{52441}{75000}>0.69.}\tag{X.12}
\]
This is a **conditional necessary test**, not a newly established positive satellite Schur lower bound. Current global residual bounds are too coarse to discharge it.

## 6. Exact second Schur gate and custody

For a genuine physical orthogonal decomposition \(F_s=Y\oplus Z_s\) understood algebraically with \(Z_s=F_s\cap Y^\perp\), the **finite** matrix \(S_Y=(\mathfrak s(y_i,y_j))\) on \(Y\) is positive by T1D-W. Let \(b(z)=(\mathfrak s(y_i,z))_i\). Completing the square gives the **exact missing obligation**
\[
\boxed{\mathfrak s(z,z)-b(z)^*S_Y^{-1}b(z)\ge0
\quad\text{for all } z\in Z_s.}\tag{X.13}
\]
Our witnesses \(\widehat z_\pm\in Z_s\) show why simply combining a positive central Ritz matrix with an unrelated high-frequency lower estimate is logically insufficient: the nonzero mixed block must be paid for. **No estimate on all of \(Z_s\)** is established, and the exact positive \(q\)-Gram in (X.10) does not by itself control the \(q-\mathcal R\) correction on the two-mode satellite span.

The accompanying `src/audits/canonical_return_t1d_x_satellite_coupling.py` verifies: complete 24-prime-power census; exact central `2,3,4` and wing `7,8,9,11,13` geometry; rational autocorrelation polynomials; mollifier transport; directed rational log and root enclosures; the new \(H(s)\) Gamma-kernel majorants; strict positive Weil Gram determinants; conditional nonzero Schur mixed bounds; and arithmetic of the necessary satellite thresholds. Its **analytical imports** are the original Suzuki Gamma kernel, T1D-W's diagonal Gamma and central residual enclosures, and (only for X.11–X.12) T1D-I's global residual bound. Re-run W independently; the script does not formalize infinite-dimensional form closure, or compute T1D-T's \(\eta_N\).

**Open:** exact true Schur cross-vector for the full central basis, \(\mathcal R[\widehat z]\) without T1D-I, the second Schur form inequality on the entire \(Z_s\), a finite Gamma spectral cutoff with usable \(\eta_N\), completed Suzuki–Connes form/metric/domain intertwining, global Weil positivity and RH.

**Next CR-INT-T1D-Y:** use the explicit original prime `7,8,9,11,13` wing compensation to bound the **selected-to-satellite** Weil functional and the actual two-mode **Schur-corrected** determinant without importing T1D-I. The second Schur test must then extend from one satellite witness to the rest of \(Z_s\) through a genuine spectral-tail or operator estimate; no finite positive example is a whole-complement certificate.