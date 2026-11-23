# CR-INT-T1D-I — Two-sided finite-width cross-block enclosures without quadrature

**Canonical Return | independent analytical research | 2026-11-23 PST.**  
**Status:** ANALYTICAL TWO-SIDED ENCLOSURE, exact-rational **AUDIT-PASS** of the finite geometry, rational inequalities, and imported T1D-G/H enclosures. **Not** independently interval-integrated residual norms, not Lean-checked, not an estimate for the infinite-dimensional complementary Weil form, and **not RH**. This note extends [T1D-E](CR-INT-T1D-E.md), [T1D-F](CR-INT-T1D-F.md), [T1D-G](CR-INT-T1D-G.md), [T1D-H](CR-INT-T1D-H.md), and [the original T1C profile](../archive/space/10-t1c.md). The source definition of the Weil form is [Suzuki, arXiv:2606.09096v4, §1.1](https://arxiv.org/html/2606.09096v4); its exact Gamma Fourier multiplier is [§2.4](https://arxiv.org/html/2606.09096v4). An elementary analytic tool below is the [DLMF digamma Binet integral, §5.9, (5.9.15)](https://dlmf.nist.gov/5.9.E15).

## 1. Exact fixed source, physical projection and goal

Fix \(I=(-2,2)\), \(\ell=10^{-4}\), and the actual real even nonnegative smooth unit bump \(v=\psi_0/\|\psi_0\|_2\), where \(\psi_0=b_{49/100}*\rho_{1/100}\) was frozen in T1C and has support \([-1/2,1/2]\). For \(d=\log 3\), \(x_j=(j-3/2)d\), and \(U_{\ell,x}v(t)=\ell^{-1/2}v((t-x)/\ell)\), retain the actual even/odd two-Mellin-pole-null selected vectors
\[
 u_{\sigma}=N_\sigma^{-1/2}\sum_{j=0}^3 c_{\sigma,j}U_{\ell,x_j}v,
 \qquad c_+=(1,1-s,1-s,1),\quad c_-=(1,-1-s,1+s,-1),\quad s=4/\sqrt3,
\]
\[N_\pm=(44\mp16\sqrt3)/3.\]
Let \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\subset L^2(I)\), \(E=\operatorname{span}\{u_+,u_-\}\perp P\), and \(F=(P\oplus E)^\perp\). Use the **source-defined** full Weil mixed representer \(\mathcal K_2\) from T1D-E, *not* Suzuki's inverse Neumann Laplacian. The actual quantities are
\[
\beta_\sigma^2=\|\Pi_F\mathcal K_2u_\sigma\|_{L^2(I)}^2,\qquad
\|\mathcal R\|=\max\{\beta_+^2/\mu_+,\ \beta_-^2/\mu_-\}.
\]
The 2026 T1D-H paper-like note had strict **lower** finite-width estimates but no two-sided interval. We produce **coarse but genuine two-sided bounds without replacing the finite-width norms by their asymptotic limits**.

## 2. A quantitative digamma multiplier remainder valid at the actual width

The exact Gamma-only whole-line Fourier multiplier from Suzuki is
\[
 m_\Gamma(\xi)=\operatorname{Re}\psi_{\rm dg}\bigl(\tfrac14+i\xi/2\bigr)-\log\pi.
\]
Write \(L=\log|D|-\log(2\pi)\), \(\Delta(\xi)=m_\Gamma(\xi)-(\log|\xi|-\log(2\pi))\) for \(\xi\ne0\). **Two explicit multiplier inequalities** are
\[
 |\Delta(\xi)|<11+|\log|\xi||\quad(0<|\xi|\le1),\qquad
 |\Delta(\xi)|<\frac{3}{2|\xi|}\quad(|\xi|\ge1). \tag{I.1}
\]
For the second, DLMF's Binet formula for \(\psi(z)\), \(z=1/4+i\xi/2\), gives
\[
 \psi(z)=\log z-\frac1{2z}-2\int_0^\infty \frac{t\,dt}{(t^2+z^2)(e^{2\pi t}-1)}.
\]
Since \(|\operatorname{Im}z^2|=|\xi|/4\), \(\int_0^\infty t/(e^{2\pi t}-1)dt=1/24\), \(|\Re(1/(2z))|\le1/|\xi|\), and \(|\log|z|-\log(|\xi|/2)|\le1/(8\xi^2)\), the difference is bounded by \(1/(8\xi^2)+1/|\xi|+1/(3|\xi|)<3/(2|\xi|)\) for \(|\xi|\ge1\). For the first, \(\psi(z)=-\gamma-1/z+\sum_{n\ge1}z/[n(n+z)]\) and \(\Re z>0\) imply \(|\psi(z)|<1+4+(3/5)\sum n^{-2}<7\), so \(|m_\Gamma|<9\) and \(|\log(2\pi)|<2\).

By compact support of \(v\) (unit length) and the unitary Fourier convention, \(|\widehat v(\eta)|^2\le 1/(2\pi)\). Splitting the Fourier integral at \(|\eta|=\ell\) and using \(\int_0^1(11-\log t)^2dt=145\), obtain the fully explicit one-cell error
\[
\boxed{\|[m_\Gamma(D/\ell)-\log(1/\ell)-L]v\|_2^2
 <\frac{589\ell}{4\pi}<\frac{589\ell}{12}<\left(\frac{71}{1000}\right)^2.}\tag{I.2}
\]
No Fourier quadrature or asymptotic remainder \(o(1)\) has been used. For the four normalized coefficient vectors, Cauchy–Schwarz yields \(\|c_\sigma\|_1/\sqrt{N_\sigma}\le2\); the total Gamma-multiplier replacement error has norm **\(<0.142\)** in \(L^2(I)\).

## 3. Localize the logarithmic Fourier shape with a second explicit error

Let \(a_v=\langle Lv,v\rangle\), \(h=(L-a_v)v\), so \(\langle h,v\rangle=0\) and \(\|h\|_2^2=V_v\). The real-space convolution kernel of \(\log|D|\) **off the support** of a smooth compact test is \(-1/(2|x|)\) with this Fourier convention; thus for \(|x|>1/2\),
\[
 h(x)=-\frac12\int_{-1/2}^{1/2}\frac{v(y)}{|x-y|}dy,
 \qquad |h(x)|\le\frac{\|v\|_1}{2(|x|-1/2)}\le\frac1{2(|x|-1/2)}.
\]
Choose physical radius \(r=1/100\) and truncate \(h\) to \(|x|<r/\ell\) in local variables. The omitted whole-line \(L^2\) mass obeys
\[
 \boxed{\tau^2\le\frac{\ell}{2(r-\ell/2)}<\left(\frac{71}{1000}\right)^2.}\tag{I.3}
\]
On the four normalized cells the truncation error is \(<0.142\). Combining (I.2) and (I.3), the true projected Weil image differs in norm by **\(<71/250=0.284\)** from the projection of the compactly supported proxy below.

## 4. An exactly orthogonal finite-translation proxy

Let \(P_3\) be the original genuine \(3,9,27\)-power symmetric prime Toeplitz matrix of T1D-G, \(\lambda_{3,\sigma}=c_\sigma^TP_3c_\sigma/N_\sigma\), and \(T_{\sigma,\mathrm{non3}}\) the complete original prime translation image from non-3 prime powers \(n\le54\), restricted to \(I\). Define
\[
\boxed{\begin{aligned}
 t_\sigma={}&\frac1{\sqrt{N_\sigma}}\sum_{j=0}^3 c_{\sigma,j}\,U_{\ell,x_j}(h\mathbf1_{|\cdot|<r/\ell})\\
 &+\frac1{\sqrt{N_\sigma}}\sum_{j=0}^3[(P_3-\lambda_{3,\sigma}I)c_\sigma]_j U_{\ell,x_j}v
 +T_{\sigma,\mathrm{non3}}.
\end{aligned}}\tag{I.4}
\]
The finite exact arithmetic census has **48** in-window non-ternary translated centers, with the **smallest multiplicative separation ratio \(37/36\)** between distinct centers, including the four original centers. Therefore the smallest logarithmic center distance is \(\log(37/36)>2/73>2r+\ell\). The closest translated center lies more than \(r\) inside the physical boundary \((-2,2)\). The 4 base neighborhoods of radius \(r\) and 48 non-ternary bumps are consequently mutually disjoint. This is certified *by exhaustive rational comparison* in the attached script; it does not depend on machine floating-point comparisons of center positions.

Because \(\langle h,v\rangle=0\), the second term has coefficient vector orthogonal to \(c_\sigma\), and the non-ternary centers are disjoint from the four original centers, **\(t_\sigma\perp E\)** exactly. Since the original profile is even and the Weil translation operator respects reflection, \(t_+\) is even and \(t_-\) odd. The three terms in (I.4) are mutually orthogonal by source support or \(h\perp v\), and the exact norm is
\[
\boxed{\|t_\sigma\|_2^2=V_v-\tau_{\rm actual}^2+V_{3,\sigma}+C_\sigma,\quad0\le\tau_{\rm actual}^2<\ell/[2(r-\ell/2)].}\tag{I.5}
\]
**Important:** this equality belongs to the *constructed finite-width proxy*, not to the actual untruncated residual and not to a substitution of the T1D-G limiting norm into the finite-width formula.

The T1D-G exact-rational census gives \(V_{3,\pm}+C_\pm\), and T1D-H independently enclosed the true smooth \(V_v\). Reusing their certified outward fractions yields
\[
 \begin{array}{lcl}
 4.82035748<T_+:=V_v+V_{3,+}+C_+<4.82993688,\\
 3.82338534<T_-:=V_v+V_{3,-}+C_-<3.83296474.
 \end{array}\tag{I.6}
\]
These are *reference-mass* enclosures for (I.5), not finite-width residual values.

## 5. The pole projection costs less than 5% in norm

Each \(t_\sigma\) is supported on a union of length at most
\[
D=8r+48\ell=\frac{53}{625}=0.0848.
\]
For the even pole direction \(p_e(x)=\cosh(x/2)\), \(\|p_e\|_2^2=\sinh2+2>5\) and \(\|p_e\|_\infty<8/5\). For the odd direction \(p_o(x)=\sinh(x/2)\), \(\|p_o\|_2^2=\sinh2-2>3/2\) and \(\|p_o\|_\infty<6/5\). Parity rules out the opposite pole, so Cauchy–Schwarz on the support gives
\[
\frac{\|\Pi_Pt_\pm\|_2^2}{\|t_\pm\|_2^2}
 \le D\max\!\left\{\frac{(8/5)^2}{5},\frac{(6/5)^2}{3/2}\right\}<\frac{39}{400}.
\]
Together with \(t_\pm\perp E\), this proves
\[
 \boxed{\|\Pi_Ft_\pm\|_2>\frac{19}{20}\|t_\pm\|_2.}\tag{I.7}
\]
The true projected residual obeys
\[
 \boxed{\left|\beta_\pm-\|\Pi_Ft_\pm\|_2\right|<\frac{71}{250}.}\tag{I.8}
\]
Every inequality applies **at the actual fixed width** \(\ell=10^{-4}\).

## 6. Main result: two-sided actual finite-width enclosures

From (I.3)–(I.8), entirely with rational comparisons, one obtains
\[
\boxed{\frac{16}{5}<\beta_+^2(10^{-4})<\frac{31}{5},\qquad
\frac{12}{5}<\beta_-^2(10^{-4})<\frac{51}{10}.}\tag{I.9}
\]
That is,
\[
\boxed{3.2<\beta_+^2<6.2,\qquad2.4<\beta_-^2<5.1.}
\]
This **does answer the former two-sided-enclosure question** with a coarse analytical certificate. It does not yet give *high-precision* enclosures or certify the actual scalar integrals by quadrature. It is much sharper than the earlier overall crude \(\|\mathcal R\|<1785\) and independently stronger than H's one-sided lower witnesses.

Using the independently certified selected Rayleigh bounds from T1C, \(7.28<\mu_+<13.93\) and \(7<\mu_-\), and exact parity diagonalization, we conclude the corresponding genuinely finite-width Schur norm satisfies
\[
 \boxed{\frac{320}{1393}<\|\mathcal R_{10^{-4}}\|<\frac{155}{182}.}\tag{I.10}
\]
Numerically, this is \(0.2297\ldots<\|\mathcal R\|<0.8517\ldots\). The norm is **neither zero nor greater than one** for the frozen selected source. These inequalities do **not** control the complementary quadratic form.

## 7. What this says about the missing Suzuki–Connes interface

A **sufficient norm-based condition** for closing the pole-null source Schur criterion is still a correct *whole-complement* bound
\[
 q(y,y)\ge d\|y\|_2^2\quad\text{for all admissible }y\in F,
 \qquad d\ge\|\mathcal R_{10^{-4}}\|.
\]
Because (I.10) gives \(\|\mathcal R\|<155/182\), a separately proved complement margin \(d\ge155/182\) would suffice. However no such positive \(d\)—nor even \(d>0\)—has been proved on this **actual physical complement**. Conversely a proposed lower estimate \(d\le320/1393\) **cannot by itself verify** the norm domination route, since \(\|\mathcal R\|>320/1393\); it does **not** disprove the stronger direction-dependent Schur inequality itself. No identification of this physical complement with Connes's Sonin/prolate-conditioned complement has been made, nor of the completed Suzuki and Connes operator domains. The original **two pole directions** are not covered by the pole-null smooth criterion. See [Andrews–Connes–La Gatta 2026, §8.6](https://arxiv.org/html/2611.13129v1) for the separately conditional semilocal Schur program.

**Disposition:** finite physical mixed-residual two-sided intervals **ANALYTIC PASS**, independent exact rational propagation and prime geometry **AUDIT-PASS**; actual high-precision residual quadrature, complement coercivity, semilocal (H2′)/R1/R2, completed-domain intertwiner, and full RH remain **OPEN**. The implementation `src/audits/canonical_return_t1d_i_residual_bracket.py` verifies the rational bounds (and imports the previously certified outward fractions from G/H, which are run separately in the regression runner). It does **not** machine-verify Binet, Fourier distribution identities, or all infinite-dimensional analytic steps.

**Next CR-INT-T1D-J:** decide the sign of the complementary form in a **source-defined, explicitly testable subspace**, or construct a genuine coercivity mechanism for an appropriately conditioned complement; as a secondary numerical improvement, replace the conservative projection and Binet bounds with certified quadrature for the four scalar residual observables. Do not infer a full-space theorem merely from the now controlled finite-rank penalty.