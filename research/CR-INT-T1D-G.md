# CR-INT-T1D-G — Exact sharp-rate residual law, non-ternary arithmetic floor, and complement-frequency obstruction

**Canonical Return · independent analytic source investigation · 2026-11-22 PST / 2026-11-23 UTC.**  
**Status:** **ANALYTIC DERIVATION** of an asymptotic law for the fixed-shape four-cell pole-null family; **AUDIT-PASS** for its exact finite prime-power combinatorics, symbolic ternary variance, and outward rational numerical bounds. **NOT** an independently certified numerical evaluation of the *actual finite-width* projected residual Gram of T1C, an infinite-complement positivity theorem, a completed Suzuki–Connes operator equivalence, Lean certification, or RH. This is a source-level refinement of [T1D-E](CR-INT-T1D-E.md) and [T1D-F](CR-INT-T1D-F.md); the draft PR stays unmerged.

## 1. Source pin and exact question

Use the **full original signed Weil explicit formula** as written in Suzuki, *Weil's quadratic form via the screw function*, [arXiv:2606.09096v4 §1.1](https://arxiv.org/html/2606.09096v4), and its verified Fourier/digamma reduction **§2.4, eqs. (2.7)–(2.8)**. In the present notations, the **Gamma-only** whole-line Fourier multiplier is exactly

\[
\boxed{m_\Gamma(\xi)=\Re\psi_{\rm dg}\!\left(\frac14+\frac{i\xi}{2}\right)-\log\pi,\quad
m_\Gamma(\xi)=\log|\xi|-\log(2\pi)+O(|\xi|^{-2})\ \ (|\xi|\to\infty),}
\tag{G.1}
\]

where \(\psi_{\rm dg}\) is the *digamma function*, not the test profile. These are **canonical published identities**, not a new discovery. The remaining original signed prime part on \(I=(-2,2)\) is the finite, exact translation sum \(-\sum_{n\le54}\Lambda(n)n^{-1/2}[T_{\log n}+T_{-\log n}]\) for the 24 genuine prime powers. This follows because \(54<e^{4}<55\), and the zero extension of any test in \(C_c^\infty(I)\) has no prime correlations at shifts outside \((-4,4)\). Both pole terms vanish *for the two-moment-null selected sources*, as already verified in T1B–F.

The previously proved bound \(\beta_\pm(\ell)\le136\) shows \(\|\mathcal R_\ell\|=O(1/\log(1/\ell))\) but does not say whether the **mixed residual norms** \(\beta_\pm\) themselves tend to zero. Here we derive their **positive exact limits** and identify the arithmetic part of those limits without any numerical integration.

## 2. Freeze the rescaled T1C carrier and its native metric

Let \(\psi_0=b_{49/100}*\rho_{1/100}\) be the even nonnegative compactly supported smooth test profile of T1C, with support in \([-1/2,1/2]\), and let \(v=\psi_0/\|\psi_0\|_2\). For \(0<\ell\le1/16\), set \(\phi_\ell(x)=\psi_0(x/\ell)\). Put \(d=\log3\), \(x_j=(j-3/2)d\) for \(j=0,1,2,3\), and define the unitary localization \(U_{\ell,x}v(t)=\ell^{-1/2}v((t-x)/\ell)\). Our real, physically normalized even/odd test vectors are exactly

\[
\boxed{u_{\sigma,\ell}=N_\sigma^{-1/2}\sum_{j=0}^3 c_{\sigma,j}U_{\ell,x_j}v,\quad
c_+=(1,1-s,1-s,1),\quad c_-=(1,-1-s,1+s,-1),\ s=4/\sqrt3,}
\tag{G.2}
\]

with \(N_\pm=(44\mp16\sqrt3)/3\). Their two Mellin moments vanish, they are orthonormal in physical \(L^2(I)\), and their complete Weil Rayleigh values are \(\mu_\pm(\ell)>0\). Write \(E_\ell=\operatorname{span}\{u_{+,\ell},u_{-,\ell}\}\), \(P=\operatorname{span}\{e^{\pm x/2}\mathbf 1_I\}\), \(\mathcal F_\ell=(E_\ell\oplus P)^\perp\), and \(r_{\sigma,\ell}=\Pi_{\mathcal F_\ell}\mathcal K_2u_{\sigma,\ell}\). All projections use **the native physical \(L^2(-2,2)\) metric**, not the Weil form or Suzuki's inverse Neumann Laplacian.

Parity remains exact: \(r_{+,\ell}\) is even, \(r_{-,\ell}\) is odd, so the mixed Gram and Schur penalty remain diagonal and \(\|\mathcal R_\ell\|=\max_\sigma\|r_{\sigma,\ell}\|_2^2/\mu_\sigma(\ell)\).

## 3. The shape-intrinsic logarithmic variance (an exactly defined finite number)

On whole-line \(L^2(\mathbb R)\), define \(Lv\) by its unitary Fourier transform

\[
\widehat{Lv}(\eta)=\bigl(\log|\eta|-\log(2\pi)\bigr)\widehat v(\eta).
\]

Because \(v\) is smooth and compactly supported, \(Lv\in L^2\) and \(v\) belongs to the domain of \(L^2\). Put \(a_v=\langle Lv,v\rangle\) and

\[
\boxed{V_v=\|Lv\|_2^2-a_v^2
=\int\left(\log|\eta|-\int\log|\zeta||\widehat v(\zeta)|^2d\zeta\right)^2|\widehat v(\eta)|^2d\eta>0.}
\tag{G.3}
\]

The Fourier transform is unitary, so \(\int|\widehat v|^2=1\). Strict positivity is immediate: \(\log|\eta|\) is not a.e. constant on the nonzero Fourier density \(|\widehat v|^2\). **We have not computed a certified numerical interval for \(V_v\).** Its appearance in the final formula is an exact analytically defined quantity, not an undisclosed quadrature approximation.

## 4. Exact base-cell ternary arithmetic variance

Let \(P_3\in\mathbb R^{4\times4}\) be the Toeplitz matrix with zero diagonal and \((P_3)_{ij}=-\kappa_{3^{|i-j|}}\) for \(i\ne j\), where \(\kappa_{3^k}=\log3/3^{k/2}\). For \(\sigma=+,-\), define the normalized selected prime offset \(\lambda_{3,\sigma}=c_\sigma^TP_3c_\sigma/N_\sigma\), and the **base-center prime residual**

\[
V_{3,\sigma}=\frac{\|(P_3-\lambda_{3,\sigma}I)c_\sigma\|_{\ell^2}^2}{N_\sigma}.
\]

Exact SymPy algebra establishes

\[
\boxed{\begin{aligned}
\lambda_{3,\pm}&=\frac{(41\mp5\sqrt3)\log3}{73},\\
V_{3,+}&=\frac{507+264\sqrt3}{5329}(\log3)^2,\\
V_{3,-}&=\frac{507-264\sqrt3}{5329}(\log3)^2.
\end{aligned}}
\tag{G.4}
\]

These variances are positive, and arise from three true prime powers 3, 9 and 27. They are *not* the selected generalized eigenvalue gap itself. The full non-ternary shifted support must also be counted.

## 5. All non-ternary prime shifts have an exact discrete support census

For every genuine prime power \(n\le54\) whose base prime \(p\ne3\), translation by \(\pm\log n\) moves the four narrow cells to new centers \(x_j\pm\log n\). A center contributes to the window only when it lies strictly inside \((-2,2)\). By the arithmetic identity \(x_i-x_j=(i-j)\log3\), **no two such non-ternary translated centers coincide with one another or with a base center**: an equality would require a product or ratio of non-ternary prime powers to equal a power of 3. Unique factorization rules this out. Thus their normalized squared \(L^2\) contributions add independently as \(\ell\downarrow0\), with **no unknown signs or cancellation**.

The exact finite census gives six internal translates for \(n\in\{2,4\}\), four for \(n\in\{5,7,8,11\}\), two for \(n\in\{13,16,17,19,23,25,29,31,32,37\}\), and none for \(n\ge41\). Let

\[
\begin{gathered}
A=\sum_{n\in\{2,4\}}\kappa_n^2,\qquad B=\sum_{n\in\{5,7,8,11\}}\kappa_n^2,\\
C=\sum_{n\in\{13,16,17,19,23,25,29,31,32,37\}}\kappa_n^2,\quad
\kappa_n=\frac{\log p}{\sqrt n}\ \text{for } n=p^k,\\
t_+=(s-1)^2,\quad t_-=(s+1)^2.
\end{gathered}
\]

Then the exact non-ternary contribution to the limiting **squared** residual norms is

\[
\boxed{C_\pm=\frac{1+2t_\pm}{1+t_\pm}A+B+\frac{1}{1+t_\pm}C.}
\tag{G.5}
\]

These are explicit finite sums involving only \(\log p\), rational weights and \(\sqrt3\). The separate exact rational-interval auditor (no floating point in its internal arithmetic) gives

\[
\boxed{\begin{array}{rcl}
3.38439635&\le C_+\le&3.38439636,\\
2.59455169&\le C_-\le&2.59455170,\\
0.21839262&\le V_{3,+}\le&0.21839263,\\
0.01126514&\le V_{3,-}\le&0.01126515.
\end{array}}
\tag{G.6}
\]

In particular, \(C_++V_{3,+}>18/5=3.6\) and \(C_++V_{3,+}>C_-+V_{3,-}\) rigorously. **This is a positive, *nonvanishing* arithmetic residual floor.** It contradicts any stronger inference that the *mixed representers themselves* must converge to zero as localization narrows; only their Schur penalty was shown to decay in T1D-F.

## 6. Sharp asymptotic theorem for projected residual and Schur penalty

**Theorem (smooth-core, fixed T1C profile):** In the exact setting of §2,

\[
\boxed{\begin{aligned}
\lim_{\ell\downarrow0}\|r_{+,\ell}\|_2^2&=V_v+V_{3,+}+C_+,\\
\lim_{\ell\downarrow0}\|r_{-,\ell}\|_2^2&=V_v+V_{3,-}+C_-,\\
\mu_\pm(\ell)&=\log(1/\ell)+a_v+\lambda_{3,\pm}+o(1),\\
\lim_{\ell\downarrow0}\log(1/\ell)\,\|\mathcal R_\ell\|
&=V_v+V_{3,+}+C_+\ >\frac{18}{5}.
\end{aligned}}
\tag{G.7}
\]

Thus \(\boxed{\|\mathcal R_\ell\|=\Theta(1/\log(1/\ell))}\). T1D-F's **upper** asymptotic decoupling remains correct, but its *rate* cannot be improved to \(o(1/\log(1/\ell))\) for this family. The nonzero residual is not evidence of a negative Weil direction; it is an unavoidable mixed interaction sourced by prime translations and the logarithmic Fourier shape variance.

**Analytic proof.** By the exact published digamma asymptotic (G.1), for each fixed \(\eta\ne0\),
\(m_\Gamma(\eta/\ell)-\log(1/\ell)\to\log|\eta|-\log(2\pi)\). The previously documented uniform estimate \(|m_\Gamma(\xi)-\log(1+|\xi|)|<40\) gives the dominating integrable Fourier weight \(41+|\log|\eta||\) for \(0<\ell\le1\). Since \(\int(1+|\log|\eta||)^2|\widehat v(\eta)|^2d\eta<\infty\), dominated convergence and Plancherel prove the **strong L² scaling limit** of \((\mathcal G-\log(1/\ell))U_{\ell,x}v\) to \(U_{\ell,x}Lv\) in the sense that the norm of the difference tends to zero.

The prime source is exactly the finite sum of 24 paired translations. For each fixed center inside the window, rescaled \(L^2\) profiles become asymptotically orthogonal to profiles centered at any other fixed position; translates whose centers lie outside the window vanish in \(L^2(I)\). Every ternary translate \(n=3,9,27\) lands at an original cell center, yielding \(P_3c_\sigma\); all non-ternary translates occupy mutually distinct off-base centers and contribute the separate positive sum (G.5). In the four base neighborhoods, the \(Lv\) part and ternary part yield the orthogonal decomposition \((Lv-a_vv)\) and \((P_3-\lambda_{3,\sigma}I)c_\sigma\), proving their squared contributions are \(V_v\) and \(V_{3,\sigma}\).

The residual projection onto the **two windowed pole directions** tends to zero: the renormalized image \((\mathcal K_2-\log(1/\ell))u_{\sigma,\ell}\) is a finite sum of bounded-norm rescaled \(L^2\) profiles at fixed centers and hence tends weakly to zero. Its pole projections therefore vanish. The removed \(\log(1/\ell)u_{\sigma,\ell}\) has **exactly zero** pole moments, so contributes none. Projection onto \(E_\ell\) removes the component \(\mu_\sigma(\ell)u_{\sigma,\ell}\) exactly by the earlier parity-diagonal Weil matrix. This establishes both residual limits and the Rayleigh asymptotics, with no assumption that \(\mathcal K_2\) is globally bounded on \(L^2\). Since (G.6) gives \(V_{3,+}+C_+>V_{3,-}+C_-\), the plus sector determines the limiting maximum in the rank-two penalty. **No identification of completed Suzuki/Connes operator domains has been used.**

## 7. The complement does not become a high-frequency subspace merely because E localizes

There is a separate **finite-codimension obstruction to a common shortcut**. Choose a nonzero real even nonnegative \(b\in C_c^\infty(-1/100,1/100)\), and put \(b_0(x)=b(x)\), \(b_1(x)=b(x-1/5)+b(x+1/5)\). Both are even and have disjoint supports in \((-1/4,1/4)\). Let

\[
\alpha=\frac{\int b_0(x)\cosh(x/2)dx}{\int b_1(x)\cosh(x/2)dx},\qquad y=b_0-\alpha b_1,
\]

using a nonnegative \(b\) to ensure a nonzero denominator. Then \(y\ne0\) is even and satisfies **both** Mellin pole moments exactly. It is physically disjoint from every four-cell selected source for all \(\ell\le1/16\): the nearest selected centers are \(\pm\tfrac12\log3\), and the nearest selected support begins beyond \(\tfrac12\log3-1/32>1/4\). Therefore

\[
\boxed{0\ne y\in F_\ell\quad\text{for every}\quad 0<\ell\le1/16.}
\tag{G.8}
\]

This is a **fixed**, smooth, finite logarithmic-Fourier-energy complement vector, independent of \(\ell\). Its correlation is supported in \((-1/2,1/2)\), so every prime shift \(\log n\ge\log2>1/2\) vanishes, and its own pole contribution is zero: \(q(y,y)\) is **Gamma-only and independent of \(\ell\)**. Consequently the complement's Rayleigh infimum satisfies the finite upper comparison \(\inf_{0\ne z\in F_\ell}q(z,z)/\|z\|^2\le q(y,y)/\|y\|^2<\infty\). Merely narrowing \(E_\ell\) therefore cannot force the whole complement to arbitrarily high logarithmic frequency or give a complementary coercivity margin diverging with \(\log(1/\ell)\). It **does not rule out** a finite positive complementary lower bound, which is exactly what remains missing.

## 8. Disposition, falsifiers, and next genuine certification

**ANALYTICAL PASS (source level):** the *sharp-rate* asymptotic identity (G.7), a strictly positive finite-arithmetic residual floor, parity reduction, and a fixed smooth low-frequency-energy vector that remains inside every narrow-source complement. **AUDIT-PASS:** exact 24-prime-power enumeration, all finite translated-center multiplicities, ternary variance identities, and certified rational enclosures (G.6). **OPEN:** the actual finite-\(\ell=10^{-4}\) residual integrals, certified value of \(V_v\), a **positive** lower bound for \(q|_{F_\ell}\), semilocal (H2′)/R1/R2 and Sonin-conditioning intertwiners, full Weil positivity, or RH.

**Falsifiers:** a discrepancy in the digamma Fourier normalization \(-\log2\pi\), a missed active prime power or center coincidence, invalid use of a noncompact \(Lv\) in the concentration limit, failure of the pole projections to vanish after subtracting \(\log1/\ell\), or an incorrect metric in \(\Pi_{E_\ell}\). The exact rational checker covers the combinatorial/finite cases, *not* the Fourier dominated-convergence or smooth-core proof.

**Next CR-INT-T1D-H:** evaluate \(V_v\) and the **frozen-width** two scalar \(\beta_\pm(10^{-4})^2\) by certified high-precision interval integration with explicit tail bounds. Independently test a plausible **positive complement coercivity bound on the actual \(F_\ell\)** or document a specific mismatch with the Connes conditioning constraints. The limit (G.7) is a precise cross-check and lower diagnostic for those calculations; it does not replace them.