# CR-INT-T1D-C — The 8 = 2³ collision: a two-entry perturbation and profile-uniform positivity to the prime-31 threshold

**Canonical Return · independent analytic research, November 22, 2026 PST.**  
**Status:** DERIVED / ANALYTIC with a separately executed exact-rational and SymPy audit. **Not Lean-certified, not externally reviewed, not a whole-Weil-form positivity theorem, not an operator equivalence, and not RH.**  
Source references: [Suzuki, *Weil's quadratic form via the screw function*, arXiv:2606.09096v4, §1.1 and §§8.4–8.5](https://arxiv.org/html/2606.09096v4), the preceding [T1D-A](CR-INT-T1D-A.md) and [T1D-B](CR-INT-T1D-B.md), and the [original T1B finite Gram derivation](../archive/space/09-t1b.md). The separate Connes trace-domain identification of [T1A](../archive/space/08-t1a.md) is **not re-proved** here. This note proves source-level facts only for the original complete signed Weil functional.

## 1. Fixed family and next collision census

Let \(d=\log3\), \(s=4/\sqrt3\), and fix source support radius \(a=2\). Choose any **nonzero real** \(\phi\in C_c^\infty(\mathbb R)\), possibly sign-changing and disconnected, whose support is contained in an interval of width

\[
0<\ell<\boxed{\delta_{31}:=\log(31/27)\approx0.13815033848}. 
\]

Center that interval; define \(\phi_j(x)=\phi(x-(j-3/2)d)\) for \(j=0,1,2,3\) and \(e_0=\phi_0-s\phi_1+\phi_2\), \(e_1=\phi_1-s\phi_2+\phi_3\). Then \(E_\phi=\mathrm{span}_{\mathbb C}\{e_0,e_1\}\) has both Mellin pole moments **zero**, because its coefficient columns are multiples of \((z-\sqrt3)(z-1/\sqrt3)\). The four cell supports are disjoint and lie in \((-2,2)\). With \(J=\|\phi\|_2^2>0\), the **unmodified physical Gram denominator** is

\[
\boxed{M=J\begin{pmatrix}22/3&-8\sqrt3/3\\-8\sqrt3/3&22/3\end{pmatrix}>0.}
\]

At spacing \(d\), the only prime power in the corresponding logarithmic band is \(3\); at spacing \(2d\) the possible prime powers are **\(8,9\)**; at spacing \(3d\), they are **\(25,27,29\)**. These are exhaustive for \(\ell<\delta_{31}\): the strict rational integer windows are

\[
\begin{array}{c|c|c}
\text{spacing center}&\text{open integer interval (conservative)}&\text{prime powers}\\\hline
3&(81/31,31/9)&3\\
9&(243/31,31/3)&8,9\\
27&(729/31,31)&25,27,29.
\end{array}
\]

The **next** collision beyond \(8\) is \(31\), at the excluded endpoint \(\delta_{31}\). Since \(31/27<(29/27)^2\) and \(31/27<(9/8)^2\), all three new shifts \(\delta_{29}=\log(29/27),\delta_{25}=\log(27/25),\delta_8=\log(9/8)\) satisfy \(\delta>\ell/2\) whenever they are active, so their normalized source autocorrelations obey \(|R_\phi(\delta)|\le1/2\).

## 2. Exact 8-perturbation: *both* numerator entries change

Let \(R_\phi(t)=J^{-1}\int\phi(x)\phi(x-t)\,dx\) and \(\kappa_8=\Lambda(8)/\sqrt8=\log2/(2\sqrt2)\), and retain the B-note outer-cell amplitude \(b_\phi=\kappa_{29}R_\phi(\delta_{29})+\kappa_{25}R_\phi(\delta_{25})\), where \(\kappa_{29}=\log29/\sqrt{29}\), \(\kappa_{25}=\log5/5\). Relative to the complete old T1B numerator \(H_0\), the **complete original Weil numerator** on \(E_\phi\) is

\[
\boxed{H=H_0-Jb_\phi\begin{pmatrix}0&1\\1&0\end{pmatrix}
-J\kappa_8R_\phi(\delta_8)\begin{pmatrix}2&-8\sqrt3/3\\-8\sqrt3/3&2\end{pmatrix}.}\tag{C.1}
\]

The new matrix follows by compressing the genuine \(8\)-prime-power convolution matrix \(D_2\), whose only nonzero entries are \((D_2)_{0,2}=(D_2)_{2,0}=(D_2)_{1,3}=(D_2)_{3,1}=1\), with the **same** four-to-two Mellin-null embedding \(X=\begin{psmallmatrix}1&0\\-s&1\\1&-s\\0&1\end{psmallmatrix}\):

\[
\boxed{X^TD_2X=\begin{pmatrix}2&-8\sqrt3/3\\-8\sqrt3/3&2\end{pmatrix}.}\tag{C.2}
\]

This matrix is indefinite (one eigenvalue \(2+8\sqrt3/3>0\), the other \(2-8\sqrt3/3<0\)). In particular, it is **not** the simple mixed-only perturbation of the 25/29 collisions. The physical denominator \(M=JX^TX\) and the exact two-pole nullspace do not change.

Because the pencil is simultaneously diagonal in the \((1,1)\) and \((1,-1)\) coefficient directions, the **exact new shift** is

\[
\boxed{\begin{aligned}
\Delta_8\mu_+&=\frac{15+32\sqrt3}{73}\,\kappa_8R_\phi(\delta_8),\\
\Delta_8\mu_-&=\frac{15-32\sqrt3}{73}\,\kappa_8R_\phi(\delta_8),\\
\Delta_8(\mu_- -\mu_+)&=-\frac{64\sqrt3}{73}\,\kappa_8R_\phi(\delta_8).
\end{aligned}}\tag{C.3}
\]

For a **nonnegative** profile, the \(8\)-term increases the symmetric Rayleigh value and decreases the antisymmetric one, so it **reduces** their separation. This explicitly contrasts with the 25 and 29 terms. For a sign-changing profile the signs can reverse.

## 3. A layered autocorrelation lemma, stronger than the T1D-B half-support estimate

Let \(\phi\in L^2\) be real and supported in an interval of width at most \(\ell\). For integers \(j\ge1\), for almost every \(t\in(\ell/(j+1),\ell/j)\),

\[
\boxed{|R_\phi(t)|\le \cos\frac{\pi}{j+2}.}\tag{C.4}
\]

**Proof:** Decompose the original support into residue fibers \(x\bmod t\), and consider the ordered values \(\phi(u+kt)\) on each fiber. The number of nonzero positions is at most \(j+1\) in the stated interval. The fiber autocorrelation is the quadratic form of the symmetric tridiagonal path matrix with \(1/2\) on its two neighboring diagonals; its spectral radius is \(\cos(\pi/(N+1))\le\cos(\pi/(j+2))\). Integrate over \(u\in[0,t)\) and normalize by \(J\). This elementary path bound is profile-independent and also applies to sign-changing functions.

The \(j=1\) case recovers \(|R(t)|\le1/2\) for \(t\in[\ell/2,\ell]\). The \(j=2,3\) cases improve the bound over \((\ell/3,\ell/2)\) and \((\ell/4,\ell/3)\). For \(j\ge4\), use \(\pi>3\) and the alternating cosine estimate \(1-\cos x\ge x^2/2-x^4/24\) with \(x=3/(j+2)\). No numerical trigonometric evaluation is needed.

## 4. New **profile-uniform** complete Gamma diagonal lower bound

As in T1D-A/B, write \(g_0=W_\Gamma(\phi*\widetilde\phi)\), \(C_\Gamma=\log(4\pi)+\gamma\), \(w(t)=e^{-t/2}/(1-e^{-2t})\). The actual Suzuki signed-Weil formula gives

\[
\frac{g_0}{J}=-C_\Gamma+2\int_0^\infty[e^{-t/2}-R_\phi(t)]w(t)\,dt.
\]

The tail \(t\ge\ell\) contributes exactly \(\log\coth(\ell/2)\ge\log(2/\ell)\). For the interior, \(2(e^{-t/2}-1)w(t)\ge-1\) and \(2w(t)\ge(1-\ell/2)/t\). Applying (C.4) on the eight layers \((\ell/(j+1),\ell/j)\), \(j=1,\dots,8\), proves

\[
\boxed{\frac{g_0}{J}\ge\log(2/\ell)-C_\Gamma-\ell+(1-\ell/2)
\sum_{j=1}^8\!\left(1-\cos\frac{\pi}{j+2}\right)\log\frac{j+1}{j}.}\tag{C.5}
\]

This is a genuinely stronger *profile-independent* localization-energy estimate than the single half-support layer in T1D-B. It is a consequence of the simple finite-path spectral bound, not a new full spectral realization of the Weil form. The new argument uses only the width, not derivatives, nonnegativity, or a prescribed profile shape.

## 5. Exact rational bounds and a wider restricted positivity theorem

Set \(L=7/50\). The elementary Taylor inequality \(e^L>1+L+L^2/2>31/27\) gives \(\ell<\delta_{31}<L\). From the rational proof in T1D-B, \(C_\Gamma<25/8\). The audit uses **exact Fraction computations only** for the following rational enclosures: logarithms by an atanh series with geometric remainder; \(\sqrt{2},\sqrt3,\sqrt5,\sqrt8,\sqrt{29}\) by integer-square-root brackets; the \(j=1,2,3\) layer gains exactly in radicals and the \(j\ge4\) gains by \(\tfrac12x^2-\tfrac1{24}x^4\) at \(x=3/(j+2)\). Therefore (C.5) gives a rigorous bound

\[
\boxed{g_0/J>-0.059463.}\tag{C.6}
\]

For the Gamma cross-correlations \(A_k\) of T1B, \(|A_k|/J\le\ell w(kd-\ell)\). Using \(e^{\ell/2}\le(1-L/2)^{-1}\), \(e^{2\ell}\le(1-2L)^{-1}\), the weights are bounded by the **rational/radical** expressions

\[
w(kd-\ell)\le\frac{3^{-k/2}}{1-L/2}\left(1-\frac{3^{-2k}}{1-2L}\right)^{-1},\quad k=1,2,3.
\]

All three *new* amplitude correlations satisfy \(|R_\phi(\delta)|\le1/2\). Substituting those bounds and the exact coefficients of T1B together with (C.3) produces **certified rational lower bounds**

\[
\boxed{\mu_+>0.04037538>\frac1{40},\qquad
\mu_->0.39863077>\frac14.}\tag{C.7}
\]

Thus the **profile-uniform extension across the first spacing-two prime collision** is

\[
\boxed{\forall\,\phi\in C_c^\infty(\mathbb R,\mathbb R)\setminus\{0\}\ \text{with width}(\phi)<\log(31/27),\quad
\forall\,0\ne f\in E_\phi:\quad W(f*\widetilde f)>\frac1{40}\|f\|_2^2.}\tag{C.8}
\]

All Mellin pole terms are zero **by the exact two-row annihilator**; every admissible 3/9/27/29/25/8 prime-power contribution and the complete signed Gamma term is included. A **selected two-dimensional compression** has been bounded, *not* the infinite-dimensional complement.

## 6. What remains of the phase ordering?

For \(\phi\ge0\), \(A_1,A_2,A_3\ge0\) and \(A_2<A_1/\sqrt3\), so the original T1B Gamma correction *increases* \(\mu_- -\mu_+\). Also \(b_\phi\ge0\), so 29 and 25 increase that difference, whereas 8 decreases it. Because \(R_\phi(\delta_8)\le1/2\), the remaining **nonnegative-profile** splitting satisfies

\[
\mu_- -\mu_+ > \frac{\sqrt3}{73}\big(10\log3-32\kappa_8\big)>\boxed{\frac1{16}}.\tag{C.9}
\]

The normalized \(\beta=1,i,-1\) tests therefore retain the order \(Q_W(f_1)<Q_W(f_i)<Q_W(f_{-1})\), with the smaller adjacent gap \(>1/88\). **For arbitrary sign-changing profiles**, (C.8) holds but this phase ordering is **not proved uniform** here: the Gamma and new prime autocorrelations can change sign, and their separate worst-case bounds do not control the gap's sign.

## 7. Exact 8-only activation, and why new primes are not monotone

Let \(0<\eta<\min\{\delta_{31}-\delta_8,\delta_8-\delta_{25},\delta_8-\delta_{29},\delta_8/3\}\). Choose a nonzero nonnegative smooth \(b_\eta\) supported inside \((0,\eta)\), and set \(\phi_\eta(x)=b_\eta(x)+b_\eta(x-\delta_8)\). Then the enclosing support has width below \(\delta_{31}\) and the two components are disjoint. Direct correlation expansion gives

\[
\boxed{R_{\phi_\eta}(\delta_8)=\tfrac12,\quad
R_{\phi_\eta}(\delta_{29})=R_{\phi_\eta}(\delta_{25})=0.}\tag{C.10}
\]

So the first spacing-two interaction can genuinely activate *without activating the earlier 25/29 cross-couplings*. Its exact selected eigenvalue changes, **relative to the no-8 matrix evaluated on that same source**, are

\[
\boxed{\Delta_8\mu_+=\frac{15+32\sqrt3}{146}\kappa_8>0,\quad
\Delta_8\mu_-=\frac{15-32\sqrt3}{146}\kappa_8<0,\quad
\Delta_8(\mu_- -\mu_+)=-\frac{32\sqrt3}{73}\kappa_8<0.}\tag{C.11}
\]

The **effect on the spectral splitting reverses** compared to the 25/29 collisions; this is an exact local obstruction to a naive principle that every added prime must widen the selected gap. The short-scale Gamma self energy of this two-peak family again grows as \(\eta\downarrow0\), so its new arithmetic activation cannot by itself be read as a negative full Weil direction.

## 8. Transfer/status and next noncircular target

**ANALYTICAL PASS + exact rational/SymPy AUDIT-PASS:** the unique next \(8=2^3\) collision, correct two-entry rank-two perturbation, a layered source autocorrelation energy estimate, **all-real-profile strictly positive selected 2D pencil through \(\ell<\log(31/27)\)**, and retained phase ordering for nonnegative profiles. The audit verifies algebra and rational comparisons; the source-level integral and fiber-path estimates are mathematical arguments in this note, **not Lean elaborations**.

**OPEN:** sign of the phase gap for every sign-changing source in this wider class; what happens at the **next 31 prime** threshold; any global support-independent lower bound, invariant complement estimate, (H2′)/R1/R2 or Schur-complement realization, operator-domain identification, or RH. There is no claim that our selected Rayleigh values are actual eigenvalues of the infinite-dimensional Suzuki or Connes operator.

**Next — CR-INT-T1D-D:** analyze the 31 collision at spacing \(3\log3\), and decide whether additional layered autocorrelation control and localized source families can maintain a meaningful positive 2D margin. In parallel, formulate a **source-defined projection and cross-block inequality** needed to interface this restricted compression with the semilocal Schur complement; do not call a positive selected subspace a positive whole spectrum.

**Reproduce:** `python3 src/audits/canonical_return_t1d_c_8_collision.py` (SymPy and pure rational arithmetic). The PR is kept draft until separate review.