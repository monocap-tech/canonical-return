# CR-INT-T1D-A — Profile-uniform positivity and the first non-ternary prime collision

**Canonical Return · analytic research note · November 22, 2026 (America/Los_Angeles).**  
**Status:** DERIVED / ANALYTIC, with a separately executed exact-rational algebra audit. **Not Lean-certified, not peer-reviewed, not global Weil positivity, and not RH.** This is a restricted smooth-core result; a general support-independent lower bound remains open. This research note follows [T1B](../../archive/space/09-t1b.md), [T1C](../../archive/space/10-t1c.md), and [Issue #1](https://github.com/monocap-tech/canonical-return/issues/1). The two computational source scripts in `src/audits/` check only the stated finite-dimensional identities and rational inequalities.

## 1. Exact source and carrier

Work with the original **signed** Weil functional in [Suzuki, arXiv:2606.09096v4, §1.1](https://arxiv.org/html/2606.09096v4); the comparison to Suzuki's generalized screw-kernel form uses his §8.4–8.5 smooth-core identities. The comparison to the semilocal Connes trace is *conditional on the separate mixed trace-source argument documented in [T1A](../../archive/space/08-t1a.md)*; it is **not** independently re-proved in this note.

Fix \(d=\log 3\), \(s=4/\sqrt3\), and the support radius \(a=2\). Let \(0<\ell\le 1/16\) and choose **any nonzero real** \(\phi\in C_c^\infty(\mathbb R)\) supported inside an interval of length at most \(\ell\). **No positivity, evenness, fixed shape, or normalization is required**; \(\phi\) may change sign. After centering the support interval, set

\[
\phi_j(x)=\phi(x-(j-\tfrac32)d),\quad j=0,1,2,3,
\qquad e_0=\phi_0-s\phi_1+\phi_2,\quad e_1=\phi_1-s\phi_2+\phi_3.
\]

All supports lie in \((-2,2)\). Both test functions satisfy the two exact Mellin pole-null conditions, since their coefficients are the shifts of \(q(z)=(z-\sqrt3)(z-1/\sqrt3)\). Thus \(E_\phi=\operatorname{span}_{\mathbb C}\{e_0,e_1\}\) is a two-dimensional pole-null smooth test space. Its **native physical mass matrix** is

\[
M=J\begin{pmatrix}22/3&-8\sqrt3/3\\-8\sqrt3/3&22/3\end{pmatrix},
\quad J=\|\phi\|_2^2>0.
\]

The determinant is \(292J^2/9>0\); this is a genuine Hermitian generalized pencil, not an ordinary eigenproblem for the numerator alone.

## 2. Exact arithmetic isolation for every profile of that width

Because \(e^x\le(1-x)^{-1}\) for \(0\le x<1\),

\[
 e^{1/16}\le16/15<29/27.
\]

For each candidate center \(3,9,27\), a prime-power \(n\) able to contribute at its shifted cross-correlation must satisfy \(n\in(15m/16,16m/15)\) with \(m=3,9,27\). These rational windows contain, respectively, only \(3\), only \(9\), and the integers \(26,27,28\). Since neither 26 nor 28 is a prime power, the **only possible** cross-prime-power shifts on these four bumps are \(3,9,27\). Same-bump prime shifts vanish because \(\ell<\log2\). This argument is profile-independent and does not require nonnegativity.

## 3. Profile-independent Gamma diagonal bound

Write \(R(t)=J^{-1}\int\phi(x)\phi(x-t)dx\), \(C_\Gamma=\log(4\pi)+\gamma\), and \(w(t)=e^{-t/2}/(1-e^{-2t})\). Since \(\phi\) is real, \(|R(t)|\le1\) by Cauchy–Schwarz, and \(R(t)=0\) for \(|t|\ge\ell\). The **full original Gamma diagonal** is

\[
\frac{g_0}{J}=-C_\Gamma+2\int_0^\infty(e^{-t/2}-R(t))w(t)dt.
\]

On \((0,\ell)\), \(2(1-e^{-t/2})w(t)\le1\). On \((\ell,\infty)\) the correlation vanishes, and

\[
2\int_\ell^\infty e^{-t/2}w(t)dt=\log\coth(\ell/2)\ge\log(2/\ell).
\]

Consequently, **for every nonzero real smooth bump of width at most \(\ell\)**,

\[
\boxed{\frac{g_0}{J}\ge\log(2/\ell)-C_\Gamma-\ell.}\tag{T1D.1}
\]

An explicit exact-rational proof of the loose constant \(C_\Gamma<16/5\) uses \(\pi<22/7\), the fifth-order positive Taylor partial sum \(\sum_{j=0}^{5}(13/5)^j/j!>88/7\) (so \(\log(4\pi)<13/5\)), and the standard harmonic estimate \(\gamma<H_4-\log4-1/10\) together with the atanh lower bound \(\log2>2(1/3+1/81+1/1215)=842/1215\), giving \(\gamma<2903/4860<3/5\). The same lower bound for \(\log2\) proves \(\log32>17/5\). Hence at \(\ell\le1/16\),

\[
\boxed{g_0/J>1/5-1/16=11/80.}\tag{T1D.2}
\]

These bounds do **not** use a bound on \(\|\phi'\|_2\), and are therefore uniform over highly oscillatory as well as nonnegative profiles.

## 4. Signed off-diagonal Gamma terms and uniform positive pencil

Define \(A_k=\iint\phi(x)\phi(y)w(kd+x-y)dxdy\), \(a_k=A_k/J\), for \(k=1,2,3\). The \(A_k\) need **not** be positive if \(\phi\) changes sign. Nevertheless

\[
|a_k|\le\frac{(\int|\phi|)^2}{J}\sup_{|t-kd|\le\ell}w(t)
\le\ell w(kd-\ell)<\ell
\]

because \(kd-\ell>\log2\) and \(w(t)\) is strictly decreasing with \(w(\log2)<1\). Thus the T1B generalized eigenvalue formulas apply with **\(|a_k|<\ell\)** and all arithmetic prime weights unchanged:

\[
\begin{aligned}
\mu_+&=g_0/J+\frac{41-5\sqrt3}{73}\log3
+\frac{-83+76\sqrt3}{146}a_1+\frac{15+32\sqrt3}{73}a_2-\frac{33+12\sqrt3}{146}a_3,\\
\mu_-&=g_0/J+\frac{41+5\sqrt3}{73}\log3
+\frac{83+76\sqrt3}{146}a_1+\frac{15-32\sqrt3}{73}a_2+\frac{33-12\sqrt3}{146}a_3.
\end{aligned}
\]

Use \(5/3<\sqrt3<7/4\) and \(\log3>1\). In the \(\mu_+\) row the total *absolute* Gamma coefficient is less than \(2\); in the \(\mu_-\) row it is less than \(3\). The prime offsets obey \((41-5\sqrt3)\log3/73>129/292\) and \((41+5\sqrt3)\log3/73>41/73\). Combining these with (T1D.1) gives the width-dependent bounds

\[
\boxed{\mu_+>\log(2/\ell)-16/5+129/292-3\ell,}\tag{T1D.3}
\]
\[
\boxed{\mu_->\log(2/\ell)-16/5+41/73-4\ell.}\tag{T1D.4}
\]

For \(\ell\le1/16\), the rough rational constants yield

\[
\mu_+>\frac15+\frac{129}{292}-\frac3{16}=\frac{2653}{5840}>\frac9{20},
\qquad \mu_->\frac15+\frac{41}{73}-\frac14=\frac{747}{1460}>\frac12.
\]

**Profile-uniform restricted positivity theorem:**

\[
\boxed{\forall\,0\ne f\in E_\phi,\quad W(f*\widetilde f)>\frac9{20}\|f\|_2^2.}\tag{T1D.5}
\]

This holds for **every nonzero real \(C_c^\infty\) bump of width at most \(1/16\)**, including sign-changing profiles. It strengthens T1C's single-profile interval result in *generality*, though its numerical margin is far smaller. It is still only a **selected 2D subspace** of the full test domain and supplies **no** global Weil-form lower bound.

## 5. Profile-uniform phase splitting, even without nonnegativity

The exact T1B splitting is

\[
\mu_--\mu_+=\frac{10\sqrt3\log3+83a_1-64\sqrt3a_2+33a_3}{73}.
\]

Using \(|a_k|<\ell\le1/16\), \(\sqrt3>5/3\), \(\sqrt3<7/4\), and \(\log3>1\),

\[
\boxed{\mu_--\mu_+>\frac{50/3-228/16}{73}=\frac{29}{876}>0.}\tag{T1D.6}
\]

For the normalized phases \(\beta=1,i,-1\) with \(f_\beta=(e_0+\beta e_1)/\|e_0+\beta e_1\|_2\), the complete signed Weil quadratic values **always** satisfy

\[
\boxed{W(f_1*\widetilde f_1)<W(f_i*\widetilde f_i)<W(f_{-1}*\widetilde f_{-1}).}\tag{T1D.7}
\]

Their **smaller adjacent** separation exceeds \(\frac{11-4\sqrt3}{22}(\mu_--\mu_+)>29/4818\). Every pole contribution is zero by exact annihilation; every Gamma and prime term is included. This bound is nontrivial for arbitrary real profiles, even ones with negative cross-Gamma integrals.

## 6. The first new arithmetic activation is **prime 29**, not a higher power of 3

Let \(\delta=\log(29/27)\). This is the **smallest positive distance** from any of \(\log3,\log9,\log27\) to another prime-power logarithm. In particular, \(\delta<\log(27/25)\) because \(29\cdot25=725<729=27^2\). Thus for bumps whose width is \(\ell\le\delta\), no extra arithmetic collision is possible. At \(\ell>\delta\), a nonnegative profile occupying nearly the full width can activate the \(n=29\) cross-term between the outermost two translated cells. The occurrence is **possible**, not compulsory for every narrower-than-window profile.

For \(\delta<\ell<\log(27/25)\), **29 is the only newly admissible prime power** beyond \(3,9,27\). Define the real autocorrelation at the shifted prime \(R_\phi(\delta)=J^{-1}\int\phi(x)\phi(x-\delta)dx\). Let \(\kappa_{29}=\log29/\sqrt{29}>0\). Relative to the no-29 Gram matrix, the actual full signed Weil 2D matrix changes *only* in its mixed entry by

\[
\boxed{\Delta H=-J\kappa_{29}R_\phi(\delta)
\begin{pmatrix}0&1\\1&0\end{pmatrix}.}\tag{T1D.8}
\]

The physical metric \(M\) does **not** change. Consequently

\[
\Delta\mu_+=-\frac{3\kappa_{29}R_\phi(\delta)}{22-8\sqrt3},
\qquad
\Delta\mu_- =+\frac{3\kappa_{29}R_\phi(\delta)}{22+8\sqrt3},
\]

and the selected generalized splitting changes by the exact amount

\[
\boxed{\Delta(\mu_--\mu_+)=\frac{33}{73}\kappa_{29}R_\phi(\delta).}\tag{T1D.9}
\]

For \(\phi\ge0\), \(R_\phi(\delta)\ge0\), so the **first newly activated genuine prime term further separates** the two selected generalized directions. For a sign-changing \(\phi\), this additional term can have either sign. **This formula is only valid before the next prime collision**; past \(\log(27/25)\) the 25 term must also be restored.

## 7. Pointwise continuity vs failure of uniform activation smoothness

For a *fixed* smooth \(\phi\), \(R_\phi(t)\) is continuous in \(t\), and its correlation vanishes when the supports do not overlap. However, there is **no uniform smallness** of the new prime coefficient over all test profiles as the admissible width approaches \(\delta\) from above.

For each \(0<\eta<\min(\delta/3,\log(27/25)-\delta)\), choose any nonzero real nonnegative \(b_\eta\in C_c^\infty(0,\eta)\) and set \(\phi_\eta(x)=b_\eta(x)+b_\eta(x-\delta)\). Its support lies within an interval of width \(\delta+\eta\), and the two components are disjoint. Direct expansion gives

\[
\boxed{R_{\phi_\eta}(\delta)=\frac12}\qquad\text{for every }\eta>0,
\]

while at or below the exact width threshold every allowable profile has zero correlation at \(\delta\). Thus the **supremum of activated normalized prime strength over changing profiles** has a nonzero one-sided jump, although each individual fixed smooth correlation behaves continuously. There is no contradiction with Suzuki's continuity theorem for the *lowest full* localized eigenvalue: the test functions are changing and become singularly localized.

Indeed, for the two-peak family and \(\eta<\delta/4\), its autocorrelation vanishes throughout \((\eta,\delta-\eta)\). On \((\eta,\delta/2)\), the positive Gamma diagonal integrand satisfies \(2e^{-t/2}w(t)=2e^{-t}/(1-e^{-2t})\ge e^{-\delta/2}/t\). The potentially negative correlation terms near \(0\) and near \(\delta\) contribute a bounded correction (at worst \(-\eta\) near zero and \(O(\eta)\) near \(\delta\)). Therefore

\[
\boxed{g_0(\phi_\eta)/\|\phi_\eta\|_2^2\ \longrightarrow+\infty\quad(\eta\downarrow0).}
\]

Meanwhile the normalized cross terms are uniformly bounded by their finite support width and the prime weights. Thus the selected generalized Rayleigh levels for this family grow rather than supply a low-energy counterexample as the new 29 interaction appears. **This is a statement about this explicit profile family, not a theorem bounding the entire Suzuki or Connes spectrum.**

## 8. Transfer conditions and disposition

**ANALYTICAL PASS, not Lean:** (i) uniform positive lower margin \(9/20\) for every real smooth profile of width ≤1/16 on its exact two-pole-null two-dimensional space, (ii) uniform full-form phase splitting >29/876, (iii) first additional arithmetic support threshold \(\log(29/27)\), (iv) exact rank-two perturbation of the selected Gram pencil and sign of its effect for nonnegative profiles, (v) a discontinuity of the *profile supremum* of normalized newly activated prime strength, accompanied by divergent localization cost for the explicit two-peak sequence.

**OPEN:** positivity for wider profiles after further prime collisions, any estimate for the infinite-dimensional orthogonal complement, global Suzuki/Connes spectral bottom, the semilocal H2′/R1/R2 or uniform Schur hypotheses, completed-domain intertwiner, and RH.

**Reproducibility:** `src/audits/canonical_return_t1d_uniform_profile.py` executes the finite symbolic Gram perturbation and all numerical rational comparison inequalities used in (T1D.5)–(T1D.9), and prints PASS. The integral/correlation arguments are human mathematical proofs in this note, not assertions independently established by the Python checker. Earlier T1B/T1C auditors remain applicable to their specific hypotheses.

**Next / Issue #1 remains OPEN:** quantify the next 25-prime collision and choose a broader explicit smooth family; independently test whether any source-uniform estimate can be connected to the Connes semilocal complement rather than repeatedly proving selected compression positivity.