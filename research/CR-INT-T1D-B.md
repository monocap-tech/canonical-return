# CR-INT-T1D-B — The 25=5² collision and profile-uniform positivity through the 8=2³ threshold

**Canonical Return · independent research note · November 22, 2026, America/Los_Angeles.**  
**Status: DERIVED / ANALYTIC + exact-rational and symbolic AUDIT-PASS.** The finite source calculation is not Lean-certified or peer-reviewed, and **does not** establish positivity of the full localized Weil form, the Connes/Suzuki infinite-dimensional complement, or RH. This note extends [CR-INT-T1D-A](CR-INT-T1D-A.md), which remains on the same draft research branch, and uses the actual original signed Weil functional in [Suzuki, arXiv:2606.09096v4, §1.1](https://arxiv.org/html/2606.09096v4). The inherited semilocal trace comparison is subject to the separately recorded [T1A domain argument](../archive/space/08-t1a.md), **not** re-proved here.

## 1. Frozen smooth-source family and exact admissibility

Put \(d=\log3\), \(s=4/\sqrt3\), and fix support radius \(a=2\). Choose any **nonzero real** \(\phi\in C_c^\infty(\mathbb R)\) supported in an interval of length at most \(0<\ell<\log(9/8)\); **sign changes, asymmetry and multiple support components are allowed**. Center the enclosing interval and set

\[
\phi_j(x)=\phi(x-(j-\tfrac32)d),\quad j=0,1,2,3,\qquad
 e_0=\phi_0-s\phi_1+\phi_2,\quad e_1=\phi_1-s\phi_2+\phi_3.
\]

All supports lie strictly inside \((-2,2)\); the four cells are disjoint because \(\ell<\log(9/8)<d\). The two Mellin pole moments of \(e_0,e_1\) vanish **exactly**, since \(q(z)=1-sz+z^2=(z-\sqrt3)(z-1/\sqrt3)\). Write \(E_\phi=\operatorname{span}_\mathbb C\{e_0,e_1\}\) and \(J=\|\phi\|_2^2>0\). Its physical mass is the unchanged, **nonorthogonal** matrix

\[
\boxed{M=J\begin{pmatrix}22/3&-8\sqrt3/3\\-8\sqrt3/3&22/3\end{pmatrix}>0.}
\]

This is a *fixed-dimensional, source-varying* statement. It is not a claim about all functions with this support radius.

## 2. Complete collision census up to the next threshold

The three cross-cell spacings are \(d,2d,3d\), corresponding to \(3,9,27\). An integer prime power \(n\) can contribute at a given spacing only when \(|\log n-kd|<\ell\) (at equality a compact smooth correlation is zero). Define

\[
\boxed{\delta_{29}=\log(29/27),\qquad
\delta_{25}=\log(27/25),\qquad
\delta_8=\log(9/8).}
\]

The exact order is \(\delta_{29}<\delta_{25}<\delta_8\) (equivalently \(29/27<27/25<9/8\)). For \(0<\ell<\delta_8\), inspection of the strict integer windows for each of \(3,9,27\) shows **all and only** the possible prime-power interactions are

\[
\boxed{n\in\{3,9,27,29,25\}.}
\]

In particular, **25 is the second newly admissible term** after prime 29; it is the second iterate \(5^2\), with genuine coefficient \(\Lambda(25)/\sqrt{25}=\log5/5\), **not** \(\log25/5\). The **next** new interaction is \(8=2^3\) at width \(\delta_8\). An *admissible* shift need not be active for a given bump: its coefficient depends on actual autocorrelation.

## 3. Exact additive 25/29 update of the same Hermitian pencil

Define the normalized real even autocorrelation

\[
R_\phi(t)=\frac{1}{J}\int_{\mathbb R}\phi(x)\phi(x-t)\,dx,
\]

and the two exact new arithmetic amplitudes

\[
\kappa_{29}=\frac{\log29}{\sqrt{29}},\qquad
\kappa_{25}=\frac{\log5}{5},\qquad
b_\phi=\kappa_{29}R_\phi(\delta_{29})+\kappa_{25}R_\phi(\delta_{25}).
\]

Relative to the original T1B expression \(H_0=\begin{psmallmatrix}Q&C\\C&Q\end{psmallmatrix}\), where \(Q,C\) retain **the full archimedean Gamma part and the genuine 3,9,27 weights**, the exact updated matrix is

\[
\boxed{H_\phi=H_0-Jb_\phi\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad M_\phi=M.}\tag{B.1}
\]

The two new integers affect only the outermost cells, hence only the **mixed** entry of the two-dimensional pole-null Gram numerator; the pole term remains zero on \(E_\phi\). Both H and M must be carried into the generalized spectral problem. In the symmetric/antisymmetric directions, the new generalized eigenvalues satisfy

\[
\boxed{\begin{aligned}
\mu_+(\phi)&=\mu_+^{(0)}(\phi)-\frac{3b_\phi}{22-8\sqrt3},\\
\mu_-(\phi)&=\mu_-^{(0)}(\phi)+\frac{3b_\phi}{22+8\sqrt3},\\
\bigl(\mu_--\mu_+\bigr)(\phi)&=\bigl(\mu_-^{(0)}-\mu_+^{(0)}\bigr)(\phi)+\frac{33}{73}b_\phi.
\end{aligned}}\tag{B.2}
\]

For **nonnegative** bumps, \(R_\phi(\delta_{25}),R_\phi(\delta_{29})\ge0\), so their added prime terms **reinforce the selected eigenvalue separation**, even though they lower the *symmetric* Rayleigh quotient. For sign-changing bumps, the new perturbation may have either sign. In the nonnegative case, T1B's exact Gamma comparison \(A_2<A_1/\sqrt3\) persists for these widths, and hence \(\mu_--\mu_+>10\sqrt3\log3/73\). We make **no corresponding gap-order claim for all sign-changing profiles** at this larger scale.

## 4. New width-uniform autocorrelation inequality

For any real \(\phi\) supported in an interval \([u,u+\ell]\), with \(J=\|\phi\|_2^2\), Cauchy–Schwarz gives \(|R_\phi(t)|\le1\). But for \(\ell/2\le t\le\ell\), the source domain \([u,u+\ell-t]\) and translated target domain \([u+t,u+\ell]\) are **disjoint**. Their squared \(L^2\) masses sum to at most \(J\), so the arithmetic–geometric mean inequality yields the stronger estimate

\[
\boxed{|R_\phi(t)|\le\tfrac12\qquad(\ell/2\le |t|\le\ell).}\tag{B.3}
\]

Since \(\delta_8<2\delta_{29}\) (exactly \(9/8<(29/27)^2\)), throughout \(\ell<\delta_8\) we have \(\delta_{29},\delta_{25}>\ell/2\). Therefore both new normalized amplitudes are at most \(1/2\) in modulus. Moreover the elementary rational inequalities \(\kappa_{29}<16/25\) and \(\kappa_{25}<1/3\) imply

\[
\boxed{|b_\phi|<\tfrac12.}\tag{B.4}
\]

This is a **profile-independent** bound, including disconnected and oscillatory real sources; no derivative norm is assumed.

## 5. A new Gamma diagonal lower bound with a half-support gain

The original signed archimedean Gamma self-term has the exact expression

\[
\frac{g_0}{J}=-(\log(4\pi)+\gamma)+2\int_0^\infty\bigl(e^{-t/2}-R_\phi(t)\bigr)w(t)dt,
\qquad w(t)=\frac{e^{-t/2}}{1-e^{-2t}}.
\]

As in T1D-A, the exterior \(t\ge\ell\) tail is \(\log\coth(\ell/2)\ge\log(2/\ell)\), and \(2(1-e^{-t/2})w(t)\le1\). But now (B.3) gives \(1-R_\phi(t)\ge1/2\) for \(t\in[\ell/2,\ell]\). The positive interior correction therefore obeys

\[
2\int_0^\ell(1-R_\phi(t))w(t)dt
\ge\int_{\ell/2}^\ell w(t)dt
\ge \frac{1-\ell/2}{2}\log2,
\]

using \(w(t)\ge(1-t/2)/(2t)\). Consequently, **without assuming \(\phi\ge0\) or a derivative estimate**,

\[
\boxed{\frac{g_0}{J}\ge\log(2/\ell)-(\log(4\pi)+\gamma)-\ell+
\frac{1-\ell/2}{2}\log2.}\tag{B.5}
\]

This strengthens the diagonal localization inequality used in T1D-A and is the essential mechanism permitting a wider theorem despite the two added prime-power interactions.

## 6. Rationally certified positivity through the 8 threshold

For \(0<\ell<\delta_8<3/25\), exact elementary inequalities establish

\[
\log(2/\ell)>14/5,\quad \log(4\pi)+\gamma<25/8,\quad \log2>69/100,\quad
\log3>109/100,\quad 5/3<\sqrt3<7/4.
\]

For transparent reproduction: \(\log(50/3)=4\log2+\log(25/24)>4(69/100)+1/25=14/5\); the bound on \(\log(4\pi)\) follows from \(\pi<22/7\) and the positive seventh exponential Taylor partial sum at \(127/50\); \(\gamma<29/50\) follows from \(H_{20}-\log20-1/42\) and a rational atanh lower bound. The accompanying exact-fraction audit independently checks **all** numerical comparisons, including \(\kappa_{29}<16/25\), \(\kappa_{25}<1/3\).

Eq. (B.5) now gives the **universal normalized diagonal** estimate

\[
\boxed{g_0/J>-1207/10000.}\tag{B.6}
\]

For \(k=1,2,3\), retain the complete signed Gamma cross-integrals \(A_k\). Even for sign-changing bumps they obey \(|A_k|/J\le\ell w(kd-\ell)\). Elementary bounds on \(e^{\ell/2}\) and \(e^{2\ell}\) yield

\[
w(d-\ell)<3/4,\qquad w(2d-\ell)<3/8,\qquad w(3d-\ell)<1/4.
\]

Propagating these into the **exact** T1B generalized eigenvalue formulas, with both new prime amplitudes bounded by (B.4), gives independently checked rational inequalities

\[
\begin{aligned}
\mu_+&> -\frac{1207}{10000}+
\frac{14061}{29200}-\frac{417}{584}\frac3{25}-\frac3{16}
=\frac{63989}{730000}>\frac2{25},\\
\mu_-&> -\frac{1207}{10000}+
\frac{148}{219}\frac{109}{100}-\frac{98}{73}\frac3{25}-\frac9{212}
=\frac{47864051}{116070000}>\frac25.
\end{aligned}\tag{B.7}
\]

**Profile-uniform restricted positivity, now across both the 29 and 25 thresholds:**

\[
\boxed{\forall\,0\ne f\in E_\phi:\qquad
W(f*\widetilde f)>\frac2{25}\|f\|_2^2,
\quad 0<\mathrm{width}(\phi)<\log(9/8).}\tag{B.8}
\]

The inequality holds for *every* nonzero real \(C_c^\infty\) profile satisfying the width restriction, including sign-changing functions. It is a **two-dimensional restriction only**, not a lower bound on the full Suzuki/Connes window. The bound is weaker than the narrower T1D-A \(9/20\) margin because the class of sources is substantially larger.

## 7. Isolating the 25 collision even after 29 became admissible

Let \(\eta>0\) satisfy

\[
\eta<\min\bigl\{\delta_{25}-\delta_{29},\;\delta_{29}/4,\;\delta_8-\delta_{25}\bigr\}.
\]

Choose any nonzero nonnegative \(b_\eta\in C_c^\infty(0,\eta)\) and set \(\phi_\eta(x)=b_\eta(x)+b_\eta(x-\delta_{25})\). The two component supports are disjoint and its enclosing width is less than \(\delta_{25}+\eta<\delta_8\). Direct evaluation of normalized autocorrelations gives

\[
\boxed{R_{\phi_\eta}(\delta_{25})=\tfrac12,
\qquad R_{\phi_\eta}(\delta_{29})=0.}\tag{B.9}
\]

Thus the **second** collision can be activated *in isolation* even though the nominal 29 threshold was crossed. The 25 term alone induces the exact finite spectral changes

\[
\boxed{\Delta\mu_+=-\frac{3\log5}{10(22-8\sqrt3)},\quad
\Delta\mu_- =+\frac{3\log5}{10(22+8\sqrt3)},\quad
\Delta(\mu_--\mu_+)=\frac{33\log5}{730}>0.}\tag{B.10}
\]

The normalized **supremum** of the newly activated 25 amplitude need not go to zero as width decreases to its exact threshold from above: in this changing-profile family it is always \(\kappa_{25}/2\). This does **not** contradict continuity for each **fixed** test function: as \(\eta\downarrow0\) the profiles concentrate. On this sequence the Gamma diagonal localization cost diverges at least logarithmically (because the autocorrelation vanishes for \(\eta<t<\delta_{25}/2\) and the positive kernel behaves like \(1/(2t)\)); thus the finite phase perturbation does not create a low-Rayleigh counterexample here.

## 8. What the transfer **does not** establish

The trace representation requires the separate complex source-domain analysis of T1A; this note evaluates the original **signed** Weil form and its native physical Gram pencil. It does not compute a Connes semilocal operator trace independently. The inequality (B.8) concerns **only the selected two-dimensional pole-null subspace**, not its infinite-dimensional orthogonal complement, the Suzuki global ground state or the Connes semilocal Schur estimate. Neither the semilocal hypotheses H2′/R1/R2 nor original RH STRICT/BAL are advanced.

**Next, CR-INT-T1D-C:** analyze the \(8=2^3\) interaction at \(\ell=\log(9/8)\), which enters at spacing \(2d=\log9\) instead of \(3d\), and changes the **diagonal as well as mixed** entries upon compression. Do not reuse the simple outer-cell perturbation matrix past that threshold. Independently examine whether controlling the profile cost in a *growing* source family has any relation to the missing infinite-dimensional complement inequality.

## 9. Exact validation and custody

**Audit script:** `src/audits/canonical_return_t1d_b_25_collision.py`, executed with Python 3 and SymPy. It verifies the rational inequalities, threshold order, Gamma coefficient bounds, exact projection of the outer-cell update through the Mellin annihilator, and the generalized-eigenvalue formulas. The **support-geometry and Gamma integral estimates are mathematical proofs above, not runtime-evaluated integral formulas**. Report: **PASS**, with rational endpoint margins \(63989/730000\) and \(47864051/116070000\). No Lean proof or external referee certification is claimed.