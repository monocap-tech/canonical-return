# CR-INT-T1D-D — The prime-31 activation, outer-cell phase robustness, and the exact complement Schur gate

**Canonical Return · November 22, 2026 PST / November 23 UTC.**  
**Status: DERIVED / ANALYTIC, with executable exact SymPy and rational arithmetic audit. Not Lean, not peer reviewed, no completed-domain operator transfer, no full Weil positivity, no RH.**  
**Reference basis:** Suzuki, *Weil's quadratic form via the screw function*, [arXiv:2606.09096v4](https://arxiv.org/html/2606.09096v4), §1.1 and smooth-core formulas (8.7)–(8.8); Andrews–Connes–La Gatta, *The uncertainty operator, the Weil ground state and the zeros of zeta*, [arXiv:2611.13129v1](https://arxiv.org/html/2611.13129), §8.6 Lemma 8.16, Theorem 8.14 and Proposition 8.18, **whose additional semilocal hypotheses are not verified here**. Prior Canonical Return [T1B](../archive/space/09-t1b.md), [T1C](../archive/space/10-t1c.md), [T1D-A](CR-INT-T1D-A.md), [T1D-B](CR-INT-T1D-B.md), and [T1D-C](CR-INT-T1D-C.md). This note acts only in the **original signed Weil smooth-test source**, unless a cited canonical identity is specifically invoked.

## 1. Freeze the original source and exact pole conditions

Fix support radius \(a=2\), \(d=\log3\), \(s=4/\sqrt3\). For any nonzero **real** \(\phi\in C_c^\infty(\mathbb R)\), possibly sign-changing or disconnected, suppose its support lies in an interval of enclosing width \(0<\ell<\log(11/9)\). Center the interval and take four translates \(\phi_j(x)=\phi(x-(j-3/2)d)\), \(j=0,1,2,3\), which fit inside \((-2,2)\) and have mutually disjoint cell supports. Define \(e_0=\phi_0-s\phi_1+\phi_2\) and \(e_1=\phi_1-s\phi_2+\phi_3\), \(E_\phi=\operatorname{span}_{\mathbb C}\{e_0,e_1\}\). Each \(e_j\) has both Mellin pole moments \(\int e_j(x)e^{\pm x/2}\,dx=0\), because the column polynomial has factors \((z-\sqrt3)(z-1/\sqrt3)\).

The **native physical mass** \(M\) and pole-annihilator remain independent of width and prime collision:

\[
X=\begin{pmatrix}1&0\\-4/\sqrt3&1\\1&-4/\sqrt3\\0&1\end{pmatrix},\quad
L=\begin{pmatrix}(\sqrt3)^j\\(\sqrt3)^{-j}\end{pmatrix}_{j=0}^3,\quad LX=0,\qquad
M=JX^*X=J\begin{pmatrix}22/3&-8\sqrt3/3\\-8\sqrt3/3&22/3\end{pmatrix},
\]
where \(J=\|\phi\|_2^2>0\). Both \(H\) and \(M\) have the real symmetric \(\begin{psmallmatrix}q&b\\b&q\end{psmallmatrix}\) pattern for all extra collisions in this note, hence are diagonal in the **physical orthonormal** directions \(u_\pm=(e_0\pm e_1)/\|e_0\pm e_1\|_2\). Denote their generalized Rayleigh values by \(\mu_\pm=Q_W(u_\pm)\). **Their sign beyond the T1D-C width is not assumed.**

## 2. Full arithmetic collision census, not a guessed continuation

At separation \(kd\), a prime power \(n\) can couple two cells only when \(|\log(n/3^k)|<\ell\). For \(0<\ell<\log(11/9)\), this forces the *exact* rational windows \(9\cdot3^k/11<n<11\cdot3^k/9\). Enumerating prime powers yields:

| Cell separation | All possibly active prime powers for \(\ell<\log(11/9)\) |
|---|---|
| \(d=\log3\) | \(3\) |
| \(2d=\log9\) | \(8,9\) |
| \(3d=\log27\) | \(23,25,27,29,31,32\) |

The genuine new collision order, with exact thresholds, is
\[
\boxed{\begin{gathered}
29:\log(29/27),\quad25:\log(27/25),\quad8:\log(9/8),\\
31:\log(31/27),\quad\underbrace{23:\log(27/23)}_{\text{next after 31}},\quad32:\log(32/27),\quad11:\log(11/9).
\end{gathered}}
\]

In particular, **the next after 31 is prime 23, then the prime power 32, and only afterward prime 11 enters at spacing \(2d\)**. The order is proved by *exact rational comparisons* of the ratios and the exhaustive windows; it is not a decimal approximation or a heuristic list. The collision thresholds are *possible activations*, not guaranteed nonzero correlations for each profile.

## 3. Exact prime-31 perturbation (and its 23/32 successors)

Let \(R_\phi(t)=J^{-1}\int\phi(x)\phi(x-t)\,dx\) (real and even), and put \(\delta_{31}=\log(31/27)\), \(\kappa_{31}=\Lambda(31)/\sqrt{31}=\log31/\sqrt{31}\). A prime at an outermost-cell shift acts on the four-cell coefficients by \(D_3=E_{0,3}+E_{3,0}\), so
\[
\boxed{X^TD_3X=\begin{pmatrix}0&1\\1&0\end{pmatrix}=:B_3.}
\]

Consequently, **from the same complete signed Weil source**, the additional prime-31 term is exactly
\[
\boxed{\Delta_{31}H=-J\kappa_{31}R_\phi(\delta_{31})B_3,\qquad\Delta_{31}M=0.}
\]

After diagonalizing the *generalized* (not ordinary) pencil, the two exact shifts are
\[
\boxed{\begin{aligned}
\Delta_{31}\mu_+&=-\frac{3}{22-8\sqrt3}\,\kappa_{31}R_\phi(\delta_{31}),\\
\Delta_{31}\mu_-&=+\frac{3}{22+8\sqrt3}\,\kappa_{31}R_\phi(\delta_{31}),\\
\Delta_{31}(\mu_--\mu_+)&=\frac{33}{73}\,\kappa_{31}R_\phi(\delta_{31}).
\end{aligned}}
\]

For **nonnegative** \(\phi\), \(R_\phi\ge0\), so 31 pushes \(\mu_+\) **down** and \(\mu_-\) **up**, reinforcing the selected spectral **difference** but **not automatically the positive minimum**. For sign-changing profiles the correlation may have either sign. The identical \(B_3\) law applies to the subsequent 23 and 32 collisions with \(\delta_{23}=\log(27/23)\), \(\kappa_{23}=\log23/\sqrt{23}\), and \(\delta_{32}=\log(32/27)\), \(\kappa_{32}=\log2/\sqrt{32}\). It applies **only until a different-spacing prime event appears**, first at \(\delta_{11}=\log(11/9)\).

For profiles with \(\ell\le\delta_{31}\), \(R_\phi(\delta_{31})=0\); for \(\ell>\delta_{31}\), an explicit two-peak smooth nonnegative profile at separation \(\delta_{31}\) realizes \(R_\phi(\delta_{31})=1/2\) with no other newly offset collisions when the individual peak width is sufficiently small. This shows sharpness of the *profile supremum* of activation, **not** low energy: by the T1D-A/C argument, the normalized Gamma self-energy of these increasingly narrow two-peak profiles diverges as their individual peak width tends to zero.

## 4. A stronger invariant: nonnegative-profile phase ordering survives through **11**, not just 31

The older T1B Gamma correction to \(\mu_--\mu_+\) equals \((83A_1-64\sqrt3 A_2+33A_3)/(73J)\). For nonnegative \(\phi\), the positive Gamma kernel satisfies \(w(t+d)<w(t)/\sqrt3\) for \(t>0\); hence \(A_2<A_1/\sqrt3\) and that entire correction is **strictly positive**. All additional outermost-cell prime powers \(25,29,31,23,32\) also add \((33/73)\kappa_nR_\phi(\delta_n)\ge0\) to the splitting. The only new **spacing-two** prime power before \(\ell=\log(11/9)\) is \(8\), whose contribution to the difference is \(-64\sqrt3\kappa_8R_\phi(\delta_8)/73\), \(\delta_8=\log(9/8)\), \(\kappa_8=\log2/\sqrt8\).

Because \(\log(11/9)<2\log(9/8)\) (equivalently \(11/9<81/64\)), for every such profile of width \(\ell<\log(11/9)\), the **half-support autocorrelation lemma** gives \(0\le R_\phi(\delta_8)\le1/2\). Therefore the *complete original signed Weil* generalized splitting obeys
\[
\boxed{\mu_--\mu_+>\frac{\sqrt3}{73}\bigl(10\log3-32\kappa_8\bigr)>\frac1{16}.}
\]

The last inequality has a fully rational proof: \(\sqrt3>5/3\), \(\log3>109/100\) by a rational atanh series, \(\kappa_8=\log2/\sqrt8<1/4\) from \(\log2<7/10\) and \(\sqrt8>14/5\). The lower bound is \(29/438>1/16\). Thus, throughout this **larger entire window** (including the 31, 23 and 32 events), the normalized three-phase tests still satisfy
\[
\boxed{Q_W(u_1)<Q_W(u_i)<Q_W(u_{-1}),\qquad
\min\{Q_W(u_i)-Q_W(u_1),Q_W(u_{-1})-Q_W(u_i)\}>1/88.}
\]
Here \(u_\beta=(e_0+\beta e_1)/\|e_0+\beta e_1\|_2\). **This is a phase-splitting theorem, not a positivity theorem** for \(E_\phi\) on this extended width: either (or both) absolute Rayleigh levels could conceivably cross zero, and no profile-uniform lower bound for them is established beyond the prior T1D-C regime \(\ell<\log(31/27)\). Nor is phase ordering for *arbitrary sign-changing* profiles proved in this wider range.

## 5. The missing source-defined infinite-complement theorem: exact form Schur identity

This is the more important transfer finding. Let the fixed **two-pole-null smooth test core** be
\[
\mathcal D_a^0=\{f\in C_c^\infty(-a,a;\mathbb C):\ \int f(x)e^{\pm x/2}dx=0\},\quad a=2.
\]

For a fixed profile whose selected Rayleigh levels satisfy \(\mu_\pm>0\) (proved profile-uniformly for widths \(<\log(31/27)\) in T1D-C, **not** assumed for every wider profile), let
\[
F_\phi=\{y\in\mathcal D_a^0:\langle y,u_+\rangle_{L^2}=\langle y,u_-\rangle_{L^2}=0\}.
\]
Because the orthogonal projection of a smooth pole-null test onto the **finite** smooth space \(E_\phi\) remains in \(\mathcal D_a^0\), there is an exact algebraic and physical-orthogonal decomposition \(\mathcal D_a^0=E_\phi\oplus F_\phi\) **on the test core**, with no assumption that the full Weil operator already possesses a self-adjoint extension in this metric.

Write \(q(f,g)=W(f*\widetilde g)\), **linear in the first argument**. For \(y\in F_\phi\), set the two *actual mixed Weil cross functionals* \(b_\pm(y)=q(u_\pm,y)\). They are given by the full original signed zeta explicit formula (Gamma + von Mangoldt, with both pole terms zero because the sources lie in \(\mathcal D_a^0\)); they **are not the physical inner products**, which are zero on \(F_\phi\). For \(f=z_+u_++z_-u_-+y\), direct sesquilinear expansion and completing squares yield the **exact identity**
\[
\boxed{\begin{aligned}
q(f,f)
&=\mu_+\left|z_++\frac{\overline{b_+(y)}}{\mu_+}\right|^2
 +\mu_-\left|z_-+\frac{\overline{b_-(y)}}{\mu_-}\right|^2\\
&\qquad +\underbrace{q(y,y)-\frac{|b_+(y)|^2}{\mu_+}-\frac{|b_-(y)|^2}{\mu_-}}_{\mathscr S_\phi(y)}.
\end{aligned}}
\]

**Necessary and sufficient source-core criterion:** With \(\mu_\pm>0\),
\[
\boxed{q(f,f)\ge0\ \forall f\in\mathcal D_a^0\quad\Longleftrightarrow\quad
q(y,y)\ge\frac{|b_+(y)|^2}{\mu_+}+\frac{|b_-(y)|^2}{\mu_-}\quad\forall y\in F_\phi.}
\]

More quantitatively, for any sought **uniform margin** \(c<\min\{\mu_+,\mu_-\}\),
\[
\boxed{q(f,f)\ge c\|f\|_2^2\ \forall f\in\mathcal D_a^0\quad\Longleftrightarrow\quad
q(y,y)-c\|y\|_2^2\ge\sum_{\sigma\in\{+,-\}}\frac{|b_\sigma(y)|^2}{\mu_\sigma-c}\quad\forall y\in F_\phi.}
\]

This is **standard but exact form-level Schur algebra**, now anchored to a concrete pole-null original-Weil source and **its actual physical metric**. It is not an existence theorem for the missing complement estimate. If one verifies \(\mathscr S_\phi(y)<0\) for some authentic smooth \(y\), the completed-square formula supplies an explicit negative test direction (and would require careful independent checking of any such extraordinary assertion). Merely knowing \(q(y,y)>0\) separately does not suffice: even the finite symmetric matrix \(\begin{psmallmatrix}1&2\\2&1\end{psmallmatrix}\) has positive diagonal compressions but one negative eigenvalue because of the mixed block.

## 6. Why this is not yet Connes's Schur reduction: an exact pole-carrier mismatch

The proof-bearing finite conditioning space of Andrews–Connes–La Gatta §8.6 includes the **two pole directions**, together with near-self-dual modes, and their Lemma 8.16 assumes a **positive lower bound on the infinite-dimensional complementary block** *before* eliminating it and obtaining a **finite** Schur complement. Our chosen \(E_\phi\), by construction, has zero Mellin pole moments and so is physically **orthogonal** in \(L^2(-a,a)\) to the two-dimensional pole span \(P_a=\operatorname{span}\{e^{x/2},e^{-x/2}\}\) on the window:
\[
\boxed{E_\phi\subset P_a^\perp,\qquad E_\phi\cap P_a=\{0\}.}
\]

Hence \(E_\phi\) **cannot be identified with Connes's entire finite conditioning space**, and positivity on \(E_\phi\) is **not** the positive-complement estimate their Schur strategy requires. Also, our elimination of positive \(E_\phi\) leaves an **infinite-dimensional** Schur inequality on \(F_\phi\), whereas their paper first assumes/derives a suitable coercive \(W\) and reduces to a **finite** Schur block. This orientation mismatch is a concrete diagnostic, not an objection to the published theorem. The two constructions would require a separately established operator-domain, metric and conditioning-subspace dictionary before their inequalities could be transferred.

Finally, even positivity on the full \(\mathcal D_a^0\) would not on its own cover the **two independent pole directions** in the full \(C_c^\infty(-a,a)\) Weil criterion. In the smooth source, those extra directions can be supplied by any pair of compactly supported tests with an invertible \(2\times2\) Mellin-moment matrix, but the corresponding mixed Weil blocks would also require control. Orthogonal projection onto literal \(e^{\pm x/2}\) window functions is an \(L^2\) construction; it does **not** preserve the compactly supported smooth core, so do not silently identify the two decompositions.

## 7. Disposition and immediate next obligations

**ANALYTICALLY DERIVED; exact SymPy/rational AUDIT-PASS:** prime-31 perturbation; complete rational prime-power census through the next spacing-two event (31 → 23 → 32 → 11); persistent strictly positive **phase splitting** for all nonnegative profile shapes through \(\ell<\log(11/9)\), with \(>1/16\) generalized split and \(>1/88\) adjacent normalized phase gap; a source-typed exact, **conditional** complement Schur criterion and exact orthogonality to the canonical pole directions. None is a new global operator bound or RH theorem.

**OPEN:** profile-uniform sign of \(\mu_\pm\) after 31, arbitrary sign-changing phase ordering through 11, sign or coercivity of \(\mathscr S_\phi\) on the infinite-dimensional \(F_\phi\), addition of the two source pole directions, semilocal H2′/R1/R2 and large-\(a\) uniform estimates, and equality of completed Suzuki and Connes operator domains.

**Next CR-INT-T1D-E:** pursue either (i) the precise **joint Gamma-localization / outer-prime autocorrelation** inequality needed to protect the minimum \(\mu_+\) after 31 *without falsely decoupling profile correlation and Gamma cost*, or (ii) an explicit source-defined formula and rigorous norm bound for \(y\mapsto(b_+(y),b_-(y))\) on a chosen complement energy space. **The second is more directly relevant to the cross-canonical interface; no direct implication for RH is currently established.**

**Reproducibility:** `python3 src/audits/canonical_return_t1d_d_31_schur.py` checks the collision census, compressed prime-31 arithmetic, rational gap comparison, and finite Schur algebra. The smooth-core formula identities, Gamma monotonicity and path-correlation estimate are mathematical input arguments, **not Lean compiler outcomes**. Earlier T1B/T1C/T1D-A/B/C auditors are separate regressions.