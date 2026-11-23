# CR-INT-T1D-F — Parity-diagonal mixed Gram and asymptotic cross-block decoupling

**Canonical Return · independent analytical study · 2026-11-22 PST / 2026-11-23 UTC.**  
**Status: DERIVED / ANALYTIC on the original Weil smooth core; exact finite algebra and rational inequalities AUDIT-PASS.** **NOT** an independently evaluated sharp residual Gram interval, Lean certification, whole-complement coercivity, operator equivalence, global Weil positivity, or RH. This note advances [T1D-E](CR-INT-T1D-E.md) using the exact T1C profile and its fixed rescalings. The draft PR remains for mathematical review.

## 1. Source contract and exact objects

Primary arithmetic form: Suzuki, *Weil's quadratic form via the screw function*, [arXiv:2606.09096v4, §1.1](https://arxiv.org/html/2606.09096v4), with Suzuki's smooth-core derivative-kernel reformulation in §8.4–8.5. The separate semilocal trace identities and conditional complement hypotheses: Andrews–Connes–La Gatta, [arXiv:2611.13129v1, §§8.2, 8.5–8.6](https://arxiv.org/html/2611.13129v1). No new completed-domain identification is asserted.

Fix \(I=(-2,2)\), \(\mathcal D^0=\{f\in C_c^\infty(I;\mathbb C): \int f(x)e^{\pm x/2}dx=0\}\), and the **even** nonnegative smooth triangle–mollifier profile \(\psi=b_{49/100}*\rho_{1/100}\) from [T1C](../archive/space/10-t1c.md), supported in \([-1/2,1/2]\). For \(0<\ell\le 1/16\), put \(\phi_\ell(x)=\psi(x/\ell)\). The four translated bumps have centers \((-3d/2,-d/2,d/2,3d/2)\), \(d=\log3\). The normalized symmetric/antisymmetric pole-null generators are \(u_\pm=(e_0\pm e_1)/\|e_0\pm e_1\|_2\) as in T1B–E. Both are real and physically orthonormal; the exact full signed Weil diagonal values are \(\mu_\pm=q(u_\pm,u_\pm)>0\).

Write \(k_\pm=\mathcal K_2u_\pm\in L^2(I)\), with **the complete original signed Weil representer** \(\mathcal K_2\) of T1D-E equation (E.1), **not** Suzuki's inverse Neumann Laplacian. Let \(P=\mathrm{span}\{e^{x/2},e^{-x/2}\}\), \(\mathcal F=(E\oplus P)^\perp\subset L^2(I)\), \(r_\pm=\Pi_\mathcal F k_\pm\), \(\beta_\pm=\|r_\pm\|_2\). The finite-rank Schur penalty is \(\mathcal R=\sum_{\pm}|r_\pm\rangle\langle r_\pm|/\mu_\pm\). The question is its **actual norm**, not the *selected* Rayleigh values.

## 2. Exact parity removes the off-diagonal Gram entry

Let \(Rf(x)=f(-x)\). Because \(\psi\) is even, reversal of the four bump indices exchanges \(e_0\) and \(e_1\). Therefore \(Ru_+=u_+\), \(Ru_-=-u_-\). The full source representer \(\mathcal K_2\) commutes with \(R\) (both prime translations and the paired Gamma shifts are reflection invariant), and both \(E\) and the physical pole span \(P\) are reflection invariant. Hence \(r_+\) is even and \(r_-\) is odd. **Their physical inner product vanishes identically:**

\[
\boxed{\langle r_+,r_-\rangle=0,\quad
\mathsf B=\begin{pmatrix}\beta_+^2/\mu_+&0\\0&\beta_-^2/\mu_-\end{pmatrix},\quad
\|\mathcal R\|=\max\{\beta_+^2/\mu_+,\beta_-^2/\mu_-\}.}\tag{F.1}
\]

This is an **exact mixed-entry result**, not numerical cancellation. Set \(p_e(x)=\cosh(x/2)\), \(p_o(x)=\sinh(x/2)\) on \(I\). Their exact masses are \(\|p_e\|_2^2=\sinh2+2\), \(\|p_o\|_2^2=\sinh2-2\). From physical projection onto \(E\oplus P\),

\[
\boxed{\begin{aligned}
\beta_+^2&=\|k_+\|_2^2-\mu_+^2-\frac{|\langle k_+,p_e\rangle|^2}{\sinh2+2},\\
\beta_-^2&=\|k_-\|_2^2-\mu_-^2-\frac{|\langle k_-,p_o\rangle|^2}{\sinh2-2}.
\end{aligned}}\tag{F.2}
\]

The exact *sharp* task has thus become **two scalar norm integrals and two pole projections**, not a general 2×2 complex numerical eigensolver. Their numerically tight enclosures are still **OPEN**; (F.2) involves cancellation and must be evaluated using certified interval integration, not unverified floating-point subtraction.

## 3. A first genuine but coarse interval for the frozen T1C profile

The T1C mass and Young inequalities imply, for every rescaling of its fixed profile, \(\ell\|\phi_\ell'\|_2/\|\phi_\ell\|_2\le 175/48\). Because disjoint cells have the same profile, the corresponding normalized selected sources also satisfy \(\ell\|u_\pm'\|_2\le175/48\). For \(0<\ell<1\), split the original signed Gamma integral in (E.1) at \(\ell\) and \(1\):

- On \((0,\ell)\), the central translation difference is at most \(2t\|u'\|_2\), and the compensating term is at most \(t\|u\|_2\); with \(w(t)\le 1/t\) the integral is at most \(2\ell\|u'\|_2+\ell\|u\|_2\).
- On \((\ell,1)\), the paired translated/compensating term has norm at most \(4\|u\|_2\), giving \(4\log(1/\ell)\|u\|_2\).
- On \((1,\infty)\), the T1D-E estimate gives \(12\|u\|_2\). The original Gamma constant has magnitude <4 and the 24 genuine prime-power paired translations contribute <48.

The elementary weight inequality \(w(t)\le 1/t\) for \(0<t\le1\) follows from \(1-e^{-2t}\ge2t/(1+2t)\) and \((1+2t)e^{-t/2}<2\) on this interval. It follows that

\[
\boxed{\|\mathcal K_2u\|_2<2\ell\|u'\|_2+
(64+\ell+4\log(1/\ell))\|u\|_2.}\tag{F.3}
\]

At the **actual T1C width** \(\ell=10^{-4}\), the fully rational exponential-series proof \(\log10000<10\) and \(2\ell\|u'\|_2\le175/24\) give \(\boxed{\|k_\pm\|_2<112}\). With T1C's strictly certified \(\mu_\pm>7\), equation (F.2) yields

\[
\boxed{0\le\mathsf B_{++}<1785,\quad
0\le\mathsf B_{--}<1785,\quad
\mathsf B_{+-}=\mathsf B_{-+}=0,\quad
0\le\|\mathcal R\|<1785.}\tag{F.4}
\]

This is an **actual analytic interval enclosure** for the T1C source, but deliberately very loose. It is not an independently computed residual norm or a useful numerical certificate for Schur closure. In particular, its upper endpoint is orders of magnitude above the previously certified selected \(\mu_+>7\).

## 4. Stronger structural result: shrinking the fixed shape forces the Schur penalty to zero

The crude first/second derivative bounds obscure a favorable cancellation. Let \(\mathcal G\) denote the **whole-line translation-invariant Gamma-only part** of the original Weil representer (with no primes). Its Fourier multiplier, with unitary Fourier transform, is

\[
\boxed{m_\Gamma(\xi)=-C_\Gamma+2\int_0^\infty
\big(e^{-t/2}-\cos(t\xi)\big)w(t)\,dt.}\tag{F.5}
\]

This integral converges for each real \(\xi\). Write \(m_\Gamma=c_0+I+E_{<1}+E_{\ge1}\), where \(I(\xi)=\int_0^1(1-\cos(t\xi))dt/t\). The elementary bounds \(|c_0|<14\), \(|I(\xi)-\log(1+|\xi|)|<5\), \(|E_{<1}|\le4\), and \(|E_{\ge1}|<16\) give the **uniform scalar multiplier estimate**

\[
\boxed{|m_\Gamma(\xi)-\log(1+|\xi|)|<39<40\qquad(\xi\in\mathbb R).}\tag{F.6}
\]

For completeness, the small-\(t\) remainder bound follows from \(|w(t)-1/(2t)|\le1\) for \(0<t\le1\): write \(1/(1-e^{-2t})=1/(2t)+d(t)\) and use \(0\le d(t)\le1\) from \(1-e^{-x}\in[x/(1+x),x]\). For \(I\), split at frequency 1 and integrate \(\cos s/s\) by parts to bound its tail by 3. The constants 14 and 16 use \(|C_\Gamma|<4\), \(w(t)\le2/t\) below 1 and \(w(t)\le2e^{-t/2}\) above 1.

Writing \(\hat\psi\) for the unitary Fourier transform of the **fixed** triangle–mollifier shape, let

\[
L_\psi^2:=\frac{1}{\|\psi\|_2^2}\int_\mathbb R(1+|\log|\eta||)^2|\hat\psi(\eta)|^2\,d\eta.
\]

This is finite: \(\psi\) is smooth and compactly supported. More precisely, on \(|\eta|\le1\) the Cauchy–Schwarz/support bound gives \(|\hat\psi|\le\|\psi\|_2/\sqrt{2\pi}\), yielding an integral at most \(5/\pi\); on \(|\eta|\ge1\), \((1+\log|\eta|)^2\le\eta^2\), yielding at most \(\|\psi'\|_2^2/\|\psi\|_2^2\le30625/2304\). Since \(\pi>3\),

\[
\boxed{L_\psi^2\le5/\pi+30625/2304<16,\quad L_\psi<4.}\tag{F.7}
\]

For \(0<\ell\le1\), Fourier rescaling and the pointwise inequality \(|\log(\ell+|\eta|)|\le 1+|\log|\eta||\) give

\[
\boxed{\|\mathcal G\phi_\ell-\log(1/\ell)\phi_\ell\|_2
\le(40+L_\psi)\|\phi_\ell\|_2<44\|\phi_\ell\|_2.}\tag{F.8}
\]

The normalized four-cell combinations have \(\sum|c_j|/\sqrt{\sum c_j^2}\le2\). Apply (F.8) to each translated bump, then use the original fixed-window prime-translation norm bound <48 and restrict to \(I\):

\[
\boxed{\|\mathcal K_2u_\sigma-\log(1/\ell)u_\sigma\|_2<2(44)+48=136.}\tag{F.9}
\]

The physical complement projector annihilates \(u_\sigma\), so **uniformly** \(\beta_\sigma\le136\) as the fixed bump *shape* is scaled towards zero. T1D-A's exact width-dependent lower bounds give \(\mu_\pm>\log(1/\ell)-3\) for \(\ell\le1/16\). Therefore, for \(0<\ell\le1/32\),

\[
\boxed{0\le\|\mathcal R_\ell\|=
\max_\sigma\frac{\beta_\sigma^2}{\mu_\sigma}
\le\frac{18496}{\log(1/\ell)-3}\ \xrightarrow[\ell\downarrow0]{}\ 0.}\tag{F.10}
\]

**This is the main T1D-F consequence:** with a fixed smooth rescaled shape, the *entire selected-to-complement Schur penalty*, measured as a physical \(L^2\) operator, tends to zero as localization narrows, even while the selected Rayleigh values diverge logarithmically and their arithmetic separation stays nonzero. The infinite-dimensional complement varies with \(\ell\), so the statement is about the norms of a **family of rank-two source-defined operators**, not norm convergence of a single fixed-space operator or an equivalence of Suzuki/Connes completions.

## 5. Why this still does not close the complement

A direct, source-valid lower estimate exists but is too weak: for every pole-null smooth \(y\in\mathcal D^0\), \(\Re\langle y(\cdot+t),y\rangle\le\|y\|_2^2\). With the complete Gamma scalar formula, \(|C_\Gamma|<4\), \(2(1-e^{-t/2})w(t)\le2\) for \(t<1\), and \(w(t)\le2e^{-t/2}\) for \(t\ge1\), the Gamma contribution is \(>-14\|y\|_2^2\). The 24 source-active prime powers contribute at worst \(>-48\|y\|_2^2\). Consequently

\[
\boxed{q(y,y)>-62\|y\|_2^2\qquad(0\ne y\in\mathcal D^0).}\tag{F.11}
\]

This **semiboundedness** is not the required **positive** complementary coercivity. Even a Schur penalty that tends to zero cannot dominate a complementary form for which one has not proved nonnegativity. The exact T1D-D criterion remains \(q(y,y)\ge\langle\mathcal R_\ell y,y\rangle\) for all smooth \(y\in F_\ell\); neither (F.10) nor (F.11) implies it. There is **no verified negative test direction** and no new information about RH.

There is an additional parity simplification: \(\mathcal F=\mathcal F_{\rm even}\oplus\mathcal F_{\rm odd}\), \(q\) has no even–odd mixed terms, and the Schur criterion decomposes into the separate conditions \(q(y_e)\ge|\langle r_+,y_e\rangle|^2/\mu_+\) and \(q(y_o)\ge|\langle r_-,y_o\rangle|^2/\mu_-\). This matches the *parity labels* seen in Suzuki and Connes, but does not identify their spectral subspaces or domains.

## 6. Custody / falsifiers / next

**ANALYTICALLY DERIVED:** exact zero off-diagonal projected mixed Gram, exact scalar projection formulas (F.2), fixed-T1C certified loose norm interval (F.4), explicit Fourier/Gamma multiplier estimate and **asymptotic suppression** (F.10), and only a negative complementary lower margin (F.11). **AUDIT-PASS:** reproducible SymPy and rational arithmetic for the finite identities, source mass and bounds. **OPEN:** tight numerical enclosure for the two positive residual norms in (F.2); an effective *positive* bound for \(q\) on the varying smooth complement; incorporation of the full two-pole directions and semilocal H2′/R1/R2; completed-domain transfer; RH. No numerical eigensolver result is passed off as certified integral evaluation.

**Falsifiers:** a source sign/normalization mismatch in (E.1); failure of reflection invariance; failure of T1C's derivative-energy bound; a counterexample to the uniform multiplier comparison (F.6); invalid Plancherel transform conventions; or a missing source support/prime-power cutoff. Correct all these before attempting to use (F.10) in an operator transfer theorem.

**Next CR-INT-T1D-G:** compute tighter *independently certified* intervals for the two scalar residual norms using (F.2) (not subtractive double precision), and separately identify whether Connes's actual conditioned Sonin complement provides a **uniform positive form lower bound on the same native physical space**. If not, record an explicit domain/metric mismatch rather than transferring a fictitious coercivity estimate.