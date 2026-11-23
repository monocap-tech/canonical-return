# CR-INT-T1D-R — The sharp archimedean floor and an exact compact spectral gate for the original low-frequency Schur complement

**Canonical Return · 2026-11-23 · independent source-level analytical investigation.**  
**Disposition:** **DERIVED / ANALYTIC** (original complete signed Weil smooth source) with an executable exact-rational audit for the sharp Gamma floor, the finite-window prime census, the negative zero-frequency *formal* symbol and the conditional Galerkin constant. **NOT** a computation of the low-mode top eigenvalue, a positive theorem on the entire pole-null complement, a completed Suzuki–Connes correspondence, Lean certification, or RH. This is the next entry after [T1D-Q](CR-INT-T1D-Q.md).

## 1. Preserve the source and distinguish the canonical compactness result

Fix \(I=(-2,2)\), \(H=L^2(I;\mathbb C)\), physical inner product linear in its **first** argument, and
\[
\mathcal D^0=\{f\in C_c^\infty(I):\int_I f(x)e^{x/2}dx=\int_I f(x)e^{-x/2}dx=0\}.
\]
Let \(E=\operatorname{span}\{u_+,u_-\}\subset\mathcal D^0\) be the previously frozen T1C two-dimensional selected space with \(\mu_\pm=q(u_\pm,u_\pm)>0\). Let \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\subset H\), \(\mathscr H_F=(E\oplus P)^\perp\), and \(F_s=C_c^\infty(I)\cap\mathscr H_F\). The previously derived smooth-source Schur condition is
\[
\boxed{q(f,f)\ge0\ (f\in\mathcal D^0)\iff q(y,y)-\langle\mathcal R y,y\rangle\ge0\ (y\in F_s),}\tag{R.1}
\]
where \(\mathcal R=\sum_{\sigma=\pm}|r_\sigma\rangle\langle r_\sigma|/\mu_\sigma\) is the positive, rank-at-most-two **physical** cross-block penalty of T1D-E–I. Here \(r_\sigma=\Pi_{\mathscr H_F}\mathcal K_2u_\sigma\in\mathscr H_F\), and all genuine Weil cross-pairings, not physical inner products, are retained.

**Original mathematical provenance:** Suzuki, *Weil's quadratic form via the screw function*, [arXiv:2606.09096v4](https://arxiv.org/html/2606.09096v4), §1.1, §2.4 eqs. (2.7)–(2.8), Theorem 1.1, and Proposition 4.1. Suzuki already proves a finite-window Friedrichs extension and compact embedding/discrete spectral properties. **Do not claim these pre-existing operator results as new.** This note makes a different, explicit **Schur-corrected, pole-null, parity-separated** spectral reduction preserving the selected physical subspace and source constraints. The semilocal conditioning work of Andrews–Connes–La Gatta, [arXiv:2611.13129v1](https://arxiv.org/html/2611.13129v1), has a distinct arithmetic transport and completed domain; no identification is made here.

## 2. Sharp unconditional lower bound for the original Gamma part

By Suzuki §2.4, on smooth zero-extended tests the **complete original Gamma** form is the Fourier multiplier
\[
\gamma[f]=\int_{\mathbb R}m_\Gamma(\xi)|\widehat f(\xi)|^2d\xi,\qquad m_\Gamma(\xi)=\operatorname{Re}\psi_{\rm dg}(\tfrac14+\tfrac{i\xi}{2})-\log\pi,
\]
using the **unitary** Fourier transform. The digamma series implies the exact pointwise monotonicity
\[
\operatorname{Re}\psi_{\rm dg}(a+it)-\psi_{\rm dg}(a)
=\sum_{k=0}^\infty\frac{t^2}{(k+a)((k+a)^2+t^2)}\ge0,\qquad a>0.
\]
Therefore its global multiplier infimum occurs at \(\xi=0\) **exactly**, with
\[
\boxed{\inf_\xi m_\Gamma(\xi)=m_\Gamma(0)=-\gamma_{\rm EM}-\frac\pi2-3\log2-\log\pi>-\frac{43}{8}.}\tag{R.2}
\]
The strict rational bound is audited from an alternating Machin-series enclosure for \(\pi\), an atanh-series enclosure for logarithms, and harmonic-number bounds for Euler's constant. In particular \(\gamma[f]>-43\|f\|_2^2/8\) for every nonzero smooth test; the sharper exact constant is approximately \(-5.3721834192\), and the script's certified outward interval encloses it.

Write the **true compressed arithmetic operator**
\[
A=\sum_{n=p^k\le54}\frac{\Lambda(n)}{\sqrt n}\,P_I(T_{\log n}+T_{-\log n})P_I.
\]
This is the same 24-prime-power signed source as T1D-Q: **the Weil form contains \(-A\)**. Importing **T1D-Q's independently audited** joint positive-supersolution bound \(\|A\|<203/25\), we get the improved, unconditional-to-R (subject to Q's source proof) full pole-null source semibound
\[
\boxed{q(y,y)>-\left(\frac{43}{8}+\frac{203}{25}\right)\|y\|_2^2=-\frac{2699}{200}\|y\|_2^2\qquad (0\ne y\in\mathcal D^0).}\tag{R.3}
\]
This improves the earlier coarse \(-62\) type semibound on **all** smooth pole-null directions, not only the frequency-conditioned portion. **It is negative: no complement positivity follows.**

A useful falsifier for simplistic Fourier-symbol proofs: if we artificially extend the *finite-window truncated* arithmetic sum as a whole-line translation multiplier, its value at \(\xi=0\) is
\[
 m_\Gamma(0)-2\sum_{p^k\le54}\frac{\log p}{\sqrt{p^k}}\in(-30,-29).
\tag{R.4}
\]
The audit encloses the exact sum; (R.4) establishes **failure of a naive nonnegative pointwise Fourier symbol**. It is **not** a negative source test: functions in \(C_c^\infty(-2,2)\) cannot be bandlimited to \(\xi=0\), and the two Mellin moments additionally constrain the admissible class. No negative original Weil direction is claimed.

## 3. A compact Birman–Schwinger gate for the *actual* Schur-corrected complement

Define the positive Gamma-only form on \(F_s\)
\[
\boxed{d[y]=\gamma[y]+6\|y\|_2^2>\frac58\|y\|_2^2.}\tag{R.5}
\]
The positivity comes from (R.2), not RH. Let \(\mathscr V_F\) be the completion of \(F_s\) in the norm \(d[y]^{1/2}\), with its continuous embedding into \(\mathscr H_F\). It is densely embedded: the four physical moment functionals (the two pole carriers and two selected modes) can be corrected using fixed compact smooth tests with an invertible moment matrix. This correction is bounded in the logarithmic form norm, so smooth constrained tests are a form core.

Because \(m_\Gamma(\xi)+6\) grows like \(\log(2+|\xi|)\) and is bounded below by \(5/8\), the shifted Gamma form is closed on its logarithmic Sobolev completion. The embedding \(\mathscr V_F\hookrightarrow\mathscr H_F\) is compact: low-Fourier-band truncations on a fixed finite interval are compact integral operators, and bounded form energy forces their high-frequency \(L^2\) tail uniformly to zero as the cutoff rises. This is an application of Suzuki's **canonical** compact-embedding mechanism (Proposition 4.1), specialized to the four physical constraints here.

Let \(D_F\ge(5/8)I\) be the positive self-adjoint operator associated with the closed form \(d\) on \(\mathscr H_F\). It has **compact resolvent**. Define the bounded self-adjoint operator
\[
\boxed{C_F=\Pi_{\mathscr H_F}A\Pi_{\mathscr H_F}\big|_{\mathscr H_F}+6I_{\mathscr H_F}+\mathcal R.}\tag{R.6}
\]
For every \(y\in F_s\), the original complete signed Weil Schur remainder is **exactly**
\[
q(y,y)-\langle\mathcal R y,y\rangle
=d[y]-\langle C_Fy,y\rangle.
\]
Consequently the compact self-adjoint **joint Gamma–prime–Schur operator**
\[
\boxed{\mathscr B_F=D_F^{-1/2}C_FD_F^{-1/2}}\tag{R.7}
\]
satisfies the equivalent criterion
\[
\boxed{q(y,y)\ge\langle\mathcal R y,y\rangle\ \forall y\in F_s
\iff\ \sup\sigma(\mathscr B_F)\le1.}\tag{R.8}
\]
**Proof:** for every form-domain vector \(y\), put \(z=D_F^{1/2}y\); the map is onto \(\mathscr H_F\) because \(D_F\ge5/8\). Then
\(d[y]-\langle C_Fy,y\rangle=\langle(I-\mathscr B_F)z,z\rangle\). Density of \(F_s\) in the form domain proves the equivalence; \(D_F^{-1/2}\) is compact by compact resolvent, hence so is \(\mathscr B_F\). **This is an exact criterion, not its verification.** The operator \(\mathscr B_F\) need not be positive semidefinite: \(A\) is positivity-preserving as a kernel, not necessarily a positive quadratic operator.

The parity involution \(Jf(x)=f(-x)\) commutes with the Gamma form, the true prime operator, \(\Pi_{\mathscr H_F}\), and the source-defined \(\mathcal R\). Thus
\[
\boxed{\mathscr B_F=\mathscr B_{F,+}\oplus\mathscr B_{F,-},\qquad
\sup\sigma(\mathscr B_F)=\max\{\sup\sigma(\mathscr B_{F,+}),\sup\sigma(\mathscr B_{F,-})\}.}\tag{R.9}
\]
This **does not** assert equality with Connes's semilocal conditioned Sonin operators. Its purpose is to replace a gigantic *prolate rank count* by two precise compact spectral questions about the actual signed Weil complement.

## 4. An explicit, finite-step certification/falsification rule

Let \(0<d_1\le d_2\le\cdots\uparrow\infty\) be the eigenvalues of \(D_F\), counting multiplicity, and let \(\Pi_N\) be its spectral projection onto the first \(N\) eigenvectors (similarly within each parity sector). The matrix
\[
\mathscr B_N=\Pi_N D_F^{-1/2}C_FD_F^{-1/2}\Pi_N
\]
retains **all original Gamma, von Mangoldt and Schur cross-pairings**. By splitting the cross and tail blocks, and using \(D_F\ge5/8\),
\[
\boxed{\begin{aligned}
\sup\sigma(\mathscr B_N)&\le\sup\sigma(\mathscr B_F)\\
&\le\sup\sigma(\mathscr B_N)+\varepsilon_N,\\
\varepsilon_N&\le\frac{2\|C_F\|}{\sqrt{d_1d_{N+1}}}+\frac{\|C_F\|}{d_{N+1}}.
\end{aligned}}\tag{R.10}
\]
The first inequality is the Rayleigh–Ritz principle; the second follows from \(\|\mathscr B_F-\Pi_N\mathscr B_F\Pi_N\|\le\varepsilon_N\). The bound tends to zero because \(d_n\to\infty\). **A rigorous upper enclosure** for \(\sup\sigma(\mathscr B_N)\) and a certified lower enclosure for \(d_{N+1}\) with
\(\sup\sigma(\mathscr B_N)+\varepsilon_N\le1\) would establish the entire smooth pole-null Schur inequality, not merely a prolate-conditioned portion. **Conversely**, a verified Ritz eigenvalue \(>1\) gives a negative Schur remainder. Since \(F_s\) is a form core, it would yield smooth source counterexamples after approximation and demand extraordinary independent review.

For a concrete numerical stopping rule, **conditional on T1D-I's separately derived** \(\|\mathcal R\|<155/182\), the new Q bound implies
\[
\|C_F\|<6+\frac{203}{25}+\frac{155}{182}<15,
\]
so
\[
\boxed{\varepsilon_N<30\sqrt{\frac{8}{5d_{N+1}}}+\frac{15}{d_{N+1}}.}\tag{R.11}
\]
Without this optional imported residual enclosure, \(C_F\) is still bounded by T1D-E's explicit one-sided mixed norm, but this \(<15\) stopping constant must **not** be used as unconditional evidence. Neither parity block's top eigenvalue has been enclosed; the new stopping rule is mathematical infrastructure, not a completed certification.

## 5. Source-to-source preservation and hard stopping point

| Obligation | Status |
| --- | --- |
| Sharp source Gamma multiplier infimum and \(>-43/8\) rational floor | **DERIVED / exact-rational AUDIT-PASS** |
| Combined genuine prime operator \(\|A\|<203/25\) | **Imported T1D-Q** |
| Full pole-null Weil semibound \(>-2699/200\) | **DERIVED using Q, still negative** |
| Compact-resolvent property of original finite-window Weil operator | **CANONICAL: Suzuki**, not new |
| Exact compact joint Gamma–prime–Schur operator on the physical complement | **DERIVED / ANALYTIC**, based on canonical compact embedding |
| Parity splitting and finite Ritz tail error bound | **DERIVED / ANALYTIC**; finite algebra audit only |
| Numerical top Ritz eigenvalues \(\sup\sigma(\mathscr B_{F,\pm})\) | **OPEN** |
| Entire source Schur positivity, two pole carriers, support-independent transfer, RH | **OPEN / NO TRANSFER** |

**Reproducibility:** `python3 src/audits/canonical_return_t1d_r_compact_schur.py`. The script performs exact rational intervals for \(\pi\), \(\gamma_{\rm EM}\), logarithms, square roots and all 24 prime powers; verifies the new Gamma floor, the artificial symbol diagnostic (R.4), the imported conditional \(<15\) bound and a finite symbolic Galerkin example. It does **not** certify the infinite-dimensional compactness proof or compute any low-mode spectral eigenvalue. No completed Suzuki/Connes operator equivalence or RH follows.

**Next — CR-INT-T1D-S:** construct actual smooth even/odd low-frequency Ritz test families for \(\mathscr B_{F,\pm}\), rigorously compute their full signed Weil Gram matrices and Dirichlet/Gamma masses, and seek **one-sided** certified top-eigenvalue information. A positive finite Ritz matrix alone cannot prove (R.8); the missing upper tail enclosure must be included in every proposed positive certificate. Do not restart prolate rank enumeration without evidence it improves this direct joint spectral gate.
