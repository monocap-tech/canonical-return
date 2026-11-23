# CR-SZ-003 — Prime-power threshold, full-form phase splitting, and Suzuki’s variational bound

**Canonical Return / Suzuki / independent analytic note / 2026-11-22 PST (2026-11-23 UTC).** Status: **PROVED ON THE SMOOTH CORE BY DIRECT ESTIMATES; NOT LEAN FORMALIZED; NOT A CLAIM OF NOVELTY OVER PRIOR PI027; NOT RH.** Predecessors: [CR-SZ-001](https://chatgpt.com/space/page_d744fb3bc0dc8191b12101c36536a925), [CR-SZ-002](https://chatgpt.com/space/page_8c63ddb9ebb481919002432bc72dab3b). All source repositories and proof workflow jurisdictions remain unchanged.

## 1. Sources and scope

- **Suzuki** arXiv:2606.09096**v4**, especially §1.1 signed Weil functional, §2.3 prime-translation identity (2.5), Corollary 1.2 smooth-core Rayleigh infimum, Theorem 1.3 continuity of \(\lambda_a\), §8.4 (8.7), §8.5 (8.8). [Primary HTML](https://arxiv.org/html/2606.09096v4).

- **Weil-Lab IV030**: [immutable phase construction](https://github.com/monocap-tech/weil-lab/blob/1061a331c195e30e85dd8eb35a6eff1c794cf054/docs/IV030_PHASE_IV_ORIGINAL_PRIME_PHASE_FIDELITY_NO_GO.md), with genuine \(\kappa=\Lambda(9)/\sqrt9=\log3/3\) and valid physical two-collar Schwartz tests.

- **Existing PI027** [PI registry #792](https://github.com/monocap-tech/weil-lab/issues/792): already establishes original full-form prime/Gamma phase competition and n=9 dominance for a different pole-null family, analytically only. Our result below transfers and quantifies a simpler **nonnegative bump** case with **both poles retained**. Do not present it as the first observation of full-form phase sensitivity.

## 2. The exact n=9 admissibility threshold

Write \(\delta=\log9=2\log3\), \(\kappa=\log3/3\), and fix an arbitrary Suzuki support radius **\(a>\log3\)**. Choose a positive width
\[0<\ell<\min\{1/5000,\;2(a-\log3),\;\log(9/8)\}.\]
Choose any **nonzero real nonnegative** \(\phi\in C_c^\infty(-\delta/2-\ell/2,-\delta/2+\ell/2)\); define \(\psi(x)=\phi(x-\delta)\). Both test functions lie in \(C_c^\infty(-a,a)\). If desired, the width may be reduced further to lie inside the frozen original NF543 two collars. Let \(J=\int\phi^2=\int\psi^2>0\).

Every same-bump prime-power correlation vanishes because the bump width \(\ell<\log2\). Every cross-bump prime-power correlation vanishes except \(n=9\) because the difference set lies in \((\delta-\ell,\delta+\ell)\), with \(\ell<\log(9/8)<\log(11/9)\), and 10 is not a prime power. At \(n=9\) the two orientations contribute \(\kappa J\) to the off-diagonal prime pairing.

Conversely, if \(a\le\log3\), the difference set of supports within \((-a,a)\) is contained in \((-2a,2a)\) and cannot contain \(\pm\log9\); so the \(n=9\) term vanishes for every smooth supported test. **The strict radius threshold for possible n=9 interactions is \(a=\log3\).** This does not assert that n=9 is the first possible prime interaction (n=2 occurs at smaller radii).

## 3. The FULL Weil mixed coefficient retains a definite negative sign

Use Suzuki's **linear-first** convention \(Q(v,w)=W(v*\widetilde w)\). Let \(M_+(\phi)=A=\int\phi(x)e^{x/2}dx>0\), \(M_-(\phi)=B=\int\phi(x)e^{-x/2}dx>0\). Translation by \(\delta=\log9\) gives \(M_+(\psi)=3A\), \(M_-(\psi)=B/3\). The **two poles** therefore contribute \((3+1/3)AB=(10/3)AB\) to \(c:=Q(\psi,\phi)\).

The prime contribution is \(-\kappa J\). For the Gamma term, the cross-correlation \(k=\psi*\widetilde\phi\) is real nonnegative, supported in \((\delta-\ell,\delta+\ell)\subset(0,\infty)\), and satisfies \(k(0)=0\). Suzuki's direct archimedean formula yields a **strictly negative** mixed contribution \(-I_\Gamma\), where
\[I_\Gamma=\int\!\!\int \phi(x)\phi(y)w(\delta+x-y)\,dx\,dy>0,\quad w(t)=\frac{e^{-t/2}}{1-e^{-2t}}>0.\]
Therefore the exact **full signed Weil** cross-pairing is
\[\boxed{c=\frac{10}{3}AB-\kappa J-I_\Gamma.}\]
All original contributions are included; this is *not* the isolated prime block from CR-SZ-001.

## 4. Fully profile-independent quantitative bounds

Let \(\ell\) bound the support width. Cauchy–Schwarz and \(e^{x/2}e^{-y/2}\le e^{\ell/2}\) on the support give
\[0<AB\le e^{\ell/2}\Big(\int\phi\Big)^2\le e^{\ell/2}\ell J.\]
As \(\delta-\ell>\log8\), we have \(0<w(t)<1\) on the cross-difference set, hence \(0<I_\Gamma<\ell J\). Consequently
\[\boxed{-\kappa-\ell<\frac{c}{J}< -\kappa+\frac{10}{3}e^{\ell/2}\ell.}\]
Since \(\ell<1/5000\), \(e^{\ell/2}<2\), \(\kappa=\log3/3>1/3\), this yields the **rational lower margin**
\[\boxed{-\frac{c}{J}>\frac13-\frac1{750}=\frac{83}{250}>0.}\]
The bound is uniform over **all** eligible nonnegative smooth bump profiles, however narrow; it does not rely on computing their exact integrals. In particular, \(c/J\to-\kappa\) uniformly as \(\ell\downarrow0\). The number \(83/250\) is deliberately coarse; no claim of optimality.

## 5. Consequences for exact phase tests and the 2D generalized pencil

Translation invariance of autocorrelation gives identical diagonal entries \(q:=Q(\phi,\phi)=Q(\psi,\psi)\), real and finite. The full Gram matrix on the ordered physical basis \((\phi,\psi)\) is
\[H=\begin{pmatrix}q&c\\c&q\end{pmatrix},\qquad M=J I_2.\]
For \(f_\beta=\phi+\beta\psi\), \(|\beta|=1\), the support set and physical mass \(\|f_\beta\|_2^2=2J\) are fixed and the genuine prime-interaction support graph is fixed, but
\[Q(f_\beta)=2q+2c\operatorname{Re}\beta.\]
Hence the full original signed Weil quadratic values have the **strict ordering**
\[\boxed{Q(f_1)<Q(f_i)<Q(f_{-1}),}\]
with each successive gap \(-2c>\tfrac{83}{125}J\). This strengthens the IV030 *prime-only* phase obstruction for this positive-bump family to the **complete original signed form**, using both poles and Gamma rather than deleting them. It is consistent with the distinct already retained PI027 analytic full-form phase results.

In Suzuki's §8.5 generalized quotient, the denominator matrix is \(M=JI_2\), because \(\langle K_a Df,Dg\rangle=\langle f,g\rangle\) on the smooth core. The **restricted two-dimensional generalized eigenvalues** are
\[\mu_\pm=\frac{q\pm c}{J}\quad(c<0),\qquad \mu_{\max}-\mu_{\min}=\frac{-2c}{J}>\frac{83}{125}.\]
**This is a spectral gap for the selected 2×2 compression, NOT a gap between genuine eigenvalues of Suzuki's full operator.** Its lower variational value is \((q+c)/J\); by Corollary 1.2,
\[\boxed{\lambda_a\le \frac{q+c}{J}<\frac qJ-\kappa+\frac1{750}.}\]
This is a legitimate **upper bound depending on the uncomputed full diagonal \(q\)**. It is not a sign result and not a lower bound on \(\lambda_a\).

## 6. Why this does NOT create a discontinuity at the prime threshold

Suzuki already proved continuity of \(\lambda_a\) in his Theorem 1.3. Our strict off-diagonal phase splitting appears whenever \(a>\log3\), but for \(a\downarrow\log3\), the admissible bumps must narrow. The diagonal cost grows with localization; for a fixed nonzero normalized profile \(\eta\in C_c^\infty(-1/2,1/2)\) and scaled bump \(\phi_\ell(x)=\ell^{-1/2}\eta((x+\log3)/\ell)\), the archimedean Fourier multiplier \(B(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi=\log|r|+O(1)\) at high frequency gives \(q_\ell/J_\ell=\log(1/\ell)+O_\eta(1)\). The prime diagonal vanishes and the pole diagonal is \(O(\ell J_\ell)\). The two-dimensional minimal quotient therefore grows rather than driving a negative eigenvalue when \(\ell\downarrow0\).

This is a standard high-frequency scaling consequence, recorded as an **analytic asymptotic sketch** rather than a newly formally compiled theorem. The strict n=9 interaction threshold therefore does not contradict continuous \(\lambda_a\), and it does not force a new eigenvalue crossing.

## 7. Classification and explicit next gate

**Outcome:** (1) **exact support threshold** \(a=\log3\) for the n=9 term (elementary recovery); (2) **uniform complete-Weil phase-splitting estimate** for positive smooth bump pairs, retaining both poles and Gamma (analytic consequence; related to earlier PI027), (3) **exact restricted generalized-pair spectral gap** and **full \(\lambda_a\) upper estimate**, no numerical sign; (4) **no full spectral lower bound**, since the complement of the 2D span is uncontrolled, (5) **no RH proof / no Lean acceptance**.

**First missing hypothesis for any full lower bound:** a complete form-compatible decomposition of the entire localized domain into the selected span plus its orthogonal complement, together with rigorous control of the complement block \(D\) and mixed couplings \(C\). If \(D\ge dI\), control of \(CD^{-1}C^*\) and a suitable Schur bound could yield a lower estimate; no such bounds follow from the 2D compression. The old PI027 cancellation examples also warn that full phase behavior depends on the test family.

**Next separate author:** [Morishita, Theorem 3.6](https://arxiv.org/abs/2508.15971): inspect whether phase-space/orbit transport preserves any measure, pairing, or phase information. Do not infer a Weil/Hodge pairing from orbit equivariance.
