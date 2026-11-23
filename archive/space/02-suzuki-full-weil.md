# CR-SZ-002 — Full scalar explicit-formula dictionary, Gamma normalization and the source-carrier boundary

**Independent Canonical Return | 2026-11-22 PST / 2026-11-23 UTC. Outcome: FULL SCALAR FORMULA RECOVERY ON SMOOTH CORE; Zeta23/Weil-Lab CARRIER EXTENSION UNPROVED; NOT A LEAN THEOREM; NOT A NEW RH RESULT.** Continues [CR-SZ-001](https://chatgpt.com/space/page_d744fb3bc0dc8191b12101c36536a925). This page has no effect on Weil-Lab's proofs, workflows, or canonical acceptance registry.

## Exact sources and version pins

- Suzuki, [*Weil’s quadratic form via the screw function*, arXiv:2606.09096v4 (2026-11-15)](https://arxiv.org/html/2606.09096v4), §1.1 displayed definition of \(W\), §1.2 (1.3), §2.3 (2.5), §2.4, §2.5 (2.9), §8.4 (8.7), §8.5 (8.8).

- Frozen original Weil-Lab mathematical research commit **`a431f903254000bb5ad8edbc276ccc8599da98f5`**: [`WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`](https://github.com/monocap-tech/weil-lab/blob/a431f903254000bb5ad8edbc276ccc8599da98f5/WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean), Lean definition `Zeta23.EF.literatureRHS` and `gammaBracket`; [`WeilDefect/Arithmetic/ActualZetaCorrelationTransform.lean`](https://github.com/monocap-tech/weil-lab/blob/a431f903254000bb5ad8edbc276ccc8599da98f5/WeilDefect/Arithmetic/ActualZetaCorrelationTransform.lean), `neutralWindowCorrelation`; [`WeilDefect/Arithmetic/ActualZetaLiteralFormula.lean`](https://github.com/monocap-tech/weil-lab/blob/a431f903254000bb5ad8edbc276ccc8599da98f5/WeilDefect/Arithmetic/ActualZetaLiteralFormula.lean), `neutralActualZetaGreenZeroForm_arithmetic`; [`WeilDefect/Arithmetic/ActualZetaOriginalCollarLiteralArch.lean`](https://github.com/monocap-tech/weil-lab/blob/a431f903254000bb5ad8edbc276ccc8599da98f5/WeilDefect/Arithmetic/ActualZetaOriginalCollarLiteralArch.lean), original Fourier/Gamma source identity.

- The correspondence below is a **new analytic comparison** of source definitions, not a newly elaborated Lean equality. Suzuki's operator identity (8.7) is already canonical.

## 1. Scalar signed Weil functional: exact three-part comparison

Let \(k\in C_c^\infty(\mathbb R)\) and adopt the positive-exponent transform \(\widehat k(r)=\int_\mathbb R k(x)e^{irx}\,dx\). Suzuki §1.1 defines the real-line distributional functional
\[\begin{aligned}
W_{\rm Su}(k)=&\underbrace{\int k(x)(e^{x/2}+e^{-x/2})dx}_{\text{two poles}}\\
&-\underbrace{\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}(k(\log n)+k(-\log n))}_{\text{signed prime}}\\
&-\underbrace{(\log(4\pi)+\gamma)k(0)-\int_0^\infty\bigl(k(x)+k(-x)-2e^{-x/2}k(0)\bigr)\frac{e^{x/2}}{e^x-e^{-x}}\,dx}_{\text{archimedean}}.
\end{aligned}\]
The last two lines are to be read with the displayed minus signs: the archimedean piece is \(- (\log(4\pi)+\gamma) k(0) - \int(\cdots)w(x)dx\).

The frozen Zeta23 `literatureRHS k` is
\[R_{\rm Z23}(k)=\widehat k(i/2)+\widehat k(-i/2)-\sum_{n\ge0}\kappa_n\left(k(L_n)+k(-L_n)\right)+\frac1{2\pi}\int_\mathbb R\widehat k(r)B(r)dr,\]
where \(\kappa_n=\Lambda(n)/\sqrt n\), \(L_n=\log n\), and \(B(r)=\operatorname{Re}\psi(1/4+ir/2)-\log\pi\). The \(n=0,1\) terms vanish by the arithmetic definition, so the signed prime terms agree.

## 2. Two original zeta poles are the same cross-moments

\[\widehat k(i/2)+\widehat k(-i/2)=\int_\mathbb R k(x)(e^{-x/2}+e^{x/2})dx.\]
Thus Suzuki's first term equals the literal Zeta23 pole evaluations, with no missing factor of 2 or \(\pi\). For mixed \(k=t*\widetilde s\), \(\widetilde s(x)=\overline{s(-x)}\), set \(M_\pm(f)=\int f(x)e^{\pm x/2}dx\). Then
\[\boxed{P_{\rm poles}(s,t)=M_-(t)\overline{M_+(s)}+M_+(t)\overline{M_-(s)}.}\]
This is the **cross-pole form**, not the sum of two moment squares. It matches the form of `neutralWindowCorrelation_poleCrossTerms` / `neutralLogHilbertCorrelation_poleOperator` in the frozen Lean source. Its sign is not intrinsically positive on arbitrary complex tests.

## 3. Exact Gamma/digamma normalization identity (analytic derivation)

**Lemma (smooth scalar tests).** For \(k\in C_c^\infty(\mathbb R)\), with \(w(x)=e^{-x/2}/(1-e^{-2x})=e^{x/2}/(e^x-e^{-x})\),
\[\boxed{\frac1{2\pi}\int_\mathbb R\widehat k(r)\left(\operatorname{Re}\psi(\tfrac14+\tfrac{ir}{2})-\log\pi\right)dr=-(\log(4\pi)+\gamma)k(0)-\int_0^\infty\left(k(x)+k(-x)-2e^{-x/2}k(0)\right)w(x)dx.}\]
**Proof.** For \(\Re z>0\), the standard integral difference identity for the digamma function gives
\[\operatorname{Re}\psi(\tfrac14+\tfrac{ir}{2})-\psi(\tfrac14)=\int_0^\infty\frac{e^{-u/4}}{1-e^{-u}}(1-\cos(ru/2))du.\]
Fourier inversion and \(u=2x\) yield
\[A(k)=(\psi(\tfrac14)-\log\pi)k(0)+\int_0^\infty[2k(0)-k(x)-k(-x)]w(x)dx.\]
Rewrite \(2k(0)=2e^{-x/2}k(0)+2(1-e^{-x/2})k(0)\) and evaluate the absolutely convergent correction integral:
\[\begin{aligned}
I&=\int_0^\infty 2(1-e^{-x/2})w(x)dx\\
&=4\int_0^1\frac{dt}{(1+t)(1+t^2)}=\log2+\frac\pi2.
\end{aligned}\]
The special value \(\psi(1/4)=-\gamma-\pi/2-3\log2\) (also quoted by Suzuki §2.2) gives \(\psi(1/4)-\log\pi+I=-\gamma-\log(4\pi)\). Substitution proves the identity. Smooth compact support justifies the Fourier inversion and Fubini steps by the decay of \(\widehat k\); the local cancellation at \(u=0\) is quadratic in \(u\) after Fourier inversion.

**Result:** All three scalar terms of Suzuki's \(W\) agree exactly with Zeta23 `EF.literatureRHS` on \(C_c^\infty\); in particular, **there is no missing Gamma constant**. This is a derivation of a classical explicit-formula equivalence, not a newly discovered theorem.

## 4. Exact mixed-convention dictionary and Fourier scaling

Under Weil-Lab's anti-linear-first correlation convention,
\[C_{\rm lab}(s,t;\tau)=\int t(y)\overline{s(y-\tau)}dy=(t*\widetilde s)(\tau).\]
Suzuki uses the **linear-first** pairing \(Q^{\rm Su}(t,s)=W_{\rm Su}(t*\widetilde s)\). Thus the slot map is \(Q^{\rm Su}(t,s)=R_{\rm Z23}(t*\widetilde s)\), which matches the laboratory's literal **pole − prime + Gamma** expression when the lab expression is defined on the same physical functions. This is a slot **reversal**, not a complex conjugation error. For real spectral \(r\),
\[\widehat{t*\widetilde s}(r)=\widehat t(r)\overline{\widehat s(r)}.\]
The project’s Mathlib transform uses \(\mathcal F f(\xi)=\int f(x)e^{-2\pi i\xi x}dx\), so \(\widehat f(r)=\mathcal F f(-r/(2\pi))\). Therefore \(dr/(2\pi)\mapsto d\xi\), and the Zeta23 Gamma integral becomes \(\int B(-2\pi\xi)\mathcal Ft(\xi)\overline{\mathcal Fs(\xi)}d\xi\), matching the actual **\(-2\pi\xi\)** argument in Weil-Lab's `neutralNF600LiteralArch_eq`. No undocumented \(2\pi\) factor is needed.

## 5. Operator identity and the precise limitation

Suzuki's §8.4 (8.7) already establishes for \(f\in C_c^\infty(-a,a)\):
\[W_{\rm Su}(f*\widetilde f)=\langle G_aDf,Df\rangle_{L^2}=\langle B_af,f\rangle_{L^2},\qquad D=i\tfrac d{dx}.\]
His paper extends the identity to \(H_0^1(-a,a)\) using norm-continuity, and identifies the Friedrichs extension \(A_a\) and the larger closed Weil form domain. These extensions are **Suzuki's results**. The lab's source-carrier Lean theorem `neutralActualZetaGreenZeroForm_arithmetic` applies to `NeutralActualZetaGreenCoefficients` with their constructed source admissibility; `neutralNF600LiteralArch_eq` applies to original-collar Schwartz tests under support hypotheses. Neither of those checked declarations is a universal lab-side theorem identifying every element of Suzuki's completed form domain. Hence:

**CR-SZ-002A PASS (analytic smooth-core scalar formula)**; **CR-SZ-002B OPEN (Weil-Lab source-domain ↔ Suzuki complete operator/form-domain equivalence and an independently compiled Lean transport theorem)**.

Also, Suzuki's \(G_a\) is built from the **same** Weil explicit formula; comparing its form to the lab expression does **not** create an independent spectral realization or prove positivity.

## 6. Quantitative information now admissible, not yet a positivity conclusion

Take the IV030 real bumps \(\phi,\psi\), each with \(\|\phi\|_2^2=\|\psi\|_2^2=J>0\), disjoint supports inside \((-11/10,11/10)\); for \(a\ge 11/10\) they lie in Suzuki's smooth core. Form the **full** 2×2 Hermitian Gram matrix \(H=\begin{pmatrix}h_{11}&c\\\overline c&h_{22}\end{pmatrix}\), with \(h_{11},h_{22}\in\mathbb R\) the actual full signed Weil diagonal values and \(c\) the full cross-pairing. Its prime contribution to \(c\) is \(-\kappa J\) for \(\kappa=\log3/3\), while **both pole cross-moments and Gamma cross-couplings must be added**. The physical \(L^2\) mass matrix is \(J I_2\). The Rayleigh–Ritz principle therefore gives the exact symbolic **upper bound**
\[\boxed{\lambda_a\le\frac{h_{11}+h_{22}-\sqrt{(h_{11}-h_{22})^2+4|c|^2}}{2J}.}\]
The expression on the right is the lowest eigenvalue of the **full selected** Gram pencil. This is an upper—not lower—bound on the true lowest eigenvalue. A positive selected finite matrix does not imply \(\lambda_a\ge0\). No certified bounds for \(h_{11},h_{22},c\) are supplied here, so no numerical or sign claim follows. A later **CR-SZ-003** can seek interval enclosures and, separately, a genuine complement estimate if lower bounds are desired.

## 7. Falsifiers, classification, next

**Exact falsifiers:** omitted pole, replacing mixed cross-poles by square moments, reversing the first-conjugate lab correlation, writing the wrong prime sign, using \(+2\pi\xi\) in the original mathlib Fourier density without the matching transform reversal, losing the \(1/(2\pi)\) gamma factor, or claiming the selected two-bump upper bound is a global lower bound.

**Classification:** canonical **RECOVERY** of complete scalar signed explicit formula; **new checked normalization dictionary** (analytic derivation, not Lean); **OPEN source-carrier/closed-form identification**; **OPEN** genuine full-form quantitative support-extension result; **NO new independent spectral model**, **NO RH proof**.

**Next exact obligation: CR-SZ-003.** Use the full \(2\times2\) source-pinned Gram coefficients to test phase-sensitive Rayleigh upper bounds and support-dependent coupling thresholds in Suzuki's generalized \((G_a,K_a)\) quotient; specify what quantitative complement estimate would be required to infer any lower bound. Only after completing or falsifying that direction move the separate return program to Morishita.
