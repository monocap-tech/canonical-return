# CR-SZ-001 — Arithmetic phase fidelity inside Suzuki’s screw kernel

**Status: ANALYTIC CORE IDENTITY ESTABLISHED / FULL FORM AUDIT OPEN / NOT LEAN / NO NEW RH RESULT.** Opened 2026-11-22 PST. This is the first concrete canonical-return analysis, distinct from all current Weil-Lab proof branches.

## I. Exact source conventions (checked against primary references)

**Suzuki**, *Weil’s quadratic form via the screw function*, arXiv:2606.09096**v4** (Nov 15, 2026): §1.1 defines \(Q_W(v_1,v_2)=W(v_1*\widetilde v_2)\), \(\widetilde v(x)=\overline{v(-x)}\). His \(L^2\) inner product is **linear in the first argument** (§1.1); §1.2 eq. (1.3) gives the screw function \(g\); §8.4 eq. (8.7) gives \(Q_W(v)=\langle G_a Dv,Dv\rangle\), with \(D=i\,d/dx\) on \(C_c^\infty(-a,a)\). §8.5 eq. (8.8) gives the generalized quotient relative to \(K_a=(-\Delta_N)^{-1}\) on the mean-zero space.

**Weil-Lab IV030**, original Zeta23 conventions: \(C_{\rm lab}(s,t;\tau)=\int_{\mathbb R}t(x)\overline{s(x-\tau)}dx\), anti-linear **first**. Its positive-named prime interaction \(P_{\rm lab}\) is a sum of \(\kappa_n(C_{\rm lab}(s,t;+\log n)+C_{\rm lab}(s,t;-\log n))\), where \(\kappa_n=\Lambda(n)/\sqrt n\). The Weil functional contains this interaction with a **minus sign**.

**Conjugation dictionary (exact on smooth functions):**

\[C_{\rm lab}(s,t;\tau)=(t*\widetilde s)(\tau).\]

Thus the Suzuki prime component in his linear-first convention satisfies

\[Q_{W,\rm prime}^{\rm Suzuki}(t,s)=-P_{\rm lab}(s,t).\]

The **full** Weil pairing equivalence still requires checking the original lab's pole, Gamma and normalization implementation against Suzuki's complete \(W\); that is a separate gate CR-SZ-002.

## II. Prime powers as derivative jumps of the continuous kernel

Suzuki’s eq. (1.3) contains an explicit prime-power term, for each \(n\ge2\) with \(\Lambda(n)>0\):

\[g_n(h)=\kappa_n\,(|h|-L_n)_+,\qquad L_n=\log n,\quad \kappa_n=\Lambda(n)/\sqrt n.\]

In distributions on the real line,

\[g_n''=\kappa_n(\delta_{L_n}+\delta_{-L_n}).\]

Since Suzuki’s remaining terms are smooth locally at \(h=\pm L_n\ne0\), the one-sided derivative jump of the **full** screw function at \(L_n\) is exactly \(\kappa_n\). The arithmetic prime-power schedule is therefore recoverable as the atomic part of \(g''\) away from the origin. This is **extracted from Suzuki’s existing explicit formula, not claimed as a new theorem**.

A support cell coupling is *possible* only if its difference set meets one of \(\pm L_n\); whether the actual sesquilinear coefficient vanishes depends on amplitudes/phases, not on the graph alone. On \(C_c^\infty(-a,a)\), a prime power with \(L_n\ge2a\) has no quadratic contribution.

## III. Exact smooth-core transport lemma

Let \(f\in C_c^\infty(-a,a)\), \(u=Df=if'\), and let \(G_a^{(n)}\) denote the contribution of \(g_n(x-y)\) to Suzuki's continuous integral operator, compressed to \(L_0^2(-a,a)\). Integration by parts twice has no boundary term and yields:

\[\begin{aligned}
\langle G_a^{(n)}Df,Df\rangle_{L^2}
&=\int\!\!\int g_n(x-y)f'(y)\overline{f'(x)}\,dy\,dx\\
&=-\kappa_n\int f(x-L_n)\overline{f(x)}dx
-\kappa_n\int f(x+L_n)\overline{f(x)}dx\\
&=-\kappa_n[C_{\rm lab}(f,f;+L_n)+C_{\rm lab}(f,f;-L_n)].
\end{aligned}\]

This is an exact **termwise** representation of the original prime summand on the smooth core, and its sign is fixed by Suzuki’s \(g\) together with the two derivatives. This specific finite-test identity does **not** establish an independent spectral realization or extend without proof to all operator/form domains.

## IV. IV030 three-phase calibration inside Suzuki

Use the original verified rational left interval \(I=[-1099/1000,-5494/5000]=[-1.099,-1.0988]\) strictly inside the original left collar. Choose a nonzero real \(\phi\in C_c^\infty(\operatorname{int}I)\); put \(L=\log9\), \(\psi(x)=\phi(x-L)\). Then \(\psi\) sits inside the original right collar, both bumps belong to \(C_c^\infty(-11/10,11/10)\), and \(J=\|\phi\|_2^2>0\). The original IV030 exact-rational source audit excludes all other active prime powers from these bumps. With \(\kappa=\Lambda(9)/3=\log3/3\), define \(f_\beta=\phi+\beta\psi\) for \(\beta\in\{1,-1,i\}\).

On the original Weil-Lab side:
\[P_{\rm lab}(f_\beta,f_\beta)=2\kappa J\operatorname{Re}\beta.\]

On Suzuki’s **n=9 kernel summand** after \(D\):
\[\boxed{\langle G_a^{(9)}Df_\beta,Df_\beta\rangle=-2\kappa J\operatorname{Re}\beta.}\]

The three exact values, ordered \(\beta=1,-1,i\), are \((-2\kappa J,+2\kappa J,0)\). This confirms **relative phase is present in Suzuki’s continuous kernel representation**; a support-graph/mass-only proxy cannot retain it. The \(\beta=i\) case cancels two nonzero complex oriented correlations; it is not a missing graph edge.

## V. Two-dimensional generalized-quotient calibration

In the ordered core basis \((\phi,\psi)\), the **n=9-only** Suzuki numerator represented on \((D\phi,D\psi)\) is:
\[H_9=-\kappa J\begin{pmatrix}0&1\\1&0\end{pmatrix}.\]

For the **correct** denominator \(K_a=(-\Delta_N)^{-1}\), Suzuki §8.5 gives \(\langle K_a Dv,Dw\rangle=\langle v,w\rangle_{L^2}\) on the core, so this selected two-bump subspace has
\[M=J\begin{pmatrix}1&0\\0&1\end{pmatrix}.\]

The generalized **prime-only** pencil \(H_9v=\mu Mv\) therefore has eigenvalues \(\mu=\pm\kappa\). **This is not** the generalized spectrum of full \(G_a\), not Suzuki's full lowest eigenvalue \(\lambda_a\), and not a negative direction for the full Weil form. Gamma, poles, and all remaining terms must be restored before any such conclusion.

## VI. What is proved, what is not, and next gates

**Established analytically here:** exact correlation orientation, negative prime sign, delta-kink coefficient extraction, smooth-core integration-by-parts transport, \(n=9\) phase calibration and the selected prime-only generalized 2×2 pencil. These are verifiable consequences of Suzuki’s existing formulas and the IV030 source, **not a novel RH theorem**. Original IV030's rational interval/log separation checks were previously run as a standalone Python audit; the new dictionary is presented as a mathematical derivation, not as a newly executed Lean proof.

**Not established:** (i) a full literal identification of the two implementations of the signed Weil pairing, including both pole terms, Gamma/digamma and global normalization; (ii) equivalence on completed operator/form domains; (iii) quantitative \(a\mapsto\lambda_a\) bounds from this prime-only subblock; (iv) an unconditional extension to Suzuki’s conditional spectral-convergence theorem; (v) any RH/Weil global positivity result.

**CR-SZ-002 next:** Full signed Weil dictionary, with both first/second argument conventions, normalizations, Gamma and poles checked termwise against Suzuki's exact eq. (1.3), eqs. (2.3), (8.7) and the original Lean source. Must not reuse the n9 prime-only result as if it controlled the full form.

**CR-SZ-003 next:** Support-extension Schur/complement estimates for the exact \((G_a,K_a)\) pencil, with uniform complement control. A finite positive minor alone cannot lower-bound \(\lambda_a\).

## Source custody

1. [Suzuki 2606.09096v4, November 15, 2026](https://arxiv.org/html/2606.09096): §1.1 eqs. (1.1)–(1.3); §2.3 prime shifts; §8.4 eq. (8.7); §8.5 eq. (8.8).

2. [Suzuki 2610.32704, October 2026](https://arxiv.org/abs/2610.32704): Theorem 1.1 and 1.2 explicitly **assume RH**; not input to the unconditional core lemma above.

3. [Weil-Lab IV030 immutable documentation](https://github.com/monocap-tech/weil-lab/blob/1061a331c195e30e85dd8eb35a6eff1c794cf054/docs/IV030_PHASE_IV_ORIGINAL_PRIME_PHASE_FIDELITY_NO_GO.md), exact phase source and original old collars.

4. [Original RH September 2026 handoff](https://github.com/monocap-tech/provenance/blob/main/rh/handoffs/RH_HODGE_PROJECT_HANDOFF_20260914.txt), relevant to eventual transfer, not used to claim that transfer now exists.
