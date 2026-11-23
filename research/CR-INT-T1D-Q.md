# CR-INT-T1D-Q — A rational positive supersolution for the *joint* prime operator; the prolate gate falls to 2²¹

**Canonical Return · November 23, 2026 · independent analytical source investigation.**  
**Disposition:** new **source-valid analytic operator bound**, with an executable all-rational path-event audit and explicit imported hypotheses for the Gamma/prolate consequences. **Not** a calculation of the omitted finite-dimensional low-mode Schur matrix, the entire original Weil complement, a completed Suzuki–Connes operator correspondence, Lean certification, or RH. Draft PR #2 remains unmerged pending independent mathematical review.

## 1. Exact source, metric, and question

Keep *unchanged* Suzuki's complete signed Weil functional ([arXiv:2606.09096v4 §1.1](https://arxiv.org/html/2606.09096v4)), the physical Hilbert space \(H=L^2(I,dx)\) on \(I=(-2,2)\), and the smooth doubly pole-null tests

\[
\mathcal D^0=\{f\in C_c^\infty(I;\mathbb C):\ \int f(x)e^{x/2}\,dx=\int f(x)e^{-x/2}\,dx=0\}.
\]

For the 24 genuine prime powers \(n\le54\), put \(\kappa_n=\Lambda(n)/\sqrt n>0\), \(T_tf(x)=f(x-t)\), and let \(P_I\) denote compression/zero extension to \(I\). The **auxiliary positive-kernel (but not positive-semidefinite) operator**

\[
\boxed{A=\sum_{n\le54}\kappa_n P_I(T_{\log n}+T_{-\log n})P_I\big|_H,\qquad q_{\mathrm{prime}}(f,f)=-\langle Af,f\rangle_H.}\tag{Q.1}
\]

It is bounded and self-adjoint and preserves pointwise nonnegativity. **No arithmetic sign is changed in the Weil form.** T1D-O established \(\|A\|<55/4\) through individual-translation norms; [T1D-P](CR-INT-T1D-P.md) improved this to \(\|A\|<93/10\) by a complete three-step path sweep. The next obligation was to exploit the path structure *jointly* without assuming cancellations or discretizing the physical interval.

## 2. A three-step **positive supersolution**, not a power-norm estimate

Let \({\bf1}\) denote the constant-one **auxiliary** vector in \(L^2(I)\), not a smooth admissible Weil test. Define the positive, bounded, measurable weight

\[
\boxed{h={\bf1}+20A{\bf1}+8A^2{\bf1}\ge1.}\tag{Q.2}
\]

All three terms are explicitly computable nonnegative step functions of the logarithmic position \(x\). Since \(A\) is a positive-kernel operator,

\[
Ah=A{\bf1}+20A^2{\bf1}+8A^3{\bf1}.
\]

**Theorem (strict rational positive-supersolution bound).** On almost every \(x\in(-2,2)\),

\[
\boxed{(Ah)(x)<\frac{203}{25}\,h(x).}\tag{Q.3}
\]

The proof is fully finite and source-specific; it does **not** infer (Q.3) from approximate floating-point samples:

1. Represent each signed shift by multiplication by \(n\) or \(1/n\) in \(z=e^{x+2}\in(1,e^4)\). Every path \(r_0=1,r_j=\prod_{i=1}^j n_i^{\varepsilon_i}\) is included only if *each intermediate step* remains in the original window. Its admissible \(z\)-interval is \((1/\min r_j,\,e^4/\max r_j)\).
2. For all lengths \(k=1,2,3\), enumerate all **48** original signed shifts without omitting prime powers. The exact admissible counts are \(48\), \(1,290\), and \(28,160\). There are \(6,054\) distinct step-function endpoints in their union.
3. Enclose each \(\kappa_n\) above and below by rational values of denominator \(10^{12}\), using explicit convergent log-series remainders and integer-square-root bounds. Enclose \(e^4\) by rational Taylor bounds. All proposed endpoint ordering is **separately certified by exact rational separation**; float is used only to propose a sort permutation.
4. On every open interval between adjacent endpoints, let \(L_k\le A^k{\bf1}\le U_k\) denote the exact-rational *lower* and *upper* path-weight sweeps. Verify the **integer inequality**

\[
\boxed{25\,(U_1+20U_2+8U_3)<203\,(1+20L_1+8L_2).}\tag{Q.4}
\]

Here the displayed quantities are rational values represented internally by integer numerators with exact powers of \(10^{12}\) denominators. The inequality is checked on every one of the \(6,053\) nondegenerate endpoint intervals. Endpoints form a measure-zero set and are irrelevant to the \(L^2\) operator result. The observed largest upper/lower ratio is about \(8.117878245374\), **a noncertifying diagnostic only**; the *strict rational* assertion (Q.4) is the certificate.

The auxiliary \(h\) need not belong to \(\mathcal D^0\). It is a bounded positive Schur-test weight, and using it does **not** enlarge or alter the Weil test class.

### Operator-norm deduction

For any symmetric nonnegative kernel or finite symmetric sum of weighted translations, the **weighted Schur test** states

\[
Kh\le\lambda h,\quad h>0\quad\Longrightarrow\quad\|K\|_{L^2\to L^2}\le\lambda.
\]

Indeed, use weighted Cauchy–Schwarz \(|Kf(x)|^2\le(Kh)(x)\,(K(|f|^2/h))(x)\), integrate, use symmetry of \(K\), then use \(Kh\le\lambda h\) twice. Since the step ratios in (Q.3) range over finitely many intervals with a **strict** finite maximum, they are bounded by a constant \(\lambda_*<203/25\). Thus

\[
\boxed{\|A\|_{L^2(I)\to L^2(I)}<\frac{203}{25}=8.12.}\tag{Q.5}
\]

This is a genuine improvement over \(9.3\) in T1D-P, derived from the **original joint compressed prime translation operator**. The arithmetic contribution to the *signed* Weil form consequently satisfies \(q_{\mathrm{prime}}(f,f)>-(203/25)\|f\|_2^2\) for nonzero smooth tests, without asserting that \(-A\) is positive.

## 3. Consequence: an **independently source-typed**, smaller high-frequency gate

**Imported, not reproved here:** the T1D-M/O Gamma bounds on the exact same original Weil source,
\[
m_\Gamma(\xi)>\log|\xi|-49/10\ (|\xi|\ge64),\qquad m_\Gamma(\xi)>-6\ (\xi\in\mathbb R),
\]
and the prolate concentration lemma for \(\mathsf T_B=P_I\mathcal F^{-1}\mathbf1_{[-B,B]}\mathcal F P_I\). Let \(E=\operatorname{span}\{u_+,u_-\}\subset\mathcal D^0\) and \(F=\mathcal D^0\cap E^\perp\). Let \(H_B^F\subset F\) be the smooth tests physically orthogonal to all prolate eigenfunctions with concentration eigenvalue \(>\theta=1/128\). The resulting \(H_B^F\) remains **infinite-dimensional**, is exactly pole-null, and has low-Fourier-band mass \(\le\theta\|y\|_2^2\). This is a **subspace of genuine smooth Weil tests**, not an assertion that the prolate eigenvectors themselves are allowable tests.

At the newly sufficient cutoff
\[
\boxed{B=2^{21}=2,097,152,}\tag{Q.6}
\]
combining (Q.5) with the **original** complete Gamma and arithmetic terms gives
\[
\begin{aligned}
q(y,y)&>\left[\frac{127}{128}\left(21\log2-\frac{49}{10}\right)-\frac6{128}-\frac{203}{25}\right]\|y\|_2^2\\
&>\boxed{\frac{180571}{128000}\|y\|_2^2>1.4107\|y\|_2^2},\qquad y\in H_B^F\setminus\{0\}.
\end{aligned}\tag{Q.7}
\]
The last step uses the **rationally proved** \(\log2>693/1000\). This is a positive estimate on an explicitly conditioned infinite-dimensional subspace of the genuine smooth two-pole-null core. It does **not** establish positivity on all of \(F\), let alone on the complete smooth domain with the pole carriers.

**Conditional finite-rank Schur transfer.** If T1D-I's separate analytical estimate \(\|\mathcal R\|<155/182\) is correct, then the original form-level complement Schur inequality holds on this conditioned subspace with the margin
\[
\boxed{q(y,y)-\langle\mathcal R y,y\rangle>
\frac{6511961}{11648000}\|y\|_2^2>\frac12\|y\|_2^2.}\tag{Q.8}
\]
The new norm theorem (Q.5) and Weil lower bound (Q.7) do **not** depend on T1D-I; only (Q.8) imports that separate residual-norm result. These are *strict sufficient inequalities* on \(H_B^F\), not a theorem about the omitted directions.

## 4. The remaining rank barrier still counts in the millions

Using T1D-N's trace/trace-defect rank theorem (imported as a separate analytic input), the \(B=2^{21}\), \(\theta=1/128\) excluded prolate rank has the exact rational bracket
\[
\boxed{2,669,099\le m_{2^{21}}\le2,796,699.}\tag{Q.9}
\]
After four known finite physical constraints (the two Mellin pole carriers and the two selected directions), the new constraints needed inside the physical complement are still at least \(2,669,095\). This is about a twofold further reduction from T1D-P's approximately 5.3–5.6 million directions, **not a solution** to the excluded low-mode Schur block.

A useful **methodological boundary** is visible from the same constants: at \(B=2^{20}\), even setting \(\theta=0\) in the *particular conservative combined lower-bound formula* gives at most \(20\log2-49/10-203/25<155/182\). The rational upper enclosure \(\log2<1733/2500\) verifies that inequality. Thus the same crude Gamma norm constants and the same uniform Schur domination test cannot certify the \(2^{20}\) cutoff. **This is a limitation of this certificate, not a no-go theorem for other form-adapted arguments.**

The persistent million-dimensional excluded sector motivates a different research direction: instead of repeatedly reducing a high-pass cutoff, build a form estimate using the **joint Gamma-plus-prime operator on low-frequency test directions**. A source-preserving comparison to the conditioned semilocal Sonin/prolate construction of Andrews–Connes–La Gatta ([arXiv:2611.13129v1 §8.6](https://arxiv.org/html/2611.13129v1)) still needs an explicit arithmetic transport, metric dictionary, and completed-domain intertwiner; nothing here supplies those missing hypotheses.

## 5. Reproduction, limitations, and next transfer gate

Execute `python3 src/audits/canonical_return_t1d_q_positive_supersolution.py` (standard Python library only). It audits the 24-prime-power list, \(e^4\) enclosure, rigorous log/square-root weight bounds, **48 + 1,290 + 28,160 path census**, **6,054 exact ordered endpoints**, the *two-sided* event inequality (Q.4), the finite rational high-band margin and Schur arithmetic, and the prior trace-defect rank substitution. The weighted Schur-test argument, exact Gamma multiplier bounds, and prolate concentration theorem are **human mathematical arguments, not compiler or symbolic outputs**. The audit does **not** calculate the actual joint spectrum \(\sigma(A)\), certify the missing low-mode matrix, or prove RH.

**Disposition:** original-joint-prime norm \(\|A\|<203/25\), **source-level analytic** with executable exact-rational interval certificate; conditioned high-band Weil positivity at \(B=2^{21}\), **analytic conditional on the pinned prior Gamma/prolate lemmas**; conditioned Schur margin, **additionally conditional on T1D-I**. **OPEN:** full original smooth complement \(F\), the excluded finite smooth low-mode matrix and its Weil cross terms, independent review of all imported inequalities, Suzuki/Connes semilocal domain and metric preservation, RH.

**Next CR-INT-T1D-R:** stop treating the low-mode sector as a mere rank statistic. Seek a falsifiable *joint source form bound* on explicit low-frequency, smooth pole-null test families, or a resolvent/Dirichlet comparison controlling that block without enumerating millions of prolate modes. If a spectral transport is proposed, test preservation of **the complete original signed Weil form and its native \(L^2\) metric** before any Phase IV handoff.
