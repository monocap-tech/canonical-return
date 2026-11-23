# CR-INT-T1D-P — A certified *joint* prime-translation norm via three-step paths; the prolate gate falls to 2^22

**Canonical Return · November 23, 2026 (America/Los_Angeles).**  
**Status:** a new analytic **bounded-operator norm estimate** for the original arithmetic source on the fixed physical interval, with an exact-rational exhaustive path-sweep audit. The resulting Gamma/prolate and Schur lower estimates **import** the separately stated T1D-M/O analytic hypotheses and T1D-I residual bound, respectively. Not a determination of the missing low-mode Schur matrix; not Lean, a completed Suzuki–Connes identity, global Weil positivity, or RH.

## 1. Freeze the form, metric, domain, and the arithmetic question

Keep the original signed Weil source of Suzuki, [arXiv:2606.09096v4 §1.1 and §2.3–2.4](https://arxiv.org/html/2606.09096v4): the physical interval \(I=(-2,2)\), \(H=L^2(I,dx)\), \(\mathcal D=C_c^\infty(I;\mathbb C)\), and the two-Mellin-pole-null smooth source \(\mathcal D^0\). Write \(\kappa_n=\Lambda(n)/\sqrt n>0\) for genuine prime powers, \(T_tf(x)=f(x-t)\), and \(P_I\) for zero-extension/compression. Let

\[
\boxed{A=\sum_{\substack{n\le54\\ n=p^j}}\kappa_n P_I(T_{\log n}+T_{-\log n})P_I\big|_H,\quad q_{\rm prime}(f,f)=-\langle Af,f\rangle.}\tag{P.1}
\]

There are exactly 24 prime powers in the window, as \(54<e^4<55\). This is the **same** original arithmetic term as T1D-O; it has *not* been replaced by a modified positive sum in the Weil form. The auxiliary \(A\) is a bounded, self-adjoint, **positivity-preserving** operator: if \(f\ge0\), then \(Af\ge0\). This does **not** mean \(A\) is positive semidefinite. Its norm nevertheless controls the original signed prime form:
\(q_{\rm prime}(f,f)\ge-\|A\|\|f\|_2^2\).

T1D-O bounded \(\|A\|<55/4\) using the triangle inequality across the 24 individual compressed shifts. Here we handle **all their overlaps jointly**, without numerical discretization of an arbitrary basis or unproved cancellation of coefficients.

## 2. Three-step paths retain every intermediate boundary crossing

For \(x\in I\) set \(z=e^{x+2}\in(1,Q)\), where \(Q=e^4\). Each of the 48 signed prime-power shifts acts by multiplying \(z\) by \(n\) or \(1/n\). For an ordered three-step path \(\omega=((n_1,\varepsilon_1),(n_2,\varepsilon_2),(n_3,\varepsilon_3))\), set

\[
 r_0=1,\qquad r_j=\prod_{i=1}^j n_i^{\varepsilon_i},\quad
 m_\omega=\min_{0\le j\le3}r_j,\quad
 M_\omega=\max_{0\le j\le3}r_j.
\]

Every \(r_j\) is an **exact positive rational number**. Crucially, an admissible path must stay inside the original interval after **each** of the three translations. It contributes to \((A^3\mathbf1_I)(x)\) precisely when

\[
\boxed{\frac1{m_\omega}<z<\frac{Q}{M_\omega},\qquad Q>\frac{M_\omega}{m_\omega}.}\tag{P.2}
\]

Its positive weight is \(w_\omega=\kappa_{n_1}\kappa_{n_2}\kappa_{n_3}\). Therefore

\[
\boxed{(A^3\mathbf1_I)(\log z-2)=
\sum_{\omega}w_\omega\,
\mathbf1_{(1/m_\omega,\,Q/M_\omega)}(z).}\tag{P.3}
\]

This formula respects the **complete fixed-window boundary truncation**; a path that leaves the support and later re-enters contributes zero, because every intermediate compression in (P.1) is retained.

The resulting step function has **only two types of breakpoints**: rational left endpoints \(1/m_\omega\) and \(e^4\) times rational right endpoints \(Q/M_\omega\). Comparisons between different types reduce to a rational comparison with \(e^4\); no rigorous conclusion depends on sorting approximate floating-point numbers.

## 3. An exact-rational proof that the *joint cubic row sum* is under 800

The accompanying audit exhausts all \(48^3=110,592\) ordered signed three-step candidate paths and counts **28,160 admissible paths** by (P.2). It certifies the weight intervals \(\kappa_n\) from the atanh-series rational enclosure of \(\log p\), integer-square-root rational enclosures of \(\sqrt n\), and outward rounding to a common denominator \(10^{12}\). It independently encloses \(e^4\) by a rational exponential series with tail, at scale \(10^{-24}\).

For each allowed path, the audit adds its **upper-rounded** product weight at its left endpoint and subtracts the same weight at its right endpoint. It merges identical exact endpoints, then sweeps **6,054 distinct breakpoints**. A floating approximation is used solely to propose an event order: *every consecutive comparison in that proposed ordering is independently checked using exact rational inequalities and the certified bracket for \(e^4\).* The computed maximum is a rational integer numerator divided by \(10^{36}\), and the checker asserts the **strict exact** inequality

\[
\boxed{\|A^3\mathbf1_I\|_{L^\infty(I)}<800.}\tag{P.4}
\]

The maximum of the upper-rounded step function is approximately \(791.835664102\), reported only as a reproducibility diagnostic, **not as a decimal interval for the exact norm**. All 24 arithmetic weights and all their three-step interactions are included. The checker neither approximates \(A\) by a mesh operator nor assumes a symmetric positive-semidefinite prime matrix.

### Operator deduction (the step the rational script does *not* by itself prove)

The 48 signed shifts are paired, so \(A^3\) is a self-adjoint positive-*kernel* operator in the sense of preserving pointwise nonnegativity; its kernel need not be positive semidefinite. The unweighted Schur test, or equivalently the \(L^1\)/\(L^\infty\) bounds combined with Riesz–Thorin, gives

\[
\|A^3\|_{2\to2}\le\|A^3\mathbf1_I\|_\infty<800.
\]

Self-adjointness of \(A\) implies \(\|A^3\|=\|A\|^3\). Since \(800<(93/10)^3=804.357\), we obtain the **new source-valid joint norm theorem**

\[
\boxed{\|A\|_{L^2(I)\to L^2(I)}<\frac{93}{10}=9.3.}\tag{P.5}
\]

This is a genuine strict improvement over T1D-O's \(55/4=13.75\), achieved solely by using the original **joint** support-sensitive arithmetic operator. No RH input is required.

## 4. Consequence for the existing time-band conditioning gate

Reuse, without strengthening or claiming independent new formalization, the source-Gamma bounds and prolate leakage theorem stated in T1D-M/O:

\[
m_\Gamma(\xi)>\log|\xi|-49/10\ (|\xi|\ge64),\qquad
m_\Gamma(\xi)>-6\ (\xi\in\mathbb R).
\]

Let \(\mathsf T_B=P_I\mathscr F^{-1}\mathbf1_{[-B,B]}\mathscr FP_I\) and condition the unchanged smooth pole-null complement \(F=\mathcal D^0\cap E^\perp\) against prolate eigenvectors with concentration eigenvalue \(>\theta=1/128\). The remaining smooth subspace \(H_B^F\) is infinite-dimensional and has Fourier low-band mass at most \(\theta\). By (P.5), for any \(y\in H_B^F\) and \(B=2^{22}\),

\[
\begin{aligned}
q(y,y)
&>\left[\frac{127}{128}\left(22\log2-\frac{49}{10}\right)-\frac6{128}-\frac{93}{10}\right]\|y\|_2^2\\
&>\boxed{\frac{58771}{64000}\|y\|_2^2>0.918\|y\|_2^2.}
\end{aligned}\tag{P.6}
\]

The rational value uses the explicit lower bound \(\log2>693/1000\), so **it is a strict analytically conservative margin**. It is not a claim of positivity on all of \(F\). The new threshold is \(2^7=128\) times smaller than T1D-O's \(2^{29}\).

**Conditional Schur comparison.** If the separately documented T1D-I finite-width residual estimate \(\|\mathcal R\|<155/182\) is correct, then the exact imported form-Schur identity yields

\[
\boxed{q(y,y)-\langle\mathcal R y,y\rangle>
\frac{388161}{5824000}\|y\|_2^2>\frac1{16}\|y\|_2^2,
\quad y\in H_B^F.}\tag{P.7}
\]

The new operator inequality (P.5) does **not** depend on T1D-I; only the Schur conclusion (P.7) depends on it. Nor does the positive high-band margin control the omitted low-band mixed matrix.

## 5. The remaining rank barrier is smaller, not gone

Import the same analytic prolate trace/trace-defect rank theorem of T1D-N, with \(B=2^{22}\):

\[
\boxed{5,338,202\le m_{2^{22}}\le5,592,921.}\tag{P.8}
\]

Compared with T1D-O, the upper rank bound falls from \(715,828,537\) to \(5,592,921\) (more than a 127-fold reduction); nevertheless, over **five million** prolate directions remain excluded. A finite-dimensional Schur condition on a **smooth corrected low-mode space** still has to be established. The raw windowed prolate eigenfunctions are not automatically allowable smooth compact tests, and an orthogonal projection onto nonsmooth pole modes does not preserve the test core; the previous domain caveats remain operative.

There is no reason to infer that the entire physical \(F\) is positive from (P.6), nor that the Connes Sonin/prolate subspaces are identified with the logarithmic \(\mathsf T_B\). The Andrews–Connes–La Gatta source, [arXiv:2611.13129v1 §8.6](https://arxiv.org/html/2611.13129v1), retains its distinct arithmetic transport, semilocal hypotheses, metric and completion.

## 6. Verification contract, open requirements, and next falsifier

**New independent analytical content:** the exact three-step joint path representation (P.2)–(P.3), together with a new positive-kernel Schur-test proof of the strict operator norm \(\|A\|<93/10\), and the resulting reduction of the formerly available prolate cutoff from \(2^{29}\) to \(2^{22}\). **Executed exact-rational audit:** prime-power census, rigorous outward arithmetic weights and \(e^4\) bounds, all candidate signed paths, exact interval endpoint ordering, upper-weight event sweep \(<800\), margin fractions and rank arithmetic. The audit does *not* verify the functional-analytic Schur test, Gamma multiplier, concentration theorem, or form-domain transfer; these must remain explicit mathematical inputs.

**OPEN:** exact norm/spectrum of the joint compressed arithmetic operator; further cutoff improvement by higher path powers or a form-adapted Gamma–prime estimate; an effective positive result on the excluded low-mode Schur block; complete original pole directions and large-support limit; completed Suzuki–Connes correspondence; global Weil positivity and RH.

**Next CR-INT-T1D-Q:** exploit higher-step path dynamics via a **rigorous positive supersolution** or an interval-residue transfer operator, aiming to reduce the remaining five-million-dimensional low block. More importantly, construct a form-adapted estimate on that *actual* low-frequency block rather than continuing to rely exclusively on ever-higher-frequency exclusions. Do not claim equivalence to semilocal Sonin conditioning unless both metrics and domains intertwine and the exact signed Weil form is preserved.

**Reproduction:** `python3 src/audits/canonical_return_t1d_p_joint_prime_cubic.py` (Python standard library only). Earlier T1D-M–O audits remain separate for imported analytic hypotheses. Branch remains a draft research review, not a merged mathematical certification.
