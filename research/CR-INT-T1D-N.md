# CR-INT-T1D-N — The intrinsic prolate rank barrier and the limits of smooth low-mode correction

**Canonical Return | independent mathematical research | 2026-11-23 America/Los_Angeles.**  
**Status: ANALYTIC OPERATOR-THEORETIC DERIVATION, independently checked exact-rational numerical consequences. Not a new Weil positivity, an actual prolate eigenvalue enumeration, a completed Suzuki–Connes transfer, Lean certification, or RH.** This note follows [T1D-M](CR-INT-T1D-M.md) and addresses its proposed finite smooth-correction gate. Its principal result is an obstruction *to reducing the number of moment conditions while retaining the same uniform low-band leakage guarantee*; it is **not** a no-go theorem for other spectral/coercivity methods.

## 1. Exact carrier, conditioning and arithmetic custody

Fix the **same physical source** as T1D-M: \(I=(-2,2)\), \(\mathcal H=L^2(I,dx)\), \(\mathcal D=C_c^\infty(I;\mathbb C)\), and
\[
\mathcal D^0=\{f\in\mathcal D:\int_I f(x)e^{x/2}dx=\int_I f(x)e^{-x/2}dx=0\}.
\]
The original full signed Weil form \(q\) and the native \(L^2\) metric are unchanged; no claim about a new Weil operator, form-domain completion, or Connes semilocal equivalence is needed for the present rank theorem. Let \(E=\operatorname{span}\{u_+,u_-\}\) be T1C's two physically orthonormal, smooth, doubly pole-null sources, \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\subset\mathcal H\), and \(\mathcal F=(E\oplus P)^\perp\). Then \(\operatorname{codim}_{\mathcal H}\mathcal F=4\). Its original smooth core is \(F=\mathcal D\cap\mathcal F\), dense in \(\mathcal F\).

For \(B>0\), define the ordinary **logarithmic-coordinate** time–band operator
\[
\mathsf T_B=\mathbf1_I\mathscr F^{-1}\mathbf1_{[-B,B]}\mathscr F\mathbf1_I\big|_{L^2(I)},\quad k_B(t)=\frac{\sin Bt}{\pi t},\quad 0\le\mathsf T_B\le1.
\]
Thus \(\langle\mathsf T_Bf,f\rangle=\int_{-B}^B|\widehat f(\xi)|^2d\xi\). This is **not** the semilocal uncertainty operator of Andrews–Connes–La Gatta, arXiv:2611.13129v1. The original signed Weil/Gamma and all 24 genuine prime powers remain as in T1D-M; the present theorem concerns the *conditioning geometry*, not a second computation of the Weil form.

Let \(\theta=1/128\) and \(V_B^>=\operatorname{Ran}\mathbf1_{(\theta,1]}(\mathsf T_B)\). Set \(m_B=\dim V_B^>\). T1D-M conditions its smooth source by orthogonality to this whole space, thereby ensuring \(\langle\mathsf T_By,y\rangle\le\theta\|y\|^2\).

## 2. Two-sided prolate count: the rank is near the Shannon number

The exact trace and trace defect, using the full-line orthogonal band projector, are
\[
\operatorname{Tr}\mathsf T_B=\frac{4B}{\pi},\qquad
D_B=\operatorname{Tr}(\mathsf T_B-\mathsf T_B^2)
=\frac2{\pi^2}\left(\int_0^4\frac{\sin^2 Bt}{t}dt+4\int_4^\infty\frac{\sin^2 Bt}{t^2}dt\right).
\]
For \(4B\ge1\), splitting the first integral at \(1/B\) yields the elementary **analytic** bound
\[
0\le D_B\le\frac2{\pi^2}\left(\log(4B)+\frac32\right).
\]
For every eigenvalue \(0\le\lambda\le1\), write
\[
\operatorname{Tr}\mathsf T_B-m_B
=-\sum_{\lambda>\theta}(1-\lambda)+\sum_{\lambda\le\theta}\lambda.
\]
On the first set \(1-\lambda\le\lambda(1-\lambda)/\theta\); on the second \(\lambda\le\lambda(1-\lambda)/(1-\theta)\). Both sums are absolutely convergent and separately bounded by the trace defect. Therefore the new **two-sided** estimate is
\[
\boxed{\frac{4B}{\pi}-\frac{D_B}{1-\theta}\le m_B\le\frac{4B}{\pi}+\frac{D_B}{\theta}.}\tag{N.1}
\]
In particular,
\[
\boxed{\frac{4B}{\pi}-\frac{256}{127\pi^2}(\log(4B)+3/2)
\le m_B\le\frac{4B}{\pi}+\frac{256}{\pi^2}(\log(4B)+3/2).}\tag{N.2}
\]
Thus \(m_B=4B/\pi+O(\log B)\), **not** merely \(O(B)\). This improves T1D-M's upper bound by an independent lower bound. It uses no numerical prolate eigenvalues and no RH-related hypothesis.

At the **already frozen** \(B=2^{44}\), the rational bounds \(3<\pi<22/7\), \(\log2<7/10\) give \(D_B<337/45\) and \(D_B/(1-\theta)<8\). Hence
\[
\boxed{m_{2^{44}}>\frac{14\cdot2^{44}}{11}-8,\qquad m_{2^{44}}\ge22{,}390{,}054{,}965{,}613.}\tag{N.3}
\]
For comparison, the earlier certified \(\log2<347/500\) gives \(D_B/\theta<951\) and
\(m_{2^{44}}<4\cdot2^{44}/3+951\), recovering T1D-M's upper estimate.

## 3. Optimality of spectral conditioning under the precise leakage requirement

**Min–max theorem.** Suppose \(X\subset\mathcal H\) is a closed linear subspace with
\[
\langle\mathsf T_Bx,x\rangle\le\theta\|x\|^2\qquad\forall x\in X.
\]
Then \(\operatorname{codim}_{\mathcal H}X\ge m_B\). Indeed, if the codimension were \(<m_B\), dimension counting would give a nonzero \(x\in X\cap V_B^>\), whereas every such \(x\) has strictly larger Rayleigh ratio than \(\theta\). The canonical choice \(X=(V_B^>)^\perp\) attains **exactly** \(m_B\) constraints. Thus the prolate spectral cutoff is rank-optimal among all closed subspaces obeying that *uniform low-frequency energy cap*. This is an exact spectral/min–max result.

In the actual original-Weil complement \(\mathcal F\), only four physical directions have already been removed. Consequently any closed subspace \(X\subset\mathcal F\) with the above leakage guarantee must have
\[
\boxed{\operatorname{codim}_{\mathcal F}X\ge m_B-4.}\tag{N.4}
\]
At \(B=2^{44}\), this gives **at least 22,390,054,965,609 further independent constraints**. The actual T1D-M moment map \(F\to(V_B^>)^*\) has rank between \(m_B-4\) and \(m_B\): the loss of rank is \(\dim(V_B^>\cap\mathcal F^\perp)\le4\). The same statement holds for smooth-source constraints because \(F\) is dense in \(\mathcal F\) and the finite moment functionals are \(L^2\)-continuous.

**Scope of this obstruction:** it rules out eliminating *most* of these constraints while retaining T1D-M's particular guarantee \(\int_{-B}^B|\hat y|^2\le\theta\|y\|^2\). It does **not** rule out a different inequality using the joint arithmetic/Gamma structure, a smaller cutoff or changed frequency threshold, compact operator estimates, or other conditions that do not require this uniform band-leakage cap.

## 4. Parity does not make the rank small

The reflection \(Rf(x)=f(-x)\) commutes with \(\mathsf T_B\). The selected source \(E\) removes one even and one odd physical direction; \(P\) likewise removes the even \(\cosh(x/2)\) and odd \(\sinh(x/2)\) directions. Thus \(\mathcal F\) has codimension **two in each parity sector**. The reflected trace is
\[
\operatorname{Tr}(R\mathsf T_B)=\frac1\pi\int_0^{4B}\frac{\sin u}{u}du=\frac{\operatorname{Si}(4B)}\pi.
\]
Splitting the sine integral at 1 and integrating its tail by parts yields \(|\operatorname{Si}(x)|\le4\) for all \(x>0\). Therefore \(|\operatorname{Tr}(R\mathsf T_B)|<4/3\). Applying the same trace-defect argument separately to the even and odd compressions gives, at \(B=2^{44}\),
\[
\boxed{m_{B,\mathrm{even}}>\frac{7B}{11}-9,
\qquad m_{B,\mathrm{odd}}>\frac{7B}{11}-9.}\tag{N.5}
\]
So each parity sector has **at least 11,195,027,482,802 excluded directions** and requires at least **11,195,027,482,800 additional constraints** after its two known physical directions are accounted for. Even/odd block diagonalization is worthwhile algebraically, but does not shrink either block to a small matrix.

## 5. Why arbitrary smooth corrections lack dimension-free derivative control

T1D-M correctly constructs a finite (but enormous) smooth correction system \(W_B\subset F\) by density and an invertible finite moment matrix. Such a system exists; the smooth tests are not replaced by nonsmooth prolate eigenvectors. However, **existence alone says nothing about the condition number of that matrix or the energy cost of its lifts**.

There is a quantitative universal obstruction to expecting a physically orthonormal smooth correction basis to have uniformly bounded first derivatives. If \(w_1,\ldots,w_r\in C_c^\infty(-2,2)\) are orthonormal in the native \(L^2\) metric, the min–max principle for the Dirichlet Laplacian on \((-2,2)\), whose eigenvalues are \((j\pi/4)^2\), gives
\[
\boxed{\sum_{j=1}^r\|w_j'\|_2^2\ge
\frac{\pi^2}{96}r(r+1)(2r+1).}\tag{N.6}
\]
In particular, \(\max_j\|w_j'\|_2^2\ge\pi^2(r+1)(2r+1)/96\). For the correction space with \(r\ge m_B-4\), this lower bound grows quadratically in \(B\) for at least one orthonormal basis vector. This says **a dimension-independent derivative bound for all correction basis vectors is impossible**; it does not prohibit a form-specific spectral estimate that avoids handling them individually.

## 6. Practical disposition: change what must be proved, not just how it is enumerated

**New exact result:** the prolate spectral cutoff used in T1D-M is essentially rank-optimal for its prescribed low-frequency leakage. Its remaining smooth correction block has \(\Theta(B)\) independent directions; at \(B=2^{44}\), over \(2.239\times10^{13}\) additional physical constraints are unavoidable within that method. Each parity block is still \(\Theta(B)\). A generic physically orthonormal smooth lift also entails a large derivative-energy cost. All bounds preserve the **actual** logarithmic physical \(L^2\) metric and exactly defined source-domain constraints.

**Not established:** the missing finite Schur matrix is positive, that the Weil complement is nonnegative, any H2′/R1/R2 estimate in the true Connes semilocal carrier, or a valid domain intertwiner with Suzuki's completed screw kernel. None of the new prolate rank results use unproved RH information.

**Next CR-INT-T1D-O:** do *not* try to enumerate \(m_B\sim2.24\times10^{13}\) low modes. Instead seek a new estimate using the **joint original Gamma + signed prime form** on a moderate-frequency space, where the prime operator norm \(2S_2\) is too pessimistic; or formulate a rigorous comparison/transfer with Connes's actual conditioned carrier, preserving both form and physical metric. A *candidate* structured low-frequency bound must be tested on genuine smooth pole-null tests and must dominate the finite-rank Schur penalty before claiming a complement theorem. In particular, no RH or global positivity conclusion follows from this negative rank diagnostic.

**Reproducibility:** `python3 src/audits/canonical_return_t1d_n_rank_barrier.py` checks Archimedean rational approximations, exact integer ranks, parity lower bounds, and the derivative-energy constant. The trace and min–max arguments are human analytical proofs, **not formalized** by that program; no numerical prolate eigenspectrum is claimed. Primary contextual sources: Suzuki, [arXiv:2606.09096v4](https://arxiv.org/html/2606.09096v4); Andrews–Connes–La Gatta, [arXiv:2611.13129v1](https://arxiv.org/html/2611.13129v1).
