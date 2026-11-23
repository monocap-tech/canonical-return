# CR-INT-T1D-M — Source-compatible logarithmic prolate cutoff and the surviving finite Schur gate

**Canonical Return — independent research, 2026-11-23 (America/Los_Angeles).**  
**Status: ANALYTICAL DERIVATION; exact-rational constant audit PASS; not a numerical prolate eigenspectrum certificate, Lean proof, global Weil positivity, or RH.** The source window is **fixed** at \(a=2\). This is a stronger, less extravagant finite-codimension conditioning theorem than [T1D-L](CR-INT-T1D-L.md), but **does not settle the remaining finite block**. The comparison with Andrews–Connes–La Gatta is at the level of a time–band operator, not an identified Sonin carrier or a proved completed-domain intertwiner.

## 1. Original signed source, metric, and prolate operator

The canonical arithmetic input is [Suzuki, *Weil's quadratic form via the screw function*, arXiv:2606.09096v4, §1.1 and §2.4](https://arxiv.org/html/2606.09096v4). Put \(I=(-2,2)\), \(\mathcal D=C_c^\infty(I;\mathbb C)\), and
\[
\mathcal D^0=\{f\in\mathcal D:\ \int_I e^{x/2}f(x)dx=\int_Ie^{-x/2}f(x)dx=0\}.
\]
Let \(q(f,g)=W(f*\widetilde g)\) be the **complete signed** original Weil form, pairing linear in the first slot, and let the *native* mass be \(\langle f,g\rangle=\int_I f\overline g\,dx\). For the frozen positive T1C source \(E=\operatorname{span}\{u_+,u_-\}\subset\mathcal D^0\), let \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\subset L^2(I)\), and \(\mathcal F=(E\oplus P)^\perp\), \(F=\mathcal D\cap\mathcal F\). The rank-at-most-two Schur penalty \(\mathcal R\) of T1D-E/F/I acts on \(\mathcal F\) and satisfies \(\mathcal R\ge0\). Its **previously derived but not independently peer-reviewed** sharper estimate \(\|\mathcal R\|<155/182\) is clearly marked whenever imported below.

For \(B>0\), define the **logarithmic time–band concentration operator** on \(L^2(I)\):
\[
\boxed{(\mathsf T_B f)(x)=\int_{-2}^2\frac{\sin B(x-y)}{\pi(x-y)}f(y)\,dy,\qquad \mathsf T_B=\mathbf1_I\mathscr F^{-1}\mathbf1_{[-B,B]}\mathscr F\mathbf1_I\big|_{L^2(I)}.}\tag{M.1}
\]
The diagonal is \(B/\pi\), by continuous extension. Here \(\mathscr F\) is the unitary **whole-line Fourier transform of the zero extension**. Thus \(0\le\mathsf T_B\le I\), it is positive trace-class, and
\[
\boxed{\langle\mathsf T_B f,f\rangle=\int_{-B}^B|\widehat f(\xi)|^2d\xi,\qquad \operatorname{Tr}\mathsf T_B=\frac{4B}{\pi}.}\tag{M.2}
\]
The \(L^2(I)\) eigenfunctions of this sinc-kernel operator are the usual **logarithmic-coordinate prolate spheroidal eigenfunctions**. They are **not** the semilocal Sonin/prolate vectors of [Andrews–Connes–La Gatta, arXiv:2611.13129v1, §§2.3–2.4, 8.5–8.6](https://arxiv.org/html/2611.13129v1); that paper uses a different additive Fourier carrier and arithmetic transport map.

Fix \(\theta=1/128\). Let \(V_B\subset L^2(I)\) be the finite span of prolate eigenfunctions with concentration eigenvalue **strictly greater** than \(\theta\). Instead of applying the prolate projector to a smooth source (which need not preserve compact support), impose **linear moment constraints** on the source:
\[
\boxed{H_B^F=\{y\in F:\ \langle y,v\rangle=0\quad\forall v\in V_B\}.}\tag{M.3}
\]
This is a source-compatible *smooth* linear space. Since \(V_B\) has finite dimension, it remains infinite-dimensional. Each member is exactly two-Mellin-pole-null and physically orthogonal to the selected space \(E\), and by the spectral theorem
\[
\boxed{\int_{-B}^B|\widehat y(\xi)|^2d\xi\le\theta\|y\|_2^2\qquad (y\in H_B^F).}\tag{M.4}
\]
No nonzero smooth compactly supported source is claimed to be exactly bandlimited or have zero Fourier transform on an open interval.

## 2. A rigorous *near-Shannon* bound on the number of excluded modes

Set \(L=|I|=4\) and \(k_B(t)=\sin(Bt)/(\pi t)\). The exact trace defect is the boundary leakage of the full band projection:
\[
\begin{aligned}
D_B&:=\operatorname{Tr}(\mathsf T_B-\mathsf T_B^2)
 =\int_I\int_{\mathbb R\setminus I}|k_B(x-y)|^2\,dy\,dx\\
&=\frac2{\pi^2}\left(\int_0^4\frac{\sin^2(Bt)}t\,dt
+4\int_4^\infty\frac{\sin^2(Bt)}{t^2}\,dt\right).\tag{M.5}
\end{aligned}
\]
For \(4B\ge1\), split the first integral at \(1/B\), using \(|\sin Bt|\le Bt\) before the split and \(|\sin Bt|\le1\) afterward. The second integral is at most \(1\). Consequently,
\[
\boxed{0\le D_B\le\frac2{\pi^2}\bigl(\log(4B)+\tfrac32\bigr).}\tag{M.6}
\]
If \(\lambda_j\) are the eigenvalues of \(\mathsf T_B\), all in \([0,1]\), and \(m_B=\#\{j:\lambda_j>\theta\}\), then
\[
\begin{aligned}
m_B&=\sum_{\lambda_j>\theta}\lambda_j+\sum_{\lambda_j>\theta}(1-\lambda_j)\\
&\le\operatorname{Tr}\mathsf T_B+\theta^{-1}\operatorname{Tr}(\mathsf T_B-\mathsf T_B^2),
\end{aligned}
\]
so
\[
\boxed{m_B\le\frac{4B}{\pi}+\frac{256}{\pi^2}\bigl(\log(4B)+\tfrac32\bigr).}\tag{M.7}
\]
This improves the crude trace/threshold estimate \(m_B\le(4B/\pi)/\theta\), and improves the **order-one factor** compared to T1D-L's \(N=10B\) polynomial-Mellin-jet constraints. No prolate eigenvalues need to be numerically guessed: the spectral projection is mathematically exact, while (M.7) rigorously bounds its finite rank.

## 3. Sharper original Gamma bound and full arithmetic census

T1D-L used \(m_\Gamma(\xi)>\log|\xi|-8\). A small analytical improvement matters greatly at large frequency. From Suzuki's original gamma distribution,
\[
 m_\Gamma(\xi)=m_\Gamma(0)+2\int_0^\infty(1-\cos\xi t)w(t)dt,\quad
 m_\Gamma(0)=-\gamma-\frac\pi2-3\log2-\log\pi,\quad
 w(t)\ge\frac1{2t}\ (0<t\le1).
\]
The classical cosine-integral identity is
\[
\int_0^1\frac{1-\cos(Xt)}t dt=\gamma+\log X-\operatorname{Ci}(X),\quad
\operatorname{Ci}(X)=-\int_X^\infty\frac{\cos t}{t}dt.
\]
Integration by parts shows \(|\operatorname{Ci}(X)|\le2/X\). Thus for \(X\ge64\), using \(\pi<22/7\), \(\log2<347/500\) and \(\log\pi<6/5\),
\[
\frac\pi2+3\log2+\log\pi+\frac2X
<\frac{11}7+\frac{1041}{500}+\frac65+\frac1{32}<\frac{49}{10}.
\]
Hence
\[
\boxed{m_\Gamma(\xi)>\log|\xi|-\frac{49}{10}\quad(|\xi|\ge64),\qquad m_\Gamma(\xi)>-6\quad(\xi\in\mathbb R).}\tag{M.8}
\]
This uses an **elementary high-frequency bound**, not a fabricated numerical digamma interval. The original signed Weil prime term on \((-2,2)\) retains *all* genuine prime powers \(n\le54\); a new exact-rational 24-term enclosure improves the conservative estimate to
\[
\boxed{S_2=\sum_{n\le54}\frac{\Lambda(n)}{\sqrt n}<\frac{61}{5},\qquad q_{\rm prime}(y,y)>-\frac{122}{5}\|y\|_2^2.}\tag{M.9}
\]
The original pole terms vanish only because \(y\in\mathcal D^0\). No prime/Gamma signs are silently changed, and no arithmetic test is replaced with an approximate model.

## 4. A substantially cheaper infinite-dimensional coercivity certificate

Choose
\[
\boxed{B=2^{44}=17,592,186,044,416,\quad\theta=1/128.}\tag{M.10}
\]
Using \(\log2>693/1000\), the complete original Weil form on \(H_B^F\) satisfies
\[
\begin{aligned}
\frac{q(y,y)}{\|y\|_2^2}
&>\bigl(\log B-\tfrac{49}{10}\bigr)(1-\theta)-6\theta-\tfrac{122}{5}\\
&>\frac{127}{128}\cdot44\cdot\frac{693}{1000}-\frac{293}{10}-\frac{11}{1280}
=\boxed{\frac{15123}{16000}=0.9451875},\qquad 0\ne y\in H_B^F.\tag{M.11}
\end{aligned}
\]
Using \(\pi>3\) and \(\log2<347/500\), the number of excluded prolate modes obeys the **independent finite** bound
\[
\boxed{m_B<\frac{4B}{3}+951.}\tag{M.12}
\]
This is still astronomically large (roughly \(2.35\times10^{13}\) directions). Nevertheless it is a substantial *quantified reduction* from T1D-L's \(10\cdot2^{60}\) jet conditions, with a better-than-0.9 positive margin at much lower frequency. It does **not** imply an efficient algorithm for enumerating the prolate basis at this B.

**Conditional Schur improvement:** import the distinct T1D-I analytical enclosure \(\|\mathcal R\|<155/182\). Then on the authentic smooth subspace \(H_B^F\),
\[
\boxed{q_*(y,y):=q(y,y)-\langle\mathcal R y,y\rangle
>\biggl(\frac{15123}{16000}-\frac{155}{182}\biggr)\|y\|_2^2
=\frac{136193}{1456000}\|y\|_2^2>\frac9{100}\|y\|_2^2.}\tag{M.13}
\]
The coercivity (M.11) **does not depend on T1D-I**; (M.13) explicitly does. Importing only the weaker T1D-F \(\|\mathcal R\|<1785\) is also possible using the same sharper Gamma/prime constants but a much larger \(B=2^{2640}\), as the rational audit checks. This supplies an independent route with respect to T1D-I's specialized residual enclosure, not a proof of either previous residual estimate.

## 5. The exact *remaining* finite-dimensional smooth-source gate

Define \(\mathscr H_B\) to be the **physical \(L^2(I)\) closure** of \(H_B^F\) inside \(\mathcal F\). The restrictions of the finitely many prolate functionals \(f\mapsto\langle f,v_j\rangle\) to \(F\) have finite rank \(m'\le m_B\). Choose \(m'\) fixed smooth tests \(w_1,\ldots,w_{m'}\in F\) so their moment vectors form a basis of that image. Elementary finite-dimensional linear algebra then yields the **algebraic smooth-source splitting**
\[
\boxed{F=H_B^F\oplus W_B,\qquad W_B=\operatorname{span}\{w_1,\ldots,w_{m'}\},\quad m'\le m_B.}\tag{M.14}
\]
The split is not asserted to be physically orthogonal: the actual prolate eigenfunctions themselves are analytic on the window and are generally **not compactly supported smooth Weil tests**. In particular, projecting a smooth \(f\) onto \(V_B\) and calling that projected vector a Weil test would be a domain error.

For any smooth \(w\in W_B\), the original full signed Weil image \(\mathcal K_2w\in L^2(I)\) from T1D-E is defined on the smooth pole-null core. Define the **cross-block representer**
\[
\boxed{h_w:=\Pi_{\mathscr H_B}\bigl(\mathcal K_2w-\mathcal Rw\bigr)\in\mathscr H_B.}\tag{M.15}
\]
Then \(q_*(w,z)=\langle h_w,z\rangle\) for \(z\in H_B^F\). By (M.13), putting \(d_*=136193/1456000\), a Cauchy–Schwarz completion gives
\[
\boxed{q_*(w+z,w+z)
\ge q_*(w,w)-\frac{\|h_w\|_2^2}{d_*}\qquad(w\in W_B,\ z\in H_B^F),}\tag{M.16}
\]
The non-strict inequality is stated to include the zero-vector and degenerate equality cases. Therefore **a sufficient, explicitly falsifiable remaining finite condition** is
\[
\boxed{q_*(w,w)\ge\frac1{d_*}\|h_w\|_2^2\quad\forall w\in W_B.}\tag{M.17}
\]
It is equivalent to positivity of a finite Hermitian matrix constructed from the true source quantities
\[
\boxed{\mathsf S_{ij}=q_*(w_i,w_j)
-\frac1{d_*}\langle h_{w_i},h_{w_j}\rangle.}\tag{M.18}
\]
**No numerical or analytic sign certification of this enormous matrix is claimed.** The criterion is **sufficient, not necessary**, because the scalar bound \(d_*\|z\|^2\) sacrifices information about the exact high-mode form. A negative eigenvalue of this conservative matrix would **not** by itself be a negative Weil direction. If it were positive semidefinite, however, it would establish \(q_*\ge0\) on all smooth pole-null \(F\), and the *existing selected-block Schur algebra* would then give \(q\ge0\) on \(\mathcal D^0\). The two independent original pole directions would still require a further genuine smooth-source extension, and no support-independent RH result would follow.

## 6. Comparison, limitations, and next gate

The dictionary \(F(e^x)=f(x)\) preserves \(du/u=dx\) and identifies log-frequency \(\xi\) with the Mellin spectral variable \(s=i\xi\). The operator \(\mathsf T_B\) is the **log-coordinate** time–band concentration operator; its eigenfunctions give a legitimate *source-domain conditioning constraint* but **not** an identification with the arithmetic \(\mathcal E\)-transferred Sonin space, the semilocal \(H2'\) conditioning of Andrews–Connes–La Gatta, or Suzuki's completed screw-kernel domains. Their different operators and metrics are not interchanged here.

**Established analytically (subject to external mathematical review):** exact \(\mathsf T_B\) low-band identity, \(O(B+\log B)\) excluded-mode count, improved Gamma/prime lower bound, and \(q>15123/16000\) on an explicitly defined infinite-dimensional smooth pole-null physical complement portion. **Conditioned on the separately stated T1D-I norm claim:** \(q_*>136193/1456000\) there. **Identified but not satisfied:** finite smooth-correction gate (M.17). **Still open:** its actual low-mode Gram and mixed couplings, positivity on all \(F\), the two original pole directions, support-uniform estimates, Suzuki–Connes/Sonin domain equivalence, Lean acceptance, and RH.

**Reproduce:** `python3 src/audits/canonical_return_t1d_m_prolate_cutoff.py` checks all 24 genuine prime powers, outward rational \(S_2\) enclosure, \(\log2\) constants, Gamma high-frequency constant, sinc trace-defect **numerical bound**, prolate-count arithmetic, and strict coercivity/Schur margins. The spectral theorem, \(\mathrm{Ci}\) identity, sinc trace boundary formula, and smooth-core finite-section existence are **human analytical proofs**, not formalized by the script. There is no fake certified prolate eigenspectrum or GitHub CI claim.

**Next CR-INT-T1D-N:** find a *constructive, finite-cost* smooth correction system for prolate moments and test whether the resulting finite \(\mathsf S\) can be controlled using rigorous off-diagonal decay or trace/Hilbert–Schmidt bounds. Specifically, seek a bound that does **not** scale like the entire \(m_B\sim4B/\pi\) dimension. Do not claim that removing finitely many low-band directions or positivity on \(H_B^F\) settles the complement unless (M.17) or an equally strong exact cross-block condition is actually proved.