# CR-INT-T1D-L — Infinite-dimensional Mellin-frequency conditioning on the original Weil smooth core

**Canonical Return · independent source-level analytical research · 2026-11-23 (America/Los_Angeles).**  
**Status: ANALYTIC DERIVATION + EXACT-RATIONAL FINITE AUDIT PASS.** This constructs *genuinely infinite-dimensional*, finite-codimension smooth source subspaces with a positive lower bound for the **full original signed Weil quadratic form** on a fixed support window. Under explicitly identified **prior** cross-block norm bounds, it also proves the source-defined Schur inequality *on those highly conditioned subspaces*, **not on the entire complement**. The constants and codimensions are deliberately enormous; this is an existence/compatibility diagnostic, not a practical certification of the full ground state, Connes's H2′/R1/R2, global Weil positivity, a completed operator equivalence, Lean acceptance, or RH.

## 1. Source, domain, and the precise conditioning

Freeze the support window \(I=(-2,2)\), the complex smooth source \(\mathcal D=C_c^\infty(I;\mathbb C)\), and the original Weil sesquilinear form \(q(f,g)=W(f*\widetilde g)\) with the physical \(L^2(dx)\) denominator and pairing linear in the first argument. Use the unchanged original two-pole-null domain

\[
\mathcal D^0=\left\{f\in\mathcal D:\int_I e^{\pm x/2}f(x)\,dx=0\right\}.
\]

The canonical inputs are Suzuki, [*Weil's quadratic form via the screw function*, arXiv:2606.09096v4](https://arxiv.org/html/2606.09096v4), §1.1 (the signed distribution) and §2.4 (the Gamma digamma multiplier and original Weil Fourier asymptotics). The comparison source is Andrews–Connes–La Gatta, [arXiv:2611.13129v1](https://arxiv.org/html/2611.13129v1), §2.3–2.4 and §8.5–8.6; **their Sonin and conditioning operators are not identified with the ones defined here**. Prior Canonical Return notes: [T1D-E](CR-INT-T1D-E.md) (source Riesz images), [T1D-F](CR-INT-T1D-F.md) (coarse Schur-norm bound), [T1D-I](CR-INT-T1D-I.md) (sharper, separately derived estimate), and [T1D-K](CR-INT-T1D-K.md) (six-dimensional selected test).

For an integer \(k\ge 60\), set the **finite but very large** integers \(B_k=2^k\) and \(N_k=10 B_k\). Define the **Mellin-jet-null** subspace

\[
\boxed{\mathcal H_k=\left\{y\in\mathcal D^0:\int_{-2}^2 x^j y(x)\,dx=0\quad (j=0,\ldots,N_k-1)\right\}.}\tag{L.1}
\]

This is a complex linear space of **infinite dimension**: it is the kernel of a linear map from the infinite-dimensional space \(C_c^\infty(I)\) into \(\mathbb C^{N_k+2}\). In particular, all its elements are honest smooth, compactly supported original Weil tests and satisfy **both zeta pole conditions**. These restrictions are severe but do not assert any impossible exact Fourier-frequency gap.

To compare the selected positive two-dimensional space \(E=\mathrm{span}\{u_+,u_-\}\) of T1C with its genuine physical complement, set \(F=\mathcal D^0\cap E^\perp\) and \(\boxed{\mathcal H_k^F=\mathcal H_k\cap F}\). These subspaces remain infinite-dimensional because adjoining the two continuous physical orthogonality constraints changes codimension by at most two. All complements and projections in this note use \(L^2(-2,2)\), **not** Suzuki's \(K_a=(-\Delta_N)^{-1}\) and not a Weil-induced seminorm.

## 2. The complete signed Weil form has a directly controlled high-frequency part

Use the **unitary** whole-line Fourier transform of the zero extension of \(y\in\mathcal D^0\), \(\widehat y(\xi)=(2\pi)^{-1/2}\int e^{-i\xi x}y(x)dx\). On \(\mathcal D^0\) the two original pole terms vanish exactly, and the remaining full signed Weil formula is

\[
\boxed{q(y,y)=\int_{\mathbb R}m_\Gamma(\xi)|\widehat y(\xi)|^2d\xi
-\sum_{2\le n\le54}\frac{2\Lambda(n)}{\sqrt n}\int_{\mathbb R}\cos(\xi\log n)|\widehat y(\xi)|^2d\xi,}\tag{L.2}
\]

where \(m_\Gamma(\xi)=\Re\psi_{\rm dg}(1/4+i\xi/2)-\log\pi\) and the sum is over true prime powers only (the others have \(\Lambda=0\)). The cutoff is exact: \(54<e^4<55\). Equation (L.2) is obtained by Fourier transforming the **full** original signed Weil distribution, not by replacing it with the prime-only or Gamma-only form. Note that Suzuki's eq. (2.7) includes an *auxiliary* operator correction: (L.2) is specifically the *original Weil form* from §1.1, with the two pole terms canceled, not a claim that Suzuki's distinct \(B_a\) equals \(q\).

The Gamma kernel has the positive-frequency-difference representation

\[
\begin{aligned}
 m_\Gamma(\xi)-m_\Gamma(0)
 &=2\int_0^\infty(1-\cos(\xi t))w(t)dt,\\
 w(t)&=\frac{e^{-t/2}}{1-e^{-2t}}=\frac{e^{t/2}}{2\sinh t}>0,\\
 m_\Gamma(0)&=-\gamma-\frac\pi2-3\log2-\log\pi>-6.
\end{aligned}\tag{L.3}
\]

The last bound follows from the elementary strict enclosures \(\gamma<3/5\), \(\pi<22/7\), \(\log2<7/10\), \(\log\pi<6/5\). For \(0<t\le1\), \(\sinh t/t\le e^{t^2/2}\le e^{t/2}\) termwise, so \(w(t)\ge1/(2t)\). For \(X\ge1\), integration by parts gives

\[
\int_0^1\frac{1-\cos(Xt)}t dt=\int_0^X\frac{1-\cos s}s ds
\ge\log X-2,
\]

since \(\int_1^X\cos(s)/s\,ds=\sin X/X-\sin1+\int_1^X\sin(s)/s^2\,ds<2\). Therefore **the exact full Gamma multiplier** satisfies

\[
\boxed{m_\Gamma(\xi)>-6\quad(\xi\in\mathbb R),\qquad
m_\Gamma(\xi)>\log|\xi|-8\quad(|\xi|\ge1).}\tag{L.4}
\]

An independent *exact-rational* 24-prime-power census with logarithmic and square-root outward intervals gives

\[
\boxed{S_2=\sum_{2\le n\le54}\Lambda(n)/\sqrt n<13.}\tag{L.5}
\]

Thus the **complete** prime component is bounded below by \(-2S_2\|y\|_2^2>-26\|y\|_2^2\), uniformly for all \(y\) in the smooth source. The numerical \(13\) is a verified conservative bound, **not** a claim that any of the 24 original prime powers is absent.

## 3. Exact finite-jet conditioning implies rigorous low-frequency leakage control

The \(N=N_k\) moment conditions imply, by subtraction of the first \(N\) Taylor terms of \(e^{-i\xi x}\),

\[
\widehat y(\xi)=\frac1{\sqrt{2\pi}}\int_I y(x)
\left(e^{-i\xi x}-\sum_{j=0}^{N-1}\frac{(-i\xi x)^j}{j!}\right)dx.
\]

For real \(\xi x\), the integral Taylor remainder is bounded in magnitude by \(|\xi x|^N/N!\) (the exponential appearing in the integral remainder has modulus one). As \(|x|<2\) and \(\|y\|_1\le 2\|y\|_2\), it follows that

\[
\boxed{\int_{|\xi|\le B}|\widehat y(\xi)|^2d\xi
\le\frac{4B}{\pi}\frac{(2B)^{2N}}{(N!)^2}\|y\|_2^2,\quad B=B_k.}\tag{L.6}
\]

Use \(N!\ge(N/e)^N\), \(e<3\), \(N=10B\), \(\pi>1\), and \(9/25<1/2\). For integer \(B\ge1\), \(B\le2^B\). Hence

\[
\frac{4B}{\pi}\frac{(2B)^{2N}}{(N!)^2}
<4B(9/25)^N<4B\,2^{-10B}
\le4\,2^{-9B}\le\boxed{1/128}.
\tag{L.7}
\]

This is an explicit *uniform* Fourier/Mellin-frequency leakage estimate for an **infinite-dimensional linear subspace**, not for a handful of Galerkin vectors. Exact disappearance of all Fourier mass on an interval is **not claimed**: a nonzero compactly supported smooth function has an entire Fourier transform and cannot vanish on an open band.

## 4. Explicit infinite-dimensional positivity, with two independently scoped Schur comparisons

Let \(\delta_B(y)=\|\mathbf1_{[-B,B]}\widehat y\|_2^2/\|y\|_2^2\). For every \(0\ne y\in\mathcal H_k\), (L.4)–(L.7), the native Plancherel identity, and the full prime estimate yield

\[
\boxed{\frac{q(y,y)}{\|y\|_2^2}
>\bigl(\log B-8\bigr)(1-\delta_B)-6\delta_B-26
\ge\frac{69k}{100}-34-\frac{69k/100-2}{128}.}\tag{L.8}
\]

The last inequality uses \(\log2>69/100\), \(\delta_B\le1/128\), and positivity of the derivative with respect to \(\log B\). The exact rational audit proves the two concrete instances:

\[
\boxed{\begin{array}{ll}
 k=60: & q(y,y)>\dfrac{4539}{640}\|y\|_2^2>7\|y\|_2^2,\\[5pt]
 k=2700:&q(y,y)>\dfrac{232251}{128}\|y\|_2^2>1800\|y\|_2^2.
\end{array}}\tag{L.9}
\]

**This is genuine unconditional-with-respect-to-RH, restricted *infinite-dimensional* source coercivity for the original complete Weil form** on the specified jet-conditioned spaces. It is a finite-codimension phenomenon and is consonant with previously known high-frequency/finite-codimension positivity methods, e.g. Yoshida as discussed by Suzuki. It is *not* a new global spectral lower bound.

For the full original T1C selected block, the **previously derived** positive matrix has \(\mu_\pm>7\) and the physical rank-two Schur operator \(\mathcal R\) defined in T1D-E. Restrict now to \(\mathcal H_k^F\subset F\):

* If one imports **T1D-I's sharper** (not independently peer-reviewed) source result \(\|\mathcal R\|<155/182<1\), then at \(k=60\) the **actual Schur remainder** satisfies \(q(y,y)-\langle\mathcal R y,y\rangle>6\|y\|_2^2\) for all \(0\ne y\in\mathcal H_{60}^F\).
* **Even without relying on T1D-I's sharp residual enclosure**, the substantially weaker **T1D-F** analytical upper bound \(\|\mathcal R\|<1785\) suffices at \(k=2700\):

\[
\boxed{q(y,y)-\langle\mathcal R y,y\rangle>
\left(\frac{232251}{128}-1785\right)\|y\|_2^2
=\frac{3771}{128}\|y\|_2^2>29\|y\|_2^2.}\tag{L.10}
\]

These are honest Schur inequalities **on infinite-dimensional subspaces of the actual smooth physical complement**. They do *not* prove the same inequality on all of \(F\). They also depend on the explicitly labeled T1D-F/I positive selected-block and projected-residual upper bounds; the new high-frequency coercivity (L.9) itself does not depend on those residual bounds.

## 5. What this is, and what it is *not*, in the Suzuki–Connes dictionary

Set \(u=e^x\) and \(F(u)=y(\log u)\). The physical metric transports exactly, \(\|F\|_{L^2(du/u)}=\|y\|_{L^2(dx)}\). The Mellin transform \(\mathcal MF(s)=\int F(u)u^sdu/u\) obeys \(\mathcal MF(\pm1/2)=0\) and

\[
\boxed{\mathcal MF^{(j)}(0)=\int_{-2}^2 x^j y(x)dx=0\quad(0\le j<N).}\tag{L.11}
\]

Thus the new conditioning says the entire Mellin transform has **a zero of order at least \(N\) at spectral center \(s=0\)**, in addition to the two original zeta pole cancellations. This is a precise source-preserving cross-coordinate dictionary. It is not a statement that the Fourier transform vanishes throughout a frequency interval.

The condition is loosely analogous to suppressing the low-frequency concentration measured by a prolate time–band limiting operator on the *logarithmic interval* \((-2,2)\). Such an operator has sinc kernel. The **Connes/Andrews–La Gatta prolate operator**, however, is constructed using the ordinary additive Fourier transform on a different carrier, then transferred through the nontrivial arithmetic map \(\mathcal E\), and its full Sonin space is a two-projection kernel. **Neither that Sonin space nor its uncertainty metric is identified with \(\mathcal H_k\) or \(\mathcal H_k^F\)**. The present finite Mellin-jet constraints are not Connes's H2′ conditioning hypothesis; they do not prove their R1/R2 estimates or any completed Suzuki–Connes intertwiner.

The codimension \(N_k+O(1)=10\cdot2^k+O(1)\) is astronomically large even at \(k=60\). Furthermore, eliminating \(\mathcal H_k^F\) leaves a very large finite-codimension complement with **uncontrolled Weil form and mixed blocks**; the standard Schur completion on \(E\oplus F\) cannot be replaced by positivity on \(\mathcal H_k^F\). Literal windowed polynomial vectors are generally not compactly supported smooth and should **not** be silently inserted as remaining test directions in the Weil form. A genuine finite-dimensional smooth correction or form-domain completion is required.

## 6. Claim inventory, regression, and next falsifier

**New analytical consequences:** (a) a complete Weil Gamma multiplier bound \(m_\Gamma(\xi)>\log|\xi|-8\) at high frequency and \(m_\Gamma>-6\) globally, (b) an independently audited original-window prime-operator bound \(2S_2<26\), (c) a finite-Mellin-jet *linear*, **infinite-dimensional** pole-null smooth source with low-band energy \(<1/128\), (d) source-positive coercivity \(>7\) for \(k=60\) and \(>1800\) for \(k=2700\), and (e) **conditional on stated earlier cross-block norm inputs**, the authentic Schur inequality on the correspondingly conditioned infinite-dimensional portion of \(F\).

**Exact-rational audit**: `src/audits/canonical_return_t1d_l_frequency_tail.py` independently recomputes the 24 true prime powers, exact-rational intervals for \(S_2\), harmonic/digamma baseline inequalities, the factorial and leakage comparison constants, and the displayed coercivity margins. **Analytic inputs**: the published digamma Fourier identity, the original signed Weil distribution, the Taylor remainder, the Gamma-kernel inequality, and the infinite-dimensionality argument are human mathematical proofs, **not** machine formalizations. All claims remain subject to independent mathematical review.

**Falsifiers:** a wrong sign or normalization in Suzuki's Gamma multiplier, missing pole terms, incorrect \(w(t)\) lower bound, misuse of support-limited/Fourier-band-limited simultaneity, failure of the moment-nullity kernel to remain infinite dimensional, or misuse of a source-conditional Schur norm as an independently certified operator theorem. Every numerical inequality in the audit is an outward rational comparison; no fabricated finite-mode compiler acceptance or GH CI is claimed.

**Still OPEN:** positivity/Schur domination on **all** \(F\); cross blocks to the massive finite-codimension low-Mellin-jet residual; two full source pole directions, completed Suzuki/Connes operator equivalence, H2′/R1/R2, support-independent estimates, global Weil positivity and RH.

**Next — CR-INT-T1D-M:** replace the astronomical Taylor-jet threshold by a **source-compatible log-prolate concentration cutoff** with verified eigenvalue bounds and a finite-codimension source correction. Then compute or bound the entire remaining low-frequency block *including its mixed Weil coupling*, rather than treating mere high-frequency positivity as completion of the complement theorem.
