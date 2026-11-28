# CR-INT-T1D-E — The exact original-Weil cross-block representer and finite-rank complement gate

**Canonical Return · independent analytic investigation · 2026-11-22 PST / 2026-11-23 UTC.**  
**Status:** **DERIVED / ANALYTIC** with an executable exact rational/SymPy audit of finite algebra and constants; **NOT** a proof of complement positivity, a completed-domain Suzuki–Connes equivalence, a Lean certification, or RH. The source comparison is confined to the *original signed Weil form on compact smooth logarithmic tests*. Earlier source notes: [T1B](../archive/space/09-t1b.md), [T1C](../archive/space/10-t1c.md), [T1D-D](CR-INT-T1D-D.md).

## 1. The primary sources and precise research question

- **Masatoshi Suzuki**, *Weil's quadratic form via the screw function*, [arXiv:2606.09096v4](https://arxiv.org/html/2606.09096v4), §1.1 original Weil distribution and form; §8.4 eq. (8.7), §8.5 eq. (8.8). The original signed Weil formula is used *directly* below; Suzuki's derivative-kernel equality is a separate source-level representation.
- **Andrews–Connes–La Gatta**, *The uncertainty operator, the Weil ground state and the zeros of zeta*, [arXiv:2611.13129v1](https://arxiv.org/html/2611.13129v1), §8.5 Lemma 8.13 (logarithmic Sobolev lower control); §8.6 Lemma 8.16 (requires complementary positivity), Theorem 8.14 (conditional on H2′), and discussion of unproved R1/R2 and semilocal localization.

**Question:** for a fixed positive two-dimensional pole-null compression, can one define and bound *the actual mixed Weil functionals* on the infinite-dimensional physical orthogonal complement? Can the bound close the Schur complement? The former is answered **yes** on the smooth source at fixed radius, with an explicit Riesz representer and finite-rank cross-block operator. The latter remains **open** because an infinite-dimensional complement lower estimate is missing.

## 2. Domain and two-pole nullity

Fix \(a=2\) and \(I=(-2,2)\). Let \(\mathcal D=C_c^\infty(I;\mathbb C)\) and
\[
\mathcal D^0=\left\{f\in\mathcal D:\int_I f(x)e^{x/2}\,dx=\int_I f(x)e^{-x/2}\,dx=0\right\}.
\]
The physical Hilbert inner product is \(\langle f,g\rangle=\int_I f(x)\overline{g(x)}\,dx\), **linear in the first argument**. Let \(q(f,g)=W(f*\widetilde g)\), where \(W\) includes **both genuine zeta pole terms, the negative prime-power sum and the complete archimedean Gamma distribution**; \(\widetilde g(x)=\overline{g(-x)}\). Denote \(m_\pm(x)=e^{\pm x/2}\mathbf1_I(x)\) and \(P=\operatorname{span}\{m_+,m_-\}\subset L^2(I)\).

Choose the explicit real two-Mellin-pole-null selected functions \(u_\pm\) obtained by physically normalizing the symmetric/antisymmetric combinations of the four-cell family from T1D-D. They form a physical orthonormal basis of \(E=\operatorname{span}_{\mathbb C}\{u_+,u_-\}\subset\mathcal D^0\) with \(q(u_\sigma,u_\tau)=\mu_\sigma\delta_{\sigma\tau}\). The T1D-C theorem gives \(\mu_\sigma>0\) for widths \(<\log(31/27)\); for wider profiles this sign remains a separate hypothesis. Because the pole moments of \(u_\sigma\) vanish, \(E\perp P\) **in the physical metric**.

Define the physical complement \(\mathcal F=(E\oplus P)^\perp\subset L^2(I)\) and its smooth source core \(F=\mathcal D^0\cap E^\perp=\mathcal D\cap\mathcal F\). This avoids treating the nonsmooth windowed pole functions \(m_\pm\) as test functions in \(\mathcal D\). The source \(F\) is dense in \(\mathcal F\): approximate any \(y\in\mathcal F\) by compact smooth functions, then correct their four continuous linear moments using four fixed compact smooth tests whose moment matrix is invertible (existence follows from linear independence and density). No global form-domain completion is assumed.

## 3. A genuine **one-sided** L² Riesz representer for the full signed Weil pairing

Put \(C_\Gamma=\log(4\pi)+\gamma\), \(w(t)=e^{-t/2}/(1-e^{-2t})\), and \(\kappa_n=\Lambda(n)/\sqrt n\). Extend a smooth \(u\in\mathcal D^0\) by zero outside \(I\). Define a **source-defined** function \(K_2u\) on \(I\) by
\[
\boxed{\begin{aligned}
(K_2u)(x)={}&-C_\Gamma u(x)
-\sum_{2\le n\le54}\kappa_n\bigl(u(x+\log n)+u(x-\log n)\bigr)\\[-1pt]
&-\int_0^\infty\bigl[u(x+t)+u(x-t)-2e^{-t/2}u(x)\bigr]w(t)\,dt.
\end{aligned}}\tag{E.1}
\]
The sum includes only actual prime powers (\(\Lambda(n)=0\) otherwise). The cutoff is **exact** for a window of radius \(a=2\), since \(54<e^4<55\). For every \(y\in\mathcal D\), **including tests that are not themselves pole-null**, the two-pole component of \(q(u,y)\) vanishes because \(u\in\mathcal D^0\). Substitution in Suzuki's exact Weil distribution and Fubini with the subtracted Gamma kernel yields
\[
\boxed{q(u,y)=\langle K_2u,y\rangle_{L^2(I)}.}\tag{E.2}
\]
The pairing is linear in \(u\) and conjugate-linear in \(y\); the equivalent **linear** functional of the complement variable is \(q(y,u)=\langle y,K_2u\rangle\). This is a *one-sided* extension only: it does **not** state that the entire quadratic form \(q(\cdot,\cdot)\) is bounded on \(L^2(I)\times L^2(I)\).

### Proof of convergence and an explicit (deliberately coarse) norm bound

For \(0<t\le1\), elementary convexity gives \(w(t)\le2/t\), \(2(1-e^{-t/2})\le t\), and the global central-difference estimate
\[
\|u(\cdot+t)+u(\cdot-t)-2u\|_{L^2(\mathbb R)}\le t^2\|u''\|_2.
\]
Consequently the small-\(t\) contribution to (E.1) has norm at most \(\|u''\|_2+2\|u\|_2\). For \(t\ge1\), \(w(t)\le2e^{-t/2}\), so the remaining integral has norm at most \(12\|u\|_2\). These bounds justify the Bochner integral, and hence pairing interchange, without assuming any zeta-zero statement.

Furthermore \(|C_\Gamma|<4\); the number of prime powers \(2\le n\le54\) is **24**, and every \(\kappa_{p^k}=\log p/p^{k/2}<1\) since \(\sup_{x>0}\log x/\sqrt{x}=2/e<1\). Thus the finite prime operator term has norm \(<48\|u\|_2\). Combining the pieces yields the entirely explicit fixed-window estimate
\[
\boxed{\|K_2u\|_{L^2(I)}<\|u''\|_{L^2(\mathbb R)}+66\|u\|_{L^2(I)}.}\tag{E.3}
\]
The constants are crude, but (E.3) proves that **every fixed smooth selected source has an L²-continuous mixed Weil functional on the physical complement**. Its derivation uses all of the original Gamma and prime terms and is not a formal assertion that an unbounded global operator is bounded on L².

## 4. Projection removes the **known** finite source and pole blocks exactly

Write \(k_\sigma=K_2u_\sigma\in L^2(I)\). The actual mixed functionals on \(F\) are
\[
 b_\sigma(y)=q(u_\sigma,y)=\langle k_\sigma,y\rangle
=\langle r_\sigma,y\rangle,\qquad
\boxed{r_\sigma=\Pi_\mathcal F k_\sigma=(I-\Pi_E-\Pi_P)k_\sigma.}\tag{E.4}
\]
Because \(F\) is dense in \(\mathcal F\), their **exact operator norm** as (conjugate-linear) functionals on the physical L² completion is
\[
\boxed{\beta_\sigma=\sup_{0\ne y\in F}\frac{|b_\sigma(y)|}{\|y\|_2}=\|r_\sigma\|_2.}\tag{E.5}
\]
The Riesz representer and projected norm are entirely **source-defined**, independently of the Connes semilocal target. At present, \(\beta_\pm\) have **not** been evaluated for a concrete chosen profile; (E.3) bounds them above by the (typically large) \(H^2\)-norm of the source.

For \(a=2\), the pole-direction physical Gram matrix is explicit:
\[
\boxed{G_P=\begin{pmatrix}2\sinh2&4\\4&2\sinh2\end{pmatrix}>0.}\tag{E.6}
\]
Set \(v_\sigma=(\langle k_\sigma,m_+\rangle,\langle k_\sigma,m_-\rangle)^T\). Because \(E\perp P\), \(u_\pm\) are orthonormal, and \(\langle k_\sigma,u_\tau\rangle=\mu_\sigma\delta_{\sigma\tau}\), the exact finite-projection identity is
\[
\boxed{\|r_\sigma\|_2^2=\|k_\sigma\|_2^2-v_\sigma^*G_P^{-1}v_\sigma-\mu_\sigma^2.}\tag{E.7}
\]
This is a **falsifiable mixed-block diagnostic**: once the original source profile is fixed, one can compute the smooth function \(k_\sigma\) using (E.1), integrate its L² norm and its two pole moments, subtract the already known selected eigenvalue, and obtain the *exact physical mixed-coupling norm*, **without constructing an infinite orthogonal basis of the complement**. The right-hand side must be nonnegative. \(\beta_\sigma=0\) if and only if \(k_\sigma\in E\oplus P\), i.e. the cross block vanishes exactly on \(F\).

## 5. The cross-block Schur penalty is **rank at most two**

For \(\mu_\pm>0\) define the positive bounded rank-at-most-two operator on \(\mathcal F\)
\[
\boxed{\mathcal R=\frac{|r_+\rangle\langle r_+|}{\mu_+}+\frac{|r_-\rangle\langle r_-|}{\mu_-},\qquad
\langle\mathcal R y,y\rangle=\frac{|b_+(y)|^2}{\mu_+}+\frac{|b_-(y)|^2}{\mu_-}.}\tag{E.8}
\]
(The rank-one notation is understood with the previously fixed complex inner-product convention.) Its nonzero eigenvalues coincide with those of the **explicit 2×2 cross-coupling Gram matrix**
\[
\boxed{\mathsf B_{\sigma\tau}=\frac{\langle r_\sigma,r_\tau\rangle}{\sqrt{\mu_\sigma\mu_\tau}},\qquad \|\mathcal R\|=\lambda_{\max}(\mathsf B).}\tag{E.9}
\]
The exact source-core Schur criterion derived in T1D-D is therefore
\[
\boxed{q\ge0\ \text{on }\mathcal D^0\quad\Longleftrightarrow\quad
q(y,y)\ge\langle\mathcal R y,y\rangle\quad\forall y\in F,}\tag{E.10}
\]
provided \(\mu_\pm>0\). **Sufficient but not established:** if \(q(y,y)\ge d\|y\|_2^2\) for all \(y\in F\) and \(d\ge\lambda_{\max}(\mathsf B)\), then (E.10) holds. Strict domination gives strict positivity on the pole-null core. The finite-rank penalty is exactly controlled by two source-defined functions; **the genuinely missing theorem is an effective lower bound on the complementary form**. No inequality of the form \(q|_F\ge d\) has been obtained here.

In particular, merely showing \(q(y,y)\ge0\) on \(F\) would still not suffice unless the projected cross coupling vanishes. The standard \(\begin{psmallmatrix}1&2\\2&1\end{psmallmatrix}\) example demonstrates the logical failure.

## 6. Precise status of the cross-canonical interface

- **Original Weil source:** (E.1)–(E.10) are analytic identities or estimates for \(C_c^\infty(-2,2)\) and its **physical L² complement**. The one-sided L² extension is *not* a theorem of full-bilinear L² boundedness.
- **Suzuki:** by [eqs. (8.7)–(8.8)](https://arxiv.org/html/2606.09096v4), the same smooth pairing can be represented after differentiation with the numerator \(G_a\) and denominator \(K_a\); transferring the complement requires transporting **both** and specifying the completed domain. Our \(K_2\) is a newly defined *Weil-source representer* and **must not be confused with Suzuki's \(K_a=(-\Delta_N)^{-1}\)**.
- **Connes:** the smooth pole-null semilocal trace identity, under T1A's **separate** complex domain argument, identifies the values \(q(u_\sigma,y)\) for smooth \(y\in F\). It does **not** imply that the trace operator extends to all of \(\mathcal F\), nor identify \(\mathcal F\) with Connes's high-frequency Sonin/prolate-conditioned complement. The latter relies on still-open semilocal **H2′**, and the original article's **R1/R2** inputs remain conditional.
- **Pole directions:** the two nonsmooth windowed vectors \(m_\pm\) are used **only** for physical projection, not substituted as smooth Weil test directions. Positivity on \(\mathcal D^0\) still does not cover the two additional independent pole directions of the full RH criterion.

## 7. Disposition / reproducibility / next falsifier

**DERIVED / ANALYTIC PASS:** an explicit absolutely convergent original signed Weil **mixed Riesz representer** on the smooth pole-null core at \(a=2\); uniform (for each fixed u) bound \(\|K_2u\|_2<\|u''\|_2+66\|u\|_2\); exact physical pole/projector subtraction; an exact rank-two cross-block operator whose norm is determined by a computable 2×2 Gram matrix. The auditable finite arithmetic constants, prime-power census, finite Gram/projection identities and rank-two spectral algebra are checked by `src/audits/canonical_return_t1d_e_mixed_bound.py`. The Bochner-integral/Fubini and smooth-core-density arguments are **human analytical proofs**, not consequences of symbolic assertions or Lean elaboration.

**OPEN:** actually evaluate or interval-enclose \(\beta_\pm\) for a fixed explicit smooth profile; estimate \(q\) from below on the **infinite-dimensional** \(F\); establish any usable bound \(d\ge\lambda_{\max}(\mathsf B)\); transfer \(\mathcal F\) and its form domain to Connes's actual conditioned subspace; complete Suzuki/Connes spectral equivalence; add the two pole directions to the full test class; prove RH.

**Next — CR-INT-T1D-F:** for the already frozen T1C smooth bump, construct independently certified numerical intervals for the **projected** cross-coupling Gram \(\mathsf B\), and compare them against the *actual available* complementary lower bound. If the interval is grossly larger than all candidate lower bounds, the method fails as currently scaled; if it is small, the analytic remaining task is to prove complementary coercivity under Connes's precise semilocal hypotheses. Do not substitute the crude \(H^2\) estimate (E.3) for a successful Schur comparison, and do not infer RH from any selected 2D identity.