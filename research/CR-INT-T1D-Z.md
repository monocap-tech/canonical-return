# CR-INT-T1D-Z — A quantitative second Schur penalty and a positive conditioned infinite remainder

**Canonical Return · November 23, 2026 · independent original-signed-Weil smooth-source research.**  
**Disposition:** NEW **ANALYTICAL / CONDITIONAL-ON-PREVIOUS-SOURCE-LEMMA** estimate with an executable exact-fraction audit. This note constructs the actual rank-at-most-four mixed functional operator after removing the four T1D-Y central/satellite directions, bounds that operator explicitly **without T1D-I**, and proves two positive **infinite-dimensional high-frequency conditioned** second-Schur subspaces. The second estimate uses an extremely large cutoff; neither statement proves positivity on the full remaining complement, a spectral upper bound below one, Suzuki–Connes completed equivalence, Lean, Phase IV, or RH.

## 1. Frozen source, metric, and exact second gate

Keep Suzuki's [original complete signed Weil functional](https://arxiv.org/html/2606.09096v4) on `I=(-2,2)`, with both zeta pole terms, all 24 genuine von Mangoldt prime powers `n<=54`, and the complete Gamma distribution. The pole terms vanish **only on** the unchanged smooth double-Mellin-pole-null core

\[
\mathcal D^0=\left\{f\in C_c^\infty(I;\mathbb C): \int_I f(x)e^{\pm x/2}\,dx=0\right\}.
\]

The physical metric is `H=L²(I,dx)`, linear in its first variable. Freeze the original physically orthonormal selected T1C pair `E=span{u_+,u_-}`, the pole carrier span `P=span{e^{x/2},e^{-x/2}}`, the physical complement `H_F=(E⊕P)^perp`, and the genuine smooth source `F_s=D^0 ∩ H_F`. The original complete signed Weil form is `q(f,g)`; its selected-space Schur correction of [T1D-E](CR-INT-T1D-E.md) is

\[
\mathcal R=\sum_{\sigma=\pm}\frac{|r_\sigma\rangle\langle r_\sigma|}{\mu_\sigma},\qquad
r_\sigma=\Pi_{H_F}\mathcal K_2u_\sigma,\quad \mu_\sigma=q[u_\sigma]>7,
\]

and the Schur remainder is `s(f,g)=q(f,g)-<R f,g>`. Here `K_2` is the original **Weil mixed representer**, **not** Suzuki's different inverse Neumann Laplacian. The `s` form is symmetric with first-argument linearity. It is still not known to be positive on all of `F_s`.

Freeze the **four** physical-orthonormal, real, parity-pure, smooth central/satellite vectors from [T1D-Y](CR-INT-T1D-Y.md):
\[
\begin{aligned}
\ell_{c,\sigma}&=\widehat y_\sigma,&
\ell_{s,\sigma}&=\widehat z_\sigma,\\
 y_\sigma&=(\partial_x^2-1/4)(b_\sigma*\rho_{10^{-7}}),&
 z_\sigma&=(T_{\log3}y_\sigma+T_{-\log3}y_\sigma)/\sqrt2,\\
 b_+(x)&=(1-16x^2)^3\mathbf1_{|x|<1/4},&
 b_-(x)&=x(1-16x^2)^3\mathbf1_{|x|<1/4}.
\end{aligned}
\]

The mollifier is the same *fixed* even smooth probability mollifier as in T1D-W/X/Y, not a raw polynomial in place of the genuine smooth source. Define
\[
L=\operatorname{span}\{\ell_{c,+},\ell_{s,+},\ell_{c,-},\ell_{s,-}\}\subset F_s,
\quad H_Z=H_F\cap L^\perp,\quad Z_s=F_s\cap L^\perp.
\]
The original prime self-pairing vanishes on the central support, but **not** on the central-to-satellite or satellite-to-satellite blocks. The complete original signed entries, including all prime powers active at those separations and the Gamma term, are those independently audited in T1D-W/X/Y.

By [T1D-Y §5](CR-INT-T1D-Y.md), the authentic `s`-Gram on `L` decomposes by reflection into two positive real `2×2` matrices `S_+` and `S_-`, with strictly certified
\[
\det S_+>1/20,\qquad\det S_->1/5.
\tag{Z.1}
\]
The unproved **whole remainder** is the exact second Schur condition
\[
\boxed{\quad s(w,w)-b(w)^*S_L^{-1}b(w)\ge0\quad\text{for every }w\in Z_s,\quad}
\tag{Z.2}
\]
where `b_i(w)=s(ell_i,w)` is a conjugate-linear functional of `w` in the chosen inner-product convention. Equation (Z.2) is an equivalence to positivity of `s` on all `F_s`, because `L` is smooth, physical-orthonormal, and `S_L>0`; **(Z.2) itself is not established**.

## 2. Original-source Riesz functions and a genuine rank-at-most-four penalty

[T1D-E](CR-INT-T1D-E.md) establishes that every fixed smooth double-pole-null source `f` has a source-preserving `L²` representer `K_2 f`, defined with the complete Gamma, negative prime, and pole-null-original form, such that
\[
q(f,w)=\langle K_2f,w\rangle_{L²(I)},\qquad w\in F_s.
\]
Define four **real and source-defined** mixed residuals by the **physical** projection
\[
\boxed{g_i=\Pi_{H_Z}\bigl(\mathcal K_2\ell_i-\mathcal R\ell_i\bigr)\in H_Z.}\tag{Z.3}
\]
Then `b_i(w)=<g_i,w>` for all `w∈Z_s`. The additional Schur penalty on the remaining physical Hilbert space is exactly
\[
\boxed{\mathcal T=G S_L^{-1}G^*,\qquad
G:\mathbb C^4\to H_Z,\quad G e_i=g_i,\quad
\langle\mathcal T w,w\rangle=b(w)^*S_L^{-1}b(w).}\tag{Z.4}
\]
The transpose/conjugate placement in (Z.4) is harmless for the real symmetric, parity-pure frozen basis; with arbitrary complex coordinates the operator is defined by the displayed nonnegative quadratic form. It is bounded, positive and rank at most **four** (at most **two per parity**). This is *not* the original rank-two `R`: it is a second independent penalty on `H_Z`, created when the certified central/satellite directions are eliminated. The exact objective becomes
\[
q[w]-\langle\mathcal Rw,w\rangle-\langle\mathcal Tw,w\rangle\ge0
\qquad(w\in Z_s).\tag{Z.5}
\]

## 3. Strict inverse-matrix bounds *without* T1D-I

We need a lower eigenvalue bound for `S_L`. Use source-level T1D-W's **whole central two-mode Gamma Rayleigh upper** `q[ell_{c,+}]<31/20` and `q[ell_{c,-}]<9/5`. From the complete original T1D-X satellite-wing census and Gamma far-cross estimate,
\[
q[\ell_{s,+}]<q[\ell_{c,+}]+29/100+1/1000,\qquad
q[\ell_{s,-}]<q[\ell_{c,-}]-9/100+1/1000.
\]
Since `R` is nonnegative, `s[ell_i] <= q[ell_i]`. Therefore
\[
\boxed{\operatorname{tr}S_+<2(31/20)+29/100+1/1000=3.391<7/2,}\tag{Z.6}
\]
\[
\boxed{\operatorname{tr}S_-<2(9/5)-9/100+1/1000=3.511<18/5.}\tag{Z.7}
\]
By positivity of both original `2×2` Schur Grams, the exact identities `lambda_min=det/lambda_max` and `lambda_max<trace`, and (Z.1),
\[
\boxed{\lambda_{\min}(S_+)>1/70,\qquad
\lambda_{\min}(S_-)>1/18.}\tag{Z.8}
\]
All the constants come from the **complete** earlier source matrices, not a prime-only approximation; no bound from the contested T1D-I residual-norm analysis is imported.

## 4. Explicit norms of the *true* off-block Weil representers

[T1D-F](CR-INT-T1D-F.md) gives on the original whole-line Gamma multiplier
\[
|m_\Gamma(\xi)-\log(1+|\xi|)|<40.
\]
For every smooth compact `f`, Plancherel and `log(1+t)<=t` imply
\[
\boxed{\|\Gamma f\|_2\le40\|f\|_2+\|f'\|_2.}\tag{Z.9}
\]
The source representer is `K_2 f=(Gamma-A)f` on the pole-null source with the original compressed prime operator `A`, so the triangle inequality bounds it without changing the signed form. [T1D-W/X](CR-INT-T1D-X.md) give the physical unit even/odd normalized derivative bounds `<38`, `<46`, preserved for the reflected satellite pair. [T1D-U](CR-INT-T1D-U.md) gives `||A P_J||<7/2` on the central support, while [T1D-Q](CR-INT-T1D-Q.md) gives `||A||<203/25` everywhere in the physical window.

To control the `R` term without T1D-I, use only [T1D-F's](CR-INT-T1D-F.md) **one-sided** fixed-width `||K_2u_sigma||<112`, so `||r_sigma||<112`, with `mu_sigma>7`. Because opposite parities decouple,
\[
\|\mathcal R\ell\|_2
\le16\,|q(u_\sigma,\ell)|.
\tag{Z.10}
\]
Use the independently certified local selected-to-central cross bounds of T1D-W, `<2/15` even and `<3/25` odd, and direct selected-to-satellite bounds of T1D-Y, `<27/100` even and `<23/100` odd. Contracting the physical projection in (Z.3), we obtain the *strict numerical inequalities*
\[
\begin{array}{c|cc}
 &\|g_{c,\sigma}\|_2&\|g_{s,\sigma}\|_2\\\hline
\sigma=+&<40+38+7/2+16(2/15)<\boxed{84}&<40+38+203/25+16(27/100)<\boxed{91}\\
\sigma=-&<40+46+7/2+16(3/25)<\boxed{92}&<40+46+203/25+16(23/100)<\boxed{98}.
\end{array}\tag{Z.11}
\]
The estimates are coarse and involve **actual Gamma/prime-mixed representers**, not substitute finite Gram norms. They hold as `L²` operator bounds on *all* of `H_Z`, not just their projections onto finitely many probes.

By parity, the two `S` blocks and corresponding `g`-images decouple. The elementary bound `||G_sigma||²<=sum_i||g_i||²`, together with (Z.8), proves
\[
\begin{aligned}
\|\mathcal T_+\|&<70(84^2+91^2)=\boxed{1,073,590},\\
\|\mathcal T_-\|&<18(92^2+98^2)=\boxed{325,224},\\
\boxed{\|\mathcal T\|&<1,073,590.}
\end{aligned}\tag{Z.12}
\]
This is a **nonzero, finite, explicit** norm control for the *actual second Schur operator*, **independent of T1D-I**. It is spectacularly loose. It does not imply coercivity of `q-R` on `Z_s`, or positivity of (Z.5). No actual numerical `||g_i||` or eigenvalues of `T` have been computed.

## 5. Two positive infinite-dimensional *conditioned* second-Schur remainders

A second-stage high-frequency result can now be stated without falsely extending positivity to all of `Z_s`. Use the already established [T1D-M/Q](CR-INT-T1D-Q.md) original complete Gamma lower bounds and prolate concentration theorem. For each integer `k>=6`, let `B=2^k`, `theta=1/128`. On the genuine smooth source subspace orthogonal (physically) to all eigenfunctions of the windowed prolate operator `P_I F^{-1}1_{[-B,B]} F P_I` with concentration eigenvalue `>theta`, the Fourier mass in `[-B,B]` is `<=theta` of the full mass. The original prime contribution remains bounded below by `-203/25` from T1D-Q. The **full source Weil form** therefore satisfies
\[
q[w]>
\left\{\frac{127}{128}\left(k\log2-\frac{49}{10}\right)-\frac6{128}-\frac{203}{25}\right\}\|w\|_2^2.
\tag{Z.13}
\]
The prolate eigenvector constraints are finitely many for every finite `B`, and intersecting them with the **smooth eight-constraint physical source** `Z_s` still leaves an infinite-dimensional space. We only impose **orthogonality** to prolate modes; those modes are not treated as original smooth Weil tests.

**Route A — four extra exact mixed-functional constraints, moderate but still enormous cutoff.** Impose `b_i(w)=<g_i,w>=0` for all four `i`, in addition to high-frequency conditioning. This imposes at most four additional linear constraints and annihilates `T[w]` *exactly*. Import only T1D-F's unconditional-to-this-note, coarse `||R||<1785`; from the rigorous `log2>693/1000` and (Z.13) at `k=2700`, the source second Schur form satisfies
\[
\boxed{q[w]-\mathcal R[w]-\mathcal T[w]>
\frac{187051}{3200}\|w\|_2^2>58\|w\|_2^2.}\tag{Z.14}
\]
This still leaves an infinite-dimensional smooth subspace of `Z_s`. No reliance on T1D-I is needed.

**Route B — without annihilating the mixed functionals.** On **all** smooth `w∈Z_s` meeting the prolate high-frequency condition, use (Z.12), T1D-F `||R||<1785`, and (Z.13). The deliberately enormous `k=1,570,000` then yields
\[
\boxed{q[w]-\mathcal R[w]-\mathcal T[w]>
\frac{26380117}{6400}\|w\|_2^2>4000\|w\|_2^2.}\tag{Z.15}
\]
Both rational margins are recalculated by the companion standard-library audit. **The cutoff `B=2^(1,570,000)` is mathematically finite but computationally unusable.** Neither (Z.14) nor (Z.15) covers the omitted low-frequency sector and neither settles the sign of the full complement. These results show only that the *new* finite-rank penalty can be dominated, or annihilated exactly, on explicitly conditioned infinite-dimensional portions of the original smooth source.

## 6. Claim level, checks and the remaining obligation

| Obligation | Status |
| --- | --- |
| Full-source mixed `L²` representers `g_i` and rank-at-most-four `T` | **DERIVED / ANALYTIC** using T1D-E |
| Positivity of finite `S_±`, strict inverse bounds | **DERIVED using Y/W/X**, rational `AUDIT-PASS` |
| Strict `||g_i||` and global `||T||<1,073,590` | **DERIVED using F/Q/U/W/Y**, rational `AUDIT-PASS`; deliberately loose |
| Positive second Schur on two infinite-dimensional high-pass source spaces | **DERIVED with stated previous Gamma/prolate/F inputs**, rational `AUDIT-PASS` |
| Actual sharp `||g_i||`, low-high `eta_N`, and omitted low-frequency matrix | **OPEN** |
| Positive second Schur on all `Z_s`, whole Weil positivity, Sonin transport, Phase IV, RH | **OPEN / NO TRANSFER** |

**Audit:** `python3 src/audits/canonical_return_t1d_z_second_schur_highpass.py` (standard Python library). The audit **rebuilds** the earlier Y determinant lower arithmetic, checks the W/X trace budget, strict representer norms from F/Q/U/W/Y, parity operator-norm arithmetic, an independent rational atanh lower bound for `log2`, and both high-pass margins. It does **not** compile analytic functional analysis, independently verify earlier Gamma/prolate lemmas, construct actual prolate eigenvectors, or calculate an upper spectral value below 1. All original source and metric distinctions are preserved.

**Next — CR-INT-T1D-AA:** obtain an **operator-adapted bound** for the new `g_i` residuals (their actual norms, projections, and low-frequency localization) and seek a much smaller, practical cutoff. In parallel, directly test the **excluded low-frequency `Z` block**, including the second rank-four penalty, instead of mistaking finite-codimension positivity for the missing whole-complement theorem. A completed Suzuki–Connes semilocal form/domain dictionary remains separate and unproved.
