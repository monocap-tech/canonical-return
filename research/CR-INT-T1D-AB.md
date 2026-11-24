# CR-INT-T1D-AB — A near-Shannon **nonasymptotic** upper bound for the extremely small prolate threshold, and a conditional finite Morse-index certificate

**Canonical Return · November 23, 2026 (America/Los_Angeles) · independent analytical source investigation.**  
**Disposition:** new **DERIVED / ANALYTIC** finite-rank approximation theorem for the *ordinary logarithmic-coordinate* time–band concentration operator, with exact-rational numerical certification. A consequent **explicit finite upper bound on the negative index of the original fixed-window smooth two-pole-null Weil form** imports the previous T1D-Y/AA analytic Schur and high-pass estimates. No computation of individual prolate eigenvalues, no low-mode Weil positivity, no new Connes/Sonin intertwiner, no Lean proof and no RH. This follows [T1D-AA](CR-INT-T1D-AA.md).

## 1. Frozen source and a threshold so small that the old upper count fails

Retain Suzuki's complete **original signed Weil form** ([*Weil's quadratic form via the screw function*, arXiv:2606.09096v4 §1.1, §2.4](https://arxiv.org/html/2606.09096v4)), the native physical metric `H=L²((-2,2),dx)`, and genuine smooth pole-null source

\[
\mathcal D^0=\{f\in C_c^\infty((-2,2)):\int f(x)e^{x/2}\,dx=\int f(x)e^{-x/2}\,dx=0\}.
\]

The T1D-Y selected-plus-central/satellite positive blocks have total physical dimension six (two selected modes and four new central/satellite modes). Together with the two physical pole carriers, they leave a closed physical remainder `\mathscr Z` of **codimension eight**, and its dense smooth core `Z_s=\mathcal D^0\cap L^\perp`. The second Schur remainder on `Z_s` is

\[
\mathfrak h[w]=q[w]-\mathcal R[w]-\mathcal T[w].
\]

T1D-AA proves, with its **explicit separate analytic source hypotheses**, `\mathfrak h[w]>(1233/1750)\|w\|²` for all nonzero `w\in Z_s` satisfying the band-concentration requirement at

\[
\boxed{B=2^{21},\qquad\theta=2^{-38}.}\tag{AB.1}
\]

Here `\mathsf T_B=P_I\mathcal F^{-1}\mathbf1_{[-B,B]}\mathcal F P_I` is the *ordinary* physical interval prolate operator. The excluded rank is

\[
 m_{B,\theta}=\dim\operatorname{Ran}\mathbf1_{(\theta,1]}(\mathsf T_B).
\]

The trace-defect theorem derived in [T1D-N](CR-INT-T1D-N.md) gives `m_{B,\theta}\le4B/\pi+D_B/\theta`. At the current tiny `\theta`, that upper estimate is **approximately 9.72×10¹¹**—far too loose. This is a weakness of that counting inequality, *not* evidence that almost one trillion prolate modes are genuinely excluded. The lower estimate from the same theorem remains useful.

## 2. Refined rational lower bound from the exact trace and defect

T1D-N proves analytically (not by numerically enumerating eigenvalues)

\[
\operatorname{Tr}\mathsf T_B=\frac{4B}{\pi},\quad
0\le D_B:=\operatorname{Tr}(\mathsf T_B-\mathsf T_B^2)
\le\frac2{\pi^2}\left(\log(4B)+\frac32\right),
\]

and, for `0<\theta<1`,

\[
\frac{4B}{\pi}-\frac{D_B}{1-\theta}\le m_{B,\theta}
\le\frac{4B}{\pi}+\frac{D_B}{\theta}.\tag{AB.2}
\]

At `B=2^21`, `\log(4B)=23\log2`. Exact rational Machin-series enclosures of `\pi` and atanh-series enclosures of `\log2` prove

\[
\frac{4B}{\pi}-\frac{2(23\log2+3/2)}{\pi^2(1-2^{-38})}>2{,}670{,}173.
\]

As the rank is an integer,

\[\boxed{m_{2^{21},2^{-38}}\ge2{,}670{,}174.}\tag{AB.3}\]

The trace itself is approximately `2,670,176.857720`; the rationally enclosed defect upper bound is approximately `3.534567`. Both decimals are *diagnostics*; the proof compares integer numerators of exact fractions.

## 3. A compact Fourier window gives a **far stronger upper bound**

**The analytic mechanism is a finite-rank approximation to the true time–band restriction, not a prolate eigenvalue computation.** Set

\[
 (S_B f)(x)=\frac1{\sqrt{2\pi}}\int_{-B}^B f(\xi)e^{ix\xi}\,d\xi,
 \quad S_B:L^2([-B,B])\to L^2((-2,2)).
\]

Then `\mathsf T_B=S_BS_B^*`; its eigenvalues are the squares of the singular values of `S_B`.

Fix the *auxiliary*, even window parameters

\[
\delta=\frac1{120},\quad a=2+\delta=\frac{241}{120},\quad
k=14,\quad C=\frac{2k}{\delta}=3360.
\]

Let `u` be the normalized uniform density on `[-\delta/(2k),\delta/(2k)]` and set

\[
\boxed{\chi=\mathbf1_{[-2-\delta/2,\,2+\delta/2]}*u^{*k}.}\tag{AB.4}
\]

It is even, `0\le\chi\le1`, supported on `[-a,a]`, is at least `C^{k-2}`, and equals **one identically** on `(-2,2)`. This auxiliary window is not a replacement of the original smooth Weil test core.

Its explicit ordinary Fourier transform satisfies

\[
\begin{aligned}
\widehat\chi_{\rm ordinary}(t)&=\frac{2\sin((2+\delta/2)t)}{t}
\left(\frac{\sin(\delta t/(2k))}{\delta t/(2k)}\right)^k,\\
|\widehat\chi_{\rm ordinary}(t)|&\le \frac{2C^k}{|t|^{k+1}}\qquad (t\ne0).
\end{aligned}\tag{AB.5}
\]

Let `P_N` be the *orthogonal Fourier series projection* onto `e^{i\pi nx/a}` for integers `|n|\le N` in `L²((-a,a))`. Define the finite-rank approximation

\[
\boxed{S_{B,N}f=\left[P_N(\chi S_Bf)\right]\big|_{(-2,2)}.}\tag{AB.6}
\]

Clearly `\operatorname{rank}S_{B,N}\le2N+1`. For `\Omega=\pi N/a-B>0`, Cauchy–Schwarz on each Fourier coefficient, Parseval and the monotone integral test give the **uniform operator estimate**

\[
\boxed{
\|S_B-S_{B,N}\|^2
\le\frac{4B}{\pi^2(2k+1)\Omega}
\left(\frac C\Omega\right)^{2k}.}\tag{AB.7}
\]

**Derivation of the constants.** For the orthonormal Fourier basis `e_n(x)=e^{i\pi nx/a}/\sqrt{2a}`, the coefficient of `\chi S_B f` is

\[
 c_n=\frac1{\sqrt{4\pi a}}\int_{-B}^B
 f(\xi)\widehat\chi_{\rm ordinary}(\pi n/a-\xi)\,d\xi.
\]

For `|n|>N`, `|\pi n/a-\xi|\ge\pi|n|/a-B`; applying (AB.5), squaring, summing both signs and integrating the decreasing tail yields

\[
 \sum_{|n|>N}|c_n|^2
 \le\frac{4B C^{2k}}{\pi a}\sum_{n>N}(\pi n/a-B)^{-2k-2}\|f\|^2
 \le\frac{4B C^{2k}}{\pi^2(2k+1)\Omega^{2k+1}}\|f\|^2.
\]

Since `\chi=1` on the target interval, restriction to that interval cannot increase the approximation error. This proves (AB.7), including its physical `L²` normalization.

Choose the **explicit integer** `N=1,346,600`. The exact-rational script proves `\pi>314159/100000`, hence

\[
 \Omega=\frac{\pi N}{a}-2^{21}>9300.
\]

Using only `\pi>3` for the remaining error gives the strict fully rational comparison

\[
\boxed{
\|S_B-S_{B,N}\|^2
<\frac{4\cdot2^{21}}{9\cdot29\cdot9300}
\left(\frac{3360}{9300}\right)^{28}
<2^{-38}=\theta.}\tag{AB.8}
\]

If a compact operator has a rank-`r` approximant within `<\sqrt\theta` in operator norm, its `(r+1)`-st singular value is `<\sqrt\theta` by min–max. Therefore **at most `r` singular values exceed `\sqrt\theta`**, giving

\[
\boxed{m_{2^{21},2^{-38}}\le 2N+1=2{,}693{,}201.}\tag{AB.9}
\]

Combining (AB.3) and (AB.9) yields the **nonasymptotic rank interval**

\[
\boxed{2{,}670{,}174\ \le\ m_{2^{21},2^{-38}}\ \le\ 2{,}693{,}201.}\tag{AB.10}
\]

This is a width-`23,027` bracket around a count near the Shannon number `4B/\pi`. The improvement is **on the rank upper bound**, not on the actual low-mode Weil form. It is source-independent time–band geometry; the Wiener/B-spline windowing construction is elementary harmonic analysis, not a claim of a new canonical prolate theory. The classical Landau–Widom asymptotics explain why such concentration near the Shannon count is natural, but the present inequality is a separate explicit finite-parameter proof and does not invoke their asymptotic error term.

## 4. Parity and exact smooth-source constraint burden

Because `\chi` is even, the finite Fourier approximation commutes with reflection; its target even modes are `1,\cos(\pi x/a),...,\cos(N\pi x/a)` (dimension `N+1`) and its odd modes are `\sin(\pi x/a),...,\sin(N\pi x/a)` (dimension `N`). Thus

\[
\boxed{m_{B,\theta,+}\le1{,}346{,}601,
\qquad m_{B,\theta,-}\le1{,}346{,}600.}\tag{AB.11}
\]

T1D-N's reflected trace estimate `|\operatorname{Tr}(J\mathsf T_B)|<4/3`, with `Jf(x)=f(-x)`, implies that **each** parity trace exceeds `2B/\pi-2/3`. Each parity trace defect is nonnegative and no larger than `D_B`. The same exact rational substitutions give

\[
\boxed{1{,}335{,}085\le m_{B,\theta,+},\qquad
1{,}335{,}085\le m_{B,\theta,-}.}\tag{AB.12}
\]

Let `Z_s` be T1D-Y's smooth physical complement after the original two pole carriers, two selected modes, and four central/satellite modes: its closure `\mathscr Z` has **codimension eight**, split as four in each parity. The prolate moment map on `Z_s` has rank between `m_{B,\theta}-8` and `m_{B,\theta}` by finite-dimensional dimension counting and smooth-core density. Consequently it requires **at least**

\[
\boxed{2{,}670{,}166\text{ further independent smooth constraints on }Z_s,}\tag{AB.13}
\]

and at least `1,335,081` in **each** parity sector. Its rank cannot exceed `2,693,201`. Merely adding finitely many certified central/satellite vectors does not eliminate this million-dimensional burden **if one insists on the same uniform low-band leakage cap**. This remains a limitation of the *conditioning method*, not of source-adapted approaches to the full Weil form.

## 5. A precise **conditional negative-index corollary**, not global positivity

Write `n_-(q;\mathcal D^0)` for the supremum of dimensions of subspaces of the original smooth pole-null source on which `q` is strictly negative definite. T1D-Y certifies the positive selected and four-mode complementary Schur matrices; T1D-AA establishes positivity of the **second Schur remainder** on the smooth prolate-conditioned subspace of `Z_s`.

By completing both finite Schur squares, negative directions of the original form correspond bijectively (as regards negative *index*) to negative directions of the second Schur remainder on `Z_s`. Any negative-definite subspace of `Z_s` has trivial intersection with the conditioned positive subspace. Its dimension therefore cannot exceed that subspace's codimension, itself no larger than `m_{B,\theta}`. **Conditional on the previous source-level analytical Weil/Gamma/Schur results being correct**, the new rank upper certificate hence yields

\[
\boxed{n_-(q;\mathcal D^0)\le2{,}693{,}201.}\tag{AB.14}
\]

The parity-specific negative indices are at most `1,346,601` (even) and `1,346,600` (odd), under the same hypotheses. This is a **quantified fixed-window finite Morse-index bound**, not positive semidefiniteness; a single genuine negative direction would still be significant. It is not an additional theorem of Suzuki, who already supplies broader compactness/discreteness mechanisms. This explicit numerical index bound is a consequence of **our unreviewed** source-level chain, not independently reviewed mathematics.

## 6. Audit, falsifiers, and stopping jurisdiction

**Exact reproduction:** `python3 src/audits/canonical_return_t1d_ab_prolate_rank_window.py` (standard Python library). Its integer/fraction computations certify: Machin/atanh rational enclosures; (AB.3); the auxiliary window's `a,k,C,N`; `\Omega>9300`; the squared approximation upper bound `<\theta`; the rank and parity limits; and the physical codimension subtraction. The analytic derivation of (AB.7), the singular-value approximation theorem, the original prolate trace/defect theorem, reflection trace, smooth-core density and the negative-index/inertia transfer are **written mathematical arguments, not symbolic-execution outputs**. They remain open to independent review.

**Not established:** the actual prolate eigenvalue sequence or exactly where its smallest excluded eigenvalue lies; a positive Weil matrix on the omitted `~2.7 million` low-frequency directions; the top eigenvalue of the full compact Schur operator; any bound on the *unconditioned* second Schur remainder; all physical pole carriers of the original Weil criterion beyond the two-pole-null source; Connes/Sonin completed-domain and native-metric transport; Phase IV positivity; RH.

**Next — CR-INT-T1D-AC:** the prolate rank is now bounded near its natural Shannon count even at `\theta=2^{-38}`. Do not spend the next increment merely tightening this numerical rank bracket. Instead exploit the explicit finite Fourier approximation `P_N\chi S_B` to build a **source-preserving low-mode transfer** with controlled Gamma energy, arithmetic translation couplings and a smooth moment-lift/conditioning number. Its finite-rank *existence* alone gives no useful Weil positivity on the low block. Alternatively find a robust negative direction by an interval-certified full original signed Weil test. **No proof transfer or RH declaration.**