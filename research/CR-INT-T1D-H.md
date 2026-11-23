# CR-INT-T1D-H — Certified logarithmic Fourier variance and a finite-width mixed-coupling witness

**Canonical Return · independent mathematical research note · November 22, 2026 (PST) / November 23 UTC.**  
**Status:** **ANALYTIC BOUNDS + EXACT-RATIONAL AUDIT PASS.** New results: a rigorous interval for the **actual** T1C smooth-profile Fourier-log variance, and explicit nonzero lower bounds for **each** true projected mixed residual norm at its actual width `ell=1/10000`. This does **NOT** fully enclose both finite-width residual norms, prove complementary coercivity, identify the completed Suzuki–Connes operators, supply Lean acceptance, or prove RH. This note documents **partial completion** of the T1D-H target and pinpoints the remaining high-precision quadrature obligation. See [T1D-G](CR-INT-T1D-G.md), [T1D-E](CR-INT-T1D-E.md), [T1C](../archive/space/10-t1c.md), and [Issue #1](https://github.com/monocap-tech/canonical-return/issues/1).

## 1. Source contract and frozen test

The full original signed Weil formula, including poles, negative von Mangoldt prime terms, and the complete Gamma distribution, is Suzuki's [arXiv:2606.09096v4 §1.1 and §2.4, equations (2.7)–(2.8)](https://arxiv.org/html/2606.09096v4). Its Gamma Fourier multiplier has asymptotics `log|xi| - log(2*pi)`. On the smooth pole-null two-dimensional source, [T1D-E](CR-INT-T1D-E.md) constructs the original Weil-source representer `K_2` and [T1D-F](CR-INT-T1D-F.md) diagonalizes the even/odd projected mixed Gram by parity. These canonical identities and previous analysis are assumed as inputs, **not independently verified operator identities by this script**.

Fix the original T1C **specific nonnegative, even, smooth** profile

\[
 b_h(t)=(1-|t|/h)_+,\quad h=49/100,\quad\rho_\varepsilon(t)=Z_\varepsilon^{-1}\exp\!\left[-\frac{1}{1-(t/\varepsilon)^2}\right]\mathbf1_{|t|<\varepsilon},\quad \varepsilon=1/100,\quad\psi_0=b_h*\rho_\varepsilon.
\]

Set `v=psi_0/||psi_0||_2`; at actual width `ell=1/10000` the physical unit cell is `U_{ell,x}v(t)=ell^(-1/2) v((t-x)/ell)`. For `d=log3`, `s=4/sqrt3`, and centers `x_j=(j-3/2)d`, the physical unit selected directions are

\[
 u_{\sigma}=N_\sigma^{-1/2}\sum_{j=0}^3 c_{\sigma,j}U_{\ell,x_j}v,\quad
 c_+=(1,1-s,1-s,1),\quad c_-=(1,-1-s,1+s,-1),\quad N_\pm=(44\mp16\sqrt3)/3.
\]

Both Mellin pole moments vanish. Set `P=span{e^(x/2),e^(-x/2)} in L2(-2,2)`, `F=(P+span{u_+,u_-})^perp`, `r_sigma=Pi_F K_2 u_sigma`, `beta_sigma=||r_sigma||_2`. The actual target Schur penalty is `||R|| = max(beta_+^2/mu_+, beta_-^2/mu_-)`, by the parity result of T1D-F and assuming `mu_±>0` (independently checked in T1C).

## 2. Exact *triangular* Fourier-log variance as a computable reference

Under the **unitary Fourier transform** and change of variable `x=eta*h/2`, the normalized spectral distribution of `b_h` has the universal density

\[
 p_\triangle(x)=\frac{3}{2\pi}\left(\frac{\sin x}{x}\right)^4,\qquad \int_{\mathbb R}p_\triangle(x)dx=1.
\]

The Mellin transform `M(z)=E_triangle(|x|^z)` near `z=0` and `z=2` is exactly

\[
 \boxed{M(z)=\frac{3\cdot2^{2-z}(2^{1-z}-1)}{2\,\Gamma(4-z)\cos(\pi z/2)}}\qquad (-1<\Re z<3),
\]

interpreting the removable sign changes across `z=1` by continuity. This follows by expanding `sin^4 x=(3-4cos 2x+cos4x)/8`, evaluating its Mellin transform by twice integrating by parts / analytic continuation, and then differentiating where convergence is absolute. Exact differentiation yields

\[
\begin{aligned}
E_\triangle(\log|x|)&=11/6-\gamma-3\log2=:m_1,\\
\operatorname{Var}_\triangle(\log|x|)&=\pi^2/12+49/36-2(\log2)^2=:V_\triangle,\\
E_\triangle(x^2)&=3/4,\\
E_\triangle(x^2\log^2|x|)&=\frac34\left[1+\frac{\pi^2}{12}-2(\log2)^2+(1-\gamma)^2\right]=:T_2.
\end{aligned}
\]

A separate symbolic derivative check reproduces these formulas. The **rational outward interval** from elementary series for `pi`, `log2`, Euler's constant is

\[
 \boxed{1.2226721166\le V_\triangle\le1.2226721167.}
\]

This is not yet `V_v`: the original carrier was explicitly **mollified**, so replacing it by the triangle without an error bound would be invalid.

## 3. Rigorous perturbation from the triangle to the **actual C-infinity T1C bump**

The normalized mollifier is an even probability density supported in `[-epsilon,epsilon]`. Thus its real Fourier transform satisfies

\[
 0\le1-|\widehat\rho_\varepsilon(\eta)|^2\le\varepsilon^2\eta^2
\]

(the inequality follows from `1-cos(eta t)<=eta²t²/2`, `|rhohat|<=1`). Let `w(x)=|rhohat(2x/h)|²`. The true normalized smooth spectral probability law is exactly `p_triangle(x) w(x)/E_triangle[w]`. Crucially, the ratio `epsilon/h=1/49` is **exact**; no quadrature approximation to `rho` is used.

With `D0=3/2401`, `D1=(2/2401)(T2+3/4)`, and `D2=(4/2401)T2`, it follows from `2|log|x||<=1+log²|x|` that

\[
 0\le E_\triangle(1-w)\le D_0,\quad
 |E_\triangle[(1-w)\log|x|]|\le D_1,\quad
 0\le E_\triangle[(1-w)\log^2|x|]\le D_2.
\]

Propagate these **correlated parameter-free inequalities** through the normalized first and second logarithmic moments. Since adding the constant `log(2/h)` does not affect a variance, this directly bounds `V_v` from [T1D-G, (G.3)](CR-INT-T1D-G.md). Exact-fraction interval arithmetic, with rigorous elementary series remainders for all transcendental constants, gives

\[
 \boxed{\mathbf{1.21756851\le V_v\le1.22714789}.}\tag{H.1}
\]

In particular, `6/5 < V_v < 5/4`. The new numerical enclosure is for the **actual fixed C∞ mollified profile**, not the triangle surrogate, and is independent of floating-point integration. The enclosure is still much wider than a direct high-precision quadrature could achieve, but it is fully rigorous under the precisely stated shape.

## 4. An explicit smooth, doubly pole-null, physically orthogonal mixed-coupling witness

Define the unit **three-cell** source

\[
 \boxed{y=\frac{U_{\ell,x_0+\log2}v-s U_{\ell,x_1+\log2}v+U_{\ell,x_2+\log2}v}{\sqrt{N_y}},\quad N_y=2+s^2=22/3.}\tag{H.2}
\]

Its three supports lie strictly in `(-2,2)`. It is a translated copy of the original two-Mellin-pole-annihilated generator, so both Mellin pole moments vanish **exactly**. Every shifted cell is disjoint from all four original selected cells, because their smallest center separation is `log(3/2)>ell`. Hence `||y||_2=1` and **`y` is an authentic smooth element of `F`**; the noncompact windowed-pole projection never needs to be used to construct the witness.

**Complete prime census:** For original center `x_i` and shifted center `x_j+log2`, the separation is `|log(2*3^(j-i))|`. The only exact prime-power equality is `n=2` with `j=i`; enumerating all genuine prime powers `2<=n<=54` and pairs `i in {0,1,2,3}`, `j in {0,1,2}` shows that every other prime-power displacement differs by at least `log(27/26)>ell`. Thus **all other original prime correlations vanish** for this pair, not just their leading terms.

Their exact prime mixed contribution is

\[
 q_{\rm prime}(u_\sigma,y)=-\kappa_2\frac{D_\sigma}{\sqrt{N_\sigma N_y}},\quad
 \kappa_2=\frac{\log2}{\sqrt2},\quad
 D_+=2-2s+s^2,\quad D_-=2+2s+s^2.
\]

The pole contribution is zero and the zeroth-order Gamma contribution vanishes from physical disjointness. Every surviving Gamma cross integration has separation at least `log(3/2)-ell`; since `w_gamma(t)=exp(-t/2)/(1-exp(-2t))<2` there, Cauchy–Schwarz gives the **profile-independent bound**

\[
 |q_\Gamma(u_\sigma,y)|<\frac{2\ell\,\|c_\sigma\|_{\ell^1}(2+s)}{\sqrt{N_\sigma N_y}},\quad
 \|c_+\|_1=2+2(s-1),\quad\|c_-\|_1=2+2(s+1).
\]

The exact rational interval audit gives:

| One-dimensional probe | Exact prime magnitude (outward interval) | Gamma absolute upper | Proven full mixed magnitude lower |
|---|---:|---:|---:|
| `|q(u_+,y)|` | `[0.21085901, 0.21085902]` | `<0.00063091` | `>0.21022810>21/100` |
| `|q(u_-,y)|` | `[0.44245336, 0.44245337]` | `<0.00056106` | `>0.44189231>11/25` |

Because `y` is a unit vector **inside the correct physical complement**,

\[
 \boxed{\beta_+^2>\left(\frac{21}{100}\right)^2=0.0441,\qquad
 \beta_-^2>\left(\frac{11}{25}\right)^2=0.1936.}\tag{H.3}
\]

This is an **actual frozen-width**, two-direction nonvanishing residual certificate, not a limit theorem or a guessed finite quadrature value. It uses only the original full signed Weil pairing and an explicit smooth witness. The bounds also apply to any **real** bump profile of width at most `10^-4` with this cell geometry, because the argument uses only its support, normalization and pole-annihilator, not its shape or sign.

Finally, combining (H.3) with the **separately established T1C** interval `mu_- < 14.18874262 < 71/5` yields the rigorous finite-width Schur-penalty floor

\[
 \boxed{\|\mathcal R_{10^{-4}}\|\ge\frac{\beta_-^2}{\mu_-}>\frac{121}{8875}\simeq 0.0136338.}\tag{H.4}
\]

This is a **one-sided lower bound**, not the requested full numerical enclosure. In particular, the Schur penalty is *provably not zero* at the original finite width, despite T1D-F's asymptotic vanishing as width tends to zero.

## 5. Cross-canonical meaning, remaining tests and honest stopping point

- **Certified now:** actual smooth-profile Fourier-log variance lies in (H.1); both original signed Weil mixed functionals are nonzero at width `1e-4` by (H.3); conditioned on the independent existing T1C positive Rayleigh enclosures, the physical projected Schur norm is at least (H.4). A genuine smooth pole-null complement witness makes the inequalities falsifiable without an infinite basis. Exact arithmetic controls the `2`, `3`, and `13` displacement near-collision margins; there is no silent omission of `n<=54` prime powers.
- **Not yet certified:** tight *two-sided* enclosures for the actual finite-width `beta_±²` (these still require numerical/integral evaluation of `K_2u` and its physical pole projections, with rigorous errors). The profile variance interval does **not** supply those finite-width residual values by substituting it into the *limiting* T1D-G formula. That substitution would be invalid.
- **Open:** any strictly positive complementary coercivity bound `q(y,y) >= d||y||²` on all of `F`; any Schur domination `d >= ||R||`; the complete pole-direction transfer; semilocal H2′/R1/R2; equality of completed Suzuki/Connes domains, global Weil positivity or RH. The new **positive lower bound on the Schur penalty** makes any proposed complementary coercivity constant `d` falsifiably testable: a claim `0<d<=121/8875` alone would be **insufficient** to establish the required Schur domination from its norm comparison.
- **Reproduction:** `python3 src/audits/canonical_return_t1d_h_mellin_probe.py`, using exact `fractions.Fraction` rational intervals and elementary remainder estimates for `log`, `pi`, and Euler's constant. The trig Mellin identity, Fourier probability-weight inequality, and smooth-source pairing arguments are human analytic inputs and are **not independently proved by the script**. All printed rational interval endpoints are rounded outward. No numerical quadrature/FFT is used, and there is no Lean certificate.

**Next CR-INT-T1D-I:** obtain **two-sided** certified finite-width `beta_±²` intervals via a shape-resolved Gamma multiplier integral and explicit prime-translated cells, including quantitative control of pole projections. Separately, test the actual complement spectrum or state a rigorous form-bound incompatibility with the hypotheses of the semilocal conditioning framework. Continue to maintain no-transfer status until a correct complement theorem is established.
