# CR-INT-T1D-U — Exact central-support prime image and near-threshold parity Ritz witnesses

**Canonical Return · November 23, 2026 (America/Los_Angeles).**  
**Disposition:** **NEW SOURCE-LEVEL ANALYTIC inequalities**, with independent executable exact-rational/SymPy checks of every finite arithmetic and polynomial inequality. We improve the genuine low-frequency trial **lower** bounds for the compact parity spectra and establish a support-sensitive bound on the complete original prime operator acting on the central interval. **No** actual largest eigenvalue, spectral high-projection norm `eta_N`, upper eigenvalue enclosure, complete pole-null Schur positivity, Suzuki–Connes completed-domain correspondence, Lean or RH proof is claimed. Draft PR #2 remains open and unmerged.

## 1. Unchanged original source and spectral obligation

Continue [T1D-T](CR-INT-T1D-T.md) with the original signed Weil form of Suzuki, [*Weil's quadratic form via the screw function*, arXiv:2606.09096v4, §1.1, §2.4](https://arxiv.org/html/2606.09096v4). The physical metric is `H=L²(I,dx)`, `I=(-2,2)`. The smooth exactly two-Mellin-pole-null core is

\[
\mathcal D^0=\{f\in C_c^\infty(I;\mathbb C):\ \int_I f(x)e^{x/2}dx=\int_I f(x)e^{-x/2}dx=0\}.
\]

Retain the original selected physical-orthonormal pair `E=span{u_+,u_-}` at width `ell=10^{-4}`, both physical pole carriers `P=span{e^{x/2},e^{-x/2}}`, and the physical complement `H_F=(E+P)^perp`, with `F_s=D^0 cap H_F`. The original full prime quadratic contribution is **negative**:

\[
A=\sum_{n=p^k\le54}\kappa_nP_I(T_{\log n}+T_{-\log n})P_I,\qquad
\kappa_n=\Lambda(n)/\sqrt n>0,\quad q_{\rm prime}(f,f)=-\langle Af,f\rangle.
\tag{U.1}
\]

On `H_F`, T1D-R's source-level specialization of Suzuki's canonical compactness theorem gives the positive closed shifted Gamma form `d[f]=gamma[f]+6||f||²`, operator `D_F>=5/8`, and compact self-adjoint

\[
\mathscr B_F=D_F^{-1/2}(6I+A_F+\mathcal R)D_F^{-1/2},\quad A_F=\Pi_F A\Pi_F,
\quad \mathcal R\ge0.
\tag{U.2}
\]

The unresolved full smooth complement Schur criterion is `sup sigma(B_F)<=1`, separated into even and odd parity. This construction is *not* identified with the Andrews–Connes–La Gatta semilocal Sonin operator [arXiv:2611.13129v1](https://arxiv.org/html/2611.13129v1); no transport between their domains, metrics or forms is inferred.

## 2. A genuine operator estimate on central supported functions

Put `a=51/200`, `J=(-a,a)`, and let `P_J` denote physical multiplication by `1_J`. This is exactly the support window for T1D-T's compact smooth broad parity probes. Every nonzero source `f` supported in `J` satisfies, **on the original unmodified prime operator**:

\[
\boxed{P_JAP_J=0,\qquad \|AP_J\|_{L^2(I)\to L^2(I)}^2<\frac{61}{5},\qquad \|AP_J\|<\frac72.}\tag{U.3}
\]

**Support geometry.** The support of `T_{+/-log n}f` intersects `I` only if `log n<2+a=451/200`. The exact rational exponential check gives `9<exp(451/200)<11`, so of the **24 genuine prime powers** in (U.1), exactly **seven** can contribute to `AP_J`:

\[
\boxed{n\in\{2,3,4,5,7,8,9\}.}\tag{U.4}
\]

Here `4,8,9` are retained as genuine prime powers; `6` is not included since `Lambda(6)=0`. Positive and negative translated images lie on disjoint sides of zero since `log2>a`; moreover `2a=51/100<log2` proves `P_JAP_J=0`.

Two different positive shifts `log n, log m` can overlap only when `|log(m/n)|<2a`. A *strict exact rational* exponential enclosure proves

\[
\frac85<e^{51/100}<\frac53.
\]

The complete unordered overlap graph among the seven contributing shifts therefore has **exactly eight** edges:

\[
\boxed{\mathcal E=\{(2,3),(3,4),(4,5),(5,7),(5,8),(7,8),(7,9),(8,9)\}.}\tag{U.5}
\]

All other pairs have disjoint supporting translated copies. The compressed images may additionally be truncated at the physical boundary `+/-2`, but truncation can only decrease their norms and cannot introduce new overlaps. Hence Cauchy–Schwarz applied **only to the permitted overlaps**, including both signs, yields for arbitrary `f` supported in `J`:

\[
\begin{split}
\|Af\|_2^2
&\le 2\left(\sum_{n\in S}\kappa_n^2
+2\sum_{(m,n)\in\mathcal E}\kappa_m\kappa_n\right)\|f\|_2^2,\\
2\left(\sum_{n\in S}\kappa_n^2
+2\sum_{(m,n)\in\mathcal E}\kappa_m\kappa_n\right)
&\in(12.132600917327,12.132600917328)<\frac{61}{5}<\frac{49}{4}.
\end{split}\tag{U.6}
\]

The decimal interval is rounded **outward from exact fractions**, not obtained by floating quadrature. Its logarithms and integer square roots are certified by rational series and integer arithmetic. This is a uniform infinite-dimensional **localized arithmetic-image** bound, stronger on the central source than the whole-window `||A||<203/25` of T1D-Q. It neither bounds the full complement nor changes the sign of the arithmetic part of Weil's form. For any compact smooth exactly pole-null central `f`, its **entire self prime contribution vanishes**; `q[f]=gamma[f]`.

## 3. A strong Gamma-energy upper bound from the exact digamma series

T1D-R established the canonical original source Gamma Fourier multiplier

\[
\gamma[f]=\int_{\mathbb R}m_\Gamma(\xi)|\widehat f(\xi)|^2d\xi,\quad
m_\Gamma(\xi)=\operatorname{Re}\psi(\tfrac14+\tfrac{i\xi}2)-\log\pi.
\]

Its exact minimum is `m_0=-gamma_Euler-pi/2-3 log2-logpi`. A direct rational Machin/atanh/harmonic enclosure in the new checker establishes **`m_0<-5372/1000`**. More importantly, the exact digamma-series increment is

\[
m_\Gamma(\xi)-m_0
=\sum_{n=0}^{\infty}
\frac{\xi^2}{(n+1/4)(4(n+1/4)^2+\xi^2)}.
\tag{U.7}
\]

For each `n`, the summand is **increasing and concave in** `v=xi²`. The summand is also decreasing in its continuous `n+1/4` argument. Consequently, for any integer `M>=0`, the discarded sum has the explicit monotone-integral upper bound

\[
\sum_{n=M+1}^\infty\frac{v}{(n+1/4)(4(n+1/4)^2+v)}
\le\frac12\log\left(1+\frac{v}{4(M+1/4)^2}\right).
\tag{U.8}
\]

The finite sum **and** the logarithmic tail envelope are concave increasing in `v`. Let `V=||f'||_2²/||f||_2²` for any nonzero compact `H¹` source `f`. Applying Plancherel and Jensen to (U.7)–(U.8) gives the uniform analytic inequality

\[
\boxed{\frac{\gamma[f]}{\|f\|_2^2}
\le m_0+\sum_{n=0}^{M}\frac{V}{(n+1/4)(4(n+1/4)^2+V)}
+\frac12\log\left(1+\frac{V}{4(M+1/4)^2}\right).}\tag{U.9}
\]

This is a considerably sharper *source-preserving* Gamma upper estimate for moderate-frequency trials than T1D-T's generic `log(1+|xi|)+40` envelope. It is not an equality or a numerical evaluation of the Gamma integral. We use `M=50`, entirely with rational bounds, below.

## 4. New near-threshold parity Ritz lower witnesses

Reuse **the same** T1D-T admissible parity functions; no substitution of a toy profile:

\[
\begin{aligned}
b_e(x)&=(1-16x^2)^5\mathbf1_{|x|<1/4},\\
b_o(x)&=x(1-16x^2)^5\mathbf1_{|x|<1/4},\\
y_\sigma&=(\partial_x^2-\tfrac14)(b_\sigma*\rho),
\end{aligned}\tag{U.10}
\]

where `rho` is any even smooth nonnegative probability mollifier supported in `[-1/200,1/200]`. Both `y_e,y_o` belong to `F_s`: their two Mellin moments vanish identically by integration by parts, and their central support is physically disjoint from `E`. They have opposite reflection parity. The exact SymPy integration of the *unmollified* `y_{0,sigma}=(d²-1/4)b_sigma` gives:

| Probe | `||y_0||²` | `||y_0'||²` |
| --- | ---: | ---: |
| even | `9581387776/2909907` | `3906493153280/2909907` |
| odd | `3621357056/66927861` | `100490301440/2909907` |

Young/Minkowski and `eps=1/200` give

\[
\|y_\sigma\|_2\ge\|y_{0,\sigma}\|_2-\varepsilon\|y_{0,\sigma}'\|_2,\quad
\|y_\sigma'\|_2\le\|y_{0,\sigma}'\|_2.
\]

Unlike T1D-T's loose derivative ratios `<24` and `<30`, direct exact-fraction comparison here sharpens them to

\[
\boxed{\frac{\|y_e'\|_2}{\|y_e\|_2}<23,\qquad
\frac{\|y_o'\|_2}{\|y_o\|_2}<29.}\tag{U.11}
\]

Insert `V_e<23²`, `V_o<29²` into (U.9), taking `M=50` and the rationally certified `m_0<-5.372`. The exact rational check yields

\[
\boxed{\frac{\gamma[y_e]}{\|y_e\|_2^2}<\frac{13}{10}=1.3,\qquad
\frac{\gamma[y_o]}{\|y_o\|_2^2}<\frac{77}{50}=1.54.}\tag{U.12}
\]

For diagnostics only, the rational *upper formula* (U.9), not the actual Gamma integral, lies in the strict outward intervals `(1.298212221565,1.298212221566)` (even), `(1.530310235491,1.530310235492)` (odd).

Now use `P_JAP_J=0` and the exact positive Schur correction `R>=0`. For either of these **authentic full-source** trial directions,

\[
\langle(6I+A_F+\mathcal R)y_\sigma,y_\sigma\rangle
=6\|y_\sigma\|_2^2+\langle\mathcal R y_\sigma,y_\sigma\rangle
\ge6\|y_\sigma\|_2^2.
\]

Their Birman–Schwinger Rayleigh quotient is evaluated on `z_sigma=D_F^{1/2}y_sigma`, using the true positive physical metric. From (U.12), **without** importing any of T1D-I's disputed two-sided residual norm estimates,

\[
\boxed{
\sup\sigma(\mathscr B_{F,+})>\frac{60}{73}\approx0.8219178,
\qquad
\sup\sigma(\mathscr B_{F,-})>\frac{300}{377}\approx0.7957560.
}\tag{U.13}
\]

These are rigorous one-sided **lower** witnesses for the actual compact parity spectral problem, improving T1D-S's `>0.28` and `>0.31`, and T1D-T's `>0.12`. The strict inequalities follow from the strict Gamma upper limits. They do **not** prove the top eigenvalue equals either value, or supply an **upper** bound below `1`. The rank-two `R` was retained by its positive sign, not discarded or claimed to vanish.

## 5. Separate optional numerical diagnostic — not a certificate

For the single explicit even probability mollifier `rho(x) proportional to exp(-1/(1-(200x)^2))` on `|x|<1/200`, a Fourier-grid diagnostic with lengths `16,32,48` and spacings approximately `6e-5` to `9e-5` gives stable *approximate* normalized Gamma energies:

| Broad smooth probe | Exploratory `gamma[y]/||y||²` | Exploratory `||Ay||/||y||` |
| --- | ---: | ---: |
| even | `~1.074326135` | `~1.856834` |
| odd | `~1.327400523` | `~1.844374` |

These figures are **not rigorous intervals**, are not used in any proof above and do **not** calculate `R[y]` or the full compact spectrum. They merely indicate how conservative the certified analytic upper estimates are. The separately supplied optional `tools/diagnostic_t1d_u_fft.py` requires NumPy and SciPy in addition to SymPy, and is deliberately **not** part of the exact-rational regression runner.

## 6. What remains missing and the next falsifier

- **Proven at source level:** for every central supported `L²` vector, the exact original seven-prime-image geometry and uniform `||AP_J||<7/2`; for the actual broad smooth pole-null parity probes, a direct Jensen/digamma upper bound on original Gamma energy, with genuine parity compact-spectrum **lower** witnesses `>60/73`, `>300/377`.
- **Machine checks:** all 24 genuine prime powers and seven surviving channels, exact exponential comparisons (`8/5<e^(.51)<5/3`), exact eight-edge overlap graph, exact rational von Mangoldt coefficient bound, rational Euler/π/digamma constant, exact broad-polynomial `H¹` norms and Jensen upper arithmetic. The analytic Schur/Plancherel/Jensen proof and source domain are written arguments, not Lean execution.
- **Not accomplished:** the actual Gamma eigenbasis, `eta_N=||Q_N(A_F+R)D_F^{-1/2}P_N||`, a lower enclosure for the first omitted Gamma eigenvalue, an upper enclosure of either leading Ritz spectrum below `1`, or positivity of **all** `F_s`. The central image estimate `||AP_J||<7/2` does not imply a bound on `eta_N` for arbitrary low-Gamma-energy eigenvectors, which need not have central support. In particular, no numerical sample above can be promoted to an operator theorem.

**Next CR-INT-T1D-V:** construct a **small, genuine broad-support low-energy Ritz subspace** containing the two newly successful parity probes, and compute rigorous full original signed cross-Grams (including Gamma, prime translations, the selected-space Schur correction, and the physical denominator). In parallel, derive a certified `D_F`-spectral projector comparison for that subspace, or an `A_F+R` commutator/source boundary estimate, before using T1D-T's tail stopping rule. A larger lower witness alone never proves the required upper bound `<=1`. Do not transfer a positivity theorem into Phase IV or identify the Connes/Sonin completion without the corresponding source, form and metric proof.

**Reproduction:** `python3 src/audits/canonical_return_t1d_u_central_prime_gamma.py` (SymPy dependency already in repository requirements); optional numerical diagnostic requires NumPy/SciPy and is not a certification script. **Disposition: DRAFT source-level analysis awaiting independent review.**