# CR-INT-T1D-Y — Direct selected-to-satellite residual bounds and a positive Schur-corrected central–satellite six-space

**Canonical Return · November 23, 2026 (America/Los_Angeles) · independent analytical smooth-source research.**  
**Disposition:** New **DERIVED / ANALYTIC** local source inequality, with executed exact-rational support-census and determinant audit. It removes the *T1D-I global residual estimate* from T1D-X's two-mode **Schur-corrected** positivity test. It **does not** establish positive semidefiniteness on the full physical complement, a certified upper bound on the compact Birman–Schwinger spectrum, the completed Suzuki–Connes domain/metric dictionary, Phase IV positivity, Lean, or RH. T1D-C/W/X source estimates remain explicit analytical imports; the script does not formalize their proofs.

## 1. Original complete Weil form and frozen domain

Keep the original complete signed Weil form of Suzuki, [*Weil's quadratic form via the screw function*, arXiv:2606.09096v4 §1.1, §2.3–2.4](https://arxiv.org/html/2606.09096v4), including both pole terms, every actual von Mangoldt prime power, and the complete Gamma kernel. Fix the *physical* Hilbert space `H=L²(I,dx)`, `I=(-2,2)`, physical inner product linear in its first argument, and the smooth double-Mellin-pole-null core
\[
\mathcal D^0=\{f\in C_c^\infty(I):\int_I f(x)e^{\pm x/2}\,dx=0\}.
\]
Let the original T1C selected directions be `u_+`,`u_-`, physically orthonormal, reflection even/odd, with `q(u_σ,u_σ)=μ_σ>7`, and `E=span{u_+,u_-}`. Put `P=span{e^{±x/2}1_I}`, `\mathscr H_F=(E\oplus P)^\perp`, `F_s=\mathcal D^0\cap\mathscr H_F`. The *true* selected-space Schur form on `F_s` is
\[
\mathfrak s(v,w)=q(v,w)-\langle\mathcal Rv,w\rangle,
\qquad
\mathcal R=\sum_{\sigma=\pm}\frac{|r_\sigma\rangle\langle r_\sigma|}{\mu_\sigma},
\quad r_\sigma=\Pi_{\mathscr H_F}\mathcal K_2u_\sigma.
\]
We use its **exact source identity** `\mathcal R[v]=\sum_σ |q(u_σ,v)|²/μ_σ` for actual smooth `v\in F_s`; no global bound on `||r_σ||` is imported in this note.

Retain T1D-W/X's specified central smooth mollifier: any fixed even, nonnegative, unit-mass `C∞` mollifier supported in `[-ε,ε]`, `ε=10^{-7}`. The central parity tests (k=3 from the W two-mode spaces) are
\[
\begin{split}
 b_+(x)&=(1-16x^2)^3\mathbf1_{|x|<1/4},\\
 b_-(x)&=x(1-16x^2)^3\mathbf1_{|x|<1/4},\\
 y_\sigma&=(\partial_x^2-1/4)(b_\sigma*\rho_\varepsilon),\\
 z_\sigma&=(T_{\log3}y_\sigma+T_{-\log3}y_\sigma)/\sqrt2.
\end{split}
\]
Both `y_σ` and `z_σ` are original smooth pole-null vectors in `F_s`, with opposite parity, pairwise disjoint central/satellite physical supports and `||z_σ||₂=||y_σ||₂`. A circumflex denotes physical normalization by `||y_σ||₂`. This does **not** substitute a raw polynomial for the genuine smoothed Weil source.

## 2. New exact support census: selected cells against both satellite wings

Write `d=log 3`, `ell=10^{-4}`, `s=4/sqrt 3`, `J=||φ||₂²`, `x_j=(j-3/2)d`. The T1C selected cells are `φ_j(x)=φ(x-x_j)`, with `φ` even, nonnegative, supported in an interval of width `ell`, and
\[
\begin{array}{ll}
 u_+={\sum_j c_{+,j}\phi_j}/{\sqrt{J N_+}},&c_+=(1,1-s,1-s,1),\ N_+=(44-16\sqrt3)/3,\\
 u_-={\sum_j c_{-,j}\phi_j}/{\sqrt{J N_-}},&c_-=(1,-1-s,1+s,-1),\ N_-=(44+16\sqrt3)/3.
\end{array}
\]
For a signed prime-power shift `sign*log n`, a selected cell at `x_j` can overlap the satellite wing centered at `wing*d` if and only if
\[
\left|x_j+\operatorname{sign}\log n-\operatorname{wing}\,d\right|<\frac14+\varepsilon+\frac{\ell}{2},
\qquad \operatorname{sign},\operatorname{wing}\in\{-1,1\}.
\]
All possible 384 comparisons (`4 original selected cells × 24 genuine prime powers × 2 shift signs × 2 wings`) are separated **strictly** from this threshold by directed rational logarithm enclosures. Exactly **14** overlap events survive:

| Selected cell `j` | `(n, sign, wing)` surviving |
|---|---|
| 0 | `(2,+,-)`, `(13,+,+)`, `(16,+,+)`, `(17,+,+)`, `(19,+,+)` |
| 1 | `(2,-,-)`, `(5,+,+)` |
| 2 | `(2,+,+)`, `(5,-,-)` |
| 3 | `(2,-,+)`, `(13,-,-)`, `(16,-,-)`, `(17,-,-)`, `(19,-,-)` |

Thus the **complete original** signed prime cross pairing between `u_σ` and satellite `z_σ` has active `n=2,5,13,16,17,19`, **not** merely the previously visible 2/5 terms from central coupling. In particular `16=2^4` has genuine weight `κ₁₆=log 2/4`; no composite of zero von Mangoldt weight enters. The other 18 prime powers vanish *only because their original physical supports do not overlap*. The signed original form retains the negative arithmetic contribution throughout.

## 3. Direct selected-to-satellite Weil cross bound with no T1D-I

Because the frozen `φ` has support length `ell`,
\[
\frac{\|\phi\|_1}{\sqrt J}\le\sqrt\ell=\frac1{100}.
\]
From T1D-V/W's genuine two-space smooth derivative bounds, `||y_+'||₂/||y_+||₂<38` and `||y_-'||₂/||y_-||₂<46`. The one-dimensional Sobolev inequality `||y||_∞²≤2||y||₂||y'||₂`, together with the **disjoint satellite wings**, gives
\[
\|\widehat z_\sigma\|_\infty<\sqrt{r_\sigma},\quad r_+=38,\quad r_-=46.
\]
Using every exactly enumerated event from §2, the **positive absolute majorant** for the true selected/satellite prime pairing is
\[
\boxed{ |q_{\rm prime}(u_\sigma,\widehat z_\sigma)|
\le\sqrt{\ell/N_\sigma}\,\sqrt{r_\sigma}
\sum_{(j,n,\text{sign},\text{wing})\in\mathcal E}\kappa_n|c_{\sigma,j}|. }
\tag{Y.1}
\]
The checker uses rational atanh-logarithm enclosures for every `κ_n=Λ(n)/sqrt n`, outward integer-root bounds, and the exact 14-event list to prove the separate strict majorants
\[
 |q_{\rm prime}(u_+,\widehat z_+)|<0.2286722100,
\qquad
 |q_{\rm prime}(u_-,\widehat z_-)|<0.1870364395.
\]
These printed decimals are conservative upper diagnostics from exact fractions, not floating quadrature.

For the original *complete* Gamma cross contribution, each selected/satellite support separation is at least
\[
 \frac{\log3}{2}-\Big(\frac14+\varepsilon+\frac\ell2\Big)>\frac{29}{100}.
\]
On disjoint supports, Suzuki's original archimedean cross-kernel is `-w(|x-t|)`, `w(t)=e^{-t/2}/(1-e^{-2t})`, and a separately directed alternating-exponential bound proves `w(29/100)<2`. The two satellite wings have total physical support length `1+4ε`, so by Cauchy–Schwarz,
\[
\boxed{|q_\Gamma(u_\sigma,\widehat z_\sigma)|
<2\sqrt{\ell/N_\sigma}\left(\sum_{j=0}^3|c_{\sigma,j}|\right)\sqrt{1+4\varepsilon}.}
\tag{Y.2}
\]
The exact rational script bounds these Gamma terms by `0.0396457856` and `0.0352565038`, respectively. Combining (Y.1) and (Y.2) with all original terms gives
\[
\boxed{|q(u_+,\widehat z_+)|<\frac{27}{100},\qquad
       |q(u_-,\widehat z_-)|<\frac{23}{100}.}\tag{Y.3}
\]
**No T1D-I global `||R||` or global `||r_σ||` estimate is used.** Opposite-parity selected-to-satellite pairings vanish *exactly*. Both Mellin pole contributions vanish by original pole-null source constraints.

## 4. Explicit local Schur penalty: satellite diagonal and mixed term

From T1C `μ_±>7`, the exact selected-source Schur identity and (Y.3),
\[
\boxed{\mathcal R[\widehat z_+]<\frac{729}{70000},\quad
       \mathcal R[\widehat z_-]<\frac{529}{70000}.}\tag{Y.4}
\]
Independently, T1D-W certifies **uniformly across each central parity two-space**, and in particular for `\widehat y_σ`,
\[
 |q(u_+,\widehat y_+)|<2/15,\qquad |q(u_-,\widehat y_-)|<3/25,
\]
so
\[
\boxed{|\langle\mathcal R\widehat y_+,\widehat z_+\rangle|<\frac9{1750},\qquad
|\langle\mathcal R\widehat y_-,\widehat z_-\rangle|<\frac{69}{17500}.}\tag{Y.5}
\]
These are direct local *rank-one-per-parity* bounds, not applications of T1D-I or its disputed full-complement residual enclosures.

## 5. Actual Schur-corrected positive determinant in BOTH parity blocks

Import from the independent, executed T1D-W/X **original full signed Weil** Gram enclosures on the physical unit central/satellite pair:
\[
\begin{array}{c|ccc}
\sigma&q[\widehat y]&q[\widehat z]&|q(\widehat y,\widehat z)|\\\hline
+&>89/100&>1159/1000&<\sqrt2\,(69/100)\\
-&>118/100&>1069/1000&<\sqrt2\,(72/100)
\end{array}
\]
These **already include all 24 original prime powers and the full Gamma term** (the prime support terms are exactly those surviving the support census). Define the true parity Schur `2×2` matrix
\[
\mathfrak S_\sigma=\Big(q(v_i,v_j)-\langle\mathcal Rv_i,v_j\rangle\Big)_{v_i,v_j=\widehat y_\sigma,\widehat z_\sigma}.
\]
Use (Y.4)–(Y.5), and T1D-W's `R[\widehat y_+]<4/1575`, `R[\widehat y_-]<9/4375`. With the strictly rational bound `sqrt2<283/200`, a completely outward **rational determinant lower bound** is
\[
\boxed{\begin{aligned}
\det \mathfrak S_+
&>\left(\frac{89}{100}-\frac4{1575}\right)
  \left(\frac{1159}{1000}-\frac{729}{70000}\right)
  -\left(\frac{283}{200}\frac{69}{100}+\frac9{1750}\right)^2
>\frac1{20},\\
\det \mathfrak S_-
&>\left(\frac{118}{100}-\frac9{4375}\right)
  \left(\frac{1069}{1000}-\frac{529}{70000}\right)
  -\left(\frac{283}{200}\frac{72}{100}+\frac{69}{17500}\right)^2
>\frac15.
\end{aligned}}\tag{Y.6}
\]
The script verifies both strict inequalities **exactly as rational numbers**, with extra lower margins approximately `0.055996` and `0.204316` before replacing them by the conservative simple fractions. Their diagonal entries are positive. Therefore
\[
\boxed{\mathfrak s\big|_{\operatorname{span}\{\widehat y_+,\widehat z_+,\widehat y_-,\widehat z_-\}}>0,}\tag{Y.7}
\]
and by the exact original selected-space Schur identity, the **original full signed Weil form** is strictly positive on the **six-dimensional** smooth pole-null span `E ⊕ span{\widehat y_+,\widehat z_+,\widehat y_-,\widehat z_-}`. This is a *different* six-dimensional positive source from T1D-K's narrow prime-2/4 probe span. The conclusion requires no new global residual-norm hypothesis.

This does **not** imply all of the central two-mode spaces `Y_±` together with the satellites have been checked: only `k=3` central directions in each parity have been used here. It does not imply positivity on `F_s` or on either full compact-spectrum parity sector.

## 6. Exact remaining second Schur gate, audit scope, next obligation

Let `L=span{\widehat y_+,\widehat z_+,\widehat y_-,\widehat z_-}` and `Z_s=F_s∩L^\perp` (physical orthogonal complement **within the smooth source**). The positive finite Schur matrix `S_L` has just been certified. For `w∈Z_s`, let `b(w)=(\mathfrak s(v_i,w))_i`. The **still missing** second Schur requirement is
\[
\boxed{\mathfrak s(w,w)-b(w)^*S_L^{-1}b(w)\ge0\qquad \forall w\in Z_s.}\tag{Y.8}
\]
No lower bound on the full remainder or certified high-energy projection norm is proved by this finite determinant. T1D-T's `eta_N`, the top compact Birman–Schwinger spectral values, and an upper enclosure below one remain **OPEN**. Do not promote this result into the separate Weil-Lab Phase IV accepted proof jurisdiction.

**Exact reproduction:** `python3 src/audits/canonical_return_t1d_y_satellite_schur.py` (Python standard library only); rerun T1D-W and T1D-X independently for the imported source Gamma / original signed-Weil Gram estimates. The new audit verifies all 384 original signed-prime overlaps, 14 survivor identities, rational log and root enclosures, disjoint-support Gamma bound, local selected-source cross inequalities, rank-one-parity Schur correction, and strict rational positive determinants. The mathematical Suzuki kernel and one-dimensional Sobolev estimate are *analytic arguments*, not symbolic output. **Status: analytical research, pending independent review; no Lean/CI/RH claim.**

**Next — CR-INT-T1D-Z:** exploit the now-certified positive 4D complement Schur block to bound its *genuine cross-operator* into `Z_s`, perhaps through finite-rank representers and the complete Gamma/prime kernel; derive a nontrivial lower estimate for the Schur-corrected remainder. Without **both** cross control and remainder coercivity, another finite positive matrix cannot establish the missing whole-complement theorem.