# CR-INT-T1D-W — Full mollified central Gamma matrix intervals and a certified two-mode Schur/Ritz window

**Canonical Return · November 23, 2026 (America/Los_Angeles) · original signed Weil source.**  
**Disposition:** **NEW ANALYTIC FINITE-SOURCE CERTIFICATES / exact rational and symbolic AUDIT-PASS.** This note computes **outward-rigorous intervals for every complete Gamma matrix entry**, propagates smoothing rigorously, proves the exact original prime-2/5 selected-to-central coupling census and a uniform bound on the *true* selected-source Schur correction, and therefore certifies **both Ritz eigenvalues of each genuine central two-dimensional subspace are below 1**. It does **not** certify that the top eigenvalue of either complete infinite-dimensional parity operator is below 1, nor calculate the exact two-dimensional Schur cross-Gram entries, `eta_N`, the Connes/Sonin intertwiner, RH, or Phase IV positivity.

## 1. Freeze the source and select an admissible mollifier

Preserve the *original complete signed Weil form* of [Suzuki, arXiv:2606.09096v4, §1.1 and §2.4](https://arxiv.org/html/2606.09096v4), the physical Hilbert space `H=L²((-2,2),dx)`, the original narrow T1C selected modes `u_+` and `u_-` at `ell=10^-4`, and the physical pole carriers `e^{±x/2}`. The source is

\[
\mathcal D^0=\{f\in C_c^\infty(-2,2):\int f(x)e^{\pm x/2}\,dx=0\},\quad
F_s=\mathcal D^0\cap\operatorname{span}\{u_+,u_-\}^{\perp},
\]

and the source-defined Schur remainder is `q[y]-R[y]`, with `R>=0` the rank-two actual Weil cross-block penalty from T1D-E. The compact Birman–Schwinger operator is the T1D-R operator
\[
\mathscr B_F=D_F^{-1/2}(6I+A_F+\mathcal R)D_F^{-1/2},\quad
D_F=\Gamma_F+6I\ge5I/8.
\]

This is not identified with the semilocal Sonin operator. The full complement criterion remains `sup sigma(B_F)<=1`.

Fix the *unchanged T1D-V polynomial carriers* for `k=3,4` and parity `sigma=+` (even) or `-` (odd):
\[
 b_{+,k}(x)=(1-16x^2)^k\mathbf1_{|x|<1/4},\qquad
 b_{-,k}(x)=x(1-16x^2)^k\mathbf1_{|x|<1/4},\qquad
 f_{\sigma,k}=(\partial_x^2-1/4)b_{\sigma,k}.
\]

The raw `f` are in compact `H¹` and have explicit rational piecewise-polynomial formulas. T1D-V allowed **any** even, nonnegative `C∞` probability mollifier supported in `[-10^-4,10^-4]`. This note fixes a narrower **allowed choice**
\[
\epsilon_*=10^{-7},\quad
\rho_*(x)=Z_*^{-1}\exp\bigl(-[1-(x/\epsilon_*)^2]^{-1}\bigr)\mathbf1_{|x|<\epsilon_*},\quad
 y_{\sigma,k}=f_{\sigma,k}*\rho_*.
\]

Every `y` is an authentic compact smooth test supported in `[-1/4-eps*,1/4+eps*]`; its two exponential Mellin moments vanish **exactly** by integration by parts. Its physical support is disjoint from the selected T1C cells. Thus `Y_±=span{y_{±,3},y_{±,4}}` are **two-dimensional genuine physical-complement subspaces**. Reflection gives exact even/odd separation. As in V, every self- or mutual prime translation within `Y_+⊕Y_-` vanishes: the support diameter is strictly below `log 2`. Hence `q(y,z)=gamma(y,z)` on each central subspace, retaining the original signed prime operator with **zero actual pairing**, not replacing it by a different form.

## 2. Complete Gamma matrix entries: an exact polynomial-integral representation

Write `M_ij=<f_i,f_j>`, `K_ij=<f_i',f_j'>` for the **unmollified** polynomial carriers of a fixed parity, both evaluated **exactly as rationals** by the T1D-V checker. Let
\[
R_{ij}(t)=\int_{-1/4}^{1/4-t} f_i(x)f_j(x+t)\,dx,
\quad 0\le t\le1/2.
\]

This is a **rational polynomial** of degree at most 19, with `R_ij(0)=M_ij` and `R'_ij(0)=0` (both profiles have the same parity). Its coefficients are recomputed exactly in the T1D-W audit, not inferred from FFT data. Suzuki's full signed Gamma distribution yields, with `C_Gamma=log(4 pi)+gamma_Euler` and `w(t)=e^{-t/2}/(1-e^{-2t})`,
\[
\boxed{\gamma(f_i,f_j)=\left[\log\coth(1/4)-C_\Gamma\right]M_{ij}
+2\int_0^{1/2}\bigl(e^{-t/2}M_{ij}-R_{ij}(t)\bigr)w(t)\,dt.}\tag{W.1}
\]

This is the **entire** original Gamma contribution, including its singular diagonal and outside-support tail. Expand `w(t)=sum_{m>=0}exp[-(2m+1/2)t]`. Each term in (W.1) is an elementary integral of an exact rational polynomial times an exponential with **rational exponent**. The checker computes 2,048 terms by interval recurrence for
\[
 I_0(\alpha)=\frac{1-e^{-\alpha/2}}{\alpha},\qquad
 I_j(\alpha)=\frac{j I_{j-1}(\alpha)-(1/2)^je^{-\alpha/2}}{\alpha}.
\]

Every transcendental endpoint is enclosed by outward-rational Taylor/alternating-series bounds, and arithmetic is fixed-scale directed integer interval arithmetic (`S=10^31`), with no unverified numerical integration. For the omitted tail, reflection gives the **analytic** bound
\[
 |M_{ij}-R_{ij}(t)|\le\frac{t^2}{2}\|f_i'\|_2\|f_j'\|_2
\le\frac{t^2}{4}(K_{ii}+K_{jj}).
\]

Since `1-e^{-2t}>=2t/(1+2t)`, letting `alpha=2N+1/2`, `N=2048`, and `D=(K_ii+K_jj)/2`, the remaining signed Gamma series has absolute value at most
\[
 \boxed{\frac{|M_{ij}|}{2\alpha}
 +\frac{2|M_{ij}|+D}{2\alpha^2}
 +\frac{2D}{\alpha^3}.}\tag{W.2}
\]

The exact rational `M_ij` are positive for these chosen profiles. **Independent numerical FFT sanity checks agree with the enclosed raw Gamma ranges, but are not used to establish any interval.**

## 3. Rigorously transport raw Gamma intervals to the genuine smooth tests

Put `delta_i=y_i-f_i` and `eps*=10^-7`. Young's inequality and translation estimates show
\[
 \|\delta_i\|_2\le\epsilon_*\|f_i'\|_2,\qquad
 \|\delta_i'\|_2\le2\|f_i'\|_2.
\]

Use the *same* positive shifted Gamma form `d[f]=gamma[f]+6||f||²` and T1D-U's digamma/Jensen upper estimate. Its increasing concave Fourier envelope has an increasing mass-energy perspective; the checker verifies the conservative shape-independent inequalities
\[
 \boxed{d[\delta_i]<30\epsilon_*^2K_{ii},\qquad d[f_i]<12M_{ii}.}\tag{W.3}
\]

The real positive-form Cauchy–Schwarz inequality therefore bounds the **actual** Gamma-entry perturbation by
\[
 |\gamma(y_i,y_j)-\gamma(f_i,f_j)|
 \le \sqrt{d[\delta_i]d[f_j]}+\sqrt{d[\delta_j]d[f_i]}
 +\sqrt{d[\delta_i]d[\delta_j]}+6|\langle y_i,y_j\rangle-M_{ij}|,
\]

with a separate exact mass perturbation bound
\[
 |\langle y_i,y_j\rangle-M_{ij}|
 \le\epsilon_*(\sqrt{K_{ii}M_{jj}}+\sqrt{K_{jj}M_{ii}})
 +\epsilon_*^2\sqrt{K_{ii}K_{jj}}.
\]

All terms are enclosed by *integer square-root rational bounds*, so the following table is a **real smooth-source** bound, not an unsmoothed substitution. Each number is the Gamma entry `gamma(y_i,y_j)/sqrt(M^0_ii M^0_jj)` using the exact known **raw** masses in its denominator; the small difference between these masses and the actual smooth mass matrix is included in the subsequent generalized-pencil certificates.

| Gamma entry, outward rational interval | Even | Odd |
| --- | ---: | ---: |
| `(3,3)` | `[0.897400, 0.897857]` | `[1.189019, 1.189549]` |
| `(3,4)` | `[0.852240, 0.852680]` | `[1.100258, 1.100763]` |
| `(4,4)` | `[0.982180, 0.982632]` | `[1.249386, 1.249902]` |

Denote the **fully smoothed** physical and Gamma Hermitian `2x2` matrices by `M_sigma^*` and `G_sigma^*`. By independent directed rational entry bounds and the **strict 2×2 Sylvester determinants**, the checker establishes the genuine **Loewner matrix inequalities**
\[
\boxed{\frac{17}{20}M_+^*\prec G_+^*\prec\frac{31}{20}M_+^*,\qquad
\frac{23}{20}M_-^*\prec G_-^*\prec\frac95M_-^*.}\tag{W.4}
\]

In particular, **every** nonzero smooth source vector in either `Y_±` has positive complete original Weil self-form `q[y]=gamma[y]>0`. These are **two-sided matrix bounds for the actual smooth profiles**, not simply Jensen diagonal upper bounds.

## 4. Complete selected-to-central Weil coupling: only primes 2 and 5

Let `ell=10^-4`, `d=log3`, `s=4/sqrt3`; the four original selected cells have centers `x_j=(j-3/2)d`, `j=0,1,2,3`, with narrow smooth bump `phi` supported within `ell/2` of each center. Their coefficient rows are
\[
 (c_{+,j})=(1,1-s,1-s,1),\qquad
 (c_{-,j})=(1,-1-s,1+s,-1),
\]

and `u_±=sum_j c_{±,j} phi(x-x_j)/sqrt(J N_±)`, where `J=||phi||²`, `N_±=(44∓16sqrt3)/3`.

Checking **all 24 genuine prime powers**, both shift orientations, and both source support widths gives **exactly four selected-cell/central-support arithmetic collisions**:

| Selected cell | Active original prime power | Shift |
| --- | --- | --- |
| `j=0` | `5` | `+log5` |
| `j=1` | `2` | `+log2` |
| `j=2` | `2` | `−log2` |
| `j=3` | `5` | `−log5` |

No other prime power can connect these supports, including the original 3/9/27 cells. **All** those excluded pairings are zero by a strict rational logarithmic separation test. For a fixed-parity central `y`, symmetry reduces the *exact* prime cross-form to
\[
\boxed{q_{\rm prime}(u_\sigma,y)=-\frac{2}{\sqrt{JN_\sigma}}
 \left[\kappa_5\int\phi(z)y\Bigl(\log5-\tfrac32\log3+z\Bigr)dz
 +c_{\sigma,1}\kappa_2\int\phi(z)y\Bigl(\log2-\tfrac12\log3+z\Bigr)dz\right].}\tag{W.5}
\]

This expression is exact for the *original* narrow selected bump and the *actual* central smooth source. **The integrals in (W.5) are not evaluated numerically or separately enclosed entrywise in this milestone.** Instead, the source-valid Sobolev inequality `||y||infty<=sqrt(2||y||||y'||)`, the T1D-V uniform parity derivative ratios `<38,<46`, and `int phi <=sqrt(ell J)` yield a uniform bound for the prime part.

For the **original Gamma cross part**, the supports are disjoint and their minimum distance exceeds `29/100`. The original disjoint-support kernel satisfies `w(t)<2` for `t>=29/100`, checked by rational exponential bounds. With `||y||_1<=sqrt(1/2+2epsilon*)||y||_2` and the exact selected coefficient `l1` masses, combining prime and Gamma estimates proves
\[
\boxed{|q(u_+,y)|<\frac{2}{15}\|y\|_2\quad(y\in Y_+),\qquad
|q(u_-,y)|<\frac{3}{25}\|y\|_2\quad(y\in Y_-).}\tag{W.6}
\]

The independent rational audit reports stronger internal upper bounds `0.129918` and `0.116814` respectively, but **the certified simple fractions** in (W.6) are used in every theorem below. The full signed Weil form is retained in (W.6); the source-defined positive Schur correction is consequently bounded on these specific two-spaces by
\[
\boxed{0\le\mathcal R[y]<\frac4{1575}\|y\|_2^2\ (y\in Y_+),\qquad
0\le\mathcal R[y]<\frac9{4375}\|y\|_2^2\ (y\in Y_-).}\tag{W.7}
\]

This uses only `mu_±>7` from T1C, **not** the disputed T1D-I global residual-norm enclosure. The two parity cross channels decouple exactly.

## 5. Genuine two-mode **Ritz intervals below 1**, not a complete-spectrum result

For each parity and for every nonzero `y` in the genuine smooth two-dimensional `Y_sigma`, the *actual* compact Birman–Schwinger Rayleigh quotient at `z=D_F^{1/2}y` is
\[
 \frac{\langle\mathscr B_Fz,z\rangle}{\|z\|^2}
 =\frac{6\|y\|^2+\mathcal R[y]}{6\|y\|^2+\gamma[y]},
\]

because the original **prime self-pairing is exactly zero**. Equations (W.4) and (W.7) therefore give strict **two-sided** uniform Rayleigh bounds on the entire central subspaces:
\[
\boxed{
\frac{120}{151}<\lambda_{j,+}^{\rm Ritz}<\frac{37816}{43155}<1,
\qquad
\frac{10}{13}<\lambda_{j,-}^{\rm Ritz}<\frac{105036}{125125}<1,
\quad j=1,2.}\tag{W.8}
\]

Here `lambda_{j,sigma}^Ritz` means the generalized eigenvalues of the authentic two-dimensional pencil
\[
(6M_\sigma^*+R_\sigma^*)v=\lambda(6M_\sigma^*+G_\sigma^*)v,
\]

where `R_sigma^*` is the **unknown entrywise, but rigorously norm-bounded** true Schur matrix induced by the original selected sources. The bounds certify the matrix pencil's two eigenvalues without pretending to compute that rank-one matrix entrywise. For orientation the even interval is `(0.7947019,0.8762832)` and the odd interval is `(0.7692307,0.8394486)`.

The actual **source Schur remainder** is thus strictly positive on all four central directions:
\[
\boxed{
 q[y]-\mathcal R[y]>
 \begin{cases}
  (17/20-4/1575)\|y\|^2,&y\in Y_+,\\
  (23/20-9/4375)\|y\|^2,&y\in Y_-.
 \end{cases}}
\]

The selected space E already has `mu_±>7`; completing the finite Schur square also proves **the original signed Weil form is positive on the new six-dimensional space `E⊕Y_+⊕Y_-`**. This is different from T1D-K's highly localized six-dimensional positive compression: these added directions are broad-central low-energy tests. **Nothing here bounds `sup sigma(B_F)` from above**: the infinite-dimensional orthogonal remainder and its cross-block couplings remain outside this calculation.

## 6. Verification, provenance and hard stopping condition

The companion `src/audits/canonical_return_t1d_w_full_central_gram.py` (standard Python plus the repository's SymPy) rederives every raw profile and polynomial correlation exactly, encloses rational transcendentals, sums the first `2048` complete Gamma-kernel terms via directed integer rational intervals, adds the explicit infinite-tail theorem (W.2), propagates nonzero mollifier error via (W.3), enforces both strict 2×2 Loewner/Sylvester inequalities, enumerates all 24 prime-power selected–central shifts, and verifies the exact combined Schur operator-norm fractions and two-sided finite Ritz windows. Its printed decimal intervals are **rounded outward by integer arithmetic**, not inferred from floating-point quadrature. **Analytical inputs**, not compiler outputs, include Suzuki's Gamma kernel, Plancherel/Cauchy–Schwarz, the digamma/Jensen envelope, and smooth form-domain preservation.

**Proved on these finite source spaces (analytically, pending external review):** complete Gamma `2×2` interval matrices on the actual smooth test source, exact prime-2/5 selected cross-census, source-defined rank-one-per-parity Schur bound, positive central Schur remainder, both authentic finite Ritz eigenvalues below 1 and six-dimensional broad-plus-selected positivity.

**NOT proved:** exact entrywise `R_sigma^*`, global top eigenvalues `sup sigma(B_{F,±})`, `eta_N`, spectral projector comparison, positive control of the complement outside `Y_+⊕Y_-`, global original Weil positivity, all pole directions or large support, Suzuki–Connes semilocal domain/metric transport, Phase IV/RH. Finite Ritz *upper* bounds must never be mistaken for complete-spectrum upper bounds.

**Next CR-INT-T1D-X:** certify the actual `Y`-to-remainder form Schur coupling (or a `D_F`-spectral projector comparison) rather than merely enlarging finite positive test spaces. Quantify `eta_N` on a genuinely bounded-energy projection or prove a useful form-adapted off-block estimate. A finite two-mode pencil below threshold is not an upper certificate for the full parity operator.
