# CR-INT-T1D-V-WIDE — Four-source parity Ritz pencils, exact prime-free internal block, and sharper compact spectral lower witnesses

**Canonical Return · November 23, 2026 (America/Los_Angeles) · companion to the existing [T1D-V](CR-INT-T1D-V.md) research note; original signed Weil source.**  
**Disposition:** new **source-level analytical lower spectral witnesses**, with an executed exact-rational/SymPy audit of the finite polynomial reference Grams, original-prime support cancellation, sharp smooth-energy inequalities, and rational Gamma/Jensen substitutions. A complete **symbolic source-defined** 4×4 Ritz pencil is derived in each parity sector, but the **actual smooth Gamma and selected-space cross entries are not numerically enclosed**, and the 4×4 top eigenvalues are **not certified**. An optional FFT diagnostic evaluates **the Gamma-only pencil**, explicitly omitting the nonnegative rank-one Schur term; it is not a certificate. Not global complement positivity, Suzuki–Connes operator transport, Lean, or RH.

## 1. Fixed original signed source, complementary space, and source expansion

Preserve Suzuki's complete original signed Weil functional ([arXiv:2606.09096v4 §1.1, §2.4](https://arxiv.org/html/2606.09096v4)) on \(I=(-2,2)\), with physical \(H=L^2(I,dx)\), and the original T1C source \(E=\operatorname{span}\{u_+,u_-\}\) of width \(10^{-4}\). The smooth two-pole-null core is
\[
\mathcal D^0=\left\{f\in C_c^\infty(I):\int_I f(x)e^{x/2}\,dx=\int_I f(x)e^{-x/2}\,dx=0\right\},
\qquad F_s=\mathcal D^0\cap E^\perp.
\]
The physical Hilbert complement is \(\mathscr H_F=(E\oplus P)^\perp\), \(P=\operatorname{span}\{e^{x/2},e^{-x/2}\}\) (windowed, not Weil tests). The original signed arithmetic contribution is \(-\langle Af,g\rangle\), with \(A=\sum_{n=p^k\le54}\frac{\log p}{\sqrt n}P_I(T_{\log n}+T_{-\log n})P_I\). No primes or Gamma terms are omitted from the *definition* of the source.

T1D-R established the compact *physical* Schur gate, relative to Suzuki's pre-existing compact-embedding foundation: \(D_F=\Gamma_F+6I\ge5/8\) is the positive shifted Gamma form operator, \(\mathcal R=\sum_{\sigma=\pm}|r_\sigma\rangle\langle r_\sigma|/\mu_\sigma\ge0\) is the actual selected-space penalty, and
\[
\mathscr B_F=D_F^{-1/2}(6I+A_F+\mathcal R)D_F^{-1/2},\qquad
q|_{F_s}-\mathcal R\ge0\iff\sup\sigma(\mathscr B_F)\le1.
\]
Its even and odd parts are independent. T1D-U's existing *specific* smooth broad source has support in \((-51/200,51/200)\) and produced lower witnesses \(>60/73\), \(>300/377\). This V-WIDE companion enlarges the source space and uses the original complete form on it.

## 2. Four genuine smooth pole-null vectors *per parity*, including the T1D-U direction

Put \(a=69/200\), \(\varepsilon=1/1000\) and choose the **explicit** even nonnegative \(C^\infty\) probability bump
\[
\rho_\varepsilon(t)=Z_\varepsilon^{-1}\exp\!\left(-\frac{1}{1-(t/\varepsilon)^2}\right)\mathbf1_{|t|<\varepsilon},
\quad Z_\varepsilon=\int_{-\varepsilon}^{\varepsilon}\exp\!\left(-\frac{1}{1-(t/\varepsilon)^2}\right)dt.
\]
For \(p\in\{0,1\}\) and \((a_j,k_j)=(a,3),(a,4),(a,5),(1/4,5)\), define the zero-extended polynomial
\[
b_{p,j}(x)=(x/a_j)^p\bigl(1-(x/a_j)^2\bigr)^{k_j}\mathbf1_{|x|<a_j},\quad
\boxed{f_{p,j}=(\partial_x^2-1/4)(b_{p,j}*\rho_\varepsilon)}.\tag{V.1}
\]
Each \(b_{p,j}\in H_0^3(-a_j,a_j)\) (indeed the degree-3 edge factor ensures its second weak derivative is in \(H^1\)), and its convolution is genuine compact \(C^\infty\). Integration by parts using \((\partial_x^2-1/4)e^{\pm x/2}=0\) proves both Mellin moments vanish **exactly**. Parity is \((-1)^p\). Their supports lie in \([-173/500,173/500]\) because \(a+\varepsilon=173/500=0.346\). The closest selected source center is \((\log3)/2>0.549\), with the old bump radius \(1/20000\); thus every \(f_{p,j}\in F_s\).

The fourth vector is precisely a scalar multiple of T1D-U's old parity source \((\partial_x^2-1/4)(x^p(1-16x^2)^5*\rho_\varepsilon)\), now with the explicit common mollifier (T1D-U allowed *any* even mollifier inside radius \(1/200\)). The four reference polynomials are linearly independent: the first three are independent degree-\(3,4,5\) polynomials on \((-a,a)\); the last has strictly smaller support. Convolution by the nonzero compact bump preserves their linear independence, since its Fourier transform is nonzero near zero and the original compact distributions have entire Fourier transforms. Thus
\[
W_p=\operatorname{span}_{\mathbb C}\{f_{p,1},f_{p,2},f_{p,3},f_{p,4}\}\subset F_s,
\quad\dim W_0=\dim W_1=4.\tag{V.2}
\]
Together \(W_0\oplus W_1\) is an **eight-dimensional smooth source trial space** in the actual physical complement, not a spectral projection of \(D_F\).

### Exact arithmetic cancellation of the **entire** internal Weil matrix

The common support diameter is \(2(a+\varepsilon)=173/250=0.692\). The exact atanh-series audit establishes
\[
\boxed{\frac{173}{250}<\frac{693}{1000}<\log2.}\tag{V.3}
\]
Every genuine von Mangoldt shift has magnitude at least \(\log2\). Therefore **for all pairs of basis functions**, not merely their diagonals,
\[
\boxed{\langle Af_{p,i},f_{q,j}\rangle=0\quad\text{for all }p,q,i,j.}\tag{V.4}
\]
Parity already kills all mixed even–odd entries. The original pole contributions vanish by the exact moments. The complete original Weil matrix on \(W_0\oplus W_1\) is exactly the Gamma matrix. This is a structural property of the *original* arithmetic operator, not its deletion or modification. The prime contribution between \(E\) and \(W_p\) **need not vanish**, and is retained in the Schur cross vector below.

## 3. The exact, physically typed four-dimensional Ritz pencil

For \(p=0\) (even, \(\sigma=+\)) or \(p=1\) (odd, \(\sigma=-\)), define the genuine **smooth-source matrices**
\[
\begin{aligned}
M^{(p)}_{ij}&=\langle f_{p,i},f_{p,j}\rangle_{L^2(I)},\\
G^{(p)}_{ij}&=\gamma(f_{p,i},f_{p,j})
 =\int_{\mathbb R}\!m_\Gamma(\xi)\widehat f_{p,i}(\xi)\overline{\widehat f_{p,j}(\xi)}\,d\xi,\\
b^{(p)}_i&=q(u_\sigma,f_{p,i})
 =\gamma(u_\sigma,f_{p,i})-\langle Au_\sigma,f_{p,i}\rangle,\\
\mu_\sigma&=q(u_\sigma,u_\sigma)>7.
\end{aligned}\tag{V.5}
\]
The *full* original prime-power translation sum is present in \(b_i\), even though it vanishes on \(W_p\times W_p\). All integrals and pairings in (V.5) refer to the exact fixed mollified functions, not an auxiliary polynomial approximation. The true Birman–Schwinger Ritz pencil on \(D_F^{1/2}W_p\) is
\[
\boxed{D_p=6M^{(p)}+G^{(p)},\qquad
C_p=6M^{(p)}+\frac{b^{(p)}(b^{(p)})^*}{\mu_\sigma},\quad
\lambda_p^{(4)}=\lambda_{\max}(C_p,D_p).}\tag{V.6}
\]
(Conjugate orientation follows the original first-argument-linear convention; the rank-one matrix is unambiguous as the positive form \(|\sum_i t_i b_i|^2/\mu_\sigma\).) Its corresponding exact **complement Schur form** is
\[
\boxed{S_p=G^{(p)}-\frac{b^{(p)}(b^{(p)})^*}{\mu_\sigma}.}\tag{V.7}
\]
Thus \(\sup\sigma(\mathscr B_{F,\sigma})\ge\lambda_p^{(4)}\). When \(G^{(p)}>0\), the exact *finite* conditional Schur test becomes \((b^{(p)})^*(G^{(p)})^{-1}b^{(p)}\le\mu_\sigma\). **Neither the actual smooth \(G\) nor the actual \(b\) has been interval-enclosed here**, so no top four-dimensional eigenvalue or finite-block Schur-sign certification is claimed. The checker computes full exact **unmollified reference** \(L^2/H^1\) matrices; those are used as error-controlled inputs below, and are **not misidentified as** \(M^{(p)}\) or \(G^{(p)}\).

## 4. Stronger **certified** parity spectral lower bounds from rational three-profile combinations

Even though the actual 4×4 matrices remain to be enclosed, we can certify one genuine smooth vector in each 4D span:
\[
\boxed{\begin{aligned}
y_+&=f_{0,1}+\frac35f_{0,2}+\frac14f_{0,3},\\
y_-&=f_{1,1}+\frac7{10}f_{1,2}+\frac35f_{1,3}.
\end{aligned}}\tag{V.8}
\]
Let \(y_{0,\sigma}\) denote the corresponding **unmollified** \(H^1\) compact polynomial and \(J_\sigma^0=\|y_{0,\sigma}\|_2^2\), \(K_\sigma^0=\|y_{0,\sigma}'\|_2^2\). The exact symbolic audit evaluates the entire 4×4 rational mass and derivative Grams and in particular
\[
\begin{array}{c|cc}
 & J_\sigma^0 & K_\sigma^0\\\hline
+&2787056943244432048/995761081940625&1198159303630810162688/2528436539263635\\
-&48357210028498864/30174578240625&1130114705986457973248/2528436539263635
\end{array}\tag{V.9}
\]
All entries are exact rational integrals on their respective compact polynomial supports. Young's inequality and the probability-bump first-moment bound give \(\|y_\sigma-y_{0,\sigma}\|_2\le\varepsilon\|y'_{0,\sigma}\|_2\), \(\|y'_\sigma\|_2\le\|y'_{0,\sigma}\|_2\), and hence
\[
\frac{\|y_\sigma'\|_2}{\|y_\sigma\|_2}
\le\frac{r_\sigma}{1-\varepsilon r_\sigma},\qquad
r_\sigma^2=K_\sigma^0/J_\sigma^0.
\]
An **exact rational square comparison** (not a floating test) certifies
\[
\boxed{\frac{\|y_+'\|_2}{\|y_+\|_2}<\frac{66}{5}=13.2,
\qquad \frac{\|y_-'\|_2}{\|y_-\|_2}<17.}\tag{V.10}
\]
Import the **source-valid Jensen/digamma inequality** from [T1D-U](CR-INT-T1D-U.md), whose finite rational form at \(M=50\) is re-evaluated here:
\[
\frac{\gamma[y]}{\|y\|_2^2}\le m_\Gamma(0)+\sum_{n=0}^{50}
\frac{V}{(n+1/4)(4(n+1/4)^2+V)}+
\frac12\log\left(1+\frac{V}{4(50+1/4)^2}\right),
\quad V=\frac{\|y'\|_2^2}{\|y\|_2^2},\quad m_\Gamma(0)<-\frac{5372}{1000}.
\]
The exact-fraction checker proves the convenient strict source-level conclusions
\[
\boxed{\frac{\gamma[y_+]}{\|y_+\|_2^2}<\frac34,
\qquad\frac{\gamma[y_-]}{\|y_-\|_2^2}<1.}\tag{V.11}
\]
All internal primes vanish by (V.4), and \(\mathcal R\ge0\) is retained by its sign. Thus the actual compact Birman–Schwinger Rayleigh quotients on \(D_F^{1/2}y_\pm\) give the **new rigorous one-sided** parity spectral bounds
\[
\boxed{\sup\sigma(\mathscr B_{F,+})>\frac{6}{6+3/4}=\frac89\approx0.88889,\quad
\sup\sigma(\mathscr B_{F,-})>\frac{6}{6+1}=\frac67\approx0.85714.}\tag{V.12}
\]
These strengthen the T1D-U bounds \(>60/73\approx0.82192\), \(>300/377\approx0.79576\), without importing T1D-I's norm bound or assuming an unproved spectral upper enclosure. They do **not** show the leading eigenvalue is below or above 1, so they do not settle the Schur inequality.

## 5. Explicitly **noncertifying** four-dimensional Fourier/Gamma-only Ritz diagnostic

For the **explicit** probability mollifier in §2, an optional Fourier-grid diagnostic constructs the four actual smooth source vectors for each parity and evaluates the generalized pencil
\[
\widetilde C_p=6M^{(p)},\qquad \widetilde D_p=6M^{(p)}+G^{(p)}.
\]
This **omits** \(bb^*/\mu\) *only in the diagnostic*, where it supplies a conservative matrix-relative comparator (the omitted piece is nonnegative); it is **not** substituted for the full source matrix in (V.6). The following approximate results were numerically stable over three periodic Fourier grids but are **not enclosures**:

| Parity | Exploratory eigenvalues of \(\widetilde C_p\) relative to \(\widetilde D_p\), increasing order | Approx. largest |
| --- | --- | ---: |
| even | 0.75553, 0.79372, 0.83580, 0.91420 | **0.91420** |
| odd | 0.74559, 0.78022, 0.81520, 0.87492 | **0.87492** |

The true source-defined Ritz tops may be larger because of the positive rank-one \(\mathcal R\) contribution. **No interval, no true \(\lambda_p^{(4)}\), and no infinite-dimensional spectral upper bound follow from the displayed decimals**. The numerical checker is deliberately separate from the rational regression runner and requires NumPy/SciPy.

## 6. Verification contract, stopping point, and next step

**Analytically derived and exact-rational audited:** eight smooth source vectors in the actual pole-null complement; nondegenerate four-source parity subspaces; all \(4\times4\) internal prime-power entries zero on source support; full exact symbolic original-form Ritz identity retaining Gamma, pole-null and selected-space cross interactions; rational reference polynomial mass and derivative matrices; two strictly improved actual compact-spectrum **lower** witnesses (V.12).

**Only analytically represented, not numerically certified:** the mollified physical mass matrix \(M^{(p)}\), Gamma matrix \(G^{(p)}\), original selected cross vector \(b^{(p)}\), and actual generalized Ritz eigenvalues. The optional numerical Fourier values are exploratory and omit the rank-one term by construction. The exact audit does not certify the digamma/Weil functional-analysis proof, Suzuki's compactness, or any RH conclusion.

**Still OPEN:** rigorous interval evaluation of \(G^{(p)}\) and \(b^{(p)}\), a certified four-dimensional Schur matrix/sign, actual low/high Gamma-energy \(D_F\)-spectral projection and \(\eta_N\), full complement positivity, Suzuki–Connes completed-domain identification, original pole directions, global Weil positivity, RH.

**Next — CR-INT-T1D-W:** use the finite polynomial structure and compact smooth mollifier to bound each entry of \(G^{(p)}\) and \(b^{(p)}\) with outward-rounded **analytic quadrature/interval arithmetic** (or a finite-dimensional form bound strong enough to avoid explicit quadrature). The immediate falsifiable gate is the sign of \(S_p=G^{(p)}-bb^*/\mu\) in each parity. Do not advertise a certified sign from numerical FFT. Even an eventually certified positive 4D matrix requires a separate infinite-dimensional tail bound before any Phase IV transfer.

**Reproduction:** `python3 src/audits/canonical_return_t1d_v_wide_ritz.py` (SymPy) for the exact rational/symbolic results. Optional `python3 tools/diagnostic_t1d_v_wide_ritz.py` (NumPy/SciPy) for the noncertifying Fourier comparator. Research claim remains **DRAFT** under [governance/CLAIMS.md](../governance/CLAIMS.md).