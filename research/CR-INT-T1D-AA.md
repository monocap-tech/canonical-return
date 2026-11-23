# CR-INT-T1D-AA — Fourier-localized mixed Weil representers close the *conditioned* second Schur gate at B = 2²¹

**Canonical Return · November 23, 2026 · independent smooth-source analytic research.**  
**Disposition:** **DERIVED / ANALYTIC** with exact-rational executable arithmetic audit and explicit imported source lemmas. The new result is an **infinite-dimensional, prolate-conditioned second-Schur positivity certificate without the four extra mixed-functional annihilation constraints**, at the comparatively small cutoff `B=2^21`. It is **NOT** the full low-mode complement, completed Suzuki–Connes correspondence, a Lean certification, a Phase IV transfer, or RH. The prolate concentration threshold is tightened to `theta=2^-38`, so this is **not** positivity on the larger subspace used in T1D-Q at `theta=1/128`.

## 1. Freeze the original source and the two successive Schur reductions

Preserve Suzuki's **complete original signed Weil functional** ([arXiv:2606.09096v4, §1.1, §2.4](https://arxiv.org/html/2606.09096v4)), physical `H=L²(I,dx)`, `I=(-2,2)`, and genuine smooth two-Mellin-pole-null core

\[
\mathcal D^0=\{f\in C_c^\infty(I):\int_I f(x)e^{\pm x/2}dx=0\}.
\]

Let the selected T1C parity pair be `u_+`,`u_-` (physical orthonormal, selected Weil values `mu_±>7`); let `P` be the span of the physical two-pole carriers. On `F_s=\mathcal D^0\cap (E\oplus P)^\perp`, the first exact source Schur remainder is

\[
\mathfrak s(v,w)=q(v,w)-\langle\mathcal Rv,w\rangle,
\qquad \mathcal R[y]=\sum_{\sigma=\pm}\frac{|q(u_\sigma,y)|^2}{\mu_\sigma}.
\tag{AA.1}
\]

Take T1D-Y's **original fixed smooth** central–satellite physical unit vectors `v_1,...,v_4` (two even and two odd), spanning `L\subset F_s`. The completed original Weil Schur matrix `S_L=(\mathfrak s(v_i,v_j))` is strictly positive, with the previously derived source bounds

\[
\boxed{\lambda_{\min}(S_+)>1/70,\qquad \lambda_{\min}(S_-)>1/18.}\tag{AA.2}
\]

Put `Z_s=F_s\cap L^\perp`. The second Schur penalty on `Z_s` is the original-source **rank-at-most-four** operator

\[
\mathcal T[w]=b(w)^*S_L^{-1}b(w),\qquad b_i(w)=\mathfrak s(v_i,w).
\tag{AA.3}
\]

The still-open unconditioned complement condition is `q[w]-\mathcal R[w]-\mathcal T[w]>=0` for **every** `w\in Z_s`. T1D-Z bounded `||T||` by a very loose number and could dominate it without extra conditions only at `B=2^1,570,000`. This note obtains a source-adapted bound **on a conditioned test sector** instead.

## 2. Whole-line Fourier representation: avoid the artificial window boundary

For compact smooth pole-null `v,w` supported in `I`, their **complete original** signed Weil mixed pairing is the Fourier pairing of the zero extensions, with symbol

\[
\boxed{\mathfrak k(\xi)=m_\Gamma(\xi)-a(\xi),\quad
m_\Gamma(\xi)=\Re\psi(1/4+i\xi/2)-\log\pi,\quad
 a(\xi)=2\sum_{n=p^j\le54}\frac{\Lambda(n)}{\sqrt n}\cos(\xi\log n).}\tag{AA.4}
\]

There are **24 true prime powers**. Prime powers larger than 54 cannot couple two tests in `I` because their log shifts exceed the physical diameter 4; this is an exact support elimination, not a modification of the Weil source. Let `\mathcal K_\infty v` be the *whole-line* Fourier-multiplier image `\mathfrak k\widehat v`. Then

\[
\boxed{q(v,w)=\langle\mathcal K_\infty v,w\rangle_{L^2(\mathbb R)}.}\tag{AA.5}
\]

Unlike the **physically projected** representer, `\mathcal K_\infty v` is not artificially cut off at `x=\pm2`. This distinction is valuable: the Fourier tails can be controlled from the derivative of the original smooth source, without introducing jump discontinuities through a zero-extended physical projection. We do **not** identify `\mathcal K_\infty` with Suzuki's inverse Neumann Laplacian.

**Imported from T1D-F (analytic, not reproved by the AA checker):** `|m_Gamma(xi)| <= log(1+|xi|)+40`. The checker independently encloses every original von Mangoldt coefficient and verifies

\[
\boxed{\sup_\xi|a(\xi)|\le2\sum_{n\le54}\Lambda(n)/\sqrt n<25.}\tag{AA.6}
\]

Thus `|\mathfrak k(\xi)|<\log(1+|\xi|)+65`. Since `\log(1+x)\le x` and `(\log(1+x)+65)/x` decreases for `x>0`, for any `v\in C_c^\infty(I)` and `B>0` with `\log(1+B)<15` we have

\[
\boxed{\begin{aligned}
\|\mathcal K_\infty v\|_2&\le 65\|v\|_2+\|v'\|_2,\\
\|\mathbf1_{|\xi|>B}\widehat{\mathcal K_\infty v}\|_2&\le\frac{80}{B}\|v'\|_2.
\end{aligned}}\tag{AA.7}
\]

These are **operator estimates for fixed smooth source images**, not an assertion that the full Weil quadratic form is bounded on `L²×L²`. The checker proves `e^15>1+2^21` by an exact rational positive Taylor partial sum.

## 3. Frequency conditioning bounds the true mixed functionals

Fix

\[
\boxed{B=2^{21},\qquad\theta=2^{-38}.}\tag{AA.8}
\]

As in T1D-M/Q, let `\mathsf T_B=P_I\mathcal F^{-1}\mathbf1_{[-B,B]}\mathcal F P_I` be the physical prolate concentration operator and impose orthogonality to its finite set of eigenvectors with eigenvalue `>theta`. Let `\mathcal H_{B,\theta}\subset Z_s` be the resulting **smooth** subspace. Then

\[
\int_{-B}^B|\widehat w(\xi)|^2d\xi\le\theta\|w\|_2^2
\qquad(w\in\mathcal H_{B,\theta})
\tag{AA.9}
\]

with unitary Fourier normalization. Its intersection with `Z_s` remains infinite-dimensional (finitely many continuous constraints on `C_c^\infty(I)`). There are **no further four mixed-functional null equations**.

Split the Fourier pairing in (AA.5) at `B`. From (AA.7), if `||v||_2=1` and `||v'||_2<r`,

\[
\boxed{|q(v,w)|<\left[(65+r)\,2^{-19}+\frac{80r}{2^{21}}\right]\|w\|_2.}\tag{AA.10}
\]

For the two selected functions `u_sigma`, T1C's exact bound `ell²||phi'||²/||phi||² <=30625/2304` with `ell=10^-4` and disjoint normalized cells gives `||u_sigma'||<37000`. For each of the four **actual** normalized central–satellite sources `v_i`, T1D-V/W/X establishes `||v_i'||<46`; no profile is altered here. The script checks the **strict exact rational inequalities**

\[
\boxed{\begin{aligned}
|q(u_\sigma,w)|/\|w\|_2
&<\frac{777065}{524288}<\frac32,\\
|q(v_i,w)|/\|w\|_2
&<\frac{1031}{524288}<\frac1{500}.
\end{aligned}}\tag{AA.11}
\]

The bounds apply to every `w` in the conditioned smooth sector and use **whole-line** arithmetic-plus-Gamma images; they are not obtained by replacing the complete Weil form with a positive prime operator.

## 4. Both finite-rank Schur penalties are dominated without T1D-I

For either parity, exactly one selected direction can pair nontrivially. By (AA.11) and `mu_sigma>7`,

\[
\boxed{\mathcal R[w]<\frac{9}{28}\|w\|_2^2.}\tag{AA.12}
\]

T1D-W/Y supplies `|q(v_i,u_sigma)|<27/100` for all four normalized source directions (with stricter bounds for the other channels). Since

\[
b_i(w)=q(v_i,w)-\sum_\sigma\frac{q(v_i,u_\sigma)q(u_\sigma,w)}{\mu_\sigma},
\]

parity and (AA.11) imply

\[
\boxed{|b_i(w)|<\left(\frac1{500}+\frac{27}{700}\frac32\right)\|w_\sigma\|_2
<\frac3{50}\|w_\sigma\|_2.}\tag{AA.13}
\]

Use (AA.2), the two coordinates in each parity and `\|S_\sigma^{-1}\|=1/\lambda_{\min}(S_\sigma)` to obtain

\[
\boxed{\mathcal T[w]<70\cdot2\left(\frac3{50}\right)^2\|w_+\|_2^2
+18\cdot2\left(\frac3{50}\right)^2\|w_-\|_2^2
\le\frac{63}{125}\|w\|_2^2.}\tag{AA.14}
\]

These are **new, frequency-localized bounds on the actual selected and second-Schur penalties**, not smaller global operator norms. No T1D-I residual-norm estimate is used.

## 5. Main conditioned second-Schur theorem at the original T1D-Q cutoff

Import, without claiming independent formalization, T1D-M/O/Q's original complete Gamma lower estimate

\[
m_\Gamma(\xi)>\log|\xi|-49/10\quad(|\xi|\ge64),\qquad
m_\Gamma(\xi)>-6\quad(\xi\in\mathbb R),
\]

and T1D-Q's **joint** compressed arithmetic norm `||A||<203/25`. By (AA.9), for nonzero `w\in\mathcal H_{B,\theta}`,

\[
q[w]>\left[(1-\theta)\left(21\log2-\frac{49}{10}\right)-6\theta-\frac{203}{25}\right]\|w\|_2^2
>\boxed{\frac{153}{100}\|w\|_2^2}.\tag{AA.15}
\]

The last inequality is checked exactly using `\log2>693/1000` and `theta=2^-38`. Combining (AA.12), (AA.14), and (AA.15) yields the **new infinite-dimensional conditioned second Schur certificate**

\[
\boxed{\forall\,0\ne w\in\mathcal H_{2^{21},2^{-38}}:\quad
q[w]-\mathcal R[w]-\mathcal T[w]
>\frac{1233}{1750}\|w\|_2^2
>\frac7{10}\|w\|_2^2.}\tag{AA.16}
\]

**Improvement over T1D-Z:** the no-additional-mixed-null version of the *same* second Schur test previously required the conservative cutoff `B=2^1,570,000`. The present bound uses `B=2^21` and no added mixed-null constraints, by exploiting **Fourier localization of fixed smooth Weil representers**. The trade-off is **much tighter** low-band leakage `theta=2^-38` rather than T1D-Q's `1/128`; this creates more excluded prolate vectors. No new precise excluded-rank bracket has been computed for the tighter threshold. This is a conditioned high-frequency theorem, **not** unconditioned positivity of `Z_s`.

## 6. Audit and non-transfer contract

**Reproduction:** `python3 src/audits/canonical_return_t1d_aa_fourier_localized_second_schur.py` (Python standard library). The checker recomputes all 24 prime powers, outward rational `log/sqrt` weights and `2S<25`, the rational exponential inequality `log(1+B)<15`, exact selected/profile `H¹` budgets, both rational Fourier-split constants, the exact parity-specific Schur arithmetic, the new `B=2^21` Gamma/prime margin and the final `1233/1750` lower bound. The **analytic arguments** (whole-line Fourier formula, Gamma symbol bounds, Plancherel/Cauchy–Schwarz, prolate compactness/concentration, finite Schur identity and form-core transfer) are explicitly imported or proved in the note, **not mechanically verified** by the script. Re-run T1D-Q, W, X, Y, Z for the independently stated source input estimates.

**OPEN / NO TRANSFER:** the excluded low-frequency subspace (and its mixed original Weil matrix), unconditioned positivity on all of `Z_s`, a certified complete top eigenvalue `sup sigma(B_F)<=1`, control of the physical pole carriers beyond the two-pole-null core, Suzuki–Connes semilocal arithmetic/metric/domain transport, support-independent global Weil positivity, Phase IV theorem promotion, and RH. A smaller high-frequency conditioning threshold is **not** a proof of any of those obligations.

**Next — CR-INT-T1D-AB:** exploit the exceptionally tight `theta=2^-38` gate to derive an honest quantitative **excluded prolate rank** from trace/trace-defect at `B=2^21`, and investigate source-aware low-frequency spectral or form estimates on that omitted block. Independently review the T1D-F whole-line Gamma upper multiplier and the W/X/Y mixed source normalizations before promoting analytical claims. The decisive gap remains the *unconditioned* second Schur inequality, not a further reduction of `B` alone.
