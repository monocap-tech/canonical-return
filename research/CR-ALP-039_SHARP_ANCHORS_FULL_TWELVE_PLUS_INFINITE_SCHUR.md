# CR-ALP-039 — The original signed Weil twelve-mode plus infinite-tail Schur certificate

**Research status:** 2026-11-23 PST. **DRAFT / NONPROMOTING**. Classical analytical continuum argument with two independently executed outward-integer arithmetic audits. **Not** a Lean theorem, not a positivity theorem on all projected sine modes through 15,000, not arbitrary-support positivity and **not RH**.

## 0. Exact source and predecessor state

The original pinned definition is `Zeta23.EF.literatureRHS` in `monocap-tech/weil-lab`, branch `research/iv068-original-source-kernel-gate`, path `WeilDefect/External/Zeta23/Definitions/Def_Zeta23_ExplicitFormula.lean`, blob `0034976fe2f1e5fa9d411108736160542eb94293`. At physical radius `a=11/10`, it carries both original Mellin pole evaluations, the **negative** paired authentic von Mangoldt prime terms at `n=2,3,4,5,7,8,9`, and the **complete** Gamma multiplier `Re psi(1/4+ir/2)-log(pi)` in the original paper Fourier convention. No sign or test-space substitution is made here.

Let `I=(-a,a)`, `H0` be the original two-pole-null subspace of `L²(I)`, and `e_j(x)=a^-1/2 sin(j*pi*(x+a)/(2a))`, `nu_j=j*pi/(2a)`. With `kappa(j)=1` for odd `j` and `2` for even `j`, let

\[
 w_j=e_j-\lambda_j e_{\kappa(j)},\quad
 \lambda_j=\frac{\nu_j(\nu_{\kappa(j)}^2+1/4)}{\nu_{\kappa(j)}(\nu_j^2+1/4)},
 \qquad u_1=\Pi_0e_1,\quad u_2=\Pi_0e_2.
\]

The physical finite block is `E12=span{w3,...,w12,u1,u2}`. The *actual* orthogonal complement beyond 15,000 projected sine modes is

\[
 Q_N=R_N-R_N B(B^*R_NB)^{-1}B^*R_N,\quad N=15000,\quad
 B=(\cosh(x/2),\sinh(x/2)),\quad R_N=I-S_N.
\]

CR-ALP-027 analytically proves `W(q,q)>=1/4 ||q||²` on `Q15000 H0` in the log-form domain; CR-ALP-034 gives the original signed twelve-dimensional lower physical matrix `Hmid` and `W|E12>=Hmid`. CR-ALP-038 certifies the entire *interior* ten-mode mixed sector by a rank-one-plus-diagonal upper residual matrix, but **explicitly excludes** both boundary anchors. The predecessor source matrix files are pinned by SHA256:

- `H_mid.tsv`: `b6a31e623298272fdceac788fc56626ec704e280903ff6baff1dc8bfe0d45f92`
- `initial.txt` (original signed prime and physical metric): `bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2`

The CR-ALP-034 and CR-ALP-038 analytical premises are **not native Lean results**. Source-indexed exact arithmetic is a separate status from compiler certification.

## 1. Exact high-sine Gamma coefficient for an original pole function

Let `b` be either genuine pole function, `b0=b(-a)` (which is nonzero in both cases), `c_j=<e_j,b>` in the matching reflection parity, and `beta_m=2m+1/2`. The original physical Gamma kernel admits the genuine digamma-resolvent representation

\[
 T_\Gamma=g(0)I+\sum_{m=0}^{\infty}J_{\beta_m},\qquad
 J_\beta(f,f)=\int_0^\infty e^{-\beta s}\|f-\tau_s f\|_2^2 ds.
\]

The exact true pole Fourier-sine coefficient is

\[
 c_j=\frac{2b_0\nu_j}{\sqrt a(\nu_j^2+1/4)}.
\]

Define `I_b(beta)=int_0^(2a) exp(-beta t)b(-a+t)dt` and `E_m=I_b(beta_m)-b0/beta_m`. The CR-ALP-037 physical Laplace-kernel identity implies, for every matching high `j`:

\[
\begin{aligned}
\langle e_j,T_\Gamma b\rangle
&=\Lambda(\nu_j)\left(c_j-\frac{b_0}{\sqrt a\,\nu_j}\right)
 -\frac{2\nu_j}{\sqrt a}\sum_{m\ge0}\frac{E_m}{\nu_j^2+\beta_m^2},\\
\Lambda(\nu)&=\Re\psi(1/4+i\nu/2)-\psi(1/4).
\end{aligned}\tag{39.1}
\]

This identity follows **without divergent renormalization**: each sum in (39.1) is convergent. The digamma series identity is exactly

\[
\Lambda(\nu)=\sum_{m\ge0}\frac{2\nu^2}{\beta_m(\beta_m^2+\nu^2)}.
\]

Let the absolutely convergent `E_infty=sum_(m>=0) E_m` and choose the source-dependent scalar

\[
 t_N(b)=\tfrac12\bigl(\log(\nu_N/2)-\psi(1/4)\bigr)-\frac{E_\infty}{b_0}.
\]

Subtracting this *exact multiple of the original pole function* is legal because `Q_N b=0`. Use the exact algebra

\[
 c_j-\frac{b_0}{\sqrt a\nu_j}=c_j\bigl(\tfrac12-\tfrac{1}{8\nu_j^2}\bigr),
\qquad
\frac1{\nu^2+\beta^2}=\frac1{\nu^2+1/4}
-\frac{\beta^2-1/4}{(\nu^2+1/4)(\nu^2+\beta^2)}.
\]

Writing `delta(nu)=Re psi(1/4+i nu/2)-log(nu/2)` gives the exact **centered-log decomposition**:

\[
\boxed{
\begin{aligned}
\langle e_j,T_\Gamma b-t_N(b)b\rangle
={}&\frac{c_j}{2}\log\frac{j}{N}
 +\frac{c_j}{2}\delta(\nu_j)
 -\frac{c_j\Lambda(\nu_j)}{8\nu_j^2}\\
&+\frac{2\nu_j}{\sqrt a(\nu_j^2+1/4)}
\sum_{m\ge0}E_m\frac{\beta_m^2-1/4}{\nu_j^2+\beta_m^2}.
\end{aligned}}
\tag{39.2}
\]

**Key point:** a naive bound would retain `(log j)/j`. The permitted `Q_N` pole subtraction cancels the constant `(log N)/j` contribution exactly. The surviving leading contribution is **`log(j/N)/j`**, whose full infinite square sum is `O(1/N)` with a small explicit constant.

## 2. Uniform rational Gamma bound for both exceptional boundary anchors

The elementary series-versus-integral total-variation proof in CR-ALP-026 gives

\[
 |\delta(\nu)|\le 2/\nu+1/(8\nu^2)<3/\nu\quad(\nu\ge1).
\]

Because `nu_N>20000`, the exact special value `psi(1/4)=-gamma_E-pi/2-3 log 2`, `gamma_E<1`, `log2<3/4`, `pi<22/7`, and `log 20000<11`, give the conservative bound `0<=Lambda(nu)<nu/1000` on every `nu>=nu_N`. The audit verifies the required exponential comparisons via finite rational Taylor sums, including `e^11>20000`.

For `B0>=||b||_infty`, `B1>=||b'||_infty`, integration by parts gives

\[
 |E_m|\le B_1/\beta_m^2+B_0 e^{-2a\beta_m}/\beta_m.
\]

The `m=0` summand in the last term of (39.2) vanishes identically. For `m>=1`, the genuine series obeys

\[
 \sum_{m\ge1}\frac1{\nu^2+\beta_m^2}\le \frac{\pi}{4\nu},\qquad
 \sum_{m\ge1}\beta_m e^{-2a\beta_m}<\frac1{70}.
\]

The second inequality is audited using `e^(11/2)>200` and `e^(22/5)>80`. Substitution into (39.2), with `3<pi<22/7`, `sqrt(a)<21/20`, gives the **pointwise full-series error bound**

\[
\boxed{
|\langle e_j,T_\Gamma b-t_N(b)b\rangle|
\le \frac{7B_0}{10}\frac{\log(j/N)}{j}
 +\frac{C_2(b)}{j^2}+\frac{C_3(b)}{j^3},\quad j>N,}
\]

where valid rational constants are

\[
 C_2(b)=\frac{121}{225}\left(\frac{12001}{4000}B_0+\frac{11}{7}B_1\right),
 \qquad C_3(b)=\frac{8a^3}{945}B_0.
\tag{39.3}
\]

A unimodal integral comparison gives `sum_(j>N) log²(j/N)/j² <= 2/N+2/N² <1/86²`, and the standard reciprocal-power tails are `<1/(3N³)<1/(3,000,000)²` and `<1/(5N⁵)<1/(10^10)²`. Dropping the parity restriction is conservative. Consequently

\[
 \|Q_NT_\Gamma b\|_2
 <\frac{7B_0}{860}+\frac{C_2(b)}{3,000,000}+\frac{C_3(b)}{10^{10}}.
\]

The original anchors `u_k=e_k-gamma_k b_k` satisfy `|gamma_1|<3/5`, `|gamma_2|<5/3`; CR-ALP-037's independently audited `||Q_NT_Gamma e_k||<10^-6` applies to both `k=1,2`. Set `(B0,B1)=(6/5,3/10)` for cosh and `(3/5,3/5)` for sinh. The **complete** original Gamma anchor column bounds are certified arithmetically as

\[
\boxed{\|Q_{15000}T_\Gamma u_1\|_2<1/170,\qquad
       \|Q_{15000}T_\Gamma u_2\|_2<1/120.}
\tag{39.4}
\]

They improve the CR-ALP-037 values `1/20` and `1/9` without truncating the Gamma series.

## 3. Source-signed prime translation: explicit oscillatory endpoint term

Let `U_{jb}(t)=<e_j,(tau_t+tau_-t)b>`, where `b` is genuinely zero-extended and then restricted to `I`, and `0<t<2a`. Pairing the actual two shifts and integrating by parts twice yields the *exact* first-order boundary expansion

\[
\boxed{
 U_{jb}(t)=\frac{2}{\sqrt a\nu_j}
\left[b(-a+t)+b(-a)\cos(\nu_jt)\right]+\mathcal R_{jb}(t),}
\]

with the uniform rigorous remainder

\[
 |\mathcal R_{jb}(t)|\le\frac{2B_1+aB_0}{\sqrt a\nu_j^2}.
\tag{39.5}
\]

The term `b(-a+t)/nu_j` is **nonoscillatory**. It can be cancelled exactly (up to an elementary `O(nu_j^-3)` resolvent mismatch) by a scalar multiple of the original pole coefficient `c_j`. The term `b(-a)cos(nu_j t)/nu_j` is **not** removed: it is the true internal jump of the translated source and must remain in the bound.

For the true ordinary sine anchor `e_k`, the exact source overlap is

\[
 U_{jk}(t)=\frac2a\frac{\nu_j\sin(\nu_kt)-\nu_k\sin(\nu_jt)}{\nu_j^2-\nu_k^2}
\]

in matching parity, and its leading `2 sin(nu_k t)/(a nu_j)` is likewise pole-aligned. After these legal subtractions, use `nu_j>=2 nu_k` and all **seven original negative weights** `omega_n=Lambda(n)/sqrt(n)`, bounded only for errors by `S_a=sum omega_n<177/50`. One obtains for each anchor `u_k=e_k-gamma_k b_k`, uniformly across every genuine high sine index,

\[
 \boxed{\left|\langle e_j,Pu_k-t_k^Pb_k\rangle\right|
 \le\frac{A_k^P}{j}+\frac{B_k^P}{j^2}+\frac{C_k^P}{j^3},}\tag{39.6}
\]

with explicit rational upper constants

\[
\begin{aligned}
 A_k^P&=\frac75|\gamma_k|S_a B_0,\\
 B_k^P&=\left(\frac{8S_a\overline\nu_k}{3a}
       +|\gamma_k|S_a(2B_1+aB_0)\right)(2a/3)^2,\\
 C_k^P&=\left(\frac{8S_a(\overline\nu_k^2+1/4)}{3a}
       +\frac12|\gamma_k|S_aB_0\right)(2a/3)^3,\qquad
 \overline\nu_k=k(22/7)/(2a).
\end{aligned}
\]

The infinite sine square-tail bound uses `sum j^-2<1/N<1/122²`, `sum j^-4<1/(3N³)`, `sum j^-6<1/(5N⁵)`. It proves the **full seven-original-prime** anchor estimates

\[
\boxed{\|Q_{15000}Pu_1\|_2<1/30,\qquad
       \|Q_{15000}Pu_2\|_2<1/23.}
\tag{39.7}
\]

No sign reversal, prime-dropping, positivity assumption, or finite-only frequency check is involved. A separate numerical integration regression of (39.5) at sine indices 13–16 and shift `log 2` agrees within the proven analytic remainder; this floating-point check is diagnostic **only**.

Combining (39.4), (39.7) by the triangle inequality gives the new source-signed **entire-infinite-tail** estimates

\[
\boxed{\|Q_{15000}Tu_1\|_2<1/25,\qquad
       \|Q_{15000}Tu_2\|_2<1/19.}
\tag{39.8}
\]

## 4. NEW complete E12-plus-infinite-tail mixed Schur theorem

From the existing CR-ALP-038 exact source-specific interior coefficient calculation, each parity's five interior source columns satisfy the rank-one-plus-diagonal bound

\[
 \|Q_NT(\sum_k x_kw_k)\|^2
 \le2\alpha^2(D^Tx)^2+2\mathfrak b^2\|x\|_{\ell^2}^2.
\tag{39.9}
\]

Here `D_k=nu_k-lambda_k nu_kappa`, and `alpha`, `b²` are the **same outward-rational values already proved** in CR-ALP-038; the original physical matrix is not altered. With an anchor coefficient `y`, take `eta=2` in the squared triangle inequality:

\[
 \|Q_NT(\sum x_kw_k+yu_\kappa)\|^2
 \le3\|Q_NT\sum x_kw_k\|^2+\tfrac32\epsilon_\kappa^2|y|^2,
\]

where `eps_1=1/25`, `eps_2=1/19`. Since `W(q,q)>=1/4||q||²` for any actual log-form-domain `q` in the original `Q_NH0`, and

\[
 2\|Q_NTp\|\|q\|\le8\|Q_NTp\|^2+\tfrac18\|q\|^2,
\]

the sufficient finite matrix is **exactly**

\[
\boxed{
 \widetilde H_{34}-48\alpha^2 DD^T-48\mathfrak b^2 I_{
m int}
 -12\epsilon_\kappa^2 e_{\rm anchor}e_{\rm anchor}^T
 -10^{-12}M_6\succ0
}
\tag{39.10}
\]

in each of the two six-dimensional reflection-parity sectors. `M6` is the **true physical** Gram, including all anchor cross-entries; it is never replaced by an identity matrix. For `E12` the full matrix has 72 ordered source entries.

### Exact arithmetic certification

The primary stdlib `Fraction` + `10^-25` outward-integer program reads the CR-ALP-034 lower matrix and physical Gram from their pinned files, independently derives the `39.4`–`39.8` scalar inequalities, constructs the **full** shifted 12×12 parity-block matrix, and performs outward interval `LDL^T`. The comparison is **strictly positive on all twelve pivots**. The independent program separately rounds the saved 72 interval entries to a `10^-20` grid, implements its own outward interval arithmetic and `LDL^T`, and rechecks the anchor inequalities. Both executed **PASS**.

The primary shifted interval bank has SHA256:

`a631f3cc7b99e47e4749a7931cd1307af336b5e20fe28b28946aa819ba56af77`

The exact last lower pivots of (39.10), at scale `10^-25`, are:

- odd sine indices / even functions, `[3,5,7,9,11,1]`: `71664015340934032282668/10^25 > 0`;
- even sine indices / odd functions, `[4,6,8,10,12,2]`: `49884032376053045649668/10^25 > 0`.

**The complete strict quantitative theorem, conditional only on the already identified classical continuum premises, is**

\[
\boxed{
 W_{\rm original}(p+q,p+q)>
 10^{-12}\|p\|_{L^2(I)}^2+\frac18\|q\|_{L^2(I)}^2,
 \quad 0\ne p+q,\ p\in E_{12},\ q\in Q_{15000}\mathscr H_0\cap\mathcal D_{\log}.
}
\tag{39.11}
\]

The form's complex polarization follows by Hermitian extension; the classical argument is real coefficientwise and separates the parity blocks. `p` is orthogonal to `q` in the genuine physical metric, so the displayed terms are comparable and meaningful. The pole-null test sources `u1,u2` lie in the operator domain despite their zero-extension endpoint jumps, because their Fourier transforms are `O(1/r)` and multiplication by `log(1+|r|)` remains in `L²`.

## 5. Precisely what remains OPEN

1. **Intermediate original signed finite shell:** none of the additional 14,988 projected physical directions between `E12` and `E15000`, nor their mixed interactions, is certified nonnegative by this result. (The numerical count is `15000-12=14988`.)
2. **Entire fixed-collar two-pole-null positivity:** still open because of that intervening shell.
3. **Original full unrestricted two-pole case, collar-radius extension, Suzuki–Connes transfer and original RH:** all open.
4. **Native Lean 4 certification:** not attempted; exact rational scripts are standalone computational audits, and continuum Fourier/Gamma facts have not entered the Lean kernel.
5. **Acceptance-ledger governance:** only draft Canonical-Return research artifacts are eligible. No Weil-Lab, formalization receiving gate, F41/F42, provenance, Hodge or RH promotion is authorized.

**Next CR-ALP-040:** extend the signed mixed Schur certificate to intermediate physical shells, preferably using the now-certifiable rank-one-plus-diagonal Gamma/prime and two-boundary updates. A positive `E12⊕Q15000` subspace is not a positivity theorem on its missing finite middle shell.