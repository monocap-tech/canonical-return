# CR-INT-T1D-AC — Weighted finite-Fourier transfer for the actual low-frequency signed Weil/second-Schur form

**Canonical Return · November 23, 2026 (America/Los_Angeles) · independent analytical research candidate.**

**Status:** new **DERIVED / ANALYTIC, conditional on the explicitly imported T1D-AB/F/AA/Y/Z lemmas**, with an exact-fraction audit of the constants. **Not Lean-certified, not independently reviewed, not a low-mode positivity certificate, not a finished source/Connes transfer, not Phase IV and not RH.** Unlike T1D-AB's original rank optimization, this note deliberately enlarges the finite Fourier window to control the *whole logarithmic form energy*.

## 1. Original source and the correct constrained low sector

Let `I=(-2,2)`, `H=L²(I,dx)`, and let `D⁰=C_c^∞(I)∩{<f,e^{x/2}>=<f,e^{-x/2}>=0}`. The source-original signed pairing (Suzuki, arXiv:2606.09096v4 §§1.1,2.4) is, **on those pole-null tests**, the full Fourier form

\[
q_0(f,g)=\int_{\mathbb R} \left(m_\Gamma(\xi)-2\sum_{n=p^j\le54}\frac{\Lambda(n)}{\sqrt n}\cos(\xi\log n)\right)\widehat f(\xi)\overline{\widehat g(\xi)}\,d\xi,
\quad m_\Gamma(\xi)=\Re\psi(\tfrac14+i\xi/2)-\log\pi.
\tag{AC.1}
\]

There are **24** actual von Mangoldt prime powers in the fixed source window; the two original pole terms vanish only under the displayed smooth constraints. The source-specific whole-line Fourier identity and `|m_Gamma(xi)|<=log(1+|xi|)+40` are analytical imports from T1D-F/AA, not verified by this audit. The arithmetic symbol is bounded in absolute value by `<25` (T1D-AA), so `|k(xi)|<=log(1+|xi|)+65`.

Let `C_8` be the physical span of **two pole carriers**, the physically orthonormal selected `u_+,u_-`, and the four physically orthonormal, central/satellite `v_i` from T1D-Y. These eight constraint carriers are mutually orthogonal after normalizing the even/odd pole basis. Write `J=P_{C_8^⊥}` and `Z=C_8^⊥`, and `Z_s=D⁰∩Z`. The authentic second Schur form on `Z_s` is `h=q_0-R-T`, where `R` and `T` are precisely the selected rank-2 and central/satellite rank-4 **positive** penalties. Extend `R,T` by zero on orthogonal complements when writing the form on the larger physical H; their imported crude global norms are `||R||<1785`, `||T||<1,073,590`.

At the fixed `B=2^21`, `theta=2^-38`, put `T_B=S_B S_B^*`, `E=Ran 1_(theta,1](T_B)` and `H_hi=Z∩E^⊥`. T1D-AA's high-pass estimate gives `h[w]>(1233/1750)||w||²` on the smooth core in `H_hi`, conditional on its analytical imports. The **actual finite low space of Z** is

\[
\boxed{L_{\rm low}=Z\ominus(Z\cap E^\perp)=\operatorname{Ran}(J|_E),\qquad \dim L_{\rm low}\le\dim E\le 2{,}693{,}201.}\tag{AC.2}
\]

It is *not* `E∩Z` in general. `J|E` can have small nonzero singular values; the previous rank theorem gives no lower angle/conditioning estimate.

## 2. The strengthened, explicit Fourier approximation

Retain T1D-AB's **same** auxiliary even B-spline window `chi=1` on I, `a=241/120`, `k=14`, `C=3360`. Set the slightly larger cutoff

\[
\boxed{N=1{,}380{,}000,\qquad \operatorname{rank} A_N\le 2N+1=2{,}760{,}001.}\tag{AC.3}
\]

For `e∈E` let `h_e=S_B^*(T_B|_E)^{-1}e` (this is a new symbol, **not** the Weil form `h`). Then `S_B h_e=e` on I and `||h_e||<=theta^{-1/2}||e||`. Let

\[
A_N e=\left[P_N(\chi S_B h_e)\right]|_I,\qquad F_Ne=J A_Ne.
\tag{AC.4}
\]

With `omega_n=pi n/a` and `Omega=pi N/a-B`, the imported exact Fourier-tail calculation from AB and its termwise derivative give

\[
\begin{aligned}
\|S_Bh-P_N(\chi S_Bh)\|_{L²(I)}^2
&\le U_0\|h\|²,\\
\|\partial_x(S_Bh-P_N(\chi S_Bh))\|_{L²(I)}^2
&\le U_1\|h\|²,\\
U_0&=\frac{4B}{\pi²(2k+1)\Omega}\left(\frac{C}{\Omega}\right)^{2k},\\
U_1&=\frac{4B}{\pi²}\left(\frac{C}{\Omega}\right)^{2k}
\left(\frac{1+B²}{(2k+1)\Omega}+\frac{2B}{2k}+\frac{\Omega}{2k-1}\right).
\end{aligned}\tag{AC.5}
\]

The derivative estimate follows by Parseval with weight `1+omega_n²` (and thus bounds derivative alone), expanding `(B+t)²t^{-2k-2}` and integrating each monotone tail. The exact-rational audit verifies `Omega>61000` using T1D-AB's `pi>314159/100000` and gives

\[
\boxed{U_0/\theta<2^{-80},\qquad U_1/\theta<2^{-36}.}\tag{AC.6}
\]

Consequently `z=e-A_Ne` has `||z||_2<2^-40||e||_2` and `||z'||_{L²(I)}<2^-18||e||_2`. **This is an approximation of actual restricted prolate vectors**, not an assertion that the unsmoothed trigonometric functions are Weil-core tests.

## 3. Complete Gamma/logarithmic form-energy control

For a function `v∈H¹(I)` let `v^0` denote its zero extension to the whole line. It need not belong to `H¹(R)` because of its boundary jumps. It does, however, belong to the logarithmic form domain. With `t=||v||_{L²(I)}`, `d=||v'||_{L²(I)}`, a one-dimensional trace estimate and integration by parts yield

\[
|\widehat{v^0}(\xi)|\le\frac{\alpha(t,d)}{\sqrt{2\pi}|\xi|},\qquad
\alpha(t,d)=2\sqrt{t²/4+2td}+2d \quad(\xi\ne0).
\tag{AC.7}
\]

For every real `R>0`, Parseval and exact integration of the frequency tail imply

\[
\boxed{\int\log(1+|\xi|)|\widehat{v^0}(\xi)|²d\xi
\le t²\log(1+R)+\frac{\alpha(t,d)^2}{\pi}
\left(\frac{\log(1+R)}R+\log(1+1/R)\right).}\tag{AC.8}
\]

Define `||v||_E²=\int (65+log(1+|xi|))|hat(v^0)|²d xi`. This positive Hilbert norm controls the entire signed original zero-pole Fourier form, including the exact Gamma part and all signed primes. In (AC.8), with `t=2^-40`, `d=2^-18`, `R=2^40`, an exact rational square comparison gives `alpha<10^-5`, while `log(1+2^40)<29`. Hence

\[
\boxed{\|e-A_Ne\|_E<10^{-10}\|e\|_2.}\tag{AC.9}
\]

For any `e∈E`, the bandlimited source representation gives `||e'||_{L²(I)}<=B theta^-1/2||e||=2^40||e||`. Taking `R=2^80` in (AC.8) proves `||e||_E<22||e||_2` by exact inequalities (`alpha(1,2^40)<2^42`, `log(1+2^80)<57`). These control **the complete Gamma energy at the original physical boundary**; pretending the zero extension is `H¹(R)` would be invalid.

## 4. The eight physical constraints and exact source-form transfer

The normalized two pole carriers `1_I cosh(x/2)` and `1_I sinh(x/2)` have interior normalized derivative `<1` (the odd one follows from `sinh 2>10/3`). Equations (AC.7)–(AC.8) show that each has `E`-norm `<9`. For the selected sources, T1D-AA imports `||u_sigma'||<37000`, hence `||u_sigma||_E<193`. For the four fixed central/satellite sources, `||v_i'||<46`, hence `||v_i||_E<11`. Orthogonality in the physical metric and exact integer arithmetic yield

\[
\boxed{\left(2\cdot9²+2\cdot193²+4\cdot11²\right)^{1/2}<275,
\quad\|(1-J)v\|_E\le275\|v\|_2.}\tag{AC.10}
\]

Therefore for all `e∈E`, with `f=Je`, `g=F_Ne=JA_Ne`,

\[
\boxed{\begin{aligned}
\|f-g\|_2&<2^{-40}\|e\|_2,\\
\|f-g\|_E&<2^{-31}\|e\|_2,\\
\big|\|f\|_2²-\|g\|_2²\big|&<2^{-38}\|e\|_2²,\\
|q_\mathrm{prime}[f]-q_\mathrm{prime}[g]|&<10^{-10}\|e\|_2²,\\
|q_\Gamma[f]-q_\Gamma[g]|&<\frac1{3{,}000{,}000}\|e\|_2²,\\
|q_0[f]-q_0[g]|&<\frac1{3{,}000{,}000}\|e\|_2²,\\
|\mathfrak h[f]-\mathfrak h[g]|&<\boxed{\frac1{400{,}000}\|e\|_2²}.
\end{aligned}}\tag{AC.11}
\]

The last line uses only the imported coarse global bounds `||R||<1785`, `||T||<1,073,590` and the exact rational comparison `1/3,000,000+(1785+1,073,590)(201/100)2^-40<1/400,000`. No positive sign is asserted. These are *quadratic-form differences measured against the representing `e`*, not uniform errors measured against `||Je||`.

**Source legitimacy.** The projected `f,g` are in `Z` and have finite `E`-norm, but are generally not compact smooth; this note **does not evaluate the original Weil distribution on them as literal tests**. On bounded I, the zero extensions of functions in `H¹(I)` lie in `H^s(R)` for each `0<s<1/2`. So do the two sharp pole carriers and `Je,JA_Ne`. By standard `H^s` density of `C_c^∞(I)` when `s<1/2`, then correcting the eight independent physical moments using a fixed smooth biorthogonal family, both `f` and `g` have sequences in the genuine `Z_s` converging in `H^s`, hence in `E`-norm. The continuous closed-form extension of (AC.1) to these vectors is therefore the unique limit of **complete original smooth Weil-form values**. Any strictly negative value of `q_0` on such an exact pole-null limit can be carried to a smooth pole-null negative witness; similarly for the second Schur remainder (and finite completed squares). Quantitative smooth cutoff/moment-lift condition numbers for **computing** individual witnesses remain an additional certification obligation.

## 5. Exact structured finite matrices, not replacement Gamma physics

Use the finite Fourier functions `phi_n(x)=1_I(x)e^{i omega_n x}/sqrt(2a)`, `|n|<=N`, which lie in the form closure, not in `C_c^∞(I)`. Put `Phi c=sum_n c_n phi_n`. Then the **physical** matrix is the exact sinc Gram

\[
G_{nm}=\langle\phi_n,\phi_m\rangle
=\begin{cases}\dfrac{\sin(2(\omega_n-\omega_m))}{a(\omega_n-\omega_m)},&n\ne m,\\[3pt]\dfrac2a,&n=m.\end{cases}\tag{AC.12}
\]

The *complete original signed prime* matrix has, for `s=log r` and `kappa_r=Lambda(r)/sqrt r`, the closed form

\[
\boxed{Q^{\rm prime}_{nm}=-\sum_{r=p^j\le54}
\frac{2\kappa_r}{a}\cos\!\left(\frac{(\omega_n+\omega_m)s}{2}\right)
\frac{\sin((\omega_n-\omega_m)(4-s)/2)}{\omega_n-\omega_m};}\tag{AC.13}
\]

interpret `sin(t L/2)/t=L/2` at `t=0`. This is a sum of 24 modulated sinc kernels, suitable for Toeplitz/Hankel-style matrix-free operations, retaining every original prime sign and weight. The **complete Gamma** matrix is the convergent exact digamma integral

\[
\boxed{Q^\Gamma_{nm}=\frac1{\pi a}\int_{\mathbb R}m_\Gamma(\xi)
\frac{\sin 2(\xi-\omega_n)}{\xi-\omega_n}
\frac{\sin 2(\xi-\omega_m)}{\xi-\omega_m}\,d\xi.}\tag{AC.14}
\]

Set `M_{in}=<phi_n,c_i>` for the eight physical constraint carriers and `J Phi=Phi-CM`. The authentic second Schur matrix on these projected columns is formed from `Q^Gamma+Q^prime`, the exact eight-moment correction and the **actual** source operators `R,T`, not an assumed diagonal metric or an arbitrary surrogate Gamma multiplier. The `G` and prime kernels are explicit; the Gamma integral and the six Schur residuals still need rigorous blockwise interval bounds for a spectral enclosure.

## 6. Low/high Schur gate and what AC does *not* prove

The conditioned `H_hi` has the imported second-Schur lower bound `d_*=1233/1750` on its smooth core. For a genuine low form-domain vector `f∈L_low`, define its mixed functional against `H_hi` by `mathfrak h(f,y)`. The exact low Schur form is

\[
\mathfrak S_{\rm low}[f]=\mathfrak h[f]
-\sup_{y\in H_{\rm hi}\cap\mathrm{dom}(\mathfrak h)}
\{2\Re\mathfrak h(f,y)-\mathfrak h[y]\}.
\tag{AC.15}
\]

For each projected finite Fourier vector, the complete `K_h f` is in physical L²: its Fourier transform is the original `k(xi) fhat(xi)` with logarithmic times `1/|xi|` tail, together with bounded rank-six corrections. Thus the *safe* (possibly very pessimistic) lower estimate is

\[
\boxed{\mathfrak S_{\rm low}[f]
\ge\mathfrak h[f]-d_*^{-1}\|P_{H_{\rm hi}}K_{\mathfrak h}f\|_2².}\tag{AC.16}
\]

Equations (AC.12)–(AC.16) specify a structured, complete-original-source spectral **screening** operator and a route to effective matrix-free/interval bounds. A negative value of the right-hand side alone does **not** prove negativity; a strict negative value of the complete original `q_0` on an actual pole-null smooth certificate *would*.

**Two independent gaps remain even after (AC.11):**

1. The map `J|E` needs a certified lower nonzero singular value (or a quotient/frame conditioning theorem). AC's absolute error is relative to `||e||`, not `||Je||`; near the at-most-eight bad physical directions, a relative low-block spectral comparison cannot be inferred.
2. The low/high **mixed** second-Schur block in (AC.15), including the full Gamma multiplier and rank-six penalties, needs an actual upper enclosure or an effective spectrum comparison. Positivity of the low compression alone would not suffice.

**Audit:** `python3 src/audits/canonical_return_t1d_ac_weighted_transfer.py` checks the exact rational inequalities (AC.3), (AC.6), (AC.9)–(AC.11), imported fixed norm bounds, and all conservative trace/log cutoffs. It does **not** certify the T1D-AB Fourier coefficient proof, Gamma/digamma identity, Sobolev trace/density, actual operator spectral values, the low-frame angle, or negative directions. The finite matrix kernels are exact mathematical identities, not numerically enclosed entries. Keep the existing draft PR #2 in research jurisdiction.

**Next — CR-INT-T1D-AD:** certify the eight-angle/quotient conditioning and the high/low mixed Gram, or produce an interval-enclosed **actual smooth** negative witness. Do not turn this screening estimate into a claim that the Morse index is zero.