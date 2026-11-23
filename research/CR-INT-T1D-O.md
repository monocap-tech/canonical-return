# CR-INT-T1D-O — Exact finite-window prime-translation contraction and a 2^29 prolate gate

**Canonical Return · November 23, 2026 (America/Los_Angeles).**  
**Disposition:** **NEW SOURCE-LEVEL ANALYTICAL ESTIMATE** with executable exact-rational arithmetic audit; not an independent Lean or peer-reviewed certificate, not a calculation of the remaining finite Schur matrix, not complete positivity on the physical complement, not an identification with the semilocal Sonin space, and not RH. This refines the existing [T1D-M](CR-INT-T1D-M.md) and [T1D-N](CR-INT-T1D-N.md) conditioning approach. Unlike T1D-N's abstract rank obstruction at its *frozen* cutoff, this note derives a smaller **source-valid cutoff** by incorporating the actual physical support in the arithmetic operator bound.

## 1. Immutable mathematical source and the oversight in the previous *bound*

Fix the **same original source**: \(I=(-2,2)\), \(H=L^2(I,dx)\), the smooth test core \(\mathcal D=C_c^\infty(I;\mathbb C)\), and the two-pole-null core
\[
\mathcal D^0=\{f\in\mathcal D:\ \int_I f(x)e^{x/2}dx=\int_I f(x)e^{-x/2}dx=0\}.
\]
The complete signed Weil form is \(q(f,g)=W(f*\widetilde g)\), **linear in its first variable**. The original zeta pole contributions vanish for \(f,g\in\mathcal D^0\). The exact original Gamma contribution and negative prime-power terms are those of Suzuki, [*Weil's quadratic form via the screw function*, arXiv:2606.09096v4, §1.1, §2.3 eq. (2.5), §2.4 eq. (2.8)](https://arxiv.org/html/2606.09096v4). This does not alter the form, the denominator metric or the admissible core.

Write \(P_I\) for multiplication by \(\mathbf1_I\) on \(L^2(\mathbb R)\), \((T_tf)(x)=f(x-t)\), \(\delta_n=\log n\) and \(\kappa_n=\Lambda(n)/\sqrt n\). On the exact fixed support the **genuine arithmetic part**, as a bounded self-adjoint operator, is
\[
\boxed{\mathsf P_I=-\sum_{2\le n\le54}\kappa_n\,P_I(T_{\delta_n}+T_{-\delta_n})P_I\big|_H.}\tag{O.1}
\]
Only the 24 genuine prime powers contribute, because \(54<e^4<55\); all other \(\Lambda(n)\) vanish. For every smooth core test \(q_{\rm prime}(f,f)=\langle\mathsf P_If,f\rangle\). The earlier T1D-M bound \(q_{\rm prime}(f,f)>-2S_2\|f\|_2^2\), \(S_2=\sum_{n\le54}\kappa_n<61/5\), was safe but replaced the **compressed shift** norm by its larger **whole-line** norm of \(2\). The next theorem removes precisely this overestimate. No imaginary arithmetic cancellation, nontrivial zeta estimate or unverified domain transfer is used.

## 2. An exact operator norm for one original prime translation

**Lemma (finite-window shift chain).** For \(0<t<L\), let \(I_L\) be any interval of length \(L\) and
\(\mathsf S_t=P_{I_L}(T_t+T_{-t})P_{I_L}|_{L^2(I_L)}\). If \(L/t\notin\mathbb Z\), then
\[
\boxed{\|\mathsf S_t\|_{L^2(I_L)\to L^2(I_L)}=2\cos\!\left(\frac{\pi}{\lceil L/t\rceil+1}\right).}\tag{O.2}
\]
The same formula with the *essential maximum* number of vertices describes the integer-ratio edge case after replacing \(\lceil L/t\rceil\) by \(L/t\); the present problem has no edge case. For \(t\ge L\) the shift compression is identically zero.

**Proof.** Parametrize \(x\in I_L\), up to null sets, by its residue \(r=x\bmod t\) and the finite consecutive set of integers \(j\) for which \(r+jt\in I_L\). This is a unitary \(L^2\)-direct-integral disintegration of \(L^2(I_L)\) into fibers \(\mathbb C^{m(r)}\), and \(\mathsf S_t\) becomes the nearest-neighbor path adjacency matrix \(A_{m(r)}\) on each fiber. Its eigenvalues are \(2\cos(k\pi/(m(r)+1))\), \(k=1,\ldots,m(r)\), by the elementary sine eigenbasis. The largest fiber has \(\lceil L/t\rceil\) vertices on a set of residues of positive measure whenever \(L/t\notin\mathbb Z\). Taking the essential supremum proves (O.2). This disintegration remains wholly within the physical \(L^2\) metric; a translated test's support is not *presumed* to remain inside the window.

For our \(L=4\), the strict rational window \(54<e^4<55\) establishes the chain length without logarithmic rounding: the largest integer \(m_n\) with \(n^{m_n-1}<e^4<n^{m_n}\) is
\[
\boxed{m_2=6,\quad m_3=4,\quad m_4=m_5=m_7=3,\quad m_n=2\text{ for the other 19 prime powers }8\le n\le54.}\tag{O.3}
\]
These are **exact**, including the near-boundary \(n=53\) shift. In particular, each \(n\ge8\) contributes an operator of norm **one**, not two.

## 3. A purely rational Perron bound on every compressed shift

The path matrices are real symmetric and entrywise nonnegative. If \(v>0\) and \(A_m v\le c v\) entrywise, then \(\|A_m\|\le c\): conjugate by the positive diagonal matrix with entries \(v_j\) and apply the maximum row-sum norm, using symmetry of the original matrix to identify the spectral radius with its largest eigenvalue.

The following **exact integer witnesses** require no trigonometric numerical computation:

| Path length | Positive witness \(v\) | Verified max \((A_mv)_j/v_j\) | Upper shift norm |
|---|---|---|---|
| 6 | \((20,36,45,45,36,20)\) | \(65/36\) | \(65/36\) |
| 4 | \((8,13,13,8)\) | \(13/8\) | \(13/8\) |
| 3 | \((70,99,70)\) | \(99/70\) | \(99/70\) |
| 2 | \((1,1)\) | \(1\) | \(1\) |

Combine (O.1)–(O.3) and the operator triangle inequality:
\[
\begin{aligned}
\|\mathsf P_I\|&\le \sum_{n\le54}\kappa_n\|\mathsf S_{\log n}\|\\
&\le S_2+\left(\frac{65}{36}-1\right)\kappa_2
+\left(\frac{13}{8}-1\right)\kappa_3
+\left(\frac{99}{70}-1\right)(\kappa_4+\kappa_5+\kappa_7)\\
&\boxed{<\frac{55}{4}=13.75.}
\end{aligned}\tag{O.4}
\]
The final **strict** inequality is certified using exact-rational outward intervals for \(\log p\) (the atanh series) and \(\sqrt n\) (integer-square-root bounds). The audited weighted upper sum is approximately \(13.72936838121\), but **the proved source estimate used subsequently is only \(55/4\)**. Compared with the former whole-line bound \(2S_2\approx24.38329\), the validated support contraction reduces the arithmetic worst-case norm by more than **10**. No assumption about the joint spectrum of distinct noncommuting compressed shifts is necessary. This may still be a substantial overestimate of the actual \(\|\mathsf P_I\|\).

**Why the window matters:** the full-line Fourier multiplier of the same *finite* arithmetic sum equals \(-2\sum_{n\le54}\kappa_n\cos(\xi\log n)\), which takes the value \(-2S_2\) at \(\xi=0\). So merely improving a **uniform whole-line pointwise estimate** cannot yield the support-specific gain in (O.4); the gain follows from honest restriction of each shift to \(I\).

## 4. A smaller prolate threshold for the *same* actual signed Weil form

Retain T1D-M/N's logarithmic time-band operator
\(\mathsf T_B=P_I\mathscr F^{-1}\mathbf1_{[-B,B]}\mathscr FP_I|_H\), with \(\langle\mathsf T_Bf,f\rangle=\int_{-B}^B|\widehat f(\xi)|^2d\xi\). Let \(\theta=1/128\), \(V_B^>=\operatorname{Ran}\mathbf1_{(\theta,1]}(\mathsf T_B)\) and
\[
H_B^F=\{y\in\mathcal D^0\cap E^\perp:\ \langle y,v\rangle=0\text{ for every }v\in V_B^>\},
\]
where \(E=\operatorname{span}\{u_+,u_-\}\) is the **unchanged original** two-dimensional positive test space. Every element of \(H_B^F\) is still a genuine compactly supported smooth two-pole-null source; it is **not** being replaced by a nonsmooth prolate eigenfunction. The spectral construction gives \(\int_{|\xi|\le B}|\widehat y(\xi)|^2d\xi\le\theta\|y\|_2^2\), and the subspace remains infinite-dimensional because the number of constraints is finite.

Import from the independently stated Suzuki/T1D-M Gamma analysis (not proved merely by the present finite checker):
\[
\boxed{m_\Gamma(\xi)>\log|\xi|-\frac{49}{10}\quad(|\xi|\ge64),\qquad m_\Gamma(\xi)>-6\quad(\xi\in\mathbb R).}\tag{O.5}
\]
Together with the **new** uniform *physical* prime lower bound \(\langle\mathsf P_Iy,y\rangle>-(55/4)\|y\|^2\), this yields for any \(B\ge64\):
\[
\boxed{\frac{q(y,y)}{\|y\|_2^2}>
(1-\theta)\left(\log B-\frac{49}{10}\right)-6\theta-\frac{55}{4}}
\qquad (0\ne y\in H_B^F).\tag{O.6}
\]
**New freeze:** choose \(B=2^{29}=536{,}870{,}912\), keeping \(\theta=1/128\) and the original source radius \(a=2\). The exact rational series gives \(\log2>693/1000\), hence
\[
\boxed{q(y,y)>\frac{164019}{128000}\,\|y\|_2^2>1.28\,\|y\|_2^2\quad(0\ne y\in H_{2^{29}}^F).}\tag{O.7}
\]
This is a direct improvement over the T1D-M \(B=2^{44}\) cutoff: the required frequency threshold falls by \(2^{15}=32{,}768\), while the *proved positive lower margin increases* from \(0.9451875\) to \(1.2813984375\). The improvement is derived from the true original arithmetic source; no changed Weil form, mass metric or pole-null class has been introduced.

**Conditioned Schur comparison — explicitly imported prior obligation.** If the separate T1D-I source-level bound \(\|\mathcal R\|<155/182\) is accepted after independent mathematical review, then on exactly the same \(H_{2^{29}}^F\),
\[
\boxed{q(y,y)-\langle\mathcal R y,y\rangle>
\frac{5005729}{11648000}\,\|y\|_2^2>\frac25\|y\|_2^2.}\tag{O.8}
\]
The new unconditional *restricted* lower bound (O.7) **does not depend on T1D-I**. Formula (O.8) does. Neither reaches all of \(F\), because an immense collection of legitimate low-band directions remains excluded.

## 5. Honest cost of the remaining low-band block

The old T1D-N two-sided time-band trace estimate remains exactly applicable, but at this **new** \(B\). With \(m_B=\dim V_B^>\), trace \(4B/\pi\) and defect
\(D_B\le(2/\pi^2)(\log(4B)+3/2)\), the same exact-rational comparison \(3<\pi<22/7\), \(\log2<347/500\) yields
\[
\boxed{683{,}290{,}247\le m_{2^{29}}\le715{,}828{,}537.}\tag{O.9}
\]
There are **at least 683,290,243 additional independent constraints** after the original four physical directions are accounted for. In particular, T1D-N's min-max rank obstruction is **not defeated**; its cost is reduced because the cutoff is smaller. The upper bound in (O.9) is more than **31,000 times smaller** than T1D-N's *minimum* excluded rank at \(B=2^{44}\), although hundreds of millions of low modes are still too many for an explicit brute-force Schur matrix.

This improvement has a concrete interpretation: any method that continues bounding primes by the **whole-line** \(2S_2\) is unnecessarily throwing away more than ten units of coercivity on this physical interval. But any method using only separate Gamma lower bounds, the new \(55/4\) arithmetic bound, and the same frequency-leakage threshold still cannot conclude positivity of the omitted low block. The **matrix (M.18)** and the two original pole directions remain unsatisfied obligations.

## 6. Relation to published spectral frameworks, exact claim status and next test

Suzuki already proves that the original fixed-window Weil form has a suitable form closure and discrete spectral realization (see his §§1, 3–4 and the compact embedding in Proposition 4.1). **We do not claim to rediscover this fact**, identify his specific completed operator with \(\mathsf T_B\), or prove a new global spectral theorem here. Our new point is the explicit **fixed-window arithmetic contraction bound (O.4)** and its direct substitution into the independent T1D-M prolate/Schur gate. The Andrews–Connes–La Gatta semilocal construction (arXiv:2611.13129v1, §8.6) has additional transport and conditioned-space hypotheses; the logarithmic \(\mathsf T_B\) and its physical metric have **not** been identified with their Sonin/prolate conditioning operator.

**Analytically established on the frozen original source:** shift-chain disintegration and exact chain lengths; the rational Perron upper bound \(\|\mathsf P_I\|<55/4\); a factor-\(2^{15}\) cheaper frequency threshold for an infinite-dimensional positive two-pole-null smooth subspace; and the new finite rank bracket. **Depends on imported analytic source facts:** Suzuki/T1D-M Gamma multiplier bounds and standard prolate trace/defect identities. **Depends additionally on unreviewed T1D-I:** the improved Schur remainder (O.8). **Not established:** any actual low-mode Weil Schur matrix sign, a full complement coercivity estimate, the two excluded original pole directions, transfer of completed Suzuki/Connes/Sonin domains, support-independent positivity, Lean certification, or RH.

**Verification:** `python3 src/audits/canonical_return_t1d_o_window_prime_norm.py` checks all exact prime-power and chain-length assignments, rational Perron vector inequalities, rigorous outward \(\log p/\sqrt n\) intervals, coercivity/Schur fractions and the T1D-N prolate rank arithmetic. The **operator direct-integral proof**, the original Weil-form/Gamma identity and the time-band spectral theorem are human mathematical arguments, not artifacts of the finite checker. No actual prolate eigenvalues are numerically enumerated. This note is a research result pending external independent review.

**Next CR-INT-T1D-P:** improve the \(55/4\) estimate by exploiting *joint* original prime translations rather than a triangle inequality over separate compressed shifts, while preserving the physically bounded interval. An independently validated certified operator-norm upper bound for \(\mathsf P_I\) would translate **exponentially** into a smaller cutoff in (O.6). In parallel, construct a **form-adapted** finite low-energy block or transfer theorem rather than brute-force enumeration of the still enormous \(m_B\). Any proposed joint estimate must pass the original physical metric, exact boundary truncation, all 24 \(\Lambda(n)\) weights, the two-pole core and the full Gamma term.