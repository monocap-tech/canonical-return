# CR-INT-T1D-AP — All-prime parity-circle diagonal comparison

**Disposition:** derived analytical comparison, executed integer/Fraction arithmetic audit PASS. Conditional on the original complete pole-null Gamma/Weil identity and AM arithmetic norm. No independent review, Lean, remote CI, whole-source positivity or RH. No earlier mixed certificate is subtracted.

## 1. Exact physical embedding and arithmetic

Keep I=(-2,2), D^0=C_c^infty(I) intersect the two Mellin pole-null constraints, physical L2 metric, and q0=Gamma-A. Embed by zero extension E into J=(-4,4). For reflection parity sigma=+ use periodic boundary condition at ±4; for sigma=- use antiperiodic boundary condition. On these spaces use normalized real cosine and sine bases, respectively, with frequencies omega_j=j*pi/8 for nonnegative even j and positive odd j. This is an isometry on each parity; the zero extensions of original smooth tests belong to the twisted circle smooth core. In general their range is a proper closed support-restricted subspace, not the whole parity circle.

Define twisted translations on the circumference-8 circle. Every shift s=log(n), n=p^k<=54, satisfies 0<s<4. If x,y in I and x-y=±s modulo 8, then x-y=±s literally, since |x-y|<4 and 8-s>4. Hence there is no wrapping contribution on physical tests, including antiperiodic signs. The complete original arithmetic quadratic form equals the diagonal circle form with a(omega)=2 sum_(n=p^k<=54) kappa_n cos(omega log n), kappa_n=log(p)/sqrt(n). No commensurability assumption or modified prime weight is used.

## 2. Parity-dependent complete Gamma lower comparison

For a>0 the real-line resolvent kernel is exp(-a|x-y|)/(2a). The periodic (+) and antiperiodic (-) circular kernels on the physical support are the absolutely convergent image sums with weights tau^l, tau=+1 or -1. Since |x-y|<4, the difference from the original real-line kernel equals

  [tau exp(-8a)/(a(1-tau exp(-8a)))] cosh(a(x-y)).

Indeed pair images at ±8l, l>=1, and sum their geometric series. On real or complex even f, the cosh(ax) moment survives and sinh(ax) moment vanishes, so the displayed correction is nonnegative for tau=+1. On odd f the cosh moment vanishes and the negative sinh outer product survives; with tau=-1 the correction is again nonnegative. Thus on each original parity the compressed whole-line resolvent is dominated by the corresponding circular resolvent in quadratic-form order. This does NOT assert dominance on the whole unseparated circle.

The complete digamma expansion with a_k=2k+1/2 is m_Gamma(w)=m0+sum_(k>=0)(2/a_k) w²/(a_k²+w²). Each summand is a positive identity minus a positive multiple of the resolvent. Consequently finite partial Gamma energies dominate the corresponding circular diagonal energies; monotone convergence of the nonnegative shifted energies gives the comparison for the full series. Both are finite on smooth zero extensions. Thus q0[f]>=sum_j (m_Gamma(omega_j)-a(omega_j)) |c_j|². The Mellin constraints are imposed on the original physical source, not replaced by a circle constraint.

## 3. Convex arithmetic comparator and exact tail

AM gives A<=399/50 I on the physical support. Retain 2/5 of the exact circle arithmetic symbol and use AM for the remaining 3/5. This is a convex lower comparison, not the exact original q0 or an assertion that the norm on the full circle is 7.98:

  q0[f] >= sum_j [m_Gamma(j*pi/8) - (2/5)a(j*pi/8) - (3/5)(399/50)] |c_j|².

All 24 signed prime powers contribute in the retained part, and all are covered by AM in the bounded part. The discarded arithmetic correction is (3/5)(399/50 I-A)>=0 on physical tests; the Gamma image correction is also nonnegative. Set eps=1/20000. For j>=1024, use the AN/AO digamma convex trapezoid lower bound, with x=j*pi/16, T=129/4, t_k=k+1/4:

  m_Gamma >= m0 + sum_(k<32)1/t_k - (sum_(k<32)t_k+T/2)/x² + log(x/T)+1/(2T).

It follows from x²/[t(t²+x²)] being decreasing and convex, its exact integral tail, and 1/(1+u)>=1-u. Freeze the inverse-square error at j=1024 and bound pi, m0 and logs rationally, giving m_Gamma>=log(j)+C. The audit verifies at j=35,000,000 that this lower bound exceeds (2/5)sum(upper 2*kappa_n)+(3/5)(399/50)+eps. Hence all larger indices satisfy the margin, independently of their phases. Diagnostic Gamma floor 14.58532325 versus threshold 14.54136674.

## 4. Integer-grid frequency count

Reproduce `python3 ap_circle.py` with `am_prime.py`, `am_high.py`, NumPy and Python standard library. The full AM certificate is re-executed. Count all j<1024 as exceptional, irrespective of sign. For every 1024<=j<35,000,000, lower-bound log(j) by dyadic octave and a 1024-bin rational logarithm table, each rounded downward on a 10^9 grid. This gives a rigorous Gamma lower bound.

For each n, represent exp(i*pi*log(n)/8) and its 8192nd power by integer pairs on the same grid. The two seeds are certified using rational pi/log bounds, reduction modulo 2*pi before multiplication, Taylor remainder and a Lipschitz bound. Both seed complex errors are <2/D (D=10^9). Complex multiplication with coordinate floor rounding adds <2/D. For table power k<=8192, induction gives error <=8k/D, since (1+2/D)*8k/D+4/D<=8(k+1)/D. Write j=r+8192 h; combine the two table values by integer complex multiplication. The resulting phase error is <1/4000; high digit h<8192 throughout the scan. Weight intervals for 2*kappa have width <=2/D and total upper weights <25. The upper symbol is the integer sum of floor-rounded products plus ceil(D/150), since 25/4000+24*3/D<1/150. Integer arithmetic ranges are checked to prevent signed 64-bit overflow. No floating-point trigonometric evaluation or eigenvalue computation enters the accepted inequalities.

Every index is conservatively classified by the integer test for m_Gamma_lower-(2/5)a_upper-(3/5)7.98>=eps. The executed output gives 3223 even periodic and 3198 odd antiperiodic exceptional frequencies, total 6421. Last exceptional indices 1278966 and 2204993. There are no further exceptional indices in the finite scan; the infinite tail is certified analytically as above. The total is 1043 smaller than AO's 7464.

## 5. Source conclusion and next gate

On the original smooth pole-null source, impose vanishing of all exceptional circle coefficient functionals. Their number bounds the codimension, independently of linear dependencies after restricting to I. Parseval and the comparator give q0[f]>=eps||f||² there. Consequently n_-(q0;D^0)<=6421, and the parity bounds are 3223 and 3198. This is a bound on possible negative directions, not a claim of 6421 actual negative eigenvalues. It is a different exceptional space than AO's, not a common Fourier cutoff. The original form is not certified nonnegative on the exceptional sector or in its mixed interactions. No previous low-mode certificates are subtracted.

AQ: exploit the physical support restriction and the exact positive resolvent image deficit to compress the exceptional circle space. The full-circle lower comparator ignores that the zero-extension range is a proper subspace; all exceptional Fourier coefficients cannot be varied independently while high coefficients vanish. Use this constraint, and the positive arithmetic deficit, to seek a tighter finite Schur or inertia certificate. Do not infer positivity from the diagonal comparator or its finite compression alone. Whole fixed-window/global positivity, Suzuki–Connes transfer, Phase IV and RH remain open.