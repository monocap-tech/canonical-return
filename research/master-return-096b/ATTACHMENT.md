# RETURN096-B — Original periodic comparison and normalization

This note establishes the analytic comparison for the original finite-interval
Gamma/prime/pole bilinear form used by RETURN048. The separate arithmetic
verifier must confirm the entire signed ledger; no numerical PASS is presumed.

## Form and periodization

Put I=(-5/4,5/4), D=5/2, L=5, beta_m=2m+1/2. Extend f by zero on
the rest of the circle of circumference L. For real f on I the
explicit-formula realization under audit is

    W(f)=g0 ||f||²
       + sum_{m>=0} [2||f||²/beta_m
           - integral_I integral_I exp(-beta_m|x-y|) f(x)f(y) dxdy]
       - 2 sum_{n in N} (log p)/sqrt(n) integral_R f(x)f(x+log n) dx
       + 2 (integral_I exp(x/2)f(x)dx)(integral_I exp(-x/2)f(x)dx),

where g0=psi(1/4)-log(pi)=-EulerGamma-pi/2-3log(2)-log(pi)
and N={2,3,4,5,7,8,9,11}, n a power of p. All real-line integrals use
zero extension. This is the realization underlying the authenticated
RETURN048 bilinear entries; literal attachment to the chosen external
explicit-formula definition must still be checked.

The periodic exponential convolution has multiplier
2 beta/(beta²+r²) at r=2*pi*k/L. For |x-y|<=D the exact image
difference between the periodized and nonperiodized exponential is

    2 exp(-beta L)/(1-exp(-beta L)) cosh(beta(x-y)).

Prime shifts have log(n)<D and L-log(n)>D, so no wrapped image contributes
on I. This follows from the support intersection, not an approximation.

Write c_beta=integral f cosh(beta x), s_beta=integral f sinh(beta x).
The original form equals the periodic symbol form plus

    sum_m 2 exp(-beta_m L)/(1-exp(-beta_m L))
                (c_beta_m²-s_beta_m²)
      + 2(c_(1/2)²-s_(1/2)²).

In even parity all correction terms are nonnegative. In odd parity on
the literal RETURN096-A complement, s_(1/2)=0. For m>=1,

    ||sinh(beta x)||² = sinh(beta D)/(2 beta)-D/2
                      <= exp(beta D)/(4 beta).

Hence the total negative correction is bounded by

    ||f||² exp(-(5/2)(L-D)) /
         [5(1-exp(-5L/2))(1-exp(-2(L-D)))] < .001 ||f||².

The bound follows by summing a positive geometric majorant. The
script independently verifies it with finite positive exponential
series. No negative pole or first image term is discarded without
the sinh constraint.

## Original signed Fourier symbol

The exact periodic symbol is

    s(k)=g0+sum_{m>=0} 2 r_k²/[beta_m(beta_m²+r_k²)]
              -2 sum_n (log p)/sqrt(n) cos(r_k log n).

The Gamma sum is monotone in |k|. For r>0, the summand as a function
of beta is decreasing and the stride-two tail after 32 terms is at
least its half-density integral, namely
(1/2)log(1+r²/(129/2)²). The lower g0 bound -5.373 uses
EulerGamma <= H_4096-log(4096) and directed pi/log bounds.
The alternate high-frequency lower bound
g(r) >= log(r)-1.85-2/r follows from the digamma harmonic sum:
for a=1/4,b=r/2 the function (t+a)/((t+a)²+b²) has total
variation at most 1/b, so its unit-step sum minus integral is at
most 2/r. Taking the limit gives
Re psi(a+ib) >= (1/2)log(a²+b²)-2/r >= log(r/2)-2/r.
Use log(2*pi)<1.85. The k<=1 bound -6 follows from g(r)>=g0.
The rational verifier reproduces RETURN046's alternate phase and
log arithmetic, verifies the whole finite ledger, and checks the
high-frequency cutoffs with the monotone bound minus the full
prime triangle upper bound.

## Physical normalization and restoration

Let a_k and b_j be the real Fourier coefficients of f on I in its
even-grid and odd-grid bases in the given parity, with a rephased
by (-1)^(k/2). The two systems are each complete, since the
odd complex exponentials are obtained from the even ones by
unitary modulation on I; combining positive/negative frequencies
gives the real parity bases. The nonconstant squared norms are 5/4.
The constant even-grid coefficient vanishes on the complement.

On the circle the zero-extended coefficients equal half the
corresponding interval coefficients in each system. Consequently

    W_periodic(f) = (5/8) [sum_even s(k) a_k²
                                      + sum_odd s(j) b_j²],
    ||f||² = (5/4) sum a_k² = (5/4) sum b_j².

The ledger checks s(k)>=.1+w_k on every complementary even
coordinate and s(j)>=.1+v_j-d_j+u_j on every nonexcluded odd
coordinate, including all omitted carriers with zero weights.
The five unpaid odd coordinates vanish. Here v,d are RETURN057
and u is the additional RETURN075 restoration weight. This check
must add the two positive weights if a carrier occurs twice.
After paying .001 for images,

    W(f) >= (5/4) [ sum A_k a_k²
              + (1/2) sum_j (v_j-d_j+u_j)b_j² ],
    A_k=.099+w_k/2.

By the overlap identity of RETURN096-A, b_j is, up to sign,
(4/pi) times the Cauchy row applied to a. This gives the precise
signed update coefficients 8(v_j-d_j)/pi². For odd restoration,
let h be the unscaled Cauchy functional vector on the 500 carriers,
B the rational eight-row coefficient bank and
D0=diag(8u_j/pi²)>0. With pi_+>=pi,
R=B diag(pi_+²/(8u_j)) B^T >= B D0^-1 B^T.
For every h, (Bh)^T R^-1 (Bh) <= h^T D0 h:
apply orthogonal projection onto row space of B D0^-1/2,
then inverse monotonicity. RETURN094-A verifies the literal R
and its positive definiteness. Thus the restoration is a lower
bound on the available positive energy, not an extra payment.
Combining with the signed updates gives exactly the model D
whose inverse is used in RETURN095-A, with physical factor 5/4.
This explains the 4/5 in the inverse payment, not a separate
empirical scaling.

The even model uses A=.174 on k>=32768; the odd model raises
A to .519 on 65536<k<=262144 and 1.199 on k>262144. The
three monotone cutoff inequalities in the script verify these
lower coefficients. Keeping the weaker value at a boundary
is safe. The baseline negative-update cap from RETURN095-B
continues to hold when A is strengthened.

## Form domain

On the circle the Gamma symbol is bounded above and below by
log(2+|k|) plus constants (split the positive Gamma series
at beta≈|r|, compare the tail with its integral). Prime terms
are bounded multipliers. The image correction is a bounded
quadratic form: finitely many first terms, then absolutely
summable operator norms by the same geometric estimate for cosh
and sinh. The pole term is bounded on the fixed interval.
Thus the natural semibounded zero-extension form domain is the
restriction of the logarithmically weighted Fourier space
sum_k log(2+|k|)|fhat(k)|²<infinity to functions supported in I.

Each retained basis function is smooth on the closed interval
and piecewise smooth after zero extension, with circle
Fourier coefficients O(1/|k|) by integration by parts.
It therefore belongs to this form domain. For a smooth compactly
supported f in I the orthogonal projection p into the finite
retained space belongs to the domain; so does q=f-p.
No vanishing boundary trace is imposed on p or q.
The same reasoning applies to arbitrary f in the stated
logarithmic form domain. Parseval and monotone convergence of
the nonnegative Gamma increments justify the periodic symbol
identity there. The image, prime and pole terms extend by
boundedness.

This note does not establish a native formal theorem, the
external explicit-formula identity, all mixed source identities
on this domain, density of a specified external test class, or
arbitrary-support globalization. In particular a failed
full-ledger check must be reported as a blocker, not repaired
by silently reallocating a weight.
