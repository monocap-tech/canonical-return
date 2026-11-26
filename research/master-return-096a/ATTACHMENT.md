# RETURN096-A — Literal retained/complement Fourier attachment

Work on I=(-5/4,5/4) in real L2, with Lebesgue measure.
The even-frequency real Fourier system on this interval is complete:
1, cos(2*pi*k*x/5), sin(2*pi*k*x/5), k=2,4,... .
Each nonconstant mode has squared norm 5/4; the constant has squared
norm 5/2. Reflection separates the even and odd sectors.

Let E be the nonnegative exceptional even indices, J the five positive
unpaid odd indices, and B the 44 selected positive even indices in the
authenticated RETURN059 chart. In each parity the retained space P is
spanned by the positive even indices in E union B and the five
corresponding odd modes with indices J, together with the constant in
even parity or sinh(x/2) in odd parity.

These 87 functions are linearly independent: the cosine and sine
functions have distinct complex exponential frequencies, with the
additional sinh function having the two distinct real exponents
+1/2,-1/2. A vanishing linear combination of finitely many distinct
exponentials on an open interval has every coefficient zero by the
Vandermonde system obtained from derivatives at zero. Hence P is closed
and every L2 function in the same parity has a unique orthogonal
decomposition p+q with p in P and q in P-perp.

## Literal constraint functionals

Write q=sum_k a_k sigma_k phi_k over positive even k outside E union B,
where sigma_k=(-1)^(k/2), phi_k is the cosine or sine for that parity,
and a is square summable. Parseval gives ||q||^2=(5/4) sum a_k^2.
Orthogonality to the omitted even-grid modes and (in even parity) to
the constant is already enforced by the index restriction.

For j odd positive put sigma_j=sin(pi*j/2), r_k=2*pi*k/5.
Product-to-sum and integration on I give

    <cos(r_j x), cos(r_k x)>
        = -(5/pi) sigma_j sigma_k j/(k^2-j^2)
    <sin(r_j x), sin(r_k x)>
        = -(5/pi) sigma_j sigma_k k/(k^2-j^2).

Consequently orthogonality to the five unpaid odd modes is exactly
the vanishing of the five Cauchy functionals
sum_k a_k j/(k^2-j^2) (even parity) or
sum_k a_k k/(k^2-j^2) (odd parity).
The factors -(5/pi) sigma_j are nonzero.

Direct integration gives

    <sinh(x/2), sin(r_k x)>
      = -(5/pi) sigma_k sinh(5/8) k/(k^2+25/(16*pi^2)).

The extra odd constraint is thus the same rational-kernel functional
with node -25/(16*pi^2). All these row sequences are in l2 (their
tails are O(k^-2) or O(k^-1)), so passage from finite sums to the L2
expansion is justified by Cauchy-Schwarz, not by pointwise convergence.

## Chart preservation

For pairwise distinct nodes t_1,...,t_n the identity

    1/prod_i(x-t_i) =
    sum_j 1/[(x-t_j) prod_{i!=j}(t_j-t_i)]

follows on multiplying by the product: the polynomial on the right,
of degree at most n-1, equals 1 at every t_j. Applying this identity
to every initial segment of the node list shows that the divided
difference rows used in RETURN063 are obtained from the original
Cauchy functionals by a triangular invertible row map. Its diagonal
entries are products of nonzero node differences, multiplied by the
nonzero powers of 116 and (even parity) the inverses of the positive
odd indices. The rational Z row chart is also invertible if its
checked rank is n. Thus the vanishing of the five or six archived
constraint rows is equivalent to membership in P-perp, not merely
a sufficient condition.

The exact finite index audit verifies that E union B removes precisely
the retained even-grid coordinates, that the finite complementary
coordinates below 32768 exhaust the remainder, and that the tail is
all positive even indices from 32768. It checks the common auxiliary
dimensions, signed carriers, constraint chart ranks, and a Machin
enclosure of the true sinh node.

This proves the stated L2 decomposition and identifies the constraint
kernels by elementary Fourier analysis, conditional only on the
finite receipt checks. The original zero-extension Weil form need not
be defined on every L2 function. In particular, orthogonal projection
onto P is not automatically a map of the desired form domain. The
next obligation is to establish the periodic-symbol comparison and
original form-domain decomposition, with all endpoint terms and the
5/4 physical normalization, on the original dense test core and its
completion. No fixed-aperture positivity, native Lean theorem or RH
is asserted by this note.
