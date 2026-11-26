# RETURN098-A — Explicit formula and mixed functional receiving proof

Scope: real tests supported in I=(-5/4,5/4), not arbitrary support.
This is an analytic derivation (not a native Lean certificate), using
the standard Guinand–Weil explicit formula as an external theorem.
Reference: Proposition 5 of "Bounding the log-derivative of the
zeta-function", Math. Z. (2021), DOI 10.1007/s00209-021-02820-9.
https://link.springer.com/article/10.1007/s00209-021-02820-9
The proposition, unlike several applications in that article, does
not assume RH.

## 1. Literal spectral object

For real f in C_c^infty(I), set F(z)=integral f(x) exp(i z x) dx,
h(z)=F(z)F(-z), and c(t)=integral f(x)f(x+t) dx.
The entire even h decays faster than any inverse power on bounded
horizontal strips, by integration by parts. It is admissible in the
cited explicit formula. For the reference's 2pi Fourier convention,
hat h(y)=2pi c(2pi y). On the real axis h=|F|².

The explicit formula therefore gives

    sum_rho h((rho-1/2)/i)
      = (1/(2pi)) integral |F(r)|² g(r) dr
        -2 sum_n Lambda(n)/sqrt(n) c(log n)
        +2 F(i/2)F(-i/2),

with g(r)=Re psi(1/4+ir/2)-log(pi).
No assumption that rho lies on the critical line is made.
In particular h at a complex zero is not an absolute square.

Put beta_m=2m+1/2 and g0=psi(1/4)-log(pi).
The convergent digamma series gives

    g(r)-g0 = sum_{m>=0} 2r²/[beta_m(beta_m²+r²)].

The Fourier transform (angular convention) of exp(-beta|x|)
is 2beta/(beta²+r²). Parseval and monotone convergence of the
nonnegative increments yield exactly the physical Gamma sum of
RETURN096-B. The prime list is exactly the eight prime powers
with log(n)<5/2; no higher prime-power contribution is omitted.
The pole term is twice the product of the two exponential
moments. Thus the explicit formula has the physical signs and
normalization of RETURN096-B.

## 2. Literal exponential bilinear entry

For complex exponents z,w, write
J(v)=integral_{-a}^a exp(vx) dx, a=5/4.
For beta>0 and beta±z nonzero, direct integration gives

    integral_I exp(-beta|x-y|) exp(zx) dx
      = 2beta exp(zy)/(beta²-z²)
        - exp(-za) exp(-beta(y+a))/(beta+z)
        - exp(za) exp(-beta(a-y))/(beta-z).

Integrating against exp(wy) shows that the m>=1 Gamma correction
to the bulk term consists of

    exp(-a(z+w)) /[(beta+z)(beta-w)]
    + exp(a(z+w)) /[(beta-z)(beta+w)]
    - exp(-beta D) exp(a(w-z)) /[(beta+z)(beta-w)]
    - exp(-beta D) exp(a(z-w)) /[(beta-z)(beta+w)].

Here D=2a. These series are absolutely convergent for the
actual retained and grid exponents (real parts in [-1/2,1/2]),
since beta>=5/2. The bulk increment sum is
(psi((5/2+z)/2)+psi((5/2-z)/2))/2-psi(5/4), multiplied by J(z+w).
The first two boundary sums equal the divided difference of
psi((5/2+v)/2) divided by 2 times the exponent difference;
the coinciding-node value is psi'((5/2+v)/2)/4.
The last two sums are the direct exponentially weighted products.
This is exactly RETURN048's gpart, ss, dec expression.

For m=0, the integral over x<=y is
[J(z+w)-exp(-a(w+1/2)) J(z-1/2)]/(w+1/2)
after naming the outer exponent z and integrating the inner
exponent w over the left triangle. At w=-1/2 take the limit,
J'(z-1/2)+a J(z-1/2). The other triangle interchanges z,w.
The diagonal has measure zero. The mass coefficient is 4.
This is the separate triangle term in RETURN048; no division
by a resonant denominator is performed.

For prime displacement s=log(n), changing variables to the
midpoint of the overlap gives the two translations jointly:

    [exp((w-z)s/2)+exp(-(w-z)s/2)]
       integral_{-a+s/2}^{a-s/2} exp((z+w)t) dt.

Multiply by -Lambda(n)/sqrt(n). The bilinear pole is
J(z+1/2)J(w-1/2)+J(z-1/2)J(w+1/2).
Consequently the complete authenticated RETURN048 raw_entry
formula is the symmetric complex-bilinear extension of the
physical form. This is not a sesquilinear matrix identity.
The real cos/sin/sinh combinations recover its real restriction.

## 3. Domain and mixed functional continuity

The physical Gamma form differs from its period-five Fourier
form by a bounded image operator (RETURN096-B). The prime and
pole operators are also bounded on L2(I). The real-line Gamma
multiplier has g(r)-g0 comparable, up to positive constants
and additive constants, to log(2+|r|), by splitting the
positive digamma series at beta≈|r| and integral bounds.
The identical argument applies on the circle. The image
identity, first on finite Gamma sums and then by the bounded
image limit, identifies the two log-weighted form domains on
zero-extended functions in I and gives equivalent shifted
form norms.

The zero extensions of retained basis functions are piecewise
smooth. Their circle Fourier coefficients are O(1/|k|), so
multiplying by the logarithmic Gamma symbol still produces
an l2 sequence. Thus for fixed retained p, the periodic
Gamma pairing with q extends to a bounded L2 functional,
and so does W(p,q) by the bounded corrections. In particular,
for any q in the form domain, expand q in the rephased
even-grid L2 basis of its parity. Continuity in L2 proves

    W(p,q) = sum_{k in complementary grid} a_k W(p,phi_k_rephased).

No assertion that the partial sums lie in P-perp is needed.
No assertion that their zero extensions converge in the
logarithmic form norm is needed. The coefficient sequence of
this functional is square summable (Parseval), so the series
converges by Cauchy–Schwarz. The actual constraint kernel is
imposed on the full coefficient vector, not on each partial
sum. This resolves the potential domain mismatch in taking
the mixed source as a dual on the comparator energy space.

The smooth compactly supported core is dense in the
zero-supported log form domain: contract support by the
L2-unitary real-line dilation f_lambda(x)=lambda^(-1/2)
f(x/lambda), lambda<1. In angular Fourier variables its
action is lambda^(1/2) F(lambda r). It is uniformly bounded
on the log-weighted L2 space for lambda near 1, because the
weights log(2+|r|) and log(2+|r|/lambda) are uniformly
comparable. Strong continuity follows by approximating
with compactly supported smooth functions in Fourier space,
on which it follows from dominated convergence. The
contracted support has a positive distance from the
endpoints. Convolve it with a smooth compactly supported
approximate identity of smaller support; Fourier multipliers
are uniformly bounded and converge pointwise to 1.
Dominated convergence yields form-norm convergence. The
mollifier can be chosen real and even to preserve parity.
This proves density in the physical form norm, not density
in the stronger domain of the Gamma operator.

## 4. What this does and does not close

This supplies the analytic attachment of the standard
explicit-formula quadratic object to the RETURN048
realization, including its resonant term, and the L2
receiving identity for mixed source sequences. Numerical
entry containment and the grouped/gauged model identities
still have their own obligations: RETURN097-A/B/C must
complete, or any alternative certified enclosures must be
transported through the terminal comparison. The 095-A
error estimates must be attached to the identical model
and norm; this note does not recheck every such estimate.
The existing conditional Schur floors are not promoted
to an unconditional result here.

No native Lean proof and no arbitrary-support theorem.
In particular fixed-diameter positivity, even if finally
certified, cannot be extended by scaling: the arithmetic
translations and Gamma multiplier are not dilation invariant.
