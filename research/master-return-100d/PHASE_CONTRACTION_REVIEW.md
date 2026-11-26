# RETURN100-D — Oscillatory feature contraction review

Scope: analytic attachment of the finite high-band Abel and infinite
leading-phase contractions. Reopened RETURN073_band.py, RETURN089
audit_high_band.py, RETURN090 audit_restricted.py and RETURN090-B
audit_asymptotic.py. This is a review of formulas, not an independent
execution of the archived RETURN086 direct-sum generator.

## 1. Exact phase algebra

Let k=2m, r=2*pi*k/5, x=k/K, and let w_n=Lambda(n)/sqrt(n)
for the eight prime powers n=2,3,4,5,7,8,9,11.
For a positive rational z define v_z=exp(4*pi*i*log(z)/5).
Then v_z^m=exp(i*r*log(z)). Write S=2 sum w_n sin(r log n)
and C=2 sum w_n cos(r log n). Their complex Laurent coefficients
at rational keys n and 1/n are respectively (-i w_n,+i w_n)
and (w_n,w_n). Products of phase polynomials use multiplication
of rational keys and summation of coefficients at equal keys.
This is precisely SP,CP and conv in RETURN073, not a list of
formally distinct prime-power phases. In particular the key 1
has exact mean coefficient 2 sum w_n² for S² and C², and zero
for SC. All nontrivial keys require a certified bound away from
v_z=1. Nontrivial rational key alone is not used as a quantitative
phase-gap certificate.

The ten scalar families in sums0, in literal order, are
weighted sums of 1,L,L²,S,C,LS,LC,S²,SC,C², with L=log x.
For a feature f=x^-p s(a+bL+cS+dC), multiplying two such
features and collecting coefficients gives exactly the ten-term
gram expression in RETURN073. The lookup array col of
RETURN089/090 is the symmetric multiplication table of the
four functions (1,L,S,C); its entries refer to the same ten
families. Consequently the model contraction is V^T M V
and the auxiliary cross is U0 V, with the retained moment
matrix V fixed at the actual anchor K=32768. The restricted
band 65536<k<=262144 changes summation endpoints, not the
anchor, features, D0 or retained moments.

## 2. Finite Abel identity with endpoint indexing

Let l=K/2+1 and u=U/2, a_m=f(2m/K), and let v have
absolute value one and v!=1. The first difference identity is

    (1-v) sum_(m=l)^u a_m v^m
      = a_l v^l - a_u v^(u+1)
        + sum_(m=l+1)^u Delta a_m v^m,

where Delta a_m=a_m-a_(m-1).
Iterating d times gives

    sum_(m=l)^u a_m v^m
      = sum_(j=0)^(d-1)
          [Delta^j a_(l+j) v^(l+j)
           - Delta^j a_u v^(u+1)]/(1-v)^(j+1)
        + sum_(m=l+d)^u Delta^d a_m v^m/(1-v)^d.

All bands here contain more than d points. At the left endpoint
Delta^j a_(l+j)=sum_(t=0)^j (-1)^t binom(j,t) a_(l+j-t);
at the right replace l+j by u. These are exactly ends() and
phase() in RETURN073. The remainder estimate uses

    sum_(m=l+d)^u |Delta^d a_m|
      <= h^(d-1) integral_1^(U/K) |f^(d)(x)| dx, h=2/K.

Indeed the dth difference is the d-fold integral of f^(d)
over increments in [0,h]. For fixed first d-1 increments,
the remaining integration intervals for distinct m have
disjoint interiors, and lie inside the source interval.
Integrate the resulting bound over the first d-1 increments.
This proves the factor h^(d-1), rather than h^d, in oscill().

For f(x)=x^-p log^q x and q=0,1,2, use parameter
differentiation of (-1)^d (p)_d x^(-p-d).
Write H=sum_(j=0)^(d-1) 1/(p+j), H2=sum 1/(p+j)^2.
The first two derivatives of (p)_d divided by (p)_d are
H and H²-H2. Taking absolute values and bounding the latter
by H²+H2, the integral bounds over x>=1 are

    q=0: (p)_d/a
    q=1: (p)_d (1/a²+H/a)
    q=2: (p)_d (2/a³+2H/a²+(H²+H2)/a),

where a=p+d-1. These are the actual deriv_int expressions.
Multiplying by the absolute phase coefficients and dividing
by |1-v|^d gives the directed remainder radius. Taking the
real part at the end is valid because the Laurent polynomials
are real-valued on the unit circle.

For the nonoscillatory terms the trapezoidal endpoint
correction is (f(U/K)-f(1))/2 for the grid excluding 1
and including U/K. The trapezoidal error after division by h
is bounded by h/8 integral |f''|, via the cell kernel
t(h-t)/2. This is the 1/(4K) factor in RETURN073 I(p,q);
the finite endpoint integral is evaluated by integration by
parts. No approximation of the integer grid density is
discarded without this reserve.

## 3. Auxiliary resolvents and direct-sum cross-check

Each auxiliary row on this band is an explicit power of
k^-1 times a product of (1-t_j/k²)^-1. The retained anchor
K is unchanged under restriction. Uniform absolute
convergence on x>=1 gives a positive omitted coefficient
mass bound E; the tail after n terms is at most E*x^-2n.
RETURN089/090 use the positive moment at exponent p+2n to
pay this tail, with all exponents within the direct bank's
2..46 inventory. The uniform scalar multiplier bound
|a|+3|b|+9(|c|+|d|) follows from x<=8 and the prime
triangle bound. In the restoration rows the zeroth power
coefficient vanishes because each rational row sum is zero.
The remaining moment terms and omitted powers are contracted
with the same 24-feature V as the common auxiliary rows.

The archived RETURN086 direct_sums bank has the 114688
full-band and 98304 restricted-band coordinates, respectively.
RETURN089 and RETURN090 use that bank to reconstruct all
feature TT/U and odd restoration crosses without invoking the
historical Abel contraction. These are independent contractions
of frozen directed scalar sums, not a new independent
generation of those scalar sums. Both archived terminal
variants receive the same physical scalar definitions.
The new analytic review establishes the formulas by which
the Abel intervals themselves must be justified.

## 4. Infinite leading phase

For k>K use the same rational phase keys and stride-two
grid, with x=k/K. A geometric partial sum has absolute value
at most 2/|1-v|. For a nonnegative unimodal or decreasing
weight f(x)=x^-p log^q x tending to zero, summation by parts
bounds the oscillatory sum by 2 sup(f)/|1-v| if one uses
the signed first-difference identity above: its numerator
is bounded by f(first)+total variation <=2 sup(f).
For q=0 sup<=1; for q>0 optimize in log x to get
sup=(q/(e p))^q <=(q/p)^q.
This justifies RETURN090-B osc_cap for nonconstant keys;
the constant key is handled by the ordinary power/log sum.

The true odd drift differs from log x by at most 8/R,
as derived in RETURN099-A. Its squared contribution and
crosses with the oscillatory polynomial are enclosed by
the explicit triangle terms in RETURN090-B, including
2 delta (9 I(2,0)+I(2,1))+delta² I(2,0).
The tail compensation and leading Gram cross remain as in
RETURN098-B. This is not a model change to pure logarithmic
drift, and no second error reserve is silently removed.

## 5. Residual obligations

The phase and endpoint formulas have been reviewed against
the named source implementations. This does not certify
every numerical phase-gap interval, D0 at both higher anchors,
the RETURN086 direct-sum generator or the provenance of all
primitive constant intervals. Those inputs are still shared
dependencies. In particular the frozen direct-sum intervals
cannot by themselves prove the accuracy of their generator
from containment of downstream contractions.

Permanent evidence custody completed separately at commit
87c9cd17b29da9d0a57cd426e0aa76a24eea0f65, covering all twenty
accepted archive ZIPs in 23 repository parts. No new runner
activated here. No unconditional theorem, native Lean proof
or arbitrary-support conclusion is asserted.
