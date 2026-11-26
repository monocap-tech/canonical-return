# RETURN099-A — Scalar model remainder receiving derivation

Scope: analytic closure of the scalar inequalities used in RETURN085,
not a new numerical error audit or an unconditional positivity theorem.
RETURN097-B/C remains queued. The actual RETURN070 leading evaluator
uses the *true* Gamma drift, enclosed by a log model plus drift reserve;
the infinite odd source is not silently replaced by a pure log model.

## 1. Gamma drift

Use the full multiplier g(r)=Re psi(1/4+ir/2)-log(pi),
beta_m=2m+1/2, and g0=g(0). The convergent expansion is

    g(r)-g0 = sum_{m>=0} 2r²/[beta_m(beta_m²+r²)].

On a compact positive r interval the derivative series converges
uniformly by comparison with sum beta_m^-3. Its summand is

    F_r(beta)=4r beta/(beta²+r²)².

For an absolutely continuous integrable function F on [0,infty)
of integrable derivative, a mesh of size 2 with sample 2m+1/2
in cell [2m,2m+2] satisfies

    |sum_m F(2m+1/2) - (1/2) integral_0^infty F|
       <= integral_0^infty |F'|.

Indeed the difference in each cell is the average of the
sample-minus-function difference, bounded by the variation in
that cell. Apply this to F_r. Its integral is 2/r, it increases
up to beta=r/sqrt(3) and decreases to zero, and its maximum is
9/(4 sqrt(3) r²). Therefore

    |g'(r)-1/r| <= 9/(2 sqrt(3) r²) < 3/r².

Integration on [R,r], r>=R>0, proves the inherited (weaker)
bound |g(r)-g(R)-log(r/R)| <= 8/R. In particular the true
odd leading scalar is l=(2Pc-g(r)+g(R))/r. The moment and Abel
bounds of RETURN090-B use log(r/R) *with* the 8/R error in
both the scalar square and the auxiliary crosses.

## 2. Resolvent scalar residuals

Set h=sum_{m>=1} 1/(beta_m²+r²),
e=sum_{m>=1} exp(-beta_m D)/(beta_m²+r²),
f=sum_{m>=1} beta_m exp(-beta_m D)/(beta_m²+r²),
D=5/2. Monotonicity of the integrand gives

    (1/2) integral_2.5^infty (beta²+r²)^-1 d beta <= h
      <= (1/2) integral_0.5^infty (beta²+r²)^-1 d beta.

The half-density integral from zero is pi/(4r), and the missing
lower integral is at most (5/4)/r². Thus the inherited
|r² h-pi r/4|<=2 holds. Write t=exp(-5). Exact geometric sums
give

    r² e <= exp(-25/4)/(1-t),
    r² f <= exp(-25/4)[(5/2)/(1-t)+2t/(1-t)²].

Both are <1: exp(5)>4 and exp(25/4)>16 by positive Taylor
partial sums (even 1+x+x²/2 suffices for the latter).
So |r²(h-e)-pi r/4|<=3.

Define the exact scalar factors

    alpha(r)=2r(h-e)+2Ps,
    d(r)=g(r)+g0+4-1/(r²+1/4)+2f-2Pc.

They follow by collecting RETURN098-C's four factors, with
g excluding no m=0 term in this note. At cutoff R set
alpha0=pi/2+2Ps, D0=g(R)+g0+4, x=r/R=k/K,
d0=D0+log(x)-2Pc. Then

    |alpha-alpha0| <= 6/r,
    |d-d0| <= 8/R+3/r².

The last inequality uses the triangle estimate on
-1/(r²+1/4)+2f (the weaker 3/r² is sufficient).

## 3. Vector resolvent expansion

For each raw retained mode j, write ell_j for the eigenvalue
of its interior second derivative. With fixed chart C/T,
define B(x),V(x) by applying that chart to raw b_j and v_j
weighted by (1+ell_j/r²)^-1. Define L(x) by applying the
same chart to c_j/(1+ell_j/r²)+d_j/[2(1+1/(4r²))].
The separate pole denominator is necessary, even though it
equals the retained denominator at the sinh mode.

For lambda_j=-ell_j/R² and lambda_p=-1/(4R²) the exact identity

    1/(1-lambda x^-2)
      = sum_{n=0}^{N-1} lambda^n x^-2n
          + lambda^N x^-2N/(1-lambda x^-2)

holds for x>=1 and |lambda|<1. It gives the RETURN066
vector moments directly, including the d_j/2 contribution
to L. The rational chart acts on the sums before taking
Euclidean norms. This is the source of the moment norms
and geometric remainder bound in RETURN085, not a
componentwise claim about the retained energy norm.

The exact normalized-chart mixed source is

    t_even = alpha B(x)/r + d V(x)/r² + 2 L(x)/r²
    t_odd  = -d B(x)/r + alpha V(x)/r² + 2 L(x)/r.

Replace the factors by alpha0,d0 and vectors by their first
N moments for either finite model. Let B,V,L now denote
bounds on the full vector norms and b_t,v_t,l_t on omitted
moment sums at x=1. With u=4/(5K)>1/R,
delta=8u+3u², and |D0|<10, 2|Ps|,2|Pc|<9,
log x<3 on the two finite bands, the exact triangle
inequality gives

    E_even <= 6B u² + delta V u² + 11 b_t u
                           + 23 v_t u² + 2 l_t u²
    E_odd  <= delta B u + 6V u³ + 23 b_t u
                           + 11 v_t u² + 2 l_t u.

These are the RETURN085 budgets, not the older weaker 42
coefficient in some historical generation scripts.

## 4. Infinite leading tail with actual odd compensation

For x>=1, |d|<=20+log x at K=262144 (use |D0|<10,
2|Pc|<9, and 8/R+3/R²<1). In even parity the model
is alpha0 B0/r. The first-omitted moments have an extra
x^-2, so

    E_even(x) <= 11 b_t u x^-3
       + (6B+20V+2L)u² x^-2 + V u² x^-2 log x.

For odd parity let G=D0 B0-2L0 and let
l=(2Pc-g(r)+g(R))/r, using the full g.
The archived compensation is -G r/(r²+1/4), not -G/r.
Subtract the actual model l B0-G r/(r²+1/4)
from the exact odd source. Algebra gives

    t_odd-m_odd =
       -d(B(x)-B0)/r + 2(L(x)-L0)/r
       + alpha V(x)/r²
       - e0 B0/r + G[r/(r²+1/4)-1/r],

where e0=-1/(r²+1/4)+2f. The compensation sign follows
from RETURN070 and RETURN090-B, which preserve the
negative mixed Gram cross. Because |e0|<=3/r² and
|r/(r²+1/4)-1/r|<=1/(4r³), this yields

    E_odd(x) <= (20 b_t+2 l_t)u x^-3
       + b_t u x^-3 log x + 11 V u² x^-2
       + (3||B0||+6V+||G||/4)u³ x^-3.

The terms are the inherited RETURN085 infinite-tail
envelope. Importantly no independent log-drift error is
charged again in this actual-source envelope: the model
uses the true drift, and the scalar contraction enclosure
already pays the approximation of that drift.

## 5. Sum and receiving norm

For p>1 and integer q=0,1,2, the half-density integral
of x^-p log^q x on [1,infty), at stride 2/K, is

    K q!/[2 (p-1)^(q+1)].

The difference from the positive-grid sum with x>1 is
bounded by total variation, as in Section 1 using the
right endpoints. For p>=2 and q=0,1,2 the maximum is <=q!
(by optimizing in log x), so 2q! is a safe reserve.
Expanding each squared positive envelope termwise gives
exactly the rational integrals and grid reserves used
in RETURN085. RETURN098-B proves the subsequent inverse
receiving estimate with the negative-update trace; no
band orthogonality is assumed. RETURN095-A composes the
three separately directed norm caps by triangle.

## 6. Boundary of this derivation

This attaches the scalar inequalities analytically to
the physical grouped formula and the actual historical
model. It does not recompute all numerical input norms,
the fixed D0 bounds, the negative trace receipts or the
infinite contraction enclosures. The earlier arithmetic
audits supply those values, subject to final joint source
identity/custody verification. In particular there is no
new PASS declaration, native Lean theorem, or
arbitrary-support conclusion in this note.
