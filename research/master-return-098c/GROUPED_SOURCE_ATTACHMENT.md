# RETURN098-C — Literal grouped mixed-row algebra

This note derives the RETURN065/066/068 grouped source formula
from the physical realization in RETURN098-A. It is not a numerical
receipt or native Lean theorem. The exact RETURN065 evaluator was
reopened, not inferred from the names of the stored coefficient fields.

## 1. Setup

Put a=5/4, D=2a, beta_m=2m+1/2, r=2*pi*k/5 for a positive
even complementary grid index k, and sigma=(-1)^(k/2).
Write p for a single raw retained basis function, of parity
eta=+1 or -1, p''=ell p, b=p(a), v=p'(a).
Take q=cos(rx) for even parity and q=sin(rx) for odd parity.
The mixed entry used in the normalized chart is sigma W(p,q).
The denominator r²+ell is nonzero on the complementary grid:
retained trigonometric frequencies are excluded and the sinh
value ell=1/4 is positive.

Set M=integral_I p(x) exp(x/2) dx, J0=exp(-a/2) M and
sh=sinh(a/2). For beta>0 define J_beta by the same expression
with beta in place of 1/2. For beta²!=ell the differential identity
[(beta p-p')exp(beta x)]'=(beta²-ell)p exp(beta x) gives

    J_beta = [(beta b-v)-eta exp(-beta D)(beta b+v)]/(beta²-ell).

For beta=1/2 this quotient is NOT used in the resonant sinh case:
M=sinh(a)-a by integration, so J0 is finite.

Integration of (p' q-p q')'=(ell+r²)p q gives mass

    integral p q = 2 sigma v/(r²+ell)   (even),
                   -2 sigma r b/(r²+ell) (odd).

## 2. Gamma boundary rows

For each beta, direct convolution of the grid mode with
exp(-beta|x-y|) gives the full-line multiplier 2 beta/(beta²+r²)
times q, minus an endpoint correction. In the even case that
correction is sigma beta/(beta²+r²) times the sum
exp(-beta(a-y))+exp(-beta(a+y)); in the odd case it is
sigma r/(beta²+r²) times the difference
exp(-beta(a-y))-exp(-beta(a+y)).
Pairing with p and subtracting the convolution therefore adds
2 sigma beta J_beta/(beta²+r²) or 2 sigma r J_beta/(beta²+r²).
Together with the mass, the bulk multiplier is
gfull=g0+g+4r²/(r²+1/4), where g denotes the m>=1 positive
Gamma series, not the full gfull.

For m>=1 let h=sum 1/(beta²+r²), e=sum exp(-beta D)/(beta²+r²),
f=sum beta exp(-beta D)/(beta²+r²). Define gj,hj,ej,fj at the
retained resolvent by replacing r² by -ell. In particular
gj=sum -2ell/[beta(beta²-ell)] and these four sums converge
absolutely, including the retained sinh case. Partial fractions give

    s0=(hj-h)/(r²+ell)
    s1=(g-gj)/(2(r²+ell))
    s2=(r² h+ell hj)/(r²+ell)
    e0=(ej-e)/(r²+ell)
    e1=(fj-f)/(r²+ell)
    e2=(r² e+ell ej)/(r²+ell).

They are respectively the sums of 1,beta,beta² divided by
(beta²+r²)(beta²-ell), with the e variants weighted by
exp(-beta D). For s1, subtract the two convergent bulk Gamma
increments; no divergent unregularized sum of beta/(beta²+r²)
is used. The beta=1/2 endpoint contribution is handled separately.

The resulting Gamma boundary terms before multiplication by sigma are

    even: 2 sigma [b(s2-e2)-v(s1+e1)+(1/2)J0/(r²+1/4)]
    odd:  2 sigma r[b(s1+e1)-v(s0-e0)+J0/(r²+1/4)].

## 3. Prime and pole rows

For a shift s=log(n) and weight w=Lambda(n)/sqrt(n),
use the actual overlap x in [-a,a-s] for p(x+s)q(x).
The derivative identity
[p'(x+s)q(x)-p(x+s)q'(x)]'=(ell+r²)p(x+s)q(x)
and reflection (same parity) give both translations together.
Let u=p(a-s), u'=p'(a-s). Substituting the grid endpoint
values gives the signed prime contributions

    even: 2 sigma sum w[b r sin(rs)-v cos(rs)-u']/(r²+ell)
    odd:  2 sigma sum w[v sin(rs)+r b cos(rs)+r u]/(r²+ell).

The bilinear pole moment is twice M times the grid moment
in even parity, and minus twice M times that moment in odd
parity. The grid moments are 2 sigma (1/2)sh/(r²+1/4) even
and -2 sigma r sh/(r²+1/4) odd. Thus the pole contributions
are 2 sigma M sh/(r²+1/4) and 4 sigma M r sh/(r²+1/4).

## 4. Grouped source identity

Define Ps=sum w sin(rs), Pc=sum w cos(rs). The fixed raw
coefficients (c,d) are

    c_even = ell b(hj-ej)+v(gj/2-fj)-sum w u'
    d_even = J0+2 M sh = M exp(a/2)
    c_odd  = b(-gj/2+fj)+v(-hj+ej)+sum w u
    d_odd  = 2 J0+4 M sh = 2 M exp(a/2).

Then the exact rephased entries are

    even:
    2[b(r²(h-e)+r Ps)+v(gfull-g/2+f-Pc)+c_even]/(r²+ell)
       + d_even/(r²+1/4)

    odd:
    2r[b(-gfull+g/2-f+Pc)+v(h-e+Ps/r)+c_odd]/(r²+ell)
       + r d_odd/(r²+1/4).

All sigma factors cancel because sigma²=1. These formulas
are obtained by collecting the RETURN065 row expressions
literally, and match RETURN068's four shared factors.
RETURN085 previously checked 870 raw coefficient enclosures;
RETURN097-A is the new exhaustive shared scalar audit.
The retained congruence C/T acts linearly on these complete
raw entries, including the pole term. No independent treatment
of divergent Gamma and pole subexpressions is necessary.

The d entry is the last raw coefficient; it is not divided
by r²+ell. At the sinh node ell=1/4 the two denominators
coincide but the formula remains finite. The later moment
expansion uses the exact geometric resolvent identity for
each denominator independently.

## 5. Scope

This removes the algebraic receiving gap between the
separated evaluator and the physical mixed form for all
nonresonant grid indices, with the beta=1/2 retained
resonance explicitly handled. It does not verify the
literal archived raw coefficient generation code or the
complete numerical grouped row containment by itself.
RETURN097-A/B/C and final model/error attachment remain
required. No new runner is activated.
