#!/usr/bin/env python3
"""CR-ALP-002: numerical regressions for exact ramp identities.

Not a Lean proof, an operator-domain certificate, or an RH proof.
Requires mpmath and sympy.
"""
import mpmath as mp
from sympy import primerange

mp.mp.dps = 55

def ramp(y, X):
    return max(mp.mpf(0), X*mp.exp(-abs(y)/2)-mp.exp(abs(y)/2))

def D(X):
    z = 1/X
    return 1-X*mp.atanh(z)-mp.log1p(-z*z)/2

def gamma_closed(X):
    return (-(mp.euler+mp.log(4*mp.pi))*(X-1)
            +2*mp.log((X+1)/2)+2*(X-1)*mp.atanh(1/X))

def gamma_integral(X):
    L = mp.log(X)
    near = mp.quad(lambda t: 2/(1+mp.exp(-t)), [0,L])
    far = mp.quad(lambda t: 2*(X-1)*mp.exp(-t)/(1-mp.exp(-2*t)), [L,mp.inf])
    return -(mp.euler+mp.log(4*mp.pi))*(X-1)+near+far

def weighted_prime_ramp(X):
    out = mp.mpf(0)
    for p in primerange(2,int(X)+1):
        n=int(p)
        while n<=X:
            out += mp.log(p)*(X/mp.mpf(n)-1)
            n *= int(p)
    return out

C = 2+mp.euler-mp.log(4*mp.pi)
for n in [2,3,10,100,1000]:
    X=mp.mpf(n)
    L=mp.log(X)
    # One compactly supported source f_X(u)=exp(u/2) 1_[0,L](u).
    for t in [mp.mpf(0), L/4, L/2, L]:
        corr=mp.quad(lambda u: mp.exp((u+t)/2)*mp.exp(u/2), [0,L-t])
        assert mp.almosteq(corr, ramp(t,X), rel_eps=mp.mpf('1e-40'))
    Fplus=X-1
    Fminus=L
    Kpole=(mp.quad(lambda t: ramp(t,X)*mp.exp(t/2), [0,L])
           +mp.quad(lambda t: ramp(t,X)*mp.exp(-t/2), [0,L]))
    assert mp.almosteq(Kpole,Fplus*Fminus,rel_eps=mp.mpf('1e-45'))
    gam=gamma_closed(X)
    assert mp.almosteq(gam,gamma_integral(X),rel_eps=mp.mpf('1e-45'))
    d=D(X)
    dseries=sum(X**(-2*j)/mp.mpf(2*j*(2*j+1)) for j in range(1,131))
    assert abs(d-dseries)<mp.mpf('1e-40')
    S=weighted_prime_ramp(X)
    G=S-X*L+(1+mp.euler)*X
    W=2*(X-1)*L-2*S+gam
    Wrenorm=C*(X+1)-2*G+2*mp.log(2*mp.pi)-2*d
    assert mp.almosteq(W,Wrenorm,rel_eps=mp.mpf('1e-44'))
    for tau in [mp.mpf(0),mp.mpf('0.3'),mp.mpf(7)]:
        F=2*mp.quad(lambda t: ramp(t,X)*mp.cos(tau*t),[0,L])
        Fc=(X+1-2*mp.sqrt(X)*mp.cos(tau*L))/(tau**2+mp.mpf('.25'))
        assert mp.almosteq(F,Fc,rel_eps=mp.mpf('1e-36'))
    print('PASS X={} pole={} G={} W={}'.format(n,mp.nstr(Kpole,12),mp.nstr(G,12),mp.nstr(W,12)))
print('CR-ALP-002 QUADRATURE/CROSS-IDENTITY PASS; no RH/Lean/domain certification')
