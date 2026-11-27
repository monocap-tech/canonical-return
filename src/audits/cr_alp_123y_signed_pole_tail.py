#!/usr/bin/env python3
"""Leading signed log/oscillatory pole envelope for the 123W-type receiver.

Inputs are original 046 interval namespace and actual trial intervals.
Output is an upper L2 bound on the leading expression plus its pole-column
remainder; caller must add the three centered sine-series charges separately.
No optimizer, no new continuum assumptions.
"""
from fractions import Fraction as F
from math import isqrt

def root(v):
    assert v>=0
    D=10**110;z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)

def leading(g,k,J,fb,T,H0):
    I,S,ZERO=(g[x] for x in ('I','S','ZERO'))
    rt,eh,em,L,qt=(g[x] for x in ('rt','eh','em','L','qt'))
    w,ns,ts=(g[x] for x in ('W','NS','TS'))
    b0=(eh+em)/2 if k==1 else -(eh-em)/2
    b1=-(eh-em)/4 if k==1 else (eh+em)/4
    ep,en=(b0+2*b1)/2,(b0-2*b1)/2
    D0=ep*L+en*(1-g['q'].sq())-2*b0+b1*g['S1']+b0/4*g['S2']
    EC=[-(-1)**k*ex*(b*b0-b1)/(bb-qt)
        for b,bb,ex in zip(g['B'][1:],g['B2'][1:],g['EX'][1:])]
    K=g['K'];p=F(2,5)**4
    # |EC_m| beta_m^2 <= 2 exp(-2a beta_m)(beta_m |b0|+|b1|).
    tail2=2*F(2,5)*p**K*((b0.up()*(2*K+F(1,2))+b1.up())/(1-p)
                                 +2*b0.up()*p/(1-p)**2)
    tail0=tail2/(2*K+F(1,2))**2
    D=D0-sum(EC,ZERO)+I(-tail0,tail0)
    E2=sum(e.up()*bb.up() for e,bb in zip(EC,g['B2'][1:]))+tail2
    bt=[ep*I(n).sqrt()+en/I(n).sqrt() for n in ns]
    J0=J-1
    assert J0>=2047
    # v=pi*j/(2a), hence log(v/(2pi))=log(j/(4a)).
    B=fb*b0/rt
    A=H0+B*(g['g0']+g['log'](I(F(J0,1))/(4*g['A'])))-2*fb*D/rt-2*T*sum((u*v for u,v in zip(w,bt)),ZERO)/rt
    C=[-2*T*b0*u/rt for u in w]
    th=[g['NUF']*t for t in ts]
    def denom(theta):
        si=g['sincos'](theta)[0]
        lower=max(F(0),F(si.l,S),-F(si.h,S))
        assert lower>0
        return lower
    # Sum over one parity. Integral of (A+B log(x/J0))^2/x^2
    # plus total variation quadrature reserve. A,B are signed intervals.
    mass=F((A.sq()+2*A*B+2*B.sq()).h,S)/(2*J0)
    a,b=A.up(),B.up()
    variation=(a*b+a*a + (b*b+2*a*b)/2 + b*b/2)/J0**2
    mass+=variation
    # Abel: |sum exp(i*j*theta) f(j)| <= TV(f)/|sin theta|
    # for a parity progression. TV(1/j^2), TV(log(j/J0)/j^2)
    # are both <= J0^-2. Retain the diagonal constant term separately.
    for i,ci in enumerate(C):
        u=ci.up()
        mass+=u*u/2*(F(1,2*J0)+F(1,J0**2)/denom(2*th[i]))
        mass+=2*u*(a+b)/(J0**2*denom(th[i]))
        for j in range(i):
            mass+=u*C[j].up()/J0**2*(1/denom(th[i]-th[j])+1/denom(th[i]+th[j]))
    assert mass>=0
    nu_lower=F(g['NUF'].l,S);assert nu_lower>F(15,11)
    lead=root(mass)/nu_lower
    # For v>= (15/11)J: |G-log(v/(2pi))| <= 4/v,
    # |R|<=1/v, |G|<=v, |G-g0|<=2v. Derivation in 123Y note.
    rr=F(rt.l,S);assert rr>1
    gamma2=4*b0.up()/rr+b0.up()/(2*rr)+2*(b1.up()+b0.up()/4)/rr
    gamma3=2*(D0.up()/4+E2+b1.up()+b0.up()/2)/rr
    prime2=2*b1.up()*sum(u.up() for u in w)/rr
    prime3=sum(u.up()*(v.up()+b0.up()) for u,v in zip(w,bt))/(2*rr)
    p2=fb.up()*gamma2+T.up()*prime2
    p3=fb.up()*gamma3+T.up()*prime3+H0.up()/4
    # Include only pole/leading-denominator error; other signed moments
    # remain with caller. Same parity integral bound with J0.
    remainder=sum(p*root(F(J0,2*(2*s-1)))/(nu_lower*J0)**s for s,p in ((2,p2),(3,p3)))
    return lead+remainder,{'constant':[str(F(A.l,S)),str(F(A.h,S))],
       'log_coefficient':[str(F(B.l,S)),str(F(B.h,S))],
       'cosine_coefficients':[[str(F(v.l,S)),str(F(v.h,S))] for v in C],
       'leading_norm_upper':str(lead),'pole_remainder_norm_upper':str(remainder),
       'pointwise_inverse_square_upper':str(p2),'pointwise_inverse_cube_upper':str(p3)}
