#!/usr/bin/env python3
"""CR-ALP-038: exact-rational combined original signed infinite tail bound.

This verifies consequences of source-specific analytical Gamma and seven-prime
coefficient formulas for two *explicit rational interior pole-null sources*.
The previous CR-ALP-034 12-dimensional positive-form margin and CR-ALP-027
15,000-mode positive tail are analytical premises, locally interval-audited
separately; this is NOT a Lean formalization, a 15k low-section proof, or RH.
"""
from fractions import Fraction as F
from math import isqrt

A=F(11,10)
N=15000
S=F(177,50) # original seven authentic Lambda(n)/sqrt(n), strict upper

class I:
    def __init__(self,x,y=None):
        if isinstance(x,I):self.lo,self.hi=x.lo,x.hi
        else:self.lo,self.hi=F(x), F(x if y is None else y)
        assert self.lo<=self.hi
    def __add__(self,v):
        v=I(v);return I(self.lo+v.lo,self.hi+v.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,v):return self+-I(v)
    def __rsub__(self,v):return I(v)+-self
    def __mul__(self,v):
        v=I(v);p=[a*b for a in (self.lo,self.hi) for b in (v.lo,v.hi)]
        return I(min(p),max(p))
    __rmul__=__mul__
    def inv(self):
        assert self.lo>0,self.lo
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,v):return self*I(v).inv()
    def __pow__(self,n):
        assert n in (2,3)
        if n==2:
            if self.lo<=0<=self.hi:return I(0,max(self.lo*self.lo,self.hi*self.hi))
            return I(min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
        return I(self.lo**3,self.hi**3)
    def absup(self):return max(abs(self.lo),abs(self.hi))

def atan_interval(z,num=32):
    assert num%2==0
    s=sum((F((-1)**m, (2*m+1)*z**(2*m+1)) for m in range(num)),F(0))
    extra=F(1,(2*num+1)*z**(2*num+1))
    return I(s,s+extra)
PI=16*atan_interval(5)-4*atan_interval(239)
assert 3<PI.lo<PI.hi<F(22,7)
nu=lambda j: PI*F(j,1)/(2*A)

def report(label,indices,rats):
    anchor=1 if indices[0]%2 else 2
    assert len(indices)==len(rats)==5 and all(j%2==anchor%2 for j in indices)
    x=[F(i,10**7) for i in rats]
    va=nu(anchor)
    ll={}
    for j in indices:
        v=nu(j)
        ll[j]=v*(va**2+F(1,4))/(va*(v**2+F(1,4)))
        assert ll[j].lo>0
    danchor = -sum((xi*ll[k] for k,xi in zip(indices,x)),I(0))
    d={**{k:I(xi) for k,xi in zip(indices,x)},anchor:danchor}
    D=sum((di*nu(j) for j,di in d.items()),I(0))
    C0=sum((di.absup() for di in d.values()),F(0))
    C2=sum((di.absup()*nu(j).hi**2 for j,di in d.items()),F(0))
    C3=sum((di.absup()*nu(j).hi**3 for j,di in d.items()),F(0))
    assert max(indices)<N
    # Uniform all-j pole-corrected signed PRIME, then GAMMA coefficient envelope:
    # |a_j-tc_j| <= Acoef/j^2+Bcoef/j^3, for matching parity j>N.
    # E_prime uses 1/(1-rho^2)<2 and nu_max/nu_j<1.
    # E_gamma uses sum_{m>=1}1/(nu_j^2+beta_m^2)<=pi/(4nu_j),
    # and sum_{m>=1}beta_m^(-2)<1/2; exp(-2a beta_m)<1/200.
    factor=2*A/3
    Acoef=(2*S/A+F(201,200)*PI.hi/(2*A))*D.absup()*factor**2
    Bcoef=(S*C0/(2*A)+8*S*C2/A+F(201,200)*C3/A)*factor**3
    floor4=isqrt(3*N**3)
    floor6=isqrt(5*N**5)
    assert floor4>0 and floor6>0
    bound=Acoef/floor4+Bcoef/floor6
    # The true physical Gram on interior source is I+lambda lambda^T,
    # so ||sum x_k w_k||^2 >= sum x_k^2 exactly.
    physical_mass_floor=sum((xi*xi for xi in x),F(0))
    # Previous certified twelve-mode lower spectral margin 3/10^12 and
    # infinite tail margin 1/4:
    # W(p,p) >= (3/10^12)||p||^2, 4||Q_NTp||^2 <=4 bound^2.
    assert F(3,10**12)*physical_mass_floor>4*bound**2, (label,float(bound))
    # Young with epsilon=1/8 gives a uniform quantitative mixed-tail margin:
    # W(alpha*p+q)>= [m_p-8*c²]*|alpha|²+(1/8)||q||².
    assert F(3,10**12)*physical_mass_floor-8*bound**2 > F(28,10**13)
    print('PASS',label)
    print(' rational vector:',[(k,str(v)) for k,v in zip(indices,x)])
    print(' D interval:',float(D.lo),float(D.hi),'(small = retained cancellation)')
    print(' C0,C2,C3 upper:',*[float(z) for z in (C0,C2,C3)])
    print(' Acoef,Bcoef:',float(Acoef),float(Bcoef))
    print(' full infinite signed residual upper:',float(bound))
    print(' low W lower:',float(F(3,10**12)*physical_mass_floor),
          '4 residual squared upper:',float(4*bound**2))
    print(' PASS quantitative mixed-tail Young lower >2.8e-12 on line coefficient, 1/8 on tail norm')
    # Export a manageable exact *outward* decimal rational rather than a
    # multi-thousand-digit unevaluated Fraction numerator/denominator.
    ceil=lambda v,scale:(v.numerator*scale+v.denominator-1)//v.denominator
    floor=lambda v,scale:v.numerator*scale//v.denominator
    export_scale=10**15
    return {'label':label,'indices':indices,'rational_num':rats,
            'D_lower_1e15':str(floor(D.lo,export_scale))+'/10^15',
            'D_upper_1e15':str(ceil(D.hi,export_scale))+'/10^15',
            'tail_norm_upper_1e15':str(ceil(bound,export_scale))+'/10^15',
            'low_weil_lower_1e25':str(floor(F(3,10**12)*physical_mass_floor,10**25))+'/10^25',
            'four_tail_sq_upper_1e25':str(ceil(4*bound**2,10**25))+'/10^25'}

def exact_coefficient_algebra_regression():
    """Independent Fraction identities in rational toy frequencies; not Lean."""
    va, vk, vj, beta=F(3,2),F(9,2),F(401,2),F(5,2)
    d3=F(1)
    d1=-(vk/(vk*vk+F(1,4)))/(va/(va*va+F(1,4)))
    d={1:(va,d1),3:(vk,d3)}
    assert sum((c*v/(v*v+F(1,4)) for v,c in d.values()),F(0))==0
    D=sum((c*v for v,c in d.values()),F(0))
    gam=sum((c*v/(beta*beta+v*v) for v,c in d.values()),F(0))
    assert gam==D/(beta*beta)-sum((c*v**3/(beta*beta*(beta*beta+v*v)) for v,c in d.values()),F(0))
    gam0=-2*vj/A*gam/(vj*vj+beta*beta)
    pole=-2*vj/A*gam/(vj*vj+F(1,4))
    corrected=(2*vj)/(A*(vj*vj+F(1,4)))*gam*(beta*beta-F(1,4))/(vj*vj+beta*beta)
    assert gam0-pole==corrected
    sj=F(-3,7)
    si={1:F(2,5),3:F(-4,9)}
    for k,(v,c) in d.items():
        exact=2/A*(vj*si[k]-v*sj)/(vj*vj-v*v)
        expanded=(2/A)*(si[k]/vj-v*sj/(vj*vj))+(2/A)*(vj*si[k]-v*sj)*v*v/(vj*vj*(vj*vj-v*v))
        assert exact==expanded
    print('PASS independent Fraction regressions: pole-null Gamma resolvent and paired-shift coefficient decomposition')

if __name__=='__main__':
    import json
    exact_coefficient_algebra_regression()
    # Rational coefficients copied from two independent original Hmid/M
    # generalized-eigenvalue numerical diagnostics, NOT theorem premises.
    odd=report('odd [3,5,7,9,11]',[3,5,7,9,11],
               [6872105,-6735272,2338391,-188293,-20112])
    even=report('even [4,6,8,10,12]',[4,6,8,10,12],
               [-7021868,6604165,-1934708,-26547,43236])
    result={'stage':'CR-ALP-038','status':'exact rational bounds for two original signed one-dimensional cross-tail Schur tests; no 12-mode block/full 15k positivity',
            'original_support_radius':'11/10','finite_tail_cutoff':N,
            'premises':'CR-ALP-027 tail W(q,q)>=1/4||q||²; CR-ALP-034 W(p,p)>=3e-12||p||² on E12; original seven negative primes and digamma coefficient identities',
            'parity_tests':[odd,even],
            'not_claimed':'original full E12-plus-tail Schur; intermediate 12..15000 signed shell; Lean; RH'}
    from pathlib import Path
    out=Path(__file__).with_name('CR-ALP-038_exact_two_line_schur_results.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print('CR-ALP-038 EXACT RATIONAL TWO-LINE INFINITE-TAIL SCHUR AUDIT PASS')
