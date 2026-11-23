#!/usr/bin/env python3
"""CR-INT-T1D-Y: source-preserving local selected-to-satellite Schur audit.

New checks: exact original 24-prime-power selected/satellite collision census;
rigorous physical |q(u_sigma,z_sigma)| bounds with prime and full Gamma;
strict Schur-corrected 2x2 determinants in each parity, with no T1D-I import.

Imports: prior T1C selected mu>7, T1D-W original Gamma self bounds and
central source cross bounds; T1D-X exact prime and Gamma near/far bounds.
The script checks finite arithmetic and determinants. It does NOT certify
whole-complement positivity, a compact top-eigenvalue upper bound, or RH.
"""
from fractions import Fraction as F
from math import isqrt

class Iv:
    def __init__(self,lo,hi=None):
        self.lo=F(lo); self.hi=F(lo if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self,x):
        x=x if isinstance(x,Iv) else Iv(x)
        return Iv(self.lo+x.lo,self.hi+x.hi)
    __radd__=__add__
    def __neg__(self):return Iv(-self.hi,-self.lo)
    def __sub__(self,x):return self+-x if isinstance(x,Iv) else self+Iv(-x)
    def __rsub__(self,x):return Iv(x)-self
    def __mul__(self,x):
        x=x if isinstance(x,Iv) else Iv(x)
        a=[p*q for p in (self.lo,self.hi) for q in (x.lo,x.hi)]
        return Iv(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=x if isinstance(x,Iv) else Iv(x)
        assert x.lo>0
        return self*Iv(1/x.hi,1/x.lo)
    def __rtruediv__(self,x):return Iv(x)/self

def lograt(q,N=90):
    q=F(q);assert q>=1
    a=(q-1)/(q+1);v=a;total=F(0)
    for k in range(N):
        total+=2*v/F(2*k+1);v*=a*a
    return Iv(total,total+2*v/(F(2*N+1)*(1-a*a)))

L2=lograt(2,100)
def logint(n):
    k=n.bit_length()-1
    return k*L2+lograt(F(n,2**k),100)

def root(n,scale=10**22):
    n=F(n);assert n>=0
    k=isqrt((n.numerator*scale*scale)//n.denominator)
    assert F(k*k,scale*scale)<=n<F((k+1)**2,scale*scale)
    return Iv(F(k,scale),F(k+1,scale))

def prime_base(n):
    for p in range(2,n+1):
        if n%p:continue
        if any(p%d==0 for d in range(2,isqrt(p)+1)):continue
        v=n
        while v%p==0:v//=p
        if v==1:return p
    return None

prime_powers={n:prime_base(n) for n in range(2,55) if prime_base(n)}
assert list(prime_powers)==[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53]
kap={n:logint(p)/root(n) for n,p in prime_powers.items()}
ln3=logint(3)
eps=F(1,10**7)
ell=F(1,10000)
width=F(1,4)+eps+ell/2
# For original cell j at (j-3/2)log3, signed prime shift, and
# satellite centered at wing*log3, require strict support overlap.
events=[]
for j in range(4):
    c=ln3*F(2*j-3,2)
    for n in prime_powers:
        for sign in (-1,1):
            x=c+logint(n)*sign
            for wing in (-1,1):
                displacement=x-ln3*wing
                if displacement.lo>-width and displacement.hi<width:
                    events.append((j,n,sign,wing))
                else:
                    assert displacement.lo>width or displacement.hi<-width, (j,n,sign,wing)

expected=[
    (0,2,1,-1),(0,13,1,1),(0,16,1,1),(0,17,1,1),(0,19,1,1),
    (1,2,-1,-1),(1,5,1,1),
    (2,2,1,1),(2,5,-1,-1),
    (3,2,-1,1),(3,13,-1,-1),(3,16,-1,-1),(3,17,-1,-1),(3,19,-1,-1),
]
assert events==expected, (events,expected)
assert {n for _,n,_,_ in events}=={2,5,13,16,17,19}
print('PASS: original 24-prime-power selected-to-satellite census, 384 strict comparisons, 14 directed overlaps')
print('CERTIFIED active original prime powers: 2,5,13,16,17,19; complete events:',events)

# Both translated satellite copies have mutually disjoint support,
# and ||(T_d y + T_-d y)/sqrt(2)|| = ||y||.
# Sobolev: ||y||_infty <= sqrt(2||y||||y_prime||).
# Therefore for unit z, ||z||_infty <= sqrt(r_sigma).
# Narrow selected bump phi has ||phi||_1 <= sqrt(ell ||phi||_2^2).
# Coefficients of u_+ / u_- on the four original selected cells are
# [1,1-s,1-s,1] / [1,-1-s,1+s,-1].
# Gamma cross: disjoint supports with separation > 29/100, w(t)<2.
# ||z||_1 <= sqrt(1+4eps) ||z||, and sum of |c| =2s or 2(s+2).

def exp_neg(z,N=102):
    z=F(z);assert 0<z<1
    v=F(1);total=F(0)
    for k in range(N):
        total+=v;v*=(-z)/F(k+1)
    assert N%2==0 and v>0
    return Iv(total,total+v)

gap=ln3/2-width
assert gap.lo>F(29,100)
ww=exp_neg(F(29,200))/(Iv(1)-exp_neg(F(29,50)))
assert ww.hi<2
print('PASS: exact selected/satellite support gap > 29/100; original Gamma off-diagonal kernel < 2')

s=4/root(3)
Nplus=(44-16*root(3))/3
Nminus=(44+16*root(3))/3
for parity in (+1,-1):
    label='even' if parity==1 else 'odd'
    r=38 if parity==1 else 46 # conservative source-wide T1D-V W bounds
    N=Nplus if parity==1 else Nminus
    assert N.lo>0
    c=[Iv(1), s-Iv(1),s-Iv(1),Iv(1)] if parity==1 else [Iv(1),s+1,s+1,Iv(1)]
    W=sum((kap[n]*c[j] for j,n,_,_ in events),Iv(0))
    # explicit upper is squared-root-rationally rounded OUTWARD
    factor=root(ell).hi/root(N.lo).lo
    p_bound=factor*root(r).hi*W.hi
    gamma_bound=2*factor*(2*s if parity==1 else 2*(s+2)).hi*root(1+4*eps).hi
    total=p_bound+gamma_bound
    target=F(27,100) if parity==1 else F(23,100)
    assert total<target,(label,float(total),target)
    print(f'CERTIFIED {label}: prime selected/satellite |cross| < {float(p_bound):.10f} ; Gamma < {float(gamma_bound):.10f}; combined < {target}')
    # T1D-W: central selected coupling < 2/15 even, 3/25 odd;
    # T1C: mu_sigma >7. New satellite R_zz and R_yz bounds
    central=F(2,15) if parity==1 else F(3,25)
    rzz=target**2/7
    ryz=central*target/7
    if parity==1:
        # Original signed Weil q two-mode lower bounds imported T1D-X
        qyy_lo=F(89,100)
        qzz_lo=F(89,100)+F(27,100)-F(1,1000)
        cross_upper_factor=F(68,100)+F(1,100)
        ryy=F(4,1575)
        certificate=F(1,20)
    else:
        qyy_lo=F(118,100)
        qzz_lo=F(118,100)-F(11,100)-F(1,1000)
        cross_upper_factor=F(71,100)+F(1,100)
        ryy=F(9,4375)
        certificate=F(1,5)
    # sqrt(2)<283/200, elementary exact square check.
    assert F(283,200)**2>2
    corrected_diag_y=qyy_lo-ryy
    corrected_diag_z=qzz_lo-rzz
    corrected_cross=(F(283,200)*cross_upper_factor)+ryz
    det_lower=corrected_diag_y*corrected_diag_z-corrected_cross**2
    assert corrected_diag_y>0 and corrected_diag_z>0 and det_lower>certificate,(label,det_lower,certificate)
    print(f'CERTIFIED {label}: Schur R[unit satellite] < {rzz}; |R[central,satellite]| < {ryz}')
    print(f'CERTIFIED {label}: full (q-R) 2x2 determinant > {certificate} ; rational margin {float(det_lower):.9f}')
print('PASS: both authentic central+satellite parity Schur blocks strictly positive WITHOUT T1D-I')
print('ANALYTIC IMPORTS: T1C mu>7; T1D-W central q/Gamma/R bounds; T1D-X original satellite q Gram')
print('OPEN: entire remainder after central+satellite; upper spectrum, eta_N; pole directions; RH')