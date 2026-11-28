#!/usr/bin/env python3
"""CR-ALP-028: original signed seven-prime finite-section interval audit.

Exact finite Fraction interval bounds (outward rounded to 1e-38), Machin pi,
atanh log, Taylor trigonometry. Verifies literal original paired translations
on several pole-null physical interior columns and rank-two Gram identities.
NOT a Gamma evaluation, a 15k Schur certificate, Lean, or RH.
"""
from fractions import Fraction as F
from math import factorial, isqrt

SCALE = 10**38
A = F(11, 10)
POWERS = {2:2,3:3,4:2,5:5,7:7,8:2,9:3}


def out_lo(v):
    x = F(v)
    return F((x.numerator * SCALE)//x.denominator, SCALE)


def out_hi(v):
    x = F(v)
    return F(-((-x.numerator * SCALE)//x.denominator), SCALE)

class I:
    __slots__=('lo','hi')
    def __init__(self, lo, hi=None):
        if isinstance(lo,I) and hi is None:
            self.lo,self.hi=lo.lo,lo.hi
            return
        if hi is None: hi=lo
        self.lo=out_lo(lo);self.hi=out_hi(hi)
        assert self.lo<=self.hi
    def __add__(self,other):
        b=I(other)
        return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,other):return self+-I(other)
    def __rsub__(self,other):return I(other)+-self
    def __mul__(self,other):
        b=I(other)
        vals=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return I(min(vals),max(vals))
    __rmul__=__mul__
    def __truediv__(self,other):
        b=I(other)
        assert b.lo>0 or b.hi<0, 'interval denominator contains zero'
        return self*I(1/b.hi,1/b.lo)
    def __rtruediv__(self,other):return I(other)/self
    def __repr__(self):return f'[{float(self.lo):.12g}, {float(self.hi):.12g}]'
    def strict_gt(self,x):return self.lo>F(x)
    def strict_lt(self,x):return self.hi<F(x)


def atanh_log_1_2(z, count=65):
    assert F(1)<=z<=F(2)
    u=(z-1)/(z+1)
    term=u; lo=F(0)
    for k in range(count):
        lo+=2*term/(2*k+1)
        term*=u*u
    hi=lo+2*term/((2*count+1)*(1-u*u))
    return (lo,hi)


def log_int(n):
    k=n.bit_length()-1
    al,ah=atanh_log_1_2(F(n,2**k)); bl,bh=atanh_log_1_2(F(2))
    return I(al+k*bl,ah+k*bh)


def atan_interval(z, terms=55):
    # Alternating arctan with exact remainder between adjacent partial sums.
    a=F(0)
    for k in range(terms):
        a+=(-1)**k*z**(2*k+1)/(2*k+1)
    nxt=(-1)**terms*z**(2*terms+1)/(2*terms+1)
    return I(min(a,a+nxt),max(a,a+nxt))

pi=16*atan_interval(F(1,5))-4*atan_interval(F(1,239))
assert pi.lo>F(3) and pi.hi<F(22,7)
assert log_int(9).hi<2*A<log_int(10).lo
assert set(POWERS)=={2,3,4,5,7,8,9}


def sqrt_int(n):
    k=isqrt(n*SCALE*SCALE)
    return I(F(k,SCALE),F(k+1,SCALE))


def trig(z,kind):
    z=I(z)
    # q is only a convenient integer for periodic reduction; correctness
    # does not rely on the accuracy of its floating-point selection.
    q=round(float((z.lo+z.hi)/2)/float(pi.lo+pi.hi))
    reduced=z-2*q*pi
    c=(reduced.lo+reduced.hi)/2
    assert abs(c)<F(4), 'unexpected phase reduction'
    width=max(abs(reduced.lo-c),abs(reduced.hi-c))
    n=25
    if kind=='sin':
        val=sum(((-1)**k*c**(2*k+1)/factorial(2*k+1) for k in range(n)),F(0))
        rem=abs(c)**(2*n+1)/factorial(2*n+1)
    elif kind=='cos':
        val=sum(((-1)**k*c**(2*k)/factorial(2*k) for k in range(n)),F(0))
        rem=abs(c)**(2*n)/factorial(2*n)
    else: raise ValueError(kind)
    return I(val-width-rem,val+width+rem)


nu={j:I(j)*pi/(2*I(A)) for j in range(1,11)}
t={n:log_int(n) for n in POWERS}
weights={n:log_int(p)/sqrt_int(n) for n,p in POWERS.items()}
S=sum(weights.values(),I(0))
assert S.strict_lt(F(177,50))
assert S.strict_gt(F(3))

# Pair of original shifts, U_ij(t)=<e_i,S_t e_j>+<e_i,S_-t e_j>.
# Source signed prime matrix P_ij = - sum Lambda(n)/sqrt(n) U_ij(log n).

def U(i,j,n):
    if (i+j)%2: return I(0)  # exact parity zero
    v=nu[i]; w=nu[j]; z=t[n]
    if i==j:
        return (2*I(A)-z)/I(A)*trig(v*z,'cos')+trig(v*z,'sin')/(I(A)*v)
    return 2*(v*trig(w*z,'sin')-w*trig(v*z,'sin'))/(I(A)*(v*v-w*w))

prime_cache={}
def P(i,j):
    if (i+j)%2: return I(0)
    key=tuple(sorted((i,j)))
    if key not in prime_cache:
        prime_cache[key]=-sum((weights[n]*U(i,j,n) for n in POWERS),I(0))
    return prime_cache[key]


def anchor(j):return 1 if j%2 else 2

def lam(j):
    k=anchor(j)
    return nu[j]*(nu[k]*nu[k]+I(F(1,4)))/(nu[k]*(nu[j]*nu[j]+I(F(1,4))))

def Pw(i,j):
    ki=anchor(i);kj=anchor(j)
    li=lam(i);lj=lam(j)
    return P(i,j)-li*P(ki,j)-lj*P(i,kj)+li*lj*P(ki,kj)

# The true original signed PRIME-ONLY low block is not PSD; Gamma is absent
# from these individual comparisons, so none is an RH witness.
for j,strict_lower,strict_upper in [
    (3,F(3,5),F(7,10)),(4,F(3,10),F(2,5)),
    (5,F(1,100),F(1,40)),(6,F(-1,4),F(-1,10)),
    (7,F(-2,5),F(-1,4)),(8,F(-3,5),F(-2,5))]:
    v=Pw(j,j)
    assert v.strict_gt(strict_lower) and v.strict_lt(strict_upper),(j,v)
    print(f'PASS original signed prime-only W(w_{j},w_{j}) in {v}')
for j,k in [(3,4),(4,5),(5,6),(3,6),(7,8)]:
    assert Pw(j,k).lo==Pw(j,k).hi==0
print('PASS: exact even/odd original signed prime cross pairs vanish')

# Algebraic Grammian regression: anchor e_k and paired pole functional
# h_j = nu_j/(nu_j^2+1/4), so w_j=e_j-(h_j/h_k)e_k has zero moment.
for k in (1,2):
    h_k=nu[k]/(nu[k]*nu[k]+I(F(1,4)))
    for j in range(3,11):
        if j%2!=k%2: continue
        h_j=nu[j]/(nu[j]*nu[j]+I(F(1,4)))
        # direct interval zero expression could expand and incur enclosure
        # error; the identity lambda=h_j/h_k is algebraic by definition.
        ratio=lam(j)
        test=ratio*h_k-h_j
        assert test.lo<=0<=test.hi
        assert abs(float(test.lo))<1e-28 and abs(float(test.hi))<1e-28
print('PASS rank-two anchor cancellation algebra in both physical parities')
print('CR-ALP-028 EXACT OUTWARD-FRACTION PRIME INTERVAL AUDIT PASS')
print('NOTE: does not audit Gamma, full H, cross-tail norm or Schur positivity')
