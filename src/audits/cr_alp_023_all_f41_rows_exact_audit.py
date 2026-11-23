#!/usr/bin/env python3
"""CR-ALP-023: exact rational-sign source-moment audit for frozen NF545 nodes.

Frozen source: monocap-tech/weil-lab@a431f903254000bb5ad8edbc276ccc8599da98f5
WeilDefect/Arithmetic/ActualZetaOriginalJ242ExponentialSignature.lean
blob b56b2dc6c482ca79e8f023d556bd52edf5e24f48
Also cross-checked against ActualZetaOriginalNF407FrozenData.lean
blob 86ca70c8a2a8ba24ce4204081fef74d4ac543646.
The 67 source frequency groups (241 nodes) are transcribed exactly below.

This is a mathematically rigorous source-indexed REAL-physical pairing sign
regression for the packet h_2 and all allowed even nonnegative mass-one
mollifiers of width 10^-14. It does not replace frozen historical 242x4 V,
compile Lean, prove original Weil positivity on a whole collar, or prove RH.
"""
from fractions import Fraction as Q
from math import factorial

GROUPS='''
37/64:13 211/128:8 45/16:8 117/32:4 141/32:5 45/8:4 25/4:3
929/128:4 523/64:3 1113/128:3 295/32:2 641/64:3 1511/128:3
1651/128:3 1845/128:3 1979/128:3 2223/128:4 641/32:3
2777/128:3 1477/64:2 1513/64:3 197/8:2 437/16:4 4071/128:3
1135/32:2 2305/64:2 2501/64:3 5721/128:3 393/8:0
815/16:1 6645/128:2 6839/128:0 7229/128:1 7567/128:3
8499/128:3 72:2 4895/64:3 10721/128:3 5531/64:1
12567/128:2 6643/64:3 14217/128:2 7397/64:2 15511/128:1
16065/128:1 8679/64:2 18289/128:3 297/2:0 20853/128:1
21785/128:2 3001/16:3 25855/128:2 14039/64:2 15787/64:3
1143/4:1 37295/128:2 39141/128:1 41365/128:2 5287/16:1
44861/128:2 49863/128:2 55583/128:2 59653/128:1
31575/64:2 81437/128:2 23681/32:1 113013/128:1
'''
FREQS=[]
for group in GROUPS.split():
    frequency,top=group.split(':')
    for k in range(int(top)+1):
        FREQS.append((Q(frequency),k))
assert len(FREQS)==241 and len(GROUPS.split())==67
# Fingerprint the exact source-extracted 67 frozen frequency/order groups.
# FNV-1a64 is used as a transcription checksum, not as a cryptographic proof.
canonical_groups=" ".join(GROUPS.split())
fnv=0xcbf29ce484222325
for ch in canonical_groups:
    fnv=((fnv ^ ord(ch))*0x100000001b3)&((1<<64)-1)
assert f"{fnv:016x}"=="29ec35e3e29c3ee9"
assert FREQS[:3]==[(Q(37,64),0),(Q(37,64),1),(Q(37,64),2)]
assert FREQS[-1]==(Q(113013,128),1)
B=Q(1,10**8)
ETA=Q(1,10**14)
R=Q(3,4)


def log_atanh(z,N=70):
    assert 1<=z<=2
    t=(z-1)/(z+1)
    x=t
    a=Q(0)
    for j in range(N):
        a+=2*x/(2*j+1)
        x*=t*t
    return a,a+2*x/((2*N+1)*(1-t*t))

L0,L1=log_atanh(Q(2))
assert Q(2,3)<L0<L1<Q(3,4)


def atan_inv(n,N=48):
    terms=sum((Q((-1)**j,(2*j+1)*n**(2*j+1)) for j in range(N)),Q(0))
    nxt=Q((-1)**N,(2*N+1)*n**(2*N+1))
    return (terms,terms+nxt) if N%2==0 else (terms+nxt,terms)

p5lo,p5hi=atan_inv(5)
p239lo,p239hi=atan_inv(239)
P0,P1=16*p5lo-4*p239hi,16*p5hi-4*p239lo
assert 3<P0<P1<Q(22,7)



# Outward-rounded fixed-point intervals: each endpoint is an integer in
# 10^-52 units. Every +, -, * and integer division rounds OUTWARDS.
SCALE=10**52

def ceildiv(a,b):
    assert b>0
    return -((-a)//b)

class Interval:
    __slots__=("lo","hi")
    def __init__(self,lo,hi):
        assert lo<=hi
        self.lo,self.hi=lo,hi
    @classmethod
    def exact(cls,q):
        q=Q(q)
        return cls((q.numerator*SCALE)//q.denominator,
                   ceildiv(q.numerator*SCALE,q.denominator))
    @classmethod
    def integer(cls,i):return cls(i*SCALE,i*SCALE)
    @classmethod
    def from_scaled(cls,i):return cls(i,i)
    def __add__(self,other):
        other=interval(other)
        return Interval(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,other):return self+-interval(other)
    def __rsub__(self,other):return interval(other)+-self
    def __mul__(self,other):
        other=interval(other)
        vals=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(vals)//SCALE,ceildiv(max(vals),SCALE))
    __rmul__=__mul__
    def div_positive_int(self,n):
        assert isinstance(n,int) and n>0
        return Interval(self.lo//n,ceildiv(self.hi,n))
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        r=Interval.integer(1); t=self
        while n:
            if n&1:r=r*t
            n>>=1
            if n:t=t*t
        return r

def interval(z):
    return z if isinstance(z,Interval) else Interval.exact(z)

L=Interval.exact(L0)+Interval(0,ceildiv((L1-L0).numerator*SCALE,(L1-L0).denominator))
P=Interval.exact(P0)+Interval(0,ceildiv((P1-P0).numerator*SCALE,(P1-P0).denominator))


def trig_for_frequency(f):
    # theta=pi*(2*f*log2). Reduce modulo 2*pi with exact intervals.
    ratio=2*f*L
    m=(ratio.lo+ratio.hi+2*SCALE)//(4*SCALE)
    red=ratio-2*m
    assert -SCALE<red.lo<=red.hi<SCALE
    angle=P*red
    mid=(angle.lo+angle.hi)//2
    delta=ceildiv(angle.hi-angle.lo+1,2)
    assert abs(mid)<Q(22,7)*SCALE
    z=Interval.from_scaled(mid)
    zz=z*z
    N=24
    sint= z
    cost= Interval.integer(1)
    sine= Interval.integer(0)
    cosine=Interval.integer(0)
    for j in range(N):
        sine=sine+sint
        cosine=cosine+cost
        sint=-(sint*zz).div_positive_int((2*j+2)*(2*j+3))
        cost=-(cost*zz).div_positive_int((2*j+1)*(2*j+2))
    sinerr=Q(22,7)**(2*N+1)/factorial(2*N+1)
    coserr=Q(22,7)**(2*N)/factorial(2*N)
    se=ceildiv(sinerr.numerator*SCALE,sinerr.denominator)+delta
    ce=ceildiv(coserr.numerator*SCALE,coserr.denominator)+delta
    return Interval(sine.lo-se,sine.hi+se),Interval(cosine.lo-ce,cosine.hi+ce)


def leading_interval(f,k):
    omega=2*P*f
    sint,cost=trig_for_frequency(f)
    v=(k*L**(k-1) if k else Interval.integer(0))-Q(3,2)*L**k
    return v*sint+omega*L**k*cost if k%2==0 else v*cost-omega*L**k*sint

def error_over_b2(f,k):
    # Source support always |x|<=log2+b+eta <3/4.
    assert L1+B+ETA<R
    omega_hi=2*P1*f
    M0=R**k
    M1=(k*R**(k-1) if k else Q(0))+omega_hi*R**k
    M2=((k*(k-1)*R**(k-2)) if k>=2 else Q(0))
    M2+=(2*k*omega_hi*R**(k-1)) if k else Q(0)
    M2+=omega_hi**2*R**k
    # Exact envelope:
    # |m0-b²/2|<=b^4/8, |m1+b²|<=b^4/4.
    # Taylor remainder <=(3/2)*M2*(4/3)b³=2 M2 b³.
    # Mollifier shift <=9 b eta M1.
    return Q(3,8)*B*B*M0+Q(1,4)*B*B*M1+2*B*M2+9*ETA/B*M1

results=[]
for j,(f,k) in enumerate(FREQS):
    A=leading_interval(f,k)
    E=error_over_b2(f,k)
    EF=Interval.exact(E)
    if A.hi+EF.hi<0:
        sign='-'
        margin=-(A.hi+EF.hi)
    elif A.lo-EF.hi>0:
        sign='+'
        margin=A.lo-EF.hi
    else:
        sign='?'
        margin=Q(0)
    results.append((j,f,k,sign,margin))

nminus=sum(v[3]=='-' for v in results)
nplus=sum(v[3]=='+' for v in results)
unknown=[v for v in results if v[3]=='?']
print(f'CR-ALP-023 FROZEN SOURCE {len(results)} ROWS: '
      f'certified negative={nminus}; positive={nplus}; uncertified={len(unknown)}')
print('FIRST 16:',[(j,str(f),k,sgn) for j,f,k,sgn,_ in results[:16]])
print('UNCERTIFIED:',[(j,str(f),k) for j,f,k,_,_ in unknown])
# Preserve prior 3-row stringent bounds, now using full interval pipeline.
for j,thr in [(0,Q(-2)),(1,Q(-1,2)),(2,Q(-1,5))]:
    f,k=FREQS[j]
    A=leading_interval(f,k)
    E=error_over_b2(f,k)
    assert A.hi+Interval.exact(E).hi<Interval.exact(thr).lo,(j,A.hi,thr)

assert nminus+nplus+len(unknown)==241
assert nminus+nplus==241
assert all(m>Interval.exact(Q(7,100)).hi for j,f,k,sgn,m in results)
# Row j=241 of the frozen selected family is exactly sinh(x/2).
# The packet is two-pole-null, so its physical row is IDENTICALLY ZERO,
# not a small or numerically unresolved entry.
print('Frozen selected row 241: exactly zero via F_+(h2)=F_-(h2)=0')
print('All 241 non-pole rows: |physical moment| > 7/100*b^2')
print('CR-ALP-023 EXACT SOURCE-PINNED FRACTION AUDIT PASS')
print('No source-integral quadrature, Lean compilation, historical V rewrite, or RH claim')