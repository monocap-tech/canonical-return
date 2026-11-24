#!/usr/bin/env python3
"""CR-ALP-038: rigorous fixed-point LDL for all 10 interior + Q_15000.

Source-pinned continuum premises: CR-ALP-027 full Q_15000 lower 1/4,
CR-ALP-034 complete original signed Hmid lower matrix, CR-ALP-036/037
source high-sine formulas. All arithmetic comparisons here are exact integer
with OUTWARD rounding at 10^-25. No Lean, E12 boundary anchors, entire E15000,
fixed-collar full positivity or RH are certified by this script.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json
from cr_alp_038_combined_infinite_tail_rational_audit import I,A,N,S,PI,nu
ROOT=Path(__file__).resolve().parent
BASE=10**25

HFILE=ROOT/'H_mid.tsv'
INFILE=ROOT/'initial.txt'
if not HFILE.exists(): HFILE=ROOT/'CR-ALP-034_midpoint_H_exact_matrix_intervals.tsv'
if not INFILE.exists(): INFILE=ROOT/'CR-ALP-034_initial_rational_intervals_1e25.txt'
assert hashlib.sha256(HFILE.read_bytes()).hexdigest()=='b6a31e623298272fdceac788fc56626ec704e280903ff6baff1dc8bfe0d45f92'
assert hashlib.sha256(INFILE.read_bytes()).hexdigest()=='bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2'


def ceil_div(a,b):
    assert b>0
    return -((-a)//b)

class Q:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        if isinstance(lo,Q):self.lo,self.hi=lo.lo,lo.hi
        elif isinstance(lo,I):
            self.lo=(lo.lo.numerator*BASE)//lo.lo.denominator
            self.hi=ceil_div(lo.hi.numerator*BASE,lo.hi.denominator)
        elif isinstance(lo,F):
            self.lo=(lo.numerator*BASE)//lo.denominator
            self.hi=ceil_div(lo.numerator*BASE,lo.denominator)
        else:
            z=F(lo)
            self.lo=(z.numerator*BASE)//z.denominator
            self.hi=ceil_div(z.numerator*BASE,z.denominator)
        if hi is not None: self.lo,self.hi=int(lo),int(hi)
        assert self.lo<=self.hi
    def __neg__(self):return Q(-self.hi,-self.lo)
    def __add__(self,v):
        z=Q(v);return Q(self.lo+z.lo,self.hi+z.hi)
    __radd__=__add__
    def __sub__(self,v):return self+-Q(v)
    def __mul__(self,v):
        z=Q(v)
        p=[a*b for a in (self.lo,self.hi) for b in (z.lo,z.hi)]
        return Q(min(p)//BASE,ceil_div(max(p),BASE))
    __rmul__=__mul__
    def __truediv__(self,v):
        z=Q(v);assert z.lo>0
        # Upper/lower of (x / y) with strictly positive y, exact outward bounds.
        ratios=[(x*BASE,y) for x in (self.lo,self.hi) for y in (z.lo,z.hi)]
        return Q(min(a//b for a,b in ratios),max(ceil_div(a,b) for a,b in ratios))


def check_sym(d,i,j):
    a,b=d[i,j],d[j,i]
    assert max(a.lo,b.lo)<=min(a.hi,b.hi),(i,j)
    return Q(min(a.lo,b.lo),max(a.hi,b.hi))


def load():
    h={}
    for line in HFILE.read_text().splitlines():
        i,j,l,u=line.split();h[int(i),int(j)]=Q(int(l),int(u))
    assert len(h)==72
    m={}
    for line in INFILE.read_text().splitlines():
        k,l,u=line.split()
        if k.startswith('M_'):
            _,i,j=k.split('_');m[int(i),int(j)]=Q(int(l),int(u))
    return h,m
H,M=load()


def check_block(ids):
    va=nu(1 if ids[0]%2 else 2)
    DD=[];beta=[]
    for j in ids:
        v=nu(j)
        lam=v*(va**2+F(1,4))/(va*(v**2+F(1,4)))
        assert lam.lo>0
        DD.append(Q(v-lam*va))
        c0=1+lam.hi
        c2=v.hi**2+lam.hi*va.hi**2
        c3=v.hi**3+lam.hi*va.hi**3
        B=(S*c0/(2*A)+8*S*c2/A+F(201,200)*c3/A)*(2*A/3)**3
        beta.append(B/isqrt(5*N**5))
    alpha=(2*S/A+F(201,200)*PI.hi/(2*A))*(2*A/3)**2/isqrt(3*N**3)
    # Convert one-sided scalar upper bounds to outward Q upper endpoints.
    alpha_q=Q(Q(alpha).hi,Q(alpha).hi)
    bb=sum((z*z for z in beta),F(0))
    bb_q=Q(Q(bb).hi,Q(bb).hi)
    eps=Q(F(1,10**12))
    n=len(ids)
    mat={}
    for ia,i in enumerate(ids):
        for ib,j in enumerate(ids):
            # 2||QTp||||q|| <= (1/8)||q||² +8||QTp||²
            # ||QTp||²<=2alpha²(Dx)² +2(sum beta²)||x||²
            mat[ia,ib]=check_sym(H,i,j)-16*alpha_q*alpha_q*DD[ia]*DD[ib] \
                       -(16*bb_q if ia==ib else Q(0))-eps*check_sym(M,i,j)
    L={};piv=[]
    for j in range(n):
        p=mat[j,j]
        for k in range(j):p=p-L[j,k]*L[j,k]*piv[k]
        assert p.lo>0,('LDL fail',ids,j,p.lo,p.hi)
        piv.append(p)
        for row in range(j+1,n):
            val=mat[row,j]
            for k in range(j):val=val-L[row,k]*L[j,k]*piv[k]
            L[row,j]=val/p
    lower_check=F(8,10**7) if ids[0]==3 else F(5,10**6)
    assert min(z.lo for z in piv)*F(1,BASE)>lower_check
    print('PASS source-signed mixed interior + Q15000 block:',ids)
    print(' alpha upper',float(F(alpha_q.hi,BASE)),'beta² sum upper',float(F(bb_q.hi,BASE)))
    print(' exact outward LDL lower pivots:',[float(F(p.lo,BASE)) for p in piv])
    print(' PASS: Hmid -16alpha² D D^T -16sum(beta²)I -10^-12 M is strictly positive')
    return {'ids':ids,'positive_pivots':n,
            'pivot_lower_1e25':[str(p.lo)+'/10^25' for p in piv],
            'alpha_rational_upper_1e25':str(alpha_q.hi)+'/10^25',
            'beta_squared_rational_upper_1e25':str(bb_q.hi)+'/10^25',
            'shifted_physical_coercivity':'1/10^12'}

if __name__=='__main__':
    blocks=[check_block([3,5,7,9,11]),check_block([4,6,8,10,12])]
    result={'stage':'CR-ALP-038',
            'status':'EXACT OUTWARD-INTEGER 10 interior mode + entire infinite Q15000 tail mixed Schur PASS',
            'conditional_analytic_theorem':'W(p+q,p+q) > 10^-12||p||² + (1/8)||q||², for nonzero p in span(w3,...,w12), q in Q15000 pole-null log-form domain',
            'premises':'CR-ALP-027 Q15000 coercivity, CR-ALP-034 original signed Hmid lower matrix; original seven negative primes and complete digamma Gamma analytical coefficient identities',
            'Hmid_SHA256':hashlib.sha256(HFILE.read_bytes()).hexdigest(),
            'source_input_SHA256':hashlib.sha256(INFILE.read_bytes()).hexdigest(),
            'arithmetic_scale':'1/10^25, all operations outward rounded',
            'blocks':blocks,
            'not_claimed':'two boundary anchors u1/u2, intermediate 12..15000 physical signed shell, full E15000 Schur, Lean, arbitrary-support, RH'}
    (ROOT/'CR-ALP-038_certified_full_interior_plus_infinite_tail.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CR-ALP-038 EXACT FULL 10-INTERIOR + INFINITE-TAIL PHYSICAL SCHUR AUDIT PASS')
