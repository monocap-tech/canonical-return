#!/usr/bin/env python3
"""CR-ALP-034 independent exact-Fraction verification of the midpoint Gamma lower matrix.

The source-analytic input is the genuine original signed Weil form on a fixed
pole-null collar. This script is a finite arithmetic certificate only. It
checks physical mass and the two 6x6 shifted original-form lower blocks,
and a separate *upper* error envelope from the convexity Taylor remainder.
No Lean compilation, full 15000 section, or RH claim.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

ROOT = Path(__file__).resolve().parent
S = 10**25
INPUT = ROOT/'CR-ALP-034_initial_rational_intervals_1e25.txt'
MATRIX = ROOT/'CR-ALP-034_midpoint_H_exact_matrix_intervals.tsv'
assert hashlib.sha256(INPUT.read_bytes()).hexdigest() == 'bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2'

class Intv:
    def __init__(self, lo, hi=None):
        if isinstance(lo, Intv): self.lo, self.hi = lo.lo, lo.hi
        else: self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi
    def __neg__(self):return Intv(-self.hi,-self.lo)
    def __add__(self, val):
        v=Intv(val);return Intv(self.lo+v.lo,self.hi+v.hi)
    __radd__=__add__
    def __sub__(self,val):return self+(-Intv(val))
    def __mul__(self,val):
        v=Intv(val);p=[x*y for x in (self.lo,self.hi) for y in (v.lo,v.hi)]
        return Intv(min(p),max(p))
    __rmul__=__mul__
    def __truediv__(self,val):
        v=Intv(val);assert v.lo>0
        return self*Intv(1/v.hi,1/v.lo)


def fixed(x,y):return Intv(F(int(x),S),F(int(y),S))

input_dict={}
for ln in INPUT.read_text().splitlines():
    k,l,u=ln.split();input_dict[k]=fixed(l,u)
assert len(input_dict)==177
H_dict={}
for ln in MATRIX.read_text().splitlines():
    i,j,l,u=ln.split();H_dict[int(i),int(j)]=fixed(l,u)
assert len(H_dict)==72

EVEN=(3,5,7,9,11,1)
ODD=(4,6,8,10,12,2)


def interval_matrix(dct,i,j):
    a,b=dct[i,j],dct[j,i]
    assert max(a.lo,b.lo)<=min(a.hi,b.hi),(i,j)
    return Intv(min(a.lo,b.lo),max(a.hi,b.hi))

def H(i,j): return interval_matrix(H_dict,i,j)
def M(i,j):
    a=input_dict[f'M_{i}_{j}'];b=input_dict[f'M_{j}_{i}']
    assert max(a.lo,b.lo)<=min(a.hi,b.hi),(i,j)
    return Intv(min(a.lo,b.lo),max(a.hi,b.hi))


def pivots(order,shift,typ):
    L={};ps=[]
    A=M if typ=='mass' else H
    for t,i in enumerate(order):
        p=A(i,i)-(shift*M(i,i) if typ!='mass' else Intv(0))
        for k in range(t):p=p-L[t,k]*L[t,k]*ps[k]
        assert p.lo>0,(typ,order,t,p.lo,p.hi)
        ps.append(p)
        for row in range(t+1,len(order)):
            j=order[row]
            v=A(j,i)-(shift*M(j,i) if typ!='mass' else Intv(0))
            for k in range(t):v=v-L[row,k]*L[t,k]*ps[k]
            L[row,t]=v/p
    return ps

shift=F(3,10**12)
report={}
for label,order in [('even',EVEN),('odd',ODD)]:
    mass=pivots(order,F(0),'mass')
    form=pivots(order,shift,'form')
    lower=[v.lo for v in form]
    report[label]={'physical_mass_pivots_positive':len(mass),
                   'shifted_form_pivots_positive':len(form),
                   'shifted_fifth_pivot_lower_rounded_down_1e15':str(lower[4]*10**15//1)+'/10^15',
                   'min_shifted_pivot_lower_rounded_down_1e15':str(min(lower)*10**15//1)+'/10^15'}

# The Loewner-convex Jensen property reduces every group to its mean beta.
# For f in H^1(I), J_beta''(f,f) <= 48*D(f)/beta^5+24*T(f)/beta^4.
# Taylor around the arithmetic mean beta and sum the variance of consecutive
# values beta_m=2m+1/2: sum (beta_m-beta_mean)^2=ell(ell^2-1)/3.
# Hence the group remainder upper is
#   8 ell(ell^2-1)/beta_L^5 * D + 4 ell(ell^2-1)/beta_L^4 * T.
# Round EACH positive exact rational term UP at scale 10^30.
Q=10**30

def ceil_div(a,b):
    assert a>=0 and b>0
    return (a+b-1)//b

L=30000;stop=2000000
cD=cT=0;number=0
while L<stop:
    ell=min(max(1,L//1000),stop-L)
    t=4*L+1
    cD+=ceil_div(256*ell*(ell*ell-1)*Q,t**5)
    cT+=ceil_div(64*ell*(ell*ell-1)*Q,t**4)
    L+=ell;number+=1
assert number==4219
beta_num=4*stop+1
cD+=ceil_div(32*Q,beta_num**3)+ceil_div(4*Q,beta_num**2)
cT+=ceil_div(16*Q,beta_num**2)+ceil_div(4*Q,beta_num)
assert cD*15*10**12<Q
assert cT*1999000<Q
report['remainder']={'midpoint_groups':number,'index_coverage':'[30000,2000000)',
  'cD_rational_upper':f'{cD}/{Q}', 'cT_rational_upper':f'{cT}/{Q}',
  'two_sided_form_error':'0 <= W-H_mid <= D/(15e12) + T/1999000 on finite H1 sources'}
report['stage']='CR-ALP-034'
report['status']='EXACT 12-MODE POSITIVE, FULL 15000-MODE SCHUR OPEN'
report['margin']='W(f,f) > (3/10^12)*||f||^2 for nonzero f in original physical E12'
report['source_input_sha256']=hashlib.sha256(INPUT.read_bytes()).hexdigest()
report['midpoint_matrix_sha256']=hashlib.sha256(MATRIX.read_bytes()).hexdigest()
report['not_claimed']='Lean formalization, entire fixed collar, cross-tail Schur, F41 or RH'
(ROOT/'CR-ALP-034_certified_12_mode_result.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: original physical mass positive in both parity blocks')
print('PASS: exact Fraction LDL of H_mid - (3/10^12)M gives 6/6 + 6/6 positive pivots')
print('PASS: Loewner convex midpoint grouped upper-error bound, cD<1/(15e12), cT<1/1999000')
print('SHA256 matrix:',report['midpoint_matrix_sha256'])
print('CR-ALP-034 exact arithmetic PASS; original signed fixed-collar global positivity OPEN')
