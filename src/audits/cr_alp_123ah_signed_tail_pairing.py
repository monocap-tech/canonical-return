#!/usr/bin/env python3
"""123AH: exact signed common-tail bilinear diagnostic, no new solves.
Inputs are authenticated original AC, AG and Z workflow ZIP archives.
All acceptance arithmetic uses rational intervals and upward rounding.
"""
from fractions import Fraction as F
from pathlib import Path
from math import factorial,isqrt
import hashlib,json,sys,zipfile
sys.set_int_max_str_digits(0)
G=10**100

def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def box(l,h):
    assert l<=h
    return F(floor(l*G),G),F(ceil(h*G),G)
def iv(v):return box(*map(F,v))
def add(a,b):return box(a[0]+b[0],a[1]+b[1])
def neg(a):return -a[1],-a[0]
def mul(a,b):
    p=[x*y for x in a for y in b];return box(min(p),max(p))
def scale(a,s):return mul(a,(s,s))
def absup(a):return max(map(abs,a))
def abslo(a):return max(F(0),a[0],-a[1])
def root(v):
    assert v>=0
    D=10**110;v=v*D*D
    return F(isqrt(v.numerator//v.denominator)+1,D)
def rounded(v):return F(ceil(F(v)*G),G)
def logseries(x):
    u=(x-1)/(x+1);assert 0<=u<=F(1,3)
    n=120;s=2*sum(u**(2*j+1)/ (2*j+1) for j in range(n))
    rem=2*u**(2*n+1)/((2*n+1)*(1-u*u))
    return box(s,s+rem)
def logn(n):
    k=n.bit_length()-1
    return add(scale(logseries(F(2)),F(k)),logseries(F(n,2**k)))
def atan(x):
    n=100;s=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n));err=x**(2*n+1)/(2*n+1)
    return box(s-err,s+err)
pi=add(scale(atan(F(1,5)),F(16)),scale(atan(F(1,239)),F(-4)))
nu=scale(pi,F(5,11));assert nu[0]>F(15,11)
def sinlo(x):
    m=(x[0]+x[1])/2;r=(x[1]-x[0])/2;assert abs(m)<7
    n=60;s=sum((-1)**j*m**(2*j+1)/factorial(2*j+1) for j in range(n));err=abs(m)**(2*n+1)/factorial(2*n+1)
    # Taylor through degree 120 (even coefficient zero), and sine is 1-Lipschitz.
    v=abs(s)-err-r;assert v>0
    return F(floor(v*10**20),10**20)
def kernel():
    # 123AD entrywise errors, not its Loewner upper matrix.
    J0=2047;th=[mul(nu,logn(n)) for n in (2,3,4,5,7,8,9)]
    C=[[F(0) for _ in range(9)] for _ in range(9)];R=[[F(0) for _ in range(9)] for _ in range(9)]
    C[0][0]=C[0][1]=C[1][0]=F(1,2*J0);C[1][1]=F(1,J0)
    for i in range(2):
        for j in range(2):R[i][j]=F(3,J0**2)
    for i,t in enumerate(th,2):
        C[i][i]=F(1,4*J0);R[i][i]=(3+1/sinlo(scale(t,F(2))))/J0**2
        for j in range(2):R[i][j]=R[j][i]=1/(J0**2*sinlo(t))
        for j in range(2,i):
            u=th[j-2];R[i][j]=R[j][i]=(1/sinlo(add(t,neg(u)))+1/sinlo(add(t,u)))/(2*J0**2)
    q=box(1/nu[1]**2,1/nu[0]**2)
    return [[mul(box(C[i][j]-R[i][j],C[i][j]+R[i][j]),q) for j in range(9)] for i in range(9)]
def read(path,member,expected):
    with zipfile.ZipFile(path) as z:b=z.read(member)
    assert hashlib.sha256(b).hexdigest()==expected
    return json.loads(b)
def coeff(p):
    c=[p['constant'],p['log_coefficient'],*p['cosine_coefficients']];assert len(c)==9
    return [iv(v) for v in c]
def pairgram(a,b,K):
    return sumiv(mul(mul(a[i],K[i][j]),b[j]) for i in range(9) for j in range(9))
def sumiv(seq):
    s=(F(0),F(0))
    for v in seq:s=add(s,v)
    return s
def upperratio(x,E):return str(F(ceil(x/E*10**9),10**9))
def lowerratio(x,E):return str(F(floor(x/E*10**9),10**9))
def main(ac,ag,z):
    A=read(ac,'result.json','a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    Gd=read(ag,'result.json','45719f6cd161b29534983c95f52941d434832bfffd3493d4d421513b436c1abf')
    Z=read(z,'joint-repaired.json','7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897')
    assert A['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER' and Gd['status']=='CONDITIONAL_DIRECTIONAL_INVERSE_BRACKET' and Z['status']=='CONDITIONAL_CONDITIONED_JOINT_CORRECTION'
    assert set(A['sectors'])==set(Gd['sectors'])==set(Z['sectors'])=={'1','2'}
    K=kernel();t0=1/F('7.072938456');dc=F('1.822938456');wc=(t0+1/dc)/2
    result={}
    for k,a in A['sectors'].items():
        g=Gd['sectors'][k];z=Z['sectors'][k][0];n=44-int(k)
        E,U,b2=map(F,(a['E'],a['U'],a['finite_charge_upper']))
        assert E==F(z['joint_energy_lower']) and U==F(z['weighted_inverse_error_upper']) and E>U>0 and b2>=0
        ep=z['signed_pole_tail_receipt']; ec=coeff(ep); el=rounded(F(ep['leading_norm_upper'])); er=rounded(F(ep['pole_remainder_norm_upper'])+F(ep['sine_series_norm_upper']))
        assert len(a['LDL_pivots'])==len(a['tail_receipts'])==len(g['directional_centers'])==n
        sums={key:F(0) for key in ('prefix_lo','prefix_hi','lead_tail_hi','prefix_plus_lead_lo','prefix_plus_lead_hi','remainder_hi','new_center_hi','old_center_hi')}
        for p,tr,c in zip(map(F,a['LDL_pivots']),a['tail_receipts'],g['directional_centers']):
            assert p>0
            pre=iv(c['prefix_and_pole_center_interval']); oldtail=rounded(F(c['tail_pair_abs_upper']))
            zp=tr['pole'];zc=coeff(zp);zl=rounded(F(zp['leading_norm_upper']));zr=rounded(F(zp['pole_remainder_norm_upper'])+F(tr['sine_series_norm_upper']))
            assert min(el,er,zl,zr,oldtail)>=0
            tail=scale(pairgram(zc,ec,K),wc)
            # All actual tail pairings are enclosed by signed leading part + remainder.
            rem=rounded(wc*(zl*er+zr*el+zr*er))
            total=add(pre,tail)
            new=add(total,(-rem,rem));old=add(pre,(-oldtail,oldtail))
            # Intersection of two valid interval enclosures.
            both=(max(new[0],old[0]),min(new[1],old[1]));assert both[0]<=both[1]
            for key,val in [('prefix_lo',abslo(pre)),('prefix_hi',absup(pre)),('lead_tail_hi',absup(tail)),('prefix_plus_lead_lo',abslo(total)),('prefix_plus_lead_hi',absup(total)),('remainder_hi',rem),('new_center_hi',absup(both)),('old_center_hi',absup(old))]:sums[key]+=p*val**2
        assert sums['prefix_lo']>E-U
        # Original AG bracket upper is already outward rounded on 10^-12 normalized grid.
        bracket=F(g['inverse_bracket_error2_over_E_upper'])*E
        center=sums['new_center_hi']
        bound=min(F(a['eta_upper'])*U,(root(center)+root(bracket))**2)
        loss=min(F(a['loss_upper']), U+(root(b2)+root(bound))**2)
        result[k]={key+'_over_E':(lowerratio(v,E) if key.endswith('_lo') else upperratio(v,E)) for key,v in sums.items()}
        result[k].update({'available_budget_over_E':upperratio(E-U,E),'bracket_error2_over_E_upper':upperratio(bracket,E),'directional_charge2_over_E_upper':upperratio(bound,E),'loss_over_E_upper':upperratio(loss,E),'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE','prefix_alone_exceeds_frozen_budget':True})
    print(json.dumps({'status':'CONDITIONAL_SIGNED_COMMON_TAIL_DIAGNOSTIC','sectors':result,'scope':'Same fixed trial, capped comparator and Dbar; no lower bound on full inverse action; no new response solves.'},indent=2))
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
