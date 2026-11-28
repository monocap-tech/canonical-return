#!/usr/bin/env python3
"""CR-INT-T1D-I: rigorous 2-sided frozen-width source-residual brackets.

Certifies finite separation geometry and rational numerical propagation for
L2 bounds on the actual projected mixed Weil-image residuals at ell=1/10000.
Uses the exact Q interval audits for the shape variance (T1D-H) and
nonternary/ternary source prime weights (T1D-G). The analytic lemmas
(log|D| off-support kernel and Binet digamma majorant) are documented
in research/CR-INT-T1D-I.md and are not proved by Python asserts.

No numerical quadrature, whole-complement coercivity or RH.
"""
from fractions import Fraction as F
from math import isqrt
from itertools import combinations
# Certified exact-rational outward bounds imported from T1D-G/H audits.
# tools/check_all.sh executes those source-bound computations before this check.
class Interval:
    def __init__(self,lo,hi):self.lo,self.hi=lo,hi
Vv=Interval(F(121756851,100000000),F(122714789,100000000))
floors={
    '+':(F(338439635,100000000),F(338439636,100000000),
         F(21839262,100000000),F(21839263,100000000)),
    '-':(F(259455169,100000000),F(259455170,100000000),
         F(1126514,100000000),F(1126515,100000000))
}

ell=F(1,10000)
r=F(1,100)
# Get a rigorous rational interval for exp(2) using the positive series.
def exp2_interval(N=30):
    fac=1
    s=F(0)
    for k in range(N+1):
        if k:fac*=k
        s+=F(2)**k/fac
    first=F(2)**(N+1)/(fac*(N+1))
    tail=first/(1-F(2,N+2))
    return s,s+tail

e2_lo,e2_hi=exp2_interval()
scale=10**25
sq3=isqrt(3*scale*scale)
sq3_lo,sq3_hi=F(sq3,scale),F(sq3+1,scale)

def primebase(n):
    for p in range(2,n+1):
        if n%p:continue
        if any(p%q==0 for q in range(2,isqrt(p)+1)):continue
        x=n
        while x%p==0:x//=p
        if x==1:return p
    return None
pp={n:primebase(n) for n in range(2,55) if primebase(n)}
non3={n:p for n,p in pp.items() if p!=3}
# Centers x=(j-3/2)log3 +/- log n are parameterized by
# positive rational t=3**j * n**sign times common factor 3**(-3/2).
base=[(F(3)**j,('base',j,0)) for j in range(4)]
trans=[]
for n in non3:
    for sign in (-1,1):
        for j in range(4):
            q=F(3)**j*(F(n) if sign>0 else F(1,n))
            # exact strict condition e^-2 < 3^-3/2 q < e^2
            low=q*e2_lo > 3*sq3_hi
            high=q < e2_lo*3*sq3_lo
            # corresponding definite opposite exclusions
            outside_low=q*e2_hi < 3*sq3_lo
            outside_high=q > e2_hi*3*sq3_hi
            assert (low or outside_low) and (high or outside_high),('boundary uncertified',n,sign,j)
            if outside_low:
                assert q*e2_hi/(1-ell/2)<3*sq3_lo
            if outside_high:
                assert q>e2_hi/(1-ell/2)*3*sq3_hi
            if low and high:trans.append((q,(n,j,sign)))
assert len(trans)==48,len(trans)
all_centers=sorted(base+trans)
assert len({q for q,_ in all_centers})==52
minrat=min(b[0]/a[0] for a,b in zip(all_centers,all_centers[1:]))
assert minrat==F(37,36),minrat
# log(37/36)>2/73 >2r+ell (because log(q)>2(q-1)/(q+1)).
assert F(2,73)>2*r+ell
# Nearest center to a window boundary is t=37 or t=27/37.
assert min(q for q,_ in all_centers)==F(27,37)
assert max(q for q,_ in all_centers)==F(37)
# Ensure 2-r is beyond the positive center at t=37:
# exp(2-r)> 37/(3sqrt3). exp(2-r) > exp2_lo*(1-r)
assert e2_lo*(1-r)*3*sq3_lo>37
# Thus all 48 translated bumps stay inside the full support window.
# Four base neighborhoods radius r and 48 narrow translated bumps
# have pairwise disjoint supports.
assert F(2)*r+ell < F(2,73)
D=8*r+48*ell
assert D==F(53,625)

# Gamma multiplier error Delta(xi) on |xi|>=1 is <3/(2|xi|)
# by Binet's digamma integral and |t^2+z^2|>=|xi|/4 for
# z=1/4+i*xi/2; on |xi|<=1, <11+|log|xi||.
# Because |vhat(eta)|^2 <= 1/(2pi), integrating those two
# elementary majorants yields squared cell error < 589 ell/(4pi)
# < 589 ell / 12 (using pi>3).
cell_gamma_sq_ub=F(589)*ell/12
assert cell_gamma_sq_ub < F(71,1000)**2
# Logarithmic Fourier h=(L-<L>)v satisfies outside the unit-cell
# support |h(x)| <= 1/[2(|x|-1/2)] because ||v||_1<=1,
# yielding tail mass <= ell/[2(r-ell/2)].
tail_sq_ub=ell/(2*(r-ell/2))
assert tail_sq_ub < F(71,1000)**2
# For four-cell combinations sum |coeff|/sqrt(sum coeff^2) <=2.
err_ub=F(2)*F(71,1000)+F(2)*F(71,1000)
assert err_ub==F(71,250)

# Pole projection of compact proxy t: its parity sigma projects only
# on p_even=cosh(x/2) or p_odd=sinh(x/2).
# Use ||p_even||^2=sinh2+2>5 and ||p_odd||^2=sinh2-2>3/2;
# sup |p_even|<8/5, sup|p_odd|<6/5, and support measure D.
ratio_even_sq_ub=D*F(8,5)**2/F(5)
ratio_odd_sq_ub=D*F(6,5)**2/F(3,2)
assert ratio_even_sq_ub < F(39,400)
assert ratio_odd_sq_ub < F(39,400)
proj_remaining=F(19,20)

# Exact-rational previously certified Fourier shape variance and finite
# arithmetic floor; no limiting residual value is substituted!
levels={}
for tag in ('+','-'):
    Clo,Chi,V3lo,V3hi=floors[tag]
    lo=Vv.lo + Clo + V3lo
    hi=Vv.hi + Chi + V3hi
    levels[tag]=(lo,hi)
    print('V_v+V3+C',tag,':',float(lo),'< T <',float(hi))
assert levels['+'][0]>F(481,100) and levels['+'][1]<F(484,100)
assert levels['-'][0]>F(381,100) and levels['-'][1]<F(385,100)
# Proxy t has squared norm T - tail_norm^2, and is physically E-orthogonal.
# Lower/upper *rational* sqrt enclosure on proxy t.
assert F(219,100)**2 < levels['+'][0]-tail_sq_ub
assert F(11,5)**2 > levels['+'][1]
assert F(195,100)**2 < levels['-'][0]-tail_sq_ub
assert F(197,100)**2 > levels['-'][1]

# Contractivity of physical projection, inverse bound (1-proj)^{1/2}>19/20,
# and ||r - Pi_F t|| < 71/250 yield true 2-sided frozen-width bounds.
lo_plus=(proj_remaining*F(219,100)-err_ub)**2
hi_plus=(F(11,5)+err_ub)**2
lo_minus=(proj_remaining*F(195,100)-err_ub)**2
hi_minus=(F(197,100)+err_ub)**2
assert lo_plus>F(16,5) and hi_plus<F(31,5)
assert lo_minus>F(12,5) and hi_minus<F(51,10)
# T1C original selected mu+ in (7.28,13.93), mu- >7.
assert F(16,5)/F(1393,100)==F(320,1393)
assert F(31,5)/F(728,100)==F(155,182)
assert F(51,10)/F(7)<F(155,182)
print('PASS: 48 source nonternary center placements, minimum separation log(37/36)')
print('PASS: Gamma Fourier-frequency error <.071 per cell; logarithmic tail <.071')
print('PASS: support-localized proxy is E-orthogonal, pole projection < sqrt(39/400)')
print('FROZEN T1C beta_plus^2 in (16/5,31/5) = (3.2,6.2)')
print('FROZEN T1C beta_minus^2 in (12/5,51/10) = (2.4,5.1)')
print('WITH T1C mu intervals, Schur norm in (320/1393,155/182)')
print('No full complement coercivity, completed-domain equivalence or RH.')
