#!/usr/bin/env python3
"""CR-INT-T1D-G: exact rational and SymPy audit of the sharp-rate residual floor.

Checks a=2 source prime-power census; four-bump two-Mellin-pole annihilation;
the exact in-cell 3-power residual variances; all non-ternary translated-center
multiplicities; and outward Q-intervals for the limiting arithmetic terms.

The Fourier/digamma scaling limit and pole-projection weak convergence are
ANALYTIC arguments in research/CR-INT-T1D-G.md, not verified by this script.
No certified finite-width residual-norm integration, complement positivity,
Lean certification, or RH claim.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as sp

def primepower(n):
    if n < 2: return None
    for p in range(2,n+1):
        if n % p: continue
        if any(p % q == 0 for q in range(2,isqrt(p)+1)): continue
        m=n
        while m % p == 0: m//=p
        if m==1: return p
    return None

pp={n:primepower(n) for n in range(2,55) if primepower(n)}
assert len(pp)==24 and [n for n,p in pp.items() if p==3]==[3,9,27]

def exp2_bounds(N=32):
    def fact(k):
        z=1
        for j in range(1,k+1): z*=j
        return z
    p=sum((F(2)**k/F(fact(k)) for k in range(N+1)),F(0))
    rem=(F(2)**(N+1)/F(fact(N+1)))/(1-F(2,N+2))
    return p,p+rem

e2lo,e2hi=exp2_bounds()
scale=10**24; root=isqrt(3*scale*scale)
rlo,rhi=F(root,scale),F(root+1,scale)
assert rlo*rlo<3<rhi*rhi

# For + log n from cell j, center x_j+log n is inside (-2,2)
# iff n < exp(2)*3^(3/2-j). The negative shifts are reflected.
for n in pp:
    hits=0
    for j in range(4):
        power=3-2*j
        if power>=0:
            lo,hi=e2lo*rlo**power,e2hi*rhi**power
        else:
            lo,hi=e2lo/rhi**(-power),e2hi/rlo**(-power)
        assert n<lo or n>hi,(n,j)
        hits+=int(n<lo)
    hits*=2
    assert hits==(6 if n<=4 else 4 if n<=11 else 2 if n<=37 else 0),(n,hits)

A=[2,4]
B=[5,7,8,11]
C=[13,16,17,19,23,25,29,31,32,37]
non3=[n for n,p in pp.items() if p!=3]
assert sorted(non3)==sorted(A+B+C+[41,43,47,49,53])

sqrt3=sp.sqrt(3); s=4/sqrt3
X=sp.Matrix([[1,0],[-s,1],[1,-s],[0,1]])
L=sp.Matrix([[sqrt3**j for j in range(4)],
             [sqrt3**(-j) for j in range(4)]])
assert sp.simplify(L*X)==sp.zeros(2,2)
z=sp.symbols('z',positive=True)
kappa={1:z/sqrt3,2:z/3,3:z/(3*sqrt3)}
P3=sp.Matrix([[0 if i==j else -kappa[abs(i-j)]
               for j in range(4)] for i in range(4)])
for sign,coeff,eps in [('+',sp.Matrix([1,1]),1),
                       ('-',sp.Matrix([1,-1]),-1)]:
    c=X*coeff; norm=sp.simplify(c.dot(c))
    lam=sp.simplify((c.T*P3*c)[0]/norm)
    variance=sp.simplify((((P3-lam*sp.eye(4))*c).T*
                         ((P3-lam*sp.eye(4))*c))[0]/norm)
    assert sp.simplify(lam-(41-eps*5*sqrt3)*z/73)==0
    assert sp.simplify(variance-(507+eps*264*sqrt3)*z*z/5329)==0

def log_q(q,N=75):
    q=F(q); assert 1<=q<=2
    t=(q-1)/(q+1); term=t; total=F(0)
    for j in range(N):
        total+=2*term/F(2*j+1)
        term*=t*t
    return total,total+2*term/(F(2*N+1)*(1-t*t))

def log_n(n):
    k=n.bit_length()-1
    a,b=log_q(F(n,2**k)); c,d=log_q(2)
    return a+k*c,b+k*d

def sq_kappa(n):
    a,b=log_n(pp[n])
    return a*a/n,b*b/n

def sigma(group):
    lo=hi=F(0)
    for n in group:
        a,b=sq_kappa(n);lo+=a;hi+=b
    return lo,hi

ag,bg,cg=sigma(A),sigma(B),sigma(C)
assert cg[0]>ag[1]  # C>A (strict, exact fractions)

def format_out(lo,hi,d=8):
    scale=10**d
    a=lo.numerator*scale//lo.denominator
    b=-((-hi.numerator*scale)//hi.denominator)
    def dec(v):
        sgn='-' if v<0 else ''
        v=abs(v)
        return sgn+str(v//scale)+'.'+str(v%scale).zfill(d)
    return '['+dec(a)+', '+dec(b)+']'

floors={}
for sign,eps in [('+',-1),('-',1)]:
    ts=[(F(4)/r+eps)**2 for r in [rlo,rhi]]
    tlo,thi=min(ts),max(ts)
    def ca(t): return (1+2*t)/(1+t)
    def cc(t): return 1/(1+t)
    # Both coefficients positive; ca increases, cc decreases with t.
    lo=ca(tlo)*ag[0]+bg[0]+cc(thi)*cg[0]
    hi=ca(thi)*ag[1]+bg[1]+cc(tlo)*cg[1]
    # V3 plus => +264sqrt3, V3 minus => -264sqrt3.
    sign3=1 if sign=='+' else -1
    l3,h3=log_n(3)
    vl=(F(507)+F(264*sign3)*(rlo if sign3>0 else rhi))*l3*l3/F(5329)
    vh=(F(507)+F(264*sign3)*(rhi if sign3>0 else rlo))*h3*h3/F(5329)
    assert vl>0
    floors[sign]=(lo,hi,vl,vh)
    print(sign,'C_non3',format_out(lo,hi),'V3',format_out(vl,vh),
          'arithmetic floor',format_out(lo+vl,hi+vh))

assert floors['+'][0]>F(338,100)
assert floors['-'][0]>F(259,100)
assert floors['+'][0]+floors['+'][2]>F(18,5)
assert floors['+'][0]>floors['-'][1]
assert floors['+'][2]>floors['-'][3]
print('PASS T1D-G: 24 source prime powers; exact inside-center census;')
print('two-pole annihilator; ternary variances; outward rational floors;')
print('limiting plus arithmetic residual-square coefficient >18/5.')
print('ANALYTIC ONLY: beta^2 -> V_profile + V3 + C_non3 and')
print('log(1/ell)*||R_ell|| -> V_profile + V3_plus + C_plus.')
print('OPEN: finite-width interval, complement coercivity, RH.')
