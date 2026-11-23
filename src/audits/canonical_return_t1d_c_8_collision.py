#!/usr/bin/env python3
"""CR-INT-T1D-C: exact prime-8 source-geometry and rational interval audit.

Verifies the pole-null 4-to-2 compression, all prime-power windows below
log(31/27), and conservative rational bounds for a shape-uniform positive
2D Weil-form margin. The fiber-path autocorrelation and Gamma integral
inequalities are proved in research/CR-INT-T1D-C.md, NOT checked by Lean.
"""
from fractions import Fraction as F
from math import isqrt
import sympy as sp

def log_interval(x, terms=240):
    x=F(x); assert x>=1
    z=(x-1)/(x+1)
    if not z: return F(0),F(0)
    term=z; partial=F(0)
    for j in range(terms):
        partial += 2*term/F(2*j+1)
        term *= z*z
    return partial, partial+2*term/(F(2*terms+1)*(1-z*z))

def sqrt_interval(n,digits=25):
    scale=10**digits
    p=isqrt(n*scale*scale)
    return F(p,scale),F(p+1,scale)

L=F(7,50)
assert 1+L+L*L/2>F(31,27)
assert F(29,27)<F(27,25)<F(9,8)<F(31,27)
assert F(31,27)<F(29,27)**2
assert F(31,27)<F(9,8)**2

def prime_power(n):
    if n<2:return False
    for p in range(2,n+1):
        if n%p:continue
        if any(p%q==0 for q in range(2,isqrt(p)+1)):continue
        y=n
        while y%p==0:y//=p
        if y==1:return True
    return False

windows=[[n for n in range(2,40)
          if prime_power(n) and F(27*m,31)<n<F(31*m,27)]
         for m in (3,9,27)]
assert windows==[[3],[8,9],[25,27,29]]

r=sp.sqrt(3);s=4/r
X=sp.Matrix([[1,0],[-s,1],[1,-s],[0,1]])
N=sp.Matrix([[r**j for j in range(4)], [r**(-j) for j in range(4)]])
assert sp.simplify(N*X)==sp.zeros(2,2)
D2=sp.zeros(4)
for j in (0,1):
    D2[j,j+2]=1
    D2[j+2,j]=1
B8=sp.simplify(X.T*D2*X)
M=sp.simplify(X.T*X)
assert B8==sp.Matrix([[2,-8*r/3],[-8*r/3,2]])
assert M==sp.Matrix([[sp.Rational(22,3),-8*r/3],
                    [-8*r/3,sp.Rational(22,3)]])
plus=sp.simplify(-(B8[0,0]+B8[0,1])/(M[0,0]+M[0,1]))
minus=sp.simplify(-(B8[0,0]-B8[0,1])/(M[0,0]-M[0,1]))
assert sp.simplify(plus-(15+32*r)/73)==0
assert sp.simplify(minus-(15-32*r)/73)==0
assert sp.simplify(minus-plus+64*r/73)==0

root2=sqrt_interval(2); root3=sqrt_interval(3); root5=sqrt_interval(5)
ln2=log_interval(2); ln3=log_interval(3); ln5=log_interval(5)
ln29=log_interval(29)
# Fiber-path layers: |R(t)| <= cos(pi/(j+2)) on ell/(j+1)<t<ell/j.
# Exact gains j=1,2,3; for j>=4, pi>3 and
# 1-cos(x)>=x*x/2-x**4/24 at x=3/(j+2).
gain=ln2[0]/2
gain+=(1-root2[1]/2)*log_interval(F(3,2))[0]
gain+=(1-(1+root5[1])/4)*log_interval(F(4,3))[0]
for j in range(4,9):
    x=F(3,j+2)
    gain+=(x*x/2-x**4/24)*log_interval(F(j+1,j))[0]
# The canonical T1D-B proof establishes C_gamma<25/8.
g0_lo=log_interval(F(100,7),300)[0]-F(25,8)-L+(1-L/2)*gain

# Gamma cross bounds: |A_k|/J <= ell*w(k log3-ell), k=1,2,3.
# exp(ell/2)<=1/(1-L/2), exp(2ell)<=1/(1-2L).
rlo,rhi=root3
w=[]
for k in (1,2,3):
    invroot=(F(1,3**(k//2)) if k%2==0
             else F(1,3**(k//2))/rlo)
    w.append(invroot/(1-L/2)/(1-F(1,3**(2*k))/(1-2*L)))
pcoef=[(-83+76*rhi)/146,(15+32*rhi)/73,(33+12*rhi)/146]
mcoef=[(83+76*rhi)/146,(32*rhi-15)/73,(33-12*rlo)/146]
assert all(v>0 for v in pcoef+mcoef)
errp=L*sum((pcoef[j]*w[j] for j in range(3)),F(0))
errm=L*sum((mcoef[j]*w[j] for j in range(3)),F(0))
k29=ln29[1]/sqrt_interval(29)[0]
k25=ln5[1]/5
k8=ln2[1]/sqrt_interval(8)[0]
b_bound=(k29+k25)/2
p_base=(41-5*rhi)*ln3[0]/73
m_base=(41+5*rlo)*ln3[0]/73
p_pen=F(3)*b_bound/(22-8*rhi)+(15+32*rhi)*k8/(2*73)
m_pen=F(3)*b_bound/(22+8*rlo)+(32*rhi-15)*k8/(2*73)
mu_plus_floor=g0_lo+p_base-errp-p_pen
mu_minus_floor=g0_lo+m_base-errm-m_pen
assert mu_plus_floor>F(1,40)
assert mu_minus_floor>F(1,4)

# Nonnegative bumps: original gamma and 25/29 gap terms >=0.
# |R(delta8)|<=1/2, k8<1/4, log3>109/100.
assert k8<F(1,4)
gap_floor=F(5,3)*(10*F(109,100)-32*F(1,4))/73
assert gap_floor>F(1,16)

print('PASS: collision windows',windows)
print('8-perturbation:',B8)
print('shift-plus:',plus,'shift-minus:',minus)
print('gamma diagonal lower:',float(g0_lo))
print('uniform mu+ >',float(mu_plus_floor),'> 1/40')
print('uniform mu- >',float(mu_minus_floor),'> 1/4')
print('nonnegative profile gap >',float(gap_floor),'> 1/16')
print('NO global Weil positivity, no RH, no Lean certification')
