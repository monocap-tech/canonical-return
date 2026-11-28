#!/usr/bin/env python3
"""CR-INT-T1D-B: exact rational checks of profile-uniform bounds and 25/29 Gram update.

Mathematical support/correlation and Gamma integrals are proved in
research/CR-INT-T1D-B.md; this audit checks exact algebra and rational estimates.
Not a Lean proof or a global Weil positivity certificate.
"""
from fractions import Fraction as F
from math import factorial
import sympy as sp

# log(2) = 2 atanh(1/3); rational lower and geometric-tail upper bounds.
l2_lo = 2*sum((F(1,3)**(2*k+1)/F(2*k+1) for k in range(3)),F(0))
l2_partial=2*sum((F(1,3)**(2*k+1)/F(2*k+1) for k in range(5)),F(0))
l2_hi=l2_partial+2*F(1,3)**11/(11*(1-F(1,9)))
assert l2_lo == F(842,1215) and F(69,100)<l2_lo
assert l2_hi < F(347,500)
# log(3)=2atanh(1/2)
assert 2*sum((F(1,2)**(2*k+1)/F(2*k+1) for k in range(3)),F(0)) > F(109,100)
# log(4*pi)<127/50: pi<22/7 and exp(127/50)>88/7, rational Taylor.
assert sum((F(127,50)**k/F(factorial(k)) for k in range(8)),F(0)) > F(88,7)
# gamma < H20 - log20 - 1/42, and log20 > 4log2 + 2/9.
h20=sum((F(1,k) for k in range(1,21)),F(0))
assert h20-(4*l2_lo+F(2,9))-F(1,42)<F(29,50)
assert F(127,50)+F(29,50)<F(25,8)
# log(50/3)>14/5, log(9/8)<3/25, and exact threshold ordering.
assert 4*F(69,100)+F(1,25)==F(14,5)
x=F(3,25)
assert 1+x+x*x/2+x*x*x/6>F(9,8)
assert F(29,27)<F(27,25)<F(9,8)<F(29,27)**2
# k29=log29/sqrt29<16/25 and k25=log5/5<1/3.
assert 5*F(347,500)-F(3,32)<F(17,5)
assert F(107,20)**2<29
assert F(17,5)/F(107,20)<F(16,25)
assert 2*F(347,500)+F(1,4)<F(5,3)
assert (F(16,25)+F(1,3))/2<F(1,2)
# w(k log3 - ell) upper bounds for ell<=3/25.
assert F(3,5)*F(50,47)/(1-F(25,19)/9)<F(3,4)
assert F(1,3)*F(50,47)/(1-F(25,19)/81)<F(3,8)
assert F(1,5)*F(50,47)/(1-F(25,19)/729)<F(1,4)
# Gamma diagonal and weighted off-diagonal coefficient bounds.
gamma_diag_lo=F(14,5)-F(25,8)-F(3,25)+F(47,50)*F(69,200)
gamma_p=F(25,73)*F(3,4)+F(71,73)*F(3,8)+F(27,73)*F(1,4)
gamma_m=F(108,73)*F(3,4)+F(41,73)*F(3,8)+F(13,146)*F(1,4)
assert gamma_diag_lo==-F(1207,10000)
assert gamma_p==F(417,584)
assert gamma_m==F(98,73)
mu_p=gamma_diag_lo+(41-5*F(7,4))*F(109,100)/73-gamma_p*F(3,25)-F(3,16)
mu_m=gamma_diag_lo+(41+5*F(5,3))*F(109,100)/73-gamma_m*F(3,25)-F(9,212)
assert mu_p==F(63989,730000)>F(2,25)
assert mu_m==F(47864051,116070000)>F(2,5)
# Compression of the outer-cell (3d) prime 25 and 29 samples.
s=4/sp.sqrt(3)
X=sp.Matrix([[1,0],[-s,1],[1,-s],[0,1]])
D=sp.zeros(4); D[0,3]=D[3,0]=1
assert sp.simplify(X.T*D*X)==sp.Matrix([[0,1],[1,0]])
J,b25,b29=sp.symbols("J b25 b29",real=True,positive=True)
M=sp.simplify(J*X.T*X)
dH=-J*(b25+b29)*(X.T*D*X)
q,c=dH[0,0],dH[0,1]
u,v=M[0,0],M[0,1]
dplus=sp.simplify((q+c)/(u+v))
dminus=sp.simplify((q-c)/(u-v))
assert sp.simplify(dplus+3*(b25+b29)/(22-8*sp.sqrt(3)))==0
assert sp.simplify(dminus-3*(b25+b29)/(22+8*sp.sqrt(3)))==0
assert sp.simplify(dminus-dplus-sp.Rational(33,73)*(b25+b29))==0
print("PASS: delta29 < delta25 < delta8, 25=5^2 coefficient, two-shift Gram update.")
print("PASS: outward rational bounds mu_plus >",mu_p, "; mu_minus >",mu_m)
print("Selected uniform Weil margin > 2/25 for support width < log(9/8).")
print("These are finite-pencil assertions only; no global Weil positivity/RH proof.")
