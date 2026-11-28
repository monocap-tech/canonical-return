#!/usr/bin/env python3
"""CR-ALP-028 Fraction gates for finite Fourier interval residual envelopes.

The continuum inequalities (integration by parts, Fourier transform,
trigamma series, midpoint error and Schur) are proved in the companion note.
Here we check rational constants, 15k/1e6 bounds and abstract Schur algebra.
No actual original Weil Gamma entries/eigenvalue have been certified.
"""
from fractions import Fraction as F
from math import factorial
A=F(11,10)
N=15_000
T=1_000_000

# True Gamma derivative ≤ (1/2)*(16 + sum_{n>=1} n^-2) <9;
# original prime derivative ≤2*(2a)*S, S<177/50.
assert F(9)+4*A*F(177,50)<26
# Native gammaBracket in [log(1+|r|)-8, log(1+|r|)+3], and
# genuine 7-prime multiplier absolute <=177/25, hence |m|≤log+16.
assert 8+F(177,25)<16
assert 3+F(177,25)<16

# Elementary bounds for exact physical anchor functions cosh(x/2), sinh(x/2).
z=F(11,20)
# Positive series; tail bounded by first missing term / (1 - next-ratio).
for typ,cut,desired in [('cosh',5,F(6,5)),('sinh',5,F(3,5))]:
    offset=0 if typ=='cosh' else 1
    part=sum((z**(2*k+offset)/factorial(2*k+offset)
              for k in range(cut)),F(0))
    exponent=2*cut+offset
    nextterm=z**exponent/factorial(exponent)
    ratio=z*z/F((exponent+1)*(exponent+2))
    assert part+nextterm/(1-ratio)<desired
assert A**3/6 + A**5/120 > F(23,100)  # sinh(a)-a>23/100

nu1lo=F(3)/(2*A);nu1hi=F(22,7)/(2*A)
nu2lo=2*nu1lo;nu2hi=2*nu1hi
# |gamma_1|=|<e1,cosh>|/||cosh||^2 and similarly gamma_2.
# sqrt(a)>1, ||cosh||²=a+sinh(a)>2a,
# ||sinh||²=sinh(a)-a>23/100; bound cosh(a/2),sinh(a/2).
gamma1_upper=(2*nu1hi*F(6,5)/(nu1lo*nu1lo+F(1,4)))/(2*A)
gamma2_upper=(2*nu2hi*F(3,5)/(nu2lo*nu2lo+F(1,4)))/F(23,100)
K1_upper=2+gamma1_upper*(4*F(6,5)-2)
K2_upper=4+gamma2_upper*4*F(3,5)
assert K1_upper<5 and K2_upper<9
print('PASS: exact physical anchor Fourier BV constants K1<5, K2<9')

# For w_j=e_j-lambda_j e_anchor, lambda_j*nu_anchor<=nu_j,
# so beta_j=nu_j+lambda_j*nu_anchor<=2nu_j.
# pi<22/7; sum beta_j² <=(pi/a)² sum_{j=3}^N j².
ss=F(N*(N+1)*(2*N+1),6)-5
B2_upper=(F(22,7)/A)**2*ss
# For |r|>T>=2nu_N, |w_jhat|<=8 beta_j/(3sqrt(a)*r²).
assert T >= N*F(22,7)/A
j0_w=F(64)*B2_upper/(27*A*3*T**3) # coefficient l2 norm -> Fourier L2 tail
j0_u=F(106)/(3*T)                    # K1²+K2²<25+81=106
# log(T)<14 since T=10^6, log 10<7/3, and |m(r)|<log(r)+17.
# z=log T +17 < 31 for T=10^6.
zup=F(31)
q3=zup*zup+F(2,3)*zup+F(2,9)
q1=zup*zup+2*zup+2
j2_w=j0_w*q3
j2_u=j0_u*q1
# Full basis p=p_w+p_u: |p_hat|²≤2|p_w_hat|²+2|p_u_hat|².
K0w=2*j0_w;K0u=2*j0_u
K2w=2*j2_w;K2u=2*j2_u
assert K2w<F(1,70) and K2u<F(3,40)
assert K0w<F(1,70_000) and K0u<F(1,14_000)
print('PASS at T=10^6, 15k: matrix envelope K2 interior <1/70, boundary <3/40')
print('PASS: K0 interior <1/70000, boundary <1/14000 (RAW basis coordinates)')

# Matrix-free residual gate:
# R_C=||T p-V Cx||² >=||Q_N Pi_0 1_I T p||²=B_N(x)
# (from orthogonal projections and V Cx in E_N).
# With tail coercivity mu=1/4, H-4R_C>=0 is SUFFICIENT.
# Toy 2D exact certified Schur / residual identity.
p=(F(1),F(1));q=(F(1),F(-1))
def dot(x,y):return sum((a*b for a,b in zip(x,y)),F(0))
Tdiag=(F(2),F(3))
Tp=tuple(Tdiag[i]*p[i] for i in range(2))
H=dot(Tp,p);M=dot(p,p);cross=dot(Tp,q)
C=H/M
vcp=tuple(C*x for x in p)
R=dot(tuple(Tp[i]-vcp[i] for i in range(2)),tuple(Tp[i]-vcp[i] for i in range(2)))
B=cross**2/dot(q,q)
assert R==B==F(1,2)
assert H-4*R==3>0
assert H-4*B==3>0
# Non-optimal C shows R_C>=B and may make the sufficient test fail.
Cbad=F(0)
Rb=dot(Tp,Tp)
assert Rb>R and H-4*Rb<0
print('PASS: exact physical residual upper-bound Schur logic (success and conservatism)')
print('CR-ALP-028 RESIDUAL ENVELOPE ARITHMETIC PASS; actual signed 15k H/RC OPEN')
