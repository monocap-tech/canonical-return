"""Exact finite checks of the CR-ALP-060 candidate matrix-moment inequalities.
No actual CR-ALP second moments or original Weil entries are evaluated.
"""
from fractions import Fraction as F
from itertools import combinations

def tr(A): return [list(r) for r in zip(*A)]
def mul(A,B): return [[sum((a*b for a,b in zip(r,c)),F(0)) for c in tr(B)] for r in A]
def scale(k,A): return [[k*x for x in r] for r in A]
def add(A,B): return [[x+y for x,y in zip(r,s)] for r,s in zip(A,B)]
def sub(A,B): return add(A,scale(-1,B))
def diag(*xs): return [[F(x) if i==j else F(0) for j in range(len(xs))] for i,x in enumerate(xs)]
def gram(R,D): return mul(tr(R),mul(D,R))
def inv2(A):
    a,b=A[0]; c,d=A[1]; det=a*d-b*c
    assert det!=0
    return scale(1/det,[[d,-b],[-c,a]])
def psd2(A):
    a,b=A[0]; c,d=A[1]
    return b==c and a>=0 and d>=0 and a*d>=b*c

def correction(B,A,Y):
    BY=mul(B,Y)
    return sub(add(BY,tr(BY)),gram(Y,A))
def moments(ts,R):
    return [gram(R,diag(*(t**k for t in ts))) for k in (0,1,2,-1)]

ell,u=F(51,200),F(1101,200)
# Full residual range, nonorthogonal source coordinates; exact recovery.
R=[[F(1),F(1)],[F(0),F(1)]]
G,M1,M2,J=moments((F(1),F(4)),R)
A=sub(M2,scale(ell,M1)); B=sub(M1,scale(ell,G))
C=sub(scale(u,M1),M2); D=sub(scale(u,G),M1)
Y=mul(inv2(A),B); Z=mul(inv2(C),D)
U=scale(1/ell,sub(G,correction(B,A,Y)))
V=scale(1/u,add(G,correction(D,C,Z)))
assert U==V==J
K=add(J,scale(F(1,100),G)); assert psd2(sub(K,U))
# Any common scalar parameter a would require |1-t/a|<1/10 at both
# t=1 and t=4, since its allowed squared error is strictly <1/100.
for t in (F(1),F(4)):
    assert F(1,100)*ell*t/(t-ell)<F(1,100)
assert F(11,40)<F(9,10) # disjoint necessary intervals for 1/a
print('PASS exact nonorthogonal 2D example; no common scalar parameter can pass')
# Proper residual range, nontrivial interval. Verify optimal and arbitrary trials.
R=[[F(1),F(2)],[F(3),F(-1)],[F(2),F(1)]]
G,M1,M2,J=moments((F(1,2),F(2),F(5)),R)
A=sub(M2,scale(ell,M1)); B=sub(M1,scale(ell,G))
C=sub(scale(u,M1),M2); D=sub(scale(u,G),M1)
assert psd2(A) and psd2(C)
Y=mul(inv2(A),B); Z=mul(inv2(C),D)
for y,z in [(Y,Z),(diag(F(1,3),F(1,4)),diag(F(1,2),F(1,5)))]:
    U=scale(1/ell,sub(G,correction(B,A,y)))
    V=scale(1/u,add(G,correction(D,C,z)))
    assert psd2(sub(U,J)) and psd2(sub(J,V))
print('PASS proper 2D residual range in 3D, optimal and arbitrary rational trials')
# Endpoint-singular A/C: arbitrary trial formula needs no inverse at all.
G,M1,M2,J=moments((ell,u),diag(1,1))
A=sub(M2,scale(ell,M1)); B=sub(M1,scale(ell,G))
C=sub(scale(u,M1),M2); D=sub(scale(u,G),M1)
Y=Z=diag(1/ell,1/u)
assert scale(1/ell,sub(G,correction(B,A,Y)))==J
assert scale(1/u,add(G,correction(D,C,Z)))==J
print('PASS singular endpoint cases without matrix inversion')
print('NO physical moment bank, Lean proof, or original Weil acceptance claimed.')
