"""Exact endpoint arithmetic for the CR-ALP-059 three-moment gate.
Inputs are outward decimal envelopes of CR-ALP-058 verification_result.json.
This script neither regenerates quadrature nor computes the second moment N.
"""
from fractions import Fraction as F
from itertools import product
ell,u,gamma=map(F,('0.255','5.505','2.88'))
assert ell+u==2*gamma
rows=[('odd',('0.15438056','0.15443126'),('-0.03046828','-0.03035098'),('0.12154757','0.12154758'), '3.04','0.85','0.1389','0.8064','0.00001050747311'),
      ('even',('0.13693461','0.13701749'),('-0.03281304','-0.03260290'),('0.08434342','0.08434343'), '2.83','1.20','0.0696','0.5741','0.00003200100411')]
for name,G,M,K,a,b,accept,reject,eps in rows:
    G,M,K=tuple(map(F,G)),tuple(map(F,M)),tuple(map(F,K))
    a,b,accept,reject,eps=map(F,(a,b,accept,reject,eps))
    def threshold(e,x):
        c=gamma**2-(2*x+e)*gamma+x*x+2*x*e
        d=2*gamma-2*x-e
        vals=[e*x*x*k-c*g-d*m for g,m,k in product(G,M,K)]
        return min(vals),max(vals)
    lo=threshold(ell,a)[0]; hi=threshold(u,b)[1]
    assert accept<lo and reject>hi
    # If an independently certified rational-shift proxy satisfies the same
    # 21/8 operator bound and norm <= 2/5, then this transport bound applies.
    s=F(21,8); delta=s*eps+F(1,10**13)
    transport=2*s*F(2,5)*delta+delta**2
    print(name,'accept N <=',accept,'reject comparison N >=',reject)
    print('  threshold endpoints:',lo,hi,'(decimal:',float(lo),float(hi),')')
    print('  conditional squared-norm transport:',float(transport))
print('PASS: exact endpoint thresholds; N not evaluated.')
