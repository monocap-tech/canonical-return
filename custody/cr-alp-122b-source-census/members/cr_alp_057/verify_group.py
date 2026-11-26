#!/usr/bin/env python3
"""Independent standard-library verifier of saved dyadic proxy coefficient intervals.
Exact rational group weights via raw power sums (not the generator's centered recurrence).
Checks saved boundary remainder budgets with conservative elementary source bounds.
Does not regenerate the inherited continuum coefficients or evaluate Schur moments.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
import json,hashlib
ROOT=Path(__file__).resolve().parent
raw=(ROOT/'result.json').read_bytes();a=json.loads(raw);b=json.loads((ROOT/'grouped_proxy.json').read_text());S=int(a['scale']);assert S==int(b['scale'])==10**100;assert hashlib.sha256(raw).hexdigest()==b['source_result_sha256'];K=b['K'];R=b['R'];END=b['end'];assert (K,R,END)==(64,40,62)
def up(x):return F(-((-x.numerator*10**90)//x.denominator),10**90)
def norm(row):return F(isqrt(sum(max(abs(int(l)),abs(int(h)))**2 for l,h in row))+1,S)
def radius(row):return F(isqrt(sum(max((int(l)+int(h))//2-int(l),int(h)-(int(l)+int(h))//2)**2 for l,h in row))+1,S)
def invsqrt(x):return F(10**40,isqrt(x*10**80))
def exact_centered(M,n):
 P=[M]
 for j in range(1,n+1):
  z=M**(j+1)-sum(comb(j+1,k)*P[k] for k in range(j));assert z%(j+1)==0;P.append(z//(j+1))
 C=[]
 for j in range(n+1):
  z=sum(comb(j,k)*2**k*(1-M)**(j-k)*P[k] for k in range(j+1));assert z>=0 and (j%2==0 or z==0);C.append(z)
 return C

def weights_exact(M):
 C=exact_centered(M,K+R-2);d=6*M-1;pw=[d**j for j in range(K+R+48)];rows=[]
 # All weights refer to the unscaled coefficient rows in result.json.
 for r in range(R):
  w=[]
  for p in range(1,49):
   z=sum((-1)**(k+r)*comb(p+k-1,k)*C[k+r]*2**(k+r+p)*pw[K-1-k+48-p] for k in range(K))
   w.append(z)
  rows.append((w,pw[K-1+r+48]))
 return rows
coeff={};bounds={};rounding={};errors={};checks={}
for k,d in a['sectors'].items():
 coeff[k]=[[tuple(map(int,x)) for x in row] for row in d['coefficients']];bounds[k]=[norm(row)/256**p for p,row in enumerate(d['coefficients'],1)];rounding[k]=sum(up(2*radius(row)*invsqrt(4*m+1)) for m,row in enumerate(d['explicit_low_amplitudes']));errors[k]=F(0);checks[k]=0
 assert d['explicit_low_amplitudes']==b['sectors'][k]['low_amplitudes']
for j in range(7,END):
 M=2**j;W=weights_exact(M);bc=F(6*M-1,2)
 for k,d in b['sectors'].items():
  g=d['groups'][j-7];assert int(g['M'])==M and F(g['center'])==bc
  for (w,den),row in zip(W,g['amplitudes']):
   for i,(l,h) in enumerate(row):
    lo=hi=0
    for p,x in enumerate(w):
     u,v=coeff[k][p][i];lo+=x*(u if x>=0 else v);hi+=x*(v if x>=0 else u)
    assert int(l)*den<=lo<=hi<=int(h)*den,(j,k,i)
    checks[k]+=1
  rounding[k]+=up(2*sum(radius(row) for row in g['amplitudes'])*invsqrt(6*M-1))
  A=bounds[k];v=sum(x*F(256,2*M)**p for p,x in enumerate(A,1));e=F(0)
  for p,x in enumerate(A,1):
   rho=F(p+K,3*(K+1));assert rho<1
   e+=x*(256/bc)**p*F(comb(p+K-1,K),3**K)/(1-rho)
  errors[k]+=up(2*M*invsqrt(4*M)*(v*F(1,2**R)+e*(1+F(1,2**R))))
 if j%8==7:print('independent groups through',j,flush=True)
for k,d in a['sectors'].items():
 M=2**END;A=bounds[k];tail=up(sum(2*x*F(256,2*M)**p*(1+F(M)/(F(p)-F(1,2)))*invsqrt(2*M) for p,x in enumerate(A,1)))
 s=b['summary'][k];assert rounding[k]<=F(s['midpoint_coefficient_error_upper']);assert errors[k]<=F(s['group_error_upper']);assert tail<=F(s['omitted_profile_tail_upper']);assert sum(A)<=F(s['signed_scaled_norm_sum_upper'])
 # Verify expansion remainder independently using nu_i < 10i/7, sqrt(a)>1,
 # and the archived outward source-coordinate coefficient rows.
 J=a['J'];M=128;B0=2*M
 U=sum(norm(row)*F(10*i,7)**(2*J+1) for i,row in zip(d['source_ids'],d['source_sine_coefficients']))
 fn=norm([*d['pole_coefficient']]);fac=1/(1-F(1,4*B0**2));odd=fn*F(d['pole_left_upper'])*fac/F(4**J*B0**(2*J+1));even=fn*F(d['pole_derivative_upper'])*fac/F(4**J*B0**(2*J+2))+U/B0**(2*J+2)
 def summ(p):return 2*(1+F(M)/(F(p)-F(1,2)))*invsqrt(B0)
 expansion=odd*summ(2*J+1)+even*summ(2*J+2);assert expansion<F(1,10**18)
 total=expansion+rounding[k]+errors[k]+tail;assert total<F(5,10**9)
 print(k,'coefficient intervals verified',checks[k],'independent expansion upper',float(expansion),'global proxy error upper',float(total),flush=True)
print('PASS: all saved dyadic coefficient intervals and complete error budgets; continuum source entries inherited',flush=True)
