#!/usr/bin/env python3
"""Finite exponential-polynomial proxy, no high-sine cutoff.
Dyadic power sums and all coefficient intervals computed with exact integers.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
raw=(ROOT/'result.json').read_bytes(); data=json.loads(raw);S=int(data['scale']);B=256;K=64;R=40;END=62
assert S==10**100 and data['J']==24

def ceildiv(n,d):return -((-n)//d)
def up90(x):return F(ceildiv(x.numerator*10**90,x.denominator),10**90)
def encl(n,d):return (n*S//d,ceildiv(n*S,d))
def divi(v,d):return v[0]//d,ceildiv(v[1],d)
def prod(x,y):
 z=[u*v for u in x for v in y];return min(z)//S,ceildiv(max(z),S)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def smul(x,n):return (x[0]*n,x[1]*n) if n>=0 else (x[1]*n,x[0]*n)
def norm_upper(row):return F(isqrt(sum(max(abs(l),abs(h))**2 for l,h in row))+1,S)
def radius_norm(row):return F(isqrt(sum(max((l+h)//2-l,h-(l+h)//2)**2 for l,h in row))+1,S)
def invsqrt_upper(x):
 d=10**40
 return F(d,isqrt(x*d*d))
def group_error(a,M):
 bc=F(6*M-1,2);bm=2*M
 trueA=sum(v*F(B,bm)**p for p,v in enumerate(a,1));errA=F(0)
 for p,v in enumerate(a,1):
  rho=F(p+K,3*(K+1));assert rho<1
  tail=F(comb(p+K-1,K),3**K)/(1-rho)
  errA+=v*(B/bc)**p*tail
 return 2*M*invsqrt_upper(2*bm)*(trueA*F(1,2**R)+errA*(1+F(1,2**R)))
def tail_error(a,M):
 return sum(2*v*F(B,2*M)**p*(1+F(M)/(F(p)-F(1,2)))*invsqrt_upper(2*M) for p,v in enumerate(a,1))
def moment_bounds(M):
 # Integer central moments via sinh(M t)/sinh t. Odd moments vanish.
 C=[0]*(K+R-1);C[0]=M
 for j in range(2,len(C),2):
  v=M**(j+1)-sum(comb(j+1,k)*C[k] for k in range(0,j,2));assert v%(j+1)==0;C[j]=v//(j+1);assert C[j]>=0
 den=6*M-1
 return [encl(c*2**j,den**j) for j,c in enumerate(C)]
def weights(M):
 H=moment_bounds(M);result=[]
 for r in range(R):
  row=[]
  for p in range(1,49):
   z=(0,0)
   for k in range(K):
    z=add(z,smul(H[k+r],(-1)**(k+r)*comb(p+k-1,k)))
   row.append(prod(z,encl((2*B)**p,(6*M-1)**p)))
  result.append(row)
 return result

def mr(rows):
 m=np.array([[(l+h)//2 for l,h in row] for row in rows],dtype=object)
 e=np.array([[max((l+h)//2-l,h-(l+h)//2) for l,h in row] for row in rows],dtype=object)
 return m,e
def mult(A,D):
 a,r=mr(A);b,s=mr(D);m=a@b;e=abs(a)@s+r@abs(b)+r@s
 return [[((int(x)-int(y))//S,ceildiv(int(x)+int(y),S)) for x,y in zip(row,rad)] for row,rad in zip(m,e)]
def string(rows):return [[[str(l),str(h)] for l,h in row] for row in rows]
scaled={}; out={};a={};rounding={};errors={}
for k,d in data['sectors'].items():
 scaled[k]=[[divi(tuple(map(int,v)),B**p) for v in row] for p,row in enumerate(d['coefficients'],1)]
 # use exact original coefficient norms for truncation budgets, not the rounded scaled bank
 a[k]=[norm_upper([tuple(map(int,x)) for x in row])/B**p for p,row in enumerate(d['coefficients'],1)]
 rounding[k]=F(0);errors[k]=F(0)
 low=[[tuple(map(int,x)) for x in row] for row in d['explicit_low_amplitudes']]
 for m,row in enumerate(low):rounding[k]+=up90(2*radius_norm(row)*invsqrt_upper(4*m+1))
 out[k]={'low_amplitudes':string(low),'groups':[]}
 print(k,'signed scaled coefficient norm sum',float(sum(a[k])),flush=True)
for j in range(7,END):
 M=2**j;W=weights(M)
 for k in out:
  C=mult(W,scaled[k]);out[k]['groups'].append({'M':str(M),'center':str(F(6*M-1,2)),'amplitudes':string(C)})
  rounding[k]+=up90(2*sum(radius_norm(row) for row in C)*invsqrt_upper(6*M-1))
  errors[k]+=up90(group_error(a[k],M))
 if j%8==7:print('dyadic group',j,'complete',flush=True)
summary={}
for k,d in data['sectors'].items():
 expand=F(d['bounds']['128']['error_upper']);tail=up90(tail_error(a[k],2**END));total=expand+tail+rounding[k]+errors[k];assert total<F(1,10**8)
 summary[k]={'expansion_error_upper':str(expand),'group_error_upper':str(errors[k]),'omitted_profile_tail_upper':str(tail),'midpoint_coefficient_error_upper':str(rounding[k]),'total_operator_error_upper':str(total),'signed_scaled_norm_sum_upper':str(up90(sum(a[k])))}
 print(k,'group error',float(errors[k]),'profile tail',float(tail),'midpoint rounding',float(rounding[k]),'TOTAL',float(total),flush=True)
result={'source_result_sha256':hashlib.sha256(raw).hexdigest(),'scale':str(S),'K':K,'R':R,'end':END,'basis_terms_per_parity':128+(END-7)*R,'summary':summary,'sectors':out}
(ROOT/'grouped_proxy.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
(ROOT/'grouped_summary.json').write_text(json.dumps({k:v for k,v in result.items() if k!='sectors'},indent=2)+'\n')
