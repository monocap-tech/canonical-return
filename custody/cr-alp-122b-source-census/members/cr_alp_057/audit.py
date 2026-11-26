#!/usr/bin/env python3
"""Global physical boundary profile remainder, in exact outward integers.
Analytic physical resolvent identities are premises; no Schur positivity test.
Run in predecessor archive layout. No numerical factorization.
"""
from pathlib import Path
import hashlib,json
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent
# Reuse exactly the inherited interval elementary constants, NOT its coefficient sweep.
source=PARENT/'cr_alp_056/extended.py'
s=source.read_text().split('nu={j:')[0].replace('ROOT=Path(__file__).resolve().parent', 'ROOT=Path(__file__).resolve().parent.parent/"cr_alp_056"')
exec(s)
ROOT=Path(__file__).resolve().parent
PARENT=ROOT.parent
PINS={1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}
def powi(x,n):
 y=I(1)
 while n:
  if n%2:y=y*x
  n//=2
  if n:x=x*x
 return y
def absI(x):return I(x.up())
def nr(row):return sum((absI(x).sq() for x in row),I(0)).sqrt()
def enc(x):return [str(x.l),str(x.h)]
def normbound(x):return F(x.h,S)
def seriesbound(p,M):
 # sum beta_m^(-p-1/2) <= (2M)^(-p-1/2) + integral_M^infinity (2x)^(-p-1/2) dx
 v=I(2*M);return (1+I(F(M,1)/(F(p)-F(1,2))))/(powi(v,p)*v.sqrt())
results={}
J=24
for bl in bank['blocks']:
 k=bl['parity'];ids=bl['ids'];n=len(ids);b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4;BM=(E-q)/2+(A if k==1 else -A)
 nu={i:i*NUF for i in ids}
 def c(i):v=nu[i];return 2*b0*v/(rt*(v.sq()+qt))
 raw=(PARENT/'cr_alp_049'/f'result_{k}.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==PINS[k]
 old=json.loads(raw);assert int(old['scale'])==S
 V=[[I(int(l),int(h),True) for l,h in row] for row in old['V']];ix=ids.index(k)
 f=[-c(k)/BM*x for x in V[ix]];d=[row[:] for row in V];d[ix]=[V[ix][j]-sum(c(i)/c(k)*V[ii][j] for ii,i in enumerate(ids) if i!=k) for j in range(n)]
 dn=[nr(row) for row in d];fn=nr(f)
 # The remainder is bounded AFTER retaining all the signed moment coefficients.
 # U bounds sum_i ||d_i|| nu_i^(2J+1)/sqrt(a). Use direct exact rational powers
 # for large terms; interval integer operations still retain all rounding.
 U=sum((v*powi(nu[i],2*J+1)/rt for i,v in zip(ids,dn)),I(0))
 coeff=[]
 for r in range(J):
  weights=[powi(nu[i],2*r+1)/rt for i in ids]
  moment=[sum((d[ii][j]*w for ii,w in enumerate(weights)),I(0)) for j in range(n)]
  coeff.append([v*b0*F(1,4**r) for v in f]);coeff.append([-v*b1*F(1,4**r)-(-1)**r*m for v,m in zip(f,moment)])
 bounds={}
 for M in (128,256,512,1024):
  # Avoid underflow of tiny interval powers by evaluating equivalent dimensionless expression.
  # For each beta>=2M, remainder <= C_odd (2M/beta)^(2J+1)+C_even (2M/beta)^(2J+2).
  B0=I(2*M);fac=1/(1-1/(4*B0.sq()));
  Co=fn*absI(b0)*fac/(4**J*powi(B0,2*J+1))
  Ce=fn*absI(b1)*fac/(4**J*powi(B0,2*J+2))+U/powi(B0,2*J+2)
  # sum of scaled powers times beta^-1/2; factor 2 for reflected exponential.
  def summ(p):return 2*(1+I(F(M,1)/(F(p)-F(1,2))))/B0.sqrt()
  err=Co*summ(2*J+1)+Ce*summ(2*J+2)
  bounds[str(M)]={'odd_amplitude_upper':str(normbound(Co)), 'even_amplitude_upper':str(normbound(Ce)), 'error_upper':str(normbound(err))}
  print(k,M,float(normbound(err)),flush=True)
 assert F(bounds['512']['error_upper'])<F(1,10**20)
 # Explicit signed boundary amplitudes for m<128, with m=0 ordinary-sine only.
 low=[]
 for m in range(128):
  beta=I(F(4*m+1,2));bb=beta.sq()
  w=[nu[i]/(rt*(bb+nu[i].sq())) for i in ids]
  row=[-sum((d[ii][j]*v for ii,v in enumerate(w)),I(0)) for j in range(n)]
  if m:row=[v+x*(beta*b0-b1)/(bb-qt) for v,x in zip(row,f)]
  low.append([enc(v) for v in row])
 results[str(k)]={'explicit_low_amplitudes':low, 'pole_coefficient':[enc(v) for v in f], 'source_sine_coefficients':[[enc(v) for v in row] for row in d], 'source_ids':ids, 'source_sha256':PINS[k], 'dimension':n,'pole_row_norm_upper':str(normbound(fn)),'sine_weighted_remainder_upper':str(normbound(U)), 'pole_left_upper':str(b0.up()),'pole_derivative_upper':str(b1.up()),'coefficients':[[enc(v) for v in row] for row in coeff], 'bounds':bounds}
(ROOT/'result.json').write_text(json.dumps({'scale':str(S),'J':J,'source_generator_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'sectors':results},indent=2)+'\n')
