#!/usr/bin/env python3
"""Exact rational CR-ALP-041 comparison, conditional on CR-ALP-027/038/039/040.
Run in a directory containing the two pinned CR-ALP-040 interval files.
No floating point acceptance; no new continuum or Lean certification.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
from hashlib import sha256
import json, sys
ROOT=Path(__file__).resolve().parent
S=10**25; N=15000; a=F(11,10); sa=F(177,50); sig=F(1,10**14)
HF=ROOT/'CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv'
IF=ROOT/'CR-ALP-040_initial_rational_intervals_1e25.txt'
assert sha256(HF.read_bytes()).hexdigest()=='7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c'
assert sha256(IF.read_bytes()).hexdigest()=='6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9'
H={}; M={}; params={}
for ln in HF.read_text().splitlines():
 i,j,l,h=map(int,ln.split());assert l<=h and (i,j) not in H;H[i,j]=(F(l,S),F(h,S))
for ln in IF.read_text().splitlines():
 k,l,h=ln.split(); v=(F(int(l),S),F(int(h),S));assert v[0]<=v[1] and k not in params;params[k]=v
 if k.startswith('M_'):
  _,i,j=k.split('_');M[int(i),int(j)]=v
assert len(H)==98

def sym(A,i,j):
 l=max(A[i,j][0],A[j,i][0]);h=min(A[i,j][1],A[j,i][1]);assert l<=h;return l,h
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):
 v=[i*j for i in x for j in y];return min(v),max(v)
def scale(x,t):return mul(x,(t,t))
def absup(x):return max(abs(x[0]),abs(x[1]))
def sm(xs):
 xs=list(xs);return sum((x[0] for x in xs),F(0)),sum((x[1] for x in xs),F(0))
def dot(x,y):return sum((u*v for u,v in zip(x,y)),F(0))
def matvec(A,x):return [dot(r,x) for r in A]
def quadint(A,x,y):return sm(scale(A[i][j],x[i]*y[j]) for i in range(len(x)) for j in range(len(y)))
def basis(C):
 V=[]
 for i in range(len(C)):
  x=[F(int(j==i)) for j in range(len(C))]
  for v in V:
   den=dot(v,matvec(C,v));assert den>0
   t=dot(v,matvec(C,x))/den;x=[u-t*z for u,z in zip(x,v)]
  x=[F(round(z*10**18),10**18) for z in x];V.append(x)
 return [[z/max(map(abs,v)) for z in v] for v in V]

def residual(ids,x):
 k=ids[-1];z=x[-1];d={i:(v,v) for i,v in zip(ids[:-1],x[:-1])}
 d[k]=sm(scale(params[f'lam_{i}'],-v) for i,v in zip(ids[:-1],x[:-1]))
 D=absup(sm(mul(v,params[f'nu_{i}']) for i,v in d.items()))
 C={s:sum((absup(v)*params[f'nu_{i}'][1]**s for i,v in d.items()),F(0)) for s in [0,2,3]}
 A=(2*sa/a+F(201,200)*F(22,7)/(2*a))*(2*a/3)**2*D/isqrt(3*N**3)
 B=(sa*C[0]/(2*a)+8*sa*C[2]/a+F(201,200)*C[3]/a)*(2*a/3)**3/isqrt(5*N**5)
 return A+B+abs(z)*(F(1,25) if k==1 else F(1,19))

def floorq(x):return F(x.numerator*S//x.denominator,S)
def ceilq(x):return -floorq(-x)
blocks=[]
path=ROOT/'CR-ALP-041_energy_conditioned_result.json'
prior=json.loads(path.read_text()) if '--verify' in sys.argv else None
for b,ids in enumerate([[3,5,7,9,11,13,1],[4,6,8,10,12,14,2]]):
 K=[[add(sym(H,i,j),scale(sym(M,i,j),-sig)) for j in ids] for i in ids]
 V=basis([[(l+h)/2 for l,h in row] for row in K]) if prior is None else [[F(z) for z in v] for v in prior['blocks'][b]['basis_rows']]
 # Rows are coordinate vectors; their triangular shape proves independence.
 assert len(V)==7 and all(len(v)==7 and v[i]!=0 and all(z==0 for z in v[i+1:]) for i,v in enumerate(V))
 G=[[quadint(K,x,y) for y in V] for x in V]
 assert all(G[i][j]==G[j][i] for i in range(7) for j in range(7))
 budgets=[]; rho=F(0)
 for i,v in enumerate(V):
  d=floorq(G[i][i][0]-sum(absup(G[i][j]) for j in range(7) if i!=j))
  r=ceilq(residual(ids,v));assert d>0 and r>=0
  rho+=r*r/d;budgets.append({'d':str(d),'r':str(r),'ratio_upper':str(ceilq(r*r/d))})
 upper=F(104,1000) if b==0 else F(79,1000)
 assert rho<upper<F(1,8)
 print('PASS',ids,'7 positive transformed energy budgets; rho <',upper,'< 1/8')
 blocks.append({'indices':ids,'basis_rows':[[str(z) for z in v] for v in V],'budgets':budgets,'rho_upper':str(upper)})
result={'stage':'CR-ALP-041','status':'PASS conditional exact rational E14 plus entire Q15000 mixed Schur','sigma':str(sig),'tail_margin':'1/8','blocks':blocks,'premises':'CR-ALP-027 positive tail; CR-ALP-038 complete signed interior residual formula; CR-ALP-039 anchor residual bounds; CR-ALP-040 original source lower bank. No new Lean certification.','not_claimed':'intermediate 14986 modes, full collar, all support radii, RH'}
if prior is not None:
 assert result==prior
 print('PASS stored basis and budgets independently rechecked against pinned inputs')
else:path.write_text(json.dumps(result,indent=2)+'\n')
