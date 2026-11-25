#!/usr/bin/env python3
"""Independent standard-library replay of saved moment comparisons; not an analytic generator."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
BANK={1:'65b38a4cd52c5370aca1fe46daf1aefd54846b6d57a71b0b76858330720aec43',2:'41f279b5f57197d511dca4590ae25bcc7e3d0cf9c74dc31ec6ac22d469670771'}
MOM={1:'dc9cf9cb76a5a70aa678464115a8282c312f26e9d4ed94e2848fb8d8e09cc88b',2:'deef59f11b8ed595e1283f4c7acd214f7fdca46881ccfb3bb02237f02ec124e3'}
PRE={1:'a133a5bf9a8229b16e6459aee952bdca677567ed9df425e7430c466bd87e0877',2:'ee7ae0d12d7beefd2ecb08c3825ca5f7ab9ea92da2fba2f28b836094c212c4f3'}
def floor(f):return f.numerator//f.denominator
def ceil(f):return -floor(-f)
def load(A):return [[tuple(map(int,x)) for x in row] for row in A]
def ld(A,shift,S):
 n=len(A);C=[];rs=[];sh=shift*S;assert sh.denominator==1
 for i,row in enumerate(A):
  cc=[];rr=0
  for j,(l,h) in enumerate(row):
   assert l<=h
   if i==j:l-=int(sh);h-=int(sh)
   m=(l+h)//2;cc.append(m);rr+=max(m-l,h-m)
  C.append(cc);rs.append(rr)
 for i in range(n):
  for j in range(i):assert C[i][j]==C[j][i]
  C[i][i]-=max(rs)
 prev=1
 for k in range(n):
  pivot=C[k][k];assert pivot>0,(k,pivot)
  for i in range(k+1,n):
   for j in range(i,n):
    q,r=divmod(C[i][j]*pivot-C[i][k]*C[k][j],prev);assert r==0;C[i][j]=C[j][i]=q
  prev=pivot
 return n
def sub(A,B):return [[(l-v,h-u) for (l,h),(u,v) in zip(a,b)] for a,b in zip(A,B)]
def quad(A,z,S):
 lo=hi=0
 for i,zi in enumerate(z):
  for j,zj in enumerate(z):
   l,h=A[i][j];v=zi*zj;lo+=min(l*v,h*v);hi+=max(l*v,h*v)
 return F(lo,S*sum(x*x for x in z)),F(hi,S*sum(x*x for x in z))
def upper(A,S):return F(max(sum(max(abs(l),abs(h)) for l,h in row) for row in A),S)
ell,u,gamma=F(51,200),F(1101,200),F(72,25)
result={}
for k in (1,2):
 raw=(ROOT/f'bank_{k}.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==BANK[k];b=json.loads(raw)
 raw=(ROOT/f'moments_{k}.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==MOM[k];m=json.loads(raw)
 assert b['predecessor_sha256']==PRE[k] and m['input_sha256']==BANK[k]
 assert b['trial']==600 and b['band_end']==850 and b['parity']==m['parity']==k
 n=44-k;assert len(b['R'])==382+k and len(b['P'])==382+k
 S=int(m['scale']);K,G,M=map(load,[m['K'],m['G'],m['M']]);assert len(K)==n
 z=[]
 for row in b['witness']['z']:
  l,h=map(int,row[0]);assert l==h;z.append(l)
 # Reconstruct the affine penalty independently, with outward rational rounding.
 A=[[(floor((gamma*gl-mh)/(ell*u)),ceil((gamma*gh-ml)/(ell*u))) for (gl,gh),(ml,mh) in zip(gr,mr)] for gr,mr in zip(G,M)]
 tests={}
 for name,t in m['tests'].items():
  J=A if name=='affine' else load(t['penalty'])
  C=sub(K,J);shift=F(t['shift']);assert ld(C,shift,S)==n
  lo,hi=quad(C,z,S);assert lo>0
  tests[name]={'positive_pivots':n,'shift':str(shift),'frozen_witness_lower':str(lo),'frozen_witness_upper':str(hi),'penalty_row_upper':str(upper(J,S))}
  print(k,name,'pivots',n,'shift',shift,'frozen witness lower',float(lo),flush=True)
 # Full inverse on the projected residual, and a conditional tolerance for missing residual.
 a,sigma,eps,sqrt_bound=(F(7,25),F(1,100),F(1,250),F(11,10)) if k==1 else (F(4,25),F(1,32),F(1,60),F(4,5))
 assert upper(A,S)<a and sqrt_bound**2>=a/ell
 remaining=sigma-2*sqrt_bound*eps-eps*eps/ell
 assert remaining>F(1,1000)
 result[str(k)]={'bank_sha256':BANK[k],'moments_sha256':MOM[k],'tests':tests,'conditional_missing_residual_norm':str(eps),'conditional_schur_margin_lower':str(remaining),'residual_gram_row_upper':str(upper(G,S))}
 print('conditional',k,'missing norm',eps,'margin',float(remaining),flush=True)
(ROOT/'verification_result.json').write_text(json.dumps(result,indent=2)+'\n')
