#!/usr/bin/env python3
"""Independent integer check of frozen interval comparisons, not analytic entry reconstruction."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
BANK={1:'7a6b1d08bc966c12d189d7387e8699a3e3fab516b35f887184c90d0085d746a4',2:'0133407e41b7cf6806366f2ea646acf143fde3289e6377fb14e8d46b21a38a2a'}
def ld(A,shift,scale):
 n=len(A);C=[];radius=[];sh=shift*scale;assert sh.denominator==1
 for i,row in enumerate(A):
  c=[];r=0
  for j,(lo,hi) in enumerate(row):
   l,h=int(lo),int(hi);assert l<=h
   if i==j:l-=int(sh);h-=int(sh)
   m=(l+h)//2;c.append(m);r+=max(m-l,h-m)
  C.append(c);radius.append(r)
 for i in range(n):
  for j in range(i):assert C[i][j]==C[j][i]
  C[i][i]-=max(radius)
 prev=1
 for k in range(n):
  pivot=C[k][k];assert pivot>0
  for i in range(k+1,n):
   for j in range(i,n):
    u,r=divmod(C[i][j]*pivot-C[i][k]*C[k][j],prev);assert r==0;C[i][j]=C[j][i]=u
  prev=pivot
 return n

def rawq(A,z,scale,zscale):
 low=high=0
 for i,zi in enumerate(z):
  for j,zj in enumerate(z):
   l,h=map(int,A[i][j]);p=zi*zj;low+=min(p*l,p*h);high+=max(p*l,p*h)
 return F(low,scale*zscale*zscale),F(high,scale*zscale*zscale)

x=F(192,175);T=[F(1),x]
for m in range(2,13):T.append(2*x*T[-1]-T[-2])
assert 1/T[12]**2<F(1,8000)
result={}
for parity in (1,2):
 raw=(ROOT/f'bank_{parity}.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==BANK[parity];b=json.loads(raw);bs=int(b['scale']);p=(ROOT/f'polynomial_{parity}.json');rawp=p.read_bytes();p=json.loads(rawp);s=int(p['scale']);assert p['input_sha256']==BANK[parity]
 n=44-parity;assert len(b['band'])==256+parity;assert len(b['K'])==n;d=bs//s;assert bs==s*d
 G=max(sum(max(abs(int(l)),abs(int(h))) for l,h in row) for row in b['G']);assert F(G,bs)<F(3,4)
 z=[]
 for l,h in b['witness']['z']:assert l==h;z.append(int(l))
 result[str(parity)]={'bank_sha256':BANK[parity],'polynomial_sha256':hashlib.sha256(rawp).hexdigest(),'projected_residual_row_bound':str(F(G,bs)),'tests':{}}
 for m in (4,8,12):
  t=p['tests'][str(m)];eta=1/T[m];assert F(t['eta'])==eta;f=1/(1-eta);cmp=t['comparison'];J=t['J']
  # independently check interval subtraction against original high-precision K and saved moments J
  for i in range(n):
   for j in range(n):
    kl,kh=map(int,b['K'][i][j]);jl,jh=map(int,J[i][j]);lo,hi=map(int,cmp[i][j]);assert F(lo,s)<=F(kl,bs)-f*F(jh,s) and F(hi,s)>=F(kh,bs)-f*F(jl,s)
  shift=F(1,4 if parity==1 and m>=8 else 8);assert F(t['shift'])==shift;assert ld(cmp,shift,s)==n
  lo,hi=rawq(cmp,z,s,bs);assert lo>0;assert lo>=F(t['witness_lower']) and hi<=F(t['witness_upper']);
  if m==12:assert lo>F(1,2)
  print('parity',parity,'degree',m,'pivots',n,'margin',shift,'witness lower',float(lo),'PASS')
  result[str(parity)]['tests'][str(m)]={'pivots':n,'margin':str(shift),'witness_lower':str(lo),'witness_upper':str(hi)}
 margin=F(1,4 if parity==1 else 8);target=F(1,4) if parity==1 else F(7,40);remaining=margin-F(200,51)*(target**2+F(3,32000));assert remaining==F(37,8160)
 result[str(parity)]['conditional_outside_residual_norm_target']=str(target);result[str(parity)]['conditional_remaining_schur_margin']=str(remaining)
print('PASS: both finite polynomial comparisons and explicit conditional outside-band targets; no infinite leakage computation')
(ROOT/'verification_result.json').write_text(json.dumps(result,indent=2)+'\n')
