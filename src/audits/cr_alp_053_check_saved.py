#!/usr/bin/env python3
"""Exact replay of saved finite outside-band Gram and comparison witnesses.
Does not reconstruct analytic entry or polynomial arithmetic."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
BANK={1:'7a6b1d08bc966c12d189d7387e8699a3e3fab516b35f887184c90d0085d746a4',2:'0133407e41b7cf6806366f2ea646acf143fde3289e6377fb14e8d46b21a38a2a'}
def quad(A,z,S):
 lo=hi=0
 for i,zi in enumerate(z):
  for j,zj in enumerate(z):
   l,h=map(int,A[i][j]); assert l<=h
   v=zi*zj;lo+=min(v*l,v*h);hi+=max(v*l,v*h)
 norm=sum(x*x for x in z);assert norm>0
 return F(lo,S*norm),F(hi,S*norm)
def ld(A,shift,S):
 C=[];radius=[];ss=shift*S;assert ss.denominator==1
 for i,row in enumerate(A):
  cc=[];rr=0
  for j,(lo,hi) in enumerate(row):
   l,h=int(lo),int(hi);assert l<=h
   if i==j:l-=int(ss);h-=int(ss)
   m=(l+h)//2;cc.append(m);rr+=max(m-l,h-m)
  C.append(cc);radius.append(rr)
 n=len(C)
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
out={}
for k in (1,2):
 raw=(ROOT/f'result_{k}.json').read_bytes(); b=json.loads(raw);S=int(b['scale']);n=44-k
 assert b['parity']==k and b['trial']==600 and b['band_end']==850
 for key in ('Knew','Gtrial','Kfinite','Gout','comparison'):
  assert len(b[key])==n and all(len(row)==n for row in b[key]);assert all(b[key][i][j]==b[key][j][i] for i in range(n) for j in range(i))
 for i in range(n):
  for j in range(n):
   l,h=map(int,b['Knew'][i][j]);gl,gh=map(int,b['Gtrial'][i][j]);kl,kh=map(int,b['Kfinite'][i][j]);ol,oh=map(int,b['Gout'][i][j]);cl,ch=map(int,b['comparison'][i][j]);
   assert kl<=l-F(200,51)*gh and kh>=h-F(200,51)*gl
   assert cl<=kl-F(200,51)*oh and ch>=kh-F(200,51)*ol
 shift=F(b['finite_shift']);assert b['finite_test']['all_positive'] and ld(b['Kfinite'],shift,S)==n
 gram_upper=F(max(sum(max(abs(int(l)),abs(int(h))) for l,h in row) for row in b['Gtrial']),S)
 result={'finite_shift':str(shift),'positive_pivots':n,'trial_gram_upper':str(gram_upper),'result_sha256':hashlib.sha256(raw).hexdigest()}
 for wkey,mkey in [('norm_witness','Gout'),('comparison_witness','comparison')]:
  w=b[wkey];z=[]
  for row in w['z']:
   assert len(row)==1;l,h=map(int,row[0]);assert l==h;z.append(l)
  lo,hi=quad(b[mkey],z,S);assert F(w['lo'])<=lo<=hi<=F(w['hi']);result[wkey]={'lower':str(lo),'upper':str(hi)}
 # Smaller coefficient is only a sensitivity test on the checked finite residual.
 def reduced(penalty):
  a=[]
  for row,grow in zip(b['Kfinite'],b['Gout']):
   rr=[]
   for (l,h),(gl,gh) in zip(row,grow):
    lo=F(int(l))-penalty*int(gh);hi=F(int(h))-penalty*int(gl);rr.append([lo.numerator//lo.denominator,-((-hi.numerator)//hi.denominator)])
   a.append(rr)
  return a
 assert ld(reduced(F(5,2)),F(1,1000),S)==n
 assert ld(reduced(F(2)),F(1,40),S)==n
 result['finite_band_sensitivity']={'penalty_5/2_shift':'1/1000','penalty_2_shift':'1/40'}
 result['scalar_correction_falsified']=F(result['comparison_witness']['upper'])<0
 print('parity',k,'positive pivots',n,'shift',shift,'trial Gram upper',float(gram_upper),'outside Gram witness lower',float(F(result['norm_witness']['lower'])),'comparison witness upper',float(F(result['comparison_witness']['upper'])),'scalar correction falsified',result['scalar_correction_falsified'])
 out[str(k)]=result
(ROOT/'verification_result.json').write_text(json.dumps(out,indent=2)+'\n')
