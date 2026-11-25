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
 raw=(ROOT/f'result_{k}.json').read_bytes(); b=json.loads(raw);S=int(b['scale']);n=44-k;assert b['parity']==k and b['N']==850 and b['bank_sha256']==BANK[k]
 bank=(ROOT.parent/'cr_alp_051'/f'bank_{k}.json').read_bytes();assert hashlib.sha256(bank).hexdigest()==BANK[k]
 for key in ('Kprime','Ghigh','Gwith_bridge','comparison'):
  assert len(b[key])==n and all(len(row)==n for row in b[key]); assert all(b[key][i][j]==b[key][j][i] for i in range(n) for j in range(i))
 shift=F(1,4 if k==1 else 8);assert ld(b['Kprime'],shift,S)==n
 for i in range(n):
  for j in range(n):
   kl,kh=map(int,b['Kprime'][i][j]);gl,gh=map(int,b['Gwith_bridge'][i][j]);l,h=map(int,b['comparison'][i][j]);assert F(l)<=kl-F(200,51)*gh and F(h)>=kh-F(200,51)*gl
 result={'finite_corrected_energy_margin':str(shift),'positive_pivots':n,'result_sha256':hashlib.sha256(raw).hexdigest()}
 for wkey,mkey in [('norm_witness','Gwith_bridge'),('comparison_witness','comparison')]:
  w=b[wkey];z=list(map(int,w['z']));assert len(z)==n and F(sum(x*x for x in z),10**24)==F(w['norm2']);lo,hi=quad(b[mkey],z,S);assert F(w['lo'])<=lo<=hi<=F(w['hi']);result[wkey]={'lower':str(lo),'upper':str(hi)}
 z=list(map(int,b['norm_witness']['z']));hl,hh=quad(b['Ghigh'],z,S);tl,th=quad(b['Gwith_bridge'],z,S);bl,bh=tl-hh,th-hl;fraction=bl/th;assert fraction>F(84 if k==1 else 71,100);result['high_only_witness']={'lower':str(hl),'upper':str(hh)};result['bridge_witness']={'lower':str(bl),'upper':str(bh),'fraction_lower':str(fraction)}
 target=F(1,4) if k==1 else F(7,40);result['target_squared']=str(target**2);result['target_falsified']=F(result['norm_witness']['lower'])>target**2;result['scalar_correction_falsified']=F(result['comparison_witness']['upper'])<0
 assert result['target_falsified'] and result['scalar_correction_falsified'];print('bridge fraction lower',float(fraction));print('parity',k,'positive pivots',n,'outside Gram witness lower',float(F(result['norm_witness']['lower'])),'target squared',float(target**2),'comparison witness upper',float(F(result['comparison_witness']['upper'])),'target falsified',result['target_falsified'],'scalar correction falsified',result['scalar_correction_falsified'])
 out[str(k)]=result
(ROOT/'verification_result.json').write_text(json.dumps(out,indent=2)+'\n')
