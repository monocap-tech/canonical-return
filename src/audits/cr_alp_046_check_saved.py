#!/usr/bin/env python3
"""Independent exact integer saved-bank verifier; does not rebuild analytic entries."""
import json, hashlib
from pathlib import Path
from fractions import Fraction as F
R=Path(__file__).parent
raw=(R/'CR-ALP-046_full_E85_intervals.json').read_bytes(); data=json.loads(raw); S=data['scale']
assert S==10**340
# Outward compression before independent integer factorization.
D=10**270
for b in data['blocks']:
 for key in ('H','M'):
  b[key]=[[[l//D,-((-h)//D)] for l,h in row] for row in b[key]]
S//=D

def ceil(n,d):return -((-n)//d)
def check(block,exp=None,mass=False):
 H,M=block['H'],block['M']; n=len(H)
 C=[[0]*n for _ in range(n)]; rad=[0]*n
 for i in range(n):
  for j in range(n):
   l,h=(M if mass else H)[i][j]; assert l<=h
   if exp is not None:
    ml,mh=M[i][j];l-=ceil(mh,10**exp);h-=ml//10**exp
   assert (M if mass else H)[i][j]==(M if mass else H)[j][i]
   c=(l+h)//2;C[i][j]=c;rad[i]+=max(c-l,h-c)
 assert max(rad)*10**65<S, 'insufficient interval accuracy'
 for i in range(n):C[i][i]-=max(rad)
 # Separate algorithm: fraction-free Schur elimination using all entries,
 # retaining the lower triangle and verifying exact division independently.
 last=1;pivots=[]
 for k in range(n):
  d=C[k][k];pivots.append(d)
  if d<=0:return False,k+1,max(rad)
  for i in range(k+1,n):
   for j in range(k+1,i+1):
    num=C[i][j]*d-C[i][k]*C[j][k]
    assert num%last==0;C[i][j]=num//last
  last=d
 return True,n,max(rad)
summary=[]
for b in data['blocks']:
 k=b['parity'];assert b['ids']==list(range(k+2,86,2))+[k]
 ok,n,r=check(b,mass=True);assert ok;print('physical metric parity',k,'PASS',n,flush=True)
 outcomes={}
 for e in (20,30,40,50,60,70,80,90,100,110):
  ok,n,r=check(b,e);outcomes[str(e)]=ok
  print('parity',k,'shift 1e-'+str(e),'PASS' if ok else 'INCONCLUSIVE','at pivot',n,flush=True)
  if ok:break
 assert ok
 summary.append({'parity':k,'dimension':len(b['ids']),'first_tested_passing_shift_exponent':e,'trials':outcomes,'interval_row_radius_lt':'1e-65'})
print('SHA256',hashlib.sha256(raw).hexdigest())
(R/'CR-ALP-046_independent_result.json').write_text(json.dumps({'sha256':hashlib.sha256(raw).hexdigest(),'blocks':summary},indent=2)+'\n')
