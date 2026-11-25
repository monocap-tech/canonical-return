#!/usr/bin/env python3
"""Independent stdlib interval contraction replay; no transcendental source re-evaluation."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent.parent
p=Path(__file__).with_name('RETURN048_retained_receipt.json');r=json.loads(p.read_bytes())
parent=(root/'return047/RETURN047_payment_receipt.json').read_bytes();assert hashlib.sha256(parent).hexdigest()==r['parent_sha256']=='def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f'
a=json.loads(parent);ks=sorted(k for k in set(a['even_exceptional']+a['unpaid_odd']) if k>0)
S=int(r['scale']);Q=int(r['congruence_scale']);assert S==10**55 and Q==10**40
assert [b['parity'] for b in r['blocks']]==[0,1]
for b in r['blocks']:
 n=len(b['names']);assert n==43 and b['status']=='PASS'
 assert b['names']==(['1']+['cos'+str(k) for k in ks] if b['parity']==0 else ['sin'+str(k) for k in ks]+['sinh(x/2)'])
 C=b['congruence'];A=b['energy'];M=b['mass']
 assert len(C)==len(A)==len(M)==n and all(len(x)==n for x in C+A+M)
 for mat in (A,M):
  for i in range(n):
   for j in range(n):
    l,u=mat[i][j];assert isinstance(l,int) and isinstance(u,int) and l<=u and u-l<=2 and mat[i][j]==mat[j][i]
 # Contract endpoint intervals directly, not using the generator center/radius method.
 # Two-stage interval matrix multiplication gives exactly the same interval as summing
 # independent entries directly (no intersections or numerical proposals in replay).
 def mul(l,u,c):return (l*c,u*c) if c>=0 else (u*c,l*c)
 AC=[]
 for i in range(n):
  row=[]
  for j in range(n):
   terms=[mul(*A[i][k],C[k][j]) for k in range(n)]
   row.append((sum(v[0] for v in terms),sum(v[1] for v in terms)))
  AC.append(row)
 T=[]
 for i in range(n):
  row=[]
  for j in range(n):
   terms=[mul(*AC[k][j],C[k][i]) for k in range(n)]
   row.append((sum(v[0] for v in terms),sum(v[1] for v in terms)))
  T.append(row)
 margins=[T[i][i][0]-sum(max(abs(x) for x in T[i][j]) for j in range(n) if i!=j) for i in range(n)]
 margin=min(margins);assert margin>0 and 2*margin==int(b['margin_numerator']);assert int(b['margin_denominator'])==2*S*Q*Q
 mass=max(sum(max(abs(x) for x in v) for v in row) for row in M)
 assert mass==int(b['mass_upper_numerator']);fro=sum(x*x for row in C for x in row);assert fro==int(b['frobenius_numerator'])
 floor=F(margin,mass*fro);assert floor==F(b['physical_floor'])
 assert floor> (F(250,10**23) if b['parity']==0 else F(837,10**21))
 print('PASS parity',b['parity'],'dimension',n,'transformed margin',float(F(margin,S*Q*Q)),'physical floor',float(floor))
print('receipt SHA256',hashlib.sha256(p.read_bytes()).hexdigest())
