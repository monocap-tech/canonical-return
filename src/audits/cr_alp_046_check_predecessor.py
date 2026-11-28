#!/usr/bin/env python3
"""Exact interval regression against pinned E14 lower bank (not independent entry reconstruction)."""
import json,hashlib
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).parent
f=p/'CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv';g=p/'CR-ALP-040_initial_rational_intervals_1e25.txt'
assert hashlib.sha256(f.read_bytes()).hexdigest()=='7434859d33c9b950b969474c0702c967d88df756b2291ab0b63290f1f5745b8c'
assert hashlib.sha256(g.read_bytes()).hexdigest()=='6f8517df1a00d46eac0b26496f220294b7b82f3571fe600bf6e3ae520092b7a9'
old={};oldm={}
for line in f.read_text().splitlines():
 i,j,l,h=map(int,line.split());old[i,j]=F(l,10**25),F(h,10**25)
for line in g.read_text().splitlines():
 key,l,h=line.split()
 if key.startswith('M_'):
  _,i,j=key.split('_');oldm[int(i),int(j)]=F(int(l),10**25),F(int(h),10**25)
d=json.loads((p/'CR-ALP-046_full_E85_intervals.json').read_text());S=d['scale'];n=0
for b in d['blocks']:
 ix={j:i for i,j in enumerate(b['ids'])};ids=[j for j in b['ids'] if j<=14];delta={}
 for i in ids:
  for j in ids:
   l,h=(F(v,S) for v in b['H'][ix[i]][ix[j]]);ol,oh=old[i,j];delta[i,j]=(l-oh,h-ol)
   l,h=(F(v,S) for v in b['M'][ix[i]][ix[j]]);ol,oh=oldm[i,j]
   assert ol<=h and l<=oh
 for i in ids:
  assert delta[i,i][0]>0
  for j in ids:
   # True omitted PSD Gamma has Cauchy-Schwarz; certify compatible bounds (upper diag budget).
   # A positive result here is a regression, not a full PSD test of the difference.
   l,h=delta[i,j];assert max(abs(l),abs(h))**2<=delta[i,i][1]*delta[j,j][1] or i==j
   n+=1
print('PASS',n,'physical Gram overlap and Gamma remainder entrywise compatibility checks against pinned E14')
