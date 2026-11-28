#!/usr/bin/env python3
"""CR-ALP-040 exact source-derived two-pole shell Gram audit.

For E_M=Pi0 S_M, Q_M e_{M+1} is the genuine orthogonal shell vector
and has norm²=G_{M+1}/G_M, with G_M=B*R_MB high sine pole Gram.
The two parity shells at M=12 are mutually orthogonal. All inequalities
below use rational intervals directly from the source-attached input bank.
"""
from pathlib import Path
from fractions import Fraction as F
import json
ROOT=Path(__file__).resolve().parent; S=10**25
P={}
for ln in (ROOT/'CR-ALP-040_initial_rational_intervals_1e25.txt').read_text().splitlines():
 k,l,h=ln.split();P[k]=(F(int(l),S),F(int(h),S))

def add(a,b): return a[0]+b[0], a[1]+b[1]
def neg(a): return -a[1], -a[0]
def sub(a,b): return add(a,neg(b))
def mul(a,b):
 z=[x*y for x in a for y in b]
 return min(z),max(z)
def div(a,b):
 assert b[0]>0
 return mul(a,(F(1,b[1]),F(1,b[0])))
def const(a): return F(a),F(a)

a=P['A'];b_even=div(add(P['EHALF'],P['EHI']),const(2));b_odd=div(sub(P['EHALF'],P['EHI']),const(2))
Gram={1:add(a,P['SinhA']),2:sub(P['SinhA'],a)}
for parity in (1,2):
 abs_b=b_even if parity==1 else b_odd
 for j in range(parity,15,2):
  nu=P[f'nu_{j}'];sqrtnu=mul(nu,nu)
  c=div(mul(mul(const(2),abs_b),nu),mul(P['SQRT_A'],add(sqrtnu,const(F(1,4)))))
  Gram[parity]=sub(Gram[parity],mul(c,c))
  if j==parity+10: P[f'G12_{parity}']=Gram[parity]
  if j==parity+12: P[f'G14_{parity}']=Gram[parity]

for k,lo12,hi12,lo14,hi14 in [
 (1,F(98,1000),F(100,1000),F(83,1000),F(86,1000)),
 (2,F(22,1000),F(24,1000),F(19,1000),F(21,1000)),
]:
 g12=P[f'G12_{k}'];g14=P[f'G14_{k}']
 assert lo12<g12[0]<=g12[1]<hi12,(k,g12)
 assert lo14<g14[0]<=g14[1]<hi14,(k,g14)
 ratio=div(g14,g12)
 assert F(4,5)<ratio[0]<=ratio[1]<F(9,10)
 print('PASS pole parity',k,'G12',*[float(x) for x in g12],
       'G14',*[float(x) for x in g14], 'shell norm²',*[float(x) for x in ratio])

print('CR-ALP-040 EXACT RATIONAL SHELL-GRAM PASS; canonical high-tail projector and shell identity require continuum Parseval premises')
