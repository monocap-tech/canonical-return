"""Exact rational mesh for a 96-node physical overlap rule. All shift approximations have error <10^-70."""
from fractions import Fraction as F
from pathlib import Path
from math import floor,ceil,isqrt
import json
root=Path(__file__).resolve().parent; D=10**80; d=F(11,5)
def dec(n):return ('-' if n<0 else '')+str(abs(n)//D)+'.'+str(abs(n)%D).zfill(80)
def enc(x,y=None):return dec(floor(x*D))+' '+dec(ceil((x if y is None else y)*D))
def lg(n):
 z=F(n-1,n+1);v=z;s=F(0)
 for j in range(400):s+=2*v/(2*j+1);v*=z*z
 return s,s+2*v/(801*(1-z*z))
ns=[2,3,4,5,7,8,9];ps=[2,3,2,5,7,2,3];sh=[];weights=[]
for n,p in zip(ns,ps):
 l,h=lg(n);x=F(round((l+h)*10**74/2),10**74);assert max(abs(l-x),abs(h-x))<F(1,10**70);sh.append(x)
 l,h=lg(p);s=F(isqrt(n*10**160),10**80);weights.append((l/(s+F(1,10**80)),h/s))
knots=sorted({F(0),d,*sh,*(d-x for x in sh)});cmax=F(3*2**31)-F(1,2)
from bisect import bisect_right
lines=[str(len(knots)-1)]
for l,h in zip(knots,knots[1:]):
 m=(l+h)/2;active=[(j,s) for j,x in enumerate(sh) for s in (-1,1) if 0<m+s*x<d];lines.append(str(len(active))+' '+' '.join(f'{j} {s}' for j,s in active))
for x,(l,h) in zip(sh,weights):lines.append(enc(x)+' '+enc(l,h))
counts=[]
for t in [F(0),*sh]:
 end=(d-t)/2;cuts=sorted({F(0),end,*[x for k in knots for x in (k,k-t) if 0<x<end]});out=[]
 def panel(l,r):
  m=(l+r)/2;h=(r-l)/2;R=4*h
  if h<=F(1,128) and all(min(y,d-y)>=2*R or cmax*R<=F(1,2) for y in (m,m+t)):
   i=bisect_right(knots,m)-1;j=bisect_right(knots,m+t)-1
   assert knots[i]<=l and r<=knots[i+1] and knots[j]<=l+t and r+t<=knots[j+1]
   out.append((m,h,i,j));return
  panel(l,m);panel(m,r)
 for l,h in zip(cuts,cuts[1:]):panel(l,h)
 assert sum(2*h for _,h,_,_ in out)==end
 lines.append(str(len(out)))
 for m,h,i,j in out:lines.append(enc(m)+' '+enc(h)+f' {i} {j}')
 counts.append(len(out))
(root/'mesh.txt').write_text('\n'.join(lines)+'\n');print('exact half-domain panel counts',counts)
# Cauchy Taylor remainder at R=4h, m=96. Uniform complex bound M=10^40.
e=F(4,3)*10**40/4**96;eg=F(4,3)*100/4**96
assert 2*d*(e*100+eg*10**40+e*eg)<F(1,10**12)
assert 4*d*e< F(1,10**12)
(root/'mesh_result.json').write_text(json.dumps({'panels':counts,'max_shift_error':'1e-70','uniform_complex_bound_assumed_for_analytic_gate':'1e40','raw_taylor_error':str(e),'projection_product_padding':'1e-12','note':'Mesh geometry and ideal Gauss error verified; complex coefficient bound must be checked in evaluator.'},indent=2)+'\n')
