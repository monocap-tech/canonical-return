#!/usr/bin/env python3
"""CR-ALP-033: independent, EXACT Fraction check of source-attached 12-mode lower bank.

The bank is the 30000 individual + 4219 decreasing-Loewner endpoint groups,
covering [30000, 2000000). This is a lower bound on the complete original
signed Weil form on the genuine pole-null physical section; an indefinite
LOWER matrix does NOT certify a negative direction for the complete form.
Continuum source identities are analytic inputs, not kernel-checked.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
S=10**25

class I:
 def __init__(self,lo,hi=None):
  if isinstance(lo,I): self.lo,self.hi=lo.lo,lo.hi
  else: self.lo,self.hi=F(lo),F(lo if hi is None else hi)
  assert self.lo<=self.hi
 def __neg__(self):return I(-self.hi,-self.lo)
 def __add__(self,v):
  b=I(v);return I(self.lo+b.lo,self.hi+b.hi)
 __radd__=__add__
 def __sub__(self,v):return self+-I(v)
 def __mul__(self,v):
  b=I(v);items=[x*y for x in [self.lo,self.hi] for y in [b.lo,b.hi]]
  return I(min(items),max(items))
 __rmul__=__mul__
 def __truediv__(self,v):
  b=I(v);assert b.lo>0
  return self*I(1/b.hi,1/b.lo)

def fixed(lo,hi):return I(F(int(lo),S),F(int(hi),S))
initial={}
file_i=ROOT/'CR-ALP-033_initial_rational_intervals_1e25.txt'
file_h=ROOT/'CR-ALP-033_refined_H_exact_matrix_intervals.tsv'
for line in file_i.read_text().splitlines():
 key,l,u=line.split();initial[key]=fixed(l,u)
assert len(initial)==177
raw={}
for line in file_h.read_text().splitlines():
 i,j,l,u=line.split();raw[int(i),int(j)]=fixed(l,u)
assert len(raw)==72

def symmat(A,i,j):
 a=A[(i,j)];b=A[(j,i)];assert max(a.lo,b.lo)<=min(a.hi,b.hi)
 return I(min(a.lo,b.lo),max(a.hi,b.hi))

def H(i,j):return symmat(raw,i,j)
def M(i,j):return symmat({(a,b):initial[f'M_{a}_{b}'] for (a,b) in raw},i,j)

def ldl(order,shift,expect_positive):
 L={};p=[]
 for t,i in enumerate(order):
  piv=H(i,i)-shift*M(i,i)
  for k in range(t):piv=piv-L[t,k]*L[t,k]*p[k]
  if t<expect_positive:assert piv.lo>0,(order,t,piv.lo,piv.hi)
  p.append(piv)
  if piv.lo<=0:
   print('UNRESOLVED/nonpositive LOWER interval pivot',order,t,float(piv.lo),float(piv.hi))
   break
  for u in range(t+1,len(order)):
   j=order[u];v=H(j,i)-shift*M(j,i)
   for k in range(t):v=v-L[u,k]*L[t,k]*p[k]
   L[u,t]=v/piv
 return p

odd=(4,6,8,10,12,2);even=(3,5,7,9,11,1)
assert all(v.lo>0 for v in ldl((3,5,7,9,1),F(1,10**11),5))
assert all(v.lo>0 for v in ldl((4,6,8,10,2),F(1,10**11),5))
# Validate full physical metric of odd 6x6 by exact pivots with H=M (not assumed identity).
def mass_ldl(order):
 L={};p=[]
 for t,i in enumerate(order):
  v=M(i,i)
  for k in range(t):v=v-L[t,k]*L[t,k]*p[k]
  assert v.lo>0
  p.append(v)
  for u in range(t+1,len(order)):
   tval=M(order[u],i)
   for k in range(t):tval=tval-L[u,k]*L[t,k]*p[k]
   L[u,t]=tval/v
 return p
assert len(mass_ldl(odd))==6
shift=F(1,10**11)
po=ldl(odd,shift,6)
assert len(po)==6 and all(v.lo>0 for v in po)
pe=ldl(even,F(0),4)
assert pe[4].hi<0,('even 5th exact lower pivot sign changed',pe[4].lo,pe[4].hi)
print('PASS: complete ORIGINAL signed odd 6-mode lower matrix > 1e-11 physical M')
print('PASS: previous even & odd five-mode 1e-11 lower certificates retained')
print('NOTE: even six-mode grouped lower matrix has a certified negative fifth pivot; THIS IS NOT NEGATIVE WEIL TEST')
report={'stage':'CR-ALP-033','status':'PARTIAL EXACT 11-DIMENSIONAL CERTIFICATE; EVEN 12TH DIRECTION UNRESOLVED','individual_gamma_terms':30000,'grouped_tail_groups':4219,'covered_gamma_indices':'[0,2000000)','true_original_prime_signs':'negative paired seven source von Mangoldt terms','pole_null_11_dimensional_certificate':{'basis_even':'w3,w5,w7,w9,u1','basis_odd':'w4,w6,w8,w10,w12,u2','margin':'1e-11 strict physical generalized lower bound'},'even_12_mode_lower_matrix_fifth_pivot_bounds_rounded_outward_1e15':[str(pe[4].lo*10**15//1)+'/10^15',str(-((-pe[4].hi*10**15)//1))+'/10^15'],'odd_six_shifted_pivot_lower_rounded_down_1e15':[str(v.lo*10**15//1)+'/10^15' for v in po], 'source_sha256':hashlib.sha256(file_i.read_bytes()).hexdigest(),'matrix_sha256':hashlib.sha256(file_h.read_bytes()).hexdigest(),'refined_Gamma_two_sided_upper':'W_complete - H_lower <= (1/10^12)D + (1/1700000)T on finite H1 section','not_claimed':'Negative complete form, full E12 positivity, full E15000 Schur, Lean or RH'}
(ROOT/'CR-ALP-033_partial_11_mode_result.json').write_text(json.dumps(report,indent=2)+'\n')

# New 4219-group two-sided Loewner envelope, same analytic lemma as CR-ALP-032.
start,stop=30000,2000000
groups=[]
while start<stop:
    length=max(1,start//1000)
    end=min(stop,start+length)
    groups.append((start,end))
    start=end
assert len(groups)==4219 and groups[0][0]==30000 and groups[-1][1]==2000000
assert sum(U-L for L,U in groups)==1970000
cD=sum((F(24*(U-L)*(U-L-1),1)/F(4*L+1,2)**4 for L,U in groups),F(0))
cT=sum((F(16*(U-L)*(U-L-1),1)/F(4*L+1,2)**3 for L,U in groups),F(0))
beta=F(4*stop+1,2)
cD+=4/beta**3+1/beta**2
cT+=4/beta**2+2/beta
assert cD<F(1,10**12)
assert cT<F(1,1700000)
print('PASS: two-sided refined Gamma remainder <= D/10^12 + T/1,700,000')
print('CR-ALP-033 INDEPENDENT EXACT FRACTION AUDIT: PARTIAL 11-DIMENSIONAL POSITIVITY PASS')
