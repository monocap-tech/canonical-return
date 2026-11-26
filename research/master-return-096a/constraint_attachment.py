"""RETURN096-A: exact retained/complement index and constraint-chart binding.
The analytic Fourier and divided-difference identities are in ATTACHMENT.md.
This checks all finite data used by those identities, not the original
Weil comparator or its domain continuity.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p,h):
 f=root/p;assert sha(f)==h,p;return json.loads(f.read_bytes())
w=load('return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
u=load('return057/RETURN057_weights.json','090ec051ef390b339357496b805808a846268352c43b70befb1c3ede3614a926')
raw=load('return059/RETURN059_raw_receipt.json','2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8')
b=load('return062/RETURN062_bank_receipt.json','429cb13fe51f9f92e684fffb9db6d6eb31004c06e20a81dc1877c94065f2cef2')
cr=load('return063/RETURN063_constraints_receipt.json','dd219b5a16d31acec79d142589e079d9c152cf69178fd4780d1b2550760e568e')
iv=load('return063/RETURN063_inverse_receipt.json','d1d64f3b8197fb8a91d3ad2ba1cf92f9b1536737cb7affc98a7366cb1d16bc5a')
pr=load('return077/proposal.json','803aa2124522c5565864c35623d92fce2eb1f1fab616347a6020387f65068120')
E=set(w['even_exceptional']);O=set(w['unpaid_odd']);P=set(w['paid_odd'])
assert all(S=={-j for j in S} for S in (E,O,P))
assert 0 in E and all(j%2==0 for j in E) and all(j%2==1 for j in O|P)
assert O.isdisjoint(P)
J=sorted(j for j in O if j>0);assert len(J)==5
even=sorted(j for j in E if j>0);assert len(even)==37
ks=raw['ks'];assert len(ks)==44 and ks==sorted(set(ks))
weights=dict(w['positive_even_weights']);assert all(k%2==0 and k>0 and v>0 for k,v in weights.items())
assert ks==sorted(weights)[:44] and ks[-1]==116
assert E.isdisjoint(weights)
all_even=set(range(0,32768,2))
assert (E&all_even)|set(weights)==all_even
ret=set(even)|set(ks);good=sorted(set(weights)-set(ks))
assert len(ret)==81 and len(good)==16302
assert set(good)==set(range(118,32768,2))-E
assert set(good).isdisjoint(ret)
assert all(j<32768 for j in ret)
pos=dict(u['odd_good']);neg={j:v for j,v in u['odd_deficit'] if v}
assert len(pos)==481 and len(neg)==32
assert set(pos).isdisjoint(neg) and not (set(pos)|set(neg))&set(J)
assert set(neg)<=P and all(j%2 and v>0 for j,v in pos.items()) and all(j%2 and v>0 for j,v in neg.items())
# Validate the exact rational chart ranks with elementary elimination.
def rank(A):
 a=[[F(x) for x in row] for row in A];r=0
 for j in range(len(a[0])):
  piv=next((i for i in range(r,len(a)) if a[i][j]),None)
  if piv is None:continue
  a[r],a[piv]=a[piv],a[r];v=a[r][j];a[r]=[x/v for x in a[r]]
  for i in range(r+1,len(a)):
   v=a[i][j]
   if v:a[i]=[x-v*y for x,y in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r
# Independent Machin bound for c=25/(16*pi^2), the sinh node.
def atan(q,n):
 a=sum((F((-1)**i,(2*i+1)*q**(2*i+1)) for i in range(n)),F(0));z=a+F((-1)**n,(2*n+1)*q**(2*n+1))
 return min(a,z),max(a,z)
a,z=atan(5,160);c,d=atan(239,64);pl,pu=16*a-4*d,16*z-4*c
cl,cu=map(F,b['c']);assert 0<cl<=F(25,16)/pu**2<=F(25,16)/pl**2<=cu
# Rational divided-difference identity at sample rational nodes:
# sum_j 1/[(x-t_j) prod_{v!=j}(t_j-t_v)] = 1/prod_j(x-t_j).
# The general polynomial proof is in the note, this exact control
# checks signs and factors on the actual integer nodes.
for par in (0,1):
 cs=J+(['sinh'] if par else [])
 blk=cr['blocks'][par];ib=iv['blocks'][par];assert blk['nodes']==cs and ib['nodes'][:len(cs)]==cs
 assert set(ib['nodes'][len(cs):])==set(pos)|set(neg)
 assert len(ib['nodes'])==518+par
 assert len(blk['Z'])==rank(blk['Z'])==len(cs) and blk['Zscale']>0
 assert all(len(row)==len(cs) for row in blk['Z'])
 # Check actual odd nodes and one exact point inside the sinh enclosure.
 ts=[F(j*j) if j!='sinh' else -(cl+cu)/2 for j in cs]
 assert len(set(ts))==len(ts)
 for k in (118,120,16384,32768,262144):
  x=F(k*k)
  for n in range(1,len(ts)+1):
   t=ts[:n]
   lhs=sum((1/((x-u)*math.prod(u-v for v in t if v!=u)) for u in t),F(0))
   assert lhs==1/math.prod(x-v for v in t)
 print('constraint chart and partial fractions',par,len(cs),'rank',rank(blk['Z']),flush=True)
# The retained basis has 81 positive even modes, five odd modes, plus
# constant (even parity) or sinh (odd parity). The complement consists
# precisely of good and infinite even-grid coordinates with these constraints.
assert 81+len(J)+1==87
assert len(pr['js'])==len(set(pr['js']))==500
assert all(j%2 and j>0 for j in pr['js'])
assert not set(pr['js'])&set(J)
records={'status':'PASS','scope':'literal index partition and rational constraint-chart identity; analytic Fourier proof in note; no original Weil comparator or density proof','retained_each':87,'positive_retained_even_modes':81,'unpaid_odd_constraints':J,'sinh_constraint':True,'finite_complement':len(good),'excluded_late_grid':sorted(j for j in E if j>=4544 and j<=16384),'chart_ranks':[5,6],'common_auxiliary_dimensions':[518,519],'restoration_carriers':500,'infinite_grid_start':32768}
(out/'receipt.json').write_text(json.dumps(records,indent=2)+'\n')
print('RETURN096-A ORIGINAL INDEX PARTITION AND CONSTRAINT CHART PASS',flush=True)
