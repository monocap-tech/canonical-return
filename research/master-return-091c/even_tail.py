"""RETURN091-C: independent even modeled asymptotic tail at K=262144.
Reuse the exact-rational interval helpers and phase-distance bounds of 090B,
not its odd gauge, logarithmic drift, or parent leading contraction.
"""
import sys,json,math,hashlib
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from mpmath import iv
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
K=262144; A=F(87,500); S=10**40; Q=10**55; N=6
def load(p,sha=None):
 b=(root/p).read_bytes()
 if sha: assert hashlib.sha256(b).hexdigest()==sha,p
 return json.loads(b)
def P(x):return F(x),F(x)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def isum(a):
 l=h=F(0)
 for x,y in a:l+=x;h+=y
 return l,h
def mag(a):return max(map(abs,a))
def wide(a,e):assert e>=0;return a[0]-e,a[1]+e
def enc(a,s=S):return [math.floor(a[0]*s),math.ceil(a[1]*s)]
def quant(a):return tuple(F(x,S) for x in enc(a))
def inv(a):assert a[0]>0;return 1/a[1],1/a[0]
def ivf(x):x=F(x);return iv.mpf(x.numerator)/x.denominator
def interval(a):return iv.mpf([ivf(a[0]).a,ivf(a[1]).b])
# Exact extraction of binary endpoints of the interval library, then rational rounding.
def binary(v):
 sign,man,exp,bc=v
 return F((-1 if sign else 1)*man)*F(2)**exp
def pair(v):
 a,b=v._mpi_;return quant((binary(a),binary(b)))
iv.dps=70
pi=pair(iv.pi);c=mul(P(F(2,5)),pi);ci=inv(c);delta=8*ci[1]/K
weights={n:pair(iv.ln(p)/iv.sqrt(n)) for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]}
C={}
for n,w in weights.items():C[F(n)]=w;C[F(1,n)]=w
assert sum(mag(w) for w in C.values())<9
C2={}
for x,w in C.items():
 for y,v in C.items():C2[x*y]=add(C2.get(x*y,P(0)),mul(w,v))
@lru_cache(None)
def gap(z):
 a=4*iv.pi*iv.ln(ivf(z))/5
 d=pair(2*abs(iv.sin(a/2)))
 assert d[0]>0
 return d[0]
@lru_cache(None)
def osc_cap(which,p,q):
 poly=C if which==1 else C2
 # Abel: sum over m>K/2, total variation plus first value <=2 sup.
 sup=F(1) if q==0 else F(q,p)**q
 return 2*sup*sum(mag(w)/gap(z) for z,w in poly.items() if z!=1)
@lru_cache(None)
def I(p,q):
 assert p>=2 and q in (0,1,2)
 # x=k/K; density K/2. Integral plus bounded total variation.
 m=F(K*math.factorial(q),2*(p-1)**(q+1))
 if q==0:return m-1,m
 e=2*F(q,p)**q
 return max(F(0),m-e),m+e
@lru_cache(None)
def J(p):
 # Weighted sum of C-log(x)-drift, x^-p.
 return quant(wide(neg(I(p,1)),osc_cap(1,p,0)+delta*I(p,0)[1]))

# Even source is (pi/2 + 2 sum w sin(4*pi*m*log(n)/5))/(c*k).
# Sine Fourier coefficients are imaginary with opposite signs at reciprocal
# frequencies. Their products are real, including the positive constant.
Sin={z:mul(P(-1 if z>1 else 1),w) for z,w in C.items()}
Sin2={}
for x,w in Sin.items():
 for y,v in Sin.items():Sin2[x*y]=add(Sin2.get(x*y,P(0)),neg(mul(w,v)))
def phase_cap(poly):
 # For decreasing x^-p the Abel bound is 2/gap times the first weight.
 return 2*sum(mag(w)/gap(z) for z,w in poly.items() if z!=1)
cap1=phase_cap(Sin);cap2=phase_cap(Sin2)
mean=add(mul(mul(pi,pi),P(F(1,4))),Sin2[F(1)])
D=quant(mul(P(1/(A*K*K)),mul(mul(ci,ci),wide(mul(mean,I(2,0)),cap2+mag(pi)*cap1))))
assert D[0]>0
raw=load('return059/RETURN059_raw_receipt.json')
lead=load('return060/RETURN060_leading_receipt.json','e9001e8dfad3e1a0fdedbaba9d27c258afb241134d16f98bbff7e6477007fe57')
cr=load('return063/RETURN063_constraints_receipt.json'); ib=load('return063/RETURN063_inverse_receipt.json')['blocks'][0]
cb=cr['blocks'][0];nodes=cb['nodes']
c0,c1=map(F,load('return062/RETURN062_bank_receipt.json')['c'])
def node(j):return (-c1,-c0) if j=='sinh' else P(j*j)
cols=list(zip(*raw['blocks'][0]['C']));B=lead['blocks'][0]
b=[tuple(F(x,lead['scale']*B['chart_denominator']) for x in r) for r in B['boundary']]
assert len(b)==len(cols)==87
for x,col in zip(b,cols):
 assert x==isum(mul(tuple(F(y,Q) for y in a),P(F(z,raw['T']))) for a,z in zip(B['raw_boundary'],col))
def polynomial(ts,offset,factor):
 a=[P(1)]+[P(0)]*(N-1);m=[F(1)]+[F(0)]*(N-1);total=F(1)
 for j in ts:
  z=mul(node(j),P(F(1,K*K)));r=mag(z);assert r<1
  w=[P(1)]
  for i in range(1,N):w.append(mul(w[-1],z))
  a=[isum(mul(a[i],w[t-i]) for i in range(t+1)) for t in range(N)]
  m=[sum(m[i]*r**(t-i) for i in range(t+1)) for t in range(N)]
  total/=1-r
 f=F(factor);assert f>=0
 return [P(0) if i<offset else mul(P(f),a[i-offset]) for i in range(N)],f*(total-sum(m[:max(0,N-offset)]))
def combine(row,polys):
 return [quant(isum(mul(P(z),a[n]) for z,(a,e) in zip(row,polys))) for n in range(N)],sum(abs(z)*e for z,(a,e) in zip(row,polys))+F(N,S)

# Even auxiliary rows are k^-2 times polynomial in (K/k)^2.
base=[polynomial(nodes[:i+1],i,F(116)*F(116,K)**(2*i)) for i in range(len(nodes))]
Z=[[F(z,cb['Zscale']) for z in row] for row in cb['Z']]
aux=[combine(row,base) for row in Z]
aux += [polynomial([j],0,F(ib['scales'][str(j)])*j) for j in ib['nodes'][len(nodes):]]
assert len(aux)==518
def J_even(p):return wide(mul(mul(pi,P(F(1,2))),I(p,0)),cap1)
def cross(row):
 a,e=row
 z=isum(mul(co,J_even(3+2*n)) for n,co in enumerate(a))
 err=e*(mag(pi)/2+9)*I(3,0)[1]
 return quant(mul(P(1/(A*K**3)),mul(ci,wide(z,err))))
TT=[[quant(mul(D,mul(x,y))) for y in b] for x in b]
gamma=[cross(row) for row in aux]
U=[[quant(mul(x,y)) for y in b] for x in gamma]
fd=load('return071/RETURN071_contract0.json','af515fe21b0c2186d3fb6c6e57764e821431fb5181be2522cc75263f90125e8f')['variant']
hd=load('return071/RETURN071_hybrid.json','816388b20efc5d06bd32b85d09dca0d51ae8294b7f1501761816f9f14f661746')
O=hd['scale'];hd=hd['variants'][0]
receipt={}
for key,mat in [('tt',TT),('U',U)]:
 contained=overlap=total=0
 for i,row in enumerate(mat):
  for j,(l,u) in enumerate(row):
   a,b0=hd[key][i][j];c0_,d=fd[key][i][j]
   lo,hi=F(a-d,O),F(b0-c0_,O)
   contained+=lo<=l<=u<=hi;overlap+=lo<=u and l<=hi;total+=1
 receipt[key]={'entries':total,'contained':contained,'overlap':overlap,'disjoint':total-overlap}
 print('even asymptotic',key,receipt[key],flush=True)
(out/'even_tail.json').write_text(json.dumps({'K':K,'scale':S,'matrices':{k:[[enc(x) for x in row] for row in m] for k,m in [('tt',TT),('U',U)]},'checks':receipt},separators=(',',':')))
assert all(x['contained']==x['entries'] for x in receipt.values()),receipt
print('RETURN091-C EVEN ASYMPTOTIC CONTAINMENT PASS',flush=True)
