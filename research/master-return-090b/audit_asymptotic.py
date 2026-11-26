"""RETURN090-B: independent infinite modeled-source contractions.
Shared source inputs: raw chart, normalized boundary and g(R) primitive;
shared analytic drift bound |g(ck)-g(cK)-log(k/K)| <= 8/(cK).
No parent leading/gauge contraction or terminal solver is called.
"""
import sys,json,math,hashlib
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from mpmath import iv
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
K=262144; A=F(87,500); A1=F(1199,1000); lam=A/A1; S=10**40; Q=10**55; N=6
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
D=wide(add(mul(C2[F(1)],I(2,0)),I(2,2)),osc_cap(2,2,0)+2*osc_cap(1,2,1)+2*delta*(9*I(2,0)[1]+I(2,1)[1])+delta**2*I(2,0)[1])
D=quant(mul(P(1/(A*K*K)),mul(mul(ci,ci),D)));assert D[0]>0
mom=load('return066/RETURN066_moments.json')
raw=load('return059/RETURN059_raw_receipt.json')
lead=load('return060/RETURN060_leading_receipt.json','e9001e8dfad3e1a0fdedbaba9d27c258afb241134d16f98bbff7e6477007fe57')
cr=load('return063/RETURN063_constraints_receipt.json'); ib=load('return063/RETURN063_inverse_receipt.json')['blocks'][1]
cb=cr['blocks'][1];nodes=cb['nodes'];assert nodes[-1]=='sinh'
c0,c1=map(F,load('return062/RETURN062_bank_receipt.json')['c'])
def node(j):return (-c1,-c0) if j=='sinh' else P(j*j)
pr=load('return077/proposal.json');js=pr['js'];BS=pr['scale'];rows=pr['rows'];assert len(rows)==8
cols=list(zip(*raw['blocks'][1]['C']));B=lead['blocks'][1]
b=[tuple(F(x,lead['scale']*B['chart_denominator']) for x in r) for r in B['boundary']]
assert len(b)==len(cols)==87
for x,col in zip(b,cols):
 assert x==isum(mul(tuple(F(y,Q) for y in a),P(F(z,raw['T']))) for a,z in zip(B['raw_boundary'],col))
# Shared g(R) source primitive only; the chart contraction is independent.
sys.path[:0]=[str(root/'return069'),str(root/'return063')]
from return069_grouped import fast,src
R=2*iv.pi*K/5
gR=src.G0+fast.sums((0,K))[0]+4*R*R/(R*R+iv.mpf(1)/4)
graw=[]
for row in mom['blocks'][1]['raw_coeff']:
 l,bb,der,cc,dd=[tuple(F(x,Q) for x in a) for a in row]
 graw.append(pair((interval(bb)*(gR+src.G0+4)-2*interval(cc)-interval(dd))*5/(2*iv.pi)))
v=[quant(isum(mul(a,P(F(z,raw['T']))) for a,z in zip(graw,col))) for col in cols]
# Each row is k^-1 * sum a_n x^-2n plus k^-1 times a bounded remainder.
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
base=[polynomial(nodes[:i+1],i,F(116,K)**(2*i)) for i in range(len(nodes))]
Z=[[F(z,cb['Zscale']) for z in row] for row in cb['Z']]
aux=[combine(row,base) for row in Z]
aux += [polynomial([j],0,F(ib['scales'][str(j)])) for j in ib['nodes'][len(nodes):]]
assert len(aux)==519
h=polynomial(['sinh'],0,F(1))
rest=[]
for row in rows:
 a=[F(z,BS) for z in row];z=[F(j*j,K*K) for j in js];assert max(z)<1 and sum(a)==0
 rest.append(([P(sum(y*t**n for y,t in zip(a,z))) for n in range(N)],sum(abs(y)*t**N/(1-t) for y,t in zip(a,z))))
def lcross(row):
 a,e=row
 z=isum(mul(co,J(2+2*n)) for n,co in enumerate(a))
 cap=e*((9+delta)*I(2,0)[1]+I(2,1)[1])
 return quant(mul(P(1/(A*K*K)),mul(ci,wide(z,cap))))
def gram(a,b):
 x,e=a;y,f=b
 z=isum(mul(mul(u,v),I(2+2*(i+j),0)) for i,u in enumerate(x) for j,v in enumerate(y))
 err=(sum(map(mag,x))*f+sum(map(mag,y))*e+e*f)*I(2,0)[1]
 return quant(mul(P(1/(A*K*K)),wide(z,err)))
H=gram(h,h);LH=lcross(h)
TT=[[quant(isum([mul(D,mul(x,y)),mul(H,mul(v[i],v[j])),neg(mul(LH,add(mul(x,v[j]),mul(v[i],y))))])) for j,y in enumerate(b)] for i,x in enumerate(b)]
def mixed(row):
 l=lcross(row);g=gram(row,h)
 return [quant(add(mul(l,x),neg(mul(g,y)))) for x,y in zip(b,v)]
U=[mixed(row) for row in aux];W=[mixed(row) for row in rest]
# Frozen parent difference is a comparison target only.
fd=load('return071/RETURN071_contract1.json','bf90b19349604c1de0d0d924eb19719afbe82291420d9237658505e34ed1d8ba')['variant']
hd=load('return071/RETURN071_hybrid.json','816388b20efc5d06bd32b85d09dca0d51ae8294b7f1501761816f9f14f661746')
O=hd['scale'];hd=hd['variants'][1]
receipt={}
for key,mat in [('tt',TT),('U',U)]:
 contained=overlap=total=0
 for i,row in enumerate(mat):
  for j,(l,u) in enumerate(row):
   a,b0=hd[key][i][j];c0_,d=fd[key][i][j]
   lo,hi=F(a-d,O),F(b0-c0_,O)
   contained+=lo<=l<=u<=hi;overlap+=lo<=u and l<=hi;total+=1
 receipt[key]={'entries':total,'contained':contained,'overlap':overlap,'disjoint':total-overlap}
 print(key,receipt[key],flush=True)
# Restore cross on the strengthened tail, without assuming it vanishes.
def sqrtup(x):
 assert x>=0
 return F(math.isqrt(x.numerator*S*S//x.denominator)+1,S)
cap=[[lam*sqrtup(max(F(0),gram(row,row)[1])*max(F(0),TT[j][j][1])) for j in range(87)] for row in rest]
# Both independent direct interval and Cauchy enclosure are valid; intersect.
rest_new=[]
for i,row in enumerate(W):
 rr=[]
 for j,a in enumerate(row):
  x=mul(P(lam),a);e=cap[i][j];z=max(x[0],-e),min(x[1],e)
  assert z[0]<=z[1],(i,j)
  rr.append(z)
 rest_new.append(rr)
def bank(a):return [[enc(x) for x in row] for row in a]
def scale(a,f):return [[mul(P(f),x) for x in row] for row in a]
(out/'asymptotic_bank.json').write_text(json.dumps({'K':K,'scale':S,'lambda':str(lam),'old':{'tt':bank(TT),'U':bank(U),'restoration_cross':bank(W)},'signed_reweighting':{'tt':bank(scale(TT,lam-1)),'U':bank(scale(U,lam-1))},'new_restoration_cross':bank(rest_new)},separators=(',',':')))
receipt['restoration']={'entries':696,'cauchy_intersections':696,'scope':'direct model cross and independent Cauchy cap; not parent reserve comparison'}
receipt['scope']='shared drift inequality and g(R) primitive; no terminal or actual-source approximation claim'
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
assert all(receipt[k]['contained']==receipt[k]['entries'] for k in ('tt','U')),receipt
print('RETURN090-B ASYMPTOTIC CONTAINMENT PASS',flush=True)
