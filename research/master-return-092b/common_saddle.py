"""RETURN092-B independent common auxiliary saddle Gram.
Recompute actual finite weighted rows from the original diagonal weights,
constraint nodes and frozen rational chart; bound the infinite constant-
weight tail via convergent power expansions. No parent Gram routines.
The odd eight-row restoration block is a separate subsequent obligation.
"""
import sys,json,hashlib,math
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=True)
def load(p):return json.loads((root/p).read_bytes())
def P(x):return F(x),F(x)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
 z=[x*y for x in a for y in b];return min(z),max(z)
def isum(xs):
 l=h=F(0)
 for a,b in xs:l+=a;h+=b
 return l,h
def mag(a):return max(map(abs,a))
def enc(a,s):return math.floor(a[0]*s),math.ceil(a[1]*s)
def ceil(a,b):return -((-a)//b)
def ivdiv(a,b):
 assert b[0]*b[1]>0
 return mul(a,(1/b[1],1/b[0]))
def rescale(a,old,new):
 l=a[:,:,0]*new//old;h=ceil(a[:,:,1]*new,old)
 assert np.all(l<=h)
 return np.stack([l,h],-1)
w=load('return057/RETURN057_weights.json')
u=load('return047/RETURN047_payment_receipt.json')
b=load('return062/RETURN062_bank_receipt.json')
cr=load('return063/RETURN063_constraints_receipt.json')
iv=load('return063/RETURN063_inverse_receipt.json')
good=sorted((k,v) for k,v in w['positive_even_weights'] if k>116)
assert len(good)==16302 and len({k for k,v in good})==len(good)
assert all(k%2==0 and v>0 for k,v in good)
print('finite weight ledger',len(good),good[0],good[-1],flush=True)
# Original RETURN063 direct sum stops immediately before the analytic
# k >= 32768 tail. Reject any different archived cutoff.
assert good[-1][0]==32766
K=32768;A=F(87,500);S=10**16;RS=10**50
c0,c1=map(F,b['c'])
def node(j):return (-c1,-c0) if j=='sinh' else P(j*j)
pi=tuple(map(F,b['pi']));p2=mul(pi,pi)
pos=dict(u['odd_good']);neg={j:v for j,v in u['odd_deficit'] if v}
assert set(pos).isdisjoint(neg)
# Base 2^18 signed limbs: every integer product and sum over ALL finite
# rows is bounded by 16302*(2^18-1)^2 < 2^50, hence exact binary64.
B=2**18;assert len(good)*(B-1)**2<2**52
def digits(x,n):
 sign=np.where(x<0,-1,1);y=abs(x);ds=[]
 for i in range(n):
  ds.append(np.asarray((y%B)*sign,dtype=np.float64));y=y//B
 assert not np.any(y)
 return ds
def finite(par,cb,ib):
 cs=cb['nodes'];nc=len(cs);nodes=ib['nodes'];fs=nodes[nc:]
 assert nodes[:nc]==cs and set(fs)==set(pos)|set(neg)
 z=np.asarray(cb['Z'],dtype=object);zs=cb['Zscale'];zp=np.maximum(z,0);zn=np.minimum(z,0)
 scales=[F(ib['scales'][str(j)]) for j in fs];sn=[a.numerator*j**(1-par) for j,a in zip(fs,scales)];sd=[a.denominator for a in scales]
 f2=np.asarray([j*j for j in fs],dtype=object);sn=np.asarray(sn,dtype=object);sd=np.asarray(sd,dtype=object)
 n=len(nodes);L=6
 accum={(i,j):np.zeros((n,n),dtype=np.float64) for i in range(L) for j in range(i,L)}
 ev=np.zeros(n,dtype=object);e2=0
 for start in range(0,len(good),512):
  rows=good[start:start+512];m=len(rows)
  cl=np.empty((m,nc),object);ch=cl.copy();fl=np.empty((m,len(fs)),object);fh=fl.copy()
  for idx,(k,v) in enumerate(rows):
   q=P(116**(1-par)*k**par)
   for i,j in enumerate(cs):
    t=node(j);q=ivdiv(q,(k*k-t[1],k*k-t[0]))
    cl[idx,i],ch[idx,i]=enc(q,RS);q=mul(q,P(116**2))
   den=sd*(k*k-f2);assert np.all(den!=0)
   num=sn*k**par*S;sg=np.where(den<0,-1,1);den=abs(den);num=num*sg
   fl[idx]=num//den;fh[idx]=ceil(num,den)
  lo=(cl@zp.T+ch@zn.T)*S//(zs*RS)
  hi=ceil((ch@zp.T+cl@zn.T)*S,zs*RS)
  lo=np.concatenate([lo,fl],axis=1);hi=np.concatenate([hi,fh],axis=1)
  # Weight rows by sqrt(2000000/(198000+w_k)), rounded outward.
  rl=np.asarray([math.isqrt(2000000*S*S//(198000+v)) for k,v in rows],object)[:,None]
  rh=rl+1
  products=[lo*rl,lo*rh,hi*rl,hi*rh]
  lo=np.minimum(np.minimum(products[0],products[1]),np.minimum(products[2],products[3]))//S
  hi=ceil(np.maximum(np.maximum(products[0],products[1]),np.maximum(products[2],products[3])),S)
  assert np.all(lo<=hi)
  mid=(lo+hi)//2;rad=np.maximum(mid-lo,hi-mid);e=np.max(rad,axis=1)
  ev+=np.sum(abs(mid)*e[:,None],axis=0);e2+=sum(int(x)**2 for x in e)
  ds=digits(mid,L)
  for i in range(L):
   if not np.any(ds[i]):continue
   for j in range(i,L):
    if np.any(ds[j]):accum[i,j]+=ds[i].T@ds[j]
  if start%4096==0:print('finite Gram rows',par,start+m,flush=True)
 M=np.zeros((n,n),object)
 for (i,j),a in accum.items():
  assert np.all(np.isfinite(a)) and np.all(a==np.rint(a)) and np.max(abs(a))<2**52
  t=a.astype(np.int64).astype(object)*B**(i+j)
  M+=t if i==j else t+t.T
 E=ev[:,None]+ev[None,:]+e2
 assert np.array_equal(M,M.T)
 return np.stack([(M-E)//S,ceil(M+E,S)],-1)
# Tail coefficient polynomials with common k^{-(2-par)} factor.
N=8;TS=10**40
def poly(ts,offset,f):
 a=[P(1)]+[P(0)]*(N-1);major=[F(1)]+[F(0)]*(N-1);total=F(1)
 for j in ts:
  r0=mul(node(j),P(F(1,K*K)));r=mag(r0);assert r<1
  v=[P(1)]
  for _ in range(1,N):v.append(mul(v[-1],r0))
  a=[isum(mul(a[i],v[t-i]) for i in range(t+1)) for t in range(N)]
  major=[sum(major[i]*r**(t-i) for i in range(t+1)) for t in range(N)]
  total/=1-r
 return [P(0) if i<offset else mul(P(f),a[i-offset]) for i in range(N)],f*(total-sum(major[:max(0,N-offset)]))
@lru_cache(None)
def moment(p):
 # k=K,K+2,..., convex trapezoidal error <= p/(4K).
 assert p>1
 m=F(K,2*(p-1))+F(1,2)
 return m,m+F(p,4*K)
def interval_array(a):return np.array([[enc(x,TS) for x in row] for row in a],object)
def mm(a,b):
 l=np.zeros((len(a),b.shape[1]),object);h=l.copy()
 for k in range(a.shape[1]):
  v=[a[:,k,i,None]*b[k,None,:,j] for i in (0,1) for j in (0,1)]
  l+=np.minimum(np.minimum(v[0],v[1]),np.minimum(v[2],v[3]))
  h+=np.maximum(np.maximum(v[0],v[1]),np.maximum(v[2],v[3]))
 return np.stack([l//TS,ceil(h,TS)],-1)
def tail(par,cb,ib):
 cs=cb['nodes'];z=[[F(a,cb['Zscale']) for a in row] for row in cb['Z']]
 base=[poly(cs[:i+1],i,F(116**(1-par))*F(116,K)**(2*i)) for i in range(len(cs))]
 C=[[isum(mul(P(a),v[t]) for a,(v,e) in zip(row,base)) for t in range(N)] for row in z]
 E=[sum(abs(a)*e for a,(v,e) in zip(row,base)) for row in z]
 for j in ib['nodes'][len(cs):]:
  a=F(ib['scales'][str(j)])*j**(1-par);v,e=poly([j],0,a);C.append(v);E.append(e)
 p=2*(2-par);den=A*K**p
 D=[[mul(P(1/den),moment(p+2*(i+j))) for j in range(N)] for i in range(N)]
 X=interval_array(C);G=mm(mm(X,interval_array(D)),X.transpose(1,0,2))
 er=np.array([math.ceil(e*TS) for e in E],object)
 M=np.array([[math.ceil(mag(a)*TS) for a in row] for row in C],object)
 ww=np.array([math.ceil(moment(p+2*(N+i))[1]*TS/den) for i in range(N)],object)
 e2=math.ceil(moment(p+4*N)[1]*TS/den);v=M@ww
 error=ceil(v[:,None]*er[None,:]+er[:,None]*v[None,:]+e2*er[:,None]*er[None,:],TS*TS)
 G[:,:,0]-=error;G[:,:,1]+=error
 return rescale(G,TS,S)
def check_overlap(a,l,h,sc):
 lower=a[:,:,0]*sc;upper=a[:,:,1]*sc;lo=l*S;hi=h*S
 inside=(lower>=lo)&(upper<=hi);over=(lower<=hi)&(upper>=lo)
 return {'entries':int(over.size),'contained':int(inside.sum()),'overlap':int(over.sum()),'disjoint':int((~over).sum()),'max_width_units':str(max((a[:,:,1]-a[:,:,0]).flat))}
records=[]
for par in (0,1):
 cb=cr['blocks'][par];ib=iv['blocks'][par];nc=len(cb['nodes']);n=len(ib['nodes'])
 assert cb['parity']==ib['parity']==par and n==518+par
 Gf=finite(par,cb,ib);Gt=tail(par,cb,ib);G=Gf+Gt
 for i,j in enumerate(ib['nodes'][nc:],nc):
  signed=F(10**6,8*pos[j]) if j in pos else -F(10**6,8*neg[j])
  a=F(ib['scales'][str(j)])**2*signed
  G[i,i]+=np.array(enc(mul(P(a),p2),S),object)
 G[:,:,0]=np.maximum(G[:,:,0],G[:,:,0].T);G[:,:,1]=np.minimum(G[:,:,1],G[:,:,1].T)
 assert np.all(G[:,:,0]<=G[:,:,1])
 rec=check_overlap(G,np.array(ib['L'],object),np.array(ib['U'],object),iv['scale'])
 rec['parity']=par;records.append(rec);print('independent common saddle',rec,flush=True)
 (out/f'common_{par}.json').write_text(json.dumps({'scale':S,'parity':par,'nodes':ib['nodes'],'K':G.tolist(),'checks':rec},separators=(',',':')))
 # No borrowing or intersection of parent intervals; independently
 # justified endpoints must at minimum agree with the frozen parent.
 assert rec['disjoint']==0,rec
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'common 518/519 auxiliary saddle; shared original weights, pi and node intervals; odd restoration not reconstructed','records':records},indent=2)+'\n')
print('RETURN092-B BOTH COMMON AUXILIARY SADDLES PASS',flush=True)
