"""RETURN092-C: independent odd restoration Gram and complete auxiliary
replacement in the four fixed terminal comparisons.
Shared: frozen restoration metric R and original scalar weight intervals.
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
# RETURN063_constraints.py binds w to RETURN047's payment ledger
# and u to RETURN057's signed-update weights (in that order).
w=load('return047/RETURN047_payment_receipt.json')
u=load('return057/RETURN057_weights.json')
assert 'positive_even_weights' in w and 'odd_good' in u and 'odd_deficit' in u
print('source schemas',list(w),list(u),flush=True)
b=load('return062/RETURN062_bank_receipt.json')
cr=load('return063/RETURN063_constraints_receipt.json')
iv=load('return063/RETURN063_inverse_receipt.json')
good=sorted((k,v) for k,v in w['positive_even_weights'] if k>116)
assert len(good)==16302 and len({k for k,v in good})==len(good)
assert all(k%2==0 and v>0 for k,v in good)
print('finite weight ledger',len(good),good[0],good[-1],flush=True)
# Original RETURN063 direct sum stops immediately before the analytic
# k >= 32768 tail. Reject any different archived cutoff.
assert good[-1][0]==32766,('finite/tail cutoff',good[-1])
K=32768;A=F(87,500);S=10**16;RS=10**50
c0,c1=map(F,b['c'])
def node(j):return (-c1,-c0) if j=='sinh' else P(j*j)
pi=tuple(map(F,b['pi']));p2=mul(pi,pi)
pos=dict(u['odd_good']);neg={j:v for j,v in u['odd_deficit'] if v}
assert set(pos).isdisjoint(neg)
print('signed carrier counts',len(pos),len(neg),flush=True)
# Base 2^18 signed limbs: every integer product and sum over ALL finite
# rows is bounded by 16302*(2^18-1)^2 < 2^50, hence exact binary64.
B=2**18;assert len(good)*(B-1)**2<2**52
def digits(x,n):
 sign=np.where(x<0,-1,1);y=abs(x);ds=[]
 for i in range(n):
  ds.append(np.asarray((y%B)*sign,dtype=np.float64));y=y//B
 assert not np.any(y)
 return ds

# Additional input arguments after root and output.
banks=Path(sys.argv[3]).resolve()
pr=load('return077/proposal.json')
assert len(pr['js'])==500 and len(pr['rows'])==8
assert all(len(r)==500 and sum(r)==0 for r in pr['rows'])
def finite_augmented(par,cb,ib):
 assert par==1
 cs=cb['nodes'];nc=len(cs);nodes=ib['nodes'];fs=nodes[nc:]
 assert nodes[:nc]==cs and set(fs)==set(pos)|set(neg)
 z=np.asarray(cb['Z'],dtype=object);zs=cb['Zscale'];zp=np.maximum(z,0);zn=np.minimum(z,0)
 scales=[F(ib['scales'][str(j)]) for j in fs];sn=[a.numerator*j**(1-par) for j,a in zip(fs,scales)];sd=[a.denominator for a in scales]
 f2=np.asarray([j*j for j in fs],dtype=object);sn=np.asarray(sn,dtype=object);sd=np.asarray(sd,dtype=object)
 n=len(nodes)+8;L=6
 rb=np.asarray(pr['rows'],dtype=object);rp=np.maximum(rb,0);rn=np.minimum(rb,0)
 j2=np.asarray([j*j for j in pr['js']],dtype=object);bs=pr['scale']
 accum={(i,j):np.zeros((n,n),dtype=np.float64) for i in range(L) for j in range(i,L)}
 ev=np.zeros(n,dtype=object);e2=0
 for start in range(0,len(good),512):
  rows=good[start:start+512];m=len(rows)
  cl=np.empty((m,nc),object);ch=cl.copy();fl=np.empty((m,len(fs)),object);fh=fl.copy()
  hl=np.empty((m,len(j2)),object);hh=hl.copy()
  for idx,(k,v) in enumerate(rows):
   q=P(116**(1-par)*k**par)
   for i,j in enumerate(cs):
    t=node(j);q=ivdiv(q,(k*k-t[1],k*k-t[0]))
    cl[idx,i],ch[idx,i]=enc(q,RS);q=mul(q,P(116**2))
   den=sd*(k*k-f2);assert np.all(den!=0)
   num=sn*k**par*S;sg=np.where(den<0,-1,1);den=abs(den);num=num*sg
   fl[idx]=num//den;fh[idx]=ceil(num,den)
   den=k*k-j2;assert np.all(den!=0)
   num=k*RS*np.where(den<0,-1,1);den=abs(den)
   hl[idx]=num//den;hh[idx]=ceil(num,den)
  lo=(cl@zp.T+ch@zn.T)*S//(zs*RS)
  hi=ceil((ch@zp.T+cl@zn.T)*S,zs*RS)
  rl0=(hl@rp.T+hh@rn.T)*S//(bs*RS)
  rh0=ceil((hh@rp.T+hl@rn.T)*S,bs*RS)
  lo=np.concatenate([lo,fl,rl0],axis=1);hi=np.concatenate([hi,fh,rh0],axis=1)
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

def tail_augmented(par,cb,ib):
 assert par==1
 cs=cb['nodes'];z=[[F(a,cb['Zscale']) for a in row] for row in cb['Z']]
 base=[poly(cs[:i+1],i,F(116**(1-par))*F(116,K)**(2*i)) for i in range(len(cs))]
 C=[[isum(mul(P(a),v[t]) for a,(v,e) in zip(row,base)) for t in range(N)] for row in z]
 E=[sum(abs(a)*e for a,(v,e) in zip(row,base)) for row in z]
 for j in ib['nodes'][len(cs):]:
  a=F(ib['scales'][str(j)])*j**(1-par);v,e=poly([j],0,a);C.append(v);E.append(e)
 for row in pr['rows']:
  a=[F(z,pr['scale']) for z in row];r=[F(j*j,K*K) for j in pr['js']]
  assert max(r)<1 and sum(a)==0
  C.append([P(sum(z*q**n for z,q in zip(a,r))) for n in range(N)])
  E.append(sum(abs(z)*q**N/(1-q) for z,q in zip(a,r)))
  assert C[-1][0]==P(0)
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

def arr(x):return np.array(x,dtype=object)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def overlap(a,b,label):
 assert a.shape==b.shape
 ok=(a[:,:,0]<=b[:,:,1])&(b[:,:,0]<=a[:,:,1])
 print(label,int(ok.sum()),ok.size,flush=True);assert np.all(ok),label
def contain_report(a,b,label):
 assert a.shape==b.shape
 ok=(a[:,:,0]>=b[:,:,0])&(a[:,:,1]<=b[:,:,1])
 print(label,'contained',int(ok.sum()),ok.size,flush=True)
 overlap(a,b,label+' overlap')
 return int(ok.sum()),ok.size
cb=cr['blocks'][1];ib=iv['blocks'][1]
G=finite_augmented(1,cb,ib)+tail_augmented(1,cb,ib)
assert G.shape==(527,527,2)
# Only after independent reconstruction import the archived setup.
# It supplies the exact restoration metric R and reference W/Gram.
sys.path.insert(0,str(root/'return083'));sys.argv=[str(root/'return083/upgrade.py')]
import setup as parent
assert not parent.alt
print('restoration receipt keys',list(parent.r),'W scale',parent.WS,flush=True)
for p in sorted((root/'return077').glob('*.py')):
 print('restoration source',p.name,sha(p),flush=True)
 for l in p.read_text().splitlines():
  if any(x in l for x in ("['R']","'R':","'gram':","R=","R =","ridge")):print(l,flush=True)
W=rescale(arr(parent.W),parent.WS,S)
gram=np.array([[enc(tuple(map(F,a)),S) for a in row] for row in parent.r['gram']],object)
contain_report(G[:519,519:],W,'odd restoration cross')
contain_report(G[519:,519:],gram,'odd restoration Gram')
# Import independently checked common saddle instead of recomputing its
# signed metric. The restoration cross and block come only from G.
def common(par):
 d=json.loads((banks/'092b'/f'common_{par}.json').read_bytes())
 assert d['parity']==par and d['scale']==S
 a=arr(d['K']);assert a.shape==(518+par,518+par,2)
 return a
G[:519,:519]=common(1)
for i in range(8):
 for j in range(8):
  G[519+i,519+j]+=np.array(enc(P(F(parent.r['R'][i][j])),S),object)
def sym(a):
 a[:,:,0]=np.maximum(a[:,:,0],a[:,:,0].T);a[:,:,1]=np.minimum(a[:,:,1],a[:,:,1].T)
 assert np.all(a[:,:,0]<=a[:,:,1]);return a
G=sym(G)
# Independently audited constant-weight tail Gram corrections, with
# signed interval endpoints. Retain scale S, round directed at each step.
def gbank(k):
 d=json.loads((banks/'088'/f'gram_{k}.json').read_bytes());assert d['K']==k
 return rescale(arr(d['gram']),d['scale'],S)
def subtract(a,b):return np.stack([a[:,:,0]-b[:,:,1],a[:,:,1]-b[:,:,0]],-1)
def fac(a,q):return np.stack([a[:,:,0]*q.numerator//q.denominator,ceil(a[:,:,1]*q.numerator,q.denominator)],-1)
g0=gbank(262144);g2=gbank(65536)
Kodd=sym(subtract(subtract(G,fac(g0,F(1025,1199))),fac(subtract(g2,g0),F(115,173))))
Keven=common(0)
(out/'odd_augmented.json').write_text(json.dumps({'scale':S,'K':Kodd.tolist(),'scope':'independent common and restoration Gram; exact archived R; independent tail corrections'},separators=(',',':')))
# Check all four previously reconstructed TT/U inputs at fixed proposals.
r86=root/'return086';sys.path.insert(0,str(r86))
from check_terminal import check
records=[]
for par,Knew,folder in [(0,Keven,'091c'),(1,Kodd,'092a')]:
 for suffix in ('','_alt'):
  name=f'reconstructed_matrices_{par}{suffix}.json'
  p=banks/folder/name;d=json.loads(p.read_bytes())
  assert F(d['target'])==(F(49,50) if par==0 else F(19,20))
  old=rescale(np.stack([arr(d['kl']),arr(d['ku'])],-1),d['T'],S)
  overlap(Knew,old,name+' complete auxiliary overlap')
  # Use common checker unchanged at its original or finer scale.
  d['T']=S;d['kl']=Knew[:,:,0].tolist();d['ku']=Knew[:,:,1].tolist()
  q=out/name;q.write_text(json.dumps(d,separators=(',',':')))
  rec=check(q);rec['input_sha256']=sha(p);records.append(rec)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'four terminal comparisons with independent mixed source and auxiliary Gram; shared raw H and restoration metric R','records':records},indent=2)+'\n')
print('RETURN092-C ALL FOUR RECONSTRUCTED AUXILIARY TERMINALS PASS',flush=True)
