"""RETURN097-A: alternate shared-scalar evaluation for every grouped mixed row.
Reuses the RETURN066 raw coefficient intervals already crosschecked in 085.
No import of historical Gamma, grouped-row or contraction evaluator.
"""
import sys,json,hashlib,math,time
from pathlib import Path
from fractions import Fraction as F
from mpmath import iv
import mpmath as mp
import numpy as np
sys.set_int_max_str_digits(0)
iv.dps=100;mp.mp.dps=110
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p,h=None):
 p=root/p
 if h:assert sha(p)==h,p
 return json.loads(p.read_bytes())
w=load('return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
raw=load('return059/RETURN059_raw_receipt.json','2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8')
mom=load('return066/RETURN066_moments.json','b9aa229e7147b8d47d56b0cbfec290255af557cfb99260624d1c4b576369d184')
Q=10**55;S=10**12;T=raw['T']
if 'scale' in mom:print('moment scale',mom['scale'],flush=True)
ks0=sorted(j for j in set(w['even_exceptional']+w['unpaid_odd']) if j>0)
assert len(ks0)==42 and len(raw['ks'])==44
charts=[]
for par in (0,1):
 ids=([(0,0)] if par==0 else [])+[(0,j)for j in ks0]+([(1,0)] if par else [])+[(0,j)for j in raw['ks']]
 assert len(ids)==87
 C=np.array(raw['blocks'][par]['C'],dtype=object);assert C.shape==(87,87)
 cf=np.array(mom['blocks'][par]['raw_coeff'],dtype=object);assert cf.shape==(87,5,2)
 charts.append((ids,cf[:,1:,:],np.maximum(C,0),np.minimum(C,0)))
def rational(t):
 s,m,e,b=t;return F((-1 if s else 1)*m)*F(2)**e
def enc(v,s=Q):
 a,b=map(rational,v._mpi_);return math.floor(a*s),math.ceil(b*s)
def point(x):
 x=F(x);return iv.mpf(x.numerator)/x.denominator
# Exact Bernoulli recurrence, independent of mp.bernfrac.
B=[F(1)]
for n in range(1,65):
 B.append(-sum(F(math.comb(n+1,k))*B[k] for k in range(n))/F(n+1))
N=128;P=32
bern=[point(B[2*j]) for j in range(1,P+1)]
err=point(abs(B[64])/F(64*N**64))
pad=iv.mpf([-err.b,err.b])
def psi(r):
 v=iv.mpc(iv.mpf(5)/4,r/2);z=v+N
 zinv=1/z;z2=zinv*zinv;t=z2
 acc=iv.mpc(0)
 for j,b in enumerate(bern,1):
  acc+=b*t/(2*j);t*=z2
 y=iv.ln(z)-zinv/2-acc-sum(1/(v+j) for j in range(N))
 return iv.mpc(y.real+pad,y.imag+pad)
psi0=psi(iv.mpf(0)).real
g0=-iv.euler-iv.pi/2-3*iv.ln(2)-iv.ln(iv.pi)
primes=[(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
shifts=[(iv.ln(n),iv.ln(p)/iv.sqrt(n)) for n,p in primes]
D=iv.mpf(5)/2
betas=[iv.mpf(2*m)+iv.mpf(1)/2 for m in range(1,49)]
weights=[iv.exp(-b*D) for b in betas]
bb=iv.mpf(197)/2
eerr=iv.exp(-bb*D)/(bb*bb*(1-iv.exp(-2*D)))
ferr=iv.exp(-bb*D)/(bb*(1-iv.exp(-2*D)))
cinv=enc(1/(2*iv.pi/5)**2)
def shared(k):
 r=2*iv.pi*k/5;r2=r*r;z=psi(r);g=z.real-psi0;h=z.imag/(2*r)
 e=sum(w/(b*b+r2) for b,w in zip(betas,weights))+iv.mpf([0,eerr.b])
 f=sum(w*b/(b*b+r2) for b,w in zip(betas,weights))+iv.mpf([0,ferr.b])
 gf=g0+g+4*r2/(r2+iv.mpf(1)/4)
 ps=sum(w*iv.sin(r*s)for s,w in shifts);pc=sum(w*iv.cos(r*s)for s,w in shifts)
 inv=1/(r2+iv.mpf(1)/4)
 fac=[[enc(x)for x in (2*(r2*(h-e)+r*ps),2*(gf-g/2+f-pc),iv.mpf(2),inv)],
      [enc(x)for x in (2*r*(-gf+g/2-f+pc),2*r*(h-e+ps/r),2*r,r*inv)]]
 return fac,cinv,enc(inv)
def mul(a,b):
 z=[x*y for x in a for y in b];return min(z),max(z)
def calc(par,k,fac,ci,pi):
 ids,cf,cp,cn=charts[par];L=[];U=[]
 for (typ,j),row in zip(ids,cf):
  if typ:den=pi
  else:
   d=k*k-j*j;assert d!=0
   a,b=ci if d>0 else (-ci[1],-ci[0]);d=abs(d)
   den=(a//d,-((-b)//d))
  l=h=0
  for n,(coef,scalar) in enumerate(zip(row,fac)):
   p=mul(coef,scalar)
   x,y=mul(p,den) if n<3 else (p[0]*Q,p[1]*Q)
   l+=x;h+=y
  L.append(l);U.append(h)
 L=np.array(L,dtype=object);U=np.array(U,dtype=object)
 lo=L@cp+U@cn;hi=U@cp+L@cn
 return np.stack([lo*S//(Q**3*T),-((-hi*S)//(Q**3*T))],axis=-1)
# Actual grouped tables; earlier 400 rows used different evaluation.
files=[root/f'return0{n}/RETURN0{n}_rows.json' for n in (68,69,71,79)]+sorted((root/'return081').glob('rows_*.json'))
records=[];seen=[];start=time.monotonic();total=contained=overlap=0;fails=[]
for file in files:
 bank=json.loads(file.read_bytes());assert bank['scale']==S
 assert len(bank['blocks'])==2 and all(len(b)==len(bank['ks']) for b in bank['blocks'])
 seen+=bank['ks'];local=lc=lo=0;new=[[],[]]
 for i,k in enumerate(bank['ks']):
  fac,ci,pi=shared(k)
  for par in (0,1):
   a=calc(par,k,fac[par],ci,pi);old=np.array(bank['blocks'][par][i],dtype=object)
   assert old.shape==a.shape==(87,2) and np.all(a[:,0]<=a[:,1])
   contain=(old[:,0]<=a[:,0])&(a[:,1]<=old[:,1])
   over=(old[:,0]<=a[:,1])&(a[:,0]<=old[:,1])
   local+=87;lc+=int(contain.sum());lo+=int(over.sum())
   if not np.all(over) and len(fails)<20:
    for j in np.where(~over)[0][:20-len(fails)]:fails.append([file.name,k,par,int(j),old[j].tolist(),a[j].tolist()])
   new[par].append(a.tolist())
  if (i+1)%128==0:print('rows',file.name,i+1,'seconds',round(time.monotonic()-start,2),flush=True)
 rec={'file':str(file.relative_to(root)),'sha256':sha(file),'rows':len(bank['ks']),'entries':local,'contained':lc,'overlap':lo}
 records.append(rec);total+=local;contained+=lc;overlap+=lo
 print('bank',rec,flush=True)
 (out/file.name).write_text(json.dumps({'ks':bank['ks'],'scale':S,'blocks':new},separators=(',',':'))+'\n')
expected=sorted(k for k,v in w['positive_even_weights'] if 116<k<=16384)[400:]
assert seen==expected and len(seen)==7711
receipt={'status':'PASS' if contained==total else 'OVERLAP_ONLY' if overlap==total else 'DISJOINT','scope':'alternate scalar regeneration and raw-to-chart contraction of every grouped row using frozen raw coefficient intervals; early 400 rows and compressed features excluded','rows':len(seen),'entries':total,'contained':contained,'overlap':overlap,'failures':fails,'records':records}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('RETURN097-A',receipt['status'],'containment',contained,total,'overlap',overlap,flush=True)
assert contained==total,('noncontained rows; no unconditional transfer',fails[:3])
