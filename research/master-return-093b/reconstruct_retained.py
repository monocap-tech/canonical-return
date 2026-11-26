"""RETURN093-B: alternate scalar evaluation of original retained energy.
Use the authenticated RETURN048 analytic bilinear formula, but independently
evaluate the original exponential basis, digamma/trigamma and decay sums.
No frozen energy matrices used during reconstruction. The analytic
Gamma/pole/prime identity and interval elementary library remain shared.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache as cache
import mpmath as mp
from mpmath import iv
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();bank=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
rawfile=root/'return059/RETURN059_raw_receipt.json'
assert sha(rawfile)=='2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8'
raw=load(rawfile);ledger=load(root/'return047/RETURN047_payment_receipt.json')
old=load(root/'return048/RETURN048_retained_receipt.json')
assert sha(root/'return048/return048_retained.py')=='eaaced4a5a9f5eadd5d428a057160aea717b0ed1daf77b8d3acd0a4f462f9834'
iv.dps=100;mp.mp.dps=110;S=10**75;N=128;P=32
a=iv.mpf(5)/4;D=2*a;zero=(0,0);one=(1,0)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
@cache(None)
def z(x):return iv.mpc(iv.mpf(x[0])/2,2*iv.pi*x[1]/5)
@cache(None)
def exa(x):
 # exact quarter-turn at physical boundary
 return iv.exp(iv.mpf(x[0])*a/2)*[iv.mpc(1),iv.mpc(0,1),iv.mpc(-1),iv.mpc(0,-1)][x[1]%4]
@cache(None)
def I(x):return iv.mpc(D) if x==zero else (exa(x)-exa(neg(x)))/z(x)
@cache(None)
def Ip(x):return iv.mpc(0) if x==zero else (a*(exa(x)+exa(neg(x)))-I(x))/z(x)
bern=[iv.mpf(str(mp.bernfrac(2*j)[0]))/int(mp.bernfrac(2*j)[1]) for j in range(1,P+1)]
def pad(v,e):
 t=iv.mpf([-e.b,e.b]);return iv.mpc(v.real+t,v.imag+t)
@cache(None)
def psi(x):
 v=(iv.mpf(5)/2+z(x))/2;w=v+N
 assert w.real.a>=N
 y=iv.ln(w)-1/(2*w)-sum(b/(2*j*w**(2*j)) for j,b in enumerate(bern,1))-sum(1/(v+j) for j in range(N))
 # Euler-Maclaurin periodic Bernoulli remainder, bounded by real part.
 return pad(y,abs(bern[-1])/(2*P*iv.mpf(N)**(2*P)))
@cache(None)
def tri(x):
 v=(iv.mpf(5)/2+z(x))/2;w=v+N
 assert w.real.a>=N
 y=1/w+1/(2*w*w)+sum(b/w**(2*j+1) for j,b in enumerate(bern,1))+sum(1/(v+j)**2 for j in range(N))
 return pad(y,abs(bern[-1])/iv.mpf(N)**(2*P+1))
@cache(None)
def ss(x,y):return tri(x)/4 if x==y else (psi(y)-psi(x))/(2*z(sub(y,x)))
# Direct product form, no subtraction of exponentially weighted simple sums.
betas=[iv.mpf(2*m)+iv.mpf(1)/2 for m in range(1,49)]
BB=iv.mpf(2*49)+iv.mpf(1)/2
weights=[iv.exp(-b*D) for b in betas]
errdec=iv.exp(-BB*D)/((BB-iv.mpf(1)/2)**2*(1-iv.exp(-2*D)))
@cache(None)
def denominators(x):return [1/(b+z(x)) for b in betas]
@cache(None)
def dec(x,y):
 return pad(sum(w*u*v for w,u,v in zip(weights,denominators(x),denominators(y))),errdec)
@cache(None)
def triangle(x,y):
 den=add(y,one)
 if den==zero:return Ip(sub(x,one))+a*I(sub(x,one))
 return (I(add(x,y))-exa(neg(den))*I(sub(x,one)))/z(den)
G0=-iv.euler-iv.pi/2-3*iv.ln(2)-iv.ln(iv.pi)
@cache(None)
def gpart(x):return (psi(x)+psi(neg(x)))/2-psi(zero)
shifts=[(iv.ln(n),iv.ln(p)/iv.sqrt(n)) for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]]
assert iv.ln(11).b<D.a and iv.ln(13).a>D.b
@cache(None)
def shiftI(x,j):
 s,w=shifts[j];v=z(x)*(a-s/2)
 return iv.mpc(D-s) if x==zero else (iv.exp(v)-iv.exp(-v))/z(x)
@cache(None)
def pair(x,j):
 s,w=shifts[j];v=z(x)*s/2
 return w*(iv.exp(v)+iv.exp(-v))
@cache(None)
def entry(x,y):
 u=add(x,y);v=sub(y,x);mass=I(u)
 # Symmetrize the analytic Gamma representation without relying on
 # cancellation of independently rounded evaluations.
 gamma=(G0+gpart(x)+4)*mass-triangle(x,y)-triangle(y,x)
 gamma+=exa(neg(u))*ss(x,neg(y))+exa(u)*ss(neg(x),y)
 gamma-=exa(v)*dec(x,neg(y))+exa(neg(v))*dec(neg(x),y)
 prime=sum(pair(v,j)*shiftI(u,j) for j in range(len(shifts)))
 pole=I(add(x,one))*I(sub(y,one))+I(sub(x,one))*I(add(y,one))
 return gamma-prime+pole
def endpoints(v):
 def r(t):
  s,m,e,b=t;return F((-1 if s else 1)*m)*F(2)**e
 l,h=map(r,v._mpi_);return [math.floor(l*S),math.ceil(h*S)]
def modes(k,par):
 return [((0,k),iv.mpc(iv.mpf(1)/2) if par==0 else iv.mpc(0,-iv.mpf(1)/2)),((0,-k),iv.mpc(iv.mpf(1)/2) if par==0 else iv.mpc(0,iv.mpf(1)/2))]
def arr(x):return np.asarray(x,dtype=object)
ks=sorted(set(ledger['even_exceptional']+ledger['unpaid_odd']))
assert len(ks)==85 and ks==sorted(-k for k in ks)
assert len(raw['ks'])==44
# Construct actual source basis in RETURN048 order followed by 44
# retained even modes, as specified by authenticated RETURN086 capture.
records=[];r86=root/'return086';sys.path.insert(0,str(r86))
from check_terminal import check
for par in (0,1):
 cols=([[(zero,iv.mpc(1))]] if par==0 else [])+[modes(k,par) for k in ks if k>0]
 if par:cols.append([((1,0),iv.mpc(iv.mpf(1)/2)),((-1,0),iv.mpc(-iv.mpf(1)/2))])
 assert len(cols)==43
 cols += [modes(k,par) for k in raw['ks']]
 assert len(cols)==87
 H=np.zeros((87,87,2),dtype=object)
 for i,c in enumerate(cols):
  for j in range(i+1):
   d=cols[j];v=sum(a*b*entry(x,y) for x,a in c for y,b in d)
   assert v.imag.a<=0<=v.imag.b,(par,i,j)
   H[i,j]=H[j,i]=endpoints(v.real)
  if i%8==0:print('retained energy rows',par,i+1,flush=True)
 # Reference capture only now; check entrywise containment before any
 # terminal congruence. Compare at exact original scale, no rounding.
 for suffix in ('','_alt'):
  name=f'matrices_{par}{suffix}.json';p=bank/('reconstructed_'+name)
  d=load(p);cap=load(r86/name)
  assert d['HS']==cap['HS'] and d['Hraw']==cap['Hraw']
  assert np.array_equal(arr(d['C'])*raw['T'],arr(raw['blocks'][par]['C'])*d['CS'])
  oldH=arr(d['Hraw']);lo=H[:,:,0]*d['HS'];hi=H[:,:,1]*d['HS']
  contained=(lo>=oldH[:,:,0]*S)&(hi<=oldH[:,:,1]*S)
  overlap=(lo<=oldH[:,:,1]*S)&(hi>=oldH[:,:,0]*S)
  print('retained containment',par,suffix,int(contained.sum()),contained.size,'overlap',int(overlap.sum()),flush=True)
  if not np.all(contained):
   ix=np.argwhere(~contained)[:5];print('retained mismatch samples',[(int(i),int(j),H[i,j].tolist(),oldH[i,j].tolist()) for i,j in ix],flush=True)
  assert np.all(overlap),(par,suffix)
  # Independent enclosure, never intersect with archived matrix.
  d['Hraw']=H.tolist();d['HS']=S
  q=out/('reconstructed_'+name);q.write_text(json.dumps(d,separators=(',',':')))
  rec=check(q);rec.update(parity=par,variant=suffix or 'generator',contained=int(contained.sum()),entries=contained.size,source_sha256=sha(p));records.append(rec)
 (out/f'H_{par}.json').write_text(json.dumps({'scale':S,'parity':par,'Hraw':H.tolist()},separators=(',',':')))
 # Free cached entry and decay products before second parity if necessary.
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'alternate scalar retained energy evaluation; shared analytic form identity and mpmath elementary interval operations; all four fixed terminal checks','records':records},indent=2)+'\n')
print('RETURN093-B INDEPENDENT RETAINED ENERGY AND FOUR TERMINALS PASS',flush=True)
