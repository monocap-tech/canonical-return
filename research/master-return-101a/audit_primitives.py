"""RETURN101-A: independent Arb primitive and complete direct-sum audit.
Original direct-sum receipt is only a comparison target.
Uses python-flint/Arb, not the historical mpmath phase evaluator.
"""
import sys,json,hashlib,time,math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,acb,fmpq,ctx
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
mode=sys.argv[3]
sys.set_int_max_str_digits(0)
ctx.prec=320
def A(v):
 v=F(v);return arb(fmpq(v.numerator,v.denominator))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads((root/p).read_bytes())
def inside(z,l,u):return A(l)<z and z<A(u)
def outward(z,s=10**60):
 # use exact dyadic endpoints from Arb's floor/ceil, converted via unique integer
 l=(z.lower()*s).floor().unique_fmpz()
 u=(z.upper()*s).ceil().unique_fmpz()
 assert l is not None and u is not None
 return [int(l),int(u)]
pi=arb.pi();ns=[(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
logs=[arb(n).log() for n,p in ns];weights=[arb(p).log()/arb(n).sqrt() for n,p in ns]
if mode=='constants':
 b=load('return062/RETURN062_bank_receipt.json')
 pl,pu=map(F,b['pi']);cl,cu=map(F,b['c'])
 assert inside(pi,pl,pu)
 assert inside(A(F(25,16))/pi**2,cl,cu)
 assert 2*sum(weights)<9
 keys={F(n)**s for n,p in ns for s in [-1,1]}
 keys=(keys|{x*y for x in keys for y in keys})-{F(1)}
 phases=[]
 for z in sorted(keys):
  t=A(z).log()*A(F(2,5))
  m=round(float(t.mid()))
  assert inside(t,F(m)-F(1,2),F(m)+F(1,2))
  gap=abs(t-m)
  assert gap>A(F(1,168)),(z,gap)
  phases.append({'key':str(z),'distance_units_1e60':outward(gap),'chord_units_1e60':outward(2*abs((pi*t).sin()))})
 g0=acb(A(F(1,4))).digamma().real-pi.log()
 assert inside(g0,-6,-5)
 # Reference primitives imported only after independent Arb values exist.
 anchors=[]
 for k in [16384,32768,262144]:
  r=2*pi*k/5
  g=acb(A(F(1,4)),r/2).digamma().real-pi.log()
  d=g+g0+4
  assert abs(d)<10
  anchors.append({'k':k,'D0':d,'units_1e60':outward(d)})
 sys.path[:0]=[str(root/'return065'),str(root/'return073')]
 import return065_fast as fast
 from mpmath import iv
 src=fast.src
 for item in anchors:
  k=item['k'];r=2*iv.pi*k/5
  gR=src.G0+fast.sums((0,k))[0]+4*r*r/(r*r+iv.mpf(1)/4)
  l,u=src.endpoints(gR+src.G0+4,10**55)
  assert inside(item['D0'],F(l,10**55),F(u,10**55)),(k,l,u,item['D0'])
  del item['D0'];item['historical_units_1e55']=[l,u]
  print('anchor D0 strict containment',k,flush=True)
 rec={'status':'PASS','scope':'independent Arb constants, all distinct first/quadratic phase keys and three literal historical anchor enclosures','phase_keys':len(phases),'distance_threshold':'1/168','phases':phases,'anchors':anchors,'pi_units_1e60':outward(pi),'sinh_node_units_1e60':outward(A(F(25,16))/pi**2),'g0_units_1e60':outward(g0)}
 (out/'constants.json').write_text(json.dumps(rec,indent=2)+'\n')
 print('RETURN101-A CONSTANTS PASS',len(phases),flush=True)
elif mode=='direct':
 # Independent complete directed direct summation; no use of Abel,
 # historical recurrence, direct-sum reserve or frozen scalar endpoints
 # until the final containment comparisons.
 K=32768;U=262144;split=65536
 powers=list(range(2,47))
 sums=[[[arb(0) for _ in range(10)] for p in powers] for band in range(2)]
 start=time.monotonic()
 for k in range(K+2,U+1,2):
  r=2*pi*k/5;x=A(F(k,K));l=x.log();s=arb(0);c=arb(0)
  for w,log in zip(weights,logs):
   si,co=(r*log).sin_cos();s+=2*w*si;c+=2*w*co
  vals=[arb(1),l,l*l,s,c,l*s,l*c,s*s,s*c,c*c]
  u=1/x;v=u*u;block=sums[int(k>split)]
  for row in block:
   for j,a in enumerate(vals):row[j]+=v*a
   v*=u
  if (k-K)%16384==0:print('direct scalar rows',(k-K)//2,'seconds',round(time.monotonic()-start,2),flush=True)
 oldpath=root/'return086/direct_sums.json';old=json.loads(oldpath.read_bytes())
 records=[];new=[];fails=[]
 for name,count in [('full',114688),('restricted',98304)]:
  bank=next(x for x in old['bands'] if x['name']==name)
  assert bank['count']==count and bank['p']==powers
  rows=sums[1] if name=='restricted' else [[a+b for a,b in zip(x,y)] for x,y in zip(*sums)]
  enc=[];contained=0
  for p,row,ref in zip(powers,rows,bank['sums']):
   rr=[]
   for j,(z,(l,u)) in enumerate(zip(row,ref)):
    passed=inside(z,F(l,bank['scale']),F(u,bank['scale']))
    contained+=int(passed);rr.append(outward(z))
    if not passed:fails.append([name,p,j,str(z),l,u])
   assert len(rr)==10;enc.append(rr)
  records.append({'band':name,'coordinates':count,'entries':450,'contained':contained})
  new.append({'name':name,'count':count,'p':powers,'scale':10**60,'sums':enc})
  print('direct containment',records[-1],flush=True)
 ok=all(r['contained']==r['entries'] for r in records)
 (out/'direct.json').write_text(json.dumps({'status':'PASS' if ok else 'NOT_CERTIFIED','scope':'independent Arb direct sums, both finite bands, all ten scalar families and exponents 2..46','historical_sha256':sha(oldpath),'records':records,'mismatch_samples':fails[:20],'bands':new},indent=2)+'\n')
 assert ok,fails[:3]
 print('RETURN101-A DIRECT SCALAR CONTAINMENT PASS',900,flush=True)
else:raise ValueError(mode)
