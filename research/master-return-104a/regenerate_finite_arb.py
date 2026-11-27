"""RETURN104-A: independent Arb regeneration of all finite mixed rows and
compressed scalar features. No mpmath in the generating evaluator.
Original interval files and early reference reader used only for comparison.
"""
import sys,json,hashlib,math,time,zipfile
from pathlib import Path
from fractions import Fraction as F
from flint import arb,acb,ctx
import numpy as np
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
p=Path(__file__).resolve().parents[1]/'master-return-101b/audit_raw_coefficients.py'
src=p.read_text();anchor='records=[];allrows=[];fail=[]'
assert src.count(anchor)==1
env={'__file__':str(p),'__name__':'arb_coefficients'}
exec(compile(src.split(anchor)[0],str(p),'exec'),env)
A=env['A'];pi=env['pi'];a=env['a'];psi0=env['psi0'];shifts=env['shifts'];bs=env['bs'];es=env['es'];pad=env['pad']
w=env['w'];raw=env['raw'];mom=env['mom'];ks0=env['ks0'];coefficient=env['coefficient']
ctx.prec=384
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x,S):
 l=(x.lower()*S).floor().unique_fmpz();u=(x.upper()*S).ceil().unique_fmpz()
 assert l is not None and u is not None
 return int(l),int(u)
def load(p):return json.loads(p.read_bytes())
with zipfile.ZipFile(root.parent/'archives/RETURN084_audit.zip') as z:
 names=[n for n in z.namelist() if n.endswith('/archive_manifest.json') or n=='archive_manifest.json'];assert len(names)==1
 data=z.read(names[0])
assert hashlib.sha256(data).hexdigest()=='38b4781ddac2b4f0437f75f8eb802e78a54c1fbc366f4690880aa73b0aa9d112'
manifest=json.loads(data);assert len(manifest)==383
def original(p):
 assert sha(p)==manifest[str(p.relative_to(root))],p
 return load(p)
g0=acb(A(F(1,4))).digamma().real-pi.log()
def phase(r):
 ps=arb(0);pc=arb(0)
 for s,w in shifts:
  si,co=(r*s).sin_cos();ps+=w*si;pc+=w*co
 return ps,pc
charts=[]
for par in (0,1):
 ids=([(0,0)] if not par else [])+[(0,k)for k in ks0]+([(1,0)]if par else [])+[(0,k)for k in raw['ks']]
 cf=[coefficient(par,*v) for v in ids];C=np.array(raw['blocks'][par]['C'],dtype=object)
 assert len(cf)==87 and C.shape==(87,87)
 charts.append((ids,cf,np.maximum(C,0),np.minimum(C,0)))
Q=10**55;S=10**12;T=raw['T']
def shared(k):
 r=2*pi*k/5;r2=r*r;z=acb(a,r/2).digamma();g=z.real-psi0;h=z.imag/(2*r)
 e=sum(t/(b*b+r2) for b,t in zip(bs,es))+pad
 f=sum(t*b/(b*b+r2) for b,t in zip(bs,es))+pad
 gf=g0+g+4*r2/(r2+A(F(1,4)));ps,pc=phase(r);p=1/(r2+A(F(1,4)))
 return [(2*(r2*(h-e)+r*ps),2*(gf-g/2+f-pc),arb(2),p),(2*r*(-gf+g/2-f+pc),2*r*(h-e+ps/r),2*r,r*p)],r2
def mixed(k):
 fac,r2=shared(k);ans=[]
 for par,(ids,cf,cp,cn) in enumerate(charts):
  x,y,z,t=fac[par];L=[];U=[]
  for (typ,j),(ell,b,v,c,d) in zip(ids,cf):
   assert typ or k!=j
   den=r2+ell
   val=(b*x+v*y+c*z)/den+d*t
   l,u=enc(val,Q);L.append(l);U.append(u)
  L=np.array(L,object);U=np.array(U,object)
  lo=L@cp+U@cn;hi=U@cp+L@cn
  ans.append(np.stack([lo*S//(Q*T),-((-hi*S)//(Q*T))],-1).tolist())
 return ans
records=[];fails=[];seen=[]
def compare(x,y,tag):
 assert len(x)==len(y)
 c=0
 for i,((l,u),(a,b)) in enumerate(zip(x,y)):
  assert l<=u and a<=b
  ok=a<=l<=u<=b;c+=int(ok)
  if not ok and len(fails)<30:fails.append([*tag,i,[l,u],[a,b]])
 return c
def audit_rows(label,ks,old):
 assert len(old)==2 and all(len(b)==len(ks) for b in old)
 new=[[],[]];count=0
 for i,k in enumerate(ks):
  row=mixed(k)
  for par in (0,1):
   count+=compare(row[par],old[par][i],[label,k,par]);new[par].append(row[par])
  if (i+1)%512==0:print('Arb mixed rows',label,i+1,flush=True)
 rec={'bank':label,'rows':len(ks),'entries':len(ks)*174,'contained':count}
 records.append(rec);seen.extend(ks)
 (out/(label+'.json')).write_text(json.dumps({'ks':ks,'scale':S,'blocks':new},separators=(',',':'))+'\n')
 print('Arb mixed containment',rec,flush=True)
expected=sorted(k for k,v in w['positive_even_weights'] if 116<k<=16384)
assert len(expected)==8111
# Early reference from accepted 097-B archive, itself already strictly
# contained in the original rephased historical reader. Comparing to
# its narrower saved bank also proves containment in the original.
early=load(Path(sys.argv[3])/'early/early_regenerated.json')
assert early['ks']==expected[:400] and early['scale']==S
audit_rows('early',early['ks'],early['blocks'])
files=[root/f'return0{n}/RETURN0{n}_rows.json' for n in (68,69,71,79)]+sorted((root/'return081').glob('rows_*.json'))
for f in files:
 b=original(f);assert b['scale']==S
 audit_rows(f.stem,b['ks'],b['blocks'])
assert seen==expected
# Compressed features: independent Arb g(R) and trig values.
old=original(root/'return080/features.json');K=16384;N=12;FS=10**10
ks=list(range(K+2,32770,2));assert old['ks']==ks and old['S']==FS
scales=[]
for par in (0,1):
 m=mom['blocks'][par]['cutoffs'][0]
 assert m['K']==32768
 scales.append([max(1,math.ceil(F(x)*4**n))for n,row in enumerate(m['norms'][:N])for x in row])
assert scales==old['scales']
R=2*pi*K/5;D0=acb(A(F(1,4)),R/2).digamma().real-pi.log()+g0+4
assert abs(D0)<10
X=[[],[]];counts=[0,0]
for i,k in enumerate(ks):
 r=2*pi*k/5;x=A(F(k,K));ps,pc=phase(r)
 aa=pi/2+2*ps;dd=D0+x.log()-2*pc
 fs=[(aa/r,dd/r**2,2/r**2),(-dd/r,aa/r**2,2/r)]
 pw=arb(1);powers=[]
 for n in range(N):powers.append(pw);pw/=x*x
 for par in (0,1):
  row=[list(enc(fs[par][j]*powers[n]*scales[par][3*n+j],FS)) for n in range(N)for j in range(3)]
  counts[par]+=compare(row,old['X'][par][i],['feature',k,par]);X[par].append(row)
 if (i+1)%2048==0:print('Arb features',i+1,flush=True)
for par in (0,1):
 records.append({'bank':'features','parity':par,'rows':len(ks),'entries':len(ks)*36,'contained':counts[par]})
(out/'features.json').write_text(json.dumps({'ks':ks,'S':FS,'scales':scales,'X':X},separators=(',',':'))+'\n')
total=sum(r['entries'] for r in records);contained=sum(r['contained'] for r in records)
assert total==2001138
ok=total==contained
(out/'receipt.json').write_text(json.dumps({'status':'PASS'if ok else 'NOT_CERTIFIED','scope':'Arb physical coefficient and complete finite mixed/feature scalar regeneration; early reference is authenticated strictly contained 097-B alternate bank; no new terminal solve','entries':total,'contained':contained,'records':records,'mismatch_samples':fails,'coefficient_source_sha256':sha(p)},indent=2)+'\n')
assert ok,fails[:3]
print('RETURN104-A ALL FINITE SOURCE INTERVALS PASS',contained,flush=True)
