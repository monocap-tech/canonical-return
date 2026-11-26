"""RETURN097-C: independently regenerate all RETURN080 finite-band features.

The 36 features per parity use the fixed RETURN080 analytic model,
not the exact mixed source. The error bound from RETURN085 is separate.
Reuse the alternate digamma primitives of RETURN097-A, without running
its mixed-row census or importing RETURN080's evaluator.
"""
import sys,json,hashlib,math,time
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
p=Path(__file__).resolve().parent.parent/'master-return-097a/regenerate_grouped.py'
src=p.read_text()
# Reuse only the independently implemented special-function primitives.
# Its setup also authenticates the raw moment and chart inputs.
anchor='def shared(k):'
assert src.count(anchor)==1
ns={'__file__':str(p),'__name__':'alternate_feature_primitives'}
exec(compile(src.split(anchor)[0],str(p),'exec'),ns)
iv=ns['iv'];psi=ns['psi'];g0=ns['g0'];psi0=ns['psi0'];enc=ns['enc'];shifts=ns['shifts']
mom=ns['mom'];K=16384;N=12;S=10**10
oldfile=root/'return080/features.json'
# Authenticate against the accepted literal archive manifest.
import zipfile
with zipfile.ZipFile(root.parent/'archives/RETURN084_audit.zip') as z:
 names=[n for n in z.namelist() if n.endswith('/archive_manifest.json') or n=='archive_manifest.json']
 assert len(names)==1,names
 data=z.read(names[0])
assert hashlib.sha256(data).hexdigest()=='38b4781ddac2b4f0437f75f8eb802e78a54c1fbc366f4690880aa73b0aa9d112'
manifest=json.loads(data)
assert len(manifest)==383 and sha(oldfile)==manifest['return080/features.json']
print('original feature bank authenticated',sha(oldfile),flush=True)
old=json.loads(oldfile.read_bytes())
ks=list(range(K+2,32768+1,2));assert old['ks']==ks and old['S']==S and len(ks)==8192
scales=[]
for par in (0,1):
 m=mom['blocks'][par]['cutoffs'][0]
 assert len(m['norms'])>=N and all(len(x)==3 for x in m['norms'][:N])
 scales.append([max(1,math.ceil(F(x)*4**n)) for n,row in enumerate(m['norms'][:N]) for x in row])
assert scales==old['scales']
# D0 = g(R)+g(0)+4, where g(R)=g0+gpart(R)+first resolvent increment.
R=2*iv.pi*K/5
D0=2*g0+psi(R).real-psi0+4*R*R/(R*R+iv.mpf(1)/4)+4
assert max(map(abs,enc(D0,10**55)))<30*10**55
# Each row has exactly 12 moment powers and 3 shared scalar components.
# Generate independent endpoint pairs at the original scale, without
# intersecting them with the historical enclosures.
new=[[],[]];counts=[dict(entries=0,contained=0,overlap=0) for _ in range(2)];fail=[];start=time.monotonic()
for i,k in enumerate(ks):
 r=2*iv.pi*k/5;x=iv.mpf(k)/K
 ps=sum(w*iv.sin(r*s) for s,w in shifts)
 pc=sum(w*iv.cos(r*s) for s,w in shifts)
 a=iv.pi/2+2*ps;d=D0+iv.ln(x)-2*pc
 fs=((a/r,d/(r*r),2/(r*r)),(-d/r,a/(r*r),2/r))
 for par in (0,1):
  row=[]
  assert len(old['X'][par][i])==36
  for n in range(N):
   for j in range(3):
    idx=3*n+j;lo,hi=enc(fs[par][j]*x**(-2*n)*scales[par][idx],S)
    l,h=old['X'][par][i][idx];assert l<=h and lo<=hi
    contained=l<=lo and hi<=h;overlap=l<=hi and lo<=h
    c=counts[par];c['entries']+=1;c['contained']+=int(contained);c['overlap']+=int(overlap)
    if not contained and len(fail)<20:fail.append([k,par,idx,[l,h],[lo,hi],overlap])
    row.append([lo,hi])
  new[par].append(row)
 if (i+1)%1024==0:print('feature coordinates',i+1,'elapsed',round(time.monotonic()-start,2),flush=True)
assert all(c['entries']==8192*36 for c in counts)
ok=all(c['contained']==c['entries'] for c in counts)
overlap=all(c['overlap']==c['entries'] for c in counts)
(out/'features.json').write_text(json.dumps({'ks':ks,'scales':scales,'X':new,'S':S},separators=(',',':'))+'\n')
receipt={'status':'PASS' if ok else 'OVERLAP_ONLY' if overlap else 'DISJOINT','scope':'all RETURN080 modeled features, independent alternate D0/phase evaluation, original moment norms and scaling; not exact mixed rows or a new error certificate','coordinates':8192,'features_per_parity':36,'counts':counts,'D0_interval_at_10e55':enc(D0,10**55),'historical_features_sha256':sha(oldfile),'alternate_primitives_sha256':sha(p),'regenerated_features_sha256':sha(out/'features.json'),'mismatch_samples':fail}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('RETURN097-C',receipt['status'],counts,flush=True)
assert ok,('strict containment required',fail[:3])
