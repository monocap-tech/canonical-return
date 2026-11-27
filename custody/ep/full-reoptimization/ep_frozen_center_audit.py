"""Read-only audit of frozen numerical center data; not a physical certificate."""
import hashlib,zipfile,json,io,math,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
p=Path(sys.argv[1]); z=zipfile.ZipFile(p); assert z.testzip() is None
m=json.loads(z.read('MANIFEST.json'))
assert set(z.namelist())=={x['name'] for x in m['files']}|{'MANIFEST.json'}
for f in m['files']:
 b=z.read(f['name']);assert len(b)==f['size'] and hashlib.sha256(b).hexdigest()==f['sha256']
def integerize(a):
 # Every stored binary64 is a rational with a power-of-two denominator.
 pairs=[float(x).as_integer_ratio() for x in a.flat];d=max(q for _,q in pairs)
 return np.array([n*(d//q) for n,q in pairs],object).reshape(a.shape),d
results=[]
for n in (384,500):
 name=f'EP-kernel-full-reopt-block-1-{n}.npz'
 with np.load(io.BytesIO(z.read(name))) as q:
  a,b,c,t,s=[q[k] for k in ('A','B','C','T','S')]
  assert np.array_equal(a,a.T) and np.array_equal(c,c.T)
  # Freeze T exactly, avoiding any assumption on inverse-solve accuracy.
  ai,ad=integerize(a);bi,bd=integerize(b);ci,cd=integerize(c);ti,td=integerize(t)
  den=max(cd,bd*td,ad*td*td)
  bt=bi.T@ti
  fi=ci*(den//cd)+(bt+bt.T)*(den//(bd*td))-(ti.T@ai@ti)*(den//(ad*td*td))
  ff=np.array([[float(F(int(x),den)) for x in row] for row in fi]);ev,v=np.linalg.eigh(ff)
  target=F(27,10**7) if n==384 else F(39,10**7)
  # Exact congruence of the shifted center. Quantized eigenvectors are merely proposals.
  R=2**40; vi=np.rint(v*R).astype(np.int64).astype(object)
  # Nonzero determinant via ||V^T V-I||_infty < 1, exactly.
  vg=vi.T@vi-np.eye(32,dtype=object)*R**2
  assert max(sum(abs(int(x)) for x in row) for row in vg)<R**2
  h=fi*target.denominator-np.eye(32,dtype=object)*(den*target.numerator)
  g=vi.T@h@vi; gd=den*target.denominator*R**2
  margins=[int(g[i,i])-sum(abs(int(g[i,j])) for j in range(32) if j!=i) for i in range(32)]
  assert min(margins)>0
  tg=ti.T@ti; norm2=F(max(sum(abs(int(x)) for x in row) for row in tg),td*td)
  # Rational norm bound rounded upward to 10^-9 grid.
  k=10**9; rad=math.isqrt((norm2.numerator*k*k)//norm2.denominator)+1; norm=F(rad,k);assert norm*norm>norm2
  results.append(dict(coordinates=n,center_min_eigenvalue=float(ev[0]),stored_center_max_difference=float(np.max(abs(ff-s))),exact_center_lower=str(target),congruence_min_margin=float(F(min(margins),gd)),trial_operator_norm_upper=str(norm),trial_operator_norm_upper_decimal=float(norm),error_B_factor=float(2*norm),error_A_factor=float(norm*norm)))
  print(results[-1],flush=True)
out=dict(status='exact dyadic audit of stored numerical centers only; no physical error enclosures or new positivity certification',archive_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),manifest_payloads=len(m['files']),results=results)
Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')
