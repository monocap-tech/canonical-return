"""Exact half-domain mesh for the global signed-prime second moment.
Uses the frozen 058 mesh convention, and verifies its byte regeneration.
Run with the extracted CR-ALP-058 package root as argument.
"""
from fractions import Fraction as F
from pathlib import Path
from bisect import bisect_right
import runpy, sys
p=Path(sys.argv[1]).resolve()/'cr_alp_058'
old=(p/'mesh.txt').read_bytes()
m=runpy.run_path(str(p/'mesh.py'))
assert (p/'mesh.txt').read_bytes()==old
d,sh,knots,cmax,enc=(m[x] for x in ('d','sh','knots','cmax','enc'))
assert 2*sum(h for l,h in m['weights'])<8
end=d/2
moves=[(j,s,s*x) for j,x in enumerate(sh) for s in (-1,1)]
cuts=sorted({F(0),end,*[k-z for k in knots for _,_,z in moves if 0<k-z<end]})
out=[]
def panel(l,r,active):
    mid=(l+r)/2; h=(r-l)/2; R=4*h
    if h<=F(1,128) and all(min(mid+z,d-mid-z)>=2*R or cmax*R<=F(1,2) for _,_,z,_ in active):
        assert all(knots[i]<=l+z and r+z<=knots[i+1] for _,_,z,i in active)
        out.append((mid,h,active)); return
    panel(l,mid,active); panel(mid,r,active)
for l,r in zip(cuts,cuts[1:]):
    mid=(l+r)/2
    active=[(j,s,z,bisect_right(knots,mid+z)-1) for j,s,z in moves if 0<mid+z<d]
    assert all((r+z<=0 or l+z>=d) for j,s,z in moves if not 0<mid+z<d)
    panel(l,r,active)
assert sum(2*h for _,h,_ in out)==end
# The complex disk bound for every active translated r is 1e40, inherited
# from 058's verified coefficient bounds. |P_rational r| <= 8e40.
# Degree-191 product remainder, including reflection, fits 1e-12.
assert 4*d*F(4,3)*64*10**80/4**192<F(1,10**12)
rows=[str(len(out))]
for mid,h,active in out:
    rows.append(enc(mid)+' '+enc(h)+' '+str(len(active))+' '+ ' '.join(f'{j} {s} {i}' for j,s,z,i in active))
(p/'second_mesh.txt').write_text('\n'.join(rows)+'\n')
print('PASS exact second-moment geometry:',len(out),'panels; total half-domain',end,flush=True)
