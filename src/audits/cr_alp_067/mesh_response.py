"""Exact half-domain geometry for a finite response through sine index 2048."""
from fractions import Fraction as F
from pathlib import Path
from bisect import bisect_right
import runpy,sys
p=Path(sys.argv[1]).resolve()/'cr_alp_058'
old=(p/'mesh.txt').read_bytes()
m=runpy.run_path(str(p/'mesh.py'))
assert (p/'mesh.txt').read_bytes()==old
d,sh,knots,cmax,enc=(m[x] for x in ('d','sh','knots','cmax','enc'))
assert 2*sum(h for l,h in m['weights'])<8
end=d/2
moves=[(j,s,s*x) for j,x in enumerate(sh) for s in (-1,1)]
cuts=sorted({F(0),end,*[k for k in knots if 0<k<end],*[-z for _,_,z in moves if 0<-z<end],*[d-z for _,_,z in moves if 0<d-z<end]})
panels=[]
def panel(l,r,ip,active):
    mid=(l+r)/2;h=(r-l)/2;R=4*h
    if h<=F(1,256) and (min(mid,d-mid)>=2*R or cmax*R<=F(1,2)):
        assert knots[ip]<=l<r<=knots[ip+1]
        assert all(0<=l+z<r+z<=d for _,_,z in active)
        panels.append((mid,h,ip,active));return
    panel(l,mid,ip,active);panel(mid,r,ip,active)
for l,r in zip(cuts,cuts[1:]):
    mid=(l+r)/2
    ip=bisect_right(knots,mid)-1
    active=[(j,s,z) for j,s,z in moves if 0<mid+z<d]
    assert all(r+z<=0 or l+z>=d for j,s,z in moves if not 0<mid+z<d)
    panel(l,r,ip,active)
assert sum(2*h for _,h,_,_ in panels)==end
# Complex sine magnitude exp(2048*(10/7)/64)<exp(46)<10^21.
# Trial coefficient sums <100, translated sum <800, so Cx<1e25 on disks.
# Original proxy <1e40; y and 2r-Cx <3e40. All products <1e82.
assert F(2048*10,7*64)<46
assert 4*d*F(4,3)*10**82/4**192<F(1,10**12)
# Shift displacement of trial: amplitude <100, derivative <100*3000.
# Zero extension is continuous (finite sine sum), so L2 shift error is
# at most 2*300000*1e-70 times paired weight sum <8.
assert 8*2*300000*F(1,10**70)<F(1,10**13)
rows=[str(len(panels))]
for mid,h,ip,active in panels:
    rows.append(enc(mid)+' '+enc(h)+' '+str(ip)+' '+str(len(active))+' '+' '.join(f'{j} {s}' for j,s,_ in active))
(p/'response_mesh.txt').write_text('\n'.join(rows)+'\n')
print('PASS exact response geometry:',len(panels),'panels',flush=True)
