"""Exact fifth moment overlap mesh: coalesced two-step translations at both arguments."""
from fractions import Fraction as F
from pathlib import Path
from bisect import bisect_right
from collections import defaultdict
import runpy,sys
p=Path(sys.argv[1]).resolve()/'cr_alp_058'
old=(p/'mesh.txt').read_bytes()
m=runpy.run_path(str(p/'mesh.py'))
assert (p/'mesh.txt').read_bytes()==old
d,sh,knots,cmax,enc=(m[x] for x in ('d','sh','knots','cmax','enc'))
assert 2*sum(h for l,h in m['weights'])<8
moves=[(j,s,s*x) for j,x in enumerate(sh) for s in (-1,1)]
vk=sorted({F(0),d,*[k-z for k in knots for _,_,z in moves if 0<k-z<d]})
wk=sorted({F(0),d,*[k-z for k in vk for _,_,z in moves if 0<k-z<d]})
def active(l,r,t):
    mid=(l+r)/2+t
    outer=[];groups=defaultdict(list)
    for j,s,z in moves:
        if not 0<mid+z<d:
            assert r+t+z<=0 or l+t+z>=d;continue
        assert 0<=l+t+z<=r+t+z<=d
        outer.append((j,s))
        for k,b,w in moves:
            a=z+w
            if not 0<mid+a<d:
                assert r+t+a<=0 or l+t+a>=d;continue
            i=bisect_right(knots,mid+a)-1
            assert knots[i]<=l+t+a and r+t+a<=knots[i+1]
            groups[a,i].append((j,s,k,b))
    return outer,dict(groups)
def emit(f,outer,groups):
    f.write(str(len(outer))+' '+' '.join(f'{j} {s}' for j,s in outer)+' '+str(len(groups))+'\n')
    for (z,i),paths in sorted(groups.items()):
        f.write(enc(z)+f' {i} {len(paths)} '+' '.join(f'{j} {s} {k} {t}' for j,s,k,t in paths)+'\n')
counts=[];raw=unique=0
with (p/'fifth_mesh.txt').open('w') as f:
    for t in sh:
        end=(d-t)/2
        cuts=sorted({F(0),end,*[k-z for k in wk for z in (F(0),t) if 0<k-z<end]})
        panels=[]
        def panel(l,r,L,Rt):
            mid=(l+r)/2;h=(r-l)/2;R=4*h
            args=[mid+z for z,i in L[1]]+[mid+t+z for z,i in Rt[1]]
            if h<=F(1,128) and all(min(y,d-y)>=2*R or cmax*R<=F(1,2) for y in args):
                panels.append((mid,h,L,Rt));return
            panel(l,mid,L,Rt);panel(mid,r,L,Rt)
        for l,r in zip(cuts,cuts[1:]):
            panel(l,r,active(l,r,F(0)),active(l,r,t))
        assert sum(2*h for _,h,_,_ in panels)==end
        counts.append(len(panels));f.write(str(len(panels))+'\n')
        for mid,h,L,Rt in panels:
            f.write(enc(mid)+' '+enc(h)+'\n')
            for outer,groups in (L,Rt):
                emit(f,outer,groups);raw+=sum(map(len,groups.values()));unique+=len(groups)
# Each w0 complex value <1e44, product <1e88; degree-191 Gauss error.
assert 4*d*F(4,3)*10**88/4**192<F(1,10**12)
print('PASS fifth geometry:',counts,'panels;',96*raw,'raw paths;',96*unique,'coalesced evaluations',flush=True)
