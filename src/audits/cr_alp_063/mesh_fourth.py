"""Exact half-domain mesh for A v0, A=P_rat, with coalesced double paths."""
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
end=d/2
cuts=sorted({F(0),end,*[k-z for k in vk for _,_,z in moves if 0<k-z<end]})
panels=[]
def panel(l,r,outer,groups):
    mid=(l+r)/2;h=(r-l)/2;R=4*h
    if h<=F(1,128) and all(min(mid+z,d-mid-z)>=2*R or cmax*R<=F(1,2) for z,i in groups):
        panels.append((mid,h,outer,groups));return
    panel(l,mid,outer,groups);panel(mid,r,outer,groups)
for l,r in zip(cuts,cuts[1:]):
    mid=(l+r)/2
    outer=[];groups=defaultdict(list)
    for j,s,z in moves:
        if not 0<mid+z<d:
            assert r+z<=0 or l+z>=d;continue
        assert 0<=l+z<=r+z<=d
        outer.append((j,s))
        for k,t,w in moves:
            a=z+w
            if not 0<mid+a<d:
                assert r+a<=0 or l+a>=d;continue
            i=bisect_right(knots,mid+a)-1
            assert knots[i]<=l+a and r+a<=knots[i+1]
            groups[a,i].append((j,s,k,t))
    panel(l,r,outer,dict(groups))
assert sum(2*h for _,h,_,_ in panels)==end
# |A v0| < 8*10^42 < 10^43 on each analytic disk.
assert 4*d*F(4,3)*10**86/4**192<F(1,10**12)
rows=[str(len(panels))];raw=unique=0
for mid,h,outer,groups in panels:
    rows.append(enc(mid)+' '+enc(h)+' '+str(len(outer))+' '+' '.join(f'{j} {s}' for j,s in outer)+' '+str(len(groups)))
    for (z,i),paths in sorted(groups.items()):
        # Store the rational total displacement directly, not a floating key.
        rows.append(enc(z)+f' {i} {len(paths)} '+' '.join(f'{j} {s} {k} {t}' for j,s,k,t in paths))
    raw+=sum(map(len,groups.values()));unique+=len(groups)
(p/'fourth_mesh.txt').write_text('\n'.join(rows)+'\n')
print('PASS fourth geometry:',len(panels),'panels;',96*raw,'raw paths;',96*unique,'coalesced evaluations per sector',flush=True)
