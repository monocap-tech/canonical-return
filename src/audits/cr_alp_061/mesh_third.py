"""Exact mesh for the third signed-prime moment of the frozen physical proxy.
Run after mesh_second.py. Refine at every translated knot of every active
inner translation, on each half-overlap for the seven outer translations.
"""
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
moves=[(j,s,s*x) for j,x in enumerate(sh) for s in (-1,1)]
vk=sorted({F(0),d,*[k-z for k in knots for _,_,z in moves if 0<k-z<d]})
rows=[];counts=[];evals=[]
for t in sh:
    end=(d-t)/2
    cuts=sorted({F(0),end,*[k-z for k in vk for z in (F(0),t) if 0<k-z<end]})
    panels=[]
    def active(mid,shift,l,r):
        a=[(j,s,z,bisect_right(knots,mid+shift+z)-1) for j,s,z in moves if 0<mid+shift+z<d]
        assert all(knots[i]<=l+shift+z and r+shift+z<=knots[i+1] for j,s,z,i in a)
        assert all(r+shift+z<=0 or l+shift+z>=d for j,s,z in moves if not 0<mid+shift+z<d)
        return a
    def panel(l,r,left,right):
        mid=(l+r)/2;h=(r-l)/2;R=4*h
        args=[mid+z for _,_,z,_ in left]+[mid+t+z for _,_,z,_ in right]
        if h<=F(1,128) and all(min(y,d-y)>=2*R or cmax*R<=F(1,2) for y in args):
            panels.append((mid,h,left,right));return
        panel(l,mid,left,right);panel(mid,r,left,right)
    for l,r in zip(cuts,cuts[1:]):
        mid=(l+r)/2
        panel(l,r,active(mid,F(0),l,r),active(mid,t,l,r))
    assert sum(2*h for _,h,_,_ in panels)==end
    rows.append(str(len(panels)));counts.append(len(panels))
    evals.append(96*sum(len(l)+len(r) for _,_,l,r in panels))
    def fmt(a):return str(len(a))+' '+' '.join(f'{j} {s} {i}' for j,s,z,i in a)
    for mid,h,left,right in panels:rows.append(enc(mid)+' '+enc(h)+' '+fmt(left)+' '+fmt(right))
# v0 = P_rat r_proxy minus a bounded exact low-sine/pole combination.
# The native evaluator checks those coefficient bounds; its complex disk
# bound is <= 10^42 (8*10^40 from the translations plus projection).
assert 4*d*F(4,3)*10**84/4**192<F(1,10**12)
(p/'third_mesh.txt').write_text('\n'.join(rows)+'\n')
print('PASS third-moment geometry:',counts,'panels;',sum(evals),'shifted evaluations per sector',flush=True)
