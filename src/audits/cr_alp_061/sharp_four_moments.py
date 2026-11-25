"""Exact four-scalar-moment extremal inverse comparison; no new quadrature."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
class I:
    def __init__(self,x,y=None):
        self.l=F(x);self.h=F(x if y is None else y);assert self.l<=self.h
    @staticmethod
    def cast(x):return x if isinstance(x,I) else I(x)
    def __add__(self,x):
        x=I.cast(x);return I(self.l+x.l,self.h+x.h)
    __radd__=__add__
    def __neg__(self):return I(-self.h,-self.l)
    def __sub__(self,x):return self+-I.cast(x)
    def __rsub__(self,x):return I.cast(x)+-self
    def __mul__(self,x):
        x=I.cast(x);a=[a*b for a in (self.l,self.h) for b in (x.l,x.h)];return I(min(a),max(a))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=I.cast(x);assert x.l*x.h>0;return self*I(1/x.h,1/x.l)
    def sq(self):return I(0 if self.l<=0<=self.h else min(self.l**2,self.h**2),max(self.l**2,self.h**2))
    def pair(self):return [str(self.l),str(self.h)]
ell,u,g,s=map(F,('0.255','5.505','2.88','2.625'))
def bounds(G,M,N,T):
    h=s*s*G-N;c=g*h+s*s*M-T
    m=g*G+M;d=g*M+N;e=g*N+T
    assert h.l>0 and c.l>0 and m.l>0
    schur=e-d.sq()/m;assert schur.l>0
    upper=(g*G-M-h.sq()/c)/(ell*u)
    lower=G.sq()/m+(M-d*G/m).sq()/schur
    return lower,upper
def test(nodes,weights,which):
    z=[I(sum(w*(t-g)**k for t,w in zip(nodes,weights))) for k in range(4)]
    j=sum(w/t for t,w in zip(nodes,weights));lo,hi=bounds(*z)
    assert lo.l<=j<=hi.h
    assert (lo.l if which=='lo' else hi.h)==j
if __name__=='__main__':
    test([F(1),F(4)],[F(2),F(3)],'lo')
    test([ell,F(2),u],[F(1),F(2),F(3)],'hi')
    print('PASS exact Gauss and Lobatto extremizers')
    if len(sys.argv)==1:sys.exit(0)
    old,second,third=(json.loads(Path(p).read_text()) for p in sys.argv[1:4]);out={}
    for k in ('1','2'):
        G,M,K=[I(*old[k][x]) for x in ('true_G','true_M','K')]
        N=I(*second[k]['true_N']);T=I(*third[k]['true_T'])
        lo,hi=bounds(G,M,N,T)
        status='PASS_WITNESS' if hi.h<K.l else 'FAIL_COMPARISON_WITNESS' if lo.l>K.h else 'FOUR_SCALAR_MOMENTS_INSUFFICIENT' if lo.h<K.l and K.h<hi.l else 'UNDECIDED_INTERVAL_WIDTH'
        out[k]={'sharp_lower':lo.pair(),'sharp_upper':hi.pair(),'K':K.pair(),'status':status}
        print(k,status,'lower',float(lo.l),float(lo.h),'upper',float(hi.l),float(hi.h),'K',float(K.l),float(K.h))
    if len(sys.argv)>4:Path(sys.argv[4]).write_text(json.dumps(out,indent=2)+'\n')
