#!/usr/bin/env python3
"""123AJ: exact interval test of the complete real scalar ray of each 123AI candidate.
No new solves. Inherited physical comparator and frozen prefix/sparse enclosures.
"""
import hashlib,json,sys,zipfile
from fractions import Fraction as F
from math import factorial,isqrt
from pathlib import Path
sys.set_int_max_str_digits(0)
S=10**80
class I:
    __slots__=('l','h')
    def __init__(self,l,h=None,raw=False):
        if raw:self.l,self.h=l,h
        else:
            l=F(l);h=l if h is None else F(h)
            self.l=l.numerator*S//l.denominator;self.h=-((-h.numerator*S)//h.denominator)
        assert self.l<=self.h
    def __add__(x,y):
        y=iv(y);return I(x.l+y.l,x.h+y.h,True)
    __radd__=__add__
    def __neg__(x):return I(-x.h,-x.l,True)
    def __sub__(x,y):return x+-iv(y)
    def __rsub__(x,y):return iv(y)+-x
    def __mul__(x,y):
        y=iv(y);z=[a*b for a in (x.l,x.h) for b in (y.l,y.h)]
        return I(min(z)//S,-((-max(z))//S),True)
    __rmul__=__mul__
    def __truediv__(x,y):
        y=iv(y)
        if y.h<0:return -x/-y
        assert y.l>0
        z=[(a*S,b) for a in (x.l,x.h) for b in (y.l,y.h)]
        return I(min(a//b for a,b in z),max(-((-a)//b) for a,b in z),True)
    def __rtruediv__(x,y):return iv(y)/x
    def sq(x):return I((0 if x.l<=0<=x.h else min(x.l*x.l,x.h*x.h))//S,-((-max(x.l*x.l,x.h*x.h))//S),True)
    def sqrt(x):assert x.l>=0;return I(isqrt(x.l*S),isqrt(x.h*S)+1,True)
    def mid(x):return F(x.l+x.h,2*S)
    def abslo(x):return F(max(0,x.l,-x.h),S)
    def up(x):return F(max(abs(x.l),abs(x.h)),S)
    def widen(x):return I(x.l-1,x.h+1,True)
    def pair(x):return [str(F(x.l,S)),str(F(x.h,S))]
def iv(v):return v if isinstance(v,I) else I(v)
ZERO=I(0);ONE=I(1);A=F(11,10);CAP=F('4.447938456');M=CAP+F(21,8)
assert F(4)**120/factorial(120)<F(1,S)
assert F(1,3)**241*2/(241*(1-F(1,9)))<F(1,S)
assert F(2)**120/factorial(120)/(1-F(2,121))<F(1,S)
def atan(x):
    assert 0<x<=F(1,3)
    s=sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(120));err=x**241/241
    return I(s-err,s+err)
PI=16*atan(F(1,5))-4*atan(F(1,239));NU=PI/(2*A)
def logseries(x):
    u=(x-1)/(x+1);assert u.up()<=F(1,3)+F(1,S)
    t=u;z=ZERO;v=u.sq()
    for k in range(120):z+=t/(2*k+1);t=t*v
    return (2*z).widen()
LOG2=logseries(I(2))
def logn(n):
    k=n.bit_length()-1
    return k*LOG2+logseries(I(F(n,2**k)))
def exp(x):
    assert 0<=x.l and x.up()<2
    t=ONE;z=ONE
    for k in range(1,120):t=t*x/k;z+=t
    return z.widen()
def sincos(x):
    n=round(x.mid()/(2*PI).mid());r=x-2*n*PI;assert r.up()<4
    q=-r.sq();s=r;c=ONE;ts=r;tc=ONE
    for k in range(1,60):
        ts=ts*q/((2*k)*(2*k+1));tc=tc*q/((2*k-1)*(2*k));s+=ts;c+=tc
    return s.widen(),c.widen()
def readzip(path,expected):
    with zipfile.ZipFile(path) as z:b=z.read('result.json')
    assert hashlib.sha256(b).hexdigest()==expected
    return json.loads(b)
def roots(v):
    assert v>=0;D=10**110;q=v*D*D;t=isqrt(q.numerator//q.denominator)
    return F(t,D),F(t+1,D)
def ratio(v,E,up=False):
    x=v/E*10**12;return str(F(-((-x).__floor__()) if up else x.__floor__(),10**12))
def main(ac,ag,ai,bands,out):
    ac=readzip(ac,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    ag=readzip(ag,'45719f6cd161b29534983c95f52941d434832bfffd3493d4d421513b436c1abf')
    ai=readzip(ai,'eac31a8bb64b36570021aa8bf51089d6e79dd2ca07f4f997562395720963cddb')
    raw=Path(bands).read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='72d7b5a04418dba9ac607fac28205f6c4d4e614ca0c8620d586bfb6a426a54c9'
    bands=json.loads(raw);assert bands[0][0]==86 and bands[-1]==[2049,None,4447938456]
    assert all(b[1]+1==c[0] for b,c in zip(bands,bands[1:]))
    def diag(j):return min(CAP,F(next(v for l,h,v in bands if l<=j and (h is None or j<=h)),10**9))
    ns=(2,3,4,5,7,8,9);ps=(2,3,2,5,7,2,3);TS=[logn(n) for n in ns];W=[logn(p)/I(n).sqrt() for n,p in zip(ns,ps)]
    eh=exp(I(A/2));em=1/eh;ee=exp(I(A));q=1/ee;rt=I(A).sqrt();qt=F(1,4)
    result={}
    for k in ('1','2'):
        a=ac['sectors'][k];g=ag['sectors'][k];r=ai['sectors'][k];kk=int(k)
        E,U,b2,eta=map(F,(a['E'],a['U'],a['finite_charge_upper'],a['eta_upper']));assert E>U>0 and eta>0 and b2>=0
        low=list(range(kk,86,2));hi=[j for j,_ in r['candidate_coordinates']];samples=[j for j,_,_ in r.get('sparse_defect_intervals',[])]
        assert hi==[j for j,_,_ in g['signed_prefix_intervals']]==[j for j,_ in a['prefix_mid_radius']]
        wg=int(r['candidate_grid']);wi={j:int(v) for j,v in r['candidate_coordinates']};w={j:I(F(v,wg)) for j,v in wi.items()}
        ids=low+hi+samples;nu={j:j*NU for j in ids};tr={j:[sincos(nu[j]*t) for t in TS] for j in ids}
        ap={j:sum((v*sc[0] for v,sc in zip(W,tr[j])),ZERO) for j in ids}
        b0=(eh+em)/2 if kk==1 else -(eh-em)/2;b1=-(eh-em)/4 if kk==1 else (eh+em)/4
        ep=(b0+2*b1)/2;en=(b0-2*b1)/2;B=(ee-q)/2+(A if kk==1 else -A)
        c={j:2*b0*nu[j]/(rt*(nu[j].sq()+qt)) for j in ids}
        st=B-sum((c[j].sq() for j in low),ZERO);assert st.l>0
        weighted=CAP*st+sum(((diag(j)-CAP)*c[j].sq() for j in hi),ZERO)
        def pss(i,j):
            if i==j:return -sum((v*((2*A-t)/A*co+si/(A*nu[i])) for v,t,(si,co) in zip(W,TS,tr[i])),ZERO)
            if j<i:return pss(j,i)
            return -4/PI*(j*ap[i]-i*ap[j])/(j*j-i*i)
        Pb={j:-sum((2*v*(nu[j]*(ep*I(n).sqrt()+en/I(n).sqrt()+b0*co)-b1*si)/(rt*(nu[j].sq()+qt)) for n,v,(si,co) in zip(ns,W,tr[j])),ZERO) for j in ids}
        Pbb=ZERO
        for n,t,v in zip(ns,TS,W):
            for pp,cp in [(F(1,2),ep),(F(-1,2),en)]:
                for rr,cr in [(F(1,2),ep),(F(-1,2),en)]:
                    vv=pp+rr;integ=2*A-t if vv==0 else ee.sq()/n-1 if vv==1 else 1-n*q.sq()
                    Pbb-=2*v*cp*cr*(I(n).sqrt() if rr>0 else 1/I(n).sqrt())*integ
        ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in ids}
        phh=Pbb-sum((c[j]*(Pb[j]+ps[j]) for j in low),ZERO)
        tau=sum((w[j]*c[j] for j in hi),ZERO)/st
        wp=sum((w[j]*ps[j] for j in hi),ZERO);wd=sum((w[j]*diag(j)*c[j] for j in hi),ZERO)
        # Exact integer interval accumulation of off-diagonal prime quadratic.
        # P_ij=-(4/pi)(j AP_i-i AP_j)/(j^2-i^2).  Keep coefficient integers
        # outside the interval division, then divide once by the rational grid.
        ql=qh=0
        for pos,i in enumerate(hi):
            ai0=ap[i];wi0=wi[i]
            for j in hi[pos+1:]:
                aj=ap[j];den=j*j-i*i;lo=(j*ai0.l-i*aj.h)//den;up=-((-(j*ai0.h-i*aj.l))//den)
                m=2*wi0*wi[j];ql+=min(m*lo,m*up);qh+=max(m*lo,m*up)
        quad=-4/PI*I(F(ql,S*wg*wg),F(qh,S*wg*wg))+sum((w[j].sq()*pss(j,j) for j in hi),ZERO)
        energy=quad+sum((w[j].sq()*diag(j) for j in hi),ZERO)-2*tau*(wp+wd)+tau.sq()*(phh+weighted)
        assert energy.l>0
        eg=int(g['prefix_grid']);pair=ZERO
        for j,l,h in g['signed_prefix_intervals']:pair+=w[j]*I(F(int(l),eg),F(int(h),eg))
        assert pair.l>0
        # Recompute weighted pairing lower from exact stored AI intervals.
        pg=int(r['pair_grid']);piv=list(map(F,a['LDL_pivots']));assert len(piv)==len(r['signed_pairing_intervals'])==44-kk
        zlo=sum(p*F(max(0,int(l),-int(h))**2,pg*pg) for p,(l,h) in zip(piv,r['signed_pairing_intervals']))
        assert zlo>0
        # For lambda, V>=|<w0,e>-lambda <w0,Cw0>|^2/<w0,Cw0>.
        # Also use original sparse defect d1 and Cw0, for d_lambda=d1+(1-lambda)Cw0.
        cw=[];d1=[];sg=int(r.get('sparse_grid','1'))
        for j,l,h in r.get('sparse_defect_intervals',[]):
            cp=c[j]/st
            cv=sum((w[i]*pss(i,j) for i in hi),ZERO)-cp*wp-tau*ps[j]+tau*cp*phh-cp*wd-tau*CAP*c[j]+tau*cp*weighted
            cw.append(cv);d1.append(I(F(int(l),sg),F(int(h),sg)))
        zroot=roots(zlo)[0];broot=roots(b2)[1];budget=roots(E-U)[1]
        # Rational interval cover, no floating decisions. Outside [-L,L],
        # the finite term alone exceeds the budget. For each cell, lower
        # bound both charges uniformly. No sign assumed for lambda.
        L=F(1,2) if kk==1 else F(5,4)
        assert L*zroot-broot>budget
        # Adaptive interval cover: subdivide only cells not yet rejected.
        stack=[(-L,L,0)];failed=[];minratio=None;N=0
        while stack:
            l,h,depth=stack.pop();lam=I(l,h)
            finite=max(F(0),lam.abslo()*zroot-broot)
            dual=(pair-lam*energy).abslo()**2/F(energy.h,S)
            sp=sum((x+(1-lam)*y).abslo()**2 for x,y in zip(d1,cw))/M
            charge=finite+roots(eta*max(dual,sp))[0]
            if charge<=budget and depth<12:
                m=(l+h)/2;stack.extend([(l,m,depth+1),(m,h,depth+1)]);continue
            N+=1;q0=charge*charge/E
            if minratio is None or q0<minratio:minratio=q0
            if charge<=budget:failed.append([str(l),str(h)])
        # Include the complementary unbounded real intervals via the finite term.
        minratio=min(minratio,(L*zroot-broot)**2/E)
        rec={'status':'RAY_REJECTED' if not failed else 'RAY_TEST_INCONCLUSIVE','cells':N,'unrejected_cells':failed,'lambda_extent':str(L),'energy_interval':energy.pair(),'e_pair_interval':pair.pair(),'energy_over_E_lower':ratio(F(energy.l,S),E),'energy_over_E_upper':ratio(F(energy.h,S),E,True),'min_cell_charge2_over_E_lower':str(F((minratio*10**12).__floor__(),10**12)),'available_budget_over_E_upper':ratio(E-U,E,True),'sparse_Cw_intervals':[x.pair() for x in cw],'source_candidate':r['classification']}
        result[k]=rec;print(k,rec['status'],len(failed),'energy/E',float(energy.mid()/E),'min charge2/E',float(minratio),'budget',float((E-U)/E),flush=True)
        Path(out).write_text(json.dumps({'status':'CONDITIONAL_FROZEN_REAL_RAY_DIAGNOSTIC','sectors':result,'scope':'Fixed AC eta, original E,U,b2; rejection of sufficient receiver on real scalar ray only, not pole sign or actual inverse lower bound.'},indent=2)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
