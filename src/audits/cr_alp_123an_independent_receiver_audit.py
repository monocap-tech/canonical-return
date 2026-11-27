#!/usr/bin/env python3
"""Independent 123AN audit of the admitted AK/AL/AM receiver.
No physical coefficient generation or original interval acceptance claimed.
Replay frozen rational frames, re-enclose signed moments on a separate grid,
and independently derive the kernel and strict AM budgets.
"""
import hashlib,json,sys,zipfile
from fractions import Fraction as F
from math import isqrt,factorial
from pathlib import Path
import numpy as np
sys.set_int_max_str_digits(0)
G=10**60; J=2047; DC=F('1.822938456')
PINS=['a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473','a2c465dd692b0afbc0cd1a15b9f3599864085101bb8f52f7d225c1125adf453c','b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e','7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897','cd7b591bc60a57bc41cc9657962815caa639b921cc8085d4e561704fb7b29d4e']
def read(p,h,name='result.json'):
    if zipfile.is_zipfile(p):
        with zipfile.ZipFile(p) as z:
            assert z.testzip() is None;b=z.read(name)
    else:b=Path(p).read_bytes()
    assert hashlib.sha256(b).hexdigest()==h
    return json.loads(b)
def floor(v):return v.numerator//v.denominator
def ceil(v):return -floor(-v)
def rt(v,grid=10**110):
    assert v>=0;v=v*grid*grid
    return F(isqrt(floor(v))+1,grid)
def up(v,grid=10**9):return F(ceil(v*grid),grid)
def arr(x):return np.array([[int(y) for y in r] for r in x],dtype=object)
def pd(M):
    A=[[int(x) for x in r] for r in M];n=len(A);d=1
    assert A==list(map(list,zip(*A)))
    for k in range(n):
        p=A[k][k];assert p>0
        for i in range(k+1,n):
            for j in range(i,n):
                q,r=divmod(p*A[i][j]-A[i][k]*A[k][j],d);assert r==0
                A[j][i]=A[i][j]=q
        d=p
    return True
def box(a,b):
    assert a<=b;return floor(a*G),ceil(b*G)
def trans(rows,B,Q):
    # Endpoint arithmetic on separate grid; no AL midpoint implementation.
    lo=np.array([a for a,b in rows],dtype=object);hi=np.array([b for a,b in rows],dtype=object)
    bp=np.maximum(B,0);bn=np.minimum(B,0)
    return [ (F(int(a),G*Q),F(int(b),G*Q)) for a,b in zip(bp@lo+bn@hi,bp@hi+bn@lo)]
def ab(v):return max(abs(x) for x in v)
def mul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def scale(a,b):return mul(a,(b,b))
def atan(x):
    s=sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(50));r=x**101/101
    return s-r,s+r
def lnsmall(x):
    t=(x-1)/(x+1);assert 0<=t<=F(1,3)
    s=2*sum(t**(2*i+1)/(2*i+1) for i in range(50))
    return s,s+2*t**101/(101*(1-t*t))
def ln(n):
    k=n.bit_length()-1;return plus(scale(lnsmall(F(2)),k),lnsmall(F(n,2**k)))
def sinlower(t):
    l,h=box(*t);m=F(l+h,2*G);r=F(h-l,2*G);assert abs(m)<7
    s=sum((-1)**i*m**(2*i+1)/factorial(2*i+1) for i in range(40))
    v=abs(s)-abs(m)**81/factorial(81)-r
    ans=F(floor(v*10**15),10**15);assert ans>0;return ans
def kernel():
    nu=scale(plus(scale(atan(F(1,5)),16),scale(atan(F(1,239)),-4)),F(5,11));assert nu[0]>F(15,11)
    th=[mul(nu,ln(n)) for n in (2,3,4,5,7,8,9)]
    # Same analytic quadrature/Abel inequalities, independently evaluated.
    c=[[F(0) for _ in range(9)] for _ in range(9)];r=[x[:] for x in c]
    c[0][0]=c[0][1]=c[1][0]=F(1,2*J);c[1][1]=F(1,J)
    for i in range(2):
        for j in range(2):r[i][j]=F(3,J*J)
    for i,t in enumerate(th,2):
        c[i][i]=F(1,4*J);r[i][i]=(3+1/sinlower(scale(t,2)))/(J*J)
        for j in range(2):r[i][j]=r[j][i]=1/(J*J*sinlower(t))
        for j in range(2,i):
            u=th[j-2];r[i][j]=r[j][i]=(1/sinlower(plus(t,u))+1/sinlower(plus(t,neg(u))))/(2*J*J)
    return [[tuple(F(x,G) for x in box(*mul((c[i][j]-r[i][j],c[i][j]+r[i][j]),(1/nu[1]**2,1/nu[0]**2)))) for j in range(9)] for i in range(9)]
def main(ac,ak,al,z,am,out):
    A,K,L,Z,M=[read(p,h,'joint-repaired.json' if i==3 else 'result.json') for i,(p,h) in enumerate(zip((ac,ak,al,z,am),PINS))]
    assert all(set(d['sectors'])=={'1','2'} for d in (A,K,L,Z,M))
    assert F(11,14)+(1+F(250,609))/F(15*2049,11)==F(M['gamma_coefficient_upper'])<F(4,5)
    KG=kernel();res={}
    for k,a in A['sectors'].items():
        s=K['sectors'][k];l=L['sectors'][k];z=Z['sectors'][k][0];m=M['sectors'][k];n=44-int(k)
        E,U,b2=map(F,(a['E'],a['U'],a['finite_charge_upper']))
        assert E==F(z['joint_energy_lower']) and U==F(z['weighted_inverse_error_upper']) and E>U>0 and b2>=0
        assert z['cutoff']==850 and z['finite_prefix']==2048
        T=int(s['matrix_grid']);Q=int(l['whitening_grid']);H=arr(s['corrected_energy_lower']);P=arr(s['prefix_inverse_gram_upper']);B=arr(l['whitening_matrix']);eye=np.eye(n,dtype=object)
        p=F(l['prefix_relative_upper']);assert pd(H) and pd(B@H@B.T-T*Q*Q*eye) and pd(p.numerator*H-p.denominator*P)
        w=list(map(F,s['frame_roots']));pv=list(map(F,a['LDL_pivots']));assert all(x*x>=y>0 for x,y in zip(w,pv))
        tails=a['tail_receipts'];known=[F(0)]*n;series={};old={(x,i):F(v) for x,i,v in l['moment_norm_uppers']};count=0
        for name,mult,shift in [('centered_combined',F(20,11),-1),('gamma_sine',F(80,11),0),('prime_sine',F(80,11),0)]:
            series[name]=F(0)
            for i in range(32):
                power=2*i+2-shift;fac=mult*rt(F(J,2*(2*power-1)))/F(15*J,11)**power
                rows=[]
                for j,(wi,t) in enumerate(zip(w,tails)):
                    x,y=map(F,t['moments'][name][i]);known[j]+=F(floor(ab((x,y))*fac*10**100),10**100)
                    rows.append(box(x*wi*fac,y*wi*fac))
                v=trans(rows,B,Q);bound=rt(sum(ab(x)**2 for x in v));assert bound<=old[name,i]
                series[name]+=old[name,i];count+=1
        rem=[]
        for j,t in enumerate(tails):
            r=F(t['sine_series_norm_upper'])-known[j];assert r>=0
            rem.append(ceil(w[j]*(r+F(t['pole']['pole_remainder_norm_upper']))*G))
        rembound=rt(sum(F(int(x),G*Q)**2 for x in abs(B)@np.array(rem,dtype=object)))
        assert rembound<=F(l['remainder_norm_upper'])
        coeff=[]
        for i in range(9):
            rows=[]
            for wi,t in zip(w,tails):
                pol=t['pole'];x,y=map(F,[pol['constant'],pol['log_coefficient'],*pol['cosine_coefficients']][i]);rows.append(box(x*wi,y*wi))
            coeff.append(trans(rows,B,Q))
        # Independent interval trace of the 9-function Gram.
        lead=F(0)
        for j in range(n):
            for i in range(9):
                for h in range(9):lead+=mul(mul(coeff[i][j],KG[i][h]),coeff[h][j])[1]
        lead=up(lead);assert lead>=0
        # The slightly different kernel rounding is allowed: independently
        # recomputed lead bound, not a forced reproduction of AL's bound.
        tail=(rt(lead)+series['centered_combined']+series['gamma_sine']/5+series['prime_sine']+F(l['remainder_norm_upper']))**2/DC
        rho=up(p+tail,10**6);assert rho<1
        eta=rho/(1-rho);loss=up((U+(rt(b2)+rt(eta*U))**2)/E);assert loss<1
        assert F(m['pole_schur_over_E_lower'])==1-F(m['loss_over_E_upper'])>0
        res[k]={'dimension':n,'positive_matrix_tests':3,'moment_bounds_rechecked':count,'remainder_upper_rechecked':True,'independent_leading_trace_upper':str(lead),'independent_rho_upper':str(rho),'independent_loss_over_E_upper':str(loss),'independent_margin_lower':str(1-loss),'committed_loss_over_E_upper':m['loss_over_E_upper'],'classification':'CONDITIONAL_RECEIVER_AUDIT_PASS'}
        print(k,res[k],flush=True)
    Path(out).write_text(json.dumps({'status':'CONDITIONAL_RECEIVER_AUDIT_PASS','sectors':res,'scope':'AK matrix premises and original physical/interval enclosures inherited; no physical regeneration or native certification.'},indent=2)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
