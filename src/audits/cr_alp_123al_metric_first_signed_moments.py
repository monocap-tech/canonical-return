#!/usr/bin/env python3
"""123AL: metric-first signed moment receiver, no response solves.
Inputs original AC/AK ZIPs; inherited comparator and asymptotic premises.
NumPy proposes whitening and shift only; exact rational tests accept.
"""
import hashlib,json,sys,zipfile
from fractions import Fraction as F
from math import isqrt,factorial
from pathlib import Path
import numpy as np
sys.set_int_max_str_digits(0)
T=10**40; Q=10**30; J0=2047; DC=F('1.822938456')
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def root(x):
    assert x>=0;G=10**110;z=x*G*G
    return F(isqrt(z.numerator//z.denominator)+1,G)
def read(p,h):
    with zipfile.ZipFile(p) as z:b=z.read('result.json')
    assert hashlib.sha256(b).hexdigest()==h
    return json.loads(b)
def positive(A):
    a=[[int(x) for x in r] for r in A];prev=1
    for k in range(len(a)):
        p=a[k][k]
        if p<=0:return False
        for i in range(k+1,len(a)):
            for j in range(i,len(a)):
                v,r=divmod(p*a[i][j]-a[i][k]*a[k][j],prev);assert r==0
                a[i][j]=a[j][i]=v
        prev=p
    return True
def fp(A):return np.array(A,float)/float(T)
def box(l,h):
    assert l<=h;m=round((l+h)*T/2)
    return m,ceil(max(F(m)-l*T,h*T-F(m)))
def absup(v):return max(abs(F(x)) for x in v)
def mat(rows):return tuple(np.array([[v[i] for v in row] for row in rows],dtype=object) for i in (0,1))
def mul(A,B):
    a,x=A;b,y=B
    m=a@b;r=np.abs(a)@y+x@np.abs(b)+x@y
    return m//T,(r+T-1)//T+1
def tr(A):return A[0].T,A[1].T
def ivadd(a,b):return a[0]+b[0],a[1]+b[1]
def ivscale(a,s):return min(a[0]*s,a[1]*s),max(a[0]*s,a[1]*s)
def ivmul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def logseries(x):
    u=(x-1)/(x+1);assert 0<=u<=F(1,3);n=120
    s=2*sum(u**(2*i+1)/(2*i+1) for i in range(n))
    return s,s+2*u**(2*n+1)/((2*n+1)*(1-u*u))
def logn(n):
    k=n.bit_length()-1;return ivadd(ivscale(logseries(F(2)),k),logseries(F(n,2**k)))
def atan(x):
    n=100;s=sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(n));r=x**(2*n+1)/(2*n+1);return s-r,s+r
def kernel():
    # Same nine-function entrywise 123AD/AH enclosure, before Loewner rounding.
    nu=ivscale(ivadd(ivscale(atan(F(1,5)),16),ivscale(atan(F(1,239)),-4)),F(5,11));assert nu[0]>F(15,11)
    th=[ivmul(nu,logn(n)) for n in (2,3,4,5,7,8,9)]
    def sinlo(v):
        # Round angle enclosure to keep Taylor rationals bounded.
        g=10**70;l=F(floor(v[0]*g),g);h=F(ceil(v[1]*g),g);m=(l+h)/2;r=(h-l)/2;assert abs(m)<7
        s=sum((-1)**i*m**(2*i+1)/factorial(2*i+1) for i in range(60))
        lo=F(floor((abs(s)-abs(m)**121/factorial(121)-r)*10**20),10**20);assert lo>0;return lo
    C=[[F(0) for j in range(9)] for i in range(9)];R=[[F(0) for j in range(9)] for i in range(9)]
    C[0][0]=C[0][1]=C[1][0]=F(1,2*J0);C[1][1]=F(1,J0)
    for i in range(2):
        for j in range(2):R[i][j]=F(3,J0**2)
    for i,t in enumerate(th,2):
        C[i][i]=F(1,4*J0);R[i][i]=(3+1/sinlo(ivscale(t,2)))/J0**2
        for j in range(2):R[i][j]=R[j][i]=1/(J0**2*sinlo(t))
        for j in range(2,i):
            u=th[j-2];R[i][j]=R[j][i]=(1/sinlo(ivadd(t,ivscale(u,-1)))+1/sinlo(ivadd(t,u)))/(2*J0**2)
    scale=(1/nu[1]**2,1/nu[0]**2)
    return mat([[box(*ivmul((C[i][j]-R[i][j],C[i][j]+R[i][j]),scale)) for j in range(9)] for i in range(9)])
def omat(A):return [[str(int(v)) for v in row] for row in A]
def upper(v,g=10**9):return str(F(ceil(v*g),g))
def main(ac,ak,out):
    A=read(ac,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    K=read(ak,'a2c465dd692b0afbc0cd1a15b9f3599864085101bb8f52f7d225c1125adf453c')
    assert A['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER' and K['status']=='CONDITIONAL_ACTUAL_SCHUR_METRIC_DIAGNOSTIC'
    assert set(A['sectors'])==set(K['sectors'])=={'1','2'}
    gram=kernel();outd={}
    for k,z in A['sectors'].items():
        s=K['sectors'][k];n=44-int(k);assert int(s['matrix_grid'])==T and s['dimension']==n
        H=np.array([[int(x) for x in row] for row in s['corrected_energy_lower']],dtype=object)
        P=np.array([[int(x) for x in row] for row in s['prefix_inverse_gram_upper']],dtype=object)
        assert positive(H)
        B=np.linalg.inv(np.linalg.cholesky(fp(H)))
        bi=np.array([[round(float(x)*1.000001*Q) for x in row] for row in B],dtype=object)
        eye=np.array([[int(i==j) for j in range(n)] for i in range(n)],dtype=object)
        assert positive(bi@H@bi.T-T*Q*Q*eye)
        # Verify directly in the original frame: P <= p H.
        approx=np.linalg.eigvalsh(B@fp(P)@B.T)[-1]
        p=F(ceil(F(float(approx))*10**6)+1,10**6)
        assert positive(p.numerator*H-p.denominator*P)
        w=list(map(F,s['frame_roots']));piv=list(map(F,z['LDL_pivots']))
        assert len(w)==len(piv)==n and all(x*x>=v>0 for x,v in zip(w,piv))
        ts=z['tail_receipts'];assert len(ts)==n
        Bbox=mat([[box(F(int(v),Q),F(int(v),Q)) for v in row] for row in bi])
        def transform(rows):return mul(mat(rows),tr(Bbox))
        cc=[]
        for l in range(9):
            row=[]
            for wi,t in zip(w,ts):
                pol=t['pole'];v=[pol['constant'],pol['log_coefficient'],*pol['cosine_coefficients']][l]
                row.append(box(F(v[0])*wi,F(v[1])*wi))
            cc.append(row)
        lead=transform(cc);ll=mul(mul(tr(lead),gram),lead)
        # PSD leading tail Gram upper: trace, with outward diagonal bounds.
        lead2=F(sum(int(ll[0][i,i]+ll[1][i,i]) for i in range(n)),T);assert lead2>=0
        known=[F(0)]*n;series=F(0);moment_receipts=[]
        # Each scalar kernel has L2 norm <= the exact 123AC factor below.
        # Combine signed moments in the new frame before taking absolute values.
        for label,mult,shift in [('centered_combined',F(20,11),-1),('gamma_sine',F(80,11),0),('prime_sine',F(80,11),0)]:
            assert all(len(t['moments'][label])==32 for t in ts)
            for l in range(32):
                ss=2*l+2-shift
                fac=mult*root(F(J0,2*(2*ss-1)))/F(15*J0,11)**ss
                row=[]
                for i,(wi,t) in enumerate(zip(w,ts)):
                    v=t['moments'][label][l];known[i]+=F(floor(absup(v)*fac*10**100),10**100)
                    row.append(box(F(v[0])*wi*fac,F(v[1])*wi*fac))
                m,r=transform([row]);bound=root(F(sum((abs(int(x))+int(y))**2 for x,y in zip(m[0],r[0])),T*T))
                series+=bound;moment_receipts.append([label,l,upper(bound)])
        # 123AC sum of positive charges: subtract downward rounded known charges,
        # leaving an upper bound on the omitted geometric remainder, not on actual norm.
        residual=[]
        for i,t in enumerate(ts):
            rem=F(t['sine_series_norm_upper'])-known[i];assert rem>=0
            residual.append(w[i]*(rem+F(t['pole']['pole_remainder_norm_upper'])))
        # Unknown vector-valued remainders: column triangle inequality after whitening.
        # Round input upward; use absolute B on every output coordinate.
        ru=np.array([ceil(x*T) for x in residual],dtype=object)
        rr=np.abs(bi)@ru
        rem=root(F(sum(int(x)**2 for x in rr),Q*Q*T*T))
        tail=(root(lead2)+series+rem)**2/DC
        rho=F(ceil((p+tail)*10**6),10**6)
        E,U,b2=map(F,(z['E'],z['U'],z['finite_charge_upper']));assert E>U>0 and b2>=0
        rec={'dimension':n,'whitening_grid':str(Q),'whitening_matrix':omat(bi),'prefix_relative_upper':str(p),'leading_tail_trace_upper':upper(lead2),'series_norm_upper':upper(series),'remainder_norm_upper':upper(rem),'tail_relative_upper':upper(tail),'total_relative_upper':str(rho),'moment_norm_uppers':moment_receipts}
        if rho<1:
            eta=rho/(1-rho)
            # S_C >= (1-rho) H. True inverse norm of b independently <= b2
            # by the original Dbar comparison; no need to use new metric for b.
            loss=U+(root(b2)+root(eta*U))**2
            rec.update({'schur_fraction_lower':str(1-rho),'eta_upper':str(eta),'loss_over_E_upper':upper(loss/E),'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if loss<E else 'IMPROVED_METRIC_SUFFICIENT_BOUND_INCONCLUSIVE'})
        else:rec['classification']='METRIC_RELATIVE_BOUND_INCONCLUSIVE'
        outd[k]=rec;print(k,'p',float(p),'lead',float(root(lead2)),'series',float(series),'rem',float(rem),'rho',float(rho),'eta',rec.get('eta_upper'),'loss/E',rec.get('loss_over_E_upper'),flush=True)
    Path(out).write_text(json.dumps({'status':'CONDITIONAL_METRIC_FIRST_SIGNED_MOMENT_RECEIVER','sectors':outd,'scope':'Frozen capped comparator, AC response bank, AK corrected energy; inherited physical and signed-tail premises.'},indent=2)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
