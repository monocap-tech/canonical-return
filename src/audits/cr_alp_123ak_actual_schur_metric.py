#!/usr/bin/env python3
"""123AK: reconstruct corrected Schur energy in the frozen AC frame.
No response solves. Exact integer interval matrix products and Loewner tests.
Floating generalized eigenvectors select diagnostic directions only.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil
import hashlib,json,sys
import numpy as np
from cr_alp_123ah_signed_tail_pairing import kernel
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent
T=10**40
CAP=F('4.447938456'); M=CAP+F(21,8); DC=CAP-F(21,8)
def read(p,h):
    b=Path(p).read_bytes();assert hashlib.sha256(b).hexdigest()==h
    return json.loads(b)
def blob(p,h):
    b=p.read_bytes();assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==h
    return b
def floor(x):return x.numerator//x.denominator
def up(x):return -floor(-x)
def root(v,G=10**110):
    assert v>=0;z=v*G*G
    return F(isqrt(z.numerator//z.denominator)+1,G)
def box(l,h):
    l,h=F(l),F(h);assert l<=h
    m=round((l+h)*T/2);return m,up(max(F(m)-l*T,h*T-F(m)))
def mat(rows):
    return tuple(np.array([[v[i] for v in row] for row in rows],dtype=object) for i in (0,1))
def tr(A):return A[0].T,A[1].T
def add(A,B,sgn=1):return A[0]+sgn*B[0],A[1]+B[1]
def mul(A,B):
    a,x=A;b,y=B
    m=a@b;err=np.abs(a)@y+x@np.abs(b)+x@y
    return m//T,(err+T-1)//T+1
def sym_bounds(A):
    a,r=A;n=len(a);mid=np.zeros((n,n),dtype=object);err=mid.copy()
    for i in range(n):
        for j in range(i,n):
            l=max(a[i,j]-r[i,j],a[j,i]-r[j,i]);h=min(a[i,j]+r[i,j],a[j,i]+r[j,i])
            assert l<=h
            m=(l+h)//2;v=max(m-l,h-m)
            mid[i,j]=mid[j,i]=m;err[i,j]=err[j,i]=v
    reserve=np.diag(err.sum(axis=1))
    return mid-reserve,mid+reserve
def positive(A):
    a=[[int(x) for x in row] for row in A];n=len(a);prev=1
    for k in range(n):
        p=a[k][k]
        if p<=0:return False
        for i in range(k+1,n):
            for j in range(i,n):
                x,r=divmod(p*a[i][j]-a[i][k]*a[k][j],prev)
                assert r==0;a[i][j]=a[j][i]=x
        prev=p
    return True
def fp(A):return np.array(A,dtype=float)/float(T)
def quad(A,v):return int(v@A@v)
def outmat(A):return [[str(int(x)) for x in row] for row in A]
def generalized(A,B):
    # Untrusted direction or upper-shift proposal.
    L=np.linalg.cholesky(fp(B));Y=np.linalg.solve(L,fp(A));Y=np.linalg.solve(L,Y.T).T
    eig,U=np.linalg.eigh((Y+Y.T)/2)
    return float(eig[-1]),np.linalg.solve(L.T,U[:,-1])
def namespace():
    raw=blob(ROOT/'cr_alp_046_full_finite.py','554e7ab54da4545987d89aa660a4abe6c951cb63')
    s=raw.decode()
    def replace(a,b):
        nonlocal s
        assert s.count(a)==1,a;s=s.replace(a,b)
    replace('range(1,N+1)};trig=', 'range(1,851)};trig=')
    replace('for j,v in nu.items():\n re,im,tre=','for j,v in nu.items():\n if j>85:continue\n re,im,tre=')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    replace(old,old+', Bmass, Hbb, c, Hb, Pbb')
    g={'__file__':str(ROOT/'cr_alp_046_full_finite.py'),'__name__':'cr_alp_123ak_source'}
    exec(compile(s,str(ROOT/'cr_alp_046_full_finite.py'),'exec'),g)
    return g
def main(archive,ac,out):
    archive=Path(archive)
    A=read(ac,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    frozen=read(archive/'cr_alp_046/CR-ALP-046_full_E85_intervals.json','d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22')
    bands=json.loads(blob(ROOT/'CR-ALP-068_pole_centered_gamma_result.json','1a7fe8629daa8b8e6cec185a79b6746d6b28975d'))['bands']
    assert A['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER' and set(A['sectors'])=={'1','2'}
    def dg(b):return F(floor(F(b['gamma_lower'])*10**9),10**9)
    assert next(dg(b) for b in bands if b['start']==2049)==CAP and all(dg(b)>=CAP for b in bands if b['start']>=2049)
    def diag(j):return min(CAP,dg(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end']))))
    g=namespace();I,S,ZERO=(g[x] for x in ('I','S','ZERO'))
    a,L,nu,TS,W,trig,rt,eh,em,qt=(g[x] for x in ('A','L','nu','TS','W','trig','rt','eh','em','qt'))
    AP={j:sum((w*t[0] for w,t in zip(W,trig[j])),ZERO) for j in nu}
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*ct+st/(a*u)) for w,t,(st,ct) in zip(W,TS,trig[i])),ZERO)
        return -2/a*(v*AP[i]-u*AP[j])/(v.sq()-u.sq())
    def ib(v):return box(F(v.l,S),F(v.h,S))
    ans={}
    for bk in frozen['blocks']:
        k=bk['parity'];ks=str(k);z=A['sectors'][ks];n=44-k
        order,H,Mass,B,Hbb,c0,Hb,Pbb=g['build'](k)
        low=bk['ids'];assert len(low)==n and set(low)==set(order)
        ix=[order.index(j) for j in low]
        fs=int(frozen['scale']);assert fs==S
        for mat0,key in [(H,'H'),(Mass,'M')]:
            for i,u in enumerate(ix):
                for j,v in enumerate(ix):
                    l,h=bk[key][i][j];assert max(mat0[u][v].l,int(l))<=min(mat0[u][v].h,int(h))
        perm=z['coordinate_order'];assert sorted(perm)==list(range(n))
        L0=[[F(v) for v in row] for row in z['LDL_lower']];p=list(map(F,z['LDL_pivots']))
        D=[[F(v) for v in row] for row in z['Dbar']]
        assert min(p)>0 and all(sum(L0[i][l]*p[l]*L0[j][l] for l in range(n))==D[i][j] for i in range(n) for j in range(n))
        w=[root(v,10**30) for v in p];ratio=max(v*v/t for v,t in zip(w,p))
        assert ratio>=1
        vb=read(archive/'cr_alp_049'/f'result_{k}.json',{1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k])
        vs=int(vb['scale']);VV=[[I(F(int(l),vs),F(int(h),vs)) for l,h in row] for row in vb['V']]
        Vt=mat([[ib(sum((VV[i][perm[l]]*I(L0[l][j]*w[j]) for l in range(n)),ZERO)) for j in range(n)] for i in range(n)])
        Ht=mul(mul(tr(Vt),mat([[ib(H[i][j]) for j in ix] for i in ix])),Vt)
        bank=z['bank_indices'];assert bank==list(range(87 if k==1 else 86,851,2))
        high=[j for j,row in z['prefix_mid_radius']]
        assert high==list(range(87 if k==1 else 86,2049,2))
        qg=int(z['response_coordinate_grid']);pg=int(z['prefix_grid']);assert qg==10**40 and pg==10**60
        Q=mat([[box(F(v,qg)*w[i],F(v,qg)*w[i]) for i,v in enumerate(row)] for row in z['response_transformed_coordinates']])
        zr=dict(z['prefix_mid_radius'])
        def zmat(js,weighted=False):
            return mat([[box(F(m-r,pg)*w[i]/(diag(j)-F(21,8) if weighted else 1),F(m+r,pg)*w[i]/(diag(j)-F(21,8) if weighted else 1)) for i,(m,r) in enumerate(zr[j])] for j in js])
        Z=zmat(high);ZB=zmat(bank)
        Glo,_=sym_bounds(mul(tr(Z),mat([[box(F(m-r,pg)*w[i]/M,F(m+r,pg)*w[i]/M) for i,(m,r) in enumerate(zr[j])] for j in high])))
        _,Gp=sym_bounds(mul(tr(Z),zmat(high,True)))
        b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4
        ep=(b0+2*b1)/2;en=(b0-2*b1)/2
        def cj(j):
            v=j*g['NUF'];return 2*b0*v/(rt*(v.sq()+qt))
        c={j:cj(j) for j in low+high}
        def pb(j):
            v=nu[j];val=ZERO
            for nn,ww,(si,co) in zip(g['NS'],W,trig[j]):
                sn=I(nn).sqrt();bt=ep*sn+en/sn
                val-=2*ww/(rt*(v.sq()+qt))*(v*(bt+b0*co)-b1*si)
            return val
        Pb={j:pb(j) for j in low+bank}
        st=B-sum((c[j].sq() for j in low),ZERO);assert st.l>0
        weighted=I(CAP)*st+sum((I(diag(j)-CAP)*c[j].sq() for j in high),ZERO)
        ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in low+bank}
        phh=Pbb-sum((c[i]*Pb[i] for i in low),ZERO)-sum((c[i]*ps[i] for i in low),ZERO)
        cp={j:c[j]/st for j in bank}
        def ce(i,j):
            return (pss(i,j)-cp[i]*ps[j]-cp[j]*ps[i]+cp[i]*cp[j]*phh
                    +(I(diag(i)) if i==j else ZERO)-cp[i]*I(diag(j))*c[j]
                    -cp[j]*I(diag(i))*c[i]+cp[i]*cp[j]*weighted)
        m=len(bank);cb=[[None]*m for _ in bank]
        for i,j in enumerate(bank):
            for l in range(i,m):cb[i][l]=cb[l][i]=ib(ce(j,bank[l]))
        CB=mat(cb)
        CQ=mul(CB,Q);QZ=mul(tr(Q),ZB)
        Hhat=add(add(add(Ht,mul(tr(Q),CQ),-1),QZ,-1),tr(QZ),-1)
        Hlo,Hup=sym_bounds(Hhat)
        print(ks,'corrected energy reconstructed',flush=True)
        # C^{-1} >= I/M, so Gactual >= Glo; S_C <= Hup-Glo.
        # Generalized ratio lower from this upper metric is valid for ANY
        # receiving inverse metric dominating the true S_C^{-1}.
        Sup=Hup-Glo
        assert positive(Sup)
        candidates=[np.array([int(i==j) for i in range(n)],dtype=object) for j in range(n)]
        try:
            _,v=generalized(Glo,Sup)
            v=v/max(abs(v));candidates.append(np.array([round(float(x)*10**9) for x in v],dtype=object))
        except np.linalg.LinAlgError:pass
        rays=[(F(quad(Glo,v),quad(Sup,v)),v) for v in candidates if quad(Sup,v)>0]
        lower,v=max(rays,key=lambda x:x[0]);lower=max(F(0),lower)
        E,U,b2=map(F,(z['E'],z['U'],z['finite_charge_upper']));assert E>U>0 and b2>=0
        ceiling=E/U-1
        rec={'dimension':n,'matrix_grid':str(T),'frame_roots':[str(x) for x in w],
             'corrected_energy_lower':outmat(Hlo),'corrected_energy_upper':outmat(Hup),
             'defect_gram_lower':outmat(Glo),'prefix_inverse_gram_upper':outmat(Gp),
             'actual_metric_ray':[str(x) for x in v],'actual_metric_eta_lower':str(lower),
             'necessary_eta_ceiling':str(ceiling),'exact_metric_scalar_obstructed':lower>=ceiling}
        # Construct several valid full Gram upper matrices from the same bank.
        # For column norm bounds a_i, Gram <= diag(a_i sum a_j).
        trc=z['tail_receipts'];assert len(trc)==n
        tn=[w[i]*F(t['tail_norm_upper']) for i,t in enumerate(trc)]
        rn=[w[i]*(F(t['pole']['pole_remainder_norm_upper'])+F(t['sine_series_norm_upper'])) for i,t in enumerate(trc)]
        def diagmajor(a):return np.diag(np.array([up(T*x*sum(a)/DC) for x in a],dtype=object))
        G0=Gp+diagmajor(tn)
        kk=mat([[box(l,h) for l,h in row] for row in kernel()])
        coeff=mat([[box(F(t['pole'][key][0])*w[i],F(t['pole'][key][1])*w[i]) if key in ('constant','log_coefficient') else box(F(t['pole']['cosine_coefficients'][key][0])*w[i],F(t['pole']['cosine_coefficients'][key][1])*w[i]) for i,t in enumerate(trc)] for key in ['constant','log_coefficient',*range(7)]])
        _,Lead=sym_bounds(mul(mul(tr(coeff),kk),coeff))
        # Divide by dc, rounding upper matrix by row reserve.
        def scaleupper(B,q):
            nn=len(B);return np.array([[round(F(int(B[i,j]))*q)+(nn if i==j else 0) for j in range(nn)] for i in range(nn)],dtype=object)
        Lead=scaleupper(Lead,1/DC);Rem=diagmajor(rn)
        uppers=[('column_norm',G0)]+[(str(t),Gp+scaleupper(Lead,1+t)+scaleupper(Rem,1+1/t)) for t in (F(1,4),F(1),F(4))]
        # S_C >= I in this weighted frame by Dbar >= S_C^{-1} and w_i^2>=p_i.
        # S_C >= Hlo-Gup as well. Convex mixtures retain Loewner order.
        eye=np.eye(n,dtype=object);proposals=[]
        for name,Gu in uppers:
            for q in [1,2,4,8,16,32,64,128,256,512,1024]:
                Sn=(q-1)*T*eye+Hlo-Gu
                try:
                    if np.linalg.eigvalsh(fp(Sn)/q)[0]<=.500001:continue
                    lam,_=generalized(q*Gu,Sn)
                    if lam>=0 and np.isfinite(lam):proposals.append((lam,name,q,Gu,Sn))
                except np.linalg.LinAlgError:continue
        proposals.sort(key=lambda x:x[0])
        rec['metric_candidates_proposed']=len(proposals)
        accepted=None
        for lam,name,q,Gu,Sn in proposals[:3]:
            if not positive(2*Sn-q*T*eye):continue
            # Exact rational upper shift (floating eigenvalue only a proposal).
            alpha=F(ceil((lam+max(1e-6,abs(lam)*1e-6))*10**6),10**6)
            if not positive(alpha.numerator*Sn-alpha.denominator*q*Gu):continue
            loss=U+(root(2*ratio*b2)+root(alpha*U))**2
            accepted={'gram_majorant':name,'mix_denominator':q,'gram_upper':outmat(Gu),'schur_lower_numerator':outmat(Sn),'eta_upper':str(alpha),'finite_charge_upper':str(2*ratio*b2),'loss_upper':str(loss),'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE'}
            break
        rec['accepted_metric']=accepted
        ans[ks]=rec
        print(ks,'true Schur metric eta lower',float(lower),'threshold',float(ceiling),'obstructed',rec['exact_metric_scalar_obstructed'],'new metric eta',float(F(accepted['eta_upper'])) if accepted else None,flush=True)
        Path(out).write_text(json.dumps({'status':'CONDITIONAL_ACTUAL_SCHUR_METRIC_DIAGNOSTIC','sectors':ans,'scope':'Same capped comparator, frozen AC responses and U,E. Conditional original physical contracts; no new response columns. True metric obstruction applies only to scalar-norm receiver.'},indent=2)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
