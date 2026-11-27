#!/usr/bin/env python3
"""Conditioned joint E850 correction of frozen E365 normalized trials.
Numerical block proposals only; exact rational interval evaluation afterwards.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
from cr_alp_123v_conditioned_block_proposal import candidate
# Internal trusted arithmetic may produce >4300-digit exact fractions.
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
PIN='554e7ab54da4545987d89aa660a4abe6c951cb63'
J=2048
def up_sqrt(v):
    assert v>=0
    D=10**110
    return F(isqrt((v*D*D).numerator//(v*D*D).denominator)+1,D)
def main(root,prior,quant,out):
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==PIN
    src=raw.decode()
    def replace(old,new):
        nonlocal src
        assert src.count(old)==1,old
        src=src.replace(old,new)
    replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=365')
    replace('def psi_pair(x,y=ZERO):','''def atan_wide(x):
 z=ZERO
 assert x.l>=0
 for _ in range(8):
  if x.up()<=F(1,3):return z+atan(x)
  assert x.l>0 and x.l+x.h>0
  z+=atan(I(F(1,4)));x=(x-F(1,4))/(1+x/4)
  assert x.l>=0
 raise AssertionError("atan range reduction")
def psi_pair(x,y=ZERO):''')
    replace('assert y.up()/F(xx.l,S)<F(1,3)','assert y.up()/F(xx.l,S)<3')
    replace('li=atan(y/xx)','li=atan_wide(y/xx)')
    # Reuse the original coefficient loop for the whole prefix, retaining N=365
    # for dense matrix construction. No large dense E2048 matrix is assembled.
    replace('range(1,N+1)};trig=', 'range(1,2049)};trig=')
    replace('for j in ids};Hb={}\n for j in ids:', 'for j in range(k,2049,2)};Hb={}\n for j in range(k,2049,2):')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    replace(old,old+', Bmass, Hbb, c, Hb, Pbb')
    g={'__file__':str(SRC),'__name__':'cr_alp_123i'}
    exec(compile(src,str(SRC),'exec'),g)
    I,S,ZERO=(g[x] for x in ('I','S','ZERO'))
    a,L,nu,TS,W,trig,rt,eh,em,qt=(g[x] for x in ('A','L','nu','TS','W','trig','rt','eh','em','qt'))
    assert len(nu)==J and len(g['AI'])==J
    assert sum(w.up() for w in W)<4
    def iv(v):return I(*map(F,v))
    def pair(v):return [str(F(v.l,S)),str(F(v.h,S))]
    def overlap(x,v):return max(F(x.l,S),F(v[0]))<=min(F(x.h,S),F(v[1]))
    rawbank=(root/'cr_alp_046/CR-ALP-046_full_E85_intervals.json').read_bytes()
    assert hashlib.sha256(rawbank).hexdigest()=='d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22'
    frozen=json.loads(rawbank)
    traw=(root/'cr_alp_057/trial.json').read_bytes()
    assert hashlib.sha256(traw).hexdigest()=='6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0'
    tb=json.loads(traw);ts=int(tb['scale'])
    rawprior=prior.read_bytes();assert hashlib.sha256(rawprior).hexdigest()=='09b34469b32c13ebf5f5b690d1eb9ad3a9fb5676f171a1c7f85325cedc199d7b'
    old=json.loads(rawprior);assert old['status']=='CONDITIONAL_EXPANDED_JOINT_AUXILIARY_TRIAL_BANK'
    assert set(old['sectors'])=={'1','2'} and all(len([r for r in v['rows'] if 'joint_trial' in r])==3 for v in old['sectors'].values())
    rawquant=quant.read_bytes();assert hashlib.sha256(rawquant).hexdigest()=='f5feb77462ff530ba85f3ee242e13b8e11b71144a65dfd3e16893188b05315de'
    qp=json.loads(rawquant);assert qp['status']=='CONDITIONAL_QUANTITATIVE_COLLAR_RECEIVER'
    bands=json.loads(SRC.with_name('CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    cap=next(dg(b) for b in bands if b['start']==2049)
    assert all(dg(b)>=F(72,25) for b in bands) and all(dg(b)>=cap for b in bands if b['start']>=2049)
    def diag(i):return min(cap,dg(next(b for b in bands if b['start']<=i and (b['end'] is None or i<=b['end']))))
    AP={j:sum((w*t[0] for w,t in zip(W,trig[j])),ZERO) for j in nu}
    AG={j:g['AI'][j]-AP[j] for j in nu}
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*ct+st/(a*u)) for w,t,(st,ct) in zip(W,TS,trig[i])),ZERO)
        return -2/a*(v*AP[i]-u*AP[j])/(v.sq()-u.sq())
    hss=g['Hss']
    ans={}
    for bk in frozen['blocks']:
        k=bk['parity'];ks=str(k);order,H,Mass,B,Hbb,c,Hb,Pbb=g['build'](k)
        low=bk['ids'];ix=[order.index(j) for j in low]
        for mat,key in [(H,'H'),(Mass,'M')]:
            for i,u in enumerate(ix):
                for j,v in enumerate(ix):
                    l,h=bk[key][i][j]
                    assert max(mat[u][v].l,l)<=min(mat[u][v].h,h)
        b0=(eh+em)/2 if k==1 else -(eh-em)/2
        b1=-(eh-em)/4 if k==1 else (eh+em)/4
        ep=(b0+2*b1)/2;en=(b0-2*b1)/2;Cb=2*b0/rt
        assert Cb.up()<3 and b0.up()<F(6,5) and b1.up()<F(3,5)
        def pb(j):
            v=nu[j];val=ZERO
            for nn,w,(st,ct) in zip(g['NS'],W,trig[j]):
                sn=I(nn).sqrt();bt=ep*sn+en/sn
                val-=2*w/(rt*(v.sq()+qt))*(v*(bt+b0*ct)-b1*st)
            return val
        Pb={j:pb(j) for j in c}
        st=B-sum((c[j].sq() for j in low),ZERO);assert st.l>0
        high=[j for j in c if j>85]
        weighted=I(cap)*st+sum((I(diag(j)-cap)*c[j].sq() for j in high),ZERO)
        active=[j for j in c if j<=850]
        hs={j:Hb[j]-sum((c[i]*hss(i,j) for i in low),ZERO) for j in active}
        ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in active}
        hb=Hbb-sum((c[i]*Hb[i] for i in low),ZERO)
        phb=Pbb-sum((c[i]*Pb[i] for i in low),ZERO)
        result=[]
        for row in old['sectors'][ks]['rows']:
            if row['cutoff']!=365:continue
            w=row['joint_trial'];d={int(j):iv(v) for j,v in w['source_sines'].items()}
            q={int(j):iv(v) for j,v in w['tail_precursor'].items()}
            fb=iv(w['source_pole']);tau=iv(w['tail_projection_tau'])
            t={j:d[j]+tau*c[j] for j in d};t.update(q);T=fb-tau
            assert set(d)==set(low) and all(j>85 for j in q)
            # Independently reconstruct exact frozen trial coordinates before
            # intersecting the two sets of valid interval enclosures.
            coords=list(map(F,w['trial_rational_coordinates']))
            order0=w['trial_order'];assert order0==low+list(q)+['normalized_pole']
            assert len(coords)==len(order0) and coords[-1]==1
            ww={j:I(v) for j,v in zip(order0[:-1],coords[:-1])}
            dr={j:ww[j] for j in low}
            dr[k]=dr[k]-sum((ww[j]*c[j]/c[k] for j in low if j!=k),ZERO)
            qr={j:ww[j] for j in q}
            fr=(1-ww[k]*c[k])/B
            tr=sum((qr[j]*c[j] for j in qr),ZERO)/st
            def intersect(v,w):
                return I(max(v.l,w.l),min(v.h,w.h),True)
            d={j:intersect(d[j],dr[j]) for j in d}
            q={j:intersect(q[j],qr[j]) for j in q}
            fb=intersect(fb,fr);tau=intersect(tau,tr)
            t={j:d[j]+tau*c[j] for j in d};t.update(q);T=fb-tau
            assert overlap(-tau,w['tail_full_pole'])
            def evaluate(q,tau,indices):
                t={j:d[j]+tau*c[j] for j in d};t.update(q);T=fb-tau
                alpha=(sum((d[j]*(hs[j]-ps[j]) for j in d),ZERO)+fb*(hb-phb)
                      +sum((t[j]*ps[j] for j in t),ZERO)+T*phb
                      +sum((q[j]*I(diag(j))*c[j] for j in q),ZERO)-tau*weighted)/st
                e={}
                for j in indices:
                    # Exact projected high-sine auxiliary residual.
                    e[j]=(sum((d[i]*(hss(i,j)-pss(i,j)) for i in d),ZERO)
                        +sum((t[i]*pss(i,j) for i in t),ZERO)
                        +fb*(Hb[j]-Pb[j])+T*Pb[j]+I(diag(j))*(q.get(j,ZERO)-tau*c[j])-alpha*c[j])
                return e,alpha,t,T
            bank=[j for j in high if j<=850]
            eold,_,_,_=evaluate(q,tau,bank)
            for j,v in w['available_projected_auxiliary_residual_sine_coefficients'].items():assert overlap(eold[int(j)],v)
            # Propose a correction in the well-conditioned compressed C matrix.
            phh=phb-sum((c[i]*ps[i] for i in low),ZERO)
            cp={j:c[j]/st for j in bank}
            def c_entry(i,j):
                return (pss(i,j)-cp[i]*ps[j]-cp[j]*ps[i]+cp[i]*cp[j]*phh
                   +(I(diag(i)) if i==j else ZERO)-cp[i]*I(diag(j))*c[j]
                   -cp[j]*I(diag(i))*c[i]+cp[i]*cp[j]*weighted)
            # Complete raw mixed covectors, in the exact 046 low basis.
            # No conversion to floating point before the V congruence.
            n=len(low)
            H0mat=[[H[u][v] for v in ix] for u in ix]
            Rbank=[]
            for j in bank:
                vals={i:hss(i,j)-cp[j]*hs[i] for i in low}
                pole=Hb[j]-cp[j]*hb
                Rbank.append([vals[i]-c[i]/c[k]*vals[k] if i!=k else vals[k]-c[k]/B*pole for i in low])
            polecross=[(Hb[i]-c[i]/c[k]*Hb[k])/B if i!=k else (Hb[k]-c[k]/B*Hbb)/B for i in low]
            x0={j:v for j,v in zip(order0[:-1],coords[:-1]) if j in low}
            qrat={j:v for j,v in zip(order0[:-1],coords[:-1]) if j>85}
            sold=[polecross[i]+sum((H0mat[i][l]*I(x0[j]) for l,j in enumerate(low)),ZERO)
                   +sum((r[i]*I(qrat.get(j,F(0))) for j,r in zip(bank,Rbank)),ZERO) for i in range(n)]
            vraw=(root/'cr_alp_049'/f'result_{k}.json').read_bytes()
            assert hashlib.sha256(vraw).hexdigest()=={1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k]
            vb=json.loads(vraw);vs=int(vb['scale']);assert S%vs==0
            VV=[[I(int(l)*(S//vs),int(h)*(S//vs),True) for l,h in row] for row in vb['V']]
            oldcond=[sum((VV[i][j]*sold[i] for i in range(n)),ZERO) for j in range(n)]
            assert overlap(sum((v.sq() for v in oldcond),ZERO),w['finite_conditioned_residual_norm2'])
            Cbank=[[ZERO]*len(bank) for _ in bank]
            for ii,i in enumerate(bank):
                for jj in range(ii,len(bank)):
                    Cbank[ii][jj]=Cbank[jj][ii]=c_entry(i,bank[jj])
            dz,dq,zz=candidate(H0mat,Rbank,Cbank,[eold[j] for j in bank],sold,VV,S,I,ZERO)
            del Cbank
            xnew={j:x0[j]+v for j,v in zip(low,dz)}
            qnew={j:qrat.get(j,F(0))+v for j,v in zip(bank,dq)}
            d={j:I(v) for j,v in xnew.items()}
            d[k]-=sum((I(xnew[j])*c[j]/c[k] for j in low if j!=k),ZERO)
            fb=(1-I(xnew[k])*c[k])/B
            q={j:I(v) for j,v in qnew.items()}
            tau=sum((q[j]*c[j] for j in bank),ZERO)/st
            e,alpha,t,T=evaluate(q,tau,high)
            snew=[sold[i]+sum((H0mat[i][l]*I(v) for l,v in enumerate(dz)),ZERO)
                  +sum((r[i]*I(v) for r,v in zip(Rbank,dq)),ZERO) for i in range(n)]
            cond=[sum((VV[i][j]*snew[i] for i in range(n)),ZERO) for j in range(n)]
            # Exact polarization of energy change along the full block direction.
            change=sum((I(v)*(a+b) for v,a,b in zip(dz,sold,snew)),ZERO)
            change+=sum((I(v)*(eold[j]+e[j]) for v,j in zip(dq,bank)),ZERO)
            energy=iv(w['joint_trial_energy'])+change
            defect=sum(F(e[j].sq().h,S) for j in bank)
            print('conditioned joint correction',k,'energy',float(F(energy.l,S)),float(F(energy.h,S)),
                  'bank defect squared upper',float(defect),'finite conditioned norm2 upper',
                  float(sum(F(v.sq().h,S) for v in cond)),flush=True)
            # 123H signed geometric expansion: use j>J and v_j>=15j/11.
            N0=32;V=F(15,11)*J;rho=(max(nu[j].up() for j in t)/V)**2;assert rho<1
            moments={};pows={}
            def charge(s,z):
                assert z>=0
                pows[s]=pows.get(s,F(0))+z
            def expand(label,x,mult,shift,skip):
                M=[]
                for l in range(N0):
                    # Power via outward interval multiplication, never endpoint replacement.
                    acc=ZERO
                    for j,v in x.items():
                        z=I(1)
                        for _ in range(l):z=z*nu[j].sq()
                        acc+=v*z
                    M.append(acc)
                    if l>=skip:charge(2*l+2-shift,mult*acc.up())
                H=sum(v.up()*nu[j].up()**(2*N0) for j,v in x.items())
                charge(2*N0+2-shift,mult*H/(1-rho))
                moments[label]=[pair(v) for v in M]
                return M
            Li={j:d.get(j,ZERO)*AG[j]+t[j]*AP[j] for j in t}
            ML0=sum(Li.values(),ZERO)
            expand('centered_combined',{j:v*(nu[j].sq()+qt) for j,v in Li.items()},F(20,11),-1,0)
            expand('gamma_sine',{j:v*nu[j] for j,v in d.items()},F(80,11),0,0)
            expand('prime_sine',{j:v*nu[j] for j,v in t.items()},F(80,11),0,0)
            H0=-2/a*ML0-(tau*I(cap)+alpha)*Cb
            charge(1,H0.up()+80*fb.up()+40*T.up())
            logs=2*fb.up()
            tail=sum(z*(F(11,15)**s)*up_sqrt(F(J,2*s-1))/J**s for s,z in pows.items())
            # log(3J) < 10 for J=2048, checked using original outward log.
            assert g['log'](I(3*J)).h<10*S
            tail+=logs*F(11,15)*up_sqrt(F(122,J))
            # Full weighted pole-null inverse comparison (123F), retaining a
            # safe lower bound on the weighted pole inner-product modulus.
            dc=cap-F(21,8)
            prefix=sum(F(e[j].sq().h,S)/(diag(j)-F(21,8)) for j in high)
            mass=st-sum((c[j].sq() for j in high),ZERO)
            assert mass.h>=0
            massup=F(mass.h,S)
            ip=sum((c[j]*e[j]/I(diag(j)-F(21,8)) for j in high),ZERO)
            iplo=max(F(0),F(ip.l,S),-F(ip.h,S))
            iplo=max(F(0),iplo-tail*up_sqrt(massup)/dc)
            denup=sum(F(c[j].sq().h,S)/(diag(j)-F(21,8)) for j in high)+massup/dc
            assert denup>0
            U=prefix+tail*tail/dc-iplo*iplo/denup
            assert U>=0
            E0=F(energy.l,S)
            # b_y = conditioned joint finite residual - X*e. X's complete pole
            # term is orthogonal to e; only finite precursor coefficients matter.
            tt=tb['sectors'][ks];assert tt['ids']==list(range(k,851,2))
            X=[]
            for col in range(len(low)):
                X.append(sum((I(F(int(v[0]),ts),F(int(v[1]),ts))*e[j] for j,rr in zip(tt['ids'],tt['D']) if j>85 for v in [rr[col]]),ZERO))
            X=[v-r for v,r in zip(X,cond)]
            xnorm=up_sqrt(sum(F(v.sq().h,S) for v in X))
            bnorm=xnorm
            pp=qp['sectors'][ks];m=F(pp['finite_gap_lower']);kap=F(pp['kappa_upper'])
            delta=F(1,10**5) if k==1 else F(21,10**5)
            A=1/F(pp['lambda_K_lower'])+1/delta;Bfac=kap/delta
            loss=U+(up_sqrt(A)*bnorm+up_sqrt(Bfac*U))**2
            lower=E0-loss
            result.append({'cutoff':850,'source_cutoff':row['cutoff'],'finite_prefix':J,
               'joint_energy_interval':pair(energy),'bank_defect_norm2_upper':str(defect),
               'low_trial_coordinates':{str(j):str(v) for j,v in xnew.items()},
               'conditioned_low_correction_coordinates':[str(v) for v in zz],
               'corrected_tail_coordinates':{str(j):str(v) for j,v in qnew.items()},
               'conditioned_joint_residual_vector':[pair(v) for v in cond],'alpha':pair(alpha),'leading_combined_moment':pair(H0),
               'source_pole':pair(fb),'combined_pole':pair(T),'signed_moments':moments,
               'infinite_tail_norm_upper':str(tail),'weighted_inverse_error_upper':str(U),
               'weighted_prefix_upper':str(prefix),'weighted_pole_subtraction_lower':str(iplo*iplo/denup),
               'correction_response_vector':[pair(v) for v in X],
               'conditioned_joint_residual_norm_upper':'0',
               'finite_corrected_residual_norm_upper':str(bnorm),'joint_energy_lower':str(E0),
               'response_loss_upper':str(loss),'conditional_pole_schur_lower':str(lower),
               'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if lower>0 else 'SUFFICIENT_BOUND_INCONCLUSIVE'})
            print('parity',k,'cutoff',row['cutoff'],'tail upper',float(tail),'U',float(U),'loss/energy',float(loss/E0) if E0>0 else None,result[-1]['classification'],flush=True)
        ans[ks]=result
    out.write_text(json.dumps({'status':'CONDITIONAL_CONDITIONED_JOINT_CORRECTION','sectors':ans,
        'scope':'Classical coefficient and comparator premises. Positive lower bound would certify corresponding conditional complete pole sector; inconclusive bound is not rejection.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
