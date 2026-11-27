#!/usr/bin/env python3
"""123AQ2: frozen response-column and direct corrected-energy audit.
No historical K,N, inverse-order certificate or response solve.
Independent shift 768 at original rigorous grid; new 10^70 prefix grid.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
import numpy as np
from cr_alp_123ap3_high_precision_replay import namespace,blob,read
from cr_alp_123y_signed_pole_tail import leading
from cr_alp_123ak_actual_schur_metric import box,mat,tr,mul,add,sym_bounds,T as MG
from cr_alp_123an_independent_receiver_audit import pd,kernel,trans,ab,box as endpoint_box,mul as ivmul,up,rt as metric_root
sys.set_int_max_str_digits(0)
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
J=2048
def up_sqrt(v):
    assert v>=0
    D=10**110;z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)
def main(root,acpath,akpath,alpath,apath,out):
    root=Path(root)
    A=read(acpath,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    K=read(akpath,'a2c465dd692b0afbc0cd1a15b9f3599864085101bb8f52f7d225c1125adf453c')
    AL=read(alpath,'b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e')
    P=read(apath,'4b5e879c5a19958c0192bc9f127b3983565309663bbf926d4a51f1660fbca222')
    assert all(set(v['sectors'])=={'1','2'} for v in (A,K,AL,P))
    g=namespace();KG=kernel()
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
    bands=json.loads(blob(SRC.with_name('CR-ALP-068_pole_centered_gamma_result.json'),'1a7fe8629daa8b8e6cec185a79b6746d6b28975d'))['bands']
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
        for physical_matrix,key in [(H,'H'),(Mass,'M')]:
            for i,u in enumerate(ix):
                for j,v in enumerate(ix):
                    l,h=bk[key][i][j]
                    assert max(physical_matrix[u][v].l,l)<=min(physical_matrix[u][v].h,h)
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
        active=list(c)
        hs={j:Hb[j]-sum((c[i]*hss(i,j) for i in low),ZERO) for j in active}
        ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in active}
        hb=Hbb-sum((c[i]*Hb[i] for i in low),ZERO)
        phb=Pbb-sum((c[i]*Pb[i] for i in low),ZERO)
        n=len(low);bank=[j for j in high if j<=850]
        ac=A['sectors'][ks];ak=K['sectors'][ks];al=AL['sectors'][ks];ap=P['sectors'][ks]
        assert ap['classification']=='FROZEN_TRIAL_DIRECT_REPLAY_PASS'
        assert bank==ac['bank_indices']
        perm=ac['coordinate_order'];assert sorted(perm)==list(range(n))
        LL=[[F(v) for v in row] for row in ac['LDL_lower']]
        assert all(len(row)==n and row[i]==1 and all(v==0 for v in row[i+1:]) for i,row in enumerate(LL))
        w=list(map(F,ak['frame_roots']));assert len(w)==n and min(w)>0
        QG=int(ac['response_coordinate_grid']);PG=10**70;Q=ac['response_transformed_coordinates']
        assert len(Q)==len(bank) and all(len(row)==n for row in Q)
        Qi=np.array(Q,dtype=object);absQi=np.abs(Qi)
        def ib(v):return box(F(v.l,S),F(v.h,S))
        def grid(v):
            m=round(F(v.l+v.h,2*S)*PG)
            rr=max(F(m,PG)-F(v.l,S),F(v.h,S)-F(m,PG))
            return m,-((-rr*PG).__floor__())
        vraw=(root/'cr_alp_049'/f'result_{k}.json').read_bytes()
        assert hashlib.sha256(vraw).hexdigest()=={1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k]
        vb=json.loads(vraw);vs=int(vb['scale'])
        VV=[[I(F(int(l),vs),F(int(h),vs)) for l,h in row] for row in vb['V']]
        LT=[[I(x) for x in row] for row in LL]
        assert all(v.l==v.h for row in VV for v in row),'V is not point rational'
        Vt=[[sum((VV[i][perm[l]]*LT[l][j] for l in range(n)),ZERO) for j in range(n)] for i in range(n)]
        phh=phb-sum((c[i]*ps[i] for i in low),ZERO)
        cp={j:c[j]/st for j in high}
        def centry(i,j):
            return (pss(i,j)-cp[i]*ps[j]-cp[j]*ps[i]+cp[i]*cp[j]*phh
               +(I(diag(i)) if i==j else ZERO)-cp[i]*I(diag(j))*c[j]
               -cp[j]*I(diag(i))*c[i]+cp[i]*cp[j]*weighted)
        def rrow(j):
            vals={i:hss(i,j)-cp[j]*hs[i] for i in low}
            pole=Hb[j]-cp[j]*hb
            raw=[vals[i]-c[i]/c[k]*vals[k] if i!=k else vals[k]-c[k]/B*pole for i in low]
            return [sum((v*Vt[i][l] for i,v in enumerate(raw)),ZERO) for l in range(n)]
        prefix=[F(0)]*n;bankmax=F(0);prefix_receipts=[];overlaps=0;RB=[];CB=[];oldprefix=dict(ac['prefix_mid_radius'])
        for j in high:
            cg=[grid(centry(j,l)) for l in bank]
            cm=np.array([v[0] for v in cg],dtype=object)
            cr=np.array([v[1] for v in cg],dtype=object)
            zz=cm@Qi;rad=cr@absQi
            if j in bank:CB.append([ib(centry(j,l)) for l in bank])
            row=[];rj=rrow(j)
            if j in bank:RB.append([ib(v) for v in rj])
            for l in range(n):
                z=rj[l]-I(F(int(zz[l]-rad[l]),PG*QG),F(int(zz[l]+rad[l]),PG*QG))
                m,r=grid(z);upper=F(abs(m)+r,PG)
                prefix[l]+=upper**2/(diag(j)-F(21,8))
                if j<=850:bankmax=max(bankmax,upper)
                row.append([m,r])
                om,orr=oldprefix[j][l];opg=int(ac['prefix_grid'])
                assert F(m-r,PG)<=F(om+orr,opg) and F(om-orr,opg)<=F(m+r,PG),(k,j,l,'prefix overlap')
                overlaps+=1
            prefix_receipts.append([j,row])
            if j%100 in (k%2,):print('prefix',k,j,flush=True)
        print('response bank',k,'max defect',float(bankmax),'weighted prefix trace',float(sum(prefix)),flush=True)
        # Each transformed column Z is Q(Gf+P(f+q)+Dcap q), q=-Mhat.
        # f=V L column, pole-null; keep signed pole and centered moments.
        J0=J-1;Vmin=F(15,11)*J;rho=(max(nu[j].up() for j in bank)/Vmin)**2;assert rho<1
        powers={i:[I(1)] for i in low+bank}
        for i in powers:
            for l in range(32):powers[i].append(powers[i][-1]*nu[i].sq())
        tails=[];receipts=[];remcharges=[];momentchecks=0
        frozen_prefix=dict(prefix_receipts)
        for col in range(n):
            d={j:Vt[i][col] for i,j in enumerate(low)}
            d[k]-=sum((d[j]*c[j]/c[k] for j in low if j!=k),ZERO)
            fb=-Vt[low.index(k)][col]*c[k]/B
            q={j:I(F(-Q[i][col],QG)) for i,j in enumerate(bank)}
            tau=sum((q[j]*c[j] for j in bank),ZERO)/st
            t={j:d[j]+tau*c[j] for j in low};t.update(q);T=fb-tau
            alpha=(sum((d[j]*(hs[j]-ps[j]) for j in d),ZERO)+fb*(hb-phb)
                   +sum((t[j]*ps[j] for j in t),ZERO)+T*phb
                   +sum((q[j]*I(diag(j))*c[j] for j in q),ZERO)-tau*weighted)/st
            for j in (bank[0],bank[-1],high[-1]):
                direct=(sum((d[i]*(hss(i,j)-pss(i,j)) for i in d),ZERO)
                    +sum((t[i]*pss(i,j) for i in t),ZERO)+fb*(Hb[j]-Pb[j])+T*Pb[j]
                    +I(diag(j))*(q.get(j,ZERO)-tau*c[j])-alpha*c[j])
                m0,r0=frozen_prefix[j][col]
                assert max(F(direct.l,S),F(m0-r0,PG))<=min(F(direct.h,S),F(m0+r0,PG)),(k,col,j)
            Li={j:d.get(j,ZERO)*AG[j]+t[j]*AP[j] for j in t}
            H0=-2/a*sum(Li.values(),ZERO)-(tau*I(cap)+alpha)*Cb
            pows={};moments={};remcharges_col=[]
            def charge(s,z):pows[s]=pows.get(s,F(0))+z
            def series(label,x,mult,shift):
                mm=[]
                for l in range(32):
                    v=sum((x[j]*powers[j][l] for j in x),ZERO)
                    charge(2*l+2-shift,mult*v.up());mm.append(pair(v))
                    assert overlap(v,ac['tail_receipts'][col]['moments'][label][l]),(k,col,label,l)
                H=sum(x[j].up()*powers[j][32].up() for j in x)
                rem=mult*H/(1-rho);charge(66-shift,rem);moments[label]=mm
                p=66-shift;remcharges_col.append(rem*up_sqrt(F(J0,2*(2*p-1)))/F(15*J0,11)**p)
            series('centered_combined',{j:v*(nu[j].sq()+qt) for j,v in Li.items()},F(20,11),-1)
            series('gamma_sine',{j:v*nu[j] for j,v in d.items()},F(80,11),0)
            series('prime_sine',{j:v*nu[j] for j,v in t.items()},F(80,11),0)
            ss=sum(v*up_sqrt(F(J0,2*(2*s-1)))/F(15*J0,11)**s for s,v in pows.items())
            pole,receipt=leading(g,k,J,fb,T,H0)
            tail=ss+pole;tails.append(tail);remcharges.append(sum(remcharges_col));momentchecks+=96
            for name,v in [('source_pole',fb),('combined_pole',T),('alpha',alpha),('leading_combined_moment',H0)]:assert overlap(v,ac['tail_receipts'][col][name]),(k,col,name)
            receipts.append({'moments':moments,'pole':receipt,'source_pole':pair(fb),'combined_pole':pair(T),'alpha':pair(alpha),'leading_combined_moment':pair(H0),'sine_series_norm_upper':str(ss),'tail_norm_upper':str(tail)})
            print('column',k,col,'tail',float(tail),'prefix',float(prefix[col]),flush=True)
        # Corrected energy directly from H, R_B, C_B and frozen Q:
        # Hhat = H - R_B^T Q - Q^T R_B + Q^T C_B Q.
        # Compute in the unscaled AC frame, then apply the frozen roots.
        Vmat=mat([[ib(v) for v in row] for row in Vt])
        H0=mat([[ib(H[i][j]) for j in ix] for i in ix])
        Ht=mul(mul(tr(Vmat),H0),Vmat)
        Rm=mat(RB);Cm=mat(CB);Qm=mat([[box(F(v,QG),F(v,QG)) for v in row] for row in Q])
        QR=mul(tr(Qm),Rm)
        Hhat=add(add(add(Ht,QR,-1),tr(QR),-1),mul(tr(Qm),mul(Cm,Qm)))
        Wm=mat([[box(w[i] if i==j else F(0),w[i] if i==j else F(0)) for j in range(n)] for i in range(n)])
        Hlo,Hup=sym_bounds(mul(mul(Wm,Hhat),Wm))
        Zm=mat([[box(F(m-r,PG)*wi,F(m+r,PG)*wi) for (m,r),wi in zip(row,w)] for j,row in prefix_receipts])
        Zd=mat([[box(F(m-r,PG)*wi/(diag(j)-F(21,8)),F(m+r,PG)*wi/(diag(j)-F(21,8))) for (m,r),wi in zip(row,w)] for j,row in prefix_receipts])
        _,Pup=sym_bounds(mul(tr(Zm),Zd))
        BG=int(al['whitening_grid']);BW=np.array([[int(v) for v in row] for row in al['whitening_matrix']],dtype=object)
        eye=np.eye(n,dtype=object);p=F(al['prefix_relative_upper'])
        assert pd(Hlo) and pd(BW@Hlo@BW.T-MG*BG*BG*eye) and pd(p.numerator*Hlo-p.denominator*Pup)
        print(k,'direct corrected-energy / prefix matrix tests PASS',flush=True)
        # Recompute the metric-first tail directly; no subtraction of old
        # receipt charges. Geometric remainders are explicit above.
        series=F(0);momentbounds=[]
        for label,mult,shift in [('centered_combined',F(20,11),-1),('gamma_sine',F(80,11),0),('prime_sine',F(80,11),0)]:
            for l in range(32):
                power=2*l+2-shift
                fac=mult*metric_root(F(J0,2*(2*power-1)))/F(15*J0,11)**power
                if label=='gamma_sine':fac/=5
                rows=[endpoint_box(F(t['moments'][label][l][0])*wi*fac,F(t['moments'][label][l][1])*wi*fac) for wi,t in zip(w,receipts)]
                v=trans(rows,BW,BG);bound=metric_root(sum(ab(x)**2 for x in v))
                series+=bound;momentbounds.append([label,l,str(up(bound))])
        # Bound unknown vector remainder using absolute whitening weights.
        rem0=[w[i]*(remcharges[i]+F(t['pole']['pole_remainder_norm_upper'])) for i,t in enumerate(receipts)]
        rem=metric_root(sum(sum(abs(int(BW[i,j]))*rem0[j]/BG for j in range(n))**2 for i in range(n)))
        coeff=[]
        for i in range(9):
            rows=[]
            for wi,t in zip(w,receipts):
                po=t['pole'];x,y=map(F,[po['constant'],po['log_coefficient'],*po['cosine_coefficients']][i])
                rows.append(endpoint_box(x*wi,y*wi))
            coeff.append(trans(rows,BW,BG))
        lead=F(0)
        for j in range(n):
            for i in range(9):
                for h in range(9):lead+=ivmul(ivmul(coeff[i][j],KG[i][h]),coeff[h][j])[1]
        lead=up(lead);assert lead>=0
        dc=cap-F(21,8)
        rho=up(p+(metric_root(lead)+series+rem)**2/dc,10**6);assert rho<1
        # Reuse independently reconstructed AP3 trial, not historical b2.
        E,U=map(F,(ap['E'],ap['U']));assert E>U>0
        s0=[iv(v) for v in ap['conditioned_residual_intervals']];eb=F(ap['bank_defect_norm2_upper'])
        bf=F(0)
        for j in range(n):
            v=sum((s0[perm[i]]*LT[i][j] for i in range(n)),ZERO)
            norm2=F(sum(int(row[j])**2 for row in Q),QG*QG)
            bf+=(w[j]*(v.up()+metric_root(eb*norm2)))**2
        frob2=F(sum(int(v)**2 for row in BW for v in row),BG*BG)
        bnew=frob2*bf/(1-rho);eta=rho/(1-rho)
        loss=up((U+(metric_root(bnew)+metric_root(eta*U))**2)/E)
        assert overlaps==len(high)*n and momentchecks==96*n
        ans[ks]={'dimension':n,'prefix_overlap_checks':overlaps,'moment_overlap_checks':momentchecks,'matrix_grid':str(MG),
            'corrected_energy_lower':[[str(int(v)) for v in row] for row in Hlo],
            'prefix_inverse_gram_upper':[[str(int(v)) for v in row] for row in Pup],
            'prefix_grid':PG,'prefix_mid_radius':prefix_receipts,'tail_receipts':receipts,
            'geometric_remainder_bounds':[str(v) for v in remcharges],'metric_moment_bounds':momentbounds,
            'prefix_relative_upper':str(p),'leading_trace_upper':str(lead),'series_upper':str(up(series)),
            'remainder_upper':str(up(rem)),'rho_upper':str(rho),'E':str(E),'U':str(U),
            'finite_charge_over_E_upper':str(up(bnew/E,10**30)),'loss_over_E_upper':str(loss),
            'classification':'DIRECT_RESPONSE_AUDIT_PASS' if loss<1 else 'DIRECT_RESPONSE_AUDIT_INCONCLUSIVE'}
        payload=json.dumps({'status':'CONDITIONAL_DIRECT_RESPONSE_AUDIT','sectors':ans,'scope':'Frozen response columns and AP3 trials, original analytic coefficient and comparator premises retained.'},indent=2)+chr(10)
        assert set(json.loads(payload)['sectors'])==set(ans);Path(out).write_text(payload)
        print(k,'rho',float(rho),'loss/E',float(loss),ans[ks]['classification'],flush=True)
if __name__=='__main__':main(*sys.argv[1:])
