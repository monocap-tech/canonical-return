#!/usr/bin/env python3
"""CR-ALP-123AC: exact weighted response-column defect receiver.
Numerical solve proposes rational columns in a rational LDL-transformed frame.
All acceptance arithmetic is outward rational interval evaluation.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
import numpy as np
from cr_alp_123y_signed_pole_tail import leading
from cr_alp_075_shifted_energy import lower_matrix,S as MS
from cr_alp_074_balanced_selection import checked,HASHES
from cr_alp_123j_common_denominator_solve import solve_common
from cr_alp_123l_matrix_order_receiver import positive_shift,HASH,DELTA
# Internal trusted arithmetic may produce >4300-digit exact fractions.
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
PIN='554e7ab54da4545987d89aa660a4abe6c951cb63'
J=2048
def up_sqrt(v):
    assert v>=0
    D=10**110
    return F(isqrt((v*D*D).numerator//(v*D*D).denominator)+1,D)
def main(root,custody,terminal,prior,quant,out):
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==PIN
    src=raw.decode()
    def replace(old,new):
        nonlocal src
        assert src.count(old)==1,old
        src=src.replace(old,new)
    replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=85')
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
    rawprior=prior.read_bytes();assert hashlib.sha256(rawprior).hexdigest()=='7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897'
    old=json.loads(rawprior);assert old['status']=='CONDITIONAL_CONDITIONED_JOINT_CORRECTION'
    frozen_terminal=json.loads(terminal.read_text());assert frozen_terminal['rank']==85
    permutation=json.loads((custody/'replayed-082/result.json').read_text())
    grams=json.loads((custody/'replayed-091/result.json').read_text())
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
        active=list(c)
        hs={j:Hb[j]-sum((c[i]*hss(i,j) for i in low),ZERO) for j in active}
        ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in active}
        hb=Hbb-sum((c[i]*Hb[i] for i in low),ZERO)
        phb=Pbb-sum((c[i]*Pb[i] for i in low),ZERO)
        n=len(low);bank=[j for j in high if j<=850]
        pp=qp['sectors'][ks];assert pp['dimension']==n
        perm=permutation['sectors'][ks]['coordinate_order'];assert sorted(perm)==list(range(n))
        ww=checked(root/'cr_alp_056'/f'result_extended_{ks}.json',HASHES[ks])
        HH=lower_matrix(ww,perm);gg=grams['sectors'][ks];assert int(gg['scale'])==MS
        Kmat=[[HH[i][j]+int(gg['correction_gram_lower_grid'][i][j]) for j in range(n)] for i in range(n)]
        aa=frozen_terminal['sectors'][ks]
        aa=aa['audit'] if k==1 else aa['audits']['second_batch_combined']
        Nmat=[[int(v) for v in row] for row in aa['integer_matrix']]
        assert aa['full_pass'] and aa['dimension']==n
        assert hashlib.sha256(json.dumps(Nmat,separators=(',',':')).encode()).hexdigest()==HASH[ks]==pp['terminal_matrix_sha256']
        assert positive_shift(Kmat,F(pp['lambda_K_lower'])) and positive_shift(Nmat,DELTA[ks])
        # Construct rational Dbar>=K^-1+N^-1 by exact column solves and
        # symmetric grid rounding plus a row-sum reserve.
        cols=[]
        for j in range(n):
            rhs=[int(i==j) for i in range(n)]
            x,dx=solve_common(Kmat,rhs);y,dy=solve_common(Nmat,rhs)
            cols.append([MS*(F(u,dx)+F(v,dy)) for u,v in zip(x,y)])
        assert cols==list(map(list,zip(*cols)))
        G=10**30
        Dbar=[[F(round(cols[i][j]*G)+(n if i==j else 0),G) for j in range(n)] for i in range(n)]
        LL=[[F(int(i==j)) for j in range(n)] for i in range(n)];piv=[]
        for j in range(n):
            d0=Dbar[j][j]-sum(LL[j][l]**2*piv[l] for l in range(j))
            assert d0>0;piv.append(d0)
            for i in range(j+1,n):
                LL[i][j]=(Dbar[i][j]-sum(LL[i][l]*LL[j][l]*piv[l] for l in range(j)))/d0
        assert all(sum(LL[i][l]*piv[l]*LL[j][l] for l in range(min(i,j)+1))==Dbar[i][j] for i in range(n) for j in range(n))
        vraw=(root/'cr_alp_049'/f'result_{k}.json').read_bytes()
        assert hashlib.sha256(vraw).hexdigest()=={1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k]
        vb=json.loads(vraw);vs=int(vb['scale'])
        VV=[[I(F(int(l),vs),F(int(h),vs)) for l,h in row] for row in vb['V']]
        LT=[[I(x) for x in row] for row in LL]
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
        # C_B and R_B are in the projected sine bank, transformed low frame.
        # Candidate solves are untrusted; exact rational coordinates follow.
        def mid(v):return float(F(v.l+v.h,2*S))
        CB=np.empty((len(bank),len(bank)))
        for i,j in enumerate(bank):
            for l in range(i,len(bank)):CB[i,l]=CB[l,i]=mid(centry(j,bank[l]))
        RB=np.array([[mid(v) for v in rrow(j)] for j in bank])
        sol=np.linalg.solve(CB,RB);assert np.all(np.isfinite(sol))
        QG=10**40;PG=10**60
        Q=[[round(F(float(v))*QG) for v in row] for row in sol]
        Qi=np.array(Q,dtype=object);absQi=np.abs(Qi)
        # Precompute rounded C row enclosure, then integer matrix products.
        # C midpoint grid error is charged exactly; R remains original interval.
        def grid(v):
            m=round(F(v.l+v.h,2*S)*PG)
            rr=max(F(m,PG)-F(v.l,S),F(v.h,S)-F(m,PG))
            return m,-((-rr*PG).__floor__())
        prefix=[F(0)]*n;bankmax=F(0);prefix_receipts=[]
        for j in high:
            cg=[grid(centry(j,l)) for l in bank]
            cm=np.array([v[0] for v in cg],dtype=object)
            cr=np.array([v[1] for v in cg],dtype=object)
            zz=cm@Qi;rad=cr@absQi
            row=[];rj=rrow(j)
            for l in range(n):
                z=rj[l]-I(F(int(zz[l]-rad[l]),PG*QG),F(int(zz[l]+rad[l]),PG*QG))
                m,r=grid(z);upper=F(abs(m)+r,PG)
                prefix[l]+=upper**2/(diag(j)-F(21,8))
                if j<=850:bankmax=max(bankmax,upper)
                row.append([m,r])
            prefix_receipts.append([j,row])
            if j%100 in (k%2,):print('prefix',k,j,flush=True)
        print('response bank',k,'max defect',float(bankmax),'weighted prefix trace',float(sum(p*v for p,v in zip(piv,prefix))),flush=True)
        # Each transformed column Z is Q(Gf+P(f+q)+Dcap q), q=-Mhat.
        # f=V L column, pole-null; keep signed pole and centered moments.
        J0=J-1;Vmin=F(15,11)*J;rho=(max(nu[j].up() for j in bank)/Vmin)**2;assert rho<1
        powers={i:[I(1)] for i in low+bank}
        for i in powers:
            for l in range(32):powers[i].append(powers[i][-1]*nu[i].sq())
        tails=[];receipts=[]
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
            Li={j:d.get(j,ZERO)*AG[j]+t[j]*AP[j] for j in t}
            H0=-2/a*sum(Li.values(),ZERO)-(tau*I(cap)+alpha)*Cb
            pows={};moments={}
            def charge(s,z):pows[s]=pows.get(s,F(0))+z
            def series(label,x,mult,shift):
                mm=[]
                for l in range(32):
                    v=sum((x[j]*powers[j][l] for j in x),ZERO)
                    charge(2*l+2-shift,mult*v.up());mm.append(pair(v))
                H=sum(x[j].up()*powers[j][32].up() for j in x)
                charge(66-shift,mult*H/(1-rho));moments[label]=mm
            series('centered_combined',{j:v*(nu[j].sq()+qt) for j,v in Li.items()},F(20,11),-1)
            series('gamma_sine',{j:v*nu[j] for j,v in d.items()},F(80,11),0)
            series('prime_sine',{j:v*nu[j] for j,v in t.items()},F(80,11),0)
            ss=sum(v*up_sqrt(F(J0,2*(2*s-1)))/F(15*J0,11)**s for s,v in pows.items())
            pole,receipt=leading(g,k,J,fb,T,H0)
            tail=ss+pole;tails.append(tail)
            receipts.append({'moments':moments,'pole':receipt,'source_pole':pair(fb),'combined_pole':pair(T),'alpha':pair(alpha),'leading_combined_moment':pair(H0),'sine_series_norm_upper':str(ss),'tail_norm_upper':str(tail)})
            print('column',k,col,'tail',float(tail),'prefix',float(prefix[col]),flush=True)
        dc=cap-F(21,8)
        eprefix=sum(p*v for p,v in zip(piv,prefix))
        etail=sum(p*v*v/dc for p,v in zip(piv,tails))
        eta=eprefix+etail
        row=old['sectors'][ks][0];U=F(row['weighted_inverse_error_upper']);E=F(row['joint_energy_lower'])
        assert E>0 and U>=0
        s0=[iv(v) for v in row['conditioned_joint_residual_vector']]
        eb=F(row['bank_defect_norm2_upper'])
        # b transformed by L. Bound its two terms separately, using actual
        # signed low residual and stored norm of e on E850.
        bparts=[]
        for l in range(n):
            st0=sum((s0[perm[i]]*LT[i][l] for i in range(n)),ZERO)
            norm2=sum(F(r[l],QG)**2 for r in Q)
            bparts.append(st0.up()+up_sqrt(norm2*eb))
        b2=sum(p*v*v for p,v in zip(piv,bparts))
        loss=U+(up_sqrt(b2)+up_sqrt(eta*U))**2
        ans[ks]={'coordinate_order':perm,'Dbar':[[str(v) for v in r] for r in Dbar],
            'LDL_lower':[[str(v) for v in r] for r in LL],'LDL_pivots':[str(v) for v in piv],
            'response_coordinate_grid':QG,'response_transformed_coordinates':Q,'bank_indices':bank,
            'prefix_grid':PG,'prefix_mid_radius':prefix_receipts,
            'column_prefix_inverse_uppers':[str(v) for v in prefix],'tail_receipts':receipts,
            'eta_prefix_upper':str(eprefix),'eta_tail_upper':str(etail),'eta_upper':str(eta),
            'finite_charge_upper':str(b2),'U':str(U),'E':str(E),'loss_upper':str(loss),
            'conditional_pole_schur_lower':str(E-loss),'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE'}
        print('terminal',k,'eta prefix',float(eprefix),'eta tail',float(etail),'eta',float(eta),
              'finite/E',float(b2/E),'loss/E',float(loss/E),ans[ks]['classification'],flush=True)
        # Preserve each parity independently, so a later failure does not
        # erase the first result from runner artifacts.
        out.write_text(json.dumps({'status':'CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER','sectors':ans,
            'scope':'Inherited physical comparator and signed coefficient premises; no arbitrary support conclusion.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
