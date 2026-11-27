#!/usr/bin/env python3
"""123AP3: high-precision repair of frozen Z trial audit.
Original 10^340 grid, independent recurrence shift 768; same classical 046 identities.
No trial optimization, historical compression, or grouped quadrature.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
from cr_alp_123y_signed_pole_tail import leading
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent
J=2048
def sha(b):return hashlib.sha256(b).hexdigest()
def read(p,h):
    b=Path(p).read_bytes();assert sha(b)==h,(p,sha(b));return json.loads(b)
def blob(p,h):
    b=p.read_bytes();assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==h;return b
def root(x):
    assert x>=0;D=10**110;v=x*D*D
    return F(isqrt(v.numerator//v.denominator)+1,D)
def up(x,g=10**9):return F(-((-x*g).__floor__()),g)
def namespace():
    src=ROOT/'cr_alp_046_full_finite.py'
    s=blob(src,'554e7ab54da4545987d89aa660a4abe6c951cb63').decode()
    def replace(a,b):
        nonlocal s
        assert s.count(a)==1,a;s=s.replace(a,b)
    replace('M=90;SHIFT=512;N=85','M=90;SHIFT=768;N=85')
    # Retain the original 10^-290 digamma remainder allowance at grid 10^340.
    replace('def psi_pair(x,y=ZERO):','''def atan_wide(x):
 z=ZERO
 for _ in range(8):
  if x.up()<=F(1,3):return z+atan(x)
  old=x;x=(x-F(1,4))/(1+x/4);z+=atan(I(F(1,4)))
  assert x.l>=0 and x.up()<old.up()
 raise AssertionError("atan range reduction")
def psi_pair(x,y=ZERO):''')
    replace('assert y.up()/F(xx.l,S)<F(1,3)','assert y.up()/F(xx.l,S)<3')
    replace('li=atan(y/xx)','li=atan_wide(y/xx)')
    replace('range(1,N+1)};trig=','range(1,2049)};trig=')
    replace('for j in ids};Hb={}\n for j in ids:','for j in range(k,2049,2)};Hb={}\n for j in range(k,2049,2):')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    replace(old,old+', Bmass, Hbb, c, Hb, Pbb')
    g={'__file__':str(src),'__name__':'cr_alp_123ap_independent'}
    exec(compile(s,str(src),'exec'),g)
    assert g['S']==10**340 and g['SHIFT']==768
    assert abs(g['BER'][180])<10**186
    return g
def main(archive,zpath,acpath,akpath,alpath,aopath,out):
    archive=Path(archive)
    Z=read(zpath,'7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897')
    AC=read(acpath,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    AK=read(akpath,'a2c465dd692b0afbc0cd1a15b9f3599864085101bb8f52f7d225c1125adf453c')
    AL=read(alpath,'b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e')
    N=read(Path(zpath).parent.parent/'cr-alp-123n-expanded-joint-trials/result.json','09b34469b32c13ebf5f5b690d1eb9ad3a9fb5676f171a1c7f85325cedc199d7b')
    # AO is used for the already independently audited rho, not its old E,U,b.
    AO=json.loads(Path(aopath).read_text());assert AO['status']=='CONDITIONAL_DIRECT_SCHUR_ROUTE_PASS'
    frozen=read(archive/'cr_alp_046/CR-ALP-046_full_E85_intervals.json','d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22')
    bands=json.loads(blob(ROOT/'CR-ALP-068_pole_centered_gamma_result.json','1a7fe8629daa8b8e6cec185a79b6746d6b28975d'))['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(a['end']+1==b['start'] for a,b in zip(bands,bands[1:]))
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    cap=next(dg(b) for b in bands if b['start']==2049);dc=cap-F(21,8)
    assert dc==F('1.822938456') and all(dg(b)>=cap for b in bands if b['start']>=2049)
    diagonal={j:min(cap,dg(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end'])))) for j in range(86,2049)}
    assert min(diagonal.values())>=F(72,25)
    g=namespace();I,S,O=(g[x] for x in ('I','S','ZERO'))
    a,L,nu,TS,W,trig,rt,eh,em,qt=(g[x] for x in ('A','L','nu','TS','W','trig','rt','eh','em','qt'))
    assert len(nu)==2048 and sum(w.up() for w in W)<4
    def pair(v):return [str(F(v.l,S)),str(F(v.h,S))]
    def overlap(v,p):return max(F(v.l,S),F(p[0]))<=min(F(v.h,S),F(p[1]))
    AP={j:sum((w*t[0] for w,t in zip(W,trig[j])),O) for j in nu}
    AG={j:g['AI'][j]-AP[j] for j in nu}
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*co+si/(a*u)) for w,t,(si,co) in zip(W,TS,trig[i])),O)
        return -2/a*(v*AP[i]-u*AP[j])/(v.sq()-u.sq())
    hss=g['Hss'];ans={}
    for bk in frozen['blocks']:
        k=bk['parity'];ks=str(k);z=Z['sectors'][ks][0];ac=AC['sectors'][ks];ak=AK['sectors'][ks];al=AL['sectors'][ks]
        order,H,Mass,B,Hbb,c,Hb,Pbb=g['build'](k)
        low=bk['ids'];assert order==low and z['cutoff']==850 and z['finite_prefix']==J
        high=[j for j in c if j>85];bank=[j for j in high if j<=850]
        assert bank==ac['bank_indices']
        oldS=int(frozen['scale']);checks=0
        for mat,key in [(H,'H'),(Mass,'M')]:
            for i,row in enumerate(mat):
                for j,v in enumerate(row):
                    p=[F(x,oldS) for x in bk[key][i][j]]
                    assert overlap(v,p),(k,key,i,j);checks+=1
        x={int(j):F(v) for j,v in z['low_trial_coordinates'].items()}
        qrat={int(j):F(v) for j,v in z['corrected_tail_coordinates'].items()}
        assert set(x)==set(low) and set(qrat)==set(bank)
        q={j:I(v) for j,v in qrat.items()}
        d={j:I(x[j]) for j in low}
        d[k]-=sum((I(x[j])*c[j]/c[k] for j in low if j!=k),O)
        fb=(1-I(x[k])*c[k])/B
        st=B-sum((c[j].sq() for j in low),O);assert st.l>0
        tau=sum((q[j]*c[j] for j in bank),O)/st
        t={j:d[j]+tau*c[j] for j in low};t.update(q);T=fb-tau
        b0=(eh+em)/2 if k==1 else -(eh-em)/2
        b1=-(eh-em)/4 if k==1 else (eh+em)/4
        ep=(b0+2*b1)/2;en=(b0-2*b1)/2;Cb=2*b0/rt
        assert Cb.up()<3
        Pb={}
        for j in c:
            v=nu[j];val=O
            for nn,w,(si,co) in zip(g['NS'],W,trig[j]):
                sn=I(nn).sqrt();bt=ep*sn+en/sn
                val-=2*w/(rt*(v.sq()+qt))*(v*(bt+b0*co)-b1*si)
            Pb[j]=val
        # Raw low/high sine contractions independently formed from final
        # source f and q, not the saved old energy plus a correction.
        wh={j:Hb[j]-sum((c[i]*hss(i,j) for i in low),O) for j in low}
        wp={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),O) for j in t}
        hb=Hbb-sum((c[i]*Hb[i] for i in low),O)
        phb=Pbb-sum((c[i]*Pb[i] for i in low),O)
        weighted=I(cap)*st+sum((I(diagonal[j]-cap)*c[j].sq() for j in high),O)
        alpha=(sum((d[j]*(wh[j]-wp[j]) for j in low),O)+fb*(hb-phb)
               +sum((t[j]*wp[j] for j in t),O)+T*phb
               +sum((q[j]*I(diagonal[j])*c[j] for j in bank),O)-tau*weighted)/st
        # Evaluate actual e coefficients on Q85; all finite prefix rows.
        e={}
        for j in high:
            e[j]=(sum((d[i]*(hss(i,j)-pss(i,j)) for i in low),O)
                  +sum((t[i]*pss(i,j) for i in t),O)
                  +fb*(Hb[j]-Pb[j])+T*Pb[j]+I(diagonal[j])*(q.get(j,O)-tau*c[j])-alpha*c[j])
            if j%200 in (k%2,):print('direct prefix',k,j,flush=True)
        # Low covectors W(e_i,f+q), W(b,f+q) contain full signed form
        # without its pole term. The pole term vanishes on low tests.
        raw={i:sum((t[j]*hss(i,j) for j in t),O)+T*Hb[i] for i in low}
        rawb=sum((t[j]*Hb[j] for j in t),O)+T*Hbb
        s=[raw[i]-c[i]/c[k]*raw[k] if i!=k else raw[k]-c[k]/B*rawb for i in low]
        # Euler's identity for the joint quadratic at (x,q,1), where
        # normalized pole derivative is rawb/B plus the exact +/-2.
        energy=sum((I(x[i])*v for i,v in zip(low,s)),O)+sum((q[j]*e[j] for j in bank),O)+rawb/B+(2 if k==1 else -2)
        phh=phb-sum((c[i]*wp[i] for i in low),O)
        def block_trial(xx,qq):
            qqi={j:I(qq.get(j,F(0))) for j in bank}
            dd={j:I(xx[j]) for j in low}
            dd[k]-=sum((I(xx[j])*c[j]/c[k] for j in low if j!=k),O)
            ff=(1-I(xx[k])*c[k])/B
            tt=sum((qqi[j]*c[j] for j in bank),O)/st
            fr={j:sum((dd[i]*hss(i,j) for i in low),O)+ff*Hb[j] for j in low+bank}
            frb=sum((dd[i]*Hb[i] for i in low),O)+ff*Hbb
            frh=sum((dd[i]*wh[i] for i in low),O)+ff*hb
            rj={j:fr[j]-c[j]*frh/st for j in bank}
            cq_h=sum((qqi[j]*(wp[j]+I(diagonal[j])*c[j]) for j in bank),O)-tt*(phh+weighted)
            ce={j:sum((qqi[i]*pss(i,j) for i in bank),O)-tt*wp[j]+I(diagonal[j])*(qqi[j]-tt*c[j])-c[j]*cq_h/st for j in bank}
            ee={j:rj[j]+ce[j] for j in bank}
            sr={j:fr[j]+sum((qqi[i]*hss(i,j) for i in bank),O)-tt*wh[j] for j in low}
            sb=frb+sum((qqi[i]*Hb[i] for i in bank),O)-tt*hb
            ss=[sr[i]-c[i]*sr[k]/c[k] if i!=k else sr[k]-c[k]*sb/B for i in low]
            Wff=sum((dd[i]*fr[i] for i in low),O)+ff*frb+(2 if k==1 else -2)
            mixed=sum((qqi[j]*rj[j] for j in bank),O)
            Cqq=sum((qqi[j]*ce[j] for j in bank),O)
            return Wff+2*mixed+Cqq,ss,ee,[Wff,2*mixed,Cqq]
        E1,s1,e1,parts=block_trial(x,qrat)
        assert overlap(E1,pair(energy)),(k,'new block vs Euler')
        assert all(overlap(v,pair(w)) for v,w in zip(s,s1)),(k,'new low block vs Euler')
        assert all(overlap(e[j],pair(e1[j])) for j in bank),(k,'new high block vs Euler')
        prior=next(v['joint_trial'] for v in N['sectors'][ks]['rows'] if v['cutoff']==365)
        coords=list(map(F,prior['trial_rational_coordinates']));ids=prior['trial_order']
        assert ids[:len(low)]==low and ids[-1]=='normalized_pole' and coords[-1]==1
        x0=dict(zip(low,coords[:len(low)]))
        q0=dict(zip(ids[len(low):-1],coords[len(low):-1]))
        assert set(q0)<=set(bank)
        E0,s0,e0,parts0=block_trial(x0,q0)
        change=sum((I(x[j]-x0[j])*(u+v) for j,u,v in zip(low,s0,s1)),O)+sum((I(qrat[j]-q0.get(j,F(0)))*(e0[j]+e1[j]) for j in bank),O)
        assert overlap(E0+change,pair(E1)),(k,'new energy polarization')
        def width_check(v,old,label):
            print(k,label,'interval',float(F(v.l,S)),float(F(v.h,S)),
                  'width',float(F(v.h-v.l,S)),'saved',float(F(old[0])),float(F(old[1])),flush=True)
            assert v.l>0,(k,label,'nonpositive lower endpoint',pair(v))
            assert (v.h-v.l)*10**8<=v.l,(k,label,'excessive relative width',pair(v))
            assert overlap(v,old),(k,label,'disjoint intervals',pair(v))
        width_check(energy,z['joint_energy_interval'],'Euler final')
        width_check(E1,z['joint_energy_interval'],'block final')
        width_check(E0,prior['joint_trial_energy'],'block initial')
        width_check(E0+change,z['joint_energy_interval'],'polarized final')
        # Both are valid independent enclosures of the same final quadratic.
        # Intersect only after separately checking widths and overlap.
        energy=I(max(energy.l,E1.l,(E0+change).l),min(energy.h,E1.h,(E0+change).h),True)
        E=F(energy.l,S)
        vraw=read(archive/'cr_alp_049'/f'result_{k}.json',{1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k])
        vs=int(vraw['scale']);V=[[I(F(int(l),vs),F(int(h),vs)) for l,h in row] for row in vraw['V']]
        n=len(low);assert len(V)==n and all(len(row)==n for row in V)
        cond=[sum((V[i][j]*s[i] for i in range(n)),O) for j in range(n)]
        assert all(overlap(v,p) for v,p in zip(cond,z['conditioned_joint_residual_vector']))
        Li={j:d.get(j,O)*AG[j]+t[j]*AP[j] for j in t}
        H0=-2/a*sum(Li.values(),O)-(tau*I(cap)+alpha)*Cb
        for v,key in [(fb,'source_pole'),(T,'combined_pole'),(alpha,'alpha'),(H0,'leading_combined_moment')]:
            assert overlap(v,z[key]),(k,key)
        # Independently evaluate the signed 32-order moments and use the
        # complete geometric absolute remainder, on the actual source trial.
        rho0=(max(nu[j].up() for j in t)/F(15*J,11))**2;assert rho0<1
        mom={};charges={}
        def add(p,v):charges[p]=charges.get(p,F(0))+v
        def series(label,coeff,mult,shift):
            terms=coeff.copy();rr=[]
            for l in range(32):
                v=sum(terms.values(),O);assert overlap(v,z['signed_moments'][label][l]),(k,label,l)
                rr.append(pair(v));add(2*l+2-shift,mult*v.up())
                terms={j:v*nu[j].sq() for j,v in terms.items()}
            add(66-shift,mult*sum(v.up() for v in terms.values())/(1-rho0));mom[label]=rr
        series('centered_combined',{j:v*(nu[j].sq()+qt) for j,v in Li.items()},F(20,11),-1)
        series('gamma_sine',{j:v*nu[j] for j,v in d.items()},F(80,11),0)
        series('prime_sine',{j:v*nu[j] for j,v in t.items()},F(80,11),0)
        J0=J-1
        tail_series=sum(v*root(F(J0,2*(2*p-1)))/F(15*J0,11)**p for p,v in charges.items())
        pole,prec=leading(g,k,J,fb,T,H0);tail=tail_series+pole
        prefix=sum(F(e[j].sq().h,S)/(diagonal[j]-F(21,8)) for j in high)
        eb=sum(F(e[j].sq().h,S) for j in bank)
        U=prefix+tail*tail/dc;assert E>U>0
        # AO direct finite charge, using newly evaluated signed covectors
        # and bank norm; frozen response coordinates and whitening only.
        perm=ac['coordinate_order'];assert sorted(perm)==list(range(n))
        LL=[[F(v) for v in row] for row in ac['LDL_lower']];w=list(map(F,ak['frame_roots']))
        qg=int(ac['response_coordinate_grid']);Q=ac['response_transformed_coordinates']
        assert len(Q)==len(bank) and all(len(row)==n for row in Q)
        bf=F(0)
        for j in range(n):
            v=sum((cond[perm[i]]*I(LL[i][j]) for i in range(n)),O)
            col2=F(sum(int(row[j])**2 for row in Q),qg*qg)
            bf+=(w[j]*(v.up()+root(eb*col2)))**2
        bg=int(al['whitening_grid']);frob2=F(sum(int(v)**2 for row in al['whitening_matrix'] for v in row),bg*bg)
        rho=F(AO['sectors'][ks]['rho_upper']);assert 0<rho<1
        bnew=frob2*bf/(1-rho);eta=rho/(1-rho)
        loss=up((U+(root(bnew)+root(eta*U))**2)/E)
        rec={'dimension':n,'grid':str(S),'finite_form_and_mass_overlaps':checks,
             'low_trial_coordinates':z['low_trial_coordinates'],'corrected_tail_coordinates':z['corrected_tail_coordinates'],
             'energy_interval':pair(energy),'energy_relative_width':str(F(energy.h-energy.l,energy.l)),
             'block_initial_energy_interval':pair(E0),'block_final_energy_interval':pair(E1),
             'polarization_change_interval':pair(change),'polarized_final_energy_interval':pair(E0+change),
             'conditioned_residual_intervals':[pair(v) for v in cond],
             'prefix_intervals':[[j,*pair(e[j])] for j in high],
             'source_pole':pair(fb),'combined_pole':pair(T),'alpha':pair(alpha),'leading_combined_moment':pair(H0),
             'signed_moments':mom,'pole_receipt':prec,'series_tail_upper':str(tail_series),
             'weighted_prefix_upper':str(prefix),'bank_defect_norm2_upper':str(eb),'U':str(U),'E':str(E),
             'rho_upper':str(rho),'U_over_E_upper':str(up(U/E)),'finite_charge_over_E_upper':str(up(bnew/E,10**30)),
             'loss_over_E_upper':str(loss),'classification':'FROZEN_TRIAL_DIRECT_REPLAY_PASS' if loss<1 else 'DIRECT_REPLAY_INCONCLUSIVE'}
        ans[ks]=rec
        payload=json.dumps({'status':'CONDITIONAL_FROZEN_TRIAL_DIRECT_REPLAY','sectors':ans,'scope':'Independent trial reconstruction from original 046 formulas. AC/AK response Gram and analytic comparator premises still inherited.'},indent=2)+chr(10)
        assert json.loads(payload)['sectors'].keys()==ans.keys()
        Path(out).write_text(payload)
        print(k,'E',float(E),'U/E',float(U/E),'loss/E',loss,rec['classification'],flush=True)
if __name__=='__main__':main(*sys.argv[1:])
