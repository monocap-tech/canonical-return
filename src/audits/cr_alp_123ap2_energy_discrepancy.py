#!/usr/bin/env python3
"""123AP: direct frozen Z trial audit, independent of historical energy updates.
Changed exact interval grid and digamma shift; same classical 046 identities.
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
    replace('S=10**340; a=F(11,10)','S=10**180; a=F(11,10)')
    replace('M=90;SHIFT=512;N=85','M=90;SHIFT=768;N=85')
    # Preserve original analytic remainder tests, use a deliberately wider
    # absolute digamma reserve than the certified 10^-290.
    replace('return re.widen(10**50),im.widen(10**50),tr.widen(10**50)','return re.widen(10**10),im.widen(10**10),tr.widen(10**10)')
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
        # Distinguish nonpositive energy from disjoint intervals. No gate
        # relaxation: this executable produces DIAGNOSTIC_ONLY records.
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
        inherited=I(*map(F,z['joint_energy_interval']))
        inherited0=I(*map(F,prior['joint_trial_energy']))
        vraw=read(archive/'cr_alp_049'/f'result_{k}.json',{1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k])
        vs=int(vraw['scale']);V=[[I(F(int(l),vs),F(int(h),vs)) for l,h in row] for row in vraw['V']]
        n=len(low);cond=[sum((V[i][j]*s[i] for i in range(n)),O) for j in range(n)]
        olde=F(z['bank_defect_norm2_upper']);newe=sum(F(e[j].sq().h,S) for j in bank)
        low_checks=[overlap(v,p) for v,p in zip(cond,z['conditioned_joint_residual_vector'])]
        oldprefix=F(z['weighted_prefix_upper']);newprefix=sum(F(e[j].sq().h,S)/(diagonal[j]-F(21,8)) for j in high)
        pole_checks={key:overlap(v,z[key]) for v,key in [(fb,'source_pole'),(T,'combined_pole'),(alpha,'alpha')]}
        def summary(v):return {'interval':pair(v),'approx':[float(F(v.l,S)),float(F(v.h,S))]}
        rec={'classification':'DIAGNOSTIC_ONLY','dimension':n,'finite_form_and_mass_overlaps':checks,
             'direct_final_energy':summary(E1),'saved_final_energy':summary(inherited),'final_difference':summary(E1-inherited),
             'direct_initial_energy':summary(E0),'saved_initial_energy':summary(inherited0),'initial_difference':summary(E0-inherited0),
             'final_energy_overlap':overlap(E1,z['joint_energy_interval']),'initial_energy_overlap':overlap(E0,prior['joint_trial_energy']),
             'direct_final_positive':E1.l>0,'direct_initial_positive':E0.l>0,
             'direct_final_energy_terms':[summary(v) for v in parts],'direct_initial_energy_terms':[summary(v) for v in parts0],
             'polarization_change':summary(change),'saved_initial_plus_direct_change_difference':summary(inherited0+change-inherited),
             'conditioned_residual_overlap_count':sum(low_checks),'conditioned_residual_dimension':n,
             'pole_coordinate_overlaps':pole_checks,'old_bank_defect_norm2_upper':str(olde),'new_bank_defect_norm2_upper':str(newe),
             'old_weighted_prefix_upper':str(oldprefix),'new_weighted_prefix_upper':str(newprefix)}
        ans[ks]=rec
        Path(out).write_text(json.dumps({'status':'DIRECT_TRIAL_ENERGY_DISCREPANCY_DIAGNOSTIC_ONLY','sectors':ans,'scope':'No positivity acceptance. Compare original initial/final energy, direct block and Euler evaluations.'},indent=2)+'\\n')
        print(k,'initial energy',float(E0.mid()),'saved',float(inherited0.mid()),'overlap',rec['initial_energy_overlap'], 'final energy',float(E1.mid()),'saved',float(inherited.mid()),'overlap',rec['final_energy_overlap'],'low overlaps',sum(low_checks),'bank ratio',float(newe/olde),'prefix ratio',float(newprefix/oldprefix),flush=True)
if __name__=='__main__':main(*sys.argv[1:])
