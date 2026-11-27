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
        assert overlap(energy,z['joint_energy_interval']) and energy.l>0,(k,'energy')
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
             'energy_interval':pair(energy),'conditioned_residual_intervals':[pair(v) for v in cond],
             'prefix_intervals':[[j,*pair(e[j])] for j in high],
             'source_pole':pair(fb),'combined_pole':pair(T),'alpha':pair(alpha),'leading_combined_moment':pair(H0),
             'signed_moments':mom,'pole_receipt':prec,'series_tail_upper':str(tail_series),
             'weighted_prefix_upper':str(prefix),'bank_defect_norm2_upper':str(eb),'U':str(U),'E':str(E),
             'rho_upper':str(rho),'U_over_E_upper':str(up(U/E)),'finite_charge_over_E_upper':str(up(bnew/E,10**30)),
             'loss_over_E_upper':str(loss),'classification':'FROZEN_TRIAL_DIRECT_REPLAY_PASS' if loss<1 else 'DIRECT_REPLAY_INCONCLUSIVE'}
        ans[ks]=rec
        Path(out).write_text(json.dumps({'status':'CONDITIONAL_FROZEN_TRIAL_DIRECT_REPLAY','sectors':ans,'scope':'Independent trial reconstruction from original 046 formulas. AC/AK response Gram and analytic comparator premises still inherited.'},indent=2)+'\n')
        print(k,'E',float(E),'U/E',float(U/E),'loss/E',loss,rec['classification'],flush=True)
if __name__=='__main__':main(*sys.argv[1:])
