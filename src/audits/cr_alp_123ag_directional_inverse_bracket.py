#!/usr/bin/env python3
"""Frozen 123AG: regenerate only signed auxiliary residual prefix, bound the specific inverse action.
No response solves, no changes to trial or bank. Inherited classical interval premises.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
sys.set_int_max_str_digits(0)
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
J=2048
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p,h):
    assert digest(p)==h
    return json.loads(p.read_bytes())
def root(v):
    assert v>=0
    D=10**110;z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)
def namespace():
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()=='554e7ab54da4545987d89aa660a4abe6c951cb63'
    src=raw.decode()
    def replace(a,b):
        nonlocal src
        assert src.count(a)==1,a
        src=src.replace(a,b)
    # Identical 123AC / 123W wide-angle and high-prefix source patches.
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
    replace('range(1,N+1)};trig=', 'range(1,2049)};trig=')
    replace('for j in ids};Hb={}\n for j in ids:', 'for j in range(k,2049,2)};Hb={}\n for j in range(k,2049,2):')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    replace(old,old+', Bmass, Hbb, c, Hb, Pbb')
    g={'__file__':str(SRC),'__name__':'cr_alp_123ag_source'}
    exec(compile(src,str(SRC),'exec'),g)
    return g
def main(src,ac,out):
    prior=read(src,'7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897')
    data=read(ac,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    assert prior['status']=='CONDITIONAL_CONDITIONED_JOINT_CORRECTION'
    assert data['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER'
    assert set(prior['sectors'])==set(data['sectors'])=={'1','2'}
    g=namespace()
    I,S,ZERO=(g[x] for x in ('I','S','ZERO'))
    a,L,nu,TS,W,trig,rt,eh,em,qt=(g[x] for x in ('A','L','nu','TS','W','trig','rt','eh','em','qt'))
    assert len(nu)==J and len(g['AI'])==J
    bands=json.loads(SRC.with_name('CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a0['end']+1 for a0,b in zip(bands,bands[1:]))
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    cap=next(dg(b) for b in bands if b['start']==2049)
    assert cap==F('4.447938456')
    assert all(dg(b)>=F(72,25) for b in bands) and all(dg(b)>=cap for b in bands if b['start']>=2049)
    def diag(j):return min(cap,dg(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end']))))
    dc=cap-F(21,8); M=cap+F(21,8); t0=1/M
    def dj(j):return diag(j)-F(21,8)
    AP={j:sum((w*t[0] for w,t in zip(W,trig[j])),ZERO) for j in nu}
    AG={j:g['AI'][j]-AP[j] for j in nu}
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*ct+st/(a*u)) for w,t,(st,ct) in zip(W,TS,trig[i])),ZERO)
        return -2/a*(v*AP[i]-u*AP[j])/(v.sq()-u.sq())
    hss=g['Hss']
    def overlap(x,p):return max(F(x.l,S),F(p[0]))<=min(F(x.h,S),F(p[1]))
    def pair(x):return [str(F(x.l,S)),str(F(x.h,S))]
    def upgrid(v,G=10**12):return str(F(-((-v*G).__floor__()),G))
    def logratio(x,E):return float(x/E)
    ans={}
    for k in (1,2):
        ks=str(k); rows=prior['sectors'][ks];assert len(rows)==1
        r=rows[0]; z=data['sectors'][ks];n=44-k
        assert r['cutoff']==850 and r['source_cutoff']==365 and r['finite_prefix']==J
        low=list(range(k,86,2));high=list(range(87 if k==1 else 86,J+1,2));bank=[j for j in high if j<=850]
        assert z['bank_indices']==bank and [j for j,_ in z['prefix_mid_radius']]==high
        assert len(z['tail_receipts'])==len(z['LDL_pivots'])==n
        assert len(low)==n and len(r['conditioned_joint_residual_vector'])==n
        _,_,_,B,Hbb,c,Hb,Pbb=g['build'](k)
        b0=(eh+em)/2 if k==1 else -(eh-em)/2
        b1=-(eh-em)/4 if k==1 else (eh+em)/4
        ep=(b0+2*b1)/2;en=(b0-2*b1)/2;Cb=2*b0/rt
        def pb(j):
            v=nu[j];val=ZERO
            for nn,w,(st,ct) in zip(g['NS'],W,trig[j]):
                sn=I(nn).sqrt();bt=ep*sn+en/sn
                val-=2*w/(rt*(v.sq()+qt))*(v*(bt+b0*ct)-b1*st)
            return val
        Pb={j:pb(j) for j in c}
        st=B-sum((c[j].sq() for j in low),ZERO);assert st.l>0
        weighted=I(cap)*st+sum((I(diag(j)-cap)*c[j].sq() for j in high),ZERO)
        active=low+bank
        hs={j:Hb[j]-sum((c[i]*hss(i,j) for i in low),ZERO) for j in active}
        ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in active}
        hb=Hbb-sum((c[i]*Hb[i] for i in low),ZERO)
        phb=Pbb-sum((c[i]*Pb[i] for i in low),ZERO)
        x={int(j):F(v) for j,v in r['low_trial_coordinates'].items()}
        q={int(j):I(F(v)) for j,v in r['corrected_tail_coordinates'].items()}
        assert set(x)==set(low) and set(q)==set(bank)
        d={j:I(v) for j,v in x.items()}
        d[k]-=sum((I(x[j])*c[j]/c[k] for j in low if j!=k),ZERO)
        fb=(1-I(x[k])*c[k])/B
        tau=sum((q[j]*c[j] for j in bank),ZERO)/st
        t={j:d[j]+tau*c[j] for j in low};t.update(q);T=fb-tau
        alpha=(sum((d[j]*(hs[j]-ps[j]) for j in d),ZERO)+fb*(hb-phb)
               +sum((t[j]*ps[j] for j in t),ZERO)+T*phb
               +sum((q[j]*I(diag(j))*c[j] for j in q),ZERO)-tau*weighted)/st
        H0=-2/a*sum((d.get(j,ZERO)*AG[j]+t[j]*AP[j] for j in t),ZERO)-(tau*I(cap)+alpha)*Cb
        assert all(overlap(v,r[name]) for v,name in [(fb,'source_pole'),(T,'combined_pole'),(alpha,'alpha'),(H0,'leading_combined_moment')])
        e={}
        for j in high:
            e[j]=(sum((d[i]*(hss(i,j)-pss(i,j)) for i in d),ZERO)
                  +sum((t[i]*pss(i,j) for i in t),ZERO)
                  +fb*(Hb[j]-Pb[j])+T*Pb[j]+I(diag(j))*(q.get(j,ZERO)-tau*c[j])-alpha*c[j])
            if j%200 in (k%2,):print('signed frozen residual',k,j,flush=True)
        normlo=sum(F(e[j].sq().l,S) for j in high)
        normhi=sum(F(e[j].sq().h,S) for j in high)
        pwlo=sum(F(e[j].sq().l,S)/dj(j) for j in high)
        pwhi=sum(F(e[j].sq().h,S)/dj(j) for j in high)
        assert pwlo<=F(r['weighted_prefix_upper'])
        assert sum(F(e[j].sq().l,S) for j in bank)<=F(r['bank_defect_norm2_upper'])
        # Exact pole-null pairing gives finite weighted pair and complete denominator.
        # B0=Q diag(dj) Q, t0 I <= C^-1 <= B0^-1.
        diff={j:1/dj(j)-1/dc for j in high}
        w={j:(t0+1/dj(j))/2 for j in high}
        wc=(t0+1/dc)/2
        assert all(v>=0 for v in diff.values())
        h=st/I(dc)+sum((diff[j]*c[j].sq() for j in high),ZERO)
        assert h.l>0
        ve=sum((diff[j]*c[j]*e[j] for j in high),ZERO)
        PG=int(z['prefix_grid']);assert PG==10**60
        piv=list(map(F,z['LDL_pivots']));assert min(piv)>0
        assert all(len(row)==n and all(type(m)==int and type(rad)==int and rad>=0 for m,rad in row) for j,row in z['prefix_mid_radius'])
        te=F(r['infinite_tail_norm_upper']);assert te>=0
        centers=[];center_sq=F(0);center_prefix_sq=F(0);tail_box_sq=F(0)
        for i,p in enumerate(piv):
            zz={j:I(F(row[i][0]-row[i][1],PG),F(row[i][0]+row[i][1],PG)) for j,row in z['prefix_mid_radius']}
            vz=sum((diff[j]*c[j]*zz[j] for j in high),ZERO)
            mid=sum((w[j]*zz[j]*e[j] for j in high),ZERO)-vz*ve/(2*h)
            tz=F(z['tail_receipts'][i]['tail_norm_upper']);assert tz>=0
            reserve=wc*tz*te
            center_sq+=p*(mid.up()+reserve)**2
            center_prefix_sq+=p*mid.up()**2;tail_box_sq+=p*reserve**2
            centers.append({'prefix_and_pole_center_interval':pair(mid),'tail_pair_abs_upper':str(reserve)})
        E,U,b2=map(F,(z['E'],z['U'],z['finite_charge_upper']))
        assert E==F(r['joint_energy_lower']) and U==F(r['weighted_inverse_error_upper'])
        assert E>0 and U>0 and b2>=0
        # 123AC trace upper also bounds the B0^-1 response Gram.
        eta=F(z['eta_upper']); assert eta>=0 and U-t0*normlo>=0
        # K=(B0^-1+t0 I)/2; R0=(B0^-1-t0 I)/2; -R0<=C^-1-K<=R0.
        # ||A*(C^-1-K)e|| <= sqrt(||A*R0 A|| <e,R0 e>).
        bracket_err_sq=eta*(U-t0*normlo)/4
        dirupper=min(root(eta*U),root(center_sq)+root(bracket_err_sq))
        loss=U+(root(b2)+dirupper)**2
        loss=min(loss,F(z['loss_upper']))
        rec={'dimension':n,'prefix_grid':str(S),'signed_prefix_intervals':[[j,str(e[j].l),str(e[j].h)] for j in high],
             'source_overlap_checks':4,'weighted_prefix_upper_ratio':upgrid(pwhi/E),
             'prefix_norm2_lower_ratio':str(F((normlo/E*10**12).__floor__(),10**12)),
             'prefix_norm2_upper_ratio':upgrid(normhi/E),
             'pole_pair_interval':pair(ve),'pole_denominator_interval':pair(h),'directional_centers':centers,
             'center_prefix_norm2_over_E_upper':upgrid(center_prefix_sq/E),
             'center_tail_box_norm2_over_E_upper':upgrid(tail_box_sq/E),
             'center_norm2_over_E_upper':upgrid(center_sq/E),
             'inverse_bracket_error2_over_E_upper':upgrid(bracket_err_sq/E),
             'directional_charge2_over_E_upper':upgrid(dirupper**2/E),
             'old_scalar_charge2_over_E_upper':upgrid(eta*U/E),
             'loss_over_E_upper':upgrid(loss/E),
             'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE'}
        ans[ks]=rec
        print(k,'prefix/U',float(pwhi/U),'center/E',float(center_sq/E),'bracket/E',float(bracket_err_sq/E),'loss/E',float(loss/E),rec['classification'],flush=True)
        out.write_text(json.dumps({'status':'CONDITIONAL_DIRECTIONAL_INVERSE_BRACKET','frozen_response_sha256':digest(ac),'frozen_pole_trial_sha256':digest(src),'sectors':ans,'scope':'Frozen comparator, interval and signed tail premises. No new response solutions, no pole sign without strict sufficient test.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
