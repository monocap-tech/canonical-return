#!/usr/bin/env python3
"""123AI: one frozen Galerkin candidate per parity; exact finite pairing and
sparse high-mode inverse-defect lower test. NOT a full inverse upper enclosure.
Floating proposals are untrusted; no response-column solve or bank expansion.
"""
import hashlib,json,sys
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parent
J=2048
CAP=F('4.447938456')
M=CAP+F(21,8)
NS=[2,3,4,5,7,8,9]
PS=[2,3,2,5,7,2,3]
def read(p,h):
    b=Path(p).read_bytes();assert hashlib.sha256(b).hexdigest()==h
    return json.loads(b)
def blob(p,h):
    b=p.read_bytes();assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==h
    return b
def roots(v):
    assert v>=0
    G=10**110;z=v*G*G;t=isqrt(z.numerator//z.denominator)
    return F(t,G),F(t+1,G)
def floor_ratio(v,E):
    return str(F((v/E*10**12).__floor__(),10**12))
def ceil_ratio(v,E):
    return str(F(-((-v/E*10**12).__floor__()),10**12))
def propose(k,a,g,diag):
    low=list(range(k,86,2));high=[j for j,l,h in g['signed_prefix_intervals']]
    ids=np.array(low+high);n=len(low);nu=np.pi*ids/2.2
    ts=np.log(np.array(NS));w=np.log(np.array(PS))/np.sqrt(NS)
    si=np.sin(nu[:,None]*ts);co=np.cos(nu[:,None]*ts);ap=si@w
    u=nu[:,None];v=nu[None,:];den=v*v-u*u;np.fill_diagonal(den,1)
    P=-2/1.1*(v*ap[:,None]-u*ap[None,:])/den
    np.fill_diagonal(P,-((2.2-ts)/1.1*co+si/(1.1*nu[:,None]))@w)
    eh=np.exp(.55);em=1/eh;ee=np.exp(1.1);q=1/ee
    B=(ee-q)/2+(1.1 if k==1 else -1.1)
    b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4
    ep=(b0+2*b1)/2;en=(b0-2*b1)/2;rt=np.sqrt(1.1)
    c=2*b0*nu/(rt*(nu**2+.25));bt=ep*np.sqrt(NS)+en/np.sqrt(NS)
    Pb=-2*((nu[:,None]*(bt+b0*co)-b1*si)@w)/(rt*(nu**2+.25));Pbb=0
    for nn,t,ww in zip(NS,ts,w):
        for pp,cp in [(.5,ep),(-.5,en)]:
            for rr,cr in [(.5,ep),(-.5,en)]:
                vv=pp+rr;it=2.2-t if vv==0 else ee**2/nn-1 if vv==1 else 1-nn*q*q
                Pbb-=2*ww*cp*cr*nn**rr*it
    st=B-c[:n]@c[:n];assert st>0
    ch=c[n:];cp=ch/st;ps=Pb[n:]-P[n:,:n]@c[:n]
    phh=Pbb-2*c[:n]@Pb[:n]+c[:n]@P[:n,:n]@c[:n]
    dd=np.array([float(diag(j)) for j in high]);cap=float(CAP)
    weighted=cap*st+((dd-cap)*ch)@ch
    Cb=P[n:,n:]+np.diag(dd)-np.outer(cp,ps)-np.outer(ps,cp)+phh*np.outer(cp,cp)-np.outer(cp,dd*ch)-np.outer(dd*ch,cp)+weighted*np.outer(cp,cp)
    scale=10**(19 if k==1 else 17);sg=int(g['prefix_grid'])
    er=np.array([float(F((int(l)+int(h))*scale,2*sg)) for j,l,h in g['signed_prefix_intervals']])
    sol=np.linalg.solve(Cb,er);assert np.all(np.isfinite(sol))
    wgrid=10**(34 if k==1 else 32)
    ints=[round(float(v)*wgrid/scale) for v in sol]
    rel=float(np.linalg.norm(Cb@sol-er)/np.linalg.norm(er))
    return high,wgrid,ints,rel
def namespace(samples):
    raw=blob(ROOT/'cr_alp_046_full_finite.py','554e7ab54da4545987d89aa660a4abe6c951cb63')
    src=raw.decode()
    def replace(a,b):
        nonlocal src
        assert src.count(a)==1,a
        src=src.replace(a,b)
    # 123AG wide-angle identity, same exact 046 arithmetic. Gamma only
    # evaluated where needed (low, original bank, and sparse test indices).
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
    replace('assert y.up()/F(xx.l,S)<F(1,3)','assert y.up()/F(xx.l,S)<5')
    replace('li=atan(y/xx)','li=atan_wide(y/xx)')
    replace('range(1,N+1)};trig=',f'(list(range(1,2049))+{samples!r})'+'};trig=')
    replace('for j,v in nu.items():\n re,im,tre=',f'for j,v in nu.items():\n if j>850 and j not in {samples!r}:continue\n re,im,tre=')
    # build() uses low E85 matrices and pole/sine values only on actual
    # Gamma indices; all needed prime trig values are still present.
    replace('for j in ids};Hb={}\n for j in ids:', 'for j in G if j%2==k%2};Hb={}\n for j in c:')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    replace(old,old+', Bmass, Hbb, c, Hb, Pbb')
    g={'__file__':str(ROOT/'cr_alp_046_full_finite.py'),'__name__':'cr_alp_123ai_source'}
    exec(compile(src,str(ROOT/'cr_alp_046_full_finite.py'),'exec'),g)
    return g
def sparse(k,high,wgrid,ints,r,diag,samples):
    g=namespace(samples);I,S,ZERO=(g[x] for x in ('I','S','ZERO'))
    a,L,nu,TS,W,trig,rt,eh,em,qt=(g[x] for x in ('A','L','nu','TS','W','trig','rt','eh','em','qt'))
    low=list(range(k,86,2));bank=[j for j in high if j<=850]
    _,_,_,B,Hbb,c0,Hb,Pbb=g['build'](k)
    b0=(eh+em)/2 if k==1 else -(eh-em)/2;b1=-(eh-em)/4 if k==1 else (eh+em)/4
    ep=(b0+2*b1)/2;en=(b0-2*b1)/2
    c={j:2*b0*nu[j]/(rt*(nu[j].sq()+qt)) for j in low+high+samples}
    def pb(j):
        v=nu[j];val=ZERO
        for nn,w,(si,co) in zip(NS,W,trig[j]):
            sn=I(nn).sqrt();bt=ep*sn+en/sn
            val-=2*w/(rt*(v.sq()+qt))*(v*(bt+b0*co)-b1*si)
        return val
    Pb={j:pb(j) for j in c}
    AP={j:sum((w*t[0] for w,t in zip(W,trig[j])),ZERO) for j in c}
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*co+si/(a*u)) for w,t,(si,co) in zip(W,TS,trig[i])),ZERO)
        return -2/a*(v*AP[i]-u*AP[j])/(v.sq()-u.sq())
    hss=g['Hss'];st=B-sum((c[j].sq() for j in low),ZERO);assert st.l>0
    weighted=I(CAP)*st+sum((I(diag(j)-CAP)*c[j].sq() for j in high),ZERO)
    ps={j:Pb[j]-sum((c[i]*pss(i,j) for i in low),ZERO) for j in high+samples}
    phb=Pbb-sum((c[i]*Pb[i] for i in low),ZERO)
    phh=phb-sum((c[i]*(Pb[i]-sum((c[l]*pss(l,i) for l in low),ZERO)) for i in low),ZERO)
    x={int(j):F(v) for j,v in r['low_trial_coordinates'].items()}
    q={int(j):I(F(v)) for j,v in r['corrected_tail_coordinates'].items()}
    assert set(x)==set(low) and set(q)==set(bank)
    d={j:I(v) for j,v in x.items()}
    d[k]-=sum((I(x[j])*c[j]/c[k] for j in low if j!=k),ZERO)
    fb=(1-I(x[k])*c[k])/B;tau=sum((q[j]*c[j] for j in bank),ZERO)/st
    t={j:d[j]+tau*c[j] for j in low};t.update(q);T=fb-tau
    def iv(v):return I(*map(F,v))
    def overlap(v,r):return max(v.l,r.l)<=min(v.h,r.h)
    assert overlap(fb,iv(r['source_pole'])) and overlap(T,iv(r['combined_pole']))
    # Alpha already certified on the original frozen trial. Use the frozen
    # enclosure, not a newly guessed projection constant.
    alpha=iv(r['alpha'])
    w={j:I(F(v,wgrid)) for j,v in zip(high,ints)}
    tw=sum((w[j]*c[j] for j in high),ZERO)/st
    wp=sum((w[j]*ps[j] for j in high),ZERO)
    wd=sum((w[j]*I(diag(j))*c[j] for j in high),ZERO)
    pairs=[]
    for j in samples:
        assert j>J and j%2==k%2
        e=(sum((d[i]*(hss(i,j)-pss(i,j)) for i in d),ZERO)+sum((t[i]*pss(i,j) for i in t),ZERO)
           +fb*(Hb[j]-Pb[j])+T*Pb[j]-I(CAP)*tau*c[j]-alpha*c[j])
        cp=c[j]/st
        cw=(sum((w[i]*pss(i,j) for i in high),ZERO)-cp*wp-tw*ps[j]+tw*cp*phh
            -cp*wd-tw*I(CAP)*c[j]+tw*cp*weighted)
        defect=e-cw
        pairs.append([j,str(defect.l),str(defect.h)])
    lower=sum(F(max(0,int(l),-int(h))**2,S*S) for j,l,h in pairs)/M
    return S,pairs,lower
def main(ac,ag,z,bands_path,out):
    A=read(ac,'a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473')
    G=read(ag,'45719f6cd161b29534983c95f52941d434832bfffd3493d4d421513b436c1abf')
    Z=read(z,'7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897')
    bands=json.loads(blob(Path(bands_path),'1a7fe8629daa8b8e6cec185a79b6746d6b28975d'))['bands']
    assert A['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER' and G['status']=='CONDITIONAL_DIRECTIONAL_INVERSE_BRACKET' and Z['status']=='CONDITIONAL_CONDITIONED_JOINT_CORRECTION'
    assert set(A['sectors'])==set(G['sectors'])==set(Z['sectors'])=={'1','2'}
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    assert next(dg(b) for b in bands if b['start']==2049)==CAP
    def diag(j):return min(CAP,dg(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end']))))
    result={}
    for k in (1,2):
        ks=str(k);a=A['sectors'][ks];g=G['sectors'][ks];r=Z['sectors'][ks][0]
        E,U,b2,eta=map(F,(a['E'],a['U'],a['finite_charge_upper'],a['eta_upper']))
        assert E==F(r['joint_energy_lower']) and U==F(r['weighted_inverse_error_upper']) and E>U>0 and b2>=0 and eta>0
        high,wgrid,ints,rel=propose(k,a,g,diag)
        assert high==[j for j,_ in a['prefix_mid_radius']]==list(range(87 if k==1 else 86,2049,2))
        PG=int(a['prefix_grid']);assert PG==10**60
        n=44-k;piv=list(map(F,a['LDL_pivots']));assert len(piv)==n and min(piv)>0
        lo=hi=F(0);pairings=[]
        for i,p in enumerate(piv):
            mid=sum(v*row[i][0] for v,(j,row) in zip(ints,a['prefix_mid_radius']))
            rad=sum(abs(v)*row[i][1] for v,(j,row) in zip(ints,a['prefix_mid_radius']))
            lo+=p*F(max(0,abs(mid)-rad)**2,(PG*wgrid)**2)
            hi+=p*F((abs(mid)+rad)**2,(PG*wgrid)**2)
            pairings.append([str(mid-rad),str(mid+rad)])
        # ||Dbar^(1/2)(b-Z*w)|| >= max(0, sqrt(lo)-sqrt(b2)).
        fmin=max(F(0),roots(lo)[0]-roots(b2)[1]);budget=roots(E-U)[1]
        rec={'candidate_grid':str(wgrid),'candidate_coordinates':[[j,str(v)] for j,v in zip(high,ints)],
             'pair_grid':str(PG*wgrid),'signed_pairing_intervals':pairings,'floating_galerkin_relative_residual_diagnostic':rel,
             'pair_norm2_over_E_lower':floor_ratio(lo,E),'pair_norm2_over_E_upper':ceil_ratio(hi,E),
             'finite_charge2_over_E_lower':floor_ratio(fmin*fmin,E),'available_budget_over_E_upper':ceil_ratio(E-U,E),
             'eta_policy':'frozen AC trace upper (fixed sufficient receiver, not true inverse action)'}
        if fmin>budget:
            rec['classification']='REJECT_SELECTED_CANDIDATE_FINITE_PAIRING'
        else:
            # The largest V permitted by this lower finite charge, for this
            # particular fixed eta receiver; use an outward upper budget.
            vmax=(budget-fmin)**2/eta
            rec['permitted_V_over_E_upper']=ceil_ratio(vmax,E)
            # Fixed parity samples, no adaptive selection or response expansion.
            samples=list(range(2049 if k==1 else 2050,2177,2))
            S,rows,vlo=sparse(k,high,wgrid,ints,r,diag,samples)
            rec.update({'sparse_grid':str(S),'sparse_defect_intervals':rows,'sparse_V_over_E_lower':floor_ratio(vlo,E),
                        'classification':'REJECT_SELECTED_CANDIDATE_SPARSE_DEFECT' if vlo>vmax else 'SPARSE_TEST_INCONCLUSIVE_FULL_UPPER_NOT_COMPUTED'})
        result[ks]=rec
        print(ks,rec['classification'],'pair2/E',float(lo/E), 'budget/E',float((E-U)/E),'allowed V/E',rec.get('permitted_V_over_E_upper'),'sparse V/E',rec.get('sparse_V_over_E_lower'),flush=True)
        Path(out).write_text(json.dumps({'status':'CONDITIONAL_ONE_VECTOR_CANDIDATE_DIAGNOSTIC','sectors':result,'scope':'Fixed prefix Galerkin candidate per parity. Lower tests reject sufficient certificate for a candidate only; no full inverse upper or pole sign.'},indent=2)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
