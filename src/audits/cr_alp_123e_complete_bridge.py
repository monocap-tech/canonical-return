#!/usr/bin/env python3
"""Complete Q85 projection, including the infinite pole bridge, at E125/E165.

The capped comparator is exact on the chosen finite projected subspace, including
its infinite pole tail. Positive finite Schur values are not infinite lower bounds.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,sys
from cr_alp_123a_pole_schur import endpoints,schur,enc
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
PIN='554e7ab54da4545987d89aa660a4abe6c951cb63'

def main(root,previous,out):
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==PIN
    bankraw=(root/'cr_alp_046/CR-ALP-046_full_E85_intervals.json').read_bytes()
    assert hashlib.sha256(bankraw).hexdigest()=='d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22'
    frozen=json.loads(bankraw)
    src=raw.decode()
    assert src.count('M=90;SHIFT=512;N=85')==1
    src=src.replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=165')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    assert src.count(old)==1
    src=src.replace(old,old+', Bmass, Hbb, c, Hb, Pbb')
    g={'__file__':str(SRC),'__name__':'cr_alp_123e_bridge'}
    exec(compile(src,str(SRC),'exec'),g)
    S,I,ZERO=g['S'],g['I'],g['ZERO'];assert S==frozen['scale']
    prior=json.loads(previous.read_text());assert prior['status']=='CONDITIONAL_POLE_RESIDUAL_CENSUS'
    bands=json.loads(SRC.with_name('CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    cap=next(dg(b) for b in bands if b['start']==2049)
    assert all(dg(b)>=F(72,25) for b in bands) and all(dg(b)>=cap for b in bands if b['start']>=2049)
    def diag(i):return min(cap,dg(next(b for b in bands if b['start']<=i and (b['end'] is None or i<=b['end']))))
    a,L,nu,TS,W,trig,rt,eh,em,qt=(g[x] for x in ('A','L','nu','TS','W','trig','rt','eh','em','qt'))
    theta={j:sum((w*t[0] for w,t in zip(W,trig[j])),ZERO) for j in nu}
    @lru_cache(None)
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*ct+st/(a*u)) for w,t,(st,ct) in zip(W,TS,trig[i])),ZERO)
        return -2/a*(v*theta[i]-u*theta[j])/(v.sq()-u.sq())
    hss=lru_cache(None)(g['Hss'])
    result={}
    for oldbank in frozen['blocks']:
        k=oldbank['parity'];order,H,Mass,B,Hbb,c,Hb,Pbb=g['build'](k)
        low=oldbank['ids'];ix=[order.index(j) for j in low]
        for mat,key in [(H,'H'),(Mass,'M')]:
            for i,u in enumerate(ix):
                for j,v in enumerate(ix):
                    assert [max(mat[u][v].l,mat[v][u].l),min(mat[u][v].h,mat[v][u].h)]==oldbank[key][i][j]
        high=[j for j in order if j>85];nl=len(low);nh=len(high)
        b0=(eh+em)/2 if k==1 else -(eh-em)/2
        b1=-(eh-em)/4 if k==1 else (eh+em)/4
        ep=(b0+2*b1)/2;en=(b0-2*b1)/2
        def pb(j):
            v=nu[j];den=v.sq()+qt;val=ZERO
            for nn,w,(st,ct) in zip(g['NS'],W,trig[j]):
                sn=I(nn).sqrt();bt=ep*sn+en/sn
                val-=2*w/(rt*den)*(v*(bt+b0*ct)-b1*st)
            return val
        Pb={j:pb(j) for j in order}
        st=B-sum((c[j].sq() for j in low),ZERO);assert st.l>0
        # h=(I-S85)b. Compute its full capped diagonal energy by Parseval;
        # beyond 2048 every coefficient carries exactly the same cap.
        def cj(j):
            v=j*g['NUF'];return 2*b0*v/(rt*(v.sq()+qt))
        weighted=I(cap)*st
        for j in range(87 if k==1 else 86,2049,2):
            weighted+=I(diag(j)-cap)*cj(j).sq()
        def full_h_data(ss,sb,bb):
            hs={j:sb[j]-sum((c[i]*ss(i,j) for i in low),ZERO) for j in order}
            hb=bb-sum((c[i]*sb[i] for i in low),ZERO)
            hh=hb-sum((c[i]*hs[i] for i in low),ZERO)
            return hs,hb,hh
        hs,hb,hh=full_h_data(hss,Hb,Hbb)
        ps,phb,phh=full_h_data(pss,Pb,Pbb)
        # Low vectors as in 046, high vectors Q85 e_j=e_j-c_j h/st.
        # hp=b/B. Use structured bilinear formulas to avoid a dense repeated
        # expansion of all complete-pole coefficients.
        lowv=[({j:I(1),k:-c[j]/c[k]},ZERO) if j!=k else ({k:I(1)},-c[k]/B) for j in low]
        def on_e(v,j):
            d,b=v;return sum((x*hss(i,j) for i,x in d.items()),ZERO)+b*Hb[j]
        def on_b(v):
            d,b=v;return sum((x*Hb[i] for i,x in d.items()),ZERO)+b*Hbb
        def on_h(v):
            d,b=v;return sum((x*hs[i] for i,x in d.items()),ZERO)+b*hb
        n=nl+nh+1;base=[[ZERO]*n for _ in range(n)]
        for i,u in enumerate(ix):
            for j,v in enumerate(ix):base[i][j]=H[u][v]
        C=[[ZERO]*nh for _ in range(nh)]
        for i,v in enumerate(lowv):
            base[i][-1]=base[-1][i]=on_b(v)/B
            hv=on_h(v)
            for j,p in enumerate(high):
                base[i][nl+j]=base[nl+j][i]=on_e(v,p)-c[p]/st*hv
        base[-1][-1]=Hbb/B.sq()+(2 if k==1 else -2)
        for i,p in enumerate(high):
            cp=c[p]/st
            base[nl+i][-1]=base[-1][nl+i]=(Hb[p]-cp*hb)/B
            for j in range(i,nh):
                q=high[j];cq=c[q]/st
                base[nl+i][nl+j]=base[nl+j][nl+i]=hss(p,q)-cp*hs[q]-cq*hs[p]+cp*cq*hh
                d=(I(diag(p)) if p==q else ZERO)-cp*I(diag(q))*c[q]-cq*I(diag(p))*c[p]+cp*cq*weighted
                C[i][j]=C[j][i]=pss(p,q)-cp*ps[q]-cq*ps[p]+cp*cq*phh+d
        rows=[]
        for N in (125,165):
            keep=[i for i,j in enumerate(high) if j<=N]
            J=list(range(nl))+[nl+i for i in keep]+[n-1]
            for t in (F(0),F(1)):
                A=[[base[i][j] for j in J] for i in J]
                for i,p in enumerate(keep):
                    for j,q in enumerate(keep):A[nl+i][nl+j]=(1-t)*C[p][q]+t*base[nl+p][nl+q]
                lo=schur(endpoints(A,-1,S));hi=schur(endpoints(A,1,S));assert lo<=hi
                if N==125 and t==1:
                    l,h=map(F,prior['sectors'][str(k)]['expanded_finite_schur_interval'])
                    assert max(l,lo)<=min(h,hi),'Original E125 span consistency failed'
                label='NEGATIVE_COMPARATOR_WITNESS' if hi<0 and t==0 else 'NEGATIVE_ORIGINAL_WITNESS' if hi<0 else 'POSITIVE_FINITE_ONLY' if lo>0 else 'UNDECIDED'
                row={'cutoff':N,'theta':str(t),'dimension':len(J),'schur_interval':[enc(lo)[0],enc(hi)[1]],'classification':label}
                rows.append(row);print('parity',k,'cutoff',N,'theta',t,'Schur',float(lo),float(hi),label,flush=True)
        result[str(k)]={'rows':rows,'complete_bridge':True,'tail_pole_mass_interval':[str(F(st.l,S)),str(F(st.h,S))],'weighted_complete_tail_pole_interval':[str(F(weighted.l,S)),str(F(weighted.h,S))]}
    out.write_text(json.dumps({'status':'CONDITIONAL_COMPLETE_BRIDGE_COMPARATOR_CENSUS','sectors':result,
      'scope':'Finite projected subspaces of actual Q85 including complete infinite pole bridge, exact capped diagonal tail by Parseval. No infinite positivity claim. Original continuum and interval premises inherited.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
