#!/usr/bin/env python3
"""Test whether the exact capped Q85 comparator can receive a pole direction.

Replace ONLY the high-high block by the historical comparator, retaining the
original low energy and original mixed interactions. A negative comparator
witness obstructs this lower-bound architecture, NOT the original Weil form.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
from cr_alp_123a_pole_schur import endpoints,schur,enc
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
PIN='554e7ab54da4545987d89aa660a4abe6c951cb63'

def main(root,out):
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==PIN
    bankraw=(root/'cr_alp_046/CR-ALP-046_full_E85_intervals.json').read_bytes()
    assert hashlib.sha256(bankraw).hexdigest()=='d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22'
    frozen=json.loads(bankraw)
    src=raw.decode()
    assert src.count('M=90;SHIFT=512;N=85')==1
    src=src.replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=125')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    assert src.count(old)==1
    src=src.replace(old,old+', Bmass, Hbb, c, Hb')
    g={'__file__':str(SRC),'__name__':'cr_alp_123d_joint'}
    exec(compile(src,str(SRC),'exec'),g)
    S,I,ZERO=g['S'],g['I'],g['ZERO'];assert S==frozen['scale']
    bands=json.loads(SRC.with_name('CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert bands[0]['start']==86 and bands[-1]['end'] is None
    assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    cap=next(dg(b) for b in bands if b['start']==2049)
    assert all(dg(b)>=F(72,25) for b in bands) and all(dg(b)>=cap for b in bands if b['start']>=2049)
    def diag(i):return min(cap,dg(next(b for b in bands if b['start']<=i and (b['end'] is None or i<=b['end']))))
    a,L,nu,TS,W,trig=(g[x] for x in ('A','L','nu','TS','W','trig'))
    theta={j:sum((w*t[0] for w,t in zip(W,trig[j])),ZERO) for j in nu}
    def pss(i,j):
        u,v=nu[i],nu[j]
        if i==j:return -sum((w*((L-t)/a*ct+st/(a*u)) for w,t,(st,ct) in zip(W,TS,trig[i])),ZERO)
        return -2/a*(v*theta[i]-u*theta[j])/(v.sq()-u.sq())
    result={}
    for oldbank in frozen['blocks']:
        k=oldbank['parity'];order,H,Mass,B,Hbb,c,Hb=g['build'](k)
        low=[order.index(j) for j in oldbank['ids']]
        for mat,key in [(H,'H'),(Mass,'M')]:
            for i,u in enumerate(low):
                for j,v in enumerate(low):
                    assert [max(mat[u][v].l,mat[v][u].l),min(mat[u][v].h,mat[v][u].h)]==oldbank[key][i][j]
        high=[j for j in order if j>85];anchor=high[0]
        # Coordinate tuples (finite sine coefficients, complete pole coefficient).
        lowv=[({j:I(1),k:-c[j]/c[k]},ZERO) if j!=k else ({k:I(1)},-c[k]/B) for j in oldbank['ids']]
        highv=[({j:I(1),anchor:-c[j]/c[anchor]},ZERO) for j in high[1:]]
        hp=({},1/B)
        def cross(v,w):
            d,b=v;e,f=w
            return (sum((x*y*g['Hss'](i,j) for i,x in d.items() for j,y in e.items()),ZERO)
                 +b*sum((y*Hb[j] for j,y in e.items()),ZERO)
                 +f*sum((x*Hb[i] for i,x in d.items()),ZERO)+b*f*Hbb)
        vectors=lowv+highv+[hp];n=len(vectors);nl=len(lowv);nh=len(highv)
        base=[[ZERO]*n for _ in range(n)];C=[[ZERO]*nh for _ in range(nh)]
        for i in range(n):
            for j in range(i,n):base[i][j]=base[j][i]=cross(vectors[i],vectors[j])
        base[-1][-1]+=2 if k==1 else -2
        for i in range(nl):
            for j in range(nl):
                v=base[i][j];l,h=oldbank['H'][i][j];assert max(v.l,l)<=min(v.h,h)
        for i in range(nh):
            for j in range(i,nh):
                v=sum((x*y*(pss(p,q)+(I(diag(p)) if p==q else ZERO)) for p,x in highv[i][0].items() for q,y in highv[j][0].items()),ZERO)
                C[i][j]=C[j][i]=v
        rows=[]
        # Interpolate the exact high-high forms, on this finite subspace only.
        # At theta=0 this is the capped comparator; at theta=1 the actual form.
        for t in [F(0),F(1,2),F(9,10),F(99,100),F(1)]:
            A=[r[:] for r in base]
            for i in range(nh):
                for j in range(nh):
                    A[nl+i][nl+j]=(1-t)*C[i][j]+t*base[nl+i][nl+j]
            try:
                lo=schur(endpoints(A,-1,S));hi=schur(endpoints(A,1,S));assert lo<=hi
                status='NEGATIVE_COMPARATOR_WITNESS' if hi<0 and t<1 else 'NEGATIVE_ORIGINAL_WITNESS' if hi<0 else 'POSITIVE_FINITE_ONLY' if lo>0 else 'UNDECIDED'
                row={'theta':str(t),'schur_interval':[enc(lo)[0],enc(hi)[1]],'classification':status}
            except AssertionError as exc:
                # Do not reinterpret an unproved leading positive block as a
                # negative pole Schur scalar.
                row={'theta':str(t),'classification':'LEADING_BLOCK_NOT_CERTIFIED','detail':str(exc)}
            rows.append(row);print('parity',k,'theta',t,row['classification'],row.get('schur_interval'),flush=True)
        result[str(k)]={'low_dimension':nl,'high_pole_null_dimension':nh,'normalized_pole_dimension':1,
          'high_ids':high,'high_anchor':anchor,'rows':rows}
    out.write_text(json.dumps({'status':'CONDITIONAL_JOINT_COMPARATOR_FEASIBILITY',
      'sectors':result,'scope':'Finite restriction of lower comparator, not original negativity unless theta=1 negative. Positive finite results are not infinite positivity. Original analytic/interval premises inherited.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
