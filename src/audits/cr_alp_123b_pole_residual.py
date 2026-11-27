#!/usr/bin/env python3
"""Two-pole Galerkin trial: extend finite upper bounds and test residual obstruction.

Only original 046 analytical coefficient formulas are used. High-band
pole-null projection is a LOWER bound on the complete response norm.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
from cr_alp_123a_pole_schur import endpoints,schur,enc
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
PIN='554e7ab54da4545987d89aa660a4abe6c951cb63'
D=10**90

def solve(A,b):
    # Exact fraction-free solve with augmented integer column.
    n=len(A);B=[list(r)+[v] for r,v in zip(A,b)];p0=1
    for k in range(n-1):
        p=B[k][k];assert p>0
        for i in range(k+1,n):
            for j in range(k+1,n+1):
                v,r=divmod(p*B[i][j]-B[i][k]*B[k][j],p0);assert r==0;B[i][j]=v
            B[i][k]=0
        p0=p
    assert B[-1][-2]>0
    x=[F(0)]*n
    for i in reversed(range(n)):
        x[i]=(B[i][-1]-sum(B[i][j]*x[j] for j in range(i+1,n)))/F(B[i][i])
    assert all(sum(a*v for a,v in zip(r,x))==v for r,v in zip(A,b))
    return x

def main(root,out):
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==PIN
    bankraw=(root/'cr_alp_046/CR-ALP-046_full_E85_intervals.json').read_bytes()
    assert hashlib.sha256(bankraw).hexdigest()=='d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22'
    frozen=json.loads(bankraw)
    source=raw.decode()
    assert source.count('M=90;SHIFT=512;N=85')==1
    source=source.replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=125')
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    assert source.count(old)==1
    source=source.replace(old,old+', Bmass, Hbb, c, Hb')
    g={'__file__':str(SRC),'__name__':'cr_alp_123b_extended'}
    exec(compile(source,str(SRC),'exec'),g)
    S,I,ZERO=g['S'],g['I'],g['ZERO'];assert S==frozen['scale']
    result={}
    for oldbank in frozen['blocks']:
        k=oldbank['parity'];order,H,M,B,Hbb,c,Hb=g['build'](k)
        low=[order.index(j) for j in oldbank['ids']]
        high=[j for j in order if j>85]
        for mat,key in [(H,'H'),(M,'M')]:
            for i,u in enumerate(low):
                for j,v in enumerate(low):
                    assert [max(mat[u][v].l,mat[v][u].l),min(mat[u][v].h,mat[v][u].h)]==oldbank[key][i][j]
        s=1 if k==1 else -1;h0=Hbb/B.sq()+2*s
        u={j:Hb[j]/B-c[j]*Hbb/B.sq() for j in order}
        cross=[u[j] if j==k else u[j]-c[j]/c[k]*u[k] for j in order]
        complete=[row+[cross[i]] for i,row in enumerate(H)]+[cross+[h0]]
        lo=schur(endpoints(complete,-1,S));hi=schur(endpoints(complete,1,S));assert lo<=hi
        # Solve a symmetric rational midpoint system; the rounded trial is
        # an exact vector, not a certified exact Galerkin solution.
        A=[[(H[i][j].l+H[i][j].h+H[j][i].l+H[j][i].h)//4 for j in low] for i in low]
        assert A==list(map(list,zip(*A)))
        rhs=[-(cross[i].l+cross[i].h)//2 for i in low]
        x=[F(round(v*D),D) for v in solve(A,rhs)]
        # Actual normalized source f=h+sum x_j w_j, h=b/B.
        coef={j:I(v) for j,v in zip(oldbank['ids'],x)}
        d={j:coef[j] for j in oldbank['ids']}
        d[k]=d[k]-sum((coef[j]*c[j]/c[k] for j in d if j!=k),ZERO)
        f=(1-coef[k]*c[k])/B
        trial_energy=h0+2*sum((I(v)*cross[i] for i,v in zip(low,x)),ZERO)+sum((I(v)*H[i][j]*I(w) for i,v in zip(low,x) for j,w in zip(low,x)),ZERO)
        # Raw sine residual, excluding the pole term. The projection into
        # the finite high-sine subspace orthogonal to b kills that term.
        r=[sum((d[i]*g['Hss'](i,j) for i in d),ZERO)+f*Hb[j] for j in high]
        ch=[c[j] for j in high]
        polemass=sum((v.sq() for v in ch),ZERO);assert polemass.l>0
        rn=sum((v.sq() for v in r),ZERO)-sum((a*b for a,b in zip(ch,r)),ZERO).sq()/polemass
        # Tight residual lower bound is useful even when it cannot certify
        # the sign. The simple global epsilon norm test would require
        # ||projected residual||^2 < epsilon * W(f,f).
        epsilon=F(1,10**38) if k==1 else F(1,10**34)
        thresh=epsilon*F(trial_energy.h,S)
        decisive=F(rn.l,S)>max(F(0),thresh)
        result[str(k)]={'expanded_dimension':len(order),'original_matrix_match':True,
            'expanded_finite_schur_interval':[enc(lo)[0],enc(hi)[1]],
            'E85_trial_energy_interval':[str(F(trial_energy.l,S)),str(F(trial_energy.h,S))],
            'E85_trial_coefficients':[str(v) for v in x],'E85_trial_order':oldbank['ids'],
            'high_band_ids':high,'high_band_projected_residual_norm2_interval':[str(F(rn.l,S)),str(F(rn.h,S))],
            'global_coercivity_residual_threshold_upper':str(thresh),
            'simple_global_norm_test_insufficient':decisive,
            'classification':'NEGATIVE_FINITE_WITNESS' if hi<0 else 'POSITIVE_FINITE_ONLY' if lo>0 else 'UNDECIDED'}
        print('parity',k,'expanded dimension',len(order),'Schur',float(lo),float(hi),
              'E85 trial energy',float(F(trial_energy.h,S)),'high residual squared lower',float(F(rn.l,S)),
              'global threshold',float(thresh),'simple norm test insufficient',decisive,flush=True)
    out.write_text(json.dumps({'status':'CONDITIONAL_POLE_RESIDUAL_CENSUS','sectors':result,
      'scope':'Original continuum/interval premises. Finite Schur upper bounds and finite high-band projected residual LOWER bounds. A failed global norm sufficient test is not a negative witness. No infinite positivity.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
