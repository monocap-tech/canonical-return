#!/usr/bin/env python3
"""Enclose the finite pole Schur scalars using original 046 interval arithmetic.

The finite minima are upper bounds on the infinite-domain Schur scalars.
Positive finite minima do not establish positivity of the infinite problem.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,runpy,sys
SRC=Path(__file__).with_name('cr_alp_046_full_finite.py')
PIN='554e7ab54da4545987d89aa660a4abe6c951cb63'
GRID=10**110
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def endpoints(a,b,scale):
    n=len(a);C=[[0]*n for _ in range(n)];r=[0]*n
    for i in range(n):
        for j in range(n):
            l=max(a[i][j].l,a[j][i].l);h=min(a[i][j].h,a[j][i].h);assert l<=h
            lo=floor(F(l*GRID,scale));hi=ceil(F(h*GRID,scale))
            m=(lo+hi)//2;C[i][j]=m;r[i]+=max(m-lo,hi-m)
    assert C==list(map(list,zip(*C)))
    return [[C[i][j]+(b*r[i] if i==j else 0) for j in range(n)] for i in range(n)]
def schur(A):
    B=[row[:] for row in A];prev=1
    for k in range(len(B)-1):
        p=B[k][k];assert p>0,'Finite pole-null block not positive at working grid'
        for i in range(k+1,len(B)):
            for j in range(i,len(B)):
                v,r=divmod(p*B[i][j]-B[i][k]*B[k][j],prev)
                assert r==0;B[i][j]=B[j][i]=v
        prev=p
    return F(B[-1][-1],prev*GRID)
def enc(x):
    D=10**90
    return [str(F(floor(x*D),D)),str(F(ceil(x*D),D))]
def main(root,out):
    assert hashlib.sha1(b'blob '+str(len(SRC.read_bytes())).encode()+b'\0'+SRC.read_bytes()).hexdigest()==PIN
    raw=(root/'cr_alp_046/CR-ALP-046_full_E85_intervals.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='d877bf95dd203c9492174afcb48c91a90bb24d9e2f5bbb6ddce5138c83c52a22'
    frozen=json.loads(raw)
    # Execute the original elementary and coefficient generator unchanged, then
    # expose its existing local pole coefficients without changing their formulas.
    g=runpy.run_path(str(SRC),run_name='cr_alp_046_import')
    src=SRC.read_text().split('def build(k):',1)[1].split('def bareiss_prefix',1)[0]
    old=' return order,[[convert(H,i,j) for j in order] for i in order],[[convert(Mass,i,j) for j in order] for i in order]'
    assert src.count(old)==1
    new=old+', Bmass, Hbb, c, Hb'
    exec('def build(k):'+src.replace(old,new),g)
    S=g['S'];I=g['I'];assert S==frozen['scale']
    result={}
    for b in frozen['blocks']:
        k=b['parity'];order,H,M,B,Hbb,c,Hb=g['build'](k);assert order==b['ids']
        for mat,key in [(H,'H'),(M,'M')]:
            for i in range(len(mat)):
                for j in range(len(mat)):
                    assert [max(mat[i][j].l,mat[j][i].l),min(mat[i][j].h,mat[j][i].h)]==b[key][i][j]
        # h=b/B, original pole term equals +2 for cosh and -2 for sinh.
        # Projected sines have zero pole term, so cross entries use signed
        # Gamma+prime Hb,Hbb only, followed by the same 046 basis transform.
        s=1 if k==1 else -1
        h0=Hbb/B.sq()+2*s
        u={j:Hb[j]/B-c[j]*Hbb/B.sq() for j in order}
        cross=[u[j] if j==k else u[j]-c[j]/c[k]*u[k] for j in order]
        complete=[row+[cross[i]] for i,row in enumerate(H)]+[cross+[h0]]
        rows=[]
        for count in sorted({2,4,8,16,32,len(order)}):
            keep=[i for i,j in enumerate(order) if j<=k+2*(count-1)]
            assert len(keep)==count
            J=keep+[len(order)];A=[[complete[i][j] for j in J] for i in J]
            lo=schur(endpoints(A,-1,S));hi=schur(endpoints(A,1,S));assert lo<=hi
            rows.append({'dimension':count,'finite_schur_lower':enc(lo)[0],'finite_schur_upper':enc(hi)[1],
                         'classification':'POSITIVE_FINITE_ONLY' if lo>0 else 'NEGATIVE_FINITE_WITNESS' if hi<0 else 'UNDECIDED'})
            print('parity',k,'dimension',count,'finite Schur',float(lo),float(hi),rows[-1]['classification'],flush=True)
        result[str(k)]={'pole_energy_interval':[str(F(h0.l,S)),str(F(h0.h,S))],
            'original_projected_matrix_match':True,'rows':rows,
            'infinite_schur_upper':min((F(x['finite_schur_upper']) for x in rows))}
        result[str(k)]['infinite_schur_upper']=str(result[str(k)]['infinite_schur_upper'])
    out.write_text(json.dumps({'status':'CONDITIONAL_FINITE_POLE_SCHUR_CENSUS','sectors':result,
      'scope':'Original analytic coefficient and interval identities inherited. Finite Schur upper bounds transport to infinite variational scalar. No infinite lower bound or full positivity claim.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
