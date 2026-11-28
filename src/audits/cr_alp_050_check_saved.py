#!/usr/bin/env python3
"""Independent integer replay of saved CR-ALP-050 matrices; not an entry reconstruction."""
from fractions import Fraction as F
from pathlib import Path
import json, hashlib
ROOT=Path(__file__).resolve().parent

def positive(A,s,shift):
    d=s//10**65; C=[];rad=[]; sh=shift*10**65; assert sh.denominator==1
    for i,row in enumerate(A):
        R=[]; r=0
        for j,(lo,hi) in enumerate(row):
            l=int(lo)//d; h=-((-int(hi))//d)
            if i==j:l-=int(sh);h-=int(sh)
            m=(l+h)//2;R.append(m);r+=max(m-l,h-m)
        C.append(R);rad.append(r)
    n=len(C)
    for i in range(n):
        for j in range(i):assert C[i][j]==C[j][i]
        C[i][i]-=max(rad)
    prev=1
    for k in range(n):
        pivot=C[k][k];assert pivot>0
        for i in range(k+1,n):
            for j in range(i,n):
                q,r=divmod(C[i][j]*pivot-C[i][k]*C[k][j],prev);assert r==0;C[i][j]=C[j][i]=q
        prev=pivot
    return n

result={}
for parity in (1,2):
    p=ROOT/f'result_{parity}.json'; raw=p.read_bytes();assert hashlib.sha256(raw).hexdigest()=={1:'1a2d6451fcd0535e900fc2351a1d7e80ec392b87f49a7bc000661e9fdd310f00',2:'de44485825fa8902c42f35d60f07520137b812da73eafe13702f2cf31a7c2d5b'}[parity];d=json.loads(raw);s=int(d['scale']); n=44-parity
    assert d['trial_end']==410 and d['checked_end']==600
    assert len(d['Knew'])==n
    margin=F(1,4 if parity==1 else 8)
    assert positive(d['Knew'],s,margin)==n
    # independently compute a lower-test upper endpoint from Knew and Gfull,
    # not from the saved difference matrix.
    z=[]
    for lo,hi in d['witness']['z']:
        assert lo==hi;z.append(int(lo))
    num=0
    for i in range(n):
        for j in range(n):
            a=z[i]*z[j];kl,kh=map(int,d['Knew'][i][j]);gl,gh=map(int,d['Gfull'][i][j]);
            num+=a*(51*kh-200*gl if a>=0 else 51*kl-200*gh)
    val=F(num,51*s**3); assert val < (-F(2,5) if parity==1 else -1)
    upper=max(sum(max(abs(int(lo)),abs(int(hi))) for lo,hi in row) for row in d['Gtrial'])
    assert F(upper,s)<F(1,10**26)
    # Check the interval subtraction encloses the result of direct entry-wise arithmetic.
    for i in range(n):
        for j in range(n):
            kl,kh=map(int,d['Knew'][i][j]);gl,gh=map(int,d['Gfull'][i][j]);lo,hi=map(int,d['comparison'][i][j]);
            assert 51*lo<=51*kl-200*gh and 51*hi>=51*kh-200*gl
    print('parity',parity,'finite margin',margin,'pivots',n,'trial residual <1e-26; witness upper',float(val),'PASS')
    result[str(parity)]={'matrix_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'finite_margin':str(margin),'positive_pivots':n,'negative_witness_upper':str(val),'negative_witness_decimal':float(val)}
# Polynomial contraction calculation: t_m(x)=T_m(x), x=gamma/delta=192/175.
# Exact integer fractions bound the entire spectral interval, not a finite sample.
x=F(192,175); a=F(51,200); b=F(1101,200)
assert (b-a)/2==F(21,8) and (b+a)/2==F(72,25)
T0,T1=F(1),x
for m in range(1,17):
    if m==1:t=T1
    else:T0,T1=T1,2*x*T1-T0;t=T1
    if m in (4,8,12,16):
        factor=1/t**2
        print('degree',m,'squared residual contraction',float(factor),'correction penalty coefficient',float(factor/a))
        result[f'chebyshev_{m}']={'squared_contraction':str(factor),'scalar_penalty':str(factor/a)}
assert F(result['chebyshev_8']['squared_contraction'])<F(1,250)
assert F(result['chebyshev_12']['squared_contraction'])<F(1,8000)
(ROOT/'replay_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: saved matrix acceptance and exact Chebyshev constants; continuum entries not independently reconstructed')
