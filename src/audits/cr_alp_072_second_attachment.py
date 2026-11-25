#!/usr/bin/env python3
"""CR-ALP-072: attach second projected residual via inherited 057/058/067 envelope.
Exact arithmetic audit, not a new proof of the continuum or floating-point premises.
Usage: python cr_alp_072_second_attachment.py ARTIFACT_DIRECTORY
"""
import hashlib,json,sys,zipfile
from pathlib import Path
from fractions import Fraction as F
from math import isqrt

def up(x):
    assert x>=0
    S=10**30;t=x*S*S
    return F(isqrt(t.numerator//t.denominator)+1,S)
def open_zip(p,h):
    assert hashlib.sha256(p.read_bytes()).hexdigest()==h
    z=zipfile.ZipFile(p);assert z.testzip() is None
    return z
def read_log(raw):
    rows={}
    for l in raw.decode().splitlines():
        n,a,b,slash,scale=l.split()
        assert slash=='/' and scale=='1e15' and n not in rows
        rows[n]=F(int(a),10**15),F(int(b),10**15)
        assert rows[n][0]<=rows[n][1]
    return rows
def run(p):
    z=open_zip(p/'CR-ALP-071-common-trial-bank.zip','ebbc51e56e67e7438483928b3570a320b77df06b92a574ad2a0e998fdb4b0ead')
    manifest={n:h for h,n in (l.split() for l in z.read('SHA256SUMS').decode().splitlines())}
    assert len(manifest)==11 and set(z.namelist())==set(manifest)|{'SHA256SUMS'}
    assert all(hashlib.sha256(z.read(n)).hexdigest()==h for n,h in manifest.items())
    raw=z.read('bank.json')
    assert hashlib.sha256(raw).hexdigest()=='b2d5f187293de20f05a9e510c1f54aeb1828fc2a20f55a0e4127546a7d99f71e'
    bank=json.loads(raw)
    assert bank['second_coordinate_zero_based']==37 and F(bank['banks']['second']['norm2'])==1
    assert F(bank['euclidean_gram_determinant'])>0
    raw=z.read('input_2_second.txt')
    assert hashlib.sha256(raw).hexdigest()=='185d9a51218712afa45169b6069856dc360b000f4b66bfd67da08b79270bf3a0'
    t=iter(raw.decode().split());assert next(t)=='2'
    def rd():
        l,h=F(next(t)),F(next(t));assert l<=h;return max(abs(l),abs(h))
    f,db=rd(),rd();assert f<100 and db<100
    assert int(next(t))==425;s=x=F(0)
    for i in range(425):
        assert int(next(t))==2*(i+1);s+=rd();x+=rd()
    assert s+x<10**17
    assert int(next(t))==128;low=F(0)
    for _ in range(128):
        assert rd()<256;low+=rd()
    assert low<10**14
    assert int(next(t))==25;group=F(0)
    for _ in range(25):
        assert rd()<=3*2**31-F(1,2)
        assert int(next(t))==40
        group+=sum(rd() for _ in range(40))
    assert group<1000 and next(t,None) is None
    # Original validate_geometry.py checks exactly these coefficient budgets.
    # 067 response.cpp reuses the original midpoint evaluator, projects on
    # mesh[0], then integrates on the exact 323-panel response mesh.
    a=open_zip(p/'CR-ALP-071-physical-second.zip','c4d06b7b19e26d708fb739746ded6627231c2ca75b2ff9213a16ac1c2b0b9eac')
    raw=a.read('response_second.log')
    assert hashlib.sha256(raw).hexdigest()=='ac4d653c09fe4363c54e54b4783c5001e4bc955297e8c57d25cad0e195bd60f7'
    r=read_log(raw);assert len(r)==2054
    direct=read_log(z.read('second_coeff.log'))
    assert r['G']==direct['G'] and 0<r['G'][0]<=r['G'][1]<1
    tail=F(bank['banks']['second']['omitted_group_midpoint_norm_upper'])
    assert tail<F(131,10**15)
    # Inherited uniform CR-ALP-057 operator proxy error, unit source,
    # omitted grouped terms, and conservative 058 rational shift allowance.
    eps=F(5,10**9)+tail+F(1,10**13)
    rho=up(r['G'][1])+eps
    assert rho<F(6,1000)
    # Independently sufficient pair test using only the older coarse correction
    # bounds, avoiding dependence on the new capped correction Gram checker.
    Acap=F(555992307,125000000)-F(72,25);dmin=F(26535863,100000000)
    aa={}
    for tag in ('first','second'):
        coeff=bank['banks'][tag]['correction_exact_intervals']
        assert set(map(int,coeff))==set(range(2,851,2))
        n2=F(0)
        for j,(l,h) in coeff.items():
            l,h=F(l),F(h);assert l<=h
            if int(j)>85:n2+=max(abs(l),abs(h))**2
        aa[tag]=Acap*n2
    B2=(rho*up(1/dmin)+up(aa['second']))**2
    B1=F(16059528681514303,62500000000000000)
    m1=F(61214171350249519,10**18)
    m2=F(bank['K_L_22'][0])-B2
    cross=max(abs(F(v)) for v in bank['K_L_12'])+up(aa['first']*aa['second'])+up(B1*B2)
    det=m1*m2-cross*cross
    assert m2>0 and det>0
    return {'status':'PASS_UNDER_INHERITED_ANALYTICAL_PREMISES',
            'source_sine_sum':str(s),'trial_sine_sum':str(x),
            'low_sum':str(low),'group_sum':str(group),
            'proxy_G_upper':str(r['G'][1]),'total_proxy_error_upper':str(eps),
            'actual_residual_norm_upper':str(rho),
            'second_inverse_upper':str(B2),'second_diagonal_lower':str(m2),
            'cross_abs_upper':str(cross),'determinant_lower':str(det)}
if __name__=='__main__':print(json.dumps(run(Path(sys.argv[1])),indent=2))
