#!/usr/bin/env python3
"""Exact sufficient 071 two-column receiver; no cross integration required.

Conditional on 070 first receiver, common comparator, original correction
construction, and second-source physical proxy norm/error envelope.
Usage: python3 cr_alp_071_coarse_receiver.py CR-ALP-071-common-trial-bank.zip
"""
import hashlib,json,sys,zipfile
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

def up(x):
    assert x>=0
    s=10**30; y=x*s*s
    return F(isqrt(y.numerator//y.denominator)+1,s)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def lo(x):
    return str(F(x*10**18//1,10**18))
def hi(x):
    return str(-F((-x)*10**18//1,10**18))
def run(path):
    assert sha(Path(path).read_bytes())=='ebbc51e56e67e7438483928b3570a320b77df06b92a574ad2a0e998fdb4b0ead'
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
        manifest={name:h for h,name in (l.split() for l in z.read('SHA256SUMS').decode().splitlines())}
        assert len(manifest)==11 and set(z.namelist())==set(manifest)|{'SHA256SUMS'}
        for name,h in manifest.items():assert sha(z.read(name))==h
        assert manifest['bank.json']=='b2d5f187293de20f05a9e510c1f54aeb1828fc2a20f55a0e4127546a7d99f71e'
        assert manifest['second_coeff.log']=='1e507cfb7ba8f9726f16eb4471359e628660d1fada62163073b16a7a1da35f50'
        b=json.loads(z.read('bank.json'))
        assert b['dimension']==2 and b['first_input_byte_match'] and F(b['euclidean_gram_determinant'])>0
        assert b['second_coordinate_zero_based']==37 and F(b['banks']['second']['norm2'])==1
        rows=[l.split() for l in z.read('second_coeff.log').decode().splitlines() if l.startswith('G ')]
        assert len(rows)==1 and rows[0][3:]==['/','1e15']
        Glo,Ghi=[F(int(x),10**15) for x in rows[0][1:3]]
        assert 0<Glo<=Ghi
    dmin=F(26535863,100000000);Acap=F(555992307,125000000)-F(72,25)
    assert dmin>0 and Acap>0
    # U_i is the projection of its exact finite sine precursor; Q85 kills
    # the pole and low sines, and projection contractivity handles the full tail.
    aa={}
    for tag in ('first','second'):
        coeff=b['banks'][tag]['correction_exact_intervals']
        assert set(map(int,coeff))==set(range(2,851,2))
        n2=F(0)
        for j,lh in coeff.items():
            l,h=map(F,lh);assert l<=h
            if int(j)>85:n2+=max(abs(l),abs(h))**2
        aa[tag]=Acap*n2
    # Deliberately loose total L2 proxy error premise for second direction.
    tail=F(b['banks']['second']['omitted_group_midpoint_norm_upper'])
    eps=F(1,1000)
    assert F(5,10**9)+tail<eps # arithmetic only; does not prove the envelope
    b2=(up(Ghi)+eps)**2/dmin
    B2=(up(b2)+up(aa['second']))**2
    B1=F(16059528681514303,62500000000000000)
    m1=F(61214171350249519,10**18)
    m2=F(b['K_L_22'][0])-B2
    c=max(abs(F(x)) for x in b['K_L_12'])+up(aa['first']*aa['second'])+up(B1*B2)
    det=m1*m2-c*c
    return {'status':'PASS_CONDITIONAL_PAIR' if m1>0 and m2>0 and det>0 else 'UNDECIDED',
            'premise':'second proxy squared norm <= Ghi and total L2 proxy error <= 1/1000; inherited 070 premises and original correction attachment',
            'authenticated_files':len(manifest),'Ghi':str(Ghi),'omitted_group_norm_upper':hi(tail),'proxy_error_assumed':str(eps),
            'a1_upper':hi(aa['first']),'a2_upper':hi(aa['second']),
            'second_old_inverse_upper':hi(b2),'second_corrected_inverse_upper':hi(B2),
            'first_pivot_lower':lo(m1),'second_diagonal_lower':lo(m2),'cross_abs_upper':hi(c),'determinant_lower':lo(det)}
if __name__=='__main__':print(json.dumps(run(sys.argv[1]),indent=2))
