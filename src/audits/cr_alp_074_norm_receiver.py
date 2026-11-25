#!/usr/bin/env python3
"""Exact norm-sum endpoint for the 074 balanced rank-20 source bank.

Requires an authenticated, independently regenerated selection and source bank,
plus native logs from the norm-only evaluator with verified geometry and runtime
checks. The result remains conditional on the inherited analytical/IEEE/MPFR
premises. No global RH or acceptance-ledger claim.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib,json,sys
S=10**30
ELL=F(51,200)
def up(x):
    assert x>=0;t=x*S*S
    return F(isqrt(t.numerator//t.denominator)+1,S)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(p):
    rows=p.read_text().splitlines()
    assert len(rows)==1
    n,l,h,s,d=rows[0].split()
    assert n=='G' and s=='/' and d=='1e15'
    l,h=F(int(l),10**15),F(int(h),10**15)
    assert 0<=l<=h<100
    return l,h
def main(bank,logs,out):
    sel=json.loads((bank/'selection.json').read_text())
    source=json.loads((bank/'source_bank.json').read_text())
    assert sel['rank']==20 and sel['conditional_gate_usable']
    assert F(sel['coercivity_lower'])==ELL
    assert len(source['entries'])==20
    manifest={n:h for h,n in (line.split() for line in (bank/'SHA256SUMS').read_text().splitlines())}
    assert set(manifest)=={'selection.json','source_bank.json',*(e['file'] for e in source['entries'])}
    assert all(sha(bank/n)==h for n,h in manifest.items())
    report={}
    for k in ('1','2'):
        s=sel['sectors'][k];J=s['indices_zero_based'];assert len(set(J))==len(J)==10
        E=[e for e in source['entries'] if e['parity']==k]
        assert [e['coordinate_zero_based'] for e in E]==J
        sum2=F(0);rows=[]
        for e in E:
            j=e['coordinate_zero_based']
            assert e['file']==f'input_{k}_{j}.txt'
            assert e['sha256']==manifest[e['file']]
            assert F(e['norm2'])==1
            G=parse(logs/f'norm_{k}_{j}.log')
            tail=F(e['omitted_group_midpoint_norm_upper']);assert tail>=0
            eps=F(5,10**9)+tail+F(1,10**13)
            q=up(G[1])+eps
            sum2+=q*q
            rows.append({'coordinate_zero_based':j,'proxy_G_upper':str(G[1]),
                         'proxy_error_upper':str(eps),'actual_residual_norm_upper':str(q),
                         'log_sha256':sha(logs/f'norm_{k}_{j}.log')})
        lam=F(s['K_lower']);assert lam>0
        margin=lam-sum2/ELL
        report[k]={'rank':10,'K_lower':str(lam),'residual_norm_squared_sum_upper':str(sum2),
                   'inverse_penalty_upper':str(sum2/ELL),'Schur_lower':str(margin),
                   'status':'PASS_CONDITIONAL' if margin>0 else 'UNDECIDED','columns':rows}
    ok=all(r['status']=='PASS_CONDITIONAL' for r in report.values())
    result={'status':'PASS_CONDITIONAL_RANK20' if ok else 'UNDECIDED',
            'sectors':report,'conditional_index_cap_if_received':65 if ok else 66,
            'premises':'Authenticated source selection, inherited 057 uniform operator proxy error, original 058 numerical and classical continuum premises, L coercivity, parity orthogonality. No native Lean or ledger promotion.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main(*(Path(p) for p in sys.argv[1:]))
