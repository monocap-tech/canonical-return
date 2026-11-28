#!/usr/bin/env python3
"""CR-ALP-075 nested prefix receiver, from authenticated norm logs.
Require all source inputs to have passed original geometry and runtime checks.
Missing columns stop a prefix; they are never assigned zero residual.
Usage: script.py BANK LOGS OUTPUT
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys,hashlib
from cr_alp_074_norm_receiver import parse,up
from cr_alp_074_balanced_selection import ELL
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(bank,logs,out):
    sel=json.loads((bank/'selection.json').read_text())
    entries=json.loads((bank/'source_bank.json').read_text())['entries']
    assert sel['rank']==70 and len(entries)==70 and F(sel['coercivity_lower'])==ELL
    manifest={n:h for h,n in (l.split() for l in (bank/'SHA256SUMS').read_text().splitlines())}
    assert set(manifest)=={'selection.json','source_bank.json',*(e['file'] for e in entries)}
    assert all(sha(bank/n)==h for n,h in manifest.items())
    # Caller must verify origin and numerical runtime checks of each log.
    report={};total=0
    for k in ('1','2'):
        s=sel['sectors'][k];J=s['indices_zero_based'];assert len(J)==len(set(J))==35
        E=[e for e in entries if e['parity']==k];assert [e['coordinate_zero_based'] for e in E]==J
        rows=[];levels=[];sum2=F(0);best=0
        for j,e in zip(J,E):
            assert F(e['norm2'])==1 and e['file']==f'input_{k}_{j}.txt'
            assert e['sha256']==manifest[e['file']]
            log=logs/f'norm_{k}_{j}.log'
            if not log.is_file():break
            _,G=parse(log)
            tail=F(e['omitted_group_midpoint_norm_upper']);assert tail>=0
            eps=F(5,10**9)+tail+F(1,10**13)
            b=up(G)+eps;sum2+=b*b
            rows.append({'coordinate':j,'norm_upper':str(b),'log_sha256':sha(log),'input_sha256':e['sha256']})
            r=len(rows)
            if r<10:continue
            lam=F(s['prefix_levels'][str(r)]['K_lower'])
            gap=lam-sum2/ELL
            if gap>0:best=r
            levels.append({'rank':r,'energy_lower':str(lam),'residual_sum_upper':str(sum2),'margin':str(gap),'pass_conditional':gap>0})
        total+=best
        report[k]={'received_prefix_length':len(rows),'best_conditional_rank':best,'columns':rows,'levels':levels}
    result={'status':'CONDITIONAL_PREFIX_RECEIVER','rank':total,'sectors':report,
        'conditional_index_upper':min(66,85-total),
        'premises':'Original continuum identities, operator proxy bound, L coercivity, parity orthogonality, authenticated source and numerical log provenance. No native Lean or unrestricted RH claim.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run(*(Path(x) for x in sys.argv[1:]))
