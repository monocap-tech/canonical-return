#!/usr/bin/env python3
"""Generate exact coordinate-source proxies for the nested rank-60 candidate.
Source only: no physical integration, positivity or index decrement.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json,sys
from cr_alp_075_nested_selection import main as select
from cr_alp_074_balanced_selection import checked
D=10**80
Q=10**40
def fmt(n):
    return ('-' if n<0 else '')+str(abs(n)//D)+'.'+str(abs(n)%D).zfill(80)
def enc(v):
    x=F(v)*D;return fmt(x.numerator//x.denominator)+' '+fmt(-((-x.numerator)//x.denominator))
def mid(pair,S):
    l,h=map(int,pair);assert l<=h;return F((l+h)//2,S)
def run(p,out):
    out.mkdir(parents=True,exist_ok=True)
    select(p,out/'selection.json')
    sel=json.loads((out/'selection.json').read_text())
    assert sel['conditional_gate_usable'], 'No usable rank-60 threshold; do not prepare physical jobs'
    a=checked(p/'cr_alp_057/result.json','e0dbb8cd24322cd5f0a96a7d8c1b188c38b7f5d2f706d644b8268cbb6533a809')
    g=checked(p/'cr_alp_057/grouped_proxy.json','12c7fbc679f244f0016cd71acd7256cf00830bb40d14fb2fedf452cab7706c0c')
    t=checked(p/'cr_alp_057/trial.json','6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0')
    S=int(a['scale']);assert S==int(g['scale'])==int(t['scale'])==10**100
    report=[]
    for k in ('1','2'):
        aa,gg,tt=(x['sectors'][k] for x in (a,g,t))
        n=43 if k=='1' else 42
        assert tt['ids']==list(range(int(k),851,2))
        kept=[b for b in gg['groups'] if int(b['M'])<2**32]
        omitted=[b for b in gg['groups'] if int(b['M'])>=2**32]
        assert len(kept)==25 and len(gg['low_amplitudes'])==128
        for j in sel['sectors'][k]['indices_zero_based']:
            def v(row):
                assert len(row)==n;return mid(row[j],S)
            f,db=v(aa['pole_coefficient']),v(tt['Db'][0])
            low={i:v(row) for i,row in zip(aa['source_ids'],aa['source_sine_coefficients'])}
            trial={i:v(row) for i,row in zip(tt['ids'],tt['D'])}
            lines=[k,enc(f),enc(db),str(len(trial))]
            for i in sorted(trial):lines.append(str(i)+' '+enc(low.get(i,F(0)))+' '+enc(trial[i]))
            lines.append(str(len(gg['low_amplitudes'])))
            for m,row in enumerate(gg['low_amplitudes']):
                lines.append(enc(F(4*m+1,2))+' '+enc(v(row)))
            lines.append(str(len(kept)))
            for b in kept:
                assert len(b['amplitudes'])==40
                lines.append(enc(F(b['center']))+' 40')
                lines.extend(enc(v(row)) for row in b['amplitudes'])
            raw=('\n'.join(lines)+'\n').encode()
            name=f'input_{k}_{j}.txt';(out/name).write_bytes(raw)
            tail=0
            for b in omitted:
                c=F(b['center']);x=2*c*Q**2
                root=isqrt(x.numerator//x.denominator);assert root>0
                val=F(2*Q,root)*sum(abs(v(row)) for row in b['amplitudes'])
                q=val*Q;tail+=-((-q.numerator)//q.denominator)
            report.append({'parity':k,'coordinate_zero_based':j,'norm2':'1','file':name,
              'sha256':hashlib.sha256(raw).hexdigest(),'omitted_group_midpoint_norm_upper':str(F(tail,Q)),
              'residual_norm_ceiling':sel['sectors'][k]['common_actual_residual_norm_ceiling']})
    assert len(report)==60
    (out/'source_bank.json').write_text(json.dumps({'status':'SOURCE_ONLY','entries':report},indent=2)+'\n')
    print('PASS source construction, 60 exact coordinate midpoints; no physical integration')
if __name__=='__main__':run(Path(sys.argv[1]),Path(sys.argv[2]))
