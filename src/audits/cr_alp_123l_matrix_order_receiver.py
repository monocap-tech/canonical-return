#!/usr/bin/env python3
"""Exact post-receiver for 123I: 123K matrix-order comparison, no new physics.

Input 123I result must come from its successful frozen workflow, with authenticated
manifest. This script is a conditional receiver, not a validator of continuum premises.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json,hashlib,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
P='f5feb77462ff530ba85f3ee242e13b8e11b71144a65dfd3e16893188b05315de'
HASH={'1':'ef6a2ce2b37654f5785c37d4d793030e843eeb075a1d0380bf6242cc283aa402',
      '2':'42c495d86072e538250e00069268a726555ad3fb17c803d5db7651904da32c3d'}
DELTA={'1':F(1,10**5),'2':F(21,10**5)}
S=10**18
def sha(b):return hashlib.sha256(b).hexdigest()
def positive_shift(rows,delta):
    n=len(rows);assert (delta*S).denominator==1
    A=[[int(v)-(int(delta*S) if i==j else 0) for j,v in enumerate(row)] for i,row in enumerate(rows)]
    assert A==list(map(list,zip(*A)))
    prev=1
    for k in range(n):
        p=A[k][k]
        if p<=0:return False
        for i in range(k+1,n):
            for j in range(i,n):
                t,r=divmod(p*A[i][j]-A[i][k]*A[k][j],prev)
                assert r==0;A[i][j]=A[j][i]=t
        prev=p
    return True
def up_sqrt(v):
    assert v>=0
    D=10**150
    z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)
def receive(E,U,b,A,B):
    assert E>0 and U>=0 and b>=0 and A>0 and B>0
    # Exact strict comparison, including cross term, with no square-root rounding.
    margin=E-(1+B)*U-A*b*b
    cross2=4*A*B*b*b*U
    passed=margin>0 and margin*margin>cross2
    lo=margin-up_sqrt(cross2)
    assert not (lo>0 and not passed)
    return {'strict_sign_test_pass':passed,'margin_before_cross':str(margin),
            'cross_term_squared':str(cross2),'outward_schur_lower':str(lo),
            'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if passed else 'SUFFICIENT_BOUND_INCONCLUSIVE',
            'zero_finite_residual_U_budget':str(E/(1+B)),
            'zero_inverse_error_b2_budget':str(E/A),
            'inverse_error_budget_ratio':str(U*(1+B)/E),
            'finite_residual_budget_ratio':str(A*b*b/E)}
def main(residual,quant,terminal,out):
    raw=quant.read_bytes();assert sha(raw)==P
    q=json.loads(raw);assert q['status']=='CONDITIONAL_QUANTITATIVE_COLLAR_RECEIVER'
    f=json.loads(terminal.read_text());assert f['status']=='CONDITIONAL_FROZEN_TERMINAL_AUDIT_PASS' and f['rank']==85
    rraw=residual.read_bytes();r=json.loads(rraw);assert r['status']=='CONDITIONAL_FULL_AUXILIARY_RESIDUAL_ENVELOPE'
    assert set(r['sectors'])=={'1','2'}
    result={}
    for k in ('1','2'):
        a=f['sectors'][k]['audit'] if k=='1' else f['sectors'][k]['audits']['second_batch_combined']
        N=[[int(v) for v in row] for row in a['integer_matrix']]
        assert a['full_pass'] and a['dimension']==len(N)==q['sectors'][k]['dimension']
        assert sha(json.dumps(N,separators=(',',':')).encode())==a['matrix_sha256']==HASH[k]==q['sectors'][k]['terminal_matrix_sha256']
        assert positive_shift(N,DELTA[k])
        lam=F(q['sectors'][k]['lambda_K_lower']);kap=F(q['sectors'][k]['kappa_upper'])
        A=1/lam+1/DELTA[k];B=kap/DELTA[k]
        rows=r['sectors'][k];assert [x['cutoff'] for x in rows]==[125,165]
        rec=[]
        for row in rows:
            assert row['finite_prefix']==2048
            E=F(row['joint_energy_lower']);U=F(row['weighted_inverse_error_upper']);b=F(row['finite_corrected_residual_norm_upper'])
            x=receive(E,U,b,A,B);x['cutoff']=row['cutoff'];rec.append(x)
            print('parity',k,'cutoff',row['cutoff'],x['classification'],
                  'inverse budget ratio',float(F(x['inverse_error_budget_ratio'])),
                  'finite budget ratio',float(F(x['finite_residual_budget_ratio'])),flush=True)
        result[k]={'delta':str(DELTA[k]),'finite_factor':str(A),'inverse_factor':str(B),'rows':rec}
    out.write_text(json.dumps({'status':'CONDITIONAL_MATRIX_ORDER_POLE_RECEIVER',
        'residual_sha256':sha(rraw),'sectors':result,
        'scope':'Exact post-receiver of authenticated 123I enclosures under original common comparator and 121/122P physical premises. Strict positive comparison does not establish arbitrary-support positivity.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
