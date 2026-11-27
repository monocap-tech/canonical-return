#!/usr/bin/env python3
"""Direction-specific exact matrix-order post-receiver for 123W signed conditioned joint correction.
No new special-function evaluation. Inherited continuum/comparator premises.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
from cr_alp_075_shifted_energy import lower_matrix,S,positive
from cr_alp_074_balanced_selection import checked,HASHES
from cr_alp_123j_common_denominator_solve import solve_common
from cr_alp_123l_matrix_order_receiver import positive_shift,HASH,DELTA
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def up_sqrt(v):
    assert v>=0
    D=10**150;z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)
def qform(A,v):
    # Exact b^T(A/S)^-1 b via common denominator solve.
    from math import lcm
    den=lcm(*(x.denominator for x in v))
    b=[int(x*den) for x in v]
    z,d=solve_common(A,b)
    q=F(S*sum(x*y for x,y in zip(b,z)),d*den*den)
    assert q>=0
    return q
def generalized(N,K,delta,kap):
    # Find certified t with N>=t K. No eigenvalue approximation.
    D=10**9
    lo=int(delta/kap*D);assert lo>=0
    def test(t):
        return positive([[D*x-t*y for x,y in zip(r,s)] for r,s in zip(N,K)])
    assert test(lo)
    hi=max(1,2*lo)
    while test(hi):
        lo=hi;hi*=2
        assert hi<=D*10**6
    while hi-lo>1:
        m=(lo+hi)//2
        if test(m):lo=m
        else:hi=m
    assert lo>0 and test(lo) and not test(hi)
    return F(lo,D)
def main(root,custody,terminal,quant,residual,out):
    raw=quant.read_bytes()
    assert sha(raw)=='f5feb77462ff530ba85f3ee242e13b8e11b71144a65dfd3e16893188b05315de'
    q=json.loads(raw);assert q['status']=='CONDITIONAL_QUANTITATIVE_COLLAR_RECEIVER'
    f=json.loads(terminal.read_text());assert f['status']=='CONDITIONAL_FROZEN_TERMINAL_AUDIT_PASS' and f['rank']==85
    rraw=residual.read_bytes();r=json.loads(rraw)
    assert r['status']=='CONDITIONAL_CONDITIONED_JOINT_CORRECTION'
    assert set(r['sectors'])=={'1','2'} and all(len(v)==1 and v[0]['cutoff']==850 and v[0]['source_cutoff']==365 for v in r['sectors'].values())
    p=json.loads((custody/'replayed-082/result.json').read_text())
    g=json.loads((custody/'replayed-091/result.json').read_text())
    result={}
    for k in ('1','2'):
        qq=q['sectors'][k];n=qq['dimension'];J=p['sectors'][k]['coordinate_order'];assert sorted(J)==list(range(n))
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',HASHES[k])
        H=lower_matrix(w,J);gs=g['sectors'][k];assert int(gs['scale'])==S
        GL=[[int(x) for x in row] for row in gs['correction_gram_lower_grid']]
        K=[[H[i][j]+GL[i][j] for j in range(n)] for i in range(n)]
        assert K==list(map(list,zip(*K)))
        lam,kap=F(qq['lambda_K_lower']),F(qq['kappa_upper'])
        assert positive_shift(K,lam)
        assert max(F(sum(abs(x) for x in row),S) for row in K)==kap
        a=f['sectors'][k]['audit'] if k=='1' else f['sectors'][k]['audits']['second_batch_combined']
        N=[[int(v) for v in row] for row in a['integer_matrix']]
        assert a['full_pass'] and a['dimension']==n
        assert sha(json.dumps(N,separators=(',',':')).encode())==a['matrix_sha256']==HASH[k]==qq['terminal_matrix_sha256']
        delta=DELTA[k];assert positive_shift(N,delta)
        t=generalized(N,K,delta,kap);B=1/t
        A=1/lam+1/delta
        rows=[]
        for row in r['sectors'][k]:
            X=[list(map(F,z)) for z in row['correction_response_vector']]
            assert len(X)==n and all(l<=h for l,h in X)
            # X now denotes V*R(w)+X*e', in the original V coordinates.
            # The old joint residual is the only omitted signed contribution.
            # K,N use 082 permutation J; apply same permutation to the covector.
            mid=[sum(X[i])/2 for i in J]
            rad=[(X[i][1]-X[i][0])/2 for i in J]
            err=up_sqrt(sum(v*v for v in rad))+F(row['conditioned_joint_residual_norm_upper'])
            qmid=qform(K,mid)+qform(N,mid)
            finite=up_sqrt(qmid)+up_sqrt(A)*err
            U=F(row['weighted_inverse_error_upper']);E=F(row['joint_energy_lower'])
            assert U>=0
            if E<=0:
                rows.append({'cutoff':row['cutoff'],'joint_energy_lower':str(E),
                    'classification':'NONPOSITIVE_TRIAL_ENERGY_LOWER_INCONCLUSIVE'})
                print(k,row['cutoff'],'NONPOSITIVE_TRIAL_ENERGY_LOWER_INCONCLUSIVE',flush=True)
                continue
            loss=U+(finite+up_sqrt(B*U))**2
            lower=E-loss
            rec={'cutoff':row['cutoff'],'midpoint_inverse_quadratic':str(qmid),
                 'finite_response_norm_upper':str(finite),'inverse_error_upper':str(U),
                 'generalized_inverse_factor':str(B),'loss_upper':str(loss),'conditional_pole_schur_lower':str(lower),
                 'finite_charge_ratio':str(finite*finite/E),'inverse_charge_ratio':str((1+B)*U/E),
                 'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if lower>0 else 'SUFFICIENT_BOUND_INCONCLUSIVE'}
            rows.append(rec)
            print(k,row['cutoff'],rec['classification'],'finite charge/E',float(finite*finite/E),'inverse charge/E',float((1+B)*U/E),'total loss/E',float(loss/E),flush=True)
        result[k]={'generalized_shift':str(t),'scalar_shift':str(delta),'coordinate_order':J,'rows':rows}
        print('generalized',k,float(t),'inverse factor',float(B),flush=True)
    out.write_text(json.dumps({'status':'CONDITIONAL_JOINT_CORRECTION_DIRECTION_SPECIFIC_RECEIVER','source_residual_sha256':sha(rraw),
        'sectors':result,'scope':'Exact direction-specific inverse of frozen rational lower matrices K,N. Inherits all original physical comparison premises. No arbitrary support conclusion.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
