#!/usr/bin/env python3
"""Repair fixed-root-grid overestimation using authenticated saved moments.
No new special-function or physical evaluation. Preserve all original signed vectors.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def digest(b):return hashlib.sha256(b).hexdigest()
def up(v):return max(abs(F(x)) for x in v)
def root(v):
    assert v>=0
    D=10**110;z=v*D*D
    return F(isqrt(z.numerator//z.denominator)+1,D)
def main(second,joint,quant,out):
    raw=second.read_bytes();assert digest(raw)=='9152dcb759260d498c7496e579285d83e28172b6bcb3192c32e961886335253f'
    d=json.loads(raw);assert d['status']=='CONDITIONAL_SECOND_AUXILIARY_CORRECTION'
    jraw=joint.read_bytes();assert digest(jraw)=='09b34469b32c13ebf5f5b690d1eb9ad3a9fb5676f171a1c7f85325cedc199d7b'
    g=json.loads(jraw)
    qraw=quant.read_bytes();assert digest(qraw)=='f5feb77462ff530ba85f3ee242e13b8e11b71144a65dfd3e16893188b05315de'
    q=json.loads(qraw)
    J=2048;V=F(15*J,11);dc=F('1.822938456')
    for k,rows in d['sectors'].items():
        assert len(rows)==1
        r=rows[0];assert r['cutoff']==850 and r['source_cutoff']==365 and r['finite_prefix']==J
        w=next(x['joint_trial'] for x in g['sectors'][k]['rows'] if x['cutoff']==365)
        ds={int(i):up(v) for i,v in w['source_sines'].items()}
        fb=up(r['source_pole']);T=up(r['combined_pole'])
        # tau = fb - T, signed endpoints.
        tl,th=map(F,r['source_pole']);ul,uh=map(F,r['combined_pole'])
        tau=max(abs(tl-uh),abs(th-ul))
        ts={i:v+tau*F(33,15*i) for i,v in ds.items()}
        ts.update({int(i):abs(F(v)) for i,v in r['corrected_tail_coordinates'].items()})
        assert max(ts)<=850 and all(i>85 for i in map(int,r['corrected_tail_coordinates']))
        n=32;rho=(F(10*max(ts),7)/V)**2;assert rho<1
        coeff={}
        def charge(s,v):
            assert v>=0
            coeff[s]=coeff.get(s,F(0))+v
        def series(label,mult,shift,H):
            mom=r['signed_moments'][label];assert len(mom)==n
            for l,z in enumerate(mom):charge(2*l+2-shift,mult*up(z))
            charge(2*n+2-shift,mult*H/(1-rho))
        # |AG|,|AP|<4; c_j <= 33/(15j). Use independent rational
        # upper absolute moments, not the lost original remainder coefficients.
        H=sum(4*(ds.get(i,F(0))+t)*(F(10*i,7)**2+F(1,4))*F(10*i,7)**(2*n) for i,t in ts.items())
        series('centered_combined',F(20,11),-1,H)
        series('gamma_sine',F(80,11),0,sum(v*F(10*i,7)**(2*n+1) for i,v in ds.items()))
        series('prime_sine',F(80,11),0,sum(v*F(10*i,7)**(2*n+1) for i,v in ts.items()))
        charge(1,up(r['leading_combined_moment'])+80*fb+40*T)
        # sqrt(J^{1-2s}/(2s-1)) = sqrt(J/(2s-1))/J^s.
        # Root is now always O(1), independent of expansion order.
        tail=sum(z*root(F(J,2*s-1))/V**s for s,z in coeff.items())
        tail+=2*fb*F(11,15)*root(F(122,J))
        U=F(r['weighted_prefix_upper'])+tail*tail/dc # omit pole subtraction safely
        E=F(r['joint_energy_lower']);b=F(r['finite_corrected_residual_norm_upper'])
        pp=q['sectors'][k];delta=F(1,10**5) if k=='1' else F(21,10**5)
        A=1/F(pp['lambda_K_lower'])+1/delta;B=F(pp['kappa_upper'])/delta
        loss=U+(root(A)*b+root(B*U))**2
        old=F(r['infinite_tail_norm_upper'])
        r.update({'unrepaired_infinite_tail_norm_upper':str(old),'infinite_tail_norm_upper':str(tail),
                  'weighted_inverse_error_upper':str(U),'weighted_pole_subtraction_lower':'0',
                  'response_loss_upper':str(loss),'conditional_pole_schur_lower':str(E-loss),
                  'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE'})
        print(k,'repaired tail',float(tail),'prefix U',float(F(r['weighted_prefix_upper'])),
              'complete U',float(U),'U/E',float(U/E) if E>0 else None,
              'loss/E',float(loss/E) if E>0 else None,flush=True)
    d['repair']={'original_sha256':digest(raw),'method':'dimensionless integral root; independently bounded absolute remainder moments; omitted nonnegative weighted pole subtraction'}
    out.write_text(json.dumps(d,indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
