#!/usr/bin/env python3
"""123AO: remove historical inverse-metric premise from receiving chain.
Exact rational postprocessing; no new physical enclosures.
"""
import hashlib,json,sys,zipfile
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
sys.set_int_max_str_digits(0)
PINS=['a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473','a2c465dd692b0afbc0cd1a15b9f3599864085101bb8f52f7d225c1125adf453c','b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e','7c45c7f493062cc1236d3f8de6c6aec636dc474119d51676308849606a191897','cd7b591bc60a57bc41cc9657962815caa639b921cc8085d4e561704fb7b29d4e','2d6706539d46f1e5603c98845c780efb80d6062e902433b5a5bfd9b1df0b52ac']
def read(p,h,name):
    if zipfile.is_zipfile(p):
        with zipfile.ZipFile(p) as z:
            assert z.testzip() is None;b=z.read(name)
    else:b=Path(p).read_bytes()
    assert hashlib.sha256(b).hexdigest()==h,(p,hashlib.sha256(b).hexdigest())
    return json.loads(b)
def ceil(x):return -((-x).__floor__())
def up(x,g):return F(ceil(x*g),g)
def root(x):
    assert x>=0;g=10**110;v=x*g*g
    return F(isqrt(v.numerator//v.denominator)+1,g)
def main(ac,ak,al,z,am,an,out):
    A,K,L,Z,M,N=[read(p,h,'joint-repaired.json' if i==3 else 'result.json') for i,(p,h) in enumerate(zip((ac,ak,al,z,am,an),PINS))]
    assert all(set(d['sectors'])=={'1','2'} for d in (A,K,L,Z,M,N))
    ans={};g=10**100
    for k,a in A['sectors'].items():
        s=K['sectors'][k];l=L['sectors'][k];z=Z['sectors'][k][0];m=M['sectors'][k];n=N['sectors'][k];dim=44-int(k)
        E,U=map(F,(a['E'],a['U']));assert E==F(z['joint_energy_lower']) and U==F(z['weighted_inverse_error_upper']) and E>U>0
        assert n['classification']=='CONDITIONAL_RECEIVER_AUDIT_PASS'
        rho=F(n['independent_rho_upper']);assert rho==F(m['total_relative_upper']) and 0<rho<1
        perm=a['coordinate_order'];assert sorted(perm)==list(range(dim))
        L0=[[F(v) for v in row] for row in a['LDL_lower']];w=list(map(F,s['frame_roots']));assert all(v>0 for v in w)
        assert len(w)==len(L0)==dim and all(len(row)==dim and row[i]==1 and all(v==0 for v in row[i+1:]) for i,row in enumerate(L0))
        v=[tuple(map(F,z['conditioned_joint_residual_vector'][i])) for i in perm];assert all(x<=y for x,y in v)
        eb=F(z['bank_defect_norm2_upper']);assert eb>=0
        bank=a['bank_indices'];q=a['response_transformed_coordinates'];qg=int(a['response_coordinate_grid']);assert len(q)==len(bank) and all(len(row)==dim for row in q) and qg>0
        bounds=[]
        for j in range(dim):
            lo=hi=0
            for i,(x,y) in enumerate(v):
                u,t=sorted((x*L0[i][j],y*L0[i][j]));lo+=(u*g).__floor__();hi+=ceil(t*g)
            vabs=F(max(abs(lo),abs(hi)),g)
            norm2=F(sum(int(row[j])**2 for row in q),qg*qg)
            bounds.append(up(w[j]*(vabs+root(eb*norm2)),g))
        bframe=sum(x*x for x in bounds)
        Q=int(l['whitening_grid']);B=l['whitening_matrix'];assert len(B)==dim and all(len(row)==dim for row in B) and Q>0
        frob2=F(sum(int(x)**2 for row in B for x in row),Q*Q)
        bnew=up(frob2*bframe/(1-rho)/E,10**30)*E
        eta=rho/(1-rho);loss=up((U+(root(bnew)+root(eta*U))**2)/E,10**9);assert loss<1
        ans[k]={'dimension':dim,'rho_upper':str(rho),'whitening_frobenius2_upper':str(up(frob2,10**9)),'frame_charge_over_E_upper':str(up(bframe/E,10**30)),'new_finite_inverse_charge_over_E_upper':str(bnew/E),'loss_over_E_upper':str(loss),'positive_margin_over_E_lower':str(1-loss),'classification':'PASS_WITHOUT_HISTORICAL_INVERSE_METRIC_PREMISE'}
        print(k,ans[k],flush=True)
    Path(out).write_text(json.dumps({'status':'CONDITIONAL_DIRECT_SCHUR_ROUTE_PASS','sectors':ans,'scope':'AK/AL/AN coefficient, energy and Gram enclosures and original analytic comparator/domain premises inherited; historical K,N compression and Dbar inverse order not used.'},indent=2)+'\n')
if __name__=='__main__':main(*sys.argv[1:])
