#!/usr/bin/env python3
"""Exact quantitative full-collar receiver, conditional on common physical attachment."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
from cr_alp_075_shifted_energy import lower_matrix,S,positive
from cr_alp_074_balanced_selection import checked,HASHES,ELL

def sqfrob(rows,scale):
    z=0
    for row in rows:
        for l,h in row:
            l,h=int(l),int(h);assert l<=h
            z+=max(abs(l),abs(h))**2
    return F(z,scale**2)

def positive_shift(N,shift):
    assert (shift*S).denominator==1
    return positive([[x-(int(shift*S) if i==j else 0) for j,x in enumerate(row)] for i,row in enumerate(N)])

def main(root,custody,frozen,out):
    a=checked(root/'cr_alp_057/result.json','e0dbb8cd24322cd5f0a96a7d8c1b188c38b7f5d2f706d644b8268cbb6533a809')
    t=checked(root/'cr_alp_057/trial.json','6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0')
    prior=json.loads((custody/'replayed-082/result.json').read_text())
    g=json.loads((custody/'replayed-091/result.json').read_text())
    f=json.loads(frozen.read_text())
    assert f['status']=='CONDITIONAL_FROZEN_TERMINAL_AUDIT_PASS' and f['rank']==85 and f['conditional_index_upper']==0
    assert int(a['scale'])==int(t['scale'])
    D=int(a['scale']);result={}
    for k,n,delta in [('1',43,F(1,10**5)),('2',42,F(1,10**4))]:
        J=prior['sectors'][k]['coordinate_order'];assert sorted(J)==list(range(n))
        w=checked(root/'cr_alp_056'/f'result_extended_{k}.json',HASHES[k])
        H=lower_matrix(w,J)
        gs=g['sectors'][k];assert int(gs['scale'])==S
        GL=[[int(x) for x in row] for row in gs['correction_gram_lower_grid']]
        K=[[H[i][j]+GL[i][j] for j in range(n)] for i in range(n)]
        assert K==list(map(list,zip(*K)))
        lam=next((F(1,10**p) for p in range(0,13) if positive_shift(K,F(1,10**p))),None)
        assert lam is not None
        kap=max(F(sum(abs(x) for x in row),S) for row in K)
        assert kap>=lam>0
        aud=f['sectors'][k]['audit'] if k=='1' else f['sectors'][k]['audits']['second_batch_combined']
        assert aud['dimension']==n and aud['full_pass']
        N=[[int(x) for x in row] for row in aud['integer_matrix']]
        assert N==list(map(list,zip(*N)))
        assert hashlib.sha256(json.dumps(N,separators=(',',':')).encode()).hexdigest()==aud['matrix_sha256']
        assert positive_shift(N,delta)
        # Actual finite coordinate p=d sine+f pole and correction X=D sine+Db pole.
        # Sines orthonormal, pole norm <2; elementary (u+v)^2<=2u^2+2v^2.
        aa=a['sectors'][k];tt=t['sectors'][k]
        assert aa['dimension']==n and len(aa['pole_coefficient'])==n
        assert all(len(row)==n for row in aa['source_sine_coefficients']+tt['D']+tt['Db'])
        assert tt['ids']==list(range(int(k),851,2))
        v2=2*sqfrob(aa['source_sine_coefficients'],D)+8*sqfrob([aa['pole_coefficient']],D)
        x2=2*sqfrob(tt['D'],D)+8*sqfrob(tt['Db'],D)
        assert v2>0
        # B<=K implies ||C^-1 r_C||^2<=kap/ELL. No new response integral.
        m=delta*lam/(kap+delta)
        M2=2*x2+2*kap/ELL
        cstar=min(m/(v2+2*M2),ELL/2)
        # Emit conservative decimal rational floor, not floating acceptance.
        e=next(p for p in range(1,40) if F(1,10**p)<=cstar)
        result[k]={'dimension':n,'lambda_K_lower':str(lam),'kappa_upper':str(kap),
            'terminal_shift':str(delta),'physical_coordinate_norm2_upper':str(v2),
            'correction_norm2_upper':str(x2),'square_shift_norm2_upper':str(M2),
            'finite_gap_lower':str(m),'physical_coercivity_power10':f'1e-{e}',
            'terminal_matrix_sha256':aud['matrix_sha256']}
        print(k,'lambda',lam,'kappa',float(kap),'V2',float(v2),'X2',float(x2),'coercivity',f'1e-{e}',flush=True)
    out.write_text(json.dumps({'status':'CONDITIONAL_QUANTITATIVE_COLLAR_RECEIVER','sectors':result,'premises':'Common original source and physical enclosures, simultaneous B_C compression majorants and K<=K_C, full tail comparator C>=ELL I; exact rational finite-to-Hilbert transport only. No native Lean or arbitrary-support claim.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
