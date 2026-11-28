#!/usr/bin/env python3
"""Generate exact Gram-residual witnesses for all six AQ2 matrix tests.
Floating Cholesky proposes factors only; every acceptance is exact integer.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
import numpy as np
sys.set_int_max_str_digits(0)
def read(path,sha):
    raw=Path(path).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==sha
    return json.loads(raw)['sectors']
def ints(rows):
    return np.array([[int(x) for x in row] for row in rows],dtype=object)
def strings(A):
    return [[str(int(x)) for x in row] for row in A]
def main(out):
    a=read('custody/cr-alp-123aq2-frozen-response/result.json','daafe6fb527d6e4cc0a9322bbc1ea8555b502192b826bacefda4a3839b918261')
    b=read('custody/cr-alp-123al-metric-first/result.json','b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e')
    S=10**12
    result={}
    for k,s in a.items():
        n=s['dimension'];T=int(s['matrix_grid']);H=ints(s['corrected_energy_lower']);P=ints(s['prefix_inverse_gram_upper'])
        Q=int(b[k]['whitening_grid']);B=ints(b[k]['whitening_matrix'])
        p=F(s['prefix_relative_upper']);assert p==F(b[k]['prefix_relative_upper'])
        I=np.eye(n,dtype=object)
        for name,M,D,delta in [
            ('H',H,T,F(1,int(k))),
            ('W',B@H@B.T-T*Q*Q*I,T*Q*Q,F(1,2**20)),
            ('P',p.numerator*H-p.denominator*P,T*p.denominator,F(1,2**(16-int(k))))]:
            assert M.shape==(n,n) and np.array_equal(M,M.T) and D>0
            A=np.array([[float(F(int(x),D)) for x in row] for row in M])
            L=np.linalg.cholesky(A-float(delta)*np.eye(n))
            N=np.array([[int(round(x*S)) for x in row] for row in L],dtype=object)
            R=delta.denominator*S*S*M-delta.denominator*D*(N@N.T)-delta.numerator*D*S*S*I
            maxr=max(abs(int(v)) for v in R.flat)
            assert 2*n*maxr<delta.numerator*D*S*S
            ratio=F(2*n*maxr,delta.numerator*D*S*S)
            result[k+name]={'dimension':n,'matrix_denominator':str(D),'matrix_numerators':strings(M),'factor_denominator':str(S),'factor_numerators':strings(N),'delta':str(delta),'epsilon':str(delta/(2*n)), 'residual_max_numerator':str(maxr),'residual_denominator':str(delta.denominator*D*S*S),'relative_residual_upper':str(ratio),'coercivity_lower':str(delta/2)}
            print(k,name,n,'exact reserve ratio',float(ratio),'PASS',flush=True)
    assert len(result)==6
    payload=json.dumps({'status':'EXACT_INTEGER_WITNESSES_PASS','matrices':result,'scope':'External integer verification only. Native entrywise checker and source interval attachment remain separate.'},indent=2)+'\n'
    Path(out).write_text(payload)
if __name__=='__main__':main(sys.argv[1])
