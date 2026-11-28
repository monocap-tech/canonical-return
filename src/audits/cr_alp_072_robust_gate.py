#!/usr/bin/env python3
"""Robust conditional pair gate, without needing a sharp physical norm enclosure."""
from fractions import Fraction as F
from math import isqrt
import json, sys, zipfile
from cr_alp_071_coarse_receiver import run as coarse

def up(x):
    assert x>=0
    s=10**30; y=x*s*s
    return F(isqrt(y.numerator//y.denominator)+1,s)
def run(path):
    r=coarse(path)
    assert F(r['a1_upper'])<=F('3.159') and F(r['a2_upper'])<=F('0.00084')
    assert F(r['first_pivot_lower'])>=F('0.0612')
    with zipfile.ZipFile(path) as z:
        b=json.loads(z.read('bank.json'))
    assert F(b['K_L_22'][0])>=F('0.999')
    assert max(abs(F(x)) for x in b['K_L_12'])<=F('0.001561')
    assert F(16059528681514303,62500000000000000)<=F('0.257')
    assert F(26535863,100000000)>=F('0.265')
    # A single new premise: ||r_L,2|| <= 0.15, actual projected residual.
    B2=(F('0.15')*up(1/F('0.265'))+up(F('0.00084')))**2
    c=F('0.001561')+up(F('3.159')*F('0.00084'))+up(F('0.257')*B2)
    m2=F('0.999')-B2
    det=F('0.0612')*m2-c*c
    assert B2<=F(102635984,10**9) and c<=F(215485071,10**9)
    assert m2>=F(896364016,10**9) and det>=F(8423662,10**9)
    return {'status':'PASS_CONDITIONAL_ROBUST_GATE',
            'actual_second_residual_norm_required':'3/20',
            'second_inverse_upper':'0.102635984','second_diagonal_lower':'0.896364016',
            'cross_abs_upper':'0.215485071','determinant_lower':'0.008423662'}
if __name__=='__main__':print(json.dumps(run(sys.argv[1]),indent=2))
