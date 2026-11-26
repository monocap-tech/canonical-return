#!/usr/bin/env python3
"""CR-ALP-080 candidate: selected coordinate second-moment inverse receiver.

Reuse the original CR-ALP-059 physical evaluator and verified shifted mesh.
Require its geometry and runtime assertions for each actual source input.
This program only transports interval endpoints and verifies a reciprocal
majorant; it does not independently validate floating arithmetic or continuum.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json,sys
ELL=F(51,200); U=F(1101,200); GAMMA=F(72,25); PRIME=F(21,8)
SHIFT=F(1,10**13); SCALE=10**30
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sqrtup(x):
    assert x>=0
    y=x*SCALE*SCALE
    return F(isqrt(y.numerator//y.denominator)+1,SCALE)
def read(p):
    rows={}
    for line in p.read_text().splitlines():
        z=line.split()
        if len(z)!=5 or z[-2:]!=['/','1e15']:continue
        n,l,h=z[:3];assert n not in rows
        rows[n]=(F(int(l),10**15),F(int(h),10**15))
        assert rows[n][0]<=rows[n][1]
    return rows
def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def upper_term(c,x):return c*(x[1] if c>=0 else x[0])
def receive(first,second,eps):
    f=read(first);s=read(second)
    assert set(f)=={'G','M','penalty',*(f'corr{i}' for i in range(7))}
    assert set(s)=={'proxy_G','V2','sine2','bridge2','N_proxy'}
    assert overlap(f['G'],s['proxy_G'])
    assert 0<s['proxy_G'][0]<=s['proxy_G'][1]<F(4,25)
    assert 0<=s['N_proxy'][0]<=s['N_proxy'][1]<2
    v,a,b,n=(s[x] for x in ('V2','sine2','bridge2','N_proxy'))
    assert overlap((max(F(0),v[0]-a[1]-b[1]),v[1]-a[0]-b[0]),n)
    g=(max(f['G'][0],s['proxy_G'][0]),min(f['G'][1],s['proxy_G'][1]))
    # Transfer ||r-rho||<=eps. The shift allowance must be the documented
    # per-source bound for the rational translations, not assumed from G.
    ge=eps*(2*sqrtup(g[1])+eps)
    G=(max(F(0),g[0]-ge),g[1]+ge)
    me=PRIME*ge+SHIFT
    M=(f['M'][0]-me,f['M'][1]+me)
    delta=PRIME*eps+SHIFT
    ne=delta*(2*sqrtup(n[1])+delta)
    N=(max(F(0),n[0]-ne),n[1]+ne)
    affine=(GAMMA*G[1]-M[0])/(ELL*U)
    scalar=G[1]/ELL
    # p_a(t) = [t^2-(2a+ell)t+a^2+2a ell]/(ell a^2).
    # p_a(t)-1/t = (t-ell)(t-a)^2/(ell a^2 t) >= 0.
    # Expand around gamma, then take sign-aware endpoint upper bounds.
    def q(a):
        c=GAMMA**2-(2*a+ELL)*GAMMA+a*a+2*a*ELL
        d=2*GAMMA-2*a-ELL
        return (N[1]+upper_term(d,M)+upper_term(c,G))/(ELL*a*a)
    trials={ELL,U,*(F(i,100) for i in range(26,551,5))}
    gm,mm,nm=((x[0]+x[1])/2 for x in (G,M,N))
    m1=GAMMA*gm+mm;m2=GAMMA**2*gm+2*GAMMA*mm+nm
    if m1>ELL*gm and m2>ELL*m1:
        a=(m2-ELL*m1)/(m1-ELL*gm)
        if ELL<=a<=U:trials.add(a)
    z,a=min((q(a),a) for a in trials)
    assert min(z,affine,scalar)>=0
    return {'true_G':list(map(str,G)),'true_M':list(map(str,M)),
        'true_N':list(map(str,N)),'transport_N':str(ne),
        'affine_upper':str(affine),'scalar_upper':str(scalar),
        'quadratic_upper':str(z),'interpolation_node':str(a),
        'inverse_upper':str(min(z,affine,scalar)),
        'first_log_sha256':sha(first),'second_log_sha256':sha(second)}
def selftest():
    for a in (ELL,F(1),F(3),U):
        for t in (ELL,F(1),F(2),F(4),U):
            p=(t*t-(2*a+ELL)*t+a*a+2*a*ELL)/(ELL*a*a)
            assert p-1/t==(t-ELL)*(t-a)**2/(ELL*a*a*t)>=0
    print('PASS exact quadratic identity')
def main(bank,logs,out):
    entries=json.loads((bank/'source_bank.json').read_text())['entries']
    E={(e['parity'],e['coordinate_zero_based']):e for e in entries}
    assert len(E)==len(entries)==85
    rows={}
    for second in sorted(logs.glob('second_*.log')):
        _,k,j=second.stem.split('_');j=int(j);e=E[k,j]
        assert F(e['norm2'])==1 and sha(bank/e['file'])==e['sha256']
        eps=F(5,10**9)+F(e['omitted_group_midpoint_norm_upper'])+SHIFT
        r=receive(logs/f'mode_{k}_{j}.log',second,eps)
        r['input_sha256']=e['sha256'];rows[f'{k}_{j}']=r
    assert rows,'No complete second-moment logs'
    out.write_text(json.dumps({'status':'CONDITIONAL_SELECTED_SECOND_MOMENTS',
        'columns':rows,'premises':'Source, geometry, shift and numerical enclosures inherited from original 057/058/059; no positivity or index promotion alone.'},indent=2)+'\n')
    print(out.read_text())
if __name__=='__main__':
    if sys.argv[1:]==['--selftest']:selftest()
    else:main(*map(Path,sys.argv[1:]))
