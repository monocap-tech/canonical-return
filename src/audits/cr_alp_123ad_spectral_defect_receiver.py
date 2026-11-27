#!/usr/bin/env python3
"""CR-ALP-123AD: exact spectral post-receiver for the 123AC defect bank.
No new response solves. Keep prefix signs and common tail oscillation geometry.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,sys
from cr_alp_075_shifted_energy import positive
from cr_alp_123z_signed_tail_readback import namespace
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
G=10**30
def ceil(x):return -((-x).__floor__())
def roots(v):
    assert v>=0
    z=v*G*G;t=isqrt(z.numerator//z.denominator)
    return F(t,G),F(t+1,G)
def up(v):return roots(v)[1]
def interval_product(a,b,c,d):
    z=[a*c,a*d,b*c,b*d];return min(z),max(z)
def box(a,b):
    assert a<=b
    m=round((a+b)*G/2)
    r=ceil(max(F(m,G)-a,b-F(m,G))*G)
    return m,r
def spectral(A):
    # Exact Loewner upper rounding then strict positive-shift bisection.
    n=len(A);assert A==list(map(list,zip(*A)))
    scale=10**14
    B=[[round(v*scale)+(n if i==j else 0) for j,v in enumerate(r)] for i,r in enumerate(A)]
    hi=max(sum(abs(v) for v in r) for r in B)+1;lo=0
    def test(t):
        return positive([[t*int(i==j)-v for j,v in enumerate(r)] for i,r in enumerate(B)])
    assert test(hi)
    while hi-lo>1:
        m=(lo+hi)//2
        if test(m):hi=m
        else:lo=m
    assert test(hi)
    return F(hi,scale)
def gram_int(rows,n):
    A=[[0]*n for _ in range(n)]
    for r in rows:
        for i in range(n):
            for j in range(i,n):A[i][j]+=r[i]*r[j]
    for i in range(n):
        for j in range(i):A[i][j]=A[j][i]
    return [[F(v,G*G) for v in r] for r in A]
def tail_kernel(g,J0):
    # Functions 1/j, log(j/J0)/j and cos(j theta_n)/j, one parity.
    # Main continuous Gram plus a diagonal upper reserve for all discrete
    # quadrature and Abel differences, uniform in the parity offset.
    theta=[g['NUF']*t for t in g['TS']]
    S=g['S'];m=2+len(theta)
    def sinlo(t):
        s=g['sincos'](t)[0]
        v=max(F(0),F(s.l,S),-F(s.h,S));assert v>0;return v
    B=[[F(0)]*m for _ in range(m)]
    R=[[F(0)]*m for _ in range(m)]
    B[0][0]=B[0][1]=B[1][0]=F(1,2*J0);B[1][1]=F(1,J0)
    for i in range(2):
        for j in range(2):R[i][j]=F(3,J0**2)
    for i,t in enumerate(theta,2):
        B[i][i]=F(1,4*J0)
        R[i][i]=(F(3)+1/sinlo(2*t))/J0**2
        for j in range(2):R[i][j]=R[j][i]=F(1,J0**2)/sinlo(t)
        for j in range(2,i):
            u=theta[j-2]
            R[i][j]=R[j][i]=(1/sinlo(t-u)+1/sinlo(t+u))/(2*J0**2)
    for i in range(m):B[i][i]+=sum(R[i])
    return B
def main(src,expected,out):
    raw=src.read_bytes();assert hashlib.sha256(raw).hexdigest()==expected
    d=json.loads(raw);assert d['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER'
    assert set(d['sectors'])=={'1','2'}
    g=namespace();J=2048;J0=J-1
    band=json.loads(Path(__file__).with_name('CR-ALP-068_pole_centered_gamma_result.json').read_text())['bands']
    assert band[0]['start']==86 and band[-1]['end'] is None
    assert all(a['end']+1==b['start'] for a,b in zip(band,band[1:]))
    def dg(b):return F((F(b['gamma_lower'])*10**9).__floor__(),10**9)
    cap=next(dg(b) for b in band if b['start']==2049);dc=cap-F(21,8)
    def di(j):return min(cap,dg(next(b for b in band if b['start']<=j and (b['end'] is None or j<=b['end']))))-F(21,8)
    B=tail_kernel(g,J0);bm=len(B);knorm=max(sum(abs(v) for v in r) for r in B)
    nu=F(g['NUF'].l,g['S']);assert nu>F(15,11)
    result={}
    for k,r in d['sectors'].items():
        n=43 if k=='1' else 42
        p=list(map(F,r['LDL_pivots']));assert len(p)==n and min(p)>0
        assert r['direct_prefix_overlap_checks']==3*n
        PG=int(r['prefix_grid']);assert PG>0
        assert [j for j,_ in r['prefix_mid_radius']]==list(range(85+(1 if k=='2' else 2),J+1,2))
        mid=[];err2=0
        for j,row in r['prefix_mid_radius']:
            assert len(row)==n and di(j)>0
            rr=[]
            for (m,a),v in zip(row,p):
                assert a>=0
                l,h=roots(v/di(j))
                x,y=interval_product(F(m-a,PG),F(m+a,PG),l,h)
                c,e=box(x,y);rr.append(c);err2+=e*e
            mid.append(rr)
        P=gram_int(mid,n);p0=spectral(P)
        pnew=(up(p0)+up(F(err2,G*G)))**2
        pold=F(r['eta_prefix_upper'])
        # Tail common basis: retain cross-column signs in each of nine
        # analytic leading functions, not independent columnwise norms.
        coeff=[[F(0)]*n for _ in range(bm)];ce2=0;rem2=F(0)
        for i,(v,t) in enumerate(zip(p,r['tail_receipts'])):
            pole=t['pole']
            seq=[pole['constant'],pole['log_coefficient']]+pole['cosine_coefficients']
            assert len(seq)==bm
            l,h=roots(v/dc)
            for j,(a,b) in enumerate(seq):
                x,y=interval_product(F(a),F(b),l,h)
                m,e=box(x,y);coeff[j][i]=F(m,G);ce2+=e*e
            rem=F(pole['pole_remainder_norm_upper'])+F(t['sine_series_norm_upper'])
            rem2+=v*rem*rem/dc
        BC=[[sum(B[i][j]*coeff[j][l] for j in range(bm)) for l in range(n)] for i in range(bm)]
        Q=[[sum(coeff[l][i]*BC[l][j] for l in range(bm))/(nu*nu) for j in range(n)] for i in range(n)]
        t0=spectral(Q)
        terr=up(knorm*F(ce2,G*G))/nu
        tnew=(up(t0)+terr+up(rem2))**2
        told=F(r['eta_tail_upper'])
        eta=min(pold,pnew)+min(told,tnew)
        E,U,b2=map(F,(r['E'],r['U'],r['finite_charge_upper']))
        assert E>0 and U>=0 and b2>=0
        loss=U+(up(b2)+up(eta*U))**2
        # For extremely small b2 or U, fixed grid root precision here is
        # 1e-30; the old 123AC receiver has a finer reserve. Keep best bound.
        loss=min(loss,F(r['loss_upper']))
        result[k]={'prefix_midpoint_gram':[[str(v) for v in row] for row in P],
            'prefix_midpoint_spectral_upper':str(p0),'prefix_error_frobenius2':str(F(err2,G*G)),
            'prefix_trace_upper':str(pold),'prefix_spectral_upper':str(pnew),
            'tail_leading_gram':[[str(v) for v in row] for row in Q],
            'tail_leading_spectral_upper':str(t0),'tail_coefficient_error_norm_upper':str(terr),
            'tail_remainder_frobenius2_upper':str(rem2),'tail_trace_upper':str(told),
            'tail_spectral_upper':str(tnew),'eta_upper':str(eta),
            'E':str(E),'U':str(U),'finite_charge_upper':str(b2),'loss_upper':str(loss),
            'conditional_pole_schur_lower':str(E-loss),
            'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL' if E>loss else 'SUFFICIENT_BOUND_INCONCLUSIVE'}
        print(k,'prefix trace/spectral',float(pold),float(pnew),'tail trace/spectral',float(told),float(tnew),'eta',float(eta),'loss/E',float(loss/E),result[k]['classification'],flush=True)
    out.write_text(json.dumps({'status':'CONDITIONAL_SPECTRAL_RESPONSE_DEFECT_RECEIVER','source_sha256':expected,
        'tail_kernel_upper':[[str(v) for v in row] for row in B],'sectors':result,
        'scope':'Inherited complete comparator and signed tail premises; no global support conclusion.'},indent=2)+'\n')
def selftest():
    for A,lam in [([[F(1),F(1)],[F(1),F(1)]],F(2)),([[F(1),F(0)],[F(0),F(4)]],F(4))]:
        b=spectral(A);assert lam<b<lam+F(1,10**10)
    assert box(F(-1,3),F(1,3))[1]>0
    print('spectral arithmetic synthetic preflight PASS',flush=True)
if __name__=='__main__':
    if sys.argv[1:]==['--selftest']:selftest()
    else:main(Path(sys.argv[1]),sys.argv[2],Path(sys.argv[3]))
