#!/usr/bin/env python3
"""CR-ALP-123AM: high-frequency Gamma kernel tightening of admitted AL.
No new moments, inverse solves, or response bank. Exact rational arithmetic.
Accept raw JSON or ZIP replay wrappers with result.json.
"""
import hashlib,json,sys,zipfile
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
sys.set_int_max_str_digits(0)
AC='a7a334e6eca63543d7f322aef124c4491c986ba878bd1434f484e86743856473'
AL='b3d984dd48dcae906a8582afbffc03b047d48275c780534086cbd2d671a87e2e'
def read(p,sha):
    if zipfile.is_zipfile(p):
        with zipfile.ZipFile(p) as z:
            assert z.testzip() is None
            b=z.read('result.json')
    else:b=Path(p).read_bytes()
    assert hashlib.sha256(b).hexdigest()==sha
    return json.loads(b)
def ceil(x):return -((-x).__floor__())
def up(x,g=10**9):return F(ceil(x*g),g)
def root(x):
    assert x>=0
    G=10**110;v=x*G*G
    return F(isqrt(v.numerator//v.denominator)+1,G)
def main(ac,al,out):
    a=read(ac,AC);l=read(al,AL)
    assert a['status']=='CONDITIONAL_RESPONSE_COLUMN_DEFECT_RECEIVER'
    assert l['status']=='CONDITIONAL_METRIC_FIRST_SIGNED_MOMENT_RECEIVER'
    assert set(a['sectors'])==set(l['sectors'])=={'1','2'}
    # beta_m=1/2+2m. For decreasing f(x)=v/(x^2+v^2),
    # sum f(beta_m)<=1/v+(1/2)*integral_0^infty f=1/v+pi/4.
    # Reflected exponential sum <= (2/5)/(1-(2/5)^4)/v <1/v.
    # The exact 046 Gamma-sine coefficient is their signed difference.
    V=F(15*2049,11);exp_sum=F(2,5)/(1-F(2,5)**4)
    assert exp_sum<1
    gamma_upper=F(11,14)+(1+exp_sum)/V # pi<22/7
    assert gamma_upper<F(4,5)
    # AL used |A_j^G|<=4. On j>2048 all 32 explicit gamma_sine
    # moment kernel norm bounds may be multiplied by (4/5)/4=1/5.
    # The entire unknown remainder is retained unchanged.
    factor=F(1,5);dc=F('1.822938456');ans={}
    for k,s in l['sectors'].items():
        z=a['sectors'][k];assert s['dimension']==44-int(k)
        assert s['classification']=='IMPROVED_METRIC_SUFFICIENT_BOUND_INCONCLUSIVE'
        rows=s['moment_norm_uppers'];labels=['centered_combined','gamma_sine','prime_sine']
        assert [(x,i) for x,i,v in rows]==[(x,i) for x in labels for i in range(32)]
        assert all(F(v)>=0 for x,i,v in rows)
        sums={x:sum((F(v) for y,i,v in rows if y==x),F(0)) for x in labels}
        assert sum(sums.values())>=F(s['series_norm_upper'])
        series=sums[labels[0]]+factor*sums[labels[1]]+sums[labels[2]]
        lead=F(s['leading_tail_trace_upper']);rem=F(s['remainder_norm_upper']);p=F(s['prefix_relative_upper'])
        assert min(lead,rem,p)>=0
        tail=(root(lead)+series+rem)**2/dc
        rho=up(p+tail,10**6);assert 0<rho<1
        eta=rho/(1-rho)
        E,U,b2=map(F,(z['E'],z['U'],z['finite_charge_upper']));assert E>U>0 and b2>=0
        loss=U+(root(b2)+root(eta*U))**2
        lossratio=up(loss/E);margin=1-lossratio
        assert margin>0
        ans[k]={'dimension':s['dimension'],'old_loss_over_E_upper':s['loss_over_E_upper'],
          'old_moment_sum_uppers':{x:str(v) for x,v in sums.items()},
          'gamma_sine_rescaled_norm_upper':str(factor*sums['gamma_sine']),
          'series_norm_upper':str(series),'leading_tail_trace_upper':str(lead),'remainder_norm_upper':str(rem),
          'prefix_relative_upper':str(p),'tail_relative_upper':str(up(tail)),
          'total_relative_upper':str(rho),'schur_fraction_lower':str(1-rho),'eta_upper':str(eta),
          'loss_over_E_upper':str(lossratio),'pole_schur_over_E_lower':str(margin),
          'classification':'PASS_COMPLETE_POLE_SECTOR_CONDITIONAL'}
    result={'status':'CONDITIONAL_HIGH_FREQUENCY_GAMMA_TIGHTENING','source_ac_sha256':AC,'source_al_sha256':AL,
      'gamma_coefficient_upper':str(gamma_upper),'used_gamma_coefficient_upper':'4/5','gamma_exp_sum_upper':str(exp_sum),
      'explicit_gamma_moment_rescale':str(factor),'sectors':ans,
      'scope':'Inherited 122R/123AC/123AK/123AL physical and interval contracts; fixed collar only. Entire old unknown remainder retained; no new solves.'}
    Path(out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main(*sys.argv[1:])
