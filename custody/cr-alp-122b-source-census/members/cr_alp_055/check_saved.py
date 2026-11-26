#!/usr/bin/env python3
"""Exact independent saved-matrix tests; no continuum coefficient reconstruction."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parent
# Reuse only elementary exact verifier functions, not prior test execution.
src=(ROOT.parent/'cr_alp_054'/'check_saved.py').read_text().split('result={}')[0]
# Define functions in this namespace, retaining this audit's path.
exec(src.replace('ROOT=Path(__file__).resolve().parent',''))
NEWBANK={1:'1ec38f55beb2d59e1467be94908fe1ca33d45d24e7cfab310084775e65f34bdb',2:'47a10b15f4bdd68e5b84fb98303455df3bd3b638eeed06c76b2f88bfc184a373'}
NEWMOM={1:'261bc15cede57ca4d95fe04c550e2e6f7c74d3a7c1d40c664b2e94ff6d1d0d42',2:'904a304e810afce13f627f83203516dade09ec97665830e2071a94041ab0ceee'}
result={}
for k in (1,2):
 raw=(ROOT/f'bank_{k}.json').read_bytes(); b=json.loads(raw); mraw=(ROOT/f'moments_{k}.json').read_bytes(); m=json.loads(mraw)
 sha=hashlib.sha256(raw).hexdigest(); assert sha==NEWBANK[k] and hashlib.sha256(mraw).hexdigest()==NEWMOM[k]; assert m['input_sha256']==sha and b['predecessor_sha256']==BANK[k]
 oldraw=(ROOT.parent/'cr_alp_054'/f'bank_{k}.json').read_bytes(); assert hashlib.sha256(oldraw).hexdigest()==BANK[k]; ob=json.loads(oldraw)
 n=44-k; assert len(b['R'])==len(ob['R'])+1 and len(b['P'])==len(ob['P'])+1
 assert b['R'][:-1]==ob['R']; assert [row[:-1] for row in b['P'][:-1]]==ob['P']; assert b['R'][-1:]==b['rbridge']
 # Independently calculate norm of the single unit-bridge residual row using exact endpoints.
 scale0=int(b['scale']); lo=hi=F(0)
 for ll,hh in b['rbridge'][0]:
  l,h=F(int(ll),scale0),F(int(hh),scale0);lo+=0 if l<=0<=h else min(l*l,h*h);hi+=max(l*l,h*h)
 target=F(1,250) if k==1 else F(1,60); print(k,'global bridge squared norm',float(lo),float(hi),'prior tolerance squared',float(target**2),flush=True)
 S=int(m['scale']);K,G,M=map(load,[m['K'],m['G'],m['M']]);assert len(K)==n
 A=[[(floor((gamma*gl-mh)/(ell*u)),ceil((gamma*gh-ml)/(ell*u))) for (gl,gh),(ml,mh) in zip(gr,mr)] for gr,mr in zip(G,M)]
 tests={}
 for name,t in m['tests'].items():
  J=A if name=='affine' else load(t['penalty']); C=sub(K,J); sh=F(t['shift']); status='inconclusive'; witness=None
  if t['positive_pivots']==n:
   assert ld(C,sh,S)==n; status='PASS'
  elif t['rejection_witness']:
   z=list(map(int,t['rejection_witness']['z'])); l,h=quad(C,z,S);witness=str(h)
   if h<0: status='FALSIFIED'
  tests[name]={'status':status,'shift':str(sh),'negative_witness_upper':witness,'penalty_row_upper':str(upper(J,S))}
  print(k,name,status,'shift',sh,'negative witness upper',float(F(witness)) if witness else None,flush=True)
 assert lo>target**2
 assert tests['chebyshev_12']['status']=='PASS' and F(tests['chebyshev_12']['shift'])==F(1,16)
 if k==2:assert tests['affine']['status']=='FALSIFIED' and F(tests['affine']['negative_witness_upper'])<-F(1049,100000)
 # Conditional complete inverse residual after the degree-12 correction on full E850 minus E85.
 assert upper(G,S)<F(1,5)
 x=F(192,175);t0,t1=F(1),x
 for j in range(2,13):t0,t1=t1,2*x*t1-t0
 assert 1/t1**2<F(1,8000)
 remaining=F(1,16)-F(200,51)*(F(1,40000)+F(1,100))
 assert remaining>F(1,50)
 print(k,'conditional Q850 norm <= 1/10; margin >',str(remaining),flush=True)
 result[str(k)]={'conditional_Q850_norm':'1/10','conditional_schur_margin_lower':str(remaining),'residual_gram_row_upper':str(upper(G,S)),'bank_sha256':sha,'moments_sha256':hashlib.sha256(mraw).hexdigest(),'global_bridge_squared_norm_lower':str(lo),'global_bridge_squared_norm_upper':str(hi),'prior_target_squared':str(target**2),'prior_target_falsified':lo>target**2,'tests':tests}
(ROOT/'verification_result.json').write_text(json.dumps(result,indent=2)+'\n')
