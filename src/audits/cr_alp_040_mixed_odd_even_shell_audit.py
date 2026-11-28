#!/usr/bin/env python3
"""CR-ALP-040: source-signed first shell mixed infinite-tail Schur tests.

Extends the CR-ALP-038 / 039 rank-one-plus-diagonal residual envelope to
w13,w14 (both have index <= 15000), retains both original boundary anchors,
reads the independent ten-million-index Gamma Jensen lower bank.
Tests the *even-index* physical 7-mode block plus Q15000 at margin 1e-14.
Also supplies an explicit negative test of the CONSERVATIVE odd-indexed
sufficient Schur comparison, NOT a negative original Weil direction.
All accepted comparisons use exact Fraction and integer Bareiss.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from hashlib import sha256
import json
from cr_alp_040_stable_14_mode_sylvester_audit import sym,det_bareiss,H,M,ROOT,S

A=F(11,10); N=15000; S_A=F(177,50)
SIGMA=F(1,10**14)
THETA=F(1,10**16)

# Independently reconstruct outward source input intervals from the pinned file.
params={}
for ln in (ROOT/'CR-ALP-040_initial_rational_intervals_1e25.txt').read_text().splitlines():
 key,lo,hi=ln.split()
 assert key not in params
 params[key]=(F(int(lo),S),F(int(hi),S))

def nu(j):return params[f'nu_{j}']
def lam(j):return params[f'lam_{j}']

def outer(z):
 lo=z.numerator*S//z.denominator
 hi=-((-z.numerator*S)//z.denominator)
 return lo,hi

def iproduct(a,b):
 v=[x*y for x in a for y in b]
 return min(v),max(v)

def compute_block(ids, anchor, ep):
 inner=ids[:-1]
 assert ids[-1]==anchor and len(inner)==6
    # D_j = nu_j - lambda_j nu_anchor, with strict outward enclosures.
 va=nu(anchor)
 Ds=[]; betas=[]
 for j in inner:
  v=nu(j);l=lam(j)
  z=iproduct(l,va)
  d=(v[0]-z[1],v[1]-z[0]);Ds.append(d)
  assert d[0]>0
  c0=1+l[1]
  c2=v[1]**2+l[1]*va[1]**2
  c3=v[1]**3+l[1]*va[1]**3
  b=(S_A*c0/(2*A)+8*S_A*c2/A+F(201,200)*c3/A)*(2*A/3)**3/F(isqrt(5*N**5))
  betas.append(b)
 alpha=(2*S_A/A+F(201,200)*F(22,7)/(2*A))*(2*A/3)**2/F(isqrt(3*N**3))
 bb=sum((b*b for b in betas),F(0))
 Z={}
 for r,i in enumerate(ids):
  for c,j in enumerate(ids):
   hlo,hhi=sym(H,i,j);mlo,mhi=sym(M,i,j)
   lo=F(hlo,S)-SIGMA*F(mhi,S)
   hi=F(hhi,S)-SIGMA*F(mlo,S)
   if r<6 and c<6:
    prod=iproduct(Ds[r],Ds[c])
    lo-=48*alpha**2*prod[1]
    hi-=48*alpha**2*prod[0]
    if r==c:
     lo-=48*bb
     hi-=48*bb
   if r==c==6:
    lo-=12*ep**2
    hi-=12*ep**2
   assert lo<=hi
   Z[r,c]=outer(lo)[0],outer(hi)[1]
 for r in range(7):
  for c in range(7):
   # Tight symmetric intersection of two independently calculated directions.
   lo=max(Z[r,c][0],Z[c,r][0]);hi=min(Z[r,c][1],Z[c,r][1]); assert lo<=hi
   Z[r,c]=(lo,hi)
 delta=max((sum((F(Z[r,c][1]-Z[r,c][0],2*S) for c in range(7)),F(0)) for r in range(7)))
 print('SOURCE MIXED', ids, 'max interval matrix row perturbation',float(delta),'alpha upper',float(alpha),'bb upper',float(bb))
 result={'indices':ids,'alpha_exact_rational_upper':str(alpha),'beta_squared_exact_upper':str(bb),
         'uncertainty_row_sum':str(delta),'sigma':str(SIGMA),'anchor_eps':str(ep)}
 # Stable congruence: Z_true lies within theta of its midpoint in operator norm.
 if delta<THETA:
  cm=[[Z[r,c][0]+Z[r,c][1]-(2*S//10**16 if r==c else 0) for c in range(7)] for r in range(7)]
  minors=[det_bareiss([row[:k] for row in cm[:k]]) for k in range(1,8)]
  result['positive_sylvester_minors']=sum(1 for d in minors if d>0)
  result['all_stable_minors_positive']=all(d>0 for d in minors)
  result['minors_sha256']=sha256(('|'.join(map(str,minors))).encode()).hexdigest()
  print('  Sylvester signs',[int(d>0)-int(d<0) for d in minors])
  if all(d>0 for d in minors):
   print('  PASS original source mixed sector + Q15000, sigma',SIGMA)
  else:
   # Construct a rational negative test of the sufficient *lower comparison*
   # using floating point only to SELECT a candidate; the final sign is checked
   # by exact integer arithmetic from outward interval entry bounds.
   # Deterministic candidate from a SEPARATE prior numerical exploration;
   # no floating-point is executed or relied upon in this exact certificate.
   z=[-658928912,698188837,-278385975,28637417,5087748,-546722,-1965]
   assert any(z)
   upper=0
   for r in range(7):
    for c in range(7):
     u=z[r]*z[c]
     upper+=u*(Z[r,c][1] if u>=0 else Z[r,c][0])
   if upper<0:
    print('  CERTIFIED: current conservative sufficient matrix is NOT PSD; NOT a negative original Weil witness')
    result['sufficient_comparison_negative_witness_integer_vector']=z
    result['sufficient_comparison_quadratic_upper_in_1e25']=str(upper)
    assert upper<0
   else:
    print('  insufficient to prove nonpositive conservative matrix at this precision')
 else:
  print('  interval error exceeds chosen theta; inconclusive')
 return result

# The even-index seven-dimensional block is the prior odd-function six-mode
# certificate enlarged by w14; the odd-index block is enlarged by w13.
res=[compute_block([4,6,8,10,12,14,2],2,F(1,19)),
     compute_block([3,5,7,9,11,13,1],1,F(1,25))]
assert res[0]['all_stable_minors_positive'], 'failed new even-indexed mixed test'
assert res[1].get('sufficient_comparison_quadratic_upper_in_1e25'), 'expected documented conservative obstacle'
result={
 'stage':'CR-ALP-040',
 'mixed_status':'PASS even-indexed physical E14 parity block + its genuine infinite Q15000 parity tail; odd-indexed enlarged conservative comparison FAILS as a sufficient test',
 'consequence':'Source-signed positivity on (E12 + span{w14}) ⊕ Q15000, with sigma=1e-14 on enlarged even-indexed sector, 1e-12 on retained odd-indexed sector, and 1/8 on both high parity tails',
 'supporting_analytic_premises':'CR-ALP-027 1/4 positive tail, CR-ALP-038 rank-one residue estimate, CR-ALP-039 boundary-anchor eps; new original signed Jensen lower H14',
 'warning':'The negative conservative odd-indexed sufficient comparison is NOT a negative complete original Weil witness; cross-terms bounded too loosely',
 'blocks':res,
 'not_claimed':'E14 ⊕ entire Q15000 mixed positivity; all 14988 intermediate shell directions; original collar positivity; RH; Lean'
}
(ROOT/'CR-ALP-040_mixed_shell_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('CR-ALP-040 MIXED SHELL EXACT INTEGER + ANALYTICAL-PREMISE AUDIT PASS')
