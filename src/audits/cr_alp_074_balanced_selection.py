#!/usr/bin/env python3
"""CR-ALP-074: balanced parity rank-20 receiver against original L.

Source-only candidate. For each parity use K_L - r_L*L^-1 r_L,
with L >= 51/200 I. There is no comparator change A and no cross-parity
interaction. Exact physical residual certificates remain required.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json,sys
ELL=F(51,200)
S=10**12
HASHES={
 '1':'5340455da04a29619f943a4360caed34a1b9bdda56d68068b9f9c523af784325',
 '2':'4017c7a656e43fe9fa97da91c37f35a9eda5a78dd3efaa2863b40de71177ac36'}
def checked(p,h):
    raw=p.read_bytes();assert hashlib.sha256(raw).hexdigest()==h,p
    return json.loads(raw)
def absmax(q):return max(map(abs,q))
def rootdown(x):
    assert x>=0;y=x*S*S
    return F(isqrt(y.numerator//y.denominator),S)
def select(w,k):
    K=w['K'];scale=int(w['scale']);n=len(K)
    assert n== (43 if k=='1' else 42)
    assert all(len(row)==n for row in K)
    M=[]
    for i in range(n):
        row=[]
        for j in range(n):
            l,h=map(int,K[i][j]);ll,hh=map(int,K[j][i])
            l,h=max(l,ll),min(h,hh);assert l<=h
            row.append((F(l,scale),F(h,scale)))
        M.append(row)
    def lam(J):
        return min(M[i][i][0]-sum(absmax(M[i][j]) for j in J if i!=j) for i in J)
    J=[]
    for _ in range(10):
        J.append(max((lam(J+[i]),-i,i) for i in range(n) if i not in J)[2])
    assert len(set(J))==10
    lower=lam(J)
    # Reserve half of the Gershgorin energy for all residual penalties.
    q=rootdown(ELL*lower/20) if lower>0 else F(0)
    margin=lower-10*q*q/ELL
    return {'parity':k,'indices_zero_based':J,'rank':len(J),
      'K_lower':str(lower),'common_actual_residual_norm_ceiling':str(q),
      'inverse_penalty_upper_at_ceiling':str(10*q*q/ELL),
      'Schur_lower_at_ceiling':str(margin),
      'conditional_gate_usable':q>0 and margin>0}
def main(p,out):
    results={k:select(checked(p/'cr_alp_056'/f'result_extended_{k}.json',h),k) for k,h in HASHES.items()}
    r={'status':'SOURCE_SELECTION_ONLY','rank':sum(x['rank'] for x in results.values()),
       'comparator':'L=(72/25)I+Q85 P Q85','coercivity_lower':str(ELL),
       'sectors':results,'conditional_gate_usable':all(x['conditional_gate_usable'] for x in results.values()),
       'remaining':'Actual projected original residual norms for every selected coordinate; source and numerical attachment; inherited analytical premises.',
       'index_if_received':65,'current_index_cap':66}
    out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
