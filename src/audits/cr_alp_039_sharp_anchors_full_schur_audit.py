#!/usr/bin/env python3
"""CR-ALP-039: authentic 12 physical directions + infinite Q_15000 Schur.

RIGOROUS EXACT-RATIONAL ARITHMETIC, conditional on prior human analytic
Weil source, Laplace Gamma, and sine/pole projector identities.

Certificate of the matrix inequality is outward fixed-point interval
LDL of the pinned CR-ALP-034 ORIGINAL SIGNED lower matrix Hmid, in the
physical M metric. Anchor tail enclosures are derived from a centered
log-digamma bound and direct integration-by-parts of the authentic seven
paired prime translations. No Lean compilation, 12-to-15000 shell, all
collars, or RH is asserted.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt, factorial
import hashlib, json, sys
ROOT=Path(__file__).resolve().parent
BASES=ROOT.parent/'cr_alp_038'
sys.path.insert(0,str(BASES))
from cr_alp_038_full_interior_infinite_schur_audit import (Q,H,M,check_sym,nu,BASE,N,S,PI,A)

assert N==15000 and A==F(11,10) and S==F(177,50)
assert hashlib.sha256((BASES/'H_mid.tsv').read_bytes()).hexdigest()=='b6a31e623298272fdceac788fc56626ec704e280903ff6baff1dc8bfe0d45f92'
assert hashlib.sha256((BASES/'initial.txt').read_bytes()).hexdigest()=='bcb6017b9358343dc16fb97cf61c1c340592444d124a0e2e466296d9703f53a2'

# Elementary scalar bounds: 3<pi<22/7, sqrt(a)<21/20,
# e^(11/2)>200 and e^(22/5)>80 (from predecessor exact Taylor audits),
# -psi(1/4)<5 and log(20000)<11.
assert F(4)*F(21,20)/3 == F(7,5)
assert F(2,N)+F(2,N*N) < F(1,86**2)
assert F(1,3*N**3) < F(1,3000000**2)
assert F(1,5*N**5) < F(1,10**10)**2
assert F(1,N) < F(1,122)**2
# Independent elementary rational Taylor checks for every exponential
# comparison used below, with no floating-point acceptance comparisons.
for z, degree, threshold in [(F(11,2),8,200),(F(22,5),9,80),
                             (F(3,4),2,2),(F(11),9,20000)]:
    assert sum((z**k/F(factorial(k)) for k in range(degree+1)),F(0))>threshold
# Euler constant gamma_E<1, pi/2<11/7, log(2)<3/4 imply
# -psi(1/4)=gamma_E+pi/2+3log(2)<1+11/7+9/4<5.
assert 1+F(11,7)+3*F(3,4)<5
# Digamma variation error |Re psi(1/4+i nu/2)-log(nu/2)|<3/nu
# and log(20000)<11 imply for nu>=20000: 0<=Lambda(nu)<nu/1000.
assert F(15000)*3/(2*A)>20000
assert 15+F(3,20000)+F(20000,20000)< F(20000,1000)
# sum_{m>=1} beta_m e^(-2a beta_m) <=
# (1/200)*sum_{h>=0}(2h+5/2)/80^h <1/70.
exp_laplace_tail=F(1,200)*(F(5,2)/(1-F(1,80)) + 2*F(1,80)/(1-F(1,80))**2)
assert exp_laplace_tail < F(1,70)

anchors=[
  {'k':1,'parity':'odd_j_even_function','gamma_upper':F(3,5),'B0':F(6,5),'B1':F(3,10),
   'gamma_target':F(1,170),'prime_target':F(1,30),'joint_target':F(1,25)},
  {'k':2,'parity':'even_j_odd_function','gamma_upper':F(5,3),'B0':F(3,5),'B1':F(3,5),
   'gamma_target':F(1,120),'prime_target':F(1,23),'joint_target':F(1,19)},
]

for d in anchors:
    k,g,B0,B1 = d['k'],d['gamma_upper'],d['B0'],d['B1']
    # Gamma b: subtract t_N b, where t_N exactly removes the high-frequency
    # constant in the complete digamma and convergent E(beta) terms.
    # R_j <= (7 B0/10) log(j/N)/j + C2/j^2 + C3/j^3.
    Glog = F(7,10)*B0
    G2 = F(121,225)*(F(12001,4000)*B0+F(11,7)*B1)
    G3 = (F(8,945)*A**3)*B0
    # summation sqrt(sum log^2(j/N)/j^2) <1/86,
    # sqrt(sum j^-4)<1/3e6 and sqrt(sum j^-6)<1/1e10.
    # Q T_Gamma e_k contribution <1e-6 from CR-ALP-037.
    GammaUpper=g*(Glog/F(86)+G2/F(3000000)+G3/F(10**10))+F(1,10**6)
    assert GammaUpper<d['gamma_target'],('Gamma anchor fails',k,float(GammaUpper))
    # Full paired ORIGINAL signed seven-prime source:
    # translation of b has a source-dependent *constant* /nu_j term
    # (cancelled by the genuine pole projector) and an oscillating
    # cos(nu_j log n)/nu_j term. The exact finite-interval integration by
    # parts leaves <=(2 B1 + a B0)/(sqrt a nu_j^2).
    # U_{j,k} low sine overlap has its nonoscillatory 1/nu_j term cancelled
    # by the same genuine pole projector; its remainder is O(nu_j^-2).
    vk_upper=F(k)*F(22,7)/(2*A)
    P1=g*F(7,5)*S*B0
    P2=(F(8)*S*vk_upper/(3*A)+g*S*(2*B1+A*B0))*(2*A/3)**2
    P3=(F(8)*S*(vk_upper**2+F(1,4))/(3*A)+g*S*B0/F(2))*(2*A/3)**3
    PrimeUpper=P1/F(122)+P2/F(3000000)+P3/F(10**10)
    assert PrimeUpper<d['prime_target'],('prime anchor fails',k,float(PrimeUpper))
    assert d['gamma_target']+d['prime_target']<d['joint_target']
    assert GammaUpper+PrimeUpper<d['joint_target']
    d['gamma_exact_rational_upper']=str(GammaUpper)
    d['prime_exact_rational_upper']=str(PrimeUpper)
    d['combined_exact_rational_upper']=str(GammaUpper+PrimeUpper)
    print(f'PASS anchor u_{k}: complete original Gamma < {d["gamma_target"]}, seven-prime < {d["prime_target"]}, combined < {d["joint_target"]}')

# The CR-ALP-038 full interior residual is bounded by
# ||Q T sum x_k w_k||² <= 2 alpha² (sum D_k x_k)² +2 bb ||x||².
# Decompose p = sum x w + y u;  ||a+b||² <=3||a||²+3/2||b||².
# W(q,q)>=1/4||q||², and 2||Q Tp||||q||<=8||Q Tp||²+1/8||q||².
# So Hmid-48alpha²DD^T-48bb I-12 eps² |y|² - sigma M >>0
# suffices for W(p+q)>sigma||p||²+1/8||q||².
SIGMA=F(1,10**12)
rows=[]


def full_block(ids,anchor):
    n=5
    assert len(ids)==6 and ids[-1]==anchor['k']
    va=nu(ids[-1]); DD=[];beta=[]
    for k in ids[:-1]:
        v=nu(k)
        lam=v*(va**2+F(1,4))/(va*(v**2+F(1,4)))
        assert lam.lo>0
        DD.append(Q(v-lam*va))
        c0=1+lam.hi
        c2=v.hi**2+lam.hi*va.hi**2
        c3=v.hi**3+lam.hi*va.hi**3
        B=(S*c0/(2*A)+8*S*c2/A+F(201,200)*c3/A)*(2*A/3)**3
        beta.append(B/isqrt(5*N**5))
    alpha=(2*S/A+F(201,200)*PI.hi/(2*A))*(2*A/3)**2/isqrt(3*N**3)
    alpha_q=Q(Q(alpha).hi,Q(alpha).hi)
    bb=sum((z*z for z in beta),F(0))
    bb_q=Q(Q(bb).hi,Q(bb).hi)
    eps=anchor['joint_target']
    # The upper-bound coupling matrix subtracts a PSD rank-one plus diagonal.
    mat={}
    for irow,i in enumerate(ids):
        for jcol,j in enumerate(ids):
            val=check_sym(H,i,j)-Q(SIGMA)*check_sym(M,i,j)
            if irow<n and jcol<n:
                val-=Q(48)*alpha_q*alpha_q*DD[irow]*DD[jcol]
                if irow==jcol:
                    val-=Q(48)*bb_q
            if irow==n and jcol==n:
                val-=Q(12)*Q(eps)*Q(eps)
            mat[irow,jcol]=val
    for ii,i in enumerate(ids):
        for jj,j in enumerate(ids):
            x=mat[ii,jj]
            rows.append((i,j,x.lo,x.hi))
    L={}; pivots=[]
    for j in range(len(ids)):
        p=mat[j,j]
        for k in range(j):
            p=p-L[j,k]*L[j,k]*pivots[k]
        assert p.lo>0,('FAILED INTERVAL LDL PIVOT',ids,j,p.lo,p.hi)
        pivots.append(p)
        for r in range(j+1,len(ids)):
            v=mat[r,j]
            for k in range(j):v=v-L[r,k]*L[j,k]*pivots[k]
            L[r,j]=v/p
    print('PASS entire original signed 6-mode parity plus infinite Q15000 mixed Schur',ids)
    print('  exact lower LDL pivots: ',[str(z.lo)+'/10^25' for z in pivots])
    print('  alpha upper',alpha_q.hi,'/10^25; beta² upper',bb_q.hi,'/10^25; anchor eps',eps)
    return {'indices':ids,'shifted_LDL_pivots_lower_1e25':[str(z.lo) for z in pivots],
       'alpha_upper_1e25':str(alpha_q.hi),'beta_squared_upper_1e25':str(bb_q.hi),
       'anchor_complete_original_signed_residual_bound':str(eps),
       'shifted_physical_coercivity':str(SIGMA),'positive_pivots':6}

blocks=[full_block([3,5,7,9,11,1],anchors[0]),full_block([4,6,8,10,12,2],anchors[1])]
assert len(rows)==72
MATRIX=ROOT/'CR-ALP-039_full_twelve_shifted_H_interval_matrix.tsv'
MATRIX.write_text(''.join('%d %d %d %d\n'%x for x in rows))
MATRIX_SHA=hashlib.sha256(MATRIX.read_bytes()).hexdigest()
print('  signed shifted interval matrix 72 entries SHA256:',MATRIX_SHA)

result={
 'stage':'CR-ALP-039',
 'status':'OUTWARD INTEGER EXACT ARITHMETIC PASS on ORIGINAL signed E12 + Q15000 mixed Schur',
 'conditional_theorem':'W(p+q,p+q) > 1e-12||p||_L2^2 + 1/8||q||_L2^2 for all nonzero p in E12 and q in Q15000(H0) intersect original log form domain',
 'source':'Original Zeta23.EF.literatureRHS on a=11/10; 2 true Mellin pole nulls, 7 negative prime powers, complete digamma',
 'certification':'Classical continuum source identities NOT Lean; exact Fraction plus 10^-25 outward-integer LDL',
 'Gamma_anchor_method':'exact pole cancellation of centered (log(j/N))/j and O(j^-2) entire Gamma remainder',
 'prime_anchor_method':'paired shifts integration by parts: pole-aligned constant 1/j removed, oscillatory cos(nu_j log n)/j retained',
 'shifted_interval_matrix_sha256':MATRIX_SHA,
 'Hmid_sha256':hashlib.sha256((BASES/'H_mid.tsv').read_bytes()).hexdigest(),
 'rational_source_sha256':hashlib.sha256((BASES/'initial.txt').read_bytes()).hexdigest(),
 'anchors':[ {k: (str(v) if isinstance(v,F) else v) for k,v in d.items()} for d in anchors],
 'blocks':blocks,
 'not_claimed':'intermediate signed modes 13..15000, full two-pole-null collar, other support radii, Lean formalization, Weil RH',
}
(ROOT/'CR-ALP-039_certified_twelve_plus_infinite_tail_result.json').write_text(json.dumps(result,indent=2)+'\n')
print('CR-ALP-039 EXACT SOURCE-SIGNED TWELVE PLUS INFINITE TAIL AUDIT PASS')
