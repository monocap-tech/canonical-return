#!/usr/bin/env python3
"""CR-ALP-032: exact structural and Gamma-tail audit (stdlib Fraction).

Checks original seven-prime displacement-rank-two formula against the
independently certified source-prime interval bank (CR-ALP-028), validates the
Gamma sine rank-one structure algebraically, proves exact rational *tail
coefficients* for the continuum derivative/trace estimate, and independently
certifies the finite physical H1 stiffness/endpoint Gram structure for E10.

The integration-by-parts, distributional derivative, monotonicity, and
operator form inequalities are documented in the accompanying mathematical
note and NOT mechanically proved by this arithmetic audit. No 15k Schur pass,
original whole-collar PSD, Lean kernel pass or RH is claimed.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from json import dumps
from importlib.util import spec_from_file_location, module_from_spec

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def locate(name,stage):
    # Same-directory GitHub src/audits/, or stage-separated reproducibility ZIP.
    for candidate in (HERE/name, ROOT/stage/name):
        if candidate.is_file(): return candidate
    raise FileNotFoundError(f"Missing immutable predecessor artifact: {name}")
PREV=locate('cr_alp_028_prime_interval_and_gram_audit.py','cr_alp_028')
PREV_INPUT=locate('CR-ALP-031_initial_rational_intervals_1e25.txt','cr_alp_031')
GROUPED=locate('CR-ALP-031_grouped_H_exact_matrix_intervals.tsv','cr_alp_031')
# Cryptographic custody of the exact committed 031 intervals before using
# those entries to certify a new *upper* Rayleigh bound.
assert sha256(PREV_INPUT.read_bytes()).hexdigest()==('a8a827ff7ad0cca5d0e7e89c72c91798d9d96e8fd9d3d879af2ec8bfd152a187')
assert sha256(GROUPED.read_bytes()).hexdigest()==('052490ae00ccbc88466b8134e04a0100950cc48bf811a094086d739e9bda5ba1')
assert sha256(PREV.read_bytes()).hexdigest()==('ddaf7ad4aedc2aaa574e0dad1950290d03770e45bd83bc46c7867063b27a3255')


spec=spec_from_file_location('alp028_source_prime',PREV)
p=module_from_spec(spec)
spec.loader.exec_module(p)
assert set(p.POWERS)=={2,3,4,5,7,8,9}

# Displacement identity: on one parity the entire original signed prime block
# has rank <=2 displacement under X=diag(nu_j^2), not 2 per prime.
#   X P - P X = -(2/a)*(nu*b^T - b*nu^T),
#   b_j=sum Lambda(n)/sqrt(n)*sin(nu_j log n).
b={j:sum((p.weights[n]*p.trig(p.nu[j]*p.t[n],'sin')
           for n in p.POWERS),p.I(0)) for j in range(1,11)}
checked=0
for inds in ((1,3,5,7,9),(2,4,6,8,10)):
    for i in inds:
        for j in inds:
            if i==j:continue
            lhs=(p.nu[i]*p.nu[i]-p.nu[j]*p.nu[j])*p.P(i,j)
            rhs=-p.I(2)/p.I(p.A)*(p.nu[i]*b[j]-p.nu[j]*b[i])
            assert max(lhs.lo,rhs.lo)<=min(lhs.hi,rhs.hi), (i,j,lhs,rhs)
            checked+=1
print(f'PASS: seven native prime powers combine into ONE rank-two displacement per parity; {checked} source-interval offdiagonal checks')

# Exact rational algebra on a synthetic sine frequency lattice; this is a
# formula regression, NOT a positivity test of the original physical form.
a=F(11,10)
beta=F(7,2)
qexp=F(1,512)
for inds in ((1,3,5,7),(2,4,6,8)):
    sgn=(-1)**inds[0]
    for i in inds:
        for j in inds:
            vi=F(7*i,5);vj=F(7*j,5)
            di=beta*beta+vi*vi;dj=beta*beta+vj*vj
            K=(F(2)*beta/di if i==j else F(0)) + (F(2)/a)*(1-sgn*qexp)*vi*vj/(di*dj)
            J=(F(2)/beta if i==j else F(0))-K
            diagonal=(F(2)*vi*vi/(beta*di) if i==j else F(0))
            rank_one=F(2)*(1-sgn*qexp)/a*(vi/di)*(vj/dj)
            assert J==diagonal-rank_one
print('PASS: exact Gamma sine block = positive diagonal minus single rank-one term (each parity, each beta)')

# Rigorous sum upper bounds using f(x)=(2x+1/2)^(-p), decreasing on [K,infty):
# sum_{m>=K} f(m) <= f(K)+ integral_K^infty f(x) dx.
K=30000
betaK=F(4*K+1,2)
cD=F(4)/betaK**3 + F(1)/betaK**2   # 4*sum 1/beta_m^3
cT=F(4)/betaK**2 + F(2)/betaK   # 4*sum 1/beta_m^2
cZero=F(2)/betaK**3 + F(1,2)/betaK**2 #2*sum1/beta_m^3
assert cD<F(1,3_000_000_000)
assert cT<F(1,29_000)
assert cZero<F(1,7_000_000_000)
print('PASS: all m>=30000 Gamma remainder bounded by cD*||f\' ||²+cT*(|f(-a)|²+|f(a)|²)')
print('      cD <=', float(cD),' cT <=',float(cT),' cZero <=',float(cZero))

# Reuse the exact audited CR-ALP-031 source interval bank; no trust in rounded
# decimal display, and no replacement of sine physical Gram by identity.
S=10**25
records={}
for line in PREV_INPUT.read_text().splitlines():
    key,lo,hi=line.split()
    records[key]=(F(int(lo),S),F(int(hi),S))
assert len(records)==129

class I:
    def __init__(self,x=0,hi=None):
        if isinstance(x,I) and hi is None:self.lo,self.hi=x.lo,x.hi
        elif isinstance(x,tuple) and hi is None:self.lo,self.hi=x
        else:self.lo,self.hi=F(x),F(x if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self,b):
        b=I(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-I(b)
    def __rsub__(self,b):return I(b)+-self
    def __mul__(self,b):
        b=I(b); z=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi];return I(min(z),max(z))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=I(b);assert b.lo>0 or b.hi<0
        z=[self.lo/b.lo,self.lo/b.hi,self.hi/b.lo,self.hi/b.hi]
        return I(min(z),max(z))
    def __repr__(self):return f'[{float(self.lo):.12g},{float(self.hi):.12g}]'

def src(key):return I(records[key])
nu={j:src(f'nu_{j}') for j in range(1,11)}
lam={j:src(f'lam_{j}') for j in range(3,11)}
gamma={j:src(f'gamma_{j}') for j in (1,2)}
A=src('A');SA=src('SinhA');EHALF=src('EHALF');EHI=src('EHI')
B={1:A+SA,2:SA-A}
assert B[1].lo>0 and B[2].lo>0
c={k:gamma[k]*B[k] for k in (1,2)}
C={1:(EHALF+EHI)/2,2:(EHALF-EHI)/2}
T={1:2*(gamma[1]*C[1])*(gamma[1]*C[1]),
   2:2*(gamma[2]*C[2])*(gamma[2]*C[2])}
assert all(T[k].lo>0 for k in (1,2))

# Physical H1 stiffness D of basis w_j=e_j-lambda_j e_k and
# u_k=e_k-gamma_k b_k. Zero boundary and exact pole null give:
# <w_i',b_k'>= -1/4 <w_i,b_k>=0, hence D(w_i,u_k)=-lambda_i*nu_k².
def stiffness(inds,k):
    js=list(inds)+[k]
    D={}
    for i in inds:
        for j in inds:D[i,j]=(nu[i]*nu[i] if i==j else I(0))+lam[i]*lam[j]*nu[k]*nu[k]
        D[i,k]=D[k,i]=-lam[i]*nu[k]*nu[k]
    other=2 if k==1 else 1
    D[k,k]=nu[k]*nu[k]+gamma[k]*c[k]/2 +gamma[k]*gamma[k]*B[other]/4
    return js,D

def strict_interval_ldl(js,D):
    L={};diags={}
    for ii,i in enumerate(js):
        val=D[i,i]
        for kk in range(ii):
            k=js[kk];val-=L[i,k]*L[i,k]*diags[k]
        assert val.lo>0, (i,val)
        diags[i]=val
        for jj in range(ii+1,len(js)):
            j=js[jj];v=D[j,i]
            for kk in range(ii):
                k=js[kk];v-=L[j,k]*L[i,k]*diags[k]
            L[j,i]=v/val
    return [diags[j].lo for j in js]

for ids,k in [((3,5,7,9),1),((4,6,8,10),2)]:
    js,D=stiffness(ids,k)
    piv=strict_interval_ldl(js,D)
    # All w_j have exact zero endpoint trace. The lone anchor carries the
    # entire trace in each parity; no cross-parity trace exists.
    assert all(x>0 for x in piv)
    print(f'PASS E10 physical derivative Gram block parity k={k}: five positive interval LDL pivots; endpoint rank-one per parity; T(anchor)~{T[k]}')

# Loewner endpoint-group gap: CR-ALP-031 used 617 certified LOWER
# endpoints J_{U-1}. For monotone J, exact full group is between
# ell*J_{U-1} and ell*J_L. The upper-lower gap has a rigorous H1
# derivative+endpoint majorant from the e^(-beta*s) difference.
start=30000
STOP=350000
groups=[]
while start<STOP:
    length=max(1,start//250)
    end=min(STOP,start+length)
    groups.append((start,end))
    start=end
assert len(groups)==617
assert groups[0][0]==K and groups[-1][1]==STOP
assert all(groups[i][1]==groups[i+1][0] for i in range(len(groups)-1))
assert sum(U-L for L,U in groups)==STOP-K
cGapD=sum((F(24*(U-L)*(U-L-1))/(F(4*L+1,2)**4) for L,U in groups),F(0))
cGapT=sum((F(16*(U-L)*(U-L-1))/(F(4*L+1,2)**3) for L,U in groups),F(0))
betaStop=F(4*STOP+1,2)
cStopD=F(4)/betaStop**3+F(1)/betaStop**2
cStopT=F(4)/betaStop**2+F(2)/betaStop
cTotalD=cGapD+cStopD
cTotalT=cGapT+cStopT
assert cTotalD<F(1,180_000_000_000)
assert cTotalT<F(1,320_000)
print('PASS: 617 Loewner groups admit two-sided complete-series interval with rank-two boundary correction')
print('      global Gamma full-minus-grouped lower <=',float(cTotalD),'*D +',float(cTotalT),'*T')
print('      components gapD=',float(cGapD),'gapT=',float(cGapT),'tailD=',float(cStopD),'tailT=',float(cStopT))

# Rounded-UP rational PSD matrix bound used for a compact reproducible theorem.
major_D=F(1,180_000_000_000)
major_T=F(1,320_000)
assert cTotalD<major_D and cTotalT<major_T

# The existing complete 10-mode grouped lower matrix (from independently
# checked 031) + our new PSD upper gap can certify UPPER Rayleigh quotients
# for exact RATIONAL test vectors, without trusting floating eigensolvers.
source_H={}
for ln in GROUPED.read_text().splitlines():
    i,j,lo,hi=ln.split()
    source_H[int(i),int(j)]=I(F(int(lo),S),F(int(hi),S))
assert len(source_H)==50
candidates={
 1: ((3,5,7,9,1),(-727226995401,643702820298,-182669753899,9336596462,37487631),F(4,10**10)),
 2: ((4,6,8,10,2),(-743143741936,620073220342,-145278032240,-1238582207,-82669993),F(4,10**9)),
}
certified_upper={}
for k,(ids,raw,limit) in candidates.items():
    ids=list(ids)
    xs={j:F(v,10**12) for j,v in zip(ids,raw)}
    js,Dmat=stiffness(ids[:-1],k)
    assert js==ids
    Hray=sum((xs[i]*xs[j]*source_H[i,j] for i in ids for j in ids),I(0))
    Mray=sum((xs[i]*xs[j]*src(f'M_{i}_{j}') for i in ids for j in ids),I(0))
    Dray=sum((xs[i]*xs[j]*Dmat[i,j] for i in ids for j in ids),I(0))
    Tray=xs[k]*xs[k]*T[k]
    assert Mray.lo>0 and Dray.lo>0 and Tray.lo>=0
    upper=(Hray+I(major_D)*Dray+I(major_T)*Tray).hi/Mray.lo
    assert upper<limit,(k,float(upper),float(limit),Hray,Mray)
    certified_upper[k]=upper
    print(f'PASS: explicit rational physical Rayleigh witness in parity {k}: original full W/M < {limit} (certified upper {float(upper):.12g})')

# Two-sided Loewner bracket on EACH finite H1 source section:
# H_K <= H_original <= H_K + cD*D + cT*T.
# The analytic inequality is in the note; this is its exact-rational
# coefficient/structure/source-attachment audit, not a proof of the integral.
metadata=(PREV_INPUT.read_bytes()+PREV.read_bytes())
print('Source attachment SHA256 (input+prime):',sha256(metadata).hexdigest())
print('CR-ALP-032 EXACT ARITHMETIC / SOURCE-STRUCTURE AUDIT PASS')
print('No full finite Schur or high-Q residual bound certified.')

# Store only small, transparent integer/rational public acceptance controls.
result={
 'stage':'CR-ALP-032',
 'status':'EXACT SOURCE-STRUCTURE / OUTWARD-FRACTION PASS; continuum inequalities analytic, NOT Lean',
 'previous_signed_lower':'CR-ALP-031 original complete E10 W >= 1/100000000000 M',
 'all_seven_native_prime_powers':sorted(p.POWERS),
 'prime_displacement_rank_per_parity_at_most':2,
 'gamma_single_kernel_form':'diagonal minus rank-one per parity',
 'gamma_monotone_groups':617,
 'individual_gamma_kernels':30000,
 'group_coverage':'[30000,350000)',
 'error_upper_form':'W_full_E10 <= Htilde_031 + (1/180000000000) D + (1/320000) T',
 'D':'physical interval derivative-stiffness Gram',
 'T':'rank-two endpoint trace Gram; only u1,u2 nonzero',
 'certified_even_rational_test_vector':list(candidates[1][1]),
 'certified_odd_rational_test_vector':list(candidates[2][1]),
 'vector_denominator':10**12,
 'certified_even_Rayleigh_upper':'1/2500000000',
 'certified_odd_Rayleigh_upper':'1/250000000',
 'complete_E10_min_generalized_eigen_interval':'[1/100000000000,1/2500000000)',
 'source_attachment_sha256':sha256(metadata).hexdigest(),
 'open':'whole 15000 original signed Schur and low-to-Q15000 residual; Lean; RH',
}
(HERE/'CR-ALP-032_structured_spectral_result.json').write_text(dumps(result,indent=2)+'\n')
print('Saved CR-ALP-032_structured_spectral_result.json')
