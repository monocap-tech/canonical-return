#!/usr/bin/env python3
"""T1D-V: certified four-dimensional (2 per parity) compact Schur spectral
LOWER bounds from genuine broad, smooth, double-Mellin-pole-null sources.

All matrix entries below are integrated exactly as SymPy rational numbers.
Gamma upper bound imports the DIGAMMA/JENSEN analytical lemma and rational
m_Gamma(0)<-5372/1000 audited in T1D-U; the local script checks the remaining
all-rational finite series inequalities and polynomial Gram definiteness.

This is NOT a full Weil Gamma/Schur matrix evaluation, an upper bound for the
leading eigenvalue, a low-high tail certificate, a Lean proof or RH.
"""
from fractions import Fraction as F
from sympy import symbols, diff, expand, integrate, Matrix, Rational, simplify

X = symbols('x', real=True)
a = Rational(1,4)
epsilon = F(1,10000)  # choose one fixed even, nonnegative, smooth unit mollifier of this support

def toF(q):
    q = simplify(q)
    assert q.is_Rational, q
    return F(int(q.p),int(q.q))

def log_upper(r, N=90):
    """True upper bound for log(r), rational r>=1."""
    r=F(r); assert r>=1
    t=(r-1)/(r+1); t2=t*t; power=t; total=F(0)
    for k in range(N):
        total += 2*power/F(2*k+1)
        power *= t2
    return total + 2*power/(F(2*N+1)*(1-t2))

def gamma_bound_for_variance(V,M=50):
    """U.9 with its rationally proved m0 < -5.372, M=50."""
    V=F(V); quarter=F(1,4)
    finite=sum((V/((F(n)+quarter)*(4*(F(n)+quarter)**2+V)) for n in range(M+1)),F(0))
    tail=F(1,2)*log_upper(1+V/(4*(F(M)+quarter)**2))
    return -F(5372,1000)+finite+tail

# All prime-power self and cross-pairings vanish within each central
# probe subspace because support diameter <=1/2+2eps < log 2.
# Source positivity of Schur correction R>=0 is imported from T1D-E.
assert F(1,2)+2*epsilon < F(2,3)

results=[]
for parity, homogeneous_r, allowed_r, second_gamma_max, mixture_t, mixture_raw_r, mixed_r, top_gamma_max in [
    (0,37,38,F(181,100),F(13,12),F(18),F(181,10),F(53,50)),
    (1,45,46,F(2),F(9,5),F(578,25),F(116,5),F(131,100)),
]:
    ys=[]
    for k in (3,4):
        b = X**parity*(1-16*X*X)**k
        ys.append(expand(diff(b,X,2)-b/4))
    M=Matrix(2,2,lambda i,j:integrate(ys[i]*ys[j],(X,-a,a)))
    K=Matrix(2,2,lambda i,j:integrate(diff(ys[i],X)*diff(ys[j],X),(X,-a,a)))
    assert all(e.is_Rational for e in M) and all(e.is_Rational for e in K)
    assert M[0,0]>0 and M.det()>0
    pos=homogeneous_r**2*M-K
    assert pos[0,0]>0 and pos.det()>0, (parity,pos)
    # All linear combinations of the raw y0,k have bounded H1 ratio:
    # K <= (homogeneous_r)^2 M. After smooth mollification, Young and
    # Minkowski imply the ratio <= r/(1-eps*r).
    r_smoothed=F(homogeneous_r)/(1-epsilon*homogeneous_r)
    assert r_smoothed<allowed_r
    bulk_upper=gamma_bound_for_variance(F(allowed_r)**2)
    assert bulk_upper<second_gamma_max,(parity,bulk_upper,second_gamma_max)
    second_witness=6/(6+second_gamma_max)

    opt=Matrix([Rational(mixture_t.numerator,mixture_t.denominator)*0+1,
                Rational(mixture_t.numerator,mixture_t.denominator)])
    m_opt=toF((opt.T*M*opt)[0]);k_opt=toF((opt.T*K*opt)[0]);
    assert m_opt>0 and k_opt>0
    assert k_opt < mixture_raw_r**2*m_opt,(parity, k_opt/m_opt)
    opt_smoothed=mixture_raw_r/(1-epsilon*mixture_raw_r)
    assert opt_smoothed<mixed_r
    opt_gamma=gamma_bound_for_variance(mixed_r**2)
    assert opt_gamma<top_gamma_max,(parity,opt_gamma,top_gamma_max)
    top_witness=6/(6+top_gamma_max)
    assert top_witness>second_witness
    # Convert exact 2x2 rational Gram entries to reproducible textual certificate.
    results.append((parity,M,K,pos,mixture_t,m_opt,k_opt,r_smoothed,bulk_upper,
                    second_witness,opt_smoothed,opt_gamma,top_witness))
    label='even' if parity==0 else 'odd'
    print(f'PASS {label}: two independent compact H1 raw polynomial probes, physical Gram SPD')
    print(f'  exact M = {M.tolist()}')
    print(f'  exact K = {K.tolist()}')
    print(f'  exact upper matrix {homogeneous_r}^2 M-K SPD, det={pos.det()}')
    print(f'  gamma[f]/||f||^2 < {second_gamma_max} for ALL 2D smooth combinations; rational Jensen actual bound ~{float(bulk_upper):.12f}')
    print(f'  SECOND positive Birman-Schwinger eigenvalue > {second_witness} ({float(second_witness):.9f})')
    print(f'  improved single witness y3 + ({mixture_t}) y4, raw ratio <{mixture_raw_r}, smoothed ratio <{mixed_r}')
    print(f'  gamma[optimal]/||optimal||^2 < {top_gamma_max}; exact Jensen bound ~{float(opt_gamma):.12f}')
    print(f'  FIRST positive Birman-Schwinger eigenvalue > {top_witness} ({float(top_witness):.9f})')

assert len(results)==2
print('CERTIFIED: at least TWO eigenvalues >600/781 in even sector, TWO >3/4 in odd sector')
print('CERTIFIED: top even >300/353, top odd >600/731 (actual original Schur operator)')
print('CHECKER SCOPE: exact rational polynomial Gram; finite Jensen inequalities; mGamma(0) bound imported from U')
print('OPEN: exact mollified Gamma matrix, true Schur cross matrix, eta_N, top spectral upper bound, RH')
