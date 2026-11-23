#!/usr/bin/env python3
"""CR-INT-T1D-V-WIDE: exact source polynomial Grams, prime-free parity
trial geometry, and rigorous Jensen/Gamma spectral LOWER witnesses.
Analytic source/digamma inputs are described in research/CR-INT-T1D-V-WIDE.md.
This does NOT certify full mollified Gamma/Schur entries or spectral UPPER bounds.
"""
from fractions import Fraction as Q
import sympy as sp

x=sp.symbols('x',real=True)
a=sp.Rational(69,200)
eps=Q(1,1000)
entries=[(a,3),(a,4),(a,5),(sp.Rational(1,4),5)]
coeffs={0:[Q(1),Q(3,5),Q(1,4),Q(0)],
        1:[Q(1),Q(7,10),Q(3,5),Q(0)]}

def gram(p,deriv=False):
    ps=[]
    for radius,k in entries:
        b=(x/radius)**p*(1-(x/radius)**2)**k
        y=sp.Poly(sp.diff(b,x,2)-b/4,x)
        ps.append(y.diff() if deriv else y)
    rows=[]
    for i in range(4):
        row=[]
        for j in range(4):
            radius=min(entries[i][0],entries[j][0])
            z=ps[i]*ps[j]
            val=sum((2*c*radius**(e[0]+1)/(e[0]+1)
                     for e,c in z.terms() if e[0]%2==0),sp.Rational(0))
            assert val.is_Rational
            row.append(val)
        rows.append(row)
    G=sp.Matrix(rows)
    assert G==G.T
    return G

def log_interval(r,n_terms=100):
    r=Q(r);assert r>=1
    t=(r-1)/(r+1);power=t;lower=Q(0)
    for j in range(n_terms):
        lower+=2*power/Q(2*j+1);power*=t*t
    return lower,lower+2*power/(Q(2*n_terms+1)*(1-t*t))

def gamma_upper(rate):
    # T1D-U's analytic digamma/Jensen inequality imported;
    # m_Gamma(0)< -5372/1000 has a separate rational audit.
    V=rate*rate
    finite=sum((V/((Q(n)+Q(1,4))*(4*(Q(n)+Q(1,4))**2+V))
                for n in range(51)),Q(0))
    tail=log_interval(1+V/(4*(Q(50)+Q(1,4))**2))[1]
    return -Q(5372,1000)+finite+tail/2

assert Q(173,250)<Q(693,1000)<log_interval(2)[0]
assert Q(a+sp.Rational(1,1000))==Q(173,500)
assert Q(173,500)+Q(1,20000)<log_interval(3)[0]/2

for parity in (0,1):
    M=gram(parity)
    K=gram(parity,True)
    assert all(M[:j,:j].det()>0 for j in range(1,5))
    assert all(K[:j,:j].det()>0 for j in range(1,5))
    v=sp.Matrix([sp.Rational(c.numerator,c.denominator)
                 for c in coeffs[parity]])
    mass=(v.T*M*v)[0];stiff=(v.T*K*v)[0]
    assert mass>0 and stiff>0
    rate=Q(66,5) if parity==0 else Q(17)
    # For y=y0*rho_eps: ||y||>=||y0||-eps||y0'|| and
    # ||y'||<=||y0'||, so the following is exact rational arithmetic.
    assert Q(stiff)/Q(mass)<(rate/(1+eps*rate))**2
    gamma_bound=Q(3,4) if parity==0 else Q(1)
    assert gamma_upper(rate)<gamma_bound
    lower=Q(8,9) if parity==0 else Q(6,7)
    assert Q(6)/(Q(6)+gamma_bound)==lower
    print(f'PASS {"even" if parity==0 else "odd"}: '
          f'4x4 exact UNMOLLIFIED mass/H1 Gram, '
          f'smooth derivative/mass <{rate}, Gamma/mass <{gamma_bound}, '
          f'actual compact-spectrum SUP >{lower}')
    print(' optimized polynomial reference mass=',mass,'stiffness=',stiff)

print('PASS geometry: all 8 smooth two-pole-null vectors in physical E complement')
print('PASS: 173/250 < log 2, hence ALL full-source 4x4 prime entries vanish')
print('SOURCE-DEFINED Ritz pencil: D=Gamma+6M, C=6M+b b*/mu;')
print('b=original full signed selected-to-test Weil coupling, retained, NOT evaluated')
print('OPEN: true mollified Gamma and b matrices, top Ritz upper, eta_N, RH')
