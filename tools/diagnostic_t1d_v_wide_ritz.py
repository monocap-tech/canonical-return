#!/usr/bin/env python3
"""NONCERTIFYING CR-INT-T1D-V-WIDE Fourier-grid Gamma-only Ritz diagnostic.

Uses the original Gamma multiplier, authentic polynomial-based smooth tests
and the L2 metric. Deliberately EXCLUDES positive selected-Schur rank-one b b*/mu.
No spectral upper bound, rigorous interval, Suzuki-Connes transfer or RH.
"""
from math import pi
import numpy as np
import sympy as sp
from scipy.special import digamma
from scipy.linalg import eigh

x=sp.symbols('x',real=True)
a=sp.Rational(69,200)
profiles={}
for parity in (0,1):
    tests=[]
    for radius,k in ((a,3),(a,4),(a,5),(sp.Rational(1,4),5)):
        b=(x/radius)**parity*(1-(x/radius)**2)**k
        fn=sp.lambdify(x,sp.diff(b,x,2)-b/4,'numpy')
        tests.append((float(radius),fn))
    profiles[parity]=tests

def grid(period,N):
    dx=period/N
    xvalues=(np.arange(N)-N//2)*dx
    frequencies=2*pi*np.fft.fftfreq(N,dx)
    eps=1/1000
    rho=np.zeros(N)
    on=np.abs(xvalues)<eps
    rho[on]=np.exp(-1/(1-(xvalues[on]/eps)**2))
    rho_hat=np.fft.fft(np.fft.ifftshift(rho))
    rho_hat/=rho_hat[0]
    gamma=np.real(digamma(.25+.5j*frequencies))-np.log(pi)
    for parity,tests in profiles.items():
        source=np.zeros((N,4))
        for j,(radius,fn) in enumerate(tests):
            inside=np.abs(xvalues)<radius
            source[inside,j]=fn(xvalues[inside])
        fhat=np.fft.fft(source,axis=0)*rho_hat[:,None]
        mass=(fhat.conj().T@fhat).real/N
        G=(fhat.conj().T@(gamma[:,None]*fhat)).real/N
        assert np.linalg.eigvalsh(mass).min()>0
        eigs=eigh(6*mass,6*mass+G,eigvals_only=True)
        print(f'period={period}, grid={N}, '
              f'{"even" if parity==0 else "odd"}: '
              f'Gamma-only Ritz {np.round(eigs,8)}')

if __name__=='__main__':
    for period,N in ((12,2**17),(16,2**18),(24,2**18)):
        grid(period,N)
    print('NUMERICAL ONLY: no b-vector, no true full Ritz matrix, no interval or RH.')
