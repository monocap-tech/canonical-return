#!/usr/bin/env python3
"""Numerical DIAGNOSTIC only: original complete signed Weil 6x6 interior block.

Uses exact mathematical physical kernel identity but numerical mpmath quadrature.
This is NOT an interval proof: neither entries nor generalized eigenvalues
are certified, and it cannot establish 15k or entire-collar positivity.
"""
import csv
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy.linalg as la

mp.mp.dps=45
A=mp.mpf(11)/10
P={2:2,3:3,4:2,5:5,7:7,8:2,9:3}
JS=(3,4,5,6,7,8)
nu=lambda j:j*mp.pi/(2*A)


def U(i,j,t):
    """Exactly paired original physical shift <e_i,T_t e_j>+<e_i,T_-t e_j>."""
    if (i+j)%2: return mp.mpf(0)
    if i==j:
        return (2*A-t)/A*mp.cos(nu(i)*t)+mp.sin(nu(i)*t)/(A*nu(i))
    return 2*(nu(i)*mp.sin(nu(j)*t)-nu(j)*mp.sin(nu(i)*t))/(A*(nu(i)**2-nu(j)**2))


def lam(j):
    k=1 if j%2 else 2
    return nu(j)*(nu(k)**2+mp.mpf('0.25'))/(nu(k)*(nu(j)**2+mp.mpf('0.25')))


def Uw(i,j,t):
    ai=1 if i%2 else 2
    aj=1 if j%2 else 2
    li,lj=lam(i),lam(j)
    return U(i,j,t)-li*U(ai,j,t)-lj*U(i,aj,t)+li*lj*U(ai,aj,t)


def mass(i,j):
    return mp.mpf(i==j)+(lam(i)*lam(j) if i%2==j%2 else 0)


def rho(s):
    return mp.exp(-s/2)/(-mp.expm1(-2*s))


def gamma_const():
    q=mp.exp(-A)
    return mp.digamma(mp.mpf('0.25'))-mp.log(mp.pi)+2*(mp.atanh(q)+mp.atan(q))


def H_entry(i,j):
    if (i+j)%2: return mp.mpf(0),mp.mpf(0)
    m=mass(i,j)
    def fn(s):
        if not s: return mp.mpf(0)
        return rho(s)*(2*m-Uw(i,j,s))
    g=gamma_const()*m+mp.quad(fn,[0,mp.mpf('0.001'),mp.mpf('0.1'),1,2*A])
    prime=-sum((mp.log(p)/mp.sqrt(n)*Uw(i,j,mp.log(n)) for n,p in P.items()),mp.mpf(0))
    return g,prime


def main():
    h=np.zeros((len(JS),len(JS)))
    mm=np.zeros_like(h)
    gd=np.zeros_like(h)
    pp=np.zeros_like(h)
    for a,i in enumerate(JS):
        for b,j in enumerate(JS[:a+1]):
            g,p=H_entry(i,j)
            gd[a,b]=gd[b,a]=float(g)
            pp[a,b]=pp[b,a]=float(p)
            h[a,b]=h[b,a]=float(g+p)
            mm[a,b]=mm[b,a]=float(mass(i,j))
    assert np.allclose(h,h.T) and np.allclose(mm,mm.T)
    eig=la.eigh(h,mm,check_finite=True)[0]
    print('Original full SIGNED 6x6 pole-null w_3,...,w_8: NUMERICAL ONLY')
    for j in range(6):
        print('w_%d: Gamma=% .13g prime=% .13g complete=% .13g'%(JS[j],gd[j,j],pp[j,j],h[j,j]))
    print('Generalized 6x6 eig(H,M), DIAGNOSTIC NOT CERTIFIED:')
    print(np.array2string(eig,precision=11))
    out=Path(__file__).with_name('CR-ALP-028_original_signed_6x6_diagnostic.csv')
    with out.open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['entry','basis_i','basis_j','physical_M','complete_signed_H','Gamma','original_prime','certification'])
        for a,i in enumerate(JS):
            for b,j in enumerate(JS[:a+1]):
                writer.writerow(['matrix',i,j,f'{mm[a,b]:.16e}',f'{h[a,b]:.16e}',f'{gd[a,b]:.16e}',f'{pp[a,b]:.16e}','NUMERICAL DIAGNOSTIC ONLY'])
        for v in eig: writer.writerow(['generalized_eigenvalue','','','','','','','NUMERICAL DIAGNOSTIC %.16e'%v])
    print('Wrote',out)

if __name__=='__main__':main()
