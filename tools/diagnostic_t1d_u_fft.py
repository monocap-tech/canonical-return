#!/usr/bin/env python3
"""NONCERTIFYING CR-INT-T1D-U numeric Fourier diagnostic.

Requires numpy, scipy, sympy. Uses an explicitly fixed C^infinity even
probability mollifier and periodic Fourier sums, with fixed support safely
inside the period. Reports approximate Gamma Rayleigh values and true
24-prime-power compressed arithmetic-image L2 norms for the broad T probes.
This provides neither interval bounds nor any upper compact eigenvalue.
Run separately; it is NOT part of the audited regression runner.
"""
from math import isqrt, log, pi, sqrt
import numpy as np
from scipy.special import digamma
import sympy as sp


def prime_power_base(n):
    for p in range(2, n + 1):
        if n % p:
            continue
        if any(p % q == 0 for q in range(2, isqrt(p) + 1)):
            continue
        k = n
        while k % p == 0:
            k //= p
        if k == 1:
            return p
    return None


PP = {n: prime_power_base(n) for n in range(2, 55) if prime_power_base(n)}
assert len(PP) == 24
X = sp.symbols('X', real=True)
BASE = {'even': (1 - 16 * X**2)**5,
        'odd': X * (1 - 16 * X**2)**5}
FUN = {k: sp.lambdify(X, sp.diff(v, X, 2) - v / 4, 'numpy')
       for k, v in BASE.items()}


def run(period, grid):
    dx = period / grid
    x = (np.arange(grid) - grid // 2) * dx
    xi = 2 * pi * np.fft.fftfreq(grid, d=dx)
    # Fixed nonnegative even normalized C^infinity probability mollifier.
    eps = 1 / 200
    rho = np.zeros(grid)
    mask = np.abs(x) < eps
    rho[mask] = np.exp(-1 / (1 - (x[mask] / eps)**2))
    # Shift the physical origin to array index 0 before convolution.
    rhohat = np.fft.fft(np.fft.ifftshift(rho))
    rhohat /= rhohat[0]
    gamma_symbol = np.real(digamma(.25 + .5j * xi)) - np.log(pi)
    # Complete original 24 prime powers, not just central contributors.
    prime_symbol = sum(2 * log(p) / sqrt(n) * np.cos(xi * log(n))
                       for n, p in PP.items())
    inner = np.abs(x) < .25
    physical = np.abs(x) < 2
    for label, f in FUN.items():
        y0 = np.zeros(grid)
        y0[inner] = f(x[inner])
        yhat = np.fft.fft(y0) * rhohat
        y = np.fft.ifft(yhat).real
        A_y = np.fft.ifft(yhat * prime_symbol).real
        A_y[~physical] = 0
        w = np.abs(yhat)**2
        gam = float(np.dot(gamma_symbol, w) / np.sum(w))
        relative_arithmetic = float(sqrt(np.sum(A_y**2) / np.sum(y**2)))
        prime_self = float(np.dot(y, A_y) / np.dot(y, y))
        print(f'L={period:2} N={grid:7} {label:4} gamma={gam:.12f} '
              f'||Ay||/||y||={relative_arithmetic:.9f} '
              f'<Ay,y>/||y||^2={prime_self:.3e} '
              f'6/(gamma+6)={6/(6+gam):.9f}')


if __name__ == '__main__':
    for period, grid in ((16, 2**18), (32, 2**19), (48, 2**19)):
        run(period, grid)
    print('NUMERICAL ONLY: no interval, no top eigenvalue, no eta_N, no RH.')
