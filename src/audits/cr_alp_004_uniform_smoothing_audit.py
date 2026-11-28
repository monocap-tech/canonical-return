#!/usr/bin/env python3
"""CR-ALP-004 numerical/symbolic regression. NOT a kernel certificate or RH proof.

Requires numpy, scipy, sympy. Numerically checks smoothing error sizes and
prime/pole/Gamma decomposition for an explicitly normalized C-infinity mollifier.
"""
import math
import numpy as np
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss
from sympy import primerange

NODES = 56
nodes, weights = leggauss(NODES)
profile = np.exp(-1.0/(1-nodes**2))
weights = weights*profile
weights /= weights.sum()
diff = (nodes[:,None]+nodes[None,:]).reshape(-1)
w_pair = (weights[:,None]*weights[None,:]).reshape(-1)

def ramp(t, X):
    t = np.asarray(t, dtype=float)
    return np.maximum(0.0, X*np.exp(-np.abs(t)/2)-np.exp(np.abs(t)/2))

def sm_ramp(t, X, width):
    return float(np.dot(w_pair,ramp(t-width*diff,X)))

def gamma_weight(t):
    return math.exp(-t/2)/(-math.expm1(-2*t))

def gamma_change(X, width):
    L=math.log(X)
    d0=sm_ramp(0,X,width)-(X-1)
    def f(t):
        if t==0:return 0.0
        dt=sm_ramp(t,X,width)-float(ramp(t,X))
        return (d0*math.exp(-t/2)-dt)*gamma_weight(t)
    pts=[0,2*width,1,L,L+2*width, max(30,L+30)]
    ints=[quad(f, a,b,epsabs=1e-6,epsrel=1e-5, limit=150)[0]
          for a,b in zip(pts[:-1],pts[1:]) if b>a]
    T=pts[-1]
    tail=d0*math.atanh(math.exp(-T))
    return -(0.5772156649015328606065+math.log(4*math.pi))*d0 + 2*(sum(ints)+tail)

def prime_change(X,width):
    end=int(math.ceil(X*math.exp(2*width)))
    out=0.0
    for p in primerange(2,end+1):
        k=int(p)
        while k<=end:
            t=math.log(k)
            out+=math.log(int(p))/math.sqrt(k)*(sm_ramp(t,X,width)-float(ramp(t,X)))
            k*=int(p)
    return -2*out

def pole_change(X,width):
    m=float(np.dot(weights,np.exp(width*nodes/2)))
    m2=float(np.dot(weights,np.exp(-width*nodes/2)))
    assert abs(m-m2)<1e-13
    return 2*(X-1)*math.log(X)*(m*m2-1)

def test():
    for X in [4.,8.,20.,64.]:
        width=X**-2
        p=prime_change(X,width)
        a=pole_change(X,width)
        g=gamma_change(X,width)
        error=abs(p+a+g)
        scale=math.log(math.e*X)/math.sqrt(X)
        # Conservative observational threshold, not proof of the O-bound.
        assert error < 300*scale,(X,error,scale)
        print(f'PASS X={X:g} width={width:.7g} prime={p:.8g} poles={a:.8g} gamma={g:.8g} '
              f'abs_error={error:.8g} scaled={error/scale:.6g}')
    from fractions import Fraction as Q
    assert (Q(1)+Q(1,2)-Q(9,2))==-3
    assert (-Q(1,16)-Q(10))==-Q(161,16)
    assert (Q(9,2)-2*Q(2))==Q(1,2)
    assert Q(10)-2*Q(8)==-6
    print('CR-ALP-004 NUMERICAL/EXPONENT REGRESSION PASS; not Lean, RH, or uniform theorem proof')

if __name__=='__main__':test()
