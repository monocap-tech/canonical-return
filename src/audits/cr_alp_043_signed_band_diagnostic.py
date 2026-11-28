#!/usr/bin/env python3
"""Floating diagnostic only: actual original signed coefficients on finite high-sine bands.
No interval certification, no negative Weil witness. Tail projection is on the finite band,
so exact band variance would lower-bound, not upper-bound, the full residual Gram.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
from fractions import Fraction as F
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
S=10**25; a=1.1; L=2*a; ns=np.array([2,3,4,5,7,8,9]); ts=np.log(ns)
w=np.log([2,3,2,5,7,2,3])/np.sqrt(ns)
nu=lambda j:np.asarray(j)*np.pi/L
H={}; params={}
for ln in (ROOT/'CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv').read_text().splitlines():
 i,j,l,h=map(int,ln.split()); H[i,j]=(l+h)/(2*S)
for ln in (ROOT/'CR-ALP-040_initial_rational_intervals_1e25.txt').read_text().splitlines():
 k,l,h=ln.split(); params[k]=(int(l)+int(h))/(2*S)
prior=json.loads((ROOT/'CR-ALP-041_energy_conditioned_result.json').read_text())
# R(v)=sum_m (beta_m^2+v^2)^(-1). Exact digamma identity.
def R(v):return digamma(.25+.5j*v).imag/(2*v)
def Lambda(v):return digamma(.25+.5j*v).real-digamma(.25)
# 50 terms suffice for exponential correction numerically; discarded part not certified here.
beta=2*np.arange(50)+.5; ex=np.exp(-L*beta)
def sine_coeff(i,j):
 v=nu(j); u=nu(i); sign=(-1)**i
 base=(R(u)-R(v))/(v*v-u*u)
 expcorr=np.sum((ex/(beta*beta+u*u))[:,None]/(beta[:,None]**2+v[None,:]**2),axis=0)
 g=-2*v*u/a*(base-sign*expcorr)
 p=-2/a*np.sum(w[:,None]*(v[None,:]*np.sin(u*ts[:,None])-u*np.sin(ts[:,None]*v[None,:]))/(v*v-u*u),axis=0)
 return g+p

def pole_coeff(k,j):
 v=nu(j); s=(-1)**(k+1); b0=np.cosh(a/2) if k==1 else -np.sinh(a/2); b1=-.5*np.sinh(a/2) if k==1 else .5*np.cosh(a/2)
 c=2*b0*v/(np.sqrt(a)*(v*v+.25))
 # E_m=I_b(beta)-b0/beta. m=0 integrated exactly; higher m split into nonexponential partial fractions and exponential correction.
 ep=(b0+2*b1)/2; em=(b0-2*b1)/2
 E0=ep*L+em*(-np.expm1(-L))-2*b0
 R1=R(v)-1/(v*v+.25)
 S1=1-np.log(2)
 S2=2*digamma(1.25)-digamma(1)-digamma(1.5)
 S3=Lambda(v)/(2*v*v)-2/(v*v+.25)
 b=beta[1:]; ee=ex[1:]
 expcorr=np.sum((ee*s*(b*b0-b1)/(b*b-.25))[:,None]/(b[:,None]**2+v[None,:]**2),axis=0)
 es=E0/(v*v+.25)+(b1*(S1-R1)+b0/4*(S2-S3))/(v*v+.25)-expcorr
 g=Lambda(v)*(c-b0/(np.sqrt(a)*v))-2*v/np.sqrt(a)*es
 bt=b0*np.cosh(ts/2)+2*b1*np.sinh(ts/2)
 # Exact prime translation coefficient: elementary primitive evaluated at both overlap endpoints.
 p=-2/np.sqrt(a)*np.sum(w[:,None]*(v[None,:]*(bt[:,None]+b0*np.cos(ts[:,None]*v[None,:]))-b1*np.sin(ts[:,None]*v[None,:]))/(v*v+.25),axis=0)
 return g+p,c

def run():
 out=[]
 for block in prior['blocks']:
  ids=block['indices']; k=ids[-1]; j=np.arange(173 if k==1 else 174,100001,2)
  a0=sine_coeff(k,j); b,c=pole_coeff(k,j)
  # g(0) and the other pole-aligned terms are killed by the band projector.
  T=np.column_stack([sine_coeff(i,j)-params[f'lam_{i}']*a0 for i in ids[:-1]]+[a0-params[f'gamma_{k}']*b])
  V=np.array([[float(F(t)) for t in v] for v in block['basis_rows']]); d=np.array([float(F(b['d'])) for b in block['budgets']]);
  U=(T@V.T)/np.sqrt(d)
  for N,mu in [(172,.25),(242,.6),(300,.75),(500,1.2),(1000,1.8),(2000,2.5),(5000,3.5),(15000,4.7)]:
   use=j>N; C=c[use]; A=U[use]; A=A-np.outer(C,(C@A)/(C@C)); B=A.T@A
   ev=np.linalg.eigvalsh(B)
   rec=dict(parity=k,N=N,band_end=int(j[-1]),trace=float(np.trace(B)),spectral=float(ev[-1]),illustrative_mu=mu,ratios=np.diag(B).tolist())
   out.append(rec); print(k,N,'trace',f'{np.trace(B):.7g}','max',f'{ev[-1]:.7g}','threshold',mu)
 (ROOT/'CR-ALP-043_signed_band_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':run()
