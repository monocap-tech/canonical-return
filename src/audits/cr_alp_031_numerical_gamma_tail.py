#!/usr/bin/env python3
"""NONCERTIFYING floating-point diagnostic for uncomputed positive Gamma tails.
Starting matrix is the rigorous CR-ALP-031 30,000-term interval midpoint;
only the incremental arithmetic in this file is floating-point.
"""
import numpy as np
import scipy.linalg as la
from pathlib import Path
from math import pi,sinh,exp,sqrt
here=Path(__file__).resolve().parent
A=1.1
NU={j:j*pi/(2*A) for j in range(1,11)}
E=exp(A); EINV=exp(-A)
SinhA=sinh(A);SqrtA=sqrt(A)
SCALE=10**25
ini={}
for line in (here/'CR-ALP-031_initial_rational_intervals_1e25.txt').read_text().splitlines():
 k,l,h=line.split();ini[k]=.5*(int(l)+int(h))/SCALE
H={}
for line in (here/'CR-ALP-031_H30000_exact_matrix_intervals.tsv').read_text().splitlines():
 i,j,l,h=line.split();H[(int(i),int(j))]=.5*(int(l)+int(h))/SCALE

def kb_sine(j,kind,beta,betasq):
 v=NU[j]; accum=np.zeros_like(beta)
 for p in (.5,-.5):
  sgn=(-1)**j
  ep=exp(p*A);em=exp(-p*A)
  mom=v*(em-sgn*ep)/(SqrtA*(v*v+.25))
  ld=em/(beta-p)
  rd=ep/(beta+p)
  val=(2*beta*mom+v/SqrtA*(ld-sgn*rd))/(betasq+v*v)
  accum+=(.5 if p==.5 else (.5 if kind==1 else -.5))*val
 return accum

def bb(kind,beta,betasq):
 vals={}
 for p in (.5,-.5):
  for r in (.5,-.5):
   z=p+r
   iz=2*A if z==0 else 2*SinhA
   ep=exp(z*A);em=1/ep
   vals[p,r]=2*beta*iz/(betasq-.25)-em/((beta+p)*(beta-r))-ep/((beta-p)*(beta+r))
 return (vals[.5,.5]+vals[-.5,-.5]+(1 if kind==1 else -1)*(vals[.5,-.5]+vals[-.5,.5]))/4

def sum_increments(start,stop,order):
 anchor=order[-1]; js=order[:-1]
 beta=2*np.arange(start,stop,dtype=float)+.5
 beta2=beta*beta
 h={j:NU[j]/(beta2+NU[j]**2) for j in (*js,anchor)}
 diag={j:2*NU[j]**2/(beta*(beta2+NU[j]**2)) for j in (*js,anchor)}
 lam={j:ini['lam_'+str(j)] for j in js}
 g=ini['gamma_'+str(anchor)]
 factor=2/A
 hk=h[anchor]
 kk=2*beta/(beta2+NU[anchor]**2)+factor*hk*hk
 kb={j:kb_sine(j,anchor,beta,beta2) for j in (*js,anchor)}
 bbb=bb(anchor,beta,beta2)
 d=np.zeros((len(order),len(order)))
 for ri,i in enumerate(js):
  vi=h[i]-lam[i]*hk
  for cj,j in enumerate(js):
   out=(diag[i] if i==j else 0) +lam[i]*lam[j]*diag[anchor]-factor*vi*(h[j]-lam[j]*hk)
   d[ri,cj]=np.sum(out)
  ki=factor*h[i]*hk
  kuw=ki-lam[i]*kk-g*(kb[i]-lam[i]*kb[anchor])
  out=2*ini[f'M_{i}_{anchor}']/beta-kuw
  d[ri,-1]=d[-1,ri]=np.sum(out)
 out=2*ini[f'M_{anchor}_{anchor}']/beta-(kk-2*g*kb[anchor]+g*g*bbb)
 d[-1,-1]=np.sum(out)
 return d

orders=[(3,5,7,9,1),(4,6,8,10,2)]
blocks=[]
for order in orders:
 mm=np.array([[ini[f'M_{i}_{j}'] for j in order] for i in order]);hh=np.array([[H[i,j] for j in order] for i in order]);blocks.append((mm,hh))

def show(tag):
 print('\n'+tag)
 for order,(mm,hh) in zip(orders,blocks):
  ev=la.eigh(hh,mm,eigvals_only=True)
  print(order,'min generalized eigenvalue',f'{ev[0]:.15g}', 'other low',f'{ev[1]:.15g}')
show('K=30000')
start=30000
for stop in [60000,100000,200000,500000,1000000,3000000]:
 for t,(order,(mm,hh)) in enumerate(zip(orders,blocks)):
  hh+=sum_increments(start,stop,order)
 show('K='+str(stop))
 start=stop
