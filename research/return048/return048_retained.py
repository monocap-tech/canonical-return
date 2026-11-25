#!/usr/bin/env python3
"""Original signed retained energy; inherited source/Gamma identities."""
import hashlib,json
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache as cache
import mpmath as mp
from mpmath import iv
import numpy as np
root=Path(__file__).resolve().parent.parent
raw=(root/'return047/RETURN047_payment_receipt.json').read_bytes();assert hashlib.sha256(raw).hexdigest()=='def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f'
r=json.loads(raw);ks=sorted(set(r['even_exceptional']+r['unpaid_odd']));assert len(ks)==85 and ks==sorted(-k for k in ks)
iv.dps=65;mp.mp.dps=75;S=10**55;Q=10**40
zero=(0,0);a=iv.mpf(5)/4;D=2*a;B=iv.mpf(5)/2
add=lambda x,y:(x[0]+y[0],x[1]+y[1]);neg=lambda x:(-x[0],-x[1]);sub=lambda x,y:add(x,neg(y))
@cache(None)
def z(x):return iv.mpc(iv.mpf(x[0])/2,2*iv.pi*x[1]/5)
@cache(None)
def exa(x):return iv.exp(iv.mpf(x[0])*a/2)*[iv.mpc(1),iv.mpc(0,1),iv.mpc(-1),iv.mpc(0,-1)][x[1]%4]
@cache(None)
def I(x):return iv.mpc(D) if x==zero else (exa(x)-exa(neg(x)))/z(x)
@cache(None)
def Ip(x):return iv.mpc(0) if x==zero else a*(exa(x)+exa(neg(x)))/z(x)-(exa(x)-exa(neg(x)))/z(x)**2
P=24;N=64
bern=[iv.mpf(str(mp.bernfrac(2*j)[0]))/int(mp.bernfrac(2*j)[1]) for j in range(1,P+1)]
errpsi=abs(bern[-1])/(2*P*iv.mpf(N)**(2*P));errtri=abs(bern[-1])/iv.mpf(N)**(2*P+1)
def pad(v,e):
 t=iv.mpf([-e.b,e.b]);return iv.mpc(v.real+t,v.imag+t)
@cache(None)
def psi(x):
 v=(B+z(x))/2;w=v+N
 y=iv.ln(w)-1/(2*w)-sum(b/(2*j*w**(2*j)) for j,b in enumerate(bern,1))-sum(1/(v+j) for j in range(N))
 return pad(y,errpsi)
@cache(None)
def tri(x):
 v=(B+z(x))/2;w=v+N
 y=1/w+1/(2*w*w)+sum(b/w**(2*j+1) for j,b in enumerate(bern,1))+sum(1/(v+j)**2 for j in range(N))
 return pad(y,errtri)
@cache(None)
def ss(x,y):return tri(x)/4 if x==y else (psi(y)-psi(x))/(2*z(sub(y,x)))
betas=[iv.mpf(2*m)+iv.mpf(1)/2 for m in range(1,25)];BB=iv.mpf(101)/2
errdec=iv.exp(-BB*D)/((BB-iv.mpf(1)/2)**2*(1-iv.exp(-2*D)))
# use divided difference of exponentially weighted simple sums, with direct repeated-denominator case
@cache(None)
def ds(x):return sum(iv.exp(-b*D)/(b+z(x)) for b in betas)
@cache(None)
def dec(x,y):
 v=sum(iv.exp(-b*D)/(b+z(x))**2 for b in betas) if x==y else (ds(x)-ds(y))/z(sub(y,x))
 return pad(v,errdec)
@cache(None)
def triangle(x,y):
 b=(1,0);den=add(y,b)
 if den==zero:return Ip(sub(x,b))+a*I(sub(x,b))
 return (I(add(x,y))-exa(neg(den))*I(sub(x,b)))/z(den)
G0=-iv.euler-iv.pi/2-3*iv.ln(2)-iv.ln(iv.pi)
@cache(None)
def gpart(x):return (psi(x)+psi(neg(x)))/2-psi(zero)
primes=[(2,2),(4,2),(8,2),(3,3),(9,3),(5,5),(7,7),(11,11)]
assert iv.ln(11).b<D.a and iv.ln(13).a>D.b
shifts=[(iv.ln(n),iv.ln(p)/iv.sqrt(n)) for n,p in primes]
@cache(None)
def shifted_I(x,j):
 s,w=shifts[j];b=a-s/2
 if x==zero:return iv.mpc(2*b)
 v=z(x)*b;return (iv.exp(v)-iv.exp(-v))/z(x)
@cache(None)
def pair_factor(x,j):
 s,w=shifts[j];v=z(x)*s/2;return w*(iv.exp(v)+iv.exp(-v))
@cache(None)
def raw_entry(x,y):
 # B is symmetric bilinear here, not sesquilinear, on exponentials.
 u=add(x,y);v=sub(y,x);mass=I(u)
 gamma=(G0+gpart(x))*mass+4*mass-triangle(x,y)-triangle(y,x)
 gamma+=exa(neg(u))*ss(x,neg(y))+exa(u)*ss(neg(x),y)
 gamma-=exa(v)*dec(x,neg(y))+exa(neg(v))*dec(neg(x),y)
 prime=sum(pair_factor(v,j)*shifted_I(u,j) for j in range(len(shifts)))
 b=(1,0);pole=I(add(x,b))*I(sub(y,b))+I(sub(x,b))*I(add(y,b))
 return gamma-prime+pole

def endpoints(v,scale=S):
 def rational(t):
  s,m,e,b=t;return F((-1 if s else 1)*m)*F(2)**e
 l,u=map(rational,v._mpi_);return [l.numerator*scale//l.denominator,-((-u.numerator*scale)//u.denominator)]
def matrix(cols,fun):
 out=[]
 for i,c in enumerate(cols):
  row=[]
  for j,d in enumerate(cols):
   if j<i:row.append(out[j][i]);continue
   x=sum(v*w*fun(a,b) for a,v in c for b,w in d)
   assert x.imag.a<=0<=x.imag.b
   row.append(endpoints(x.real))
  out.append(row)
 return out

def cert(A):
 n=len(A);lo=np.array([[v[0] for v in row] for row in A],dtype=object);hi=np.array([[v[1] for v in row] for row in A],dtype=object)
 mid=mp.matrix([[mp.mpf(int(lo[i,j]+hi[i,j]))/(2*S) for j in range(n)] for i in range(n)])
 try:L=mp.cholesky(mid);K=(L.T)**-1
 except ValueError:
  E,U=mp.eigsy(mid);print('midpoint min',mp.nstr(E[0],25),flush=True)
  x=np.array([int(mp.nint(U[i,0]*Q)) for i in range(n)],dtype=object);up=int(x@(lo+hi)@x+abs(x)@(hi-lo)@abs(x));return {'status':'NEGATIVE' if up<0 else 'NOT CERTIFIED','witness':x.tolist(),'upper_numerator':str(up)}
 C=np.array([[int(mp.nint(K[i,j]*Q)) for j in range(n)] for i in range(n)],dtype=object)
 A0=C.T@(lo+hi)@C;R=abs(C).T@(hi-lo)@abs(C)
 margin=min(int(A0[i,i]-R[i,i]-sum(abs(A0[i,j])+R[i,j] for j in range(n) if j!=i)) for i in range(n))
 return {'status':'PASS' if margin>0 else 'NOT CERTIFIED','congruence':C.tolist(),'margin_numerator':str(margin),'margin_denominator':str(2*S*Q*Q),'frobenius_numerator':str(sum(int(x)**2 for x in C.flat))}
def run():
 out={'scope':'actual original retained signed form; no mixed complement gate','parent_sha256':hashlib.sha256(raw).hexdigest(),'scale':str(S),'congruence_scale':str(Q),'psi_error':endpoints(errpsi),'trigamma_error':endpoints(errtri),'decay_error':endpoints(errdec),'blocks':[]}
 for parity in (0,1):
  cols=[];names=[]
  if parity==0:cols.append([(zero,iv.mpc(1))]);names.append('1')
  for k in ks:
   if k<=0:continue
   cols.append([((0,k),iv.mpc(iv.mpf(1)/2) if parity==0 else iv.mpc(0,-iv.mpf(1)/2)),((0,-k),iv.mpc(iv.mpf(1)/2) if parity==0 else iv.mpc(0,iv.mpf(1)/2))]);names.append(('cos' if parity==0 else 'sin')+str(k))
  if parity==1:cols.append([((1,0),iv.mpc(iv.mpf(1)/2)),((-1,0),iv.mpc(-iv.mpf(1)/2))]);names.append('sinh(x/2)')
  mass=matrix(cols,lambda a,b:I(add(a,b)));print('mass',parity,flush=True)
  H=matrix(cols,raw_entry);print('energy built',parity,'width',max(v[1]-v[0] for row in H for v in row),flush=True)
  mass_upper=max(sum(max(abs(v[0]),abs(v[1])) for v in row) for row in mass)
  c=cert(H);print('parity',parity,len(cols),c['status'],flush=True)
  if c['status']=='PASS':
   lower=F(int(c['margin_numerator']),2*mass_upper*int(c['frobenius_numerator']));c['physical_floor']=str(lower);print('physical floor',float(lower),flush=True)
  out['blocks'].append({'parity':parity,'names':names,'mass':mass,'energy':H,'mass_upper_numerator':str(mass_upper),**c})
 p=Path(__file__).with_name('RETURN048_retained_receipt.json');p.write_text(json.dumps(out,indent=2)+'\n');print('receipt SHA256',hashlib.sha256(p.read_bytes()).hexdigest(),flush=True)
if __name__=='__main__':run()
