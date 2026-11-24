#!/usr/bin/env python3
"""Extend pinned CR-ALP-044 fitted residual Gram downward, not refit.
Exact intervals; inherited continuum coefficient and remainder formulas.
"""
from cr_alp_044_fitted_matrix_audit import *
from cr_alp_045_layered_gamma_audit import lower
CUTS=(85,90,100,110,111,112,114,116,118,123,150,172)
BASE=172
# The inherited digamma routine's shifted arctangent requires |y/x|<=1/3.
# For new lower frequencies use the exact recurrence to increase Re(z).
# This also improves its existing real-shift remainder bound.
def psi_safe(x,y):
 if y.l>=120*S:return psi(x,y)
 n=ceildiv(3*y.h,S)
 re,im=psi(x+n,y)
 for m in range(n):
  rr,ii=inv(x+m,y);re-=rr;im-=ii
 return re,im
saved=json.loads((ROOT/'CR-ALP-044_exact_result.json').read_text())
assert saved['stage']=='CR-ALP-044'
assert hashlib.sha256((ROOT/'CR-ALP-044_exact_result.json').read_bytes()).hexdigest()=='8d135f8ff80c9e82293b415217be6e09e151f1ca003497fef89cabfaef24dffc'
assert hashlib.sha256((ROOT/'CR-ALP-041_energy_conditioned_result.json').read_bytes()).hexdigest()=='6cc721da03a5dedadfb277566cc73845f08085b5ef794a6815bb9350e7c4b752'
assert hashlib.sha256((ROOT/'cr_alp_043_exact_band_audit.py').read_bytes()).hexdigest()=='a56ef90efd6fdd8d23b532def2bde92b300406e3598f9b827f5cab1145416de0'
# Use prior pinned matrices/energy basis, and stored actual fitted Gram.
# Independently recheck prior fit acceptance via CR-ALP-044 checker.
results=[]
for bl,old in zip(prior['blocks'],saved['blocks']):
 ids=bl['indices'];k=ids[-1];sgn=(-1)**k;s=-sgn
 assert old['parity']==k and old['N']==BASE and old['K']==30000
 fit=[I(F(x)) for x in old['fit']]
 G=[[I(int(t[0]),int(t[1]),True) for t in row] for row in old['fitted_band_gram_intervals']]
 tail=F(old['tail_operator_upper'])
 b0=(P['EHALF']+P['EHI'])/2 if k==1 else -(P['EHALF']-P['EHI'])/2
 b1=-(P['EHALF']-P['EHI'])/4 if k==1 else (P['EHALF']+P['EHI'])/4
 ep=(b0+2*b1)/2;em=(b0-2*b1)/2
 E0=ep*(2*a)+em*(1-P['EINV'].sq())-2*b0
 EC=[s*e*(b*b0-b1)/(bb-QUARTER) for b,bb,e in zip(BETA[1:],B2[1:],EX[1:])]
 Einf=E0+b1*S1+b0/4*S2-sum(EC).widen(1)
 t0=(log(NUF*BASE/2)-psi0)/2-Einf/b0
 coords=[];ygs=[]
 for v,bu in zip(bl['basis_rows'],bl['budgets']):
  d=I(F(bu['d'])).sqrt();assert d.l>0
  x=[I(F(z))/d for z in v];y=x[-1];co={i:z for i,z in zip(ids,x)}
  co[k]-=sum(co[i]*P[f'lam_{i}'] for i in ids[:-1]);coords.append(co);ygs.append(y*P[f'gamma_{k}'])
 for N in reversed(CUTS):
  upper=BASE if N==BASE else last
  for j in range(N+1,upper+1):
   if j%2!=k%2:continue
   v=j*NUF;v2=v.sq();den=v2+QUARTER
   psir,psii=psi_safe(QUARTER,v/2);rv=psii/(2*v);Lam=psir-psi0
   epart=sum(e/(bb+v2) for e,bb in zip(EX,B2)).widen(1);Rj=rv-sgn*epart
   trig=[sincos(v*t) for t in TS];ss=sum(w*t[0] for w,t in zip(W,trig));cc=sum(w*t[1] for w,t in zip(W,trig))
   CS={i:2/I(a)*(v*AI[i]/den-(v*AI[i]-NUS[i]*(v*Rj+ss))/(v2-NUS[i].sq())) for i in ids}
   cj=2*b0*v/(N0*den);R1=rv-1/den;S3=Lam/(2*v2)-2/den
   Es=(E0+b1*(S1-R1)+b0/4*(S2-S3))/den-sum(e/(bb+v2) for e,bb in zip(EC,B2[1:])).widen(1)
   Cb=Lam*(cj-b0/(N0*v))-2*v/N0*Es-t0*cj-2/(N0*den)*(v*b0*cc-b1*ss)
   u=[sum(co[i]*CS[i] for i in ids)-yg*Cb-ti*cj for co,yg,ti in zip(coords,ygs,fit)]
   for i in range(7):
    for l in range(i+1):G[i][l]=G[l][i]=G[i][l]+(u[i]*u[l] if i!=l else u[i].sq())
  last=N
  # search rational grid for a strictly verified upper bound, no floating eigensolver.
  for m in range(1,201):
   band=F(m,100)
   try:rad,piv=pd([[I(band if i==l else 0)-G[i][l] for l in range(7)] for i in range(7)])
   except AssertionError:continue
   break
  else:raise RuntimeError('no upper bound <= 2')
  mu,_=lower(N);rho=band+tail
  print('parity',k,'N',N,'band',band,'tail',tail,'rho',float(rho),'mu',float(mu),'PASS' if rho<mu else 'inconclusive',flush=True)
  results.append(dict(N=N,parity=k,band_upper=str(band),tail_upper=str(tail),rho_upper=str(rho),mu_lower=str(mu),pass_sufficient=rho<mu,gram_intervals=[[[str(t.l),str(t.h)] for t in row] for row in G]))
Path(__file__).with_name('CR-ALP-045_fitted_extension_result.json').write_text(json.dumps(results,indent=2)+'\n')
