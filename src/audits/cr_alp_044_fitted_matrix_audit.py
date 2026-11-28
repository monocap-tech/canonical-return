#!/usr/bin/env python3
"""Exact fitted residual Gram and parity-resolved infinite remainder.
Uses CR-ALP-043 interval primitives and analytical coefficient identities.
All acceptance arithmetic integer/Fraction. No new Lean certification.
"""
from cr_alp_043_exact_band_audit import *
N=172
K=30000
J=[F(factorial(m),2*K)+F(2,K*K) for m in range(3)]
def rup(x,den=10**6):return F(ceildiv(x.numerator*den,x.denominator),den)
def pd(A):
 n=len(A);C=[[F(x.l+x.h,2*S) for x in row] for row in A]
 err=max(sum(F(x.h-x.l,2*S) for x in row) for row in A)
 for i in range(n):C[i][i]-=err
 piv=[]
 for k in range(n):
  p=C[k][k];assert p>0,(k,p);piv.append(p)
  for i in range(k+1,n):
   for j in range(k+1,n):C[i][j]-=C[i][k]*C[k][j]/p
 return err,min(piv)
def run():
 results=[]
 for bl in prior['blocks']:
  ids=bl['indices'];k=ids[-1];sgn=(-1)**k;s=-sgn
  b0=(P['EHALF']+P['EHI'])/2 if k==1 else -(P['EHALF']-P['EHI'])/2
  b1=-(P['EHALF']-P['EHI'])/4 if k==1 else (P['EHALF']+P['EHI'])/4
  ep=(b0+2*b1)/2;em=(b0-2*b1)/2
  E0=ep*(2*a)+em*(1-P['EINV'].sq())-2*b0
  EC=[s*e*(b*b0-b1)/(bb-QUARTER) for b,bb,e in zip(BETA[1:],B2[1:],EX[1:])]
  Einf=E0+b1*S1+b0/4*S2-sum(EC).widen(1)
  t0=(log(NUF*N/2)-psi0)/2-Einf/b0
  coords=[];ygs=[]
  for v,bu in zip(bl['basis_rows'],bl['budgets']):
   d=I(F(bu['d'])).sqrt();assert d.l>0
   x=[I(F(z))/d for z in v];y=x[-1];co={i:z for i,z in zip(ids,x)}
   co[k]-=sum(co[i]*P[f'lam_{i}'] for i in ids[:-1]);coords.append(co);ygs.append(y*P[f'gamma_{k}'])
  G=[[ZERO for _ in ids] for _ in ids];cross=[ZERO for _ in ids];cc2=ZERO
  trigstep=[sincos(2*NUF*t) for t in TS];trig=None;nrows=0
  for j in range(N+1+(k==2),K+1,2):
   v=j*NUF;v2=v.sq();den=v2+QUARTER
   psir,psii=psi(QUARTER,v/2);rv=psii/(2*v);Lam=psir-psi0
   epart=sum(e/(bb+v2) for e,bb in zip(EX,B2)).widen(1);Rj=rv-sgn*epart
   if nrows%16==0:trig=[sincos(v*t) for t in TS]
   else:trig=[(u*c+v0*s,v0*c-u*s) for (u,v0),(s,c) in zip(trig,trigstep)]
   ss=sum(w*t[0] for w,t in zip(W,trig));cc=sum(w*t[1] for w,t in zip(W,trig))
   CS={i:2/I(a)*(v*AI[i]/den-(v*AI[i]-NUS[i]*(v*Rj+ss))/(v2-NUS[i].sq())) for i in ids}
   cj=2*b0*v/(N0*den);R1=rv-1/den;S3=Lam/(2*v2)-2/den
   Es=(E0+b1*(S1-R1)+b0/4*(S2-S3))/den-sum(e/(bb+v2) for e,bb in zip(EC,B2[1:])).widen(1)
   Cb=Lam*(cj-b0/(N0*v))-2*v/N0*Es-t0*cj-2/(N0*den)*(v*b0*cc-b1*ss)
   u=[sum(co[i]*CS[i] for i in ids)-yg*Cb for co,yg in zip(coords,ygs)]
   cc2+=cj.sq()
   for i in range(7):
    cross[i]+=cj*u[i]
    for l in range(i+1):G[i][l]+=u[i]*u[l] if i!=l else u[i].sq()
   nrows+=1
   if nrows%2000==0:print('parity',k,'rows',nrows,flush=True)
  assert nrows==(K-N)//2 and cc2.l>0
  # Any rational fit is legal on the whole infinite tail.
  fit=[F(round((c/cc2).mid()*10**12),10**12) for c in cross]
  FG=[[ZERO for _ in ids] for _ in ids]
  for i in range(7):
   for l in range(i+1):
    FG[i][l]=FG[l][i]=G[i][l]-I(fit[i])*cross[l]-I(fit[l])*cross[i]+I(fit[i]*fit[l])*cc2
  band=F(18,100) if k==1 else F(16,100)
  rad,piv=pd([[I(band if i==l else 0)-FG[i][l] for l in range(7)] for i in range(7)])
  print('parity',k,'band bound',band,'min pivot',float(piv),'radius',float(rad),flush=True)
  C=(2*b0/(N0*NUF)).up();logKN=log(I(F(K,N)))
  B0,B1=(F(6,5),F(3,10)) if k==1 else (F(3,5),F(3,5))
  C2=F(121,225)*(F(12001,4000)*B0+F(11,7)*B1)+2*B1*sa*(2*a/3)**2
  C3=8*a**3*B0/945
  logtrace=F(0);primetrace=F(0);errtrace=F(0)
  for co,yg,f in zip(coords,ygs,fit):
   A=-yg/2; B=A*logKN-I(f)
   logtrace+=C*C*(B.up()**2*J[0]+2*B.up()*A.up()*J[1]+A.up()**2*J[2])
   primetrace+=(C*yg.up()*sa)**2*J[0]
   D=sum(co[i]*NUS[i] for i in ids).up()
   Cs={t:sum(co[i].up()*NUS[i].up()**t for i in ids) for t in (0,2,3)}
   AA=(2*sa/a+F(201,200)*F(22,7)/(2*a))*(2*a/3)**2*D
   BB=(sa*Cs[0]/(2*a)+8*sa*Cs[2]/a+F(201,200)*Cs[3]/a)*(2*a/3)**3
   r=(AA+yg.up()*C2)/isqrt(3*K**3)+(BB+yg.up()*C3)/isqrt(5*K**5)
   errtrace+=r*r
  norms=[I(t).sqrt().up() for t in (logtrace,primetrace,errtrace)]
  tail=rup(sum(norms)**2);rho=band+tail
  print('parity',k,'norms',list(map(float,norms)),'tail',tail,'total',rho,flush=True)
  assert rho<F(23,100)<F(1,4),(k,rho)
  results.append(dict(parity=k,N=N,K=K,rows=nrows,fit=list(map(str,fit)),band_operator_upper=str(band),band_radius_upper=str(rup(rad,10**30)),midpoint_pivot_lower=str(F(piv.numerator*10**12//piv.denominator,10**12)),log_norm_upper=str(rup(norms[0])),prime_norm_upper=str(rup(norms[1])),error_norm_upper=str(rup(norms[2])),tail_operator_upper=str(tail),rho_upper=str(rho),fitted_band_gram_intervals=[[[str(t.l),str(t.h)] for t in row] for row in FG]))
 result=dict(stage='CR-ALP-044',status='PASS conditional on predecessor continuum and interval inputs',sigma='1/100000000000000',eta='23/100',tail_margin='1/50',blocks=results,not_claimed='E172 finite positivity, intervening 158 dimensions, Lean, RH')
 (ROOT/'CR-ALP-044_exact_result.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS E14 plus Q172 mixed certificate, conditional on analytical premises',flush=True)
if __name__=='__main__':run()
