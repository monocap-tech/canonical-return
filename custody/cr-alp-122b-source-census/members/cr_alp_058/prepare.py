"""Exact witness contraction, discarded group norm, verified Gauss-Legendre root/weight brackets."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math
import mpmath as mp
from math import comb,isqrt
root=Path(__file__).resolve().parent; pre=root.parent/'cr_alp_057'; old=root.parent/'cr_alp_056'
files={'result':'e0dbb8cd24322cd5f0a96a7d8c1b188c38b7f5d2f706d644b8268cbb6533a809','grouped_proxy':'12c7fbc679f244f0016cd71acd7256cf00830bb40d14fb2fedf452cab7706c0c','trial':'6576fb333c7b2fdec6a4c236d70b3257b4bdedcb9208ae77a10e3e1931f512b0'}
d={}
for n,h in files.items():
 raw=(pre/(n+'.json')).read_bytes();assert hashlib.sha256(raw).hexdigest()==h;d[n]=json.loads(raw)
S=int(d['result']['scale']); assert S==10**100
# Exact decimal enclosures at 80 places; all downstream endpoints parsed directed.
D=10**80
def dec(n):return ('-' if n<0 else '')+str(abs(n)//D)+'.'+str(abs(n)%D).zfill(80)
def enc(x):return dec(math.floor(x*D))+' '+dec(math.ceil(x*D))
def encint(x):return enc(x[0]).split()[0]+' '+enc(x[1]).split()[1]
def mid(row,z):return F(sum(((int(l)+int(h))//2)*v for (l,h),v in zip(row,z)),S*10**12)
def interval_quad(A,z,sc):
 lo=hi=0
 for i,row in enumerate(A):
  for j,(a,b) in enumerate(row):
   v=z[i]*z[j];x,y=int(a)*v,int(b)*v;lo+=min(x,y);hi+=max(x,y)
 return F(lo,sc*10**24),F(hi,sc*10**24)
# Gauss roots proposed numerically, certified exactly by opposite signs of P_N.
N=96; mp.mp.dps=100
from scipy.special import roots_legendre
nodes=[]
for xx in roots_legendre(N)[0]:
 x=mp.mpf(float(xx))
 for _ in range(5):
  a,b=mp.mpf(1),x
  for j in range(2,N+1):a,b=b,((2*j-1)*x*b-(j-1)*a)/j
  x-=b*(1-x*x)/(N*(a-x*b))
 nodes.append(x)
# P_N coefficients: denominator 2^N, exact integers.
coeff=[0]*(N+1)
for j in range(N//2+1):coeff[N-2*j]=(-1)**j*comb(N,j)*comb(2*N-2*j,N)
def poly(co,x):
 # coefficients all have common denominator 2^N; integer Horner with homogeneous powers.
 m,q=x.numerator,x.denominator; v=0; qp=1
 for a in reversed(co):
  v=v*m+int(a*2**N)*qp; qp*=q
 return F(v,qp//q*2**N)
p=[F(a,2**N) for a in coeff]; dp=[(i+1)*p[i+1] for i in range(N)]; dd=sum(abs(i*dp[i]) for i in range(1,N))
Q=10**75; roots=[]
for x in nodes:
 m=int(mp.floor(x*Q)); l,h=F(m-1,Q),F(m+2,Q);assert -1<l<h<1 and poly(p,l)*poly(p,h)<0
 midx=(l+h)/2;der=abs(poly(dp,midx));er=dd*(h-l)/2;assert der>er
 ma=max(abs(l),abs(h));mi=min(abs(l),abs(h)) if l*h>0 else F(0)
 wl=2/((1-mi**2)*(der+er)**2);wh=2/((1-ma**2)*(der-er)**2);assert 0<wl<wh
 roots.append((l,h,F(math.floor(wl*D),D),F(math.ceil(wh*D),D)))
print('certified',len(roots),'Gauss nodes',flush=True)
assert all(roots[i][1]<roots[i+1][0] for i in range(N-1));assert sum(r[2] for r in roots)<2<sum(r[3] for r in roots)
(root/'gauss.txt').write_text(str(N)+'\n'+'\n'.join(encint((l,h))+' '+encint((wl,wh)) for l,h,wl,wh in roots)+'\n')
summary={}
for k in ['1','2']:
 r=d['result']['sectors'][k]; t=d['trial']['sectors'][k]; g=d['grouped_proxy']['sectors'][k]; raw=(old/f'result_extended_{k}.json').read_bytes();w=json.loads(raw);assert hashlib.sha256(raw).hexdigest()=={'1':'5340455da04a29619f943a4360caed34a1b9bdda56d68068b9f9c523af784325','2':'4017c7a656e43fe9fa97da91c37f35a9eda5a78dd3efaa2863b40de71177ac36'}[k]
 z=w['comparison_witness']['z'];norm2=F(sum(x*x for x in z),10**24);assert norm2<F(1000001,1000000)
 low={i:mid(a,z) for i,a in zip(r['source_ids'],r['source_sine_coefficients'])};trial={i:mid(a,z) for i,a in zip(t['ids'],t['D'])};f=mid(r['pole_coefficient'],z);db=mid(t['Db'][0],z)
 lines=[str(k),enc(f),enc(db),str(len(trial))]
 for i in sorted(trial):lines.append(str(i)+' '+enc(low.get(i,F(0)))+' '+enc(trial[i]))
 lines.append(str(len(g['low_amplitudes'])))
 for m,a in enumerate(g['low_amplitudes']):lines.append(enc(F(4*m+1,2))+' '+enc(mid(a,z)))
 groups=[x for x in g['groups'] if int(x['M'])<2**32]; omitted=[x for x in g['groups'] if int(x['M'])>=2**32]; lines.append(str(len(groups)))
 tail=F(0)
 for b in omitted:
  c=F(b['center']); # 2/sqrt(2c), upper via floor sqrt rational at 10^50
  inv=F(10**50,isqrt((2*c*10**100).numerator//(2*c*10**100).denominator));tail+=2*inv*sum(abs(mid(a,z)) for a in b['amplitudes'])
 for b in groups:
  lines.append(enc(F(b['center']))+' '+str(len(b['amplitudes'])))
  lines.extend(enc(mid(a,z)) for a in b['amplitudes'])
 (root/f'input_{k}.txt').write_text('\n'.join(lines)+'\n')
 K=interval_quad(w['K'],z,int(w['scale']));summary[k]={'norm2':str(norm2),'K_lower':str(K[0]),'K_upper':str(K[1]),'discarded_group_norm_upper':str(math.ceil(tail*10**12)/10**12),'discarded_group_norm_rational_upper':str(F(math.ceil(tail*10**12),10**12)),'witness_sha256':hashlib.sha256(raw).hexdigest()}; print(k,'K',float(K[0]),float(K[1]),'discarded',float(tail),flush=True)
(root/'prepare_result.json').write_text(json.dumps(summary,indent=2)+'\n')
