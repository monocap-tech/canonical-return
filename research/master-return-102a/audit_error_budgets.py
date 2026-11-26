"""RETURN102-A independent exact receiving error-budget recomputation.
All norms evaluated after the signed retained chart. No inherited error
norm or remainder value used in constructing the six bounds.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
Q=10**80;N=32
def load(p,h):
 b=(root/p).read_bytes();assert hashlib.sha256(b).hexdigest()==h,p
 return json.loads(b)
raw=load('return059/RETURN059_raw_receipt.json','2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8')
mom=load('return066/RETURN066_moments.json','b9aa229e7147b8d47d56b0cbfec290255af557cfb99260624d1c4b576369d184')
w=load('return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
signed=load('return057/RETURN057_weights.json','090ec051ef390b339357496b805808a846268352c43b70befb1c3ede3614a926')
b=load('return062/RETURN062_bank_receipt.json','429cb13fe51f9f92e684fffb9db6d6eb31004c06e20a81dc1877c94065f2cef2')
pl,pu=map(F,b['pi']);assert F(25,8)<pl<pu<F(22,7)
def ceil(x):return math.ceil(x)
def E(l,u=None):
 if u is None:u=l
 return math.floor(l*Q),ceil(u*Q)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
 t=[x*y for x in a for y in b];return min(t)//Q,-((-max(t))//Q)
def mag(a):return F(max(map(abs,a)),Q)
def sqrtup(x):
 assert x>=0
 S=10**40;n=math.isqrt(x.numerator*S*S//x.denominator)
 if F(n*n,S*S)<x:n+=1
 return F(n,S)
def norm(rows,C,T):
 z=[max(abs(sum(min(l*c,u*c) for (l,u),c in zip(rows,col))),abs(sum(max(l*c,u*c) for (l,u),c in zip(rows,col)))) for col in zip(*C)]
 return sqrtup(F(sum(x*x for x in z),Q*Q*T*T))
good=dict(w['positive_even_weights']);good[32768]=150000
assert all(k in good for k in range(16386,32769,2))
# Upward-rounded reciprocals avoid an enormous common denominator.
weight_sum=sum(F(ceil(F(2000000*10**30,198000+good[k])),10**30) for k in range(16386,32769,2))
def moment(p,q,K):return F(K*math.factorial(q),2*(p-1)**(q+1))+2*math.factorial(q)
def trace(par,K,A):
 total=F(0);p=4-2*par
 for j,d in signed['odd_deficit']:
  if not d:continue
  assert 0<j<K
  # Strict k>K decreasing power integral, geometric denominator bound.
  total+=F(8*d,10**6)/pl**2/A*F(j*j if par==0 else 1,2*(p-1)*K**(p-1))/(1-F(j*j,K*K))**2
 return total
records=[];results=[]
for par in (0,1):
 C=raw['blocks'][par]['C'];T=raw['T'];floor=F(raw['blocks'][par]['margin']);assert floor>0
 rn=[sqrtup(F(sum(x*x for x in row),T*T)) for row in C]
 cf=[[E(F(l,10**55),F(u,10**55)) for l,u in row] for row in mom['blocks'][par]['raw_coeff']]
 assert len(cf)==len(C)==87
 caps=[]
 for K,n0 in [(16384,12),(32768,8),(262144,1)]:
  u=F(4,5*K);inv=E(F(5,2*K)**2/pu**2,F(5,2*K)**2/pl**2)
  rat=[mul((-r[0][1],-r[0][0]),inv) for r in cf];pr=mul(E(-F(1,4)),inv)
  rho=max(map(mag,rat+[pr]));assert rho<1
  pw=[E(F(1))]*87;pp=E(F(1));half=E(F(1,2));ns=[]
  for n in range(N):
   vv=[[mul(row[j],p) if j!=3 else add(mul(row[j],p),mul(mul(row[4],pp),half)) for row,p in zip(cf,pw)] for j in (1,2,3)]
   ns.append([norm(a,C,T) for a in vv])
   pw=[mul(a,b) for a,b in zip(pw,rat)];pp=mul(pp,pr)
  rem=[]
  for j in (1,2,3):
   e=sum(r*mag(mul(row[j],p)) for row,p,r in zip(cf,pw,rn))
   if j==3:e+=sum(r*mag(mul(mul(row[4],pp),half)) for row,r in zip(cf,rn))
   rem.append(e/(1-rho))
  B,V,L=[sum(r[j] for r in ns)+rem[j] for j in range(3)]
  bt,vt,lt=[sum(r[j] for r in ns[n0:])+rem[j] for j in range(3)]
  delta=8*u+3*u*u
  if K<262144:
   err=(6*B*u*u+delta*V*u*u+11*bt*u+23*vt*u*u+2*lt*u*u) if not par else (delta*B*u+6*V*u**3+23*bt*u+11*vt*u*u+2*lt*u)
   W=weight_sum if K==16384 else F((262144-32768)//2,F(87,500))
   square=F(4,5)*F(25,12)*W*err**2/floor
   tau=None
  else:
   A=F(87,500) if not par else F(1199,1000)
   if not par:terms=[(11*bt*u,3,0),((6*B+20*V+2*L)*u*u,2,0),(V*u*u,2,1)]
   else:
    B0,V0,L0=ns[0];G=10*B0+2*L0
    terms=[((20*bt+2*lt)*u,3,0),(bt*u,3,1),(11*V*u*u,2,0),((3*B0+6*V+G/4)*u**3,3,0)]
   tau=trace(par,K,A)
   square=F(4,5)*(1+tau/F(12,25))/A/floor*sum(a*b*moment(p+r,q+s,K) for a,p,q in terms for b,r,s in terms)
  cap=sqrtup(square);limit=F(23,10**6) if K==16384 else F(1,10000) if K==32768 else F(1,200) if not par else F(24,10000)
  assert cap<limit,(par,K,float(cap),float(limit));caps.append(cap)
  records.append({'parity':par,'K':K,'N':n0,'squared_error_upper':str(F(ceil(square*10**50),10**50)),'norm_upper':str(cap),'target':str(limit),'trace_upper':None if tau is None else str(F(ceil(tau*10**40),10**40))})
  print('independent error budget',par,K,float(cap),'target',float(limit),flush=True)
 r=(F(99,100) if not par else F(39,40))+sum(caps);assert r<1
 masscap=F(5*87,2)*F(sum(x*x for row in C for x in row),T*T)
 physical=(1-r)*min(floor/masscap,F(12,25)*F(99,1000))
 assert physical>=F(3,10**59),(par,float(physical))
 results.append({'parity':par,'error_norm_sum':str(sum(caps)),'true_squared_upper':str(r*r),'relative_floor':str(1-r),'conditional_physical_floor':str(physical)})
 print('independent conditional physical floor',par,float(physical),flush=True)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'independent six model-error budgets and common conditional physical floor; original raw coefficients and comparator; no unconditional proof promotion','records':records,'results':results},indent=2)+'\n')
print('RETURN102-A ALL SIX ERROR BUDGETS AND CONDITIONAL FLOOR PASS',flush=True)
