"""RETURN095-B: independently reconstruct unprojected negative-update
Gram for the original baseline complement and test the inherited .52
norm threshold. No new source-form identification is assumed.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def load(p,h):
 b=(root/p).read_bytes();assert hashlib.sha256(b).hexdigest()==h,p;return json.loads(b)
w=load('return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
u=load('return057/RETURN057_weights.json','090ec051ef390b339357496b805808a846268352c43b70befb1c3ede3614a926')
good=sorted((k,v) for k,v in w['positive_even_weights'] if k>116)
K=32768;A=F(87,500);S=10**50
assert len(good)==16302 and good[-1][0]==K-2
assert len({k for k,v in good})==len(good) and all(k%2==0 and k>116 and v>0 for k,v in good)
d={j:F(v,10**6) for j,v in u['odd_deficit'] if v}
js=sorted(d);assert len(js)==32 and all(j%2 and 0<j<K and d[j]>0 for j in js)
def atan(q,n):
 a=sum((F((-1)**i,(2*i+1)*q**(2*i+1)) for i in range(n)),F(0));b=a+F((-1)**n,(2*n+1)*q**(2*n+1))
 return min(a,b),max(a,b)
a,b=atan(5,80);c,e=atan(239,32);pi=16*a-4*e;assert pi>3
def ceil(q):return -((-q.numerator)//q.denominator)
def enc(q):return q.numerator*S//q.denominator,ceil(q*S)
def floor(q):return q.numerator//q.denominator
# H_j=sum A_k^-1/(k²-j²), J_j=sum A_k^-1/(k²-j²)^2.
# Exact per-term directed rounding, no preexisting scalar bank.
H={};J={}
for j in js:
 hl=hu=jl=ju=0
 for k,v in good:
  q=k*k-j*j;assert q
  a,b=enc(F(2000000,(198000+v)*q));hl+=a;hu+=b
  a,b=enc(F(2000000,(198000+v)*q*q));jl+=a;ju+=b
 H[j]=(F(hl,S),F(hu,S));J[j]=(F(jl,S),F(ju,S))
print('finite signed scalar rows',len(good)*len(js),flush=True)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def scale(x,c):return (x[0]*c,x[1]*c) if c>=0 else (x[1]*c,x[0]*c)
def sub(x,y):return add(x,scale(y,-1))
def tail(i,j,par):
 # all tail rows positive: h_i h_j = c k^-p /[(1-i²/k²)(1-j²/k²)].
 # decreasing power sum bounded by integral and first cell, half density.
 p=4-2*par;c=F(i*j if par==0 else 1)/A
 q=F(1,2*(p-1)*K**(p-1))
 f=1/(1-F(i*i,K*K))/(1-F(j*j,K*K))
 return c*q,c*f*(q+F(1,K**p))
def gram(i,j,par):
 if i==j:
  a=scale(J[i],i*i)
  if par:a=add(H[i],a)
 elif par:a=scale(sub(scale(H[i],i*i),scale(H[j],j*j)),F(1,i*i-j*j))
 else:a=scale(sub(H[i],H[j]),F(i*j,i*i-j*j))
 return add(a,tail(i,j,par))
def pivots(A):
 n=len(A);a=[row[:] for row in A];prev=1
 for k in range(n):
  t=a[k][k]
  if t<=0:return k
  for i in range(k+1,n):
   for j in range(i,n):
    q,r=divmod(a[i][j]*t-a[i][k]*a[k][j],prev);assert r==0
    a[i][j]=a[j][i]=q
  prev=t
 return n
records=[]
for par in (0,1):
 G=[[gram(i,j,par) for j in js] for i in js]
 assert G==[list(row) for row in zip(*G)]
 # Certify rho diag(pi_lower²/(8d)) - G positive by exact
 # integer lower matrix: midpoint minus max row radius.
 # This is equivalent to unprojected normalized negative Gram < rho.
 M=[[floor((lo+hi)*S) for lo,hi in row] for row in G]
 rad=[[ceil((hi-lo)*S)+1 for lo,hi in row] for row in G]
 rho=F(13,25)
 for attempt in (rho,F(13,20)):
  N=[[-x for x in row] for row in M]
  for a,j in enumerate(js):
   N[a][a]+=floor(2*S*attempt*pi*pi/(8*d[j]))-sum(rad[a])
  n=pivots(N);print('negative-update bound',par,str(attempt),'positive pivots',n,flush=True)
  if n==len(js):break
 records.append({'parity':par,'bound':str(attempt),'positive_pivots':n,'target_pass':attempt==rho and n==32,'gram':[[[str(v) for v in a] for a in row] for row in G]})
(out/'receipt.json').write_text(json.dumps({'scope':'unprojected baseline negative-update bound on weighted grid; source-form comparator attachment remains open','status':'PASS' if all(r['target_pass'] for r in records) else 'TARGET_NOT_CERTIFIED','finite_rows':len(good),'negative_rows':len(js),'records':records},indent=2)+'\n')
assert all(r['target_pass'] for r in records),'inherited .52 threshold not independently certified'
print('RETURN095-B BOTH NEGATIVE UPDATE NORM BOUNDS PASS',flush=True)
