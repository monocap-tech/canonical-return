"""RETURN096-B: check full original signed-symbol ledger against
period-five comparison. Alternate rational endpoint arithmetic from
RETURN046 replay and the analytic lower bound in ATTACHMENT.md.
No assumption of the finite saddle's positivity.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p,h):
 f=root/p;assert sha(f)==h,p;return json.loads(f.read_bytes())
w=load('return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
u=load('return057/RETURN057_weights.json','090ec051ef390b339357496b805808a846268352c43b70befb1c3ede3614a926')
raw=load('return059/RETURN059_raw_receipt.json','2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8')
old=load('return046/RETURN046_periodic_receipt.json','34a0f40975e097c50fa7205d90715565a01972edcaf1614c9f42c4d6b7752d42')
pr=load('return077/proposal.json','803aa2124522c5565864c35623d92fce2eb1f1fab616347a6020387f65068120')
restfile=root/'return075/RETURN075_receipt.json'
rest=json.loads(restfile.read_bytes());assert all(sha(root/p)==h for p,h in rest['pins'].items())
S=10**20

def fl(x):return x.numerator*S//x.denominator
def ce(x):return -fl(-x)
def atan(x):
 z=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(32)),F(0));return z,z+x**65/65
l,u=atan(F(1,5));a,b=atan(F(1,239));pl=fl(16*l-4*b);pu=ce(16*u-4*a)
assert 3*S<pl<pu<22*S//7

def log_series(x):
 y=(x-1)/(x+1);t=y;v=F(0)
 for j in range(36):v+=2*t/(2*j+1);t*=y*y
 return v,v+2*t/(73*(1-y*y))
l2,u2=log_series(F(2))
def log_bounds(x):
 x=F(x);k=0
 while x>=2:x/=2;k+=1
 assert x>=1
 a,b=log_series(x);return a+k*l2,b+k*u2
assert log_bounds(11)[1]<F(5,2)<log_bounds(13)[0]
# Reconstruct actual prime powers by integer trial division.
W=[]
for n in range(2,13):
 p=next((i for i in range(2,n+1) if n%i==0));m=n
 while m%p==0:m//=p
 if m!=1:continue
 l,u=log_bounds(n);a,b=log_bounds(p);z=math.isqrt(n*S*S)
 W.append((n,fl(l),ce(u),fl(a/F(z+1,S)),ce(b/F(z,S))))
assert [x[0] for x in W]==[2,3,4,5,7,8,9,11]
tri=2*sum(w[-1] for w in W)
# Independent lower log: fixed point positive atanh using recurrence of powers.
# Every truncation downward is a lower bound, no remainder required.
def log_lower(a,b):
 k=a.bit_length()-b.bit_length();k=max(k,0)
 if a<b*2**k:k-=1
 b*=2**k;assert b<=a<2*b
 t=(a-b)*S//(a+b);y=t;v=0
 for j in range(20):
  v+=2*t//(2*j+1);t=t*y*y//(S*S)
 return v+k*fl(l2)
# Forward alternating Taylor with rounded-down absolute terms.
# At |x|<4, sum of term errors < 20*cosh(4) fixed-point units; allow 1000.
rem=ce(F(4**42,math.factorial(42)))
def cosine(x,err):
 assert abs(x)<4*S
 t=S;z=S
 for j in range(1,21):
  t=t*x*x//(S*S*(2*j-1)*(2*j));z+=(-1)**j*t
 e=1000+rem+err
 return max(-S,z-e),min(S,z+e)
def pupper(k):
 v=0
 for n,l,u,wl,wu in W:
  low=k*l//5;high=-((-k*u)//5);q=(low+S//2)//S
  a=low-q*S;b=high-q*S
  # 2pi times the interval; round down and up respectively.
  corners=[2*p*t for p in (pl,pu) for t in (a,b)]
  x=min(corners)//S;upper=-((-max(corners))//S)
  cl,cu=cosine(x,upper-x)
  v+=max(wl*cl,wl*cu,wu*cl,wu*cu)
 return -((-2*v)//S)
def glo(k):
 if k<=1:return -6*S
 a=2*pl*k;b=5*S
 return max(-6*S,log_lower(a,b)-185*S//100-(-(-2*S*b//a)))

# Recompute a tighter small/medium frequency Gamma lower bound:
# g0 from Euler constant upper and 32 positive terms plus integral.
Hup=F(sum(-(-S//i)for i in range(1,4097)),S)
g0lo=-(Hup-log_bounds(4096)[0])-F(pu,2*S)-3*log_bounds(2)[1]-log_bounds(F(pu,S))[1]
assert g0lo>F(-5373,1000)
def gamma(k):
 x=(F(2*pl*k,5*S))**2
 v=F(-5373,1000);B=F(129,2)
 for j in range(32):
  beta=F(4*j+1,2);t=2*x/(beta*(beta*beta+x));v+=F(fl(t),S)
 v+=F(log_lower((1+x/B**2).numerator,(1+x/B**2).denominator),2*S)
 return max(v,F(glo(k),S))
def symbol(k):return gamma(k)-F(pupper(k),S)
def exp_lower(x):
 t=s=F(1)
 for j in range(1,61):t*=x/j;s+=t
 return s
image=1/exp_lower(F(25,4))/(5*(1-1/exp_lower(F(25,2)))*(1-1/exp_lower(F(5))))
assert image<F(1,1000)
assert len(W)==8
# Finite diagonal ledger: use original weights (not upgraded retained weights).
good={k:v for k,v in w['positive_even_weights'] if k>116}
assert len(good)==16302
min_even=None
for k,v in sorted(good.items()):
 lower=F(glo(k)-pupper(k),S)-F(1,10)-F(v,10**6)
 assert lower>=0,(k,str(lower))
 min_even=lower if min_even is None else min(min_even,lower)
print('finite even weights',len(good),'min slack',float(min_even),flush=True)
# All odd modes, not merely the positive selected rows. Sum signed and
# restoration weights where they overlap; no double spending allowed.
pos=dict(u['odd_good']);neg=dict(u['odd_deficit']);rv=dict(rest['weights'])
assert len(rv)==len(rest['weights'])==500 and set(rv)==set(pr['js'])
assert all(v>0 for v in rv.values()) and all(v>0 for v in pos.values()) and all(v>=0 for v in neg.values())
J={j for j in w['unpaid_odd'] if j>0}
assert not (set(pos)|set(neg)|set(rv))&J
assert all(j%2 and 0<j<32768 for j in set(pos)|set(neg)|set(rv))
bad=[];min_odd=None
for j in range(1,32768,2):
 if j in J:continue
 v=pos.get(j,0)-neg.get(j,0)+rv.get(j,0)
 # Old lower bound is sufficient for most of the omitted odd rows.
 # Upgraded lower bound is used on selected low/extra carriers.
 lower=F(glo(j)-pupper(j),S)-F(1,10)-F(v,10**6)
 if lower<0:lower=symbol(j)-F(1,10)-F(v,10**6)
 if lower<0:bad.append((j,v,str(lower)))
 min_odd=lower if min_odd is None else min(min_odd,lower)
print('odd full ledger',16384-len(J),'restoration overlap',len(set(rv)&set(pos)),'bad',len(bad),'first',bad[:10],flush=True)
# Monotonic g lower minus the fixed prime triangle upper gives each
# infinite-band bound. A = .099 + half of the even surplus.
tails=[]
for k,A in [(32768,F(87,500)),(65536,F(519,1000)),(262144,F(1199,1000))]:
 l=F(glo(k)-tri,S);target=F(1,10)+2*(A-F(99,1000));assert l>=target,(k,str(l),str(target))
 tails.append({'cutoff':k,'A':str(A),'symbol_lower':str(l),'required':str(target)})
 print('monotone tail weight',k,str(A),float(l),flush=True)
assert F(glo(32768)-tri,S)>=F(1,10)
rec={'status':'PASS' if not bad else 'REFUSED','scope':'exact signed-symbol ledger check on periodic zero-extension comparison; analytic derivation in ATTACHMENT.md; not Lean or arbitrary support','even_count':len(good),'odd_count':16384-len(J),'restoration_overlap':len(set(rv)&set(pos)),'odd_failures':bad,'min_even_slack':str(min_even),'min_odd_slack':str(min_odd),'gamma_image_upper':str(image),'tail':tails,'restoration_receipt_sha256':sha(restfile)}
(out/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
assert not bad,'actual signed-symbol ledger does not certify all odd weights'
print('RETURN096-B COMPLETE SIGNED SYMBOL LEDGER PASS',flush=True)
