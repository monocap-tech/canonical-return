"""RETURN094-A: exact restoration metric factor audit, staged independently.
No terminal or scalar Gram computation. Run only after RETURN093-B closes.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
proposal=root/'return077/proposal.json'
assert sha(proposal)=='803aa2124522c5565864c35623d92fce2eb1f1fab616347a6020387f65068120'
cert=root/'return077/cert.py'
assert sha(cert)=='cc9ef1bc6478875ed304a323769a7ced90f80f0da9ae05ac73ab008bfdc6379b'
print('authenticated pi and weight definitions',flush=True)
for line in cert.read_text().splitlines():
 if any(t in line for t in ('pi1','weights','R=','oldfile','pins=')):print(line,flush=True)
p=load(proposal);weights_file=root/'return075/RETURN075_receipt.json';weights=dict(load(weights_file)['weights'])
js=p['js'];B=p['rows'];W=p['scale']
assert len(js)==len(set(js))==500 and len(B)==8 and all(len(x)==500 and sum(x)==0 for x in B)
assert all(weights[j]>0 for j in js)
# Compute exact signed weighted products without reading the old metric.
M=[[sum((F(int(x)*int(y)*10**6,W*W*8*weights[k]) for x,y,k in zip(b,c,js)),F(0)) for c in B] for b in B]
assert M[0][0]>0
# Locate the literal RETURN077 certificate by its authenticated proposal
# hash and exact shape; no guessed receipt filename.
candidates=[]
for f in sorted((root/'return077').glob('*.json')):
 d=load(f)
 if d.get('proposal_sha256')==sha(proposal) and 'R' in d:candidates.append((f,d))
assert len(candidates)==1,[f.name for f,d in candidates]
f,d=candidates[0];R=[[F(v) for v in row] for row in d['R']]
assert len(R)==8 and all(len(row)==8 for row in R)
assert d['rank']==8
# Recover the scalar from one exact entry; prove it is a positive
# rational square, then audit the other entries independently.
factor=R[0][0]/M[0][0]
n,m=math.isqrt(factor.numerator),math.isqrt(factor.denominator)
assert n*n==factor.numerator and m*m==factor.denominator
pi1=F(n,m);assert pi1>0
# Match the factor to an authenticated original pi endpoint, if present.
pi_bank=root/'return062/RETURN062_bank_receipt.json'
pl,pu=map(F,load(pi_bank)['pi'])
print('restoration pi scalar equals RETURN062 endpoint',pi1==pl,pi1==pu,flush=True)
# Independent rigorous Machin bound. Alternating series tails are
# enclosed by consecutive partial sums.
def atan(q,N):
 s=sum((F((-1)**i,(2*i+1)*q**(2*i+1)) for i in range(N)),F(0))
 t=s+F((-1)**N,(2*N+1)*q**(2*N+1))
 return min(s,t),max(s,t)
a,b=atan(5,160);c,e=atan(239,64)
lo,hi=16*a-4*e,16*b-4*c
assert lo<hi
assert pi1>=hi,('metric pi not independently verified as upper bound',str(pi1))
# Exact rational factorization using 500 weights. Positive weights
# and full row rank imply strict positive definiteness.
def rank(A):
 A=[[F(v) for v in row] for row in A];r=0
 for j in range(len(A[0])):
  piv=next((i for i in range(r,len(A)) if A[i][j]),None)
  if piv is None:continue
  A[r],A[piv]=A[piv],A[r];v=A[r][j];A[r]=[x/v for x in A[r]]
  for i in range(r+1,len(A)):
   v=A[i][j]
   if v:A[i]=[x-v*y for x,y in zip(A[i],A[r])]
  r+=1
  if r==len(A):break
 return r
assert rank(B)==8
assert all(R[i][j]==pi1*pi1*M[i][j] for i in range(8) for j in range(8))
record={'status':'PASS','scope':'exact weighted restoration metric factorization and pi upper bound; factor recovered from one archived entry, remaining 63 checked; independent row rank','entries':64,'independent_cross_checks':63,'rank':8,'positive_weights':500,'pi_factor_matches_062_upper':pi1==pu,'pi_factor':str(pi1),'pins':{str(v.relative_to(root)):sha(v) for v in (proposal,cert,weights_file,f,pi_bank)}}
(out/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
print('RETURN094-A RESTORATION METRIC FACTORIZATION PASS',flush=True)
