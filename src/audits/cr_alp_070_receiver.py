"""Independent exact Fraction replay of CR-ALP-069; no numerical solve or quadrature."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json,sys
root=Path(__file__).resolve().parent
ev=root/'cr_alp_067/evidence'
hashes={
'weighted_trial.txt':'27670ee547ec2f5fc30fd3adb155378267e82df6f49bbaf90701c31f1fe4c15e',
'weighted_diag.txt':'38024951a1daf15f1a19c0298d793443458e9033e5b5620af3178a8378728126',
'weighted_coeff.log':'40b6f178ae1b5d12aab9276252bff6b2c63d673c7ceeaac36a8c66d73a9ef206',
'weighted_response.log':'d8c54082f7135d1c7192b3b46bc8df00e5e1482b66d82fe0c8140bfa10cfe93b',
'weighted_response_result.json':'80fb1d55980de033e53bd2a2ed2bb16c222d6639b331914242d5f283ecf53393'}
for file,sha in hashes.items():
    assert hashlib.sha256((ev/file).read_bytes()).hexdigest()==sha,file
print('PASS all five original weighted artifact hashes')
def load(p):return json.loads(p.read_text())
bands=load(root/'CR-ALP-068_pole_centered_gamma_result.json')['bands']
oldbands=load(root/'CR-ALP-066_gamma_diagonal_result.json')['bands']
old=load(root/'cr_alp_058_recovery/verification_result.json')['2']
saved=load(root/'CR-ALP-069_transfer_response_result.json')
N=2048; ids=list(range(86,N+1,2))
def P(x):return (F(x),F(x))
Z=P(0)
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[1],a[1]-b[0])
def mul(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)
def div(a,b):
    assert b[0]>0;return mul(a,(1/b[1],1/b[0]))
def sq(a):
    return (F(0) if a[0]<=0<=a[1] else min(x*x for x in a),max(x*x for x in a))
def total(xs):
    lo=hi=F(0)
    for l,h in xs:lo+=l;hi+=h
    return lo,hi
def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def rootup(x):
    assert x>=0
    scale=10**30
    return F(isqrt((x*scale**2).numerator//(x*scale**2).denominator)+1,scale)
def bank(name):
    out={}
    for j,l,h in (row.split() for row in (ev/name).read_text().splitlines()):
        j=int(j);assert j not in out and F(l)==F(h)
        out[j]=P(l)
    assert list(out)==ids
    return out
x0=bank('weighted_trial.txt');d0=bank('weighted_diag.txt')
r={}
for name,l,h,slash,sc in (row.split() for row in (ev/'weighted_response.log').read_text().splitlines()):
    assert slash=='/' and sc=='1e15' and name not in r
    r[name]=(F(int(l),10**15),F(int(h),10**15));assert r[name][0]<=r[name][1]
assert set(r)=={'G','X2','E','Y2','Yb','pole_st',*(f'y{j}' for j in range(2,N+1,2)),*(f'pole{j}' for j in range(2,N+1,2))}
assert overlap(r['G'],tuple(map(F,old['proxy_G'])))
assert bands[0]['start']==86 and bands[-1]['end'] is None
assert all(b['start']==a['end']+1 for a,b in zip(bands,bands[1:]))
def floor9(x):return F(x*10**9//1,10**9)
def dnew(b):return floor9(F(b['gamma_lower']))
s=F(21,8);gamma=F(72,25);theta=F(1,10)
dmin=min(dnew(b)-s for b in bands)
dtail=min(dnew(b)-s for b in bands if b['end'] is None or b['end']>N)
assert dmin>0 and dtail>=dmin
cap=next(dnew(b) for b in bands if b['start']==2049)
assert cap-s==dtail
assert all(dnew(b)>=gamma for b in bands)
d1={j:P(min(cap,dnew(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end']))))) for j in ids}
for j in ids:
    b=next(b for b in reversed(oldbands) if b['j']<=j)
    assert d0[j]==P(floor9((1-theta)*gamma+theta*F(b['gamma_diagonal_lower'])))
p={j:r[f'pole{j}'] for j in range(2,N+1,2)}
y={j:r[f'y{j}'] for j in range(2,N+1,2)}
c=div(total(mul(x0[j],p[j]) for j in ids),total(sq(p[j]) for j in ids))
x={j:sub(x0[j],mul(c,p[j])) for j in ids}
X2=total(sq(x[j]) for j in ids)
assert overlap(X2,r['X2']) and overlap(total(mul(x[j],p[j]) for j in ids),Z)

# The archived input carries the contracted source and correction coefficients.
# Treat the archived correction as an interval enclosure for a proxy, not
# automatically the exact original X'. The final source attachment is separate.
package=Path(sys.argv[1]).resolve()
tokens=iter((package/'cr_alp_058/input_2.txt').read_text().split())
assert int(next(tokens))==2
def rd():
    l,h=F(next(tokens)),F(next(tokens));assert l<=h
    return l,h
f=rd();db=rd();n=int(next(tokens));assert n==425
t={}
for i in range(n):
    j=int(next(tokens));assert j==2*(i+1)
    source=rd();t[j]=rd()
for j in ids:
    if j not in t:t[j]=Z
st=r['pole_st'];assert st[0]>F(3,1000)
# Q85 is high-sine projection followed by removal of the remaining pole.
# Project the archived finite sine bank, which kills its explicit db*b term.
tau=div(total(mul(t[j],r[f'pole{j}']) for j in ids),st)
# The archived coefficient db should equal -tau if its correction is pole-null.
assert overlap(db,sub(Z,tau)),('pole bridge',db,tau)
U={j:sub(t[j],mul(tau,r[f'pole{j}'])) for j in ids}
# A = Q85 diag(D_new-gamma) Q85, with a constant cap for all j>2048.
# Q85 A_raw X_proxy = Q85 z, where z is finite and
# z_j=A_j t_j -(A_j-A_cap) tau c_j.
Acap=cap-gamma
assert Acap>0
A={j:sub(d1[j],P(gamma)) for j in ids}
assert all(0<=v[0]<=v[1]<=Acap for v in A.values())
z={j:sub(mul(A[j],t[j]),mul(sub(A[j],P(Acap)),mul(tau,r[f'pole{j}']))) for j in ids}
AU=total(mul(U[j],z[j]) for j in ids)
assert AU[1]>=0
# Archived correction norm <= norm of its finite high sine precursor.
Unorm=rootup(total(sq(t[j]) for j in ids)[1])
delta={j:mul(x[j],sub(d1[j],d0[j])) for j in ids}
change={j:add(delta[j],z[j]) for j in ids}
E=sub(sub(r['E'],total(mul(x[j],delta[j]) for j in ids)),mul(P(2),total(mul(x[j],z[j]) for j in ids)))
# y_new = y_old - delta - z before the final Q85 projection.
Y2=add(sub(r['Y2'],mul(P(2),total(mul(r[f'y{j}'],change[j]) for j in ids))),total(sq(change[j]) for j in ids))
Yb=sub(r['Yb'],total(mul(change[j],r[f'pole{j}']) for j in ids))
B=sub(Yb,total(mul(r[f'y{j}'],r[f'pole{j}']) for j in range(2,86,2)))
lam=div(B,st)
Q=sub(sub(Y2,total(sq(r[f'y{j}']) for j in range(2,86,2))),div(sq(B),st))
v={j:sub(sub(r[f'y{j}'],change[j]),mul(lam,r[f'pole{j}'])) for j in ids}
tail=sub(Q,total(sq(v[j]) for j in ids))
assert tail[1]>=0
penalty=sum((sq(v[j])[1]/(d1[j][0]-s) for j in ids),F(0))+tail[1]/dtail
eta=F(1,10**13)
# Test the explicit conditional attachment tolerance ||X'-X_proxy|| <= 1e-10.
# No assumption that this attachment has already been established.
deltaX=F(1,10**10)
eps=F(old['proxy_error_upper'])+Acap*deltaX
K=F(old['K'][0])+max(F(0),AU[0]-Acap*deltaX*(2*Unorm+deltaX))
proxy=E[1]+rootup(X2[1])*eta+(rootup(penalty)+eta*rootup(1/dmin))**2
assert proxy>=0
bound=(rootup(proxy)+eps*rootup(1/dmin))**2
margin=K-bound
# Output compact outward rational endpoints to 1e-18, avoiding giant Fraction serialization.
scale=10**18
def lo(x):return str(F(x*scale//1,scale))
def hi(x):return str(-F((-x)*scale//1,scale))
result={'source_attachment':'CONDITIONAL: ||X_original-X_archived_projected|| <= 1e-10','status':'PASS_CONDITIONAL_RECEIVER' if margin>0 else 'UNDECIDED','cap':str(cap),'dmin':str(dmin),'dtail':str(dtail),'pole_bridge_overlap':True,'correction_norm_upper':hi(Unorm),'correction_energy_lower':lo(AU[0]),'correction_energy_upper':hi(AU[1]),'K_C_lower':lo(K),'energy_upper':hi(E[1]),'Q_residual_norm2_upper':hi(Q[1]),'tail_norm2_upper':hi(tail[1]),'weighted_residual_upper':hi(penalty),'inverse_upper':hi(bound),'margin_lower':lo(margin)}
print(json.dumps(result,indent=2))
