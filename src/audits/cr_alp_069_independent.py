"""Independent exact Fraction replay of CR-ALP-069; no numerical solve or quadrature."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib,json
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
d1={j:P(dnew(next(b for b in bands if b['start']<=j and (b['end'] is None or j<=b['end'])))) for j in ids}
for j in ids:
    b=next(b for b in reversed(oldbands) if b['j']<=j)
    assert d0[j]==P(floor9((1-theta)*gamma+theta*F(b['gamma_diagonal_lower'])))
p={j:r[f'pole{j}'] for j in range(2,N+1,2)}
y={j:r[f'y{j}'] for j in range(2,N+1,2)}
c=div(total(mul(x0[j],p[j]) for j in ids),total(sq(p[j]) for j in ids))
x={j:sub(x0[j],mul(c,p[j])) for j in ids}
X2=total(sq(x[j]) for j in ids)
assert overlap(X2,r['X2']) and overlap(total(mul(x[j],p[j]) for j in ids),Z)
delta={j:mul(x[j],sub(d1[j],d0[j])) for j in ids}
E=sub(r['E'],total(mul(x[j],delta[j]) for j in ids))
# Independent organization: first compute the old projected residual norm and
# projected finite coordinates, then apply the projected diagonal difference.
low=range(2,86,2);st=r['pole_st'];assert st[0]>F(3,1000)
B=sub(r['Yb'],total(mul(y[j],p[j]) for j in low))
lam=div(B,st)
Qold=sub(sub(r['Y2'],total(sq(y[j]) for j in low)),div(sq(B),st))
bd=total(mul(delta[j],p[j]) for j in ids)
lam_delta=div(bd,st)
v={j:sub(y[j],mul(lam,p[j])) for j in ids}
# ||Q delta||²=||delta||²-|<delta,b>|²/st;
# <Q y,Q delta>=<Q y,delta> by orthogonal projection.
Qnew=add(sub(Qold,mul(P(2),total(mul(v[j],delta[j]) for j in ids))),sub(total(sq(delta[j]) for j in ids),div(sq(bd),st)))
vnew={j:sub(v[j],sub(delta[j],mul(lam_delta,p[j]))) for j in ids}
tail=sub(Qnew,total(sq(vnew[j]) for j in ids))
assert tail[1]>=0
A=sum((sq(vnew[j])[1]/(d1[j][0]-s) for j in ids),F(0))+tail[1]/dtail
eta=F(1,10**13);eps=F(old['proxy_error_upper'])
proxy=E[1]+rootup(X2[1])*eta+(rootup(A)+eta*rootup(1/dmin))**2
bound=(rootup(proxy)+eps*rootup(1/dmin))**2
K=F(old['K'][0]);margin=K-bound
assert bound<K and margin>F(2,100)
# Reorganized interval arithmetic may yield different widths; both must pass.
assert F(saved['true_inverse_upper'])<K
print('PASS independent projected-difference replay')
print('inverse upper',float(bound),'K',float(K),'margin',float(margin))
print(json.dumps({'method':'projected-difference Fraction replay','inverse_upper':str(bound),'K':str(K),'margin':str(margin),'status':'PASS_WITNESS'},indent=2))
