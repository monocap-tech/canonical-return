"""Exact rational endpoint audit of full physical response, including infinite sine tail."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json,sys
p=Path(sys.argv[1]); bands=json.loads(Path(sys.argv[2]).read_text())['bands']
old=json.loads(Path(sys.argv[3]).read_text())['2']
theta=F(sys.argv[4]);assert theta in (F(0),F(1,10));N=2048
def parse(path):
    out={}
    for line in path.read_text().splitlines():
        v=line.split()
        if len(v)==5 and v[-2:]==['/','1e15']:
            assert v[0] not in out
            out[v[0]]=tuple(F(int(t),10**15) for t in v[1:3])
    assert all(l<=h for l,h in out.values())
    return out
r=parse(p/'response.log')
assert set(r)=={'X2','G','E','Y2','Yb','pole_st',*(f'y{j}' for j in range(2,N+1,2)),*(f'pole{j}' for j in range(2,N+1,2))}
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    z=[x*y for x in a for y in b];return min(z),max(z)
def div(a,b):
    assert b[0]>0;return mul(a,(1/b[1],1/b[0]))
def sq(a):return (F(0) if a[0]<=0<=a[1] else min(x*x for x in a),max(x*x for x in a))
def rootup(x):
    assert x>=0;scale=10**18
    return F(isqrt((x*scale**2).numerator//(x*scale**2).denominator)+1,scale)
def dg(j):
    b=max((b for b in bands if b['j']<=j),key=lambda b:b['j'])
    return F(((1-theta)*F(72,25)+theta*F(b['gamma_diagonal_lower']))*10**9//1,10**9)
di=(p/'diag.txt').read_text().splitlines()
assert len(di)==(N-84)//2
for row,j in zip(di,range(86,N+1,2)):
    i,l,h=row.split();assert int(i)==j and F(l)==F(h)==dg(j)
s=F(21,8);dmin=dg(86)-s;dtail=dg(N+1)-s
assert dmin>0 and dtail>=dmin
assert all(dg(j+2)>=dg(j) for j in range(86,N,2))
assert r['pole_st'][0]>F(3,1000)
# Q85 projection of the unprojected residual y: remove low sines and orthogonal pole bridge.
Y=r['Y2'];B=r['Yb']
for j in range(2,86,2):
    Y=sub(Y,sq(r[f'y{j}']));B=sub(B,mul(r[f'y{j}'],r[f'pole{j}']))
lam=div(B,r['pole_st'])
Q=sub(Y,div(sq(B),r['pole_st']))
assert Q[1]>=0
weighted=F(0);tail=Q
for j in range(86,N+1,2):
    c=sub(r[f'y{j}'],mul(lam,r[f'pole{j}']))
    v=sq(c);weighted+=v[1]/(dg(j)-s);tail=sub(tail,v)
assert tail[1]>=0
A=weighted+tail[1]/dtail
# True proxy vs rational-shift response: ||error|| <= 1e-13.
shift=F(1,10**13)
assert 0<=r['X2'][1]<100
# Use an upper bound on 1/sqrt(dmin), NOT reciprocal of upward sqrt.
invroot=rootup(1/dmin)
proxy_bound=r['E'][1]+rootup(r['X2'][1])*shift+(rootup(A)+shift*invroot)**2
assert proxy_bound>=0
eps=F(old['proxy_error_upper'])
true_bound=(rootup(proxy_bound)+eps*invroot)**2
K=F(old['K'][0]);status='PASS_WITNESS' if true_bound<K else 'UNDECIDED'
out={'theta':str(theta),'N':N,'dmin':str(dmin),'dtail':str(dtail),'Q_residual_norm2':list(map(str,Q)),'finite_weighted_upper':str(weighted),'tail_norm2_upper':str(tail[1]),'weighted_residual_upper':str(A),'proxy_inverse_upper':str(proxy_bound),'true_inverse_upper':str(true_bound),'K_lower':str(K),'margin':str(K-true_bound),'status':status}
print('theta',theta,'inverse upper',float(true_bound),'K',float(K),status)
(p/'response_result.json').write_text(json.dumps(out,indent=2)+'\n')
