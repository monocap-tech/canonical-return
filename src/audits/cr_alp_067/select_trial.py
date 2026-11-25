"""Untrusted floating trial selection; acceptance does not rely on this solve."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,numpy as np
p=Path(sys.argv[1]); bands=json.loads(Path(sys.argv[2]).read_text())['bands']
theta=F(sys.argv[3]); assert theta in (F(0),F(1,10))
N=2048; ids=np.arange(86,N+1,2); a=1.1; v=ids*np.pi/(2*a)
def parse(path):
    out={}
    for line in path.read_text().splitlines():
        x=line.split()
        if len(x)==5 and x[-2:]==['/','1e15']:
            assert x[0] not in out
            out[x[0]]=tuple(F(int(t),10**15) for t in x[1:3])
    return out
rows=parse(p/'coeff.log')
assert set(rows)=={'G',*(f'r{j}' for j in range(2,N+1,2))}
r=np.array([float(sum(rows[f'r{j}'])/2) for j in ids])
def lower(j):
    b=max((b for b in bands if b['j']<=j),key=lambda b:b['j'])
    d=(1-theta)*F(72,25)+theta*F(b['gamma_diagonal_lower'])
    return (d*10**9//1)
di=[lower(int(j)) for j in ids]
assert min(di)>F(21,8)*10**9
def dec(n):
    return f'{n//10**9}.{n%10**9:09d}'
(p/'diag.txt').write_text(''.join(f'{j} {dec(x)} {dec(x)}\n' for j,x in zip(ids,di)))
sh=np.log([2,3,4,5,7,8,9]);wt=np.log([2,3,2,5,7,2,3])/np.sqrt([2,3,4,5,7,8,9])
S=np.sin(v[:,None]*sh); C=np.cos(v[:,None]*sh); T=S@wt
num=v[None,:]*T[:,None]-v[:,None]*T[None,:]
den=v[None,:]**2-v[:,None]**2
np.fill_diagonal(den,1)
P=-2/a*num/den
np.fill_diagonal(P,-((2*a-sh)/a*C+S/(a*v[:,None]))@wt)
A=P+np.diag(np.array(di,dtype=float)/1e9)
# Common factor in exact pole coefficients cancels in finite constraint.
c=v/(v*v+.25)
z=np.linalg.solve(A,np.column_stack([r,c]))
x=z[:,0]-z[:,1]*(c@z[:,0])/(c@z[:,1])
assert np.all(np.isfinite(x)) and np.linalg.norm(x)<10
scale=10**12; X=[int(round(float(y)*scale)) for y in x]
def signed(n):
    return ('-' if n<0 else '')+f'{abs(n)//scale}.{abs(n)%scale:012d}'
(p/'trial.txt').write_text(''.join(f'{j} {signed(x)} {signed(x)}\n' for j,x in zip(ids,X)))
print('NUMERICAL ONLY theta',theta,'finite inverse response',r@x,'trial norm',np.linalg.norm(x),flush=True)
