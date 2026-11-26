"""RETURN090-A: independent restricted-band contractions; shared RETURN089 source primitives and D0. Not an asymptotic or terminal certification."""
import json, math, sys
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
root=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=True)
K=32768; A=F(87,500); Q=10**55; S=10**50; N=8; O=10**12
def load(p):return json.loads((root/p).read_bytes())
def P(x):return F(x),F(x)
def mul(a,b):
    t=[x*y for x in a for y in b];return min(t),max(t)
def isum(xs):
    l=h=F(0)
    for a,b in xs:l+=a;h+=b
    return l,h
def neg(a):return -a[1],-a[0]
def mag(a):return max(map(abs,a))
def widen(a,e):return a[0]-e,a[1]+e
def enc(a):return math.floor(a[0]*S),math.ceil(a[1]*S)
def quant(a):return tuple(F(x,S) for x in enc(a))
def arr(a):return np.array([[enc(x) for x in row] for row in a],object)
def mm(a,b):
    lo=np.zeros((len(a),b.shape[1]),object);hi=lo.copy()
    for k in range(a.shape[1]):
        v=[a[:,k,i,None]*b[k,None,:,j] for i in (0,1) for j in (0,1)]
        lo+=np.minimum(np.minimum(v[0],v[1]),np.minimum(v[2],v[3]))
        hi+=np.maximum(np.maximum(v[0],v[1]),np.maximum(v[2],v[3]))
    return np.stack([lo//S,-((-hi)//S)],-1)

mom=load('return066/RETURN066_moments.json'); raw=load('return059/RETURN059_raw_receipt.json')
cr=load('return063/RETURN063_constraints_receipt.json'); inv=load('return063/RETURN063_inverse_receipt.json')
c0,c1=map(F,load('return062/RETURN062_bank_receipt.json')['c'])
pr=load('return077/proposal.json'); js=pr['js']; BS=pr['scale']; rows=pr['rows']
assert len(rows)==8 and all(sum(r)==0 for r in rows)
direct=load('return086/direct_sums.json')
bank=next(x for x in direct['bands'] if x['name']=='restricted')
assert bank['count']==98304 and bank['p']==list(range(2,47))
moments={p:[tuple(F(y,bank['scale']) for y in a) for a in r] for p,r in zip(bank['p'],bank['sums'])}
# Only D0 is used from parent initialization, not its scalar/contraction methods.
sys.path.insert(0,str(root/'return073'))
import return073_band as ref
D0=ref.D0
pi=(F(314159265358979323846,10**20),F(314159265358979323847,10**20))
ri=(F(5,2*K)/pi[1],F(5,2*K)/pi[0])
col=[[0,1,3,4],[1,2,5,6],[3,5,7,8],[4,6,8,9]]
def feature(p,s,c):return p,s,tuple(c)
def features(par):
    out=[]
    for n in range(N):
        a=feature(2*n+1,ri,[mul(pi,P(F(1,2))),P(0),P(1),P(0)])
        b=feature(2*n+2,mul(ri,ri),[D0,P(1),P(0),P(-1)])
        c=feature(2*n+2,mul(ri,ri),[P(2),P(0),P(0),P(0)])
        if par:
            a=feature(2*n+1,ri,[neg(D0),P(-1),P(0),P(1)])
            b=feature(2*n+2,mul(ri,ri),[mul(pi,P(F(1,2))),P(0),P(1),P(0)])
            c=feature(2*n+1,ri,[P(2),P(0),P(0),P(0)])
        out.extend([a,b,c])
    return out

@lru_cache(None)
def linear(p,c):
    return isum(mul(x,moments[p][i]) for x,i in zip(c,[0,1,3,4]) if x!=P(0))
def gram(f,g):
    p,s,c=f;q,t,d=g
    v=isum(mul(mul(c[i],d[j]),moments[p+q][col[i][j]]) for i in range(4) for j in range(4) if c[i]!=P(0) and d[j]!=P(0))
    return mul(P(1/A),mul(mul(s,t),v))
def node(j):return (-c1,-c0) if j=='sinh' else P(j*j)
@lru_cache(None)
def coeff(ts,n):
    v=[P(1)]+[P(0)]*(n-1); m=[F(1)]+[F(0)]*(n-1); total=F(1)
    for j in ts:
        a=mul(node(j),P(F(1,K*K)));r=mag(a);assert r<1
        w=[P(1)]
        for i in range(1,n):w.append(quant(mul(w[-1],a)))
        v=[quant(isum(mul(v[i],w[t-i]) for i in range(t+1))) for t in range(n)]
        m=[sum(m[i]*r**(t-i) for i in range(t+1)) for t in range(n)]
        total/=1-r
    return v,total-sum(m)

def cross(f,ts,pref,deg):
    p,s,c=f;p+=deg
    # Keep the extra moment required by the remainder below p=46.
    n=min(12,(46-p)//2);assert n>=8
    v,tail=coeff(tuple(ts),n)
    z=isum(mul(a,linear(p+2*i,c)) for i,a in enumerate(v))
    bound=mag(c[0])+3*mag(c[1])+9*(mag(c[2])+mag(c[3]))
    e=tail*bound*moments[p+2*n][0][1]
    return quant(mul(P(pref/A),mul(s,widen(z,e))))

def restore(f,row):
    p,s,c=f;p+=1;n=min(12,(46-p)//2);assert n>=8
    a=[F(z,BS) for z in row];r=[F(j*j,K*K) for j in js];assert max(r)<1
    v=[sum(z*t**i for z,t in zip(a,r)) for i in range(n)]
    assert v[0]==0
    e=sum(abs(z)*t**n/(1-t) for z,t in zip(a,r))
    bound=mag(c[0])+3*mag(c[1])+9*(mag(c[2])+mag(c[3]))
    z=isum(mul(P(v[i]),linear(p+2*i,c)) for i in range(1,n))
    return quant(mul(P(F(1,K)/A),mul(s,widen(z,e*bound*moments[p+2*n][0][1]))))

def compare(g,p,scale):
    l=g[:,:,0]*scale;h=g[:,:,1]*scale
    pl=p[:,:,0]*S;ph=p[:,:,1]*S
    inside=(l>=pl)&(h<=ph);over=(l<=ph)&(h>=pl)
    return dict(entries=int(inside.size),contained=int(inside.sum()),overlap=int(over.sum()),disjoint=int((~over).sum()),max_width_units=str(max((g[:,:,1]-g[:,:,0]).flat)))


# Reconstruct the odd restricted band at the ORIGINAL anchor 32768.
# No use of the archived restricted-band phase or matrix routine here.
par=1
cb=cr['blocks'][par]; ib=inv['blocks'][par]; nodes=cb['nodes']
Z=[[F(z,cb['Zscale']) for z in row] for row in cb['Z']]
m=mom['blocks'][par]['cutoffs'][0]; assert m['K']==K
VV=[]
for n in range(N):
    for a in m['raw_moments'][n]:
        VV.append([tuple(F(sum(fn(x*z,y*z) for (x,y),z in zip(a,col)),Q*raw['T']) for fn in (min,max)) for col in zip(*raw['blocks'][par]['C'])])
V=arr(VV); fs=features(par)
M=arr([[gram(f,g) for g in fs] for f in fs])
U=[]
for f in fs:
    u=[cross(f,nodes[:i+1],F(116**(2*i),K**(2*i+1)),2*i+1) for i in range(len(nodes))]
    v=[quant(isum(mul(P(z),a) for z,a in zip(row,u))) for row in Z]
    v += [cross(f,[j],F(ib['scales'][str(j)])/K,1) for j in ib['nodes'][len(nodes):]]
    U.append(v)
mats={'tt':mm(mm(V.transpose(1,0,2),M),V),
      'U':mm(arr(list(map(list,zip(*U)))),V),
      'restoration_cross':mm(arr([[restore(f,row) for f in fs] for row in rows]),V)}
assert mats['tt'].shape==(87,87,2)
assert mats['U'].shape==(519,87,2)
assert mats['restoration_cross'].shape==(8,87,2)
for a in mats.values(): assert np.all(a[:,:,0]<=a[:,:,1])
theta=1-A/F(519,1000)
def reweight(a):
    # Directed interval for the SUBTRACTED contribution -theta * B.
    return np.stack([(-a[:,:,1]*theta.numerator)//theta.denominator,
                     -((a[:,:,0]*theta.numerator)//theta.denominator)],axis=-1)
(out/'restricted_bank.json').write_text(json.dumps({'scale':S,'anchor':K,
    'band':[65536,262144],'theta':str(theta),
    'matrices':{k:v.tolist() for k,v in mats.items()},
    'signed_reweighting':{k:reweight(v).tolist() for k,v in mats.items()}},
    separators=(',',':')))
# Only now evaluate the frozen parent restricted-band implementation, as
# a reference. Do not execute its solver or terminal comparison.
src=(root/'return083/upgrade.py').read_text()
start=src.index('from functools import lru_cache',src.index('# Restrict the same anchored'))
end=src.index('def subband(',start)
snippet=src[start:end]
receipt={}
for alt in (False,True):
    ref.alt=alt
    env=dict(b=ref,F=F,math=math,np=np,kb=K,K2=65536,K0=262144,
             nb=8,fb=ref.features(1),VB=ref.iarr(VV),cs=nodes,Z=Z,ib=ib,
             B0=rows,js=js,BS=BS,A0=A,P=ref.P,mul=ref.mul,add=ref.add,
             isum=ref.isum,alt=alt)
    exec(compile(snippet,'archived_restricted_reference','exec'),env)
    for k,key in [('tt','bandtt'),('U','bandu'),('restoration_cross','bandh')]:
        r=compare(mats[k],env[key],ref.J)
        receipt[f'{alt}_{k}']=r
        print(alt,k,r,flush=True)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
assert all(x['contained']==x['entries'] for x in receipt.values()),receipt
print('RETURN090 RESTRICTED BAND CONTAINMENT PASS',flush=True)
