"""Independent RETURN088 odd augmented tail Gram.
No parent mathematics imported during reconstruction. Parent tail routine is
executed separately, only as a containment reference, after reconstruction.
"""
import ast, json, math, sys, hashlib
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np

root=Path(sys.argv[1]).resolve()
out=Path(sys.argv[2]).resolve(); out.mkdir(parents=True,exist_ok=True)
S=10**50
A=F(87,500)
N=6
def load(p): return json.loads((root/p).read_bytes())
def pt(x): return (F(x),F(x))
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b):
    v=[x*y for x in a for y in b]; return min(v),max(v)
def isum(xs):
    l=h=F(0)
    for a,b in xs: l+=a;h+=b
    return l,h
def mag(a): return max(map(abs,a))
def enc(a): return math.floor(a[0]*S),math.ceil(a[1]*S)
def arr(a): return np.array([[enc(x) for x in row] for row in a],dtype=object)
def mm(a,b):
    # Exact endpoints, independent of parent midpoint/radius implementation.
    lo=np.zeros((len(a),b.shape[1]),object); hi=lo.copy()
    for k in range(a.shape[1]):
        v=[a[:,k,i,None]*b[k,None,:,j] for i in (0,1) for j in (0,1)]
        lo+=np.minimum(np.minimum(v[0],v[1]),np.minimum(v[2],v[3]))
        hi+=np.maximum(np.maximum(v[0],v[1]),np.maximum(v[2],v[3]))
    return np.stack([lo//S,-((-hi)//S)],-1)

cr=load('return063/RETURN063_constraints_receipt.json')['blocks'][1]
inv=load('return063/RETURN063_inverse_receipt.json')['blocks'][1]
bank=load('return062/RETURN062_bank_receipt.json')
pr=load('return077/proposal.json')
cs=cr['nodes']; Z=[[F(z,cr['Zscale']) for z in row] for row in cr['Z']]
js=pr['js']; BS=pr['scale']; B0=pr['rows']
assert len(B0)==8 and len(js)==500 and all(sum(row)==0 for row in B0)
c0,c1=map(F,bank['c']); assert 0<c0<=c1
def node(j): return (-c1,-c0) if j=='sinh' else pt(j*j)

@lru_cache(None)
def moment(K,p):
    # Sum over k=K+2,K+4,... of (K/k)^p?  Inverse power (K/k)^p.
    # Convex trapezoid correction: 0<=error<=h/8*integral f''=p/(4K).
    assert K%2==0 and p>1
    mid=F(K,2*(p-1))-F(1,2)
    return mid,mid+F(p,4*K)

def coefficients(K):
    prefix=[]; er=[]
    for l in range(len(cs)):
        v=[pt(1)]+[pt(0)]*(N-1); major=[F(1)]+[F(0)]*(N-1)
        total=F(1)
        for j in cs[:l+1]:
            a=mul(node(j),pt(F(1,K*K))); r=mag(a); assert r<1
            # Geometric polynomial coefficients via recurrence, exact rational.
            w=[]; z=pt(1)
            for n in range(N): w.append(z); z=mul(z,a)
            v=[isum(mul(v[i],w[n-i]) for i in range(n+1)) for n in range(N)]
            major=[sum(major[i]*r**(n-i) for i in range(n+1)) for n in range(N)]
            total/=1-r
        f=F(116**(2*l),K**(2*l))
        prefix.append([pt(0) if n<l else mul(pt(f),v[n-l]) for n in range(N)])
        er.append(f*(total-sum(major[:max(0,N-l)])))
    C=[[isum(mul(pt(z),v[n]) for z,v in zip(row,prefix)) for n in range(N)] for row in Z]
    E=[sum(abs(z)*e for z,e in zip(row,er)) for row in Z]
    for j in inv['nodes'][len(cs):]:
        a=F(inv['scales'][str(j)]); r=F(j*j,K*K); assert r<1
        C.append([pt(a*r**n) for n in range(N)]); E.append(abs(a)*r**N/(1-r))
    for row in B0:
        a=[F(z,BS) for z in row]; r=[F(j*j,K*K) for j in js]; assert max(r)<1
        C.append([pt(sum(z*t**n for z,t in zip(a,r))) for n in range(N)])
        E.append(sum(abs(z)*t**N/(1-t) for z,t in zip(a,r)))
        assert C[-1][0]==pt(0)
    assert len(C)==len(E)==527 and min(E)>=0
    return C,E

def reconstruct(K):
    C,E=coefficients(K)
    D=[[mul(pt(1/(A*K*K)),moment(K,2+2*(i+j))) for j in range(N)] for i in range(N)]
    X=arr(C); G=mm(mm(X,arr(D)),X.transpose(1,0,2))
    # Omitted row polynomial is bounded by E_i*(K/k)^(2N)/k.
    # Pay each product term with its actual decay, rather than the parent
    # uniform k^-2 remainder estimate.
    erow=np.array([math.ceil(e*S) for e in E],object)
    M=np.array([[math.ceil(mag(c)*S) for c in row] for row in C],object)
    w=np.array([math.ceil(moment(K,2+2*(N+i))[1]*S/(A*K*K)) for i in range(N)],object)
    e2=math.ceil(moment(K,2+4*N)[1]*S/(A*K*K))
    v=M@w
    error=v[:,None]*erow[None,:]+erow[:,None]*v[None,:]+e2*erow[:,None]*erow[None,:]
    error=-((-error)//(S*S))
    G[:,:,0]-=error; G[:,:,1]+=error
    G[:,:,0]=np.maximum(G[:,:,0],G[:,:,0].T)
    G[:,:,1]=np.minimum(G[:,:,1],G[:,:,1].T)
    assert np.all(G[:,:,0]<=G[:,:,1])
    assert all(G[i,i,0]>=0 for i in range(527))
    return G

def compare(g,p,scale):
    lo=g[:,:,0]*scale; hi=g[:,:,1]*scale
    l=p[:,:,0]*S; h=p[:,:,1]*S
    inside=(lo>=l)&(hi<=h); overlap=(lo<=h)&(hi>=l)
    return dict(entries=int(inside.size),contained=int(inside.sum()),overlap=int(overlap.sum()),disjoint=int((~overlap).sum()),max_width_units=str(max((g[:,:,1]-g[:,:,0]).flat)))

# Independent computation first, before importing the reference.
banks={K:reconstruct(K) for K in (65536,262144)}
for K,G in banks.items():
    (out/f'gram_{K}.json').write_text(json.dumps({'scale':S,'K':K,'gram':G.tolist()},separators=(',',':')))
    print('independent augmented Gram',K,G.shape,flush=True)

# Reference: compile only the archived parent function; not its terminal solve.
sys.path.insert(0,str(root/'return073'))
import return073_band as b
src=root/'return083/upgrade.py'
fn=next(x for x in ast.parse(src.read_text()).body if isinstance(x,ast.FunctionDef) and x.name=='make_tailgram')
env=dict(b=b,F=F,math=math,np=np,A0=A,cs=cs,Z=Z,ib=inv,N=4,B0=B0,js=js,BS=BS,P=b.P,mul=b.mul,add=b.add,div=b.div,isum=b.isum)
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(src),'exec'),env)
receipt={}
for alt in (False,True):
    env['alt']=alt
    for K,G in banks.items():
        r=compare(G,env['make_tailgram'](K),b.J)
        receipt[f'{K}_{"alternate" if alt else "generator"}']=r
        print(K,alt,r,flush=True)
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
assert all(x['contained']==x['entries'] for x in receipt.values()),receipt
print('RETURN088 AUGMENTED TAIL GRAM CONTAINMENT PASS',flush=True)
