"""Independent arithmetic audit of the DI comparison bank, stdlib only."""
from pathlib import Path
from fractions import Fraction as F
import json
root=Path("custody/di/generated")
out=Path("custody/dj/generated")
out.mkdir(parents=True,exist_ok=True)
def mat(r,k):
    a=[[F(x) for x in row] for row in r[k]]
    n=r["n"]; assert len(a)==n and all(len(row)==n for row in a)
    assert all(a[i][j]==a[j][i] for i in range(n) for j in range(n))
    return a
def integer_positive(a):
    a=[row[:] for row in a];prev=1;ds=[]
    for k in range(len(a)):
        d=a[k][k];assert d>0
        ds.append(d)
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                v,rem=divmod(d*a[i][j]-a[i][k]*a[k][j],prev)
                assert rem==0
                a[i][j]=v
        prev=d
    return ds
S=10**12
checks=[]
result=json.loads((root/"result.json").read_text())
for dim in (23,):
    diagonal=None
    before=len(checks)
    summary=result[str(dim)]["best"]
    assert summary["positive"]
    for path in sorted((root/f"trial{dim}").glob("cj-mode*.json")):
        r=json.loads(path.read_text());n=r["n"];assert n==dim
        if not r["positive"]:
            assert r["mode"]==0 and F(r["lower"])<=0 and "margin" not in r and "upper_triangle" not in r
            continue
        assert r["mode"]==1
        th,delta,t,eta=map(F,r["parameters"]);ws=list(map(F,r["weights"]))
        assert min(th,delta,t,eta,*ws)>0 and sum(ws)==1
        margin,qdefect,lower=map(F,[r["margin"],r["qdefect"],r["lower"]])
        assert 0<margin<lower and 0<=qdefect<1
        M,H,R,P,G,V=(mat(r,k) for k in ("matrix","H","R","P","Gu_upper","F_lower"))
        D=[[H[i][j]-R[i][j]-V[i][j] for j in range(n)] for i in range(n)]
        if diagonal is None:diagonal=D
        assert diagonal==D
        for i in range(n):
            for j in range(n):
                expected=H[i][j]-D[i][j] if r["mode"]==0 else (1+t)*G[i][j]+(1+1/t)*R[i][j]
                assert P[i][j]==expected
        shift=margin*(1+qdefect)
        T=[[(M[i][j]-(shift if i==j else 0))*S//1-(n if i==j else 0) for j in range(n)] for i in range(n)]
        assert [T[i][j] for i in range(n) for j in range(i,n)]==r["upper_triangle"]
        minors=integer_positive(T)
        # The symmetric rounding error is bounded in spectral norm by n/S.
        # Therefore positive T implies M > margin*(1+qdefect)*I.
        pivots=[float(F(d, (minors[i-1] if i else 1)*S)) for i,d in enumerate(minors)]
        checks.append(dict(dimension=dim,file=path.name,mode=r["mode"],margin=str(margin),lower=str(lower),qdefect=str(qdefect),leading_minors=list(map(str,minors)),pivots=pivots))
        print("DJ EXACT PASS",dim,path.name,margin,min(pivots),flush=True)
    bank=checks[before:]
    assert len(bank)==3
    assert max(F(x["lower"]) for x in bank)==F(summary["lower"])
    assert any(x["margin"]==summary["margin"] and x["lower"]==summary["lower"] for x in bank)
negative=[]
for dim in (23,24,25):
    if dim!=23: assert not result[str(dim)]["best"]["positive"]
    bank=[]
    for path in sorted((root/f"trial{dim}").glob("cj-mode*.json")):
        rec=json.loads(path.read_text())
        assert rec["n"]==dim
        if rec["positive"]:
            assert dim==23 and rec["mode"]==1
            continue
        assert F(rec["lower"])<=0
        assert "margin" not in rec and "upper_triangle" not in rec
        bank.append(dict(dimension=dim,file=path.name,mode=rec["mode"],lower=rec["lower"]))
    assert len(bank)==(3 if dim==23 else 6)
    if dim!=23: assert max(F(x["lower"]) for x in bank)==F(result[str(dim)]["best"]["lower"])
    else: assert all(x["mode"]==0 for x in bank)
    negative.extend(bank)
assert len(negative)==15
assert len(checks)==3 and sum(len(x["leading_minors"]) for x in checks)==69
(out/"audit.json").write_text(json.dumps(dict(checks=checks,nonpassing_controls=negative,scope="terminal rational matrix to integer positive certificate; inherited analytical enclosure not regenerated"),indent=2)+"\n")
print("DJ AUDIT PASS",len(checks),sum(len(x["leading_minors"]) for x in checks),flush=True)
