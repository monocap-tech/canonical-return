#!/usr/bin/env python3
"""Original 046 prime-window logarithms, deduplicated reduced rational traces."""
import hashlib,json,sys
from pathlib import Path
from fractions import Fraction as F
src=Path(sys.argv[1]);out=Path(sys.argv[2]);old=Path(sys.argv[3])
out.mkdir(parents=True,exist_ok=True)
raw=src.read_bytes();blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
assert blob=="554e7ab54da4545987d89aa660a4abe6c951cb63"
code=raw.decode();marker="def atan(x):";assert code.count(marker)==1
env={"__file__":str(src.resolve()),"__name__":"alp046_prime_logs"}
exec(compile(code.split(marker)[0],str(src),"exec"),env)
I,S=env["I"],env["S"]
assert S==10**340
log2=env["LOG2"]
prev=json.loads(old.read_text())
assert prev["actual_endpoints"]==[str(log2.l),str(log2.h)]
ns=(2,3,4,5,7,8,9)
assert "NS=(2,3,4,5,7,8,9)" in code
reduced={};cases={};imports=["import Log2Actual"]
blocks=[]
for n in ns:
    x=I(n);y=F(n);k=0
    while x.l+x.h>4*S:x=x/2;y/=2;k+=1
    while x.l+x.h<2*S:x=x*2;y*=2;k-=1
    assert k>=0 and y==F(n,2**k) and x.l*y.denominator<=y.numerator*S<=x.h*y.denominator
    u=(x-1)/(x+1);v=u.sq();t=u;z=I(0)
    rat=(y-1)/(y+1);p,q=rat.numerator,rat.denominator
    assert 0<=p and q>0 and u.up()<F(7,20)
    assert u.l*q<=p*S<=u.h*q
    rows=[]
    for i in range(400):
        d=2*i+1;term=t/d;num=p**d;den=q**d*d
        assert term.l*den<=num*S<=term.h*den
        rows.append([str(term.l),str(term.h)])
        z+=term;t=t*v
    a=(2*z).widen(1)
    b=env["log_series"](x)
    assert (a.l,a.h)==(b.l,b.h)
    assert a.l==2*sum(int(r[0]) for r in rows)-1 and a.h==2*sum(int(r[1]) for r in rows)+1
    entry={"y":[y.numerator,y.denominator],"p":p,"q":q,"rows":rows,"endpoints":[str(a.l),str(a.h)],"width":str(a.h-a.l)}
    key=str(y)
    if y==2:
        assert entry["endpoints"]==prev["actual_endpoints"] and rows==prev["term_bounds"]
        mod="Log2Actual";actual="actualLog2"
    elif key not in reduced:
        mod="PrimeLog"+str(n);actual="actual"
        def arr(j):return ",\n  ".join(r[j] for r in rows)
        lean=f"""import PrimeLogTrace
namespace CRALP.{mod}
open Finset
open scoped BigOperators
open CRALP.RationalTrace CRALP.DirectedEvaluation CRALP.DirectedInterval CRALP.PrimeLogTrace
set_option maxRecDepth 200000
set_option maxHeartbeats 0
def lo (i : ℕ) : ℤ := (#[{arr(0)}] : Array ℤ)[i]!
def hi (i : ℕ) : ℤ := (#[{arr(1)}] : Array ℤ)[i]!
def L : ℤ := {a.l}
def U : ℤ := {a.h}
theorem allTerms : ∀ i : Fin 400, row (termNum {p}) (termDen {q}) lo hi i := by decide
theorem lowerSum : L = 2*(∑ i ∈ range 400, lo i)-1 := by decide
theorem upperSum : U = 2*(∑ i ∈ range 400, hi i)+1 := by decide
theorem actual : Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.log ({y.numerator}/{y.denominator})) :=
  reducedTrace ({y.numerator}/{y.denominator}) {p} {q} lo hi L U
    (by norm_num) (by norm_num) (by norm_num) allTerms lowerSum upperSum
end CRALP.{mod}
"""+''.join(f"#print axioms CRALP.{mod}.{j}\n" for j in ["allTerms","lowerSum","upperSum","actual"])
        (out/(mod+".lean")).write_text(lean)
        imports.append("import "+mod)
        entry["sha256"]=hashlib.sha256(lean.encode()).hexdigest()
    else:
        raise AssertionError("unexpected non-LOG2 duplicate")
    if key in reduced:
        assert reduced[key]==entry
    else:
        reduced[key]=entry
    full=a+k*log2
    original=env["log"](I(n))
    assert (full.l,full.h)==(original.l,original.h)
    assert full.l==a.l+k*log2.l and full.h==a.h+k*log2.h
    cases[str(n)]={"k":k,"reduced":key,"endpoints":[str(full.l),str(full.h)],"width":str(full.h-full.l)}
    yr="2" if y==2 else f"({y.numerator}/{y.denominator})"
    blocks.append(f"""namespace N{n}
def L : ℤ := {full.l}
def U : ℤ := {full.h}
theorem lowerEq : L=CRALP.{mod}.L+({k}:ℤ)*CRALP.Log2Actual.L := by decide
theorem upperEq : U=CRALP.{mod}.U+({k}:ℤ)*CRALP.Log2Actual.U := by decide
theorem actual : Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.log {n}) :=
  scaleTrace {n} {yr} {k} CRALP.{mod}.L CRALP.{mod}.U
    CRALP.Log2Actual.L CRALP.Log2Actual.U L U (by norm_num) (by norm_num)
    CRALP.{mod}.{actual} CRALP.Log2Actual.actualLog2 lowerEq upperEq
end N{n}
""")
    print(n,"k",k,"reduced",y,"width",full.h-full.l)
lean="\n".join(imports+["import PrimeLogTrace","namespace CRALP.PrimeLogs","open CRALP.PrimeLogTrace CRALP.DirectedEvaluation CRALP.DirectedInterval"]+blocks+["end CRALP.PrimeLogs"]+[f"#print axioms CRALP.PrimeLogs.N{n}.{j}" for n in ns for j in ["lowerEq","upperEq","actual"]])+"\n"
(out/"PrimeLogs.lean").write_text(lean)
assert len(reduced)==5 and len(cases)==7
(out/"prime-log-traces.json").write_text(json.dumps({"source_blob":blob,"scale":str(S),"reduced":reduced,"cases":cases,"assembly_sha256":hashlib.sha256(lean.encode()).hexdigest()},indent=2)+"\n")
