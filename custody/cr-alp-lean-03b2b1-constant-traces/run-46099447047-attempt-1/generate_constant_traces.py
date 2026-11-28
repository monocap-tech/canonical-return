#!/usr/bin/env python3
"""Actual 046 fixed elementary constants. Exact term enclosures, not rounded substitutes."""
import ast, hashlib, json, sys
from pathlib import Path
from fractions import Fraction as F
from math import factorial
source=Path(sys.argv[1]); out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
raw=source.read_bytes();blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
assert blob=="554e7ab54da4545987d89aa660a4abe6c951cb63"
code=raw.decode(); marker="PI=16*atan(I(F(1,5)))-4*atan(I(F(1,239)))"
assert code.count(marker)==1
env={"__file__":str(source.resolve()),"__name__":"alp046_constants"}
exec(compile(code.split(marker)[0],str(source),"exec"),env)
node=next(n for n in ast.parse(code).body if isinstance(n,ast.FunctionDef) and n.name=="exp_pos")
exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),"exec"),env)
I,S=env["I"],env["S"];assert S==10**340
cases=[("Atan5","atan",1,5),("Atan239","atan",1,239),("ExpA","exp",11,10),("ExpHalfA","exp",11,20)]
traces={}
for key,kind,p,q in cases:
    x=I(F(p,q));rows=[];z=I(0)
    if kind=="atan":
        assert x.up()<=F(1,3)
        v=-x.sq();t=x
        for i in range(400):
            term=t/(2*i+1)
            num=(-1)**i*p**(2*i+1);den=q**(2*i+1)*(2*i+1)
            assert den>0 and term.l*den<=num*S<=term.h*den
            rows.append([term.l,term.h]);z+=term;t=t*v
        original=env["atan"](x)
    else:
        assert 0<=x.l and x.up()<3
        t=I(1)
        for i in range(400):
            if i:t=t*x/i
            num=p**i;den=q**i*factorial(i)
            assert den>0 and t.l*den<=num*S<=t.h*den
            rows.append([t.l,t.h]);z+=t
        original=env["exp_pos"](x)
    a=z.widen(1)
    assert (a.l,a.h)==(original.l,original.h)
    assert a.l==sum(l for l,h in rows)-1 and a.h==sum(h for l,h in rows)+1
    def arr(k):return ",\n  ".join(str(row[k]) for row in rows)
    fname="Const"+key
    func="arctan" if kind=="atan" else "exp"
    terminal=f"{kind}Trace {p} {q} lo hi L U (by norm_num) allTerms lowerSum upperSum"
    if kind=="atan":
        terminal="by\n  simpa only [Nat.cast_one, Nat.cast_ofNat] using ("+terminal+")"
    lean=f"""import RationalTrace
namespace CRALP.{fname}
open Finset
open scoped BigOperators
open CRALP.RationalTrace CRALP.DirectedEvaluation CRALP.DirectedInterval
set_option maxRecDepth 200000
set_option maxHeartbeats 0
def lo (i : ℕ) : ℤ := (#[{arr(0)}] : Array ℤ)[i]!
def hi (i : ℕ) : ℤ := (#[{arr(1)}] : Array ℤ)[i]!
def L : ℤ := {a.l}
def U : ℤ := {a.h}
theorem allTerms : ∀ i : Fin 400, row ({kind}Num {p}) ({kind}Den {q}) lo hi i := by decide
theorem lowerSum : L = (∑ i ∈ range 400, lo i)-1 := by decide
theorem upperSum : U = (∑ i ∈ range 400, hi i)+1 := by decide
theorem actual : Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) (Real.{func} ({p}/{q})) :=
  {terminal}
end CRALP.{fname}
"""+''.join(f"#print axioms CRALP.{fname}.{n}\n" for n in ["allTerms","lowerSum","upperSum","actual"])
    (out/(fname+".lean")).write_text(lean)
    traces[key]={"kind":kind,"p":p,"q":q,"terms":400,"rows":[list(map(str,row)) for row in rows],"endpoints":[str(a.l),str(a.h)],"width":str(a.h-a.l),"sha256":hashlib.sha256(lean.encode()).hexdigest()}
    print(key,400,"width",a.h-a.l)
a=env["atan"](I(F(1,5)));b=env["atan"](I(F(1,239)));pi=16*a-4*b
assert pi.l>3*S and pi.h*7<22*S
lean=f"""import ConstAtan5
import ConstAtan239
namespace CRALP.ConstPi
open CRALP.DirectedEvaluation CRALP.DirectedInterval
def L : ℤ := {pi.l}
def U : ℤ := {pi.h}
theorem lowerEq : L=16*CRALP.ConstAtan5.L-4*CRALP.ConstAtan239.U := by decide
theorem upperEq : U=16*CRALP.ConstAtan5.U-4*CRALP.ConstAtan239.L := by decide
theorem actual : Contains ((L:ℝ)/logScale) ((U:ℝ)/logScale) Real.pi := by
  have h := machinAttach CRALP.ConstAtan5.actual CRALP.ConstAtan239.actual
  have hL : (L:ℝ)/logScale =
      16*((CRALP.ConstAtan5.L:ℝ)/logScale)-4*((CRALP.ConstAtan239.U:ℝ)/logScale) := by
    rw [lowerEq]
    push_cast <;> ring
  have hU : (U:ℝ)/logScale =
      16*((CRALP.ConstAtan5.U:ℝ)/logScale)-4*((CRALP.ConstAtan239.L:ℝ)/logScale) := by
    rw [upperEq]
    push_cast <;> ring
  rw [hL,hU]
  exact h
end CRALP.ConstPi
#print axioms CRALP.ConstPi.lowerEq
#print axioms CRALP.ConstPi.upperEq
#print axioms CRALP.ConstPi.actual
"""
(out/"ConstPi.lean").write_text(lean)
traces["Pi"]={"endpoints":[str(pi.l),str(pi.h)],"width":str(pi.h-pi.l),"sha256":hashlib.sha256(lean.encode()).hexdigest()}
(out/"constant-traces.json").write_text(json.dumps({"source_blob":blob,"scale":str(S),"cases":traces},indent=2)+"\n")
