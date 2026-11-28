#!/usr/bin/env python3
"""Replay actual 046 LOG2 and emit an exact termwise interval certificate."""
import hashlib, json, sys
from pathlib import Path
from fractions import Fraction

source=Path(sys.argv[1])
out=Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
raw=source.read_bytes()
blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
assert blob=="554e7ab54da4545987d89aa660a4abe6c951cb63"
marker="LOG2=log_series(I(2))"
code=raw.decode()
assert code.count(marker)==1
env={"__file__":str(source.resolve()),"__name__":"alp046_log2_replay"}
exec(compile(code.split(marker)[0],str(source),"exec"),env)
I,S=env["I"],env["S"]
assert S==10**340
x=I(2); u=(x-1)/(x+1); v=u.sq(); t=u; z=I(0)
assert u.up()<Fraction(7,20)
assert 3*u.l<=S<=3*u.h
rows=[]
for i in range(400):
    term=t/(2*i+1)
    den=3**(2*i+1)*(2*i+1)
    assert term.l*den<=S<=term.h*den
    rows.append([term.l,term.h])
    z+=term
    t=t*v
actual=(2*z).widen(1)
original=env["log_series"](I(2))
assert (actual.l,actual.h)==(original.l,original.h)
assert actual.l==2*sum(row[0] for row in rows)-1
assert actual.h==2*sum(row[1] for row in rows)+1
assert actual.l<actual.h
def arr(k):
    return ",\n  ".join(str(row[k]) for row in rows)
lean=f"""import DirectedEvaluation

namespace CRALP.Log2Actual
open Finset
open scoped BigOperators
open CRALP.DirectedEvaluation
open CRALP.DirectedInterval
set_option maxRecDepth 200000
set_option maxHeartbeats 0

def lo (i : ℕ) : ℤ := (#[{arr(0)}] : Array ℤ)[i]!
def hi (i : ℕ) : ℤ := (#[{arr(1)}] : Array ℤ)[i]!
def L : ℤ := {actual.l}
def U : ℤ := {actual.h}

theorem allTerms : ∀ i : Fin 400, logRow lo hi i := by decide
theorem lowerSum : L = 2*(∑ i ∈ range 400, lo i)-1 := by decide
theorem upperSum : U = 2*(∑ i ∈ range 400, hi i)+1 := by decide
theorem actualLog2 :
    Contains ((L : ℝ)/(logScale:ℝ)) ((U : ℝ)/(logScale:ℝ)) (Real.log 2) :=
  log2Trace lo hi L U allTerms lowerSum upperSum

end CRALP.Log2Actual
#print axioms CRALP.Log2Actual.allTerms
#print axioms CRALP.Log2Actual.lowerSum
#print axioms CRALP.Log2Actual.upperSum
#print axioms CRALP.Log2Actual.actualLog2
"""
(out/"Log2Actual.lean").write_text(lean)
trace={"source_blob":blob,"scale":str(S),"terms":400,"term_bounds":[list(map(str,row)) for row in rows],"reduced_argument":[str(u.l),str(u.h)],"actual_endpoints":[str(actual.l),str(actual.h)],"width_units":str(actual.h-actual.l),"generated_sha256":hashlib.sha256(lean.encode()).hexdigest(),"scope":"Original 046 LOG2 only; each term checked against exact rational at u=1/3, total endpoint equality, analytic remainder in Lean."}
(out/"log2-trace.json").write_text(json.dumps(trace,indent=2)+"\n")
print("Actual 046 LOG2 replay:",len(rows),"terms, width units",actual.h-actual.l)
