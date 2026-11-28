#!/usr/bin/env python3
"""Literal 046 log(PI) replay, with bounded-width recurrence witnesses."""
import hashlib,json,sys
from pathlib import Path
src,out,pi_data,log_data=map(Path,sys.argv[1:5])
out.mkdir(parents=True,exist_ok=True)
raw=src.read_bytes();blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
assert blob=="554e7ab54da4545987d89aa660a4abe6c951cb63"
code=raw.decode();marker="# Euler-Maclaurin with derivative"
assert code.count(marker)==1
env={"__file__":str(src.resolve()),"__name__":"alp046_log_pi"}
exec(compile(code.split(marker)[0],str(src),"exec"),env)
I,S=env["I"],env["S"];pi=env["PI"];log2=env["LOG2"]
assert S==10**340
def ends(a):return [str(a.l),str(a.h)]
assert ends(pi)==json.loads(pi_data.read_text())["cases"]["Pi"]["endpoints"]
assert ends(log2)==json.loads(log_data.read_text())["actual_endpoints"]
x=pi;k=0
while x.l+x.h>4*S:x=x/2;k+=1
while x.l+x.h<2*S:x=x*2;k-=1
assert k==1
u=(x-1)/(x+1);v=u.sq();assert u.l>0 and u.up()<env["F"](7,20)
t=u;z=I(0);powers=[t];terms=[]
def mulok(a,b,c):return all(c.l*S<=i*j<=c.h*S for i in (a.l,a.h) for j in (b.l,b.h))
def divok(a,b,c):return b.l>0 and all(c.l*j<=i*S<=c.h*j for i in (a.l,a.h) for j in (b.l,b.h))
assert divok(pi,I(2),x) and divok(x-1,x+1,u) and mulok(u,u,v)
for i in range(400):
    d=2*i+1;term=t/d;nt=t*v
    assert divok(t,I(d),term) and mulok(t,v,nt)
    terms.append(term);powers.append(nt);z+=term;t=nt
red=(2*z).widen(1);full=red+log2
assert ends(red)==ends(env["log_series"](x))
assert ends(full)==ends(env["log"](pi))
assert red.l==2*sum(a.l for a in terms)-1 and red.h==2*sum(a.h for a in terms)+1
def box(a):return f"⟨{a.l}, {a.h}⟩"
def arr(xs):return ",\n  ".join(box(a) for a in xs)
lean=f"""import GridTrace
import ConstPi
import Log2Actual
namespace CRALP.LogPiActual
open Finset
open scoped BigOperators
open CRALP.GridTrace CRALP.DirectedEvaluation CRALP.DirectedInterval
set_option maxRecDepth 200000
set_option maxHeartbeats 0
def piBox : Grid := ⟨CRALP.ConstPi.L,CRALP.ConstPi.U⟩
def log2Box : Grid := ⟨CRALP.Log2Actual.L,CRALP.Log2Actual.U⟩
def xBox : Grid := {box(x)}
def uBox : Grid := {box(u)}
def vBox : Grid := {box(v)}
def t (i : ℕ) : Grid := (#[{arr(powers)}] : Array Grid)[i]!
def term (i : ℕ) : Grid := (#[{arr(terms)}] : Array Grid)[i]!
def reduced : Grid := {box(red)}
def actualBox : Grid := {box(full)}
theorem halfCheck : divOK piBox (point 2) xBox := by decide
theorem mobiusCheck : divOK (minus xBox (point 1)) (plus xBox (point 1)) uBox := by decide
theorem squareCheck : mulOK uBox uBox vBox := by decide
theorem guardCheck : absGuard uBox 7 20 := by decide
theorem allSteps : ∀ i : Fin 400, logStep vBox t term i := by decide
theorem lowerSum : reduced.lo=2*(∑ i ∈ range 400, (term i).lo)-1 := by decide
theorem upperSum : reduced.hi=2*(∑ i ∈ range 400, (term i).hi)+1 := by decide
theorem assemblyCheck : actualBox=plus reduced log2Box := by decide
theorem halfSound : Encl xBox (Real.pi/2) := by
  have h := divSound (show Encl piBox Real.pi from CRALP.ConstPi.actual) (pointSound 2) halfCheck
  simpa only [Int.cast_ofNat] using h
theorem mobiusSound : Encl uBox ((Real.pi/2-1)/(Real.pi/2+1)) := by
  have h := divSound (minusSound halfSound (pointSound 1)) (plusSound halfSound (pointSound 1)) mobiusCheck
  simpa only [Int.cast_one] using h
theorem squareSound : Encl vBox (((Real.pi/2-1)/(Real.pi/2+1))^2) := by
  simpa only [pow_two] using mulSound mobiusSound mobiusSound squareCheck
theorem reducedSound : Encl reduced (Real.log (Real.pi/2)) := by
  have hr := absSound mobiusSound guardCheck
  have hr' : |(Real.pi/2-1)/(Real.pi/2+1)| ≤ 7/20 := by
    simpa only [Int.cast_ofNat] using hr
  exact logTrace vBox t term reduced (by positivity) hr'
    (show Encl (t 0) ((Real.pi/2-1)/(Real.pi/2+1)) from mobiusSound)
    squareSound allSteps lowerSum upperSum
theorem actual : Encl actualBox (Real.log Real.pi) := by
  rw [assemblyCheck,ElementarySeries.logHalf Real.pi Real.pi_pos]
  exact plusSound reducedSound (show Encl log2Box (Real.log 2) from CRALP.Log2Actual.actualLog2)
end CRALP.LogPiActual
"""+''.join(f"#print axioms CRALP.LogPiActual.{n}\n" for n in ["halfCheck","mobiusCheck","squareCheck","guardCheck","allSteps","lowerSum","upperSum","assemblyCheck","halfSound","mobiusSound","squareSound","reducedSound","actual"])
(out/"LogPiActual.lean").write_text(lean)
trace={"source_blob":blob,"scale":str(S),"pi":ends(pi),"log2":ends(log2),"k":k,"x":ends(x),"u":ends(u),"v":ends(v),"powers":[ends(a) for a in powers],"terms":[ends(a) for a in terms],"reduced":ends(red),"actual":ends(full),"width":str(full.h-full.l),"sha256":hashlib.sha256(lean.encode()).hexdigest()}
(out/"log-pi-trace.json").write_text(json.dumps(trace,indent=2)+"\n")
print("Literal log(PI) recurrence",len(terms),"terms, width",full.h-full.l)
