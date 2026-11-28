#!/usr/bin/env python3
"""Original 046 high-frequency sine/cosine and shifted arctangent pilot."""
import hashlib,json,sys
from pathlib import Path
src,out,pi_data,log_data=map(Path,sys.argv[1:5]);out.mkdir(parents=True,exist_ok=True)
raw=src.read_bytes();blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
assert blob=="554e7ab54da4545987d89aa660a4abe6c951cb63"
code=raw.decode();marker="# Euler-Maclaurin with derivative";assert code.count(marker)==1
env={"__file__":str(src.resolve()),"__name__":"alp046_signed_pilot"}
exec(compile(code.split(marker)[0],str(src),"exec"),env)
I,S,F=env["I"],env["S"],env["F"];pi=env["PI"];one=I(1)
def ends(a):return [str(a.l),str(a.h)]
assert ends(pi)==json.loads(pi_data.read_text())["cases"]["Pi"]["endpoints"]
lg=env["log"](I(9));assert ends(lg)==json.loads(log_data.read_text())["cases"]["9"]["endpoints"]
collar=I(2*env["a"]);nu=85*env["NUF"];x=nu*lg
period=round(x.mid()/env["TWO_PI"].mid());r=x-period*env["TWO_PI"];q=-r.sq()
y=nu/2;shift=I(F(2049,4));w=y/shift;v=-w.sq()
assert r.up()<4 and w.up()<=F(1,3)
assert not r.l<=0<=r.h and w.l>0
def mulok(a,b,c):return all(c.l*S<=i*j<=c.h*S for i in (a.l,a.h) for j in (b.l,b.h))
def divok(a,b,c):return b.l>0 and all(c.l*j<=i*S<=c.h*j for i in (a.l,a.h) for j in (b.l,b.h))
assert divok(pi,collar,env["NUF"]) and mulok(I(85),env["NUF"],nu)
assert mulok(nu,lg,x) and mulok(I(2*period),pi,period*env["TWO_PI"])
assert mulok(r,r,-q) and divok(nu,I(2),y) and divok(y,shift,w) and mulok(w,w,-v)
traces={}
for name,initial,e in [("Sin",r,1),("Cos",one,0),("Atan",w,None)]:
    t=initial;ts=[t];ms=[];z=I(0)
    if e is None:
        for i in range(400):
            term=t/(2*i+1);nt=t*v
            assert divok(t,I(2*i+1),term) and mulok(t,v,nt)
            ms.append(term);ts.append(nt);z+=term;t=nt
    else:
        z=t
        for i in range(399):
            m=t*q;nt=m/((2*i+e+1)*(2*i+e+2))
            assert mulok(t,q,m) and divok(m,I((2*i+e+1)*(2*i+e+2)),nt)
            ms.append(m);ts.append(nt);z+=nt;t=nt
    result=z.widen(1)
    original=env["atan"](w) if e is None else env["sincos"](x)[1-e]
    assert ends(result)==ends(original)
    traces[name]={"t":[ends(a) for a in ts],"mid":[ends(a) for a in ms],"out":ends(result),"width":str(result.h-result.l)}
def box(a):return "⟨"+", ".join(ends(a))+"⟩"
def arr(rows):return ",\n  ".join("⟨"+", ".join(a)+"⟩" for a in rows)
lean=f"""import SignedGridSeries
import ConstPi
import PrimeLogs
namespace CRALP.SignedPilot
open Finset
open scoped BigOperators
open CRALP.GridTrace CRALP.SignedGridSeries CRALP.DirectedEvaluation CRALP.DirectedInterval
set_option maxRecDepth 200000
set_option maxHeartbeats 0
def piBox : Grid := ⟨CRALP.ConstPi.L,CRALP.ConstPi.U⟩
def logBox : Grid := ⟨CRALP.PrimeLogs.N9.L,CRALP.PrimeLogs.N9.U⟩
def collar : Grid := {box(collar)}
def nu0 : Grid := {box(env["NUF"])}
def nu : Grid := {box(nu)}
def arg : Grid := {box(x)}
def periodic : Grid := {box(period*env["TWO_PI"])}
def reduced : Grid := {box(r)}
def q : Grid := {box(q)}
def halfNu : Grid := {box(y)}
def shift : Grid := {box(shift)}
def atanArg : Grid := {box(w)}
def v : Grid := {box(v)}
theorem collarSound : Encl collar (11/5) := by
  have h := pointSound 11
  have g := pointSound 5
  have hh := divSound h g (show divOK (point 11) (point 5) collar by decide)
  simpa only [Int.cast_ofNat] using hh
theorem shiftSound : Encl shift (2049/4) := by
  have h := divSound (pointSound 2049) (pointSound 4) (show divOK (point 2049) (point 4) shift by decide)
  simpa only [Int.cast_ofNat] using h
theorem nu0Sound : Encl nu0 (Real.pi/(11/5)) :=
  divSound (show Encl piBox Real.pi from CRALP.ConstPi.actual) collarSound (by decide)
theorem nuSound : Encl nu (85*(Real.pi/(11/5))) := by
  simpa only [Int.cast_ofNat] using mulSound (pointSound 85) nu0Sound (show mulOK (point 85) nu0 nu by decide)
theorem argSound : Encl arg (85*(Real.pi/(11/5))*Real.log 9) :=
  mulSound nuSound (show Encl logBox (Real.log 9) from CRALP.PrimeLogs.N9.actual) (by decide)
theorem periodicSound : Encl periodic (({period}:ℤ)*(2*Real.pi)) := by
  have h := mulSound (pointSound {2*period}) (show Encl piBox Real.pi from CRALP.ConstPi.actual) (show mulOK (point {2*period}) piBox periodic by decide)
  convert h using 1 <;> norm_num <;> ring
theorem reducedSound : Encl reduced (85*(Real.pi/(11/5))*Real.log 9-({period}:ℤ)*(2*Real.pi)) := by
  have he : reduced=minus arg periodic := by decide
  rw [he]
  exact minusSound argSound periodicSound
theorem qSound : Encl q (-( (85*(Real.pi/(11/5))*Real.log 9-({period}:ℤ)*(2*Real.pi))^2)) := by
  have h := mulSound reducedSound reducedSound (show mulOK reduced reduced (negative q) by decide)
  have he : negative (negative q)=q := by decide
  simpa only [he,pow_two] using negativeSound h
theorem halfSound : Encl halfNu (85*(Real.pi/(11/5))/2) := by
  simpa only [Int.cast_ofNat] using divSound nuSound (pointSound 2) (show divOK nu (point 2) halfNu by decide)
theorem atanArgSound : Encl atanArg ((85*(Real.pi/(11/5))/2)/(2049/4)) :=
  divSound halfSound shiftSound (by decide)
theorem vSound : Encl v (-(((85*(Real.pi/(11/5))/2)/(2049/4))^2)) := by
  have h := mulSound atanArgSound atanArgSound (show mulOK atanArg atanArg (negative v) by decide)
  have he : negative (negative v)=v := by decide
  simpa only [he,pow_two] using negativeSound h
"""
names=["collarSound","shiftSound","nu0Sound","nuSound","argSound","periodicSound","reducedSound","qSound","halfSound","atanArgSound","vSound"]
for name,e in [("Sin",1),("Cos",0),("Atan",None)]:
    t=traces[name];nt=len(t["t"]);nm=len(t["mid"])
    start="reduced" if e==1 else "point 1" if e==0 else "atanArg"
    vname="q" if e is not None else "v"
    hs="reducedSound" if e==1 else "(by simpa only [Int.cast_one] using pointSound 1)" if e==0 else "atanArgSound"
    hv="qSound" if e is not None else "vSound"
    receiver="sinTrace" if e==1 else "cosTrace" if e==0 else "CRALP.SignedGridSeries.atanTrace"
    step=f"trigStep {e} {vname} t mid i" if e is not None else "logStep v t mid i"
    summand="t" if e is not None else "mid"
    guard="reduced" if e is not None else "atanArg"
    gsound="reducedSound" if e is not None else "atanArgSound"
    p,qden=(4,1) if e is not None else (1,3)
    arg="(85*(Real.pi/(11/5))*Real.log 9-("+str(period)+":ℤ)*(2*Real.pi))" if e is not None else "((85*(Real.pi/(11/5))/2)/(2049/4))"
    fn="sin" if e==1 else "cos" if e==0 else "arctan"
    lean+=f"""namespace {name}
def ts : Array Grid := #[{arr(t["t"])}]
def ms : Array Grid := #[{arr(t["mid"])}]
theorem tSize : ts.size = {nt} := by decide
theorem mSize : ms.size = {nm} := by decide
def t (i : ℕ) : Grid := if h : i < {nt} then ts[i]'(by rw [tSize]; exact h) else ⟨1,0⟩
def mid (i : ℕ) : Grid := if h : i < {nm} then ms[i]'(by rw [mSize]; exact h) else ⟨1,0⟩
def out : Grid := ⟨{t["out"][0]},{t["out"][1]}⟩
theorem initial : t 0 = {start} := by decide
theorem allSteps : ∀ i : Fin {399 if e is not None else 400}, {step} := by decide
theorem lowerSum : out.lo=(∑ i ∈ range 400, ({summand} i).lo)-1 := by decide
theorem upperSum : out.hi=(∑ i ∈ range 400, ({summand} i).hi)+1 := by decide
theorem actualReduced : Encl out (Real.{fn} {arg}) := by
  have hg := absSound {gsound} (show absGuard {guard} {p} {qden} by decide)
  have hg' : |{arg}| ≤ {p}/{qden} := by simpa only [Int.cast_one,Int.cast_ofNat] using hg
  have hg'' : |{arg}| ≤ {"4" if e is not None else "1/3"} := by simpa using hg'
  exact {receiver} {vname} t mid out hg'' (by rw [initial]; exact {hs}) {hv} allSteps lowerSum upperSum
"""
    names.extend(name+"."+n for n in ["tSize","mSize","initial","allSteps","lowerSum","upperSum","actualReduced"])
    if e is not None:
        lean+=f"""theorem actual : Encl out (Real.{fn} (85*(Real.pi/(11/5))*Real.log 9)) :=
  {fn}PeriodAttach {period} actualReduced
"""
        names.append(name+".actual")
    lean+=f"end {name}\n"
lean+="end CRALP.SignedPilot\n"+''.join(f"#print axioms CRALP.SignedPilot.{n}\n" for n in names)
(out/"SignedPilot.lean").write_text(lean)
data={"source_blob":blob,"scale":str(S),"period":period,"pi":ends(pi),"log9":ends(lg),"collar":ends(collar),"nu0":ends(env["NUF"]),"nu":ends(nu),"arg":ends(x),"periodic":ends(period*env["TWO_PI"]),"reduced":ends(r),"q":ends(q),"halfNu":ends(y),"shift":ends(shift),"atanArg":ends(w),"v":ends(v),"traces":traces,"sha256":hashlib.sha256(lean.encode()).hexdigest()}
(out/"signed-pilot.json").write_text(json.dumps(data,indent=2)+"\n")
print("period",period,"widths",{k:v["width"] for k,v in traces.items()})
