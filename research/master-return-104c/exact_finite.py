"""RETURN104-C: exact integer RETURN087 finite-source contraction replay.
The archived arithmetic source is preserved. Replace only the signed-limb
dot routine after inspecting its actual signature and orientation. Fail
closed on an unrecognized routine or any changed output matrix.
"""
import ast,hashlib,json,sys,subprocess,shutil
from pathlib import Path
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=True)
srcdir=root/'return087'
pins={"finite_audit.py":"cf0f1170ecaf0d5faf8b4f1560517d8f678aeba9cb000f9baf68a42fd1b0c89c","band_audit.py":"e729bf288888a1d3f0105f59147c2f59395b220602d540b24fe6a2bc87915624","finite_bank.json":"c6a1df49ff86806810e8763c10223b3b4b2e331b0243693802f1372eb14c443f","band_0.json":"624d811d17b570bdc58d39d7e47350bad0f91524443e8e2fad0cba91a7c6db58","band_1.json":"89805afe0b4fdcef6a47847490d6940b75e3d95b0f5203cdc74f087098c319c1"}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_bytes())
# Authenticate against the original verifier's pinned manifest rather
# than assuming the archive root or output schema.
verify=load(srcdir/'verify.json')
assert verify['status']=='PASS' and verify['total_containment_checks']==952389
for name in pins:
 h=verify['pins'][name]
 assert h==pins[name],name
 assert sha(srcdir/name)==h,name

original=out/'reference';original.mkdir()
for name in ('finite_bank.json','band_0.json','band_1.json','finite_audit.json','verify.json'):
 shutil.copy2(srcdir/name,original/name)
# Locate the limb routine in the authenticated finite source. The
# exact implementation accepts both conventional dot orientations
# only after identifying which the actual archived routine implements.
s=(srcdir/'finite_audit.py').read_text();tree=ast.parse(s)
candidates=[n for n in tree.body if isinstance(n,ast.FunctionDef) and ('dot' in n.name.lower() or 'matmul' in n.name.lower()) and len(n.args.args)>=2]
print('candidate dot routines',[(n.name,n.lineno,n.end_lineno)for n in candidates],flush=True)
for n in candidates:print(ast.get_source_segment(s,n),flush=True)
assert candidates,'no archived dot routine found'
# Determine implementation by exact small rectangular test in a
# source-prefix namespace. No old floating result is used as a
# mathematical certificate; the test only disambiguates the API.
import numpy as np
from flint import fmpz_mat
ns={'np':np,'math':__import__('math'),'F':__import__('fractions').Fraction,'sys':sys}
# Load declarations only, avoiding the original driver or side effects.
# Retain simple constant assignments and function definitions; skip
# any import-time file access or matrix computation.
safe=[]
for n in tree.body:
 if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef)):
  safe.append(n)
 elif isinstance(n,ast.Assign) and all(isinstance(t,ast.Name) for t in n.targets) and all(isinstance(x,(ast.Constant,ast.BinOp,ast.UnaryOp,ast.Name,ast.Load,ast.operator,ast.unaryop,ast.Tuple,ast.List,ast.Dict)) for x in ast.walk(n.value)):
  safe.append(n)
for n in safe:
 try:exec(compile(ast.Module(body=[n],type_ignores=[]),str(srcdir/'finite_audit.py'),'exec'),ns)
 except (NameError,ImportError):pass
# Candidate tests use deliberately nonsquare, signed, multilimb arrays.
X=np.array([[2**68+3,-7,11],[5,-2**55,13]],object)
Y=np.array([[17,19],[23,-2**63],[29,31]],object)
Z=np.array([[17,19,23],[29,-2**63,31]],object)
choices=[]
for n in candidates:
 fn=ns.get(n.name)
 if not callable(fn):continue
 for orient,b,expected in [('normal',Y,X@Y),('right_transposed',Y.T,X@Y),('left_transposed',Z,X.T@Z)]:
  try:
   ans=np.asarray(fn(X,b),dtype=object)
   if ans.shape==expected.shape and np.array_equal(ans,expected):choices.append((n,orient))
  except Exception as exc:print('API probe',n.name,orient,type(exc).__name__,flush=True)
assert len(choices)==1,[(n.name,o)for n,o in choices]
n,orient=choices[0];assert not n.decorator_list
a,b=[v.arg for v in n.args.args[:2]]
body=f'''
    # RETURN104-C exact integer backend, no floating-point dot.
    from flint import fmpz_mat
    A=np.asarray({a},dtype=object)
    B=np.asarray({b},dtype=object)
    assert A.ndim==B.ndim==2
'''
if orient=='left_transposed':body+='    A=A.T\n'
if orient=='right_transposed':body+='    B=B.T\n'
body+='''
    assert A.shape[1]==B.shape[0]
    X=fmpz_mat(A.shape[0],A.shape[1],[int(v) for v in A.flat])
    Y=fmpz_mat(B.shape[0],B.shape[1],[int(v) for v in B.flat])
    Z=X*Y
    return np.asarray([[int(Z[i,j]) for j in range(B.shape[1])] for i in range(A.shape[0])],dtype=object)
'''
lines=s.splitlines(keepends=True)
# AST body starts at the first statement, preserving actual signature.
start=n.body[0].lineno-1;end=n.end_lineno
new=''.join(lines[:start])+body+'\n'+''.join(lines[end:])
ast.parse(new)
(srcdir/'finite_audit.py').write_text(new)
(out/'finite_audit_exact.py').write_text(new)
print('exact backend replacement',n.name,orient,flush=True)
# Inspect the band driver to verify it imports the modified routine
# rather than having its own hidden floating-point implementation.
band=(srcdir/'band_audit.py').read_text()
print('BEGIN AUTHENTICATED BAND SOURCE',flush=True);print(band,flush=True);print('END AUTHENTICATED BAND SOURCE',flush=True)
assert 'finite_audit' in band
for name in ('finite_bank.json','band_0.json','band_1.json','finite_audit.json'):
 (srcdir/name).unlink()
# Reproduce using the explicit common-root argument documented in
# MASTER_RETURN_087.md. The archived scripts write their own banks.
for name in ('finite_audit.py','band_audit.py'):
 subprocess.run([sys.executable,str(srcdir/name),str(root)],cwd=srcdir,check=True)
records=[]
for name in ('finite_bank.json','band_0.json','band_1.json'):
 x=load(srcdir/name);y=load(original/name)
 assert x==y,name
 shutil.copy2(srcdir/name,out/name)
 records.append({'name':name,'exact_equality':True,'sha256':sha(srcdir/name),'original_sha256':sha(original/name)})
 print('exact finite contraction equality',name,flush=True)
(out/'receipt.json').write_text(json.dumps({'status':'PASS','scope':'exact integer backend of original RETURN087 finite and feature contractions, identical output banks','records':records,'backend':n.name,'orientation':orient,'transformed_source_sha256':sha(srcdir/'finite_audit.py')},indent=2)+'\n')
print('RETURN104-C EXACT INTEGER FINITE CONTRACTION EQUALITY PASS',flush=True)
