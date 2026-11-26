"""Merge the authenticated old 144 and CT's new 16 odd responses.
All arithmetic error and squared-tail bounds inherited by block concatenation."""
from pathlib import Path
from fractions import Fraction as F
import numpy as np,sys
old,new,out=map(Path,sys.argv[1:4]);out.mkdir(exist_ok=True,parents=True)
A=np.load(old/'odd_merged_corrector.npz');B=np.load(new/'odd_corrector.npz')
ai,bi=A['indices'],B['indices'];newi=B['new'];S=10**60
assert A['coefficient'].shape==(len(ai),144) and B['coefficient'].shape==(len(bi),16)
assert len(ai)==len(set(ai))==283 and len(bi)==len(set(bi))
assert len(newi)==len(set(newi))==16 and not np.intersect1d(ai,newi).size
assert np.array_equal(np.setdiff1d(bi,ai),np.sort(newi))
assert int(A['coefficient_scale'])==int(B['coefficient_scale'])==S
idx=np.union1d(ai,bi);assert len(idx)==299
C=np.zeros((len(idx),160),object)
C[np.searchsorted(idx,ai),:144]=np.vectorize(int,otypes=[object])(A['coefficient'])
C[np.searchsorted(idx,bi),144:]=np.vectorize(int,otypes=[object])(B['coefficient'])
err=F(str(A['coefficient_error']))+F(str(B['coefficient_error']))
assert err>=0
np.savez_compressed(out/'odd_merged_corrector.npz',indices=idx,coefficient=np.array(C,str),coefficient_scale=str(S),coefficient_error=str(err))
P=np.load(old/'odd_combined_projection.npz');U=np.load(new/'odd_U_full.npz')
assert P['Y'].shape==(8192,263) and U['Y'].shape==(8192,16)
assert int(P['d0'])==119 and int(P['Yscale'])==int(U['scale'])==10**12
assert np.array_equal(U['indices'],np.arange(1,16384,2))
assert F(str(U['error']))==sum(map(F,map(str,U['tile_errors'])),F(0))
e=F(str(P['U_error']))+F(str(U['error']));assert e>=0
np.savez_compressed(out/'odd_combined_projection.npz',Y=np.column_stack([P['Y'],U['Y']]),Yscale=P['Yscale'],head_error=P['head_error'],U_error=str(e),d0=119)
T=np.load(old/'odd_merged_tail.npz');V=np.load(new/'odd_U_tail.npz')
assert np.array_equal(T['cutoff'],V['cutoff'])
tau=[F(str(a))+F(str(b)) for a,b in zip(T['tail'],V['tail'])]
assert all(v>=0 for v in tau)
np.savez_compressed(out/'odd_merged_tail.npz',cutoff=T['cutoff'],tail=np.array(list(map(str,tau))))
print('CV MERGE PASS',len(idx),160,'coefficient error',float(err),'action error',float(e),'squared tail',float(tau[0]),flush=True)
