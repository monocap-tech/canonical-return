from pathlib import Path
import numpy as np
from fractions import Fraction as F
root=Path(__file__).resolve().parent
chunks=[np.load(root/f'odd_U_{a}_{b}.npz') for a,b in [(1,2047),(2049,6143),(6145,11263),(11265,16383)]]
ms=np.concatenate([a['indices'] for a in chunks]);Y=np.concatenate([a['Y'] for a in chunks]);sc=int(chunks[0]['scale']);errs=np.concatenate([a['tile_errors'] for a in chunks]);err=sum(map(F,map(str,errs)),F(0));assert np.array_equal(ms,np.arange(1,16385,2)) and all(int(a['scale'])==sc for a in chunks) and all(F(str(a['error']))==sum(map(F,map(str,a['tile_errors'])),F(0)) for a in chunks)
np.savez_compressed(root/'odd_U_full.npz',Y=Y,indices=ms,scale=sc,error=str(err),tile_errors=errs);print('combined full U head PASS',len(ms),float(err))
