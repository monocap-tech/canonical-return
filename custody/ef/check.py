from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import numpy as np,json
root=Path(__file__).resolve().parents[2];out=root/'custody/ef/generated';p=root/'ce-output/pristine/cr_int_t1d_ef_pilot'
old=np.load(root/'custody/dv/generated/merged/odd_combined_projection.npz')
assert old['Y'].shape==(8192,471) and int(old['Yscale'])==10**12
records=[]
for lo,hi in [(15873,16383),(16385,18431),(65025,65535)]:
    a=np.load(p/f'odd_U_{lo}_{hi}.npz');Y=a['Y'];ids=a['indices'];S=int(a['scale']);e=F(str(a['error']))
    assert np.array_equal(ids,np.arange(lo,hi+1,2)) and Y.shape==(len(ids),352) and S==10**12
    assert e==sum(map(lambda x:F(str(x)),a['tile_errors']),F(0)) and e>0
    rec=dict(lo=lo,hi=hi,rows=len(ids),columns=352,error=str(e),error_float=float(e),error_per_row=float(e/len(ids)))
    if hi==16383:
        D=Y.astype(object)-old['Y'][(ids-1)//2,119:].astype(object)
        v=sum(int(x)**2 for x in D.flat);bound=e+F(str(old['U_error']))
        assert F(v,S*S)<=bound*bound
        rec['old_overlap_frobenius_upper']=float(F(isqrt(v)+1,S));rec['overlap_pass']=True
    records.append(rec)
# This is a cost diagnostic, not an extrapolated finite-error certificate.
a=records[1];estimate=F(a['error'])*F(24576,a['rows'])
report=dict(scope='directed pilot rows only, no full higher-cutoff action or head regeneration',records=records,linear_error_extrapolation_diagnostic=float(estimate),old_physical_error=float(F(str(old['U_error']))))
(out/'pilot.json').write_text(json.dumps(report,indent=2)+'\n')
print('EF PILOT',report,flush=True)
