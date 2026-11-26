#!/usr/bin/env python3
"""Receive a simultaneous compression majorant by exact polarization of certified scalar intervals."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
from cr_alp_078_moment_receiver import sha

def main(root):
    sel=json.loads((root/'selection/result.json').read_text())
    assert sel['status']=='CONDITIONAL_BLOCK_SELECTION' and sel['rank']==84
    s=sel['sectors']['2'];Z=[list(map(F,z)) for z in s['columns']];m=len(Z);assert m==5
    names=[f'd{i}' for i in range(m)]+[f'p{i}_{j}' for i in range(m) for j in range(i+1,m)]
    assert s['cases']==names
    bounds={};reports={}
    for name in names:
        d=root/'selection'/name
        one=json.loads((d/'result.json').read_text())
        ps=json.loads((root/'cases'/name/'result.json').read_text())
        a=one['sectors']['2'];b=ps['sectors']['2']
        assert ps['status']=='CONDITIONAL_EVEN_ADAPTIVE_RESPONSE'
        assert ps['sectors']['1']['certificate']==sel['sectors']['1']['certificate']==one['sectors']['1']['certificate']
        assert ps['sectors']['1']['replayed_inertia']==[43,0,0]
        assert a['adaptive_source']['coordinates']==s['coordinates']
        assert a['adaptive_source']['input_sha256']==b['source_sha256']==sha(d/'input_2_adaptive.txt')
        if name[0]=='d':q=Z[int(name[1:])]
        else:
            i,j=map(int,name[1:].split('_'));q=[x+y for x,y in zip(Z[i],Z[j])]
        assert q==list(map(F,a['adaptive_source']['coefficients']))
        lo,hi=map(F,b['inverse_interval']);assert 0<=lo<=hi
        bounds[name]=(lo,hi);reports[name]=b
    C=[[F(0) for j in range(m)] for i in range(m)]
    R=[[F(0) for j in range(m)] for i in range(m)]
    for i in range(m):
        for j in range(i,m):
            li,ui=bounds[f'd{i}'];lj,uj=bounds[f'd{j}']
            if i==j:lo,hi=li,ui
            else:
                lp,up=bounds[f'p{i}_{j}']
                lo,hi=(lp-ui-uj)/2,(up-li-lj)/2
            assert lo<=hi
            C[i][j]=C[j][i]=(lo+hi)/2;R[i][j]=R[j][i]=(hi-lo)/2
    A=[[C[i][j]+(sum(R[i])+F(1,10**18) if i==j else 0) for j in range(m)] for i in range(m)]
    # For every symmetric error E with |Eij|<=Rij, diag(row sums R)-E is PSD.
    # The extra positive diagonal makes the majorant strict for a PSD compression.
    mat=lambda X:[[str(x) for x in row] for row in X]
    result={'status':'CONDITIONAL_BLOCK_COMPRESSION','rank':84,'conditional_index_upper':1,
      'odd_certificate':sel['sectors']['1']['certificate'],'coordinates':s['coordinates'],
      'columns':s['columns'],'K_compression':s['K_compression'],'center':mat(C),'radius':mat(R),
      'compression_majorant':mat(A),'scalar_intervals':{k:list(map(str,v)) for k,v in bounds.items()},
      'scalar_pass_count':sum(v['classification']=='PASS_SCALAR_C_COMPARISON' for v in reports.values()),
      'premises':'Inherited physical inverse interval premises. Polarization uses exact pair sums; diagonal row-sum reserve proves the simultaneous majorant. No rank increase claimed.'}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS exact polarization of',len(names),'physical intervals; block dimension',m,flush=True)
if __name__=='__main__':main(Path(sys.argv[1]))
