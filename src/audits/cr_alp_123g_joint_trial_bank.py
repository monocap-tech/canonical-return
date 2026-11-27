#!/usr/bin/env python3
"""Extract rational joint pole/capped-tail witnesses with correlated interval residuals."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,sys,json
from cr_alp_123b_pole_residual import solve
from cr_alp_123a_pole_schur import enc
SRC=Path(__file__).with_name('cr_alp_123e_complete_bridge.py')
PIN='c92c1c9e341071bc43c5e40994a1111e0fdbeaa9'
D=10**90
def witness(A,base,C,keep,nl,high,low,c,st,B,k,S,I,ZERO,root):
    n=len(A);j0=n-1
    P=[[(A[i][j].l+A[i][j].h+A[j][i].l+A[j][i].h)//4 for j in range(j0)] for i in range(j0)]
    rhs=[-(A[i][-1].l+A[i][-1].h)//2 for i in range(j0)]
    x=[F(round(v*D),D) for v in solve(P,rhs)]+[F(1)]
    w=[I(v) for v in x]
    res=[sum((v*z for v,z in zip(row,w)),ZERO) for row in A]
    en=sum((v*z for v,z in zip(res,w)),ZERO)
    assert en.l>0
    # z are ordinary historical E85 coefficients. q is the negative of the
    # auxiliary tail trial y in 123F. Joint residual finite = a-R*y.
    q={high[p]:w[nl+i] for i,p in enumerate(keep)}
    tau=sum((v*c[j] for j,v in q.items()),ZERO)/st
    f={j:w[i] for i,j in enumerate(low)}
    f[k]-=sum((f[j]*c[j]/c[k] for j in low if j!=k),ZERO)
    fb=(1-w[low.index(k)]*c[k])/B
    # q= sum q_j e_j - tau*(b-sum_low c_i e_i)
    # q has exact pole nullness by symbolic construction, not interval overlap.
    vraw=(root/'cr_alp_049'/f'result_{k}.json').read_bytes()
    assert hashlib.sha256(vraw).hexdigest()=={1:'47c63e3ac654b0992c5dc97d3aa4287d92aec45b2241debeaad16ccbe9279ee3',2:'9632a259df1194aa7fbde145081b2918437db7e0387fde365d2ae794f48b930c'}[k]
    vb=json.loads(vraw);vs=int(vb['scale']);assert S%vs==0
    V=[[I(int(l)*(S//vs),int(h)*(S//vs),True) for l,h in row] for row in vb['V']]
    assert len(V)==nl and all(len(row)==nl for row in V)
    cond=[sum((V[i][j]*res[i] for i in range(nl)),ZERO) for j in range(nl)]
    cond2=sum((v.sq() for v in cond),ZERO)
    # Compute residual functionals on the available Q85 e_j, including omitted
    # high modes, with exact complete bridge. These are actual sine coefficients
    # of the projected auxiliary error, because <e,Q85 e_j>=<e,e_j>.
    coeff=[]
    for i in range(len(high)):
        r=base[nl+i][-1]+sum((base[nl+i][j]*w[j] for j in range(nl)),ZERO)+sum((C[i][p]*w[nl+j] for j,p in enumerate(keep)),ZERO)
        coeff.append(r)
    tail2=sum((v.sq() for v in coeff),ZERO)
    # Upper bound on norm of source f and tail q is NOT a residual upper bound.
    # Preserve unrounded analytic coefficient intervals to retain correlated
    # cancellation in future tail bounds; exact rational trial is x.
    def iv(v):return [str(F(v.l,S)),str(F(v.h,S))]
    result={'trial_order':low+[high[p] for p in keep]+['normalized_pole'],
        'trial_rational_coordinates':[str(v) for v in x],
        'joint_trial_energy':iv(en),'finite_conditioned_residual_norm2':iv(cond2),
        'available_projected_auxiliary_residual_norm2':iv(tail2),
        'available_projected_auxiliary_residual_sine_coefficients':{str(j):iv(v) for j,v in zip(high,coeff)},
        'source_sines':{str(j):iv(v) for j,v in f.items()},'source_pole':iv(fb),
        'tail_precursor':{str(j):iv(v) for j,v in q.items()},'tail_projection_tau':iv(tau),
        'tail_full_pole':iv(-tau)}
    print('joint witness parity',k,'cutoff',high[keep[-1]],'energy',float(F(en.l,S)), 'conditioned finite residual norm2 upper',float(F(cond2.h,S)),'available auxiliary tail norm2 lower',float(F(tail2.l,S)),flush=True)
    return result

def main(root,previous,out):
    raw=SRC.read_bytes()
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==PIN
    source=raw.decode()
    old="                rows.append(row);print('parity',k,'cutoff',N,'theta',t,'Schur',float(lo),float(hi),label,flush=True)"
    assert source.count(old)==1
    source=source.replace(old,old+"\n                if t==0: row['joint_trial']=witness(A,base,C,keep,nl,high,low,c,st,B,k,S,I,ZERO,root)")
    source=source.replace("'CONDITIONAL_COMPLETE_BRIDGE_COMPARATOR_CENSUS'","'CONDITIONAL_JOINT_AUXILIARY_TRIAL_BANK'")
    g={'__file__':str(SRC),'__name__':'cr_alp_123g_trial_bank','witness':witness}
    exec(compile(source,str(SRC),'exec'),g)
    g['main'](root,previous,out)
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
