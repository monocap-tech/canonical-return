#!/usr/bin/env python3
"""CR-ALP-092: receive complete capped inverse Gram on the obstruction block."""
from pathlib import Path
from fractions import Fraction as F
import sys,json
import numpy as np
from cr_alp_086_terminal_moment import data,full,exact_spd
from cr_alp_085_obstruction_block import mm,transpose,inv
from cr_alp_075_shifted_energy import S,floor
from cr_alp_082_inertia_receiver import inertia
def fmat(A):return np.asarray([[float(x) for x in row] for row in A])
def common(work):
    def load(tag):
        rows=[r.split() for r in (work/tag/'trial.txt').read_text().splitlines()]
        assert all(F(l)==F(h) for _,l,h in rows)
        return {int(j):F(l) for j,l,h in rows}
    def dec(x):
        n=x*10**12;assert n.denominator==1;n=n.numerator
        return ('-' if n<0 else '')+f'{abs(n)//10**12}.{abs(n)%10**12:012d}'
    banks=[load(str(i)) for i in range(3)]
    assert all(list(x)==list(banks[0]) for x in banks)
    for tag in ('01','02','12'):
        a,b=map(int,tag);p=work/tag;p.mkdir(exist_ok=True)
        (p/'trial.txt').write_text(''.join(f'{j} {dec(banks[a][j]+banks[b][j])} {dec(banks[a][j]+banks[b][j])}\n' for j in banks[0]))
        (p/'diag.txt').write_bytes((work/'0'/'diag.txt').read_bytes())
def main(root,prior,block,previous,quadratic,out):
    g=json.loads(previous.read_text());assert g['rank']==83
    q=json.loads(quadratic.read_text());assert q['status']=='CONDITIONAL_CAPPED_BLOCK_QUADRATIC'
    report={}
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        n=len(J);assert inertia(full(H,D,W,DP,A,F(1,10)))[0]>=n-1
        z=q['sectors'][k]
        I={tag:tuple(map(F,v['majorant'])) for tag,v in z.items()}
        assert set(I)=={'0','1','2','01','02','12'}
        ints=[[None]*3 for _ in range(3)]
        for i in range(3):
            ints[i][i]=I[str(i)]
            for j in range(i+1,3):
                l,h=I[f'{i}{j}'];a,b=I[str(i)];c,d=I[str(j)]
                ints[i][j]=ints[j][i]=((l-b-d)/2,(h-a-c)/2)
        B=[[(l+h)/2 for l,h in r] for r in ints]
        for i in range(3):B[i][i]+=sum((h-l)/2 for l,h in ints[i])
        assert exact_spd(B)
        GL=[[F(int(v),S) for v in r] for r in g['sectors'][k]['correction_gram_lower_grid']]
        GU=[[F(int(v),S) for v in r] for r in g['sectors'][k]['correction_gram_upper_grid']]
        K=[[F(H[i][j],S)+GL[i][j] for j in range(n)] for i in range(n)]
        assert exact_spd(K)
        # B_C <= B_L+G. Compare on complement of Z in various
        # positive metric choices, using the actual block bound B_C|Z <= B.
        BP=mm(transpose(W),mm(A,W))
        def global_bound(t):
            if not t:return [[F(D[i],S)*(i==j)+GU[i][j] for j in range(n)] for i in range(n)]
            return [[(1+t)*(F(D[i],S)*(i==j)-DP[i][j])+(1+1/t)*BP[i][j]+GU[i][j] for j in range(n)] for i in range(n)]
        hf=fmat(K);zf=fmat(Z);bf=fmat(B)
        props=[]
        ts=[F(0),F(1,20),F(1,10),F(1,5),F(3,10)]
        for t in ts:
            M=global_bound(t);mf=fmat(M)
            if np.linalg.eigvalsh(mf)[0]<=0:continue
            for metric in ('global','energy','original_D'):
                Fm=mf if metric=='global' else hf if metric=='energy' else np.diag(np.array(D,dtype=float)/S)
                Wf=np.linalg.solve(zf.T@Fm@zf,zf.T@Fm)
                V=np.eye(n)-zf@Wf
                L=V.T@mf@V;R=Wf.T@bf@Wf
                for e in sorted({F(i,100) for i in range(1,101)}|{F(i,10) for i in range(11,51)}):
                    ef=float(e)
                    score=float(np.linalg.eigvalsh(hf-(1+ef)*L-(1+1/ef)*R)[0])
                    props.append((score,t,metric,e))
        assert props
        props.sort(reverse=True);best=n-1;trials=[]
        for score,t,metric,e in props[:8]:
            M=global_bound(t);assert exact_spd(M)
            Fm=M if metric=='global' else K if metric=='energy' else [[F(D[i],S)*(i==j) for j in range(n)] for i in range(n)]
            FZ=mm(Fm,Z);W1=mm(inv(mm(transpose(Z),FZ)),transpose(FZ))
            P=mm(Z,W1);V=[[F(i==j)-P[i][j] for j in range(n)] for i in range(n)]
            L=mm(transpose(V),mm(M,V));R=mm(transpose(W1),mm(B,W1))
            C=[[floor(S*(K[i][j]-(1+e)*L[i][j]-(1+1/e)*R[i][j])) for j in range(n)] for i in range(n)]
            assert all(C[i][j]==C[j][i] for i in range(n) for j in range(n))
            for i in range(n):C[i][i]-=n-1
            ii=inertia(C);best=max(best,ii[0]);trials.append({'global_eta':str(t),'metric':metric,'split_eta':str(e),'inertia':ii,'numeric_min_eigenvalue':score})
            if best==n:break
        report[k]={'rank':best,'pass_full_sector':best==n,'block_majorant':[[str(x) for x in r] for r in B],
                   'quadratic_intervals':z,'trials':trials}
        print(k,'RANK',best,'best numeric margin',props[0][0],flush=True)
    rank=sum(v['rank'] for v in report.values())
    ans={'status':'CONDITIONAL_CAPPED_BLOCK_RECEIVER','rank':rank,'conditional_index_upper':85-rank,
         'sectors':report,'premises':'Inherited original source and Gamma comparison, common finite trial projection, exact polarization of a single weighted residual quadratic majorant, authenticated 091 correction Gram. Exact rational congruence. No unrestricted RH claim.'}
    out.write_text(json.dumps(ans,indent=2)+'\n')
if __name__=='__main__':
    if sys.argv[1]=='common':common(Path(sys.argv[2]))
    else:main(*map(Path,sys.argv[1:]))
