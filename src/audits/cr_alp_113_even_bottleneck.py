#!/usr/bin/env python3
"""CR-ALP-113: joint envelope refined by authenticated adaptive response.

All compression bounds are simultaneous necessary constraints for B_C.
Convex inverse metric is sufficient, not claimed complete for >2 constraints.
An exactly feasible vv^T with v^T K^-1 v>1 rejects only the lower-K
relaxation; the separating source q is then selected for a physical audit.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import numpy as np
from scipy.optimize import linprog,minimize
from cr_alp_086_terminal_moment import data,full,exact_spd,entrymap,input_for
from cr_alp_085_obstruction_block import mm,transpose,inv
from cr_alp_075_shifted_energy import S,floor,ceil,positive
from cr_alp_082_inertia_receiver import inertia
from cr_alp_093_envelope_duality import outer,sub,fmat
from cr_alp_078_moment_receiver import sha
def grid(A,upper):
    n=len(A);C=[[floor(v*S) for v in row] for row in A]
    assert C==transpose(C)
    for i in range(n):C[i][i]+=n if upper else -n
    return [[F(v,S) for v in row] for row in C]
def spd(A):
    return positive([[int(v*S) for v in row] for row in grid(A,False)])
def quad(v,A):
    return sum(a*b for a,b in zip(v,[r[0] for r in mm(A,[[x] for x in v])]))
def search(K,terms):
    L=np.linalg.cholesky(fmat(K));n=len(K);m=len(terms)
    # Each term = (label,U,A); evaluate U A^-1 U^T in K metric.
    M=np.array([L.T@fmat(U)@np.linalg.solve(fmat(A),fmat(U).T)@L for _,U,A in terms])
    M=(M+M.transpose(0,2,1))/2
    cuts=[];best=None;lp=None;V0=[]
    def evaluate(w):
        nonlocal best
        e,V=np.linalg.eigh(np.einsum('i,ijk->jk',w,M))
        if best is None or e[0]>best[0]:best=(float(e[0]),w.copy(),V.copy())
        for j in range(4):
            v=V[:,j];cuts.append(np.r_[-np.einsum('i,kil,l->k',v,M,v),1.])
        return float(e[0])
    w=np.zeros(m);w[0]=1;evaluate(w);evaluate(np.ones(m)/m)
    for _ in range(350):
        r=linprog(np.r_[np.zeros(m),-1.],A_ub=np.array(cuts),b_ub=np.zeros(len(cuts)),
                  A_eq=[np.r_[np.ones(m),0.]],b_eq=[1.],bounds=[(0,1)]*m+[(None,None)],method='highs')
        if not r.success:break
        w=np.maximum(0,r.x[:-1]);w/=w.sum();lp=float(r.x[-1]);val=evaluate(w)
        if val>1.00001 or lp-best[0]<1e-7:break
    val,w,V=best
    return val,w,V,lp,M,L
def certificate(K,terms,w):
    a=[max(0,round(float(x)*10**12)) for x in w];a[0]=max(a[0],1)
    weights=[F(x,sum(a)) for x in a]
    assert sum(weights)==1 and min(weights)>=0 and weights[0]>0
    n=len(K);T=[[F(0)]*n for _ in range(n)]
    for t,(_,U,A) in zip(weights,terms):
        if not t:continue
        B=mm(U,mm(inv(A),transpose(U)))
        T=[[x+t*y for x,y in zip(r,s)] for r,s in zip(T,B)]
    C=grid(sub(mm(K,mm(T,K)),K),False)
    ii=inertia([[int(v*S) for v in row] for row in C])
    return {'weights':list(map(str,weights)),'inertia':ii,'full_pass':ii==[n,0,0]}
def witness(K,terms,M,L,V,old):
    n=len(K);rng=np.random.default_rng(9400+n)
    def normalize(u):
        q=np.einsum('i,kij,j->k',u,M,u)
        return .999*u/np.sqrt(max(q))
    starts=[V[:,j] for j in range(6)]
    starts += [np.linalg.solve(L,np.array([float(F(x)) for x in v])) for v in old]
    starts += [V[:,:6]@rng.normal(size=6) for _ in range(8)]
    best=None
    def fun(u):return -float(u@u)
    def jac(u):return -2*u
    def con(u):return 1-np.einsum('i,kij,j->k',u,M,u)
    def cj(u):return -2*np.einsum('kij,j->ki',M,u)
    for v in starts:
        u=normalize(v)
        r=minimize(fun,u,jac=jac,constraints=[{'type':'ineq','fun':con,'jac':cj}],method='SLSQP',
                   options={'maxiter':300,'ftol':1e-10})
        for p in (u,r.x):
            if not np.isfinite(p).all():continue
            p=normalize(p);gain=float(p@p)
            if best is None or gain>best[0]:best=(gain,L@p)
    if best is None:return {}
    gain,v=best;ans={'numeric_rank_one_gain':gain}
    if gain<=1:return ans
    vv=[F(round(float(x)*10**12),10**12) for x in v]
    for label,U,A in terms:
        c=[r[0] for r in mm(transpose(U),[[x] for x in vv])]
        if not spd(sub(A,outer(c))):
            ans['exact_rejection_blocked_by']=label;return ans
    # Rational separating source from K^-1 v; rounding does not enter
    # the exact rejection check. The chosen direction is normalized.
    z=np.linalg.solve(fmat(K),np.array(list(map(float,vv))))
    z=z/(4*np.linalg.norm(z));q=[F(round(float(x)*10**12),10**12) for x in z]
    energy=quad(q,K);gainq=sum(a*b for a,b in zip(q,vv))**2
    if gainq<=energy:return ans
    ans['exact_rejection']={'v':list(map(str,vv)),'separating_source':list(map(str,q)),
            'K_lower_energy':str(energy),'rank_one_energy':str(gainq),
            'strict_gap':str(gainq-energy),'all_constraints_verified':len(terms)}
    return ans
def main(root,bank,prior,block,correction,terminal,capped,received,dual,adaptive,replayed,adaptive2,replayed2,adaptive3,replayed3,adaptive4,replayed4,mixed,mixedreplay,frozen,mixed2,mixedreplay2,frozen2,mixed3,mixedreplay3,frozen3,mixed4,mixedreplay4,frozen4,even,evenreplay,accepted,latest,out):
    old=json.loads(prior.read_text());g=json.loads(correction.read_text());r=json.loads(received.read_text())
    c=json.loads(capped.read_text());p=json.loads((terminal/'terminal_selection.json').read_text());dr=json.loads(dual.read_text())
    assert old['rank']==79 and g['rank']==83 and r['rank']==83 and dr['rank']==83
    ar=json.loads((adaptive/'result.json').read_text());physical=json.loads(replayed.read_text())
    assert ar['status']=='CONDITIONAL_JOINT_ENVELOPE' and physical['status']=='CONDITIONAL_ADAPTIVE_RESPONSE'
    ar2=json.loads((adaptive2/'result.json').read_text());physical2=json.loads(replayed2.read_text())
    assert ar2['status']=='CONDITIONAL_REFINED_JOINT_ENVELOPE' and physical2['status']=='CONDITIONAL_ADAPTIVE_RESPONSE'
    ar3=json.loads((adaptive3/'result.json').read_text());physical3=json.loads(replayed3.read_text())
    assert ar3['status']=='CONDITIONAL_REFINED_JOINT_ENVELOPE' and physical3['status']=='CONDITIONAL_ADAPTIVE_RESPONSE'
    ar4=json.loads((adaptive4/'result.json').read_text());physical4=json.loads(replayed4.read_text())
    assert ar4['status']=='CONDITIONAL_REFINED_JOINT_ENVELOPE' and physical4['status']=='CONDITIONAL_ADAPTIVE_RESPONSE'
    ar5=json.loads((mixed/'result.json').read_text());physical5=json.loads(mixedreplay.read_text())
    fr=json.loads(frozen.read_text())
    assert ar5['status']=='CONDITIONAL_BOTTLENECK_SELECTION' and physical5['status']=='CONDITIONAL_ADAPTIVE_RESPONSE' and fr['rank']==83
    ar6=json.loads((mixed2/'result.json').read_text());physical6=json.loads(mixedreplay2.read_text())
    fr2=json.loads(frozen2.read_text())
    assert ar6['status']=='CONDITIONAL_BOTTLENECK_SELECTION' and physical6['status']=='CONDITIONAL_ADAPTIVE_RESPONSE' and fr2['rank']==83
    ar7=json.loads((mixed3/'result.json').read_text());physical7=json.loads(mixedreplay3.read_text())
    fr3=json.loads(frozen3.read_text())
    assert ar7['status']=='CONDITIONAL_BOTTLENECK_SELECTION' and physical7['status']=='CONDITIONAL_ADAPTIVE_RESPONSE' and fr3['rank']==83
    ar8=json.loads((mixed4/'result.json').read_text());physical8=json.loads(mixedreplay4.read_text())
    fr4=json.loads(frozen4.read_text())
    assert ar8['status']=='CONDITIONAL_BOTTLENECK_SELECTION' and physical8['status']=='CONDITIONAL_ADAPTIVE_RESPONSE' and fr4['rank']==83
    ar9=json.loads((even/'result.json').read_text());physical9=json.loads(evenreplay.read_text());fr5=json.loads(accepted.read_text())
    assert ar9['status']=='CONDITIONAL_EVEN_BOTTLENECK_SELECTION' and physical9['status']=='CONDITIONAL_EVEN_ADAPTIVE_RESPONSE' and fr5['rank']==84 and fr5['conditional_index_upper']==1
    latest_result=json.loads(latest.read_text());assert latest_result['rank']==84 and latest_result['conditional_index_upper']==1
    E=entrymap(bank);report={}
    for k,J,H,D,Z,A,W,DP in data(root,prior,block):
        n=len(J);assert inertia(full(H,D,W,DP,A,F(1,10)))[0]>=n-1
        B=[list(map(F,row)) for row in r['sectors'][k]['block_majorant']]
        GL=[[F(int(v),S) for v in row] for row in g['sectors'][k]['correction_gram_lower_grid']]
        GU=[[F(int(v),S) for v in row] for row in g['sectors'][k]['correction_gram_upper_grid']]
        K=[[F(H[i][j],S)+GL[i][j] for j in range(n)] for i in range(n)];assert spd(K)
        eye=[[F(i==j) for j in range(n)] for i in range(n)]
        BP=mm(transpose(W),mm(A,W));terms=[]
        def attach(label,U,M):
            M=grid(M,True);assert spd(M);terms.append((label,U,M))
        for t in [F(0),F(1,20),F(1,10),F(1,5),F(3,10)]:
            M=[[F(D[i],S)*(i==j)+GU[i][j] for j in range(n)] for i in range(n)] if not t else [[(1+t)*(F(D[i],S)*(i==j)-DP[i][j])+(1+1/t)*BP[i][j]+GU[i][j] for j in range(n)] for i in range(n)]
            attach('global_'+str(t),eye,M)
        attach('physical_C_block',Z,B)
        GZ=mm(transpose(Z),mm(GU,Z))
        attach('transferred_L_block',Z,[[a+b for a,b in zip(x,y)] for x,y in zip(A,GZ)])
        q=list(map(F,p[k]['q']));z=[sum(a*b for a,b in zip(row,q)) for row in Z]
        assert c['sectors'][k]['classification']=='PASS_SCALAR_C_COMPARISON'
        attach('physical_C_terminal',[[x] for x in z],[[F(c['sectors'][k]['inverse_interval'][1])]])
        aa=[F(old['sectors'][k]['inverse_bounds'][str(j)]) for j in J]
        assert sum(x*S/F(d) for x,d in zip(aa,D))<=1
        for i,a in enumerate(aa):
            attach('transferred_coordinate_'+str(J[i]),[[F(j==i)] for j in range(n)],[[a+GU[i][i]]])
        ap=ar['sectors'][k]['adaptive_source'];ps=physical['sectors'][k]
        assert ap['coordinates']==J and ps['source_sha256']==ap['input_sha256']
        q0=list(map(F,ap['coefficients']));assert len(q0)==n
        raw,eps=input_for(bank,E,k,J,eye,q0)
        f=adaptive/f'input_{k}_adaptive.txt'
        assert raw==f.read_text() and sha(f)==ap['input_sha256'] and eps==F(ap['error_upper'])
        vprev=list(map(F,ar['sectors'][k]['exact_rejection']['v']))
        previous_penalty=sum(x*y for x,y in zip(vprev,q0))**2
        bnew=F(ps['inverse_interval'][1])
        assert ps['witness_excluded'] and ps['positive_against_joint_lower_energy']
        assert 0<bnew<quad(q0,K)<previous_penalty
        attach('physical_C_adaptive_095',[[x] for x in q0],[[bnew]])
        # Check exclusion after the outward grid reserve, not before it.
        assert terms[-1][2][0][0]<previous_penalty
        ap=ar2['sectors'][k]['adaptive_source'];ps=physical2['sectors'][k]
        assert ap['coordinates']==J and ps['source_sha256']==ap['input_sha256']
        q0=list(map(F,ap['coefficients']));assert len(q0)==n
        raw,eps=input_for(bank,E,k,J,eye,q0)
        f=adaptive2/f'input_{k}_adaptive.txt'
        assert raw==f.read_text() and sha(f)==ap['input_sha256'] and eps==F(ap['error_upper'])
        vprev=list(map(F,ar2['sectors'][k]['exact_rejection']['v']))
        previous_penalty=sum(x*y for x,y in zip(vprev,q0))**2
        bnew2=F(ps['inverse_interval'][1])
        assert ps['witness_excluded'] and ps['positive_against_joint_lower_energy']
        assert 0<bnew2<quad(q0,K)<previous_penalty
        attach('physical_C_adaptive_097',[[x] for x in q0],[[bnew2]])
        # Check exclusion after the outward grid reserve, not before it.
        assert terms[-1][2][0][0]<previous_penalty
        ap=ar3['sectors'][k]['adaptive_source'];ps=physical3['sectors'][k]
        assert ap['coordinates']==J and ps['source_sha256']==ap['input_sha256']
        q0=list(map(F,ap['coefficients']));assert len(q0)==n
        raw,eps=input_for(bank,E,k,J,eye,q0)
        f=adaptive3/f'input_{k}_adaptive.txt'
        assert raw==f.read_text() and sha(f)==ap['input_sha256'] and eps==F(ap['error_upper'])
        vprev=list(map(F,ar3['sectors'][k]['exact_rejection']['v']))
        previous_penalty=sum(x*y for x,y in zip(vprev,q0))**2
        bnew3=F(ps['inverse_interval'][1])
        assert ps['witness_excluded'] and ps['positive_against_joint_lower_energy']
        assert 0<bnew3<quad(q0,K)<previous_penalty
        attach('physical_C_adaptive_099',[[x] for x in q0],[[bnew3]])
        # Check exclusion after the outward grid reserve, not before it.
        assert terms[-1][2][0][0]<previous_penalty
        ap=ar4['sectors'][k]['adaptive_source'];ps=physical4['sectors'][k]
        assert ap['coordinates']==J and ps['source_sha256']==ap['input_sha256']
        q0=list(map(F,ap['coefficients']));assert len(q0)==n
        raw,eps=input_for(bank,E,k,J,eye,q0)
        f=adaptive4/f'input_{k}_adaptive.txt'
        assert raw==f.read_text() and sha(f)==ap['input_sha256'] and eps==F(ap['error_upper'])
        vprev=list(map(F,ar4['sectors'][k]['exact_rejection']['v']))
        previous_penalty=sum(x*y for x,y in zip(vprev,q0))**2
        bnew4=F(ps['inverse_interval'][1])
        assert ps['witness_excluded'] and ps['positive_against_joint_lower_energy']
        assert 0<bnew4<quad(q0,K)<previous_penalty
        attach('physical_C_adaptive_101',[[x] for x in q0],[[bnew4]])
        # Check exclusion after the outward grid reserve, not before it.
        assert terms[-1][2][0][0]<previous_penalty
        sec=ar5['sectors'][k];ap=sec['adaptive_source'];ps=physical5['sectors'][k]
        assert sec['labels']==fr['sectors'][k]['labels']==[t[0] for t in terms]
        assert ap['coordinates']==J and ps['source_sha256']==ap['input_sha256']
        assert sec['selection_kind']==ps['selection_kind']
        q0=list(map(F,ap['coefficients']));assert len(q0)==n
        raw,eps=input_for(bank,E,k,J,eye,q0)
        f=mixed/f'input_{k}_adaptive.txt'
        assert raw==f.read_text() and sha(f)==ap['input_sha256'] and eps==F(ap['error_upper'])
        bnew5=F(ps['inverse_interval'][1]);assert bnew5>0
        target=quad(q0,K);assert target==F(ps['joint_lower_energy'])
        attach('physical_C_adaptive_103',[[x] for x in q0],[[bnew5]])
        rounded_bound=terms[-1][2][0][0]
        excluded=None
        if sec['selection_kind']=='exact_rejection':
            assert sec['exact_rejection']==fr['sectors'][k]['exact_rejection']
            vprev=list(map(F,sec['exact_rejection']['v']))
            previous_penalty=sum(x*y for x,y in zip(vprev,q0))**2
            assert previous_penalty==F(sec['exact_rejection']['rank_one_energy'])>target
            excluded=rounded_bound<previous_penalty
        else:
            assert sec['selection_kind']=='dual_bottleneck' and 'exact_rejection' not in sec
            assert target==F(sec['K_lower_energy'])
        # A valid upper bound is a necessary constraint even if the
        # individual scalar comparison or previous witness exclusion fails.
        sec6=ar6['sectors'][k];ap6=sec6['adaptive_source'];ps6=physical6['sectors'][k]
        assert sec6['labels']==fr2['sectors'][k]['labels']==[t[0] for t in terms]
        assert ap6['coordinates']==J and ps6['source_sha256']==ap6['input_sha256']
        assert sec6['selection_kind']==ps6['selection_kind']=='dual_bottleneck'
        assert 'exact_rejection' not in sec6 and ps6['witness_excluded'] is None
        q6=list(map(F,ap6['coefficients']));assert len(q6)==n
        raw6,eps6=input_for(bank,E,k,J,eye,q6)
        f6=mixed2/f'input_{k}_adaptive.txt'
        assert raw6==f6.read_text() and sha(f6)==ap6['input_sha256'] and eps6==F(ap6['error_upper'])
        bnew6=F(ps6['inverse_interval'][1]);assert bnew6>0
        target6=quad(q6,K);assert target6==F(ps6['joint_lower_energy'])==F(sec6['K_lower_energy'])
        attach('physical_C_adaptive_105',[[x] for x in q6],[[bnew6]])
        rounded6=terms[-1][2][0][0]
        sec7=ar7['sectors'][k];ap7=sec7['adaptive_source'];ps7=physical7['sectors'][k]
        assert sec7['labels']==fr3['sectors'][k]['labels']==[t[0] for t in terms]
        assert ap7['coordinates']==J and ps7['source_sha256']==ap7['input_sha256']
        assert sec7['selection_kind']==ps7['selection_kind']=='dual_bottleneck'
        assert 'exact_rejection' not in sec7 and ps7['witness_excluded'] is None
        q7=list(map(F,ap7['coefficients']));assert len(q7)==n
        raw7,eps7=input_for(bank,E,k,J,eye,q7)
        f7=mixed3/f'input_{k}_adaptive.txt'
        assert raw7==f7.read_text() and sha(f7)==ap7['input_sha256'] and eps7==F(ap7['error_upper'])
        bnew7=F(ps7['inverse_interval'][1]);assert bnew7>0
        target7=quad(q7,K);assert target7==F(ps7['joint_lower_energy'])==F(sec7['K_lower_energy'])
        attach('physical_C_adaptive_107',[[x] for x in q7],[[bnew7]])
        rounded7=terms[-1][2][0][0]
        sec8=ar8['sectors'][k];ap8=sec8['adaptive_source'];ps8=physical8['sectors'][k]
        assert sec8['labels']==fr4['sectors'][k]['labels']==[t[0] for t in terms]
        assert ap8['coordinates']==J and ps8['source_sha256']==ap8['input_sha256']
        assert sec8['selection_kind']==ps8['selection_kind']=='dual_bottleneck'
        assert 'exact_rejection' not in sec8 and ps8['witness_excluded'] is None
        q8=list(map(F,ap8['coefficients']));assert len(q8)==n
        raw8,eps8=input_for(bank,E,k,J,eye,q8)
        f8=mixed4/f'input_{k}_adaptive.txt'
        assert raw8==f8.read_text() and sha(f8)==ap8['input_sha256'] and eps8==F(ap8['error_upper'])
        bnew8=F(ps8['inverse_interval'][1]);assert bnew8>0
        target8=quad(q8,K);assert target8==F(ps8['joint_lower_energy'])==F(sec8['K_lower_energy'])
        attach('physical_C_adaptive_109',[[x] for x in q8],[[bnew8]])
        rounded8=terms[-1][2][0][0]
        fs=fr5['sectors'][k];sec9=ar9['sectors'][k];ps9=physical9['sectors'][k]
        assert fs['labels']==sec9['labels']==[t[0] for t in terms]
        if k=='1':
            assert fs['rank']==n and fs['certificate']['full_pass']
            assert sec9['selection_kind']==ps9['selection_kind']=='retained_exact_certificate'
            assert fs['certificate']==sec9['certificate']==ps9['certificate']
            weights=list(map(F,fs['certificate']['weights']))
            assert len(weights)==len(terms) and min(weights)>=0 and sum(weights)==1 and weights[0]>0
            T=[[F(0)]*n for _ in range(n)]
            for weight,(_,U,A) in zip(weights,terms):
                if not weight:continue
                B=mm(U,mm(inv(A),transpose(U)))
                T=[[x+weight*y for x,y in zip(row,col)] for row,col in zip(T,B)]
            C=grid(sub(mm(K,mm(T,K)),K),False)
            ii=inertia([[int(x*S) for x in row] for row in C])
            assert ii==fs['certificate']['inertia']==sec9['replayed_inertia']==ps9['replayed_inertia']==[n,0,0]
            ls=latest_result['sectors'][k]
            assert ls['rank']==n and ls['labels']==fs['labels'] and ls['certificate']==fs['certificate']
            report[k]={'rank':n,'labels':fs['labels'],'selection_kind':'retained_exact_certificate',
                       'certificate':fs['certificate'],'replayed_inertia':ii}
            print(k,'PASS frozen exact certificate',ii,flush=True)
            continue
        assert k=='2' and fs['rank']==n-1
        ap9=sec9['adaptive_source']
        assert ap9['coordinates']==J and ps9['source_sha256']==ap9['input_sha256']
        assert sec9['selection_kind']==ps9['selection_kind']=='dual_bottleneck'
        assert 'exact_rejection' not in sec9 and ps9['witness_excluded'] is None
        q9=list(map(F,ap9['coefficients']));assert len(q9)==n
        raw9,eps9=input_for(bank,E,k,J,eye,q9)
        f9=even/f'input_{k}_adaptive.txt'
        assert raw9==f9.read_text() and sha(f9)==ap9['input_sha256'] and eps9==F(ap9['error_upper'])
        bnew9=F(ps9['inverse_interval'][1]);assert bnew9>0
        target9=quad(q9,K);assert target9==F(ps9['joint_lower_energy'])==F(sec9['K_lower_energy'])
        attach('physical_C_adaptive_111',[[x] for x in q9],[[bnew9]])
        rounded9=terms[-1][2][0][0]
        ls=latest_result['sectors'][k];assert ls['labels']==[t[0] for t in terms]
        assert ls['rank']==n-1 and 'exact_rejection' not in ls
        val,w,V,lp,M,L=search(K,terms)
        u=V[:,0];z0=np.linalg.solve(L.T,u);z0=z0/(4*np.linalg.norm(z0))
        z=[F(round(float(t)*10**12),10**12) for t in z0]
        assert 0<sum(t*t for t in z)<F(1,10)
        raw,eps=input_for(bank,E,k,J,eye,z)
        f=out.parent/f'input_{k}_adaptive.txt';f.write_text(raw)
        ans={'rank':n-1,'labels':ls['labels'],'selection_kind':'dual_bottleneck',
             'numeric_dual_margin':val,'K_lower_energy':str(quad(z,K)),
             'adaptive_source':{'coordinates':J,'coefficients':list(map(str,z)),
                  'input_sha256':sha(f),'error_upper':str(eps)}}
        report[k]=ans
        print(k,ans['selection_kind'],sha(f),flush=True)
    rank=sum(v['rank'] for v in report.values())
    out.write_text(json.dumps({'status':'CONDITIONAL_EVEN_BOTTLENECK_SELECTION','rank':rank,'conditional_index_upper':85-rank,
     'sectors':report,'premises':'Odd exact weights replayed without numerical search. Only even bottleneck source selected numerically; no additional rank claim.'},indent=2)+'\n')
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
