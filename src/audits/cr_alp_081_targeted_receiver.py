#!/usr/bin/env python3
"""Selective affine inverse refinement, exact receiving tests, no rank assumption."""
from pathlib import Path
from fractions import Fraction as F
import sys,json
from cr_alp_074_balanced_selection import checked,HASHES,ELL
from cr_alp_075_shifted_energy import lower_matrix,positive,S,ceil
from cr_alp_078_moment_receiver import moments,sha,candidate
from cr_alp_079_inverse_gram_floor import lower
TARGET={'1':[13,10,11,14,15,7],'2':[11,13,16,24,19]}
def main(p,bank,oldpath,r70path,fullpath,oldlogs,newlogs,out):
    old=json.loads(oldpath.read_text());assert old['rank']==69
    r70=json.loads(r70path.read_text());assert r70['status']=='PASS_CONDITIONAL_RANK70' and F(r70['theta_upper'])<1
    assert F(old['sectors']['1']['levels'][-1]['margin'])>0
    full=json.loads(fullpath.read_text());assert full['status']=='CONDITIONAL_FULL85_DIAGNOSTIC'
    entries=json.loads((bank/'source_bank.json').read_text())['entries']
    E={(e['parity'],e['coordinate_zero_based']):e for e in entries}
    assert len(E)==len(entries)==85
    for e in entries:assert sha(bank/e['file'])==e['sha256'] and F(e['norm2'])==1
    report={};total=0
    for k,h in HASHES.items():
        ss=full['sectors'][k];J=ss['coordinate_order'];n=43 if k=='1' else 42
        assert sorted(J)==list(range(n))
        cols=old['sectors'][k]['columns'];assert [x['coordinate'] for x in cols]==J[:35]
        a={}
        for x in cols:
            j=x['coordinate'];assert x['input_sha256']==E[k,j]['sha256']
            a[j]=F(x['norm_upper'])**2/ELL
        if k=='2':
            for j,v in r70['refined_inverse_upper'].items():a[int(j)]=min(a[int(j)],F(v))
        assert sorted(map(int,ss['new_moments']))==sorted(J[35:])
        for j in J[35:]:
            q=moments(oldlogs/f'mode_{k}_{j}.log',E[k,j])
            assert q==ss['new_moments'][str(j)]
            a[j]=F(q['inverse_upper'])
        refined={}
        assert set(TARGET[k])<=set(J[:35])
        for j in TARGET[k]:
            row=moments(newlogs/f'mode_{k}_{j}.log',E[k,j])
            before=a[j];a[j]=min(before,F(row['inverse_upper']))
            refined[str(j)]={**row,'prior_inverse_upper':str(before),'retained_inverse_upper':str(a[j])}
        w=checked(p/'cr_alp_056'/f'result_extended_{k}.json',h)
        best=35;attempts={}
        for r in range(n,35,-1):
            I=J[:r];H=lower_matrix(w,I)
            v=candidate(H,[a[j] for j in I])
            # Retain alternate finite Gram comparison, only if direct test fails.
            vv=None
            if not v['pass_conditional']:
                H2,c=lower(w,I,a);vv=candidate(H2,c)
            attempts[str(r)]={'direct':v,'inverse_gram':vv}
            if v['pass_conditional'] or (vv and vv['pass_conditional']):
                best=r;break
        total+=best
        report[k]={'coordinate_order':J,'refinements':refined,'inverse_bounds':{str(j):str(v) for j,v in a.items()},'attempts':attempts,'best_conditional_rank':best}
    result={'status':'CONDITIONAL_TARGETED_MOMENT_RECEIVER','rank':total,'conditional_index_upper':85-total,'sectors':report,
      'premises':'Original authenticated fixed-collar source, norm and moment banks; inherited continuum and numerical enclosure premises; finite projection attachment for Gram alternative; exact integer pivots and rational sums. No unrestricted RH claim.'}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
