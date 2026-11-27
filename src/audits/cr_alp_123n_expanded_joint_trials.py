#!/usr/bin/env python3
"""Larger joint trial proposal bank: only outward interval acceptance.
mpmath supplies candidates, never proof or an inverse error enclosure.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys,inspect
import mpmath as mp
import cr_alp_123g_joint_trial_bank as G
SRC=Path(__file__).with_name('cr_alp_123e_complete_bridge.py')
PIN='c92c1c9e341071bc43c5e40994a1111e0fdbeaa9'
G_PIN='2a40d6466cd54e8af04d62c6421fc53d078eb1bf'
D=10**90
def blob(raw):return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
def candidate(A,b):
    with mp.workdps(145):
        M=mp.matrix([[mp.mpf(v)/mp.mpf(10)**340 for v in row] for row in A])
        rhs=mp.matrix([mp.mpf(v)/mp.mpf(10)**340 for v in b])
        x=mp.lu_solve(M,rhs)
        return [F(int(mp.nint(v*D)),D) for v in x]
def main(root,previous,joint,out):
    raw=SRC.read_bytes();assert blob(raw)==PIN
    assert blob(Path(G.__file__).read_bytes())==G_PIN
    jraw=joint.read_bytes()
    assert hashlib.sha256(jraw).hexdigest()=='d079ddf237af10624928245972156367d677347746623ba8e14611ea4ac64e82'
    old=json.loads(jraw);assert old['status']=='CONDITIONAL_JOINT_AUXILIARY_TRIAL_BANK'
    # Retain the original physical trial extraction and interval residual code.
    ws=inspect.getsource(G.witness)
    anchor='    x=[F(round(v*D),D) for v in solve(P,rhs)]+[F(1)]'
    assert ws.count(anchor)==1
    ws=ws.replace(anchor,'    x=candidate(P,rhs)+[F(1)]')
    assert ws.count('    assert en.l>0')==1
    ws=ws.replace('    assert en.l>0','    # A negative upper energy is a comparator witness only.')
    env=dict(vars(G),candidate=candidate)
    exec(compile(ws,'123n-witness','exec'),env)
    witness=env['witness']
    src=raw.decode()
    def replace(a,b):
        nonlocal src
        assert src.count(a)==1,a
        src=src.replace(a,b)
    replace("    src=src.replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=165')",
"""    src=src.replace('M=90;SHIFT=512;N=85','M=90;SHIFT=512;N=365')
    assert src.count('li=atan(y/xx)')==1
    assert src.count('assert y.up()/F(xx.l,S)<F(1,3)')==1
    src=src.replace('assert y.up()/F(xx.l,S)<F(1,3)','assert y.up()/F(xx.l,S)<F(3,5)')
    src=src.replace('li=atan(y/xx)','li=atan(y/xx) if (y/xx).up()<=F(1,3) else atan(I(F(1,4)))+atan((y/xx-F(1,4))/(1+y/xx/4))')
""")
    replace('        for N in (125,165):','        for N in (165,245,365):')
    replace('            for t in (F(0),F(1)):','            for t in (F(0),):')
    start=src.index('                lo=schur(endpoints(A,-1,S));')
    end=src.index('        result[str(k)]',start)
    src=src[:start]+"""                w=witness(A,base,C,keep,nl,high,low,c,st,B,k,S,I,ZERO,root)
                lo,hi=map(F,w['joint_trial_energy'])
                assert lo<=hi
                if N==165:
                    p=next(v for v in frozen_joint['sectors'][str(k)]['rows'] if v['cutoff']==165 and v['theta']=='0')
                    l,h=map(F,p['schur_interval'])
                    assert hi>=l and hi<=h+F(1,10**100)
                label='NEGATIVE_COMPARATOR_TRIAL' if hi<0 else 'POSITIVE_TRIAL_ENERGY_ONLY' if lo>0 else 'UNDECIDED_TRIAL_ENERGY'
                row={'cutoff':N,'theta':'0','dimension':len(J),'classification':label,'joint_trial':w}
                rows.append(row)
                print('candidate parity',k,'cutoff',N,'energy',float(lo),float(hi),label,flush=True)
"""+src[end:]
    replace("'CONDITIONAL_COMPLETE_BRIDGE_COMPARATOR_CENSUS'","'CONDITIONAL_EXPANDED_JOINT_AUXILIARY_TRIAL_BANK'")
    replace("'Finite projected subspaces of actual Q85 including complete infinite pole bridge, exact capped diagonal tail by Parseval. No infinite positivity claim. Original continuum and interval premises inherited.'","'Exact interval energy and finite residuals of rounded mpmath candidate trials. Not finite Schur lower bounds. Negative trial energy in lower comparator is not an original Weil witness.'")
    g={'__file__':str(SRC),'__name__':'cr_alp_123n_expanded','witness':witness,'frozen_joint':old}
    exec(compile(src,str(SRC),'exec'),g)
    g['main'](root,previous,out)
if __name__=='__main__':main(*map(Path,sys.argv[1:]))
