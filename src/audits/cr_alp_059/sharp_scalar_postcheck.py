"""CR-ALP-059: optimal scalar three-moment limitation and cubic targets.
Inputs: exact 058 and 059 endpoint JSONs, in that order. No new integral.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,sys
old,new=(json.loads(Path(p).read_text()) for p in sys.argv[1:3])
ell,u,g,s=map(F,('0.255','5.505','2.88','2.625'))
assert ell==g-s and u==g+s
result={}
for k in ('1','2'):
    G,M,K=(tuple(map(F,old[k][x])) for x in ('true_G','true_M','K'))
    N=tuple(map(F,new[k]['true_N']))
    assert K==tuple(map(F,new[k]['K'])) and new[k]['status']=='UNDECIDED'
    gl,gh=G; ml,mh=M; kl,kh=K; nl,nh=N
    B=(s*gl+ml,s*gh+mh); D=(s*gl-mh,s*gh-ml)
    A=(g*s*gl+u*ml+nl,g*s*gh+u*mh+nh)
    C=(g*s*gl-ell*mh-nh,g*s*gh-ell*ml-nl)
    assert all(0<x[0]<=x[1] for x in (G,A,B,C,D))
    U=(gl/ell-B[1]**2/(ell*A[0]),gh/ell-B[0]**2/(ell*A[1]))
    V=(gl/u+D[0]**2/(u*C[1]),gh/u+D[1]**2/(u*C[0]))
    # Every admissible exact moment triple has extremizing measures on
    # {ell,A/B} and {C/D,u}. Both sides of K remain possible.
    assert V[1]<kl<=kh<U[0]
    a,alpha,beta,Tpass,Tfail=map(F,('2.517','.671','-.589','.073','-.568') if k=='1' else ('2.134','.542','-.310','.269','-.383'))
    # Cubic reciprocal upper bound: H(t)=(t-ell)(u-t).
    # p3(t)=secant(t)+H(t)*(t/a**2-2/a)/(ell*u).
    # T=<r,P_Q^3 r>.  Upper penalty < K when T exceeds this threshold.
    def pass_threshold(gg,mm,kk,nn):
        h=s*s*gg-nn; f=g*gg-mm-ell*u*kk
        return g*h+s*s*mm+f*a*a-2*h*a
    corners=list(product(G,M,K,N))
    p=max(pass_threshold(*x) for x in corners)
    # Lower inverse comparison via the variational trial alpha*r+beta*P_Q*r.
    def fail_threshold(gg,mm,kk,nn):
        return (2*alpha*gg+2*beta*mm-alpha*alpha*(g*gg+mm)-2*alpha*beta*(g*mm+nn)-beta*beta*g*nn-kk)/(beta*beta)
    q=min(fail_threshold(*x) for x in corners)
    assert Tpass>p and Tfail<q
    result[k]={'optimal_upper_enclosure':list(map(str,U)),'optimal_lower_enclosure':list(map(str,V)),'cubic_a':str(a),'variational_alpha':str(alpha),'variational_beta':str(beta),'third_moment_accept_threshold':str(Tpass),'third_moment_reject_threshold':str(Tfail),'worst_accept_threshold':str(p),'worst_reject_threshold':str(q),'status':'THREE_SCALAR_MOMENTS_INSUFFICIENT; THIRD_MOMENT_NOT_EVALUATED'}
    print(k,'sharp lower',*[float(x) for x in V],'sharp upper',*[float(x) for x in U],'K',*[float(x) for x in K])
    print('third moment sufficient accept >=',float(Tpass),'reject <=',float(Tfail),'exact extrema',float(p),float(q))
print('PASS: optimal three-moment bounds straddle K for both witnesses; T not evaluated.')
if len(sys.argv)>3: Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n')
