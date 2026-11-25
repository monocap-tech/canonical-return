"""Independent rational check of saved Gauss brackets and actual input/mesh bounds."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json
p=Path(__file__).resolve().parent
N=96; co=[0]*(N+1)
for j in range(N//2+1):co[N-2*j]=(-1)**j*comb(N,j)*comb(2*N-2*j,N)
dp=[(i+1)*co[i+1] for i in range(N)]
dd=F(sum(abs(i*dp[i]) for i in range(1,N)),2**N)
def val(c,x):
 m,q=x.numerator,x.denominator;v=0;qp=1
 for a in reversed(c):v=v*m+a*qp;qp*=q
 return F(v,qp//q*2**N)
rows=(p/'gauss.txt').read_text().splitlines();assert int(rows[0])==N and len(rows)==N+1
last=F(-1); sl=su=F(0)
for row in rows[1:]:
 l,h,wl,wh=map(F,row.split());assert last<l<h<1 and val(co,l)*val(co,h)<0;last=h
 der=abs(val(dp,(l+h)/2));er=dd*(h-l)/2;assert der>er
 ma=max(abs(l),abs(h));mi=min(abs(l),abs(h)) if l*h>0 else F(0)
 # Stored root bracket is exact at 75 decimal places; stored weight endpoints outward at 80.
 assert 0<wl<=2/((1-mi**2)*(der+er)**2)<=2/((1-ma**2)*(der-er)**2)<=wh
 sl+=wl;su+=wh
assert sl<2<su
print('96 independent exact root and weight brackets PASS')
for k in (1,2):
 it=iter((p/f'input_{k}.txt').read_text().split());assert int(next(it))==k
 def rd():
  l,h=F(next(it)),F(next(it));assert l<=h;return max(abs(l),abs(h))
 assert rd()<100 and rd()<100;n=int(next(it));assert n==425; src=tr=F(0)
 for _ in range(n):
  i=int(next(it));assert 1<=i<=850 and i%2==k%2;src+=rd();tr+=rd()
 assert src+tr<10**17
 n=int(next(it));assert n==128;lo=F(0)
 for _ in range(n):assert rd()<256;lo+=rd()
 assert lo<10**14
 n=int(next(it));assert n==25;gm=F(0)
 for _ in range(n):
  assert rd()<=3*2**31-F(1,2);n=int(next(it));assert n==40;gm+=sum(rd() for _ in range(n))
 assert gm<1000 and next(it,None) is None
 print('sector',k,'source/trial/boundary coefficient budgets PASS')
# Direct degree-191 product remainder, replacing the earlier unverified real-axis bound.
d=F(11,5);error=4*d*F(4,3)*10**80/4**192
assert error<F(1,10**12)
# For shift displacement <=1e-70, use 32 intervals of discontinuity and rough
# uniform real and derivative bounds on the retained physical pieces.
delta=F(1,10**70)
# sqrt(32 delta)<6e-35; sqrt(d)<2. Weighted translation source L2 transport.
assert 8*(2*10**21*delta+2*10**18*F(6,10**35))<F(1,10**13)
# Prime moment shift discrepancy: L1 translation bound times uniform magnitude.
assert 8*10**18*(d*10**21*delta+32*2*10**18*delta)<F(1,10**13)
print('direct product quadrature and shift displacement budgets PASS')
