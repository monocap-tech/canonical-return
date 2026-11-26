// High-precision physical-space witness integration diagnostic.
// Gauss nodes/weights come from prepare.py's exact brackets; elementary evaluations
// and this mesh are not interval-certified. Do not use numerical sign for acceptance.
#include <quadmath.h>
#include <bits/stdc++.h>
using namespace std;using Q=__float128;
const Q a=1.1Q,d=2.2Q,gamma0=2.88Q; Q pi=acosq(-1),rt=sqrtq(a),nu=pi/d;
Q rd(istream &f){string a,b;f>>a>>b; if(!f)throw runtime_error("input");return (strtoflt128(a.c_str(),0)+strtoflt128(b.c_str(),0))/2;}
string ss(Q x){char buf[128];quadmath_snprintf(buf,sizeof(buf),"%.25Qe",x);return buf;}
struct Group{Q c;vector<Q>A;};
struct Piece{vector<Q>s,c; Q ep=0,en=0;};
int k,sgn;vector<int> ids;vector<Q> src,tr;Q f,db,bp,bm;vector<Q>lo;vector<Group> gs;vector<Q>sh,wt,knots,ga,gw;vector<Piece>pieces;vector<Q>projv;Q lam=0;long long evals=0;
Q bval(Q t){return bp*expq(t/2)+bm*expq(-t/2);}
Q group(Q t,const Group&g){Q x=g.c*t;if(x>1000)return 0; Q v=expq(-x),sum=0;for(int j=0;j<(int)g.A.size();j++){if(j)v*=x/j;sum+=v*g.A[j];}return sum;}
Q boundary(Q t){Q v=expq(-t/2),w=expq(-(d-t)/2),u=expq(-2*t),z=expq(-2*(d-t));Q val=0;for(auto A:lo){val+=A*(v+sgn*w);v*=u;w*=z;}for(auto &g:gs)val+=group(t,g)+sgn*group(d-t,g);
 Q et=expq(t/2),em=1/et;
 Q kr=bp*(et*(1+d-t)-em)+bm*(em*(t+1)-et*expq(-d));return val-f*kr;}
int locate(Q t){auto it=upper_bound(knots.begin(),knots.end(),t);return max(0,min((int)pieces.size()-1,int(it-knots.begin())-1));}
// recurrence evaluated in quad precision; restart every 16 modes.
Q eval(Q t,int p,vector<Q>*basis=nullptr){evals++;auto &v=pieces[p];Q sum=boundary(t)+v.ep*expq(t/2)+v.en*expq(-t/2);Q ca=cosq(2*nu*t),sa=sinq(2*nu*t),c=0,s=0;
 for(int j=0;j<(int)ids.size();j++){
  if(j%16==0){s=sinq(nu*ids[j]*t);c=cosq(nu*ids[j]*t);}else{Q next=s*ca+c*sa;c=c*ca-s*sa;s=next;}
  sum+=v.s[j]*s+v.c[j]*c;if(basis && ids[j]<=85)(*basis)[j]=s/rt;
 }return sum;}
void uniqueq(vector<Q>&v){sort(v.begin(),v.end());vector<Q>w;for(auto x:v)if(w.empty()||fabsq(x-w.back())>1e-28Q)w.push_back(x);v=w;}
void init(string path){ifstream in(path);int n;in>>k;sgn=k==1?1:-1;f=rd(in);db=rd(in);in>>n;for(int i=0;i<n;i++){int j;in>>j;ids.push_back(j);src.push_back(rd(in));tr.push_back(rd(in));}in>>n;for(int i=0;i<n;i++){Q beta=rd(in);(void)beta;lo.push_back(rd(in));}in>>n;for(int i=0;i<n;i++){Group g;g.c=rd(in);int r;in>>r;for(int j=0;j<r;j++)g.A.push_back(rd(in));gs.push_back(g);}
 bp=expq(-a/2)/2;bm=sgn*expq(a/2)/2;
 int nn[]={2,3,4,5,7,8,9},pp[]={2,3,2,5,7,2,3};knots={0,d};for(int j=0;j<7;j++){sh.push_back(logq(nn[j]));wt.push_back(logq(pp[j])/sqrtq(nn[j]));knots.push_back(sh.back());knots.push_back(d-sh.back());}uniqueq(knots);
 vector<vector<Q>> cs(7),sn(7);for(int j=0;j<7;j++)for(int i:ids){cs[j].push_back(cosq(nu*i*sh[j]));sn[j].push_back(sinq(nu*i*sh[j]));}
 for(int m=0;m+1<(int)knots.size();m++){Q mid=(knots[m]+knots[m+1])/2;Piece p;for(int i=0;i<(int)ids.size();i++){p.s.push_back(-gamma0*tr[i]/rt);p.c.push_back(0);}p.ep=-gamma0*db*bp;p.en=-gamma0*db*bm;
  for(int j=0;j<7;j++)for(int sign:{-1,1}){Q shift=sign*sh[j];if(mid+shift>0&&mid+shift<d){for(int i=0;i<(int)ids.size();i++){Q v=-wt[j]*(src[i]-tr[i])/rt;p.s[i]+=v*cs[j][i];p.c[i]+=v*sign*sn[j][i];}p.ep-=wt[j]*(f-db)*bp*expq(shift/2);p.en-=wt[j]*(f-db)*bm*expq(-shift/2);}}
  pieces.push_back(p);
 }
}
// recursively bisect a cell until the complex disk of radius 4h is controlled
// at both evaluation locations; small endpoint cells use cmax*4h <=1/2.
void panel(Q l,Q r,Q shift,vector<pair<Q,Q>>&out){Q h=(r-l)/2,m=(l+r)/2,cmax=gs.back().c;bool ok=h<=1.0Q/128;
 for(Q t:{m,m+shift})if(min(t,d-t)<8*h && 4*h*cmax>0.5Q)ok=false;
 if(ok){out.emplace_back(m,h);return;} if(m<=l||m>=r)throw runtime_error("mesh");panel(l,m,shift,out);panel(m,r,shift,out);}
vector<pair<Q,Q>> mesh(Q shift){vector<Q>c={0,d-shift};for(auto t:knots)for(Q x:{t,t-shift})if(x>1e-28Q&&x<d-shift-1e-28Q)c.push_back(x);uniqueq(c);vector<pair<Q,Q>>out;for(int i=0;i+1<(int)c.size();i++)panel(c[i],c[i+1],shift,out);return out;}
int main(int argc,char**argv){if(argc!=3)return 1;init(argv[1]);ifstream in(argv[2]);int N;in>>N;for(int i=0;i<N;i++){ga.push_back(rd(in));gw.push_back(rd(in));}int nl=0;for(auto j:ids)if(j<=85)nl++;projv.assign(nl,0);Q vb=0;
 auto ps=mesh(0);cerr<<"projection panels "<<ps.size()<<"\n";for(auto [m,h]:ps){int p=locate(m);for(int j=0;j<N;j++){Q t=m+h*ga[j],w=h*gw[j];vector<Q>bas(nl);Q u=eval(t,p,&bas);for(int i=0;i<nl;i++)projv[i]+=w*u*bas[i];vb+=w*u*bval(t);}}
 Q bm2=(expq(a)-expq(-a))/2+sgn*a,ss2=bm2;vector<Q>c(nl);for(int i=0;i<nl;i++){Q v=nu*ids[i];c[i]=2*(bp+bm)*v/(rt*(v*v+0.25Q));ss2-=c[i]*c[i];vb-=c[i]*projv[i];}lam=vb/ss2;
 cerr<<"projection pole denominator "<<ss(ss2)<<" lambda "<<ss(lam)<<"\n";
 for(auto &p:pieces){for(int i=0;i<nl;i++)p.s[i]+=(c[i]*lam-projv[i])/rt;p.ep-=lam*bp;p.en-=lam*bm;}
 Q G=0;for(auto [m,h]:ps){int p=locate(m);for(int j=0;j<N;j++){Q t=m+h*ga[j],w=h*gw[j];Q u=eval(t,p);G+=w*u*u;}}
 cout<<"G "<<ss(G)<<endl;Q M=0;for(int j=0;j<7;j++){auto q=mesh(sh[j]);Q C=0;for(auto [m,h]:q){int p=locate(m),r=locate(m+sh[j]);for(int i=0;i<N;i++){Q t=m+h*ga[i],w=h*gw[i];C+=w*eval(t,p)*eval(t+sh[j],r);}}M-=2*wt[j]*C;cout<<"corr "<<j<<" "<<ss(C)<<" panels "<<q.size()<<endl;}
 Q ell=51.0Q/200,u=1101.0Q/200;cout<<"M "<<ss(M)<<"\npenalty "<<ss((gamma0*G-M)/(ell*u))<<"\nevals "<<evals<<endl;}
