// Global fifth projected signed-prime moment; first reproduce second moment; includes the frozen 058 evaluator.
// Build against the extracted, SHA-verified original source; no replacement coefficients.
#define main cr_alp_058_original_main
#include "validated.cpp"
#undef main
struct Move59{int j,s,p;};
struct Panel59{I m,h;vector<Move59> active;};
I sq59(I x){Q z=x.absmax(),l=x.l<=0&&x.h>=0?0:min(fabsq(x.l),fabsq(x.h));return I(l==0?0:down(l*l),up(z*z));}
int main(int argc,char**argv){
 assert((argc==5||argc==7)&&fegetround()==FE_TONEAREST);
 mpfr_init2(mm,192);mpfr_init2(nn,192);mpfr_init2(oo,192);
 init(argv[1],argv[2],argv[3]);assert(k==2);
 // Same first Q85 projection as CR-ALP-058. Reintegrate, not a fitted finite projection.
 vector<I> proj(nl);I vb;
 for(auto &p:mesh[0])for(int j=0;j<96;j++){
  I t=p.m+p.h*gn[j],w=I(2)*p.h*gw[j];vector<I>b(nl);I x=eval(t,p.p,&b);
  for(int i=0;i<nl;i++)proj[i]=proj[i]+w*x*b[i];vb=vb+w*x*bval(t);
 }
 vb=widen(vb,1e-12Q);I st=(ex(a)-ex(-a))/I(2)+I(sgn)*a;vector<I> c(nl);
 for(int i=0;i<nl;i++){
  proj[i]=widen(proj[i],1e-12Q);I v=nu*I(ids[i]);
  c[i]=I(2)*(bp+bm)*v/(rt*(v*v+I(1)/I(4)));st=st-sq59(c[i]);vb=vb-c[i]*proj[i];
 }
 assert(st.l>0.003Q);I lam=vb/st;
 for(auto &p:pieces){for(int i=0;i<nl;i++)p.s[i]=p.s[i]+(c[i]*lam-proj[i])/rt;p.ep=p.ep-lam*bp;p.en=p.en-lam*bm;}
 midpointize();

 // N=2048. Original source coefficient projection exactly as in 065.
 constexpr int N=2048,nN=N/2;
 vector<I> poles(nN),rcoef(nN),ycoef(nN);
 for(int i=0;i<nN;i++){I v=nu*I(2*(i+1));poles[i]=I(2)*(bp+bm)*v/(rt*(v*v+I(1)/I(4)));}
 // At every 16th index restart with MPFR; between restarts use interval rotations.
 auto basis=[&](I t,vector<I>& bs,vector<I>& bc){
   auto step=sincosi(I(2)*nu*t);I sn,cs;
   for(int i=0;i<nN;i++){
     if(i%16==0){auto v=sincosi(I(2*(i+1))*nu*t);sn=v.first;cs=v.second;}
     else {I snext=sn*step.second+cs*step.first;cs=cs*step.second-sn*step.first;sn=snext;}
     bs[i]=sn/rt;bc[i]=cs/rt;
   }
 };
 ifstream mf(argv[4]);int np;mf>>np;assert(mf&&np>0);
 vector<Panel59> panels;vector<int> piece;
 for(int p=0;p<np;p++){Panel59 z;z.m=rd(mf);z.h=rd(mf);int ip,n;mf>>ip>>n;assert(mf&&ip>=0&&ip<(int)pieces.size()&&n>=0&&n<=14);piece.push_back(ip);for(int i=0;i<n;i++){Move59 m;mf>>m.j>>m.s;assert(mf&&m.j>=0&&m.j<7&&(m.s==-1||m.s==1));z.active.push_back(m);}panels.push_back(z);}
 string extra;assert(!(mf>>extra));
 vector<I> bs(nN),bc(nN);
 if(argc==5){
   I G;for(int p=0;p<np;p++){auto &z=panels[p];for(int j=0;j<96;j++){I t=z.m+z.h*gn[j],w=I(2)*z.h*gw[j],r=eval(t,piece[p]);G=G+w*sq59(r);basis(t,bs,bc);for(int i=0;i<nN;i++)rcoef[i]=rcoef[i]+w*r*bs[i];}if((p+1)%100==0)cerr<<"coeff panels "<<p+1<<"/"<<np<<"\n";}
   assert(G.l>0&&G.h<1);out("G",widen(G,1e-12Q));for(int i=0;i<nN;i++){assert(rcoef[i].absmax()<2);out("r"+to_string(2*(i+1)),widen(rcoef[i],1e-12Q));}
 }else{
   // x is exact orthogonal projection of rational trial coordinates onto the finite pole-null sine band.
   ifstream xf(argv[5]),df(argv[6]);vector<I> x(nN),diag(nN),dx(nN);I cp,pp0;
   for(int i=nl;i<nN;i++){int j;xf>>j;assert(xf&&j==2*(i+1));x[i]=rd(xf);df>>j;assert(df&&j==2*(i+1));diag[i]=rd(df);assert(diag[i].l>0);cp=cp+x[i]*poles[i];pp0=pp0+sq59(poles[i]);}
   assert(!(xf>>extra)&&!(df>>extra)&&pp0.l>0.003Q);
   I correction=cp/pp0,sx,sdx,X2;
   for(int i=nl;i<nN;i++){x[i]=x[i]-correction*poles[i];dx[i]=diag[i]*x[i];sx=sx+I(x[i].absmax());sdx=sdx+I(dx[i].absmax());X2=X2+sq59(x[i]);}
   assert(sx.h<100&&sdx.h<100&&X2.h<100);
   out("X2",X2);
   // Precompute all translation coefficients, retaining exact enclosed weights and rational shifts.
   vector<vector<I>> ps(7,vector<I>(nN)),pc(7,vector<I>(nN));
   for(int z=0;z<7;z++)for(int i=nl;i<nN;i++){auto sc=sincosi(nu*I(2*(i+1))*sh[z]);ps[z][i]=-wt[z]*x[i]*sc.second;pc[z][i]=-wt[z]*x[i]*sc.first;}
   I G,E,Y2,Yb;
   for(int p=0;p<np;p++){
     auto &z=panels[p];vector<I> vs=dx,vc(nN);
     for(auto m:z.active)for(int i=nl;i<nN;i++){vs[i]=vs[i]+ps[m.j][i];vc[i]=vc[i]+I(m.s)*pc[m.j][i];}
     for(int j=0;j<96;j++){
       I t=z.m+z.h*gn[j],w=I(2)*z.h*gw[j],r=eval(t,piece[p]),xx,v;
       basis(t,bs,bc);
       for(int i=nl;i<nN;i++){xx=xx+x[i]*bs[i];v=v+vs[i]*bs[i]+vc[i]*bc[i];}
       I y=r-v;
       G=G+w*sq59(r);E=E+w*xx*(I(2)*r-v);Y2=Y2+w*sq59(y);Yb=Yb+w*y*bval(t);
       for(int i=0;i<nN;i++)ycoef[i]=ycoef[i]+w*y*bs[i];
     }
     if((p+1)%100==0)cerr<<"response panels "<<p+1<<"/"<<np<<"\n";
   }
   assert(G.l>0&&G.h<1&&E.absmax()<100&&Y2.h<100&&Yb.absmax()<100);out("G",widen(G,1e-12Q));out("E",widen(E,1e-12Q));out("Y2",widen(Y2,1e-12Q));out("Yb",widen(Yb,1e-12Q));out("pole_st",st);
   for(int i=0;i<nN;i++){assert(ycoef[i].absmax()<100&&poles[i].absmax()<10);out("y"+to_string(2*(i+1)),widen(ycoef[i],1e-12Q));out("pole"+to_string(2*(i+1)),poles[i]);}
 }
 mpfr_clear(mm);mpfr_clear(nn);mpfr_clear(oo);
}
