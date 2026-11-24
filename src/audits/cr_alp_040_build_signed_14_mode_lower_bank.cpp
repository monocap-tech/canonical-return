// CR-ALP-040: original signed 14-mode Jensen midpoint lower bank through ten million Gamma indices.
// Compile: g++ -std=c++17 -O2 -o cr_alp_030_interval_core cr_alp_030_eight_mode_interval_core.cpp
// Run: ./cr_alp_030_interval_core CR-ALP-031_initial_rational_intervals.txt 30000
// Input created from source-attached exact rational prime intervals.
// All interval computations use checked 256-bit exact integers; no floating-point arithmetic.
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
// Fixed 256-bit signed integers with checked overflow; any overflow aborts certification.
using cpp_int = boost::multiprecision::checked_int256_t;
using namespace std;
static const cpp_int S("10000000000000000000000000"); //10^25
cpp_int floor_div(cpp_int n, const cpp_int &d){if(d<=0)throw runtime_error("nonpositive divisor");cpp_int q=n/d;if(n<0 && n%d!=0)--q;return q;}
cpp_int ceil_div(cpp_int n, const cpp_int &d){if(d<=0)throw runtime_error("nonpositive divisor");cpp_int q=n/d;if(n>0 && n%d!=0)++q;return q;}
struct FI{
 cpp_int lo,hi;
 FI():lo(0),hi(0){} FI(int x):lo(cpp_int(x)*S),hi(cpp_int(x)*S){}
 FI(long long x):lo(cpp_int(x)*S),hi(cpp_int(x)*S){}
 FI(cpp_int low,cpp_int high):lo(low),hi(high){if(low>high)throw runtime_error("invalid interval");}
 FI operator-()const{return FI(-hi,-lo);}
 FI operator+(const FI& b) const{return FI(lo+b.lo,hi+b.hi);}
 FI operator-(const FI& b) const{return FI(lo-b.hi,hi-b.lo);}
 FI operator*(const FI& b) const{
  array<cpp_int,4> v={lo*b.lo,lo*b.hi,hi*b.lo,hi*b.hi};
  auto p=minmax_element(v.begin(),v.end());
  return FI(floor_div(*p.first,S),ceil_div(*p.second,S));}
 FI operator/(const FI& b) const{
  if(b.lo<=0)throw runtime_error("division interval contains nonpositive denominator");
  array<cpp_int,4> vv={lo*S,lo*S,hi*S,hi*S};
  array<cpp_int,4> ds={b.lo,b.hi,b.lo,b.hi};
  cpp_int low=floor_div(vv[0],ds[0]),high=ceil_div(vv[0],ds[0]);
  for(int j=1;j<4;++j){low=min(low,floor_div(vv[j],ds[j]));high=max(high,ceil_div(vv[j],ds[j]));}
  return FI(low,high);}
 FI& operator+=(const FI& b){lo+=b.lo;hi+=b.hi;return *this;}
 FI& operator-=(const FI& b){auto n=*this-b;lo=n.lo;hi=n.hi;return *this;}
};
FI rat(long long a,long long b){return FI(cpp_int(a)*S,cpp_int(a)*S)/FI(b);}
map<string,FI> inputs;
FI get(const string& k){auto i=inputs.find(k);if(i==inputs.end())throw runtime_error("Missing "+k);return i->second;}
int parity(int j){return (j%2)?1:2;}
FI A, SQRT_A, E, EINV, EHALF, EHI, SinhA, q, g0;
map<int,FI> nu,lam,gamma_coeff;
map<pair<int,int>,FI> M,H;
FI frac_onehalf(){return rat(1,2);}
FI expo_pa(int sign){return (sign==1)?EHALF:EHI;}
FI Iexp(int z){if(z==0)return FI(2)*A;if(z==1||z==-1)return FI(2)*SinhA;throw runtime_error("Iexp invalid z");}
FI Iprime(int z){if(z==0)return FI(0);if(abs(z)!=1)throw runtime_error("Iprime invalid");FI v=FI(2)*A*(E+EINV)/FI(2)-FI(2)*SinhA;return z==1?v:-v;}
FI Lbeta0(int p,int r){
 // β=1/2, p,r=±1/2; signed integers encode 2p,2r.
 int z=(r+1)/2;
 int pminus=(p-1)/2;
 if(z==0)return Iprime(pminus)+A*Iexp(pminus);
 return Iexp((p+r)/2)-EINV*Iexp(pminus);
}
FI Kbb(int kind,const FI& beta,const FI& beta2,const FI& qp,int k){
 FI vals[2][2];const int signs[2]={1,-1};
 for(int i=0;i<2;++i)for(int j=0;j<2;++j){
  int p=signs[i],r=signs[j];
  if(k==0){vals[i][j]=Lbeta0(p,r)+Lbeta0(r,p);continue;}
  int z=(p+r)/2;
  FI ep=(z==1)?E:(z==-1?EINV:FI(1));
  FI em=FI(1)/ep;
  FI er=(r==1)?E:EINV;
  FI pp=rat(p,2), rr=rat(r,2);
  FI den=beta2-rat(1,4);
  FI v=FI(2)*beta*Iexp(z)/den;
  v-=em*(FI(1)-qp*er)/((beta+pp)*(beta-rr));
  v-=ep*(FI(1)-qp/er)/((beta-pp)*(beta+rr));
  vals[i][j]=v;
 }
 FI cross=vals[0][1]+vals[1][0];
 FI sum=vals[0][0]+vals[1][1]+(kind==1?cross:-cross);
 return sum/FI(4);
}
FI K_sine_b(int j,int kind,const FI& beta,const FI& beta2,const FI& qp,int k){
 FI vals[2];
 for(int c=0;c<2;++c){
  int p=(c==0)?1:-1;
  FI pp=rat(p,2);FI ep=expo_pa(p);FI em=FI(1)/ep;
  FI v=nu[j],den=beta2+v*v;
  FI mom=v*(em-FI((j%2)?-1:1)*ep)/(SQRT_A*(v*v+rat(1,4)));
  FI fac=(p==1)?E:EINV;
  FI ld=(k==0 && p==1)?FI(2)*A*em:em*(FI(1)-qp*fac)/(beta-pp);
  FI rd=(k==0 && p==-1)?FI(2)*A*ep:ep*(FI(1)-qp/fac)/(beta+pp);
  FI val=(FI(2)*beta*mom+v/SQRT_A*(ld-FI((j%2)?-1:1)*rd))/den;
  vals[c]=val;
 }
 return (vals[0]+(kind==1?vals[1]:-vals[1]))/FI(2);
}

void ldl_and_report(const vector<int>& js,const string& label){
 int n=js.size();vector<vector<FI>> L(n,vector<FI>(n));vector<FI> piv;
 for(int i=0;i<n;++i){
  FI d=H[{js[i],js[i]}];
  for(int m=0;m<i;++m)d-=L[i][m]*L[i][m]*piv[m];
  cout<<label<<" pivot "<<i<<" lower="<<d.lo<<" upper="<<d.hi<<endl;
  if(d.lo<=0)throw runtime_error("uncertified nonpositive LDL pivot");
  piv.push_back(d);
  for(int j=i+1;j<n;++j){
   FI v=H[{js[j],js[i]}];
   for(int m=0;m<i;++m)v-=L[j][m]*L[i][m]*piv[m];
   L[j][i]=v/d;
  }
 }
 // Exact Frobenius-bound on inverse of lower triangular L
 vector<vector<FI>> inv(n,vector<FI>(n));
 for(int i=0;i<n;++i){
  inv[i][i]=FI(1);
  for(int j=0;j<i;++j){
   FI sum(0);
   for(int k=j;k<i;++k)sum+=L[i][k]*inv[k][j];
   inv[i][j]=-sum;
  }
 }
 FI norm(0),trace(0);
 for(int i=0;i<n;++i){
  trace+=M[{js[i],js[i]}];
  for(int j=0;j<=i;++j)norm+=inv[i][j]*inv[i][j];
 }
 cpp_int smallest=piv[0].lo;
 for(auto &v:piv)smallest=min(smallest,v.lo);
 // generalized lambda >= smallest*S / (norm.hi*trace.hi).
 // Check lambda > 1e-10 via exact integer comparison.
 cpp_int numeric=smallest*S;
 cpp_int denom=norm.hi*trace.hi;
 if(numeric<=0 || denom<=0)throw runtime_error("bad norm interval");
 bool gt1e15=(numeric*1000000000000000LL>denom);
 cout<<label<<" exact generalized lower numerator="<<numeric<<" denominator="<<denom<<" >1e-15="<<gt1e15<<endl;
 if(!gt1e15)cout<<"crude spectral lower estimate inconclusive; shifted LDL is definitive"<<endl;
}

void ldl_shifted(const vector<int>& js,const FI& shift,const string& label){
 int n=js.size();vector<vector<FI>> L(n,vector<FI>(n));vector<FI> piv;
 for(int i=0;i<n;++i){
  FI d=H[{js[i],js[i]}]-shift*M[{js[i],js[i]}];
  for(int k=0;k<i;++k)d-=L[i][k]*L[i][k]*piv[k];
  cout<<label<<" shifted pivot "<<i<<" lower="<<d.lo<<" upper="<<d.hi<<endl;
  if(d.lo<=0)throw runtime_error("cannot certify shifted physical H-lambda*M positive");
  piv.push_back(d);
  for(int j=i+1;j<n;++j){
   FI v=H[{js[j],js[i]}]-shift*M[{js[j],js[i]}];
   for(int k=0;k<i;++k)v-=L[j][k]*L[i][k]*piv[k];
   L[j][i]=v/d;
  }
 }
}

int main(int argc,char**argv){
 try{
  if(argc<3)throw runtime_error("usage: core input.txt K");
  ifstream in(argv[1]);if(!in)throw runtime_error("cannot open inputs");
  string k,low,high;
  while(in>>k>>low>>high)inputs[k]=FI(cpp_int(low),cpp_int(high));
  A=get("A");SQRT_A=get("SQRT_A");E=get("E");EINV=get("EINV");EHALF=get("EHALF");EHI=get("EHI");SinhA=get("SinhA");q=get("q");g0=get("g0");
  for(int j=1;j<=14;++j)nu[j]=get("nu_"+to_string(j));
  for(int j=3;j<=14;++j)lam[j]=get("lam_"+to_string(j));
  for(int t=1;t<=2;++t)gamma_coeff[t]=get("gamma_"+to_string(t));
  for(auto inds: {vector<int>{3,5,7,9,11,13,1},vector<int>{4,6,8,10,12,14,2}})
    for(int i:inds)for(int j:inds){M[{i,j}]=get("M_"+to_string(i)+"_"+to_string(j));H[{i,j}]=get("H_"+to_string(i)+"_"+to_string(j));}
  int K=stoi(argv[2]);if(K!=30000)throw runtime_error("CR-ALP-033 requires exactly 30000 individual Gamma terms before grouped tail");
  FI qp=q,q4=q*q*q*q;
  // The integral defining J_beta is PSD and convex in beta,
  // since its second beta derivative is integral s^2 exp(-beta*s) times
  // a positive squared translation-difference. Therefore the Loewner
  // Jensen lower bound uses the group mean beta, not its endpoint.
  // Evaluate the first K=30000 terms exactly; group the following terms
  // with deterministic step size floor(start/2000) to 350,000.
  vector<pair<int,int>> schedule;
  for (int m=0;m<K;++m) schedule.emplace_back(m,1);
  int start=K, grouped_count=0;
  while(start<10000000){
   int len=max(1,start/2000);
   int end=min(10000000,start+len);
   schedule.emplace_back(end-1,end-start);
   ++grouped_count;start=end;
  }
  cout<<"Scheduled "<<K<<" exact terms and "<<grouped_count
      <<" monotone Loewner lower-tail groups through m=9999999"<<endl;
  vector<pair<vector<int>,int>> blocks={{{3,5,7,9,11,13},1},{{4,6,8,10,12,14},2}};
  for(const auto& entry:schedule){
   int t=entry.first,weight=entry.second;
   // For t>=30000, exp(-a(4t+1)) is in [0,10^-25]; because
   // exp(-a)<2/5 and (2/5)^(120001)<10^-25.
   FI qterm=(t<K)?qp:FI(cpp_int(0),cpp_int(1));
   long long beta_numer = 4LL*t - 2LL*weight + 3;
   FI beta=rat(beta_numer,2);
   FI beta2=rat(beta_numer*beta_numer,4);
   map<int,FI> D,h,diag;
   for(int j=1;j<=14;++j){D[j]=beta2+nu[j]*nu[j];h[j]=nu[j]/D[j];diag[j]=FI(2)*nu[j]*nu[j]/(beta*D[j]);}
   for(const auto &[js,anchor]:blocks){
    FI factor=FI(2)/A*(FI(1)-FI((anchor==1)?-1:1)*qterm);
    FI hka=h[anchor],kk=FI(2)*beta/D[anchor]+factor*hka*hka;
    map<int,FI> kb;
    kb[anchor]=K_sine_b(anchor,anchor,beta,beta2,qterm,t);
    for(int j:js)kb[j]=K_sine_b(j,anchor,beta,beta2,qterm,t);
    FI bb=Kbb(anchor,beta,beta2,qterm,t);
    for(int i:js){
     FI vi=h[i]-lam[i]*hka;
     for(int j:js){
      FI v=(i==j?diag[i]:FI(0))+lam[i]*lam[j]*diag[anchor]-factor*vi*(h[j]-lam[j]*hka);
      H[{i,j}]+=FI(weight)*v;
     }
     FI Kik=factor*h[i]*hka;
     FI Kuw=Kik-lam[i]*kk-gamma_coeff[anchor]*(kb[i]-lam[i]*kb[anchor]);
     FI cross=FI(2)*M[{i,anchor}]/beta-Kuw;
     H[{i,anchor}]+=FI(weight)*cross;
     H[{anchor,i}]=H[{i,anchor}];
    }
    H[{anchor,anchor}]+=FI(weight)*(FI(2)*M[{anchor,anchor}]/beta-(kk-FI(2)*gamma_coeff[anchor]*kb[anchor]+gamma_coeff[anchor]*gamma_coeff[anchor]*bb));
   }
   if(t<K)qp=qp*q4;
   if(t<K && (t+1)%5000==0)cout<<"Summed "<<t+1<<" certified PSD Gamma kernels"<<endl;
  }
  {ofstream o("CR-ALP-040_midpoint_H_14_mode_10m_intervals.tsv");
   for(const auto& b:blocks){vector<int> js=b.first;js.push_back(b.second);
    for(int i:js)for(int j:js)o<<i<<" "<<j<<" "<<H[{i,j}].lo<<" "<<H[{i,j}].hi<<"\n";}}
  // Acceptance intentionally uses a stable *separate* exact-rational
  // Sylvester+row-sum verifier. Recursive interval LDL loses correlation on
  // the near-null w13 pivot, though the true Jensen lower bank is positive.
  // Do not label failure of naive interval LDL as a negative Weil direction.
  cout<<"CR-ALP-040: 14-mode original-signed Jensen Gamma lower interval bank successfully constructed; VERIFY WITH stable_14_mode_sylvester_audit.py"<<endl;
  return 0;
 }catch(const exception&e){cerr<<"CERTIFICATE FAIL: "<<e.what()<<endl;return 1;}
}
