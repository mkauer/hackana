#include "dollimit.hpp"

dollimit::dollimit(){
  initialised=false;
}

void dollimit::Init(TH1 *h2, TH1 *h3,Int_t lx, Int_t rx){
  bh=h2;
  sh=h3;
  signal=0;
 
  for(Int_t i=lx;i<rx+1;i++){
    s.push_back(sh->GetBinContent(i));
    b.push_back(bh->GetBinContent(i));
    signal+=sh->GetBinContent(i);
  }
  initialised=true;

}
void dollimit::fit(){
  if(!initialised){
    std::cout<<"DOLFIT ERROR. Histograms not initialised before the fit() call."<<std::endl;
    return;
  }
  rightside=2.71;// 90% CL exclusion
  double xup=10000;

  RF.SetFunction(Eval2,this, 1E-10,xup);
  RF.Solve(100,1E-8,1E-8);
  N= RF.Root();
}

double Eval2(double x, void *param){
  dollimit *p = (dollimit *)param;
  double value=0;
  if(x<=0) return -p->rightside;
  for(size_t i=0;i<p->b.size();i++){
   if(p->b[i]>0 && p->s[i]>0) value+=2*(p->b[i]*TMath::Log(p->b[i]/(p->b[i]+x*p->s[i])) + x*p->s[i]);
   if(p->b[i]==0 ) value += 2 *x* p->s[i];
   //cout<<"Eval i="<<i<<" b="<<p->b[i]<<" s="<<x*p->s[i]<<" value="<<value<<endl;
  }
  value=value - p->rightside;
  //cout<<"Eval call x="<<x<<" value = "<<value<<" rightside="<<p->rightside<<endl;
  return value;
}

