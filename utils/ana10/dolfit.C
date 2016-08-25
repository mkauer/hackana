#include "dolfit.hpp"

dolfit::dolfit(){
  initialised=false;
}

void dolfit::Init(TH1* h1, TH1 *h2, TH1 *h3,Int_t lx, Int_t rx,Int_t rebin){
  dh=h1;
  bh=h2;
  sh=h3;
  signal=0;
  data=0;

  Double_t ss,bs,ds;
  ss=0;
  bs=0;
  ds=0;
  Int_t rs = 0;
  for(Int_t i=lx;i<rx+1;i++){
    ss += sh->GetBinContent(i);
    bs += bh->GetBinContent(i);
    ds += dh->GetBinContent(i);
    rs ++;
    if(rs == rebin){
      s.push_back(ss);
      b.push_back(bs);
      d.push_back(ds);
      signal+=ss;
      data+=ds;
      ss = 0;
      bs = 0;
      ds = 0;
      rs = 0;
    }
  }
  if(rs){
    s.push_back(ss);
    b.push_back(bs);
    d.push_back(ds);
    signal+=ss;
    data+=ds;
    ss = 0;
    bs = 0;
    ds = 0;
    rs = 0;
  }
  initialised=true;

}
void dolfit::fit(){
  if(!initialised){
    std::cout<<"DOLFIT ERROR. Histograms not initialised before the fit() call."<<std::endl;
    return;
  }
  rightside=signal;
  double xup=data/signal;

  RF.SetFunction(Eval,this, 1E-10,xup);
  RF.Solve(100,1E-8,1E-8);
  N= RF.Root();

  rightside=0;
  double log_max=LEval(N, this);
  rightside = log_max-0.5*1;
  double xlow = N-1E-10;
  xup = xlow + 30*sqrt(data)/signal;
  RF.SetFunction(LEval,this, xlow,xup);
  RF.Solve(100,1E-8,1E-8);
  NR=RF.Root();

  xup = N+1E-10;
  xlow = xup - 30*sqrt(data)/signal;
  if (xlow<0) xlow = 1E-10;
  RF.SetFunction(LEval,this, 1E-10,xup);
  RF.Solve(100,1E-8,1E-8);
  NL=RF.Root();
  
}
double Eval(double x, void *param){
  dolfit *p = (dolfit *)param;
  double value=0;
  for(size_t i=0;i<p->d.size();i++){
    if(p->d[i]!=0) value+=(p->d[i]*p->s[i])/(p->b[i]+x*p->s[i]);
    //    cout<<"Eval i="<<i<<" d="<<p->d[i]<<" b="<<p->b[i]<<" s="<<x*p->s[i]<<" value="<<value<<endl;
  }
  value=value - p->rightside;
  //cout<<"Eval call x="<<x<<" value = "<<value<<" rightside="<<p->rightside<<endl;
  return value;
}
double LEval(double x, void *param){
  dolfit *p = (dolfit *)param;
  double value=0;
  for(size_t i=0;i<p->d.size();i++){
    if(p->d[i]!=0) value+=-(x*p->s[i]) + p->d[i]*TMath::Log(p->b[i]+x*p->s[i]);
    //    cout<<"LEval i="<<i<<" d="<<p->d[i]<<" b="<<p->b[i]<<" s="<<x*p->s[i]<<" value="<<value<<endl;
  }
  value=value - p->rightside;
  //cout<<"LEval call x="<<x<<" value = "<<value<<" rightside="<<p->rightside<<endl;
  return value;
}

