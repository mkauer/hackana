#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include "../../utils/h10.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TMath.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"

// use hcalc.C to solve Helene formula and set limit on signal.

Double_t HelenLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t lx, Int_t ux){
  //init histogram errors
  Float_t rintegral[400];
  
  if(nlim !=1){
    std::cout<<"Helen formula can put limit on one component only. But "<<nlim<<" components were provided"<<endl;
    std::cout<<"No limit was set"<<endl;
    return 0;
  }
  
  Int_t n0=(Int_t)data->Geth1(name)->Integral(lx,ux);
  Double_t mb=0;
  Double_t emb=0;
  Int_t i1;
  Double_t nev;
  
  for(Int_t i=0;i<nbgr;i++){
    nev=bgr[i]->Geth1(name)->Integral(lx,ux);
    mb+=bgr[i]->Geth1(name)->Integral(lx,ux)*bgr[i]->norm;
    emb+=pow(bgr[i]->GetErrorExpected(name,lx,ux),2);
  }
  for(Int_t i=0;i<nsig;i++){
    nev=sig[i]->Geth1(name)->Integral(lx,ux);
    mb+=sig[i]->Geth1(name)->Integral(lx,ux)*sig[i]->norm;
    emb+=pow(sig[i]->GetErrorExpected(name,lx,ux),2);
  }
  emb=sqrt(emb);
  
  Float_t CL=0.9;
  Float_t limit=SolveHelen(CL,mb,n0);
  Float_t enorm=limit/lim[0]->Geth1(name)->Integral(lx,ux);
  lim[0]->setnormalization(0,enorm);
  Double_t eff=lim[0]->Geth1(name)->Integral(lx,ux)/lim[0]->totaltime;
  
  
  cout<<"Helene Limit Results: "<<endl;
  cout<<"\t Energy window  = "<<fixed<<setprecision(2)<<data->Geth1(name)->GetBinLowEdge(lx)<<" - "
	<<data->Geth1(name)->GetBinLowEdge(ux)+data->Geth1(name)->GetBinWidth(ux)<<" MeV "<<endl;
  cout<<"\t Data Events    = "<<fixed<<setprecision(0)<<n0<<endl;
  cout<<"\t Exp Bkg Events = "<<fixed<<setprecision(2)<<mb<<" +/- "<<emb<<endl;
  cout<<"\t Exp 0v Events  < "<<fixed<<setprecision(2)<<limit<<" @ "<<setprecision(0)<<CL*100<<" % CL "<<endl;
  cout<<"\t 0v Efficiency  = "<<fixed<<setprecision(2)<<eff*100<<" % "<<endl;
  
  anaout<<"Energy window  = "<<fixed<<setprecision(2)<<data->Geth1(name)->GetBinLowEdge(lx)<<" - "
	<<data->Geth1(name)->GetBinLowEdge(ux)+data->Geth1(name)->GetBinWidth(ux)<<" MeV "<<endl;
  anaout<<"Data Events    = "<<fixed<<setprecision(0)<<n0<<endl;
  anaout<<"Exp Bkg Events = "<<fixed<<setprecision(2)<<mb<<" +/- "<<emb<<endl;
  anaout<<"Exp 0v Events  < "<<fixed<<setprecision(2)<<limit<<" @ "<<setprecision(0)<<CL*100<<" % CL "<<endl;
  anaout<<"0v Efficiency  = "<<fixed<<setprecision(2)<<eff*100<<" % "<<endl;
  
  
  return eff;
}
Double_t HelenLimitf(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Float_t lx, Float_t ux){
  
  Int_t ilx,iux;
  TH1D * h = data->Geth1(name);
  Int_t nbins=h->GetNbinsX();
  ilx=0;
  iux=nbins;
  
  for(Int_t i=1;i<nbins;i++){
    if(h->GetBinLowEdge(i)<lx) ilx=i;
    if(h->GetBinLowEdge(i)+h->GetBinWidth(i)<ux) iux=i;    
  }
  
  return HelenLimit(name,data,bgr,nbgr,sig,nsig,lim,nlim,ilx,iux);
}

