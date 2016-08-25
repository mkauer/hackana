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
#include "TCanvas.h"
#include "TVectorD.h"
#include "TObjString.h"

// Set 90% limit on signal with the help of ROOT  TLimit package

TConfidenceLevel * MCLimit_basic(TH1D *d,TH1D *bg,TH1D *s,Double_t &c,Double_t &xc,Double_t es,Double_t eb){
  // for verbose output
  Bool_t verbose=false;
  
  // init histogram errors 
  Float_t rintegral[400];


  const Double_t eps=0.1;
  const Double_t bgrmin=0.0;
  const Double_t cld=0.1;
  const Int_t nmc=10000;
  const Int_t nmc_max=20000;
  const Int_t nmc_min=20000;

  Bool_t stat = true;
  Double_t xa,xb,a,b;
  Int_t i;
  TH1D *snee;
  
  TVectorD errorb(2);
  TVectorD errors(2);
  TObjArray* names = new TObjArray();
  TObjString name1("bg uncertainty");
  TObjString name2("sig uncertainty");
  names->AddLast(&name1);
  names->AddLast(&name2);
  errorb[0]=eb;  // error source 1: 5%
  errorb[1]=0;   // error source 2: 0%
  errors[0]=0;   // error source 1: 0%
  errors[1]=es;  // error source 2: 1%

  a=0;
  snee=(TH1D*)s->Clone("snee");
  snee->Scale(a);
  TLimitDataSource* mydatasource = new TLimitDataSource(snee,bg,d);//&errors,&errorb,names);
  TConfidenceLevel* myconfidence = TLimit::ComputeLimit(mydatasource,1000);
  xa=myconfidence->CLs();
  if(verbose){
    cout<<  "Norm = "<<b<<" -- CLs = "<<myconfidence->CLs()<<" -- CLb = "<<myconfidence->CLb()<<" -- CLsb = "<<myconfidence->CLsb()<<endl;
    anaout<<"Norm = "<<b<<" -- CLs = "<<myconfidence->CLs()<<" -- CLb = "<<myconfidence->CLb()<<" -- CLsb = "<<myconfidence->CLsb()<<endl;
  }
  delete myconfidence;
  delete mydatasource;
  delete snee;  
  
  b=1000/s->Integral();
  snee=(TH1D*)s->Clone("snee");
  snee->Scale(b);
  mydatasource = new TLimitDataSource(snee,bg,d);//&errors,&errorb,names);
  myconfidence = TLimit::ComputeLimit(mydatasource,1000);
  xb=myconfidence->CLs();
  if(verbose){
    cout<<  "Norm = "<<b<<" -- CLs = "<<myconfidence->CLs()<<" -- CLb = "<<myconfidence->CLb()<<" -- CLsb = "<<myconfidence->CLsb()<<endl;
    anaout<<"Norm = "<<b<<" -- CLs = "<<myconfidence->CLs()<<" -- CLb = "<<myconfidence->CLb()<<" -- CLsb = "<<myconfidence->CLsb()<<endl;
  }
  if(1-myconfidence->CLb()<bgrmin){
    c=1;
    cout<<"Background only hypothesis excluded at the level <"<<bgrmin<<endl;
    cout<<"Set signal limit to inf."<<endl;
    delete snee;
    return 0;
  }
  while(fabs(a-b)*s->Integral()*10>eps){
    delete myconfidence;
    delete mydatasource;
    delete snee;
    if(a!=0){
      //c=a+(b-a)*((xa-cld)/(xa-xb));
      c=a+(b-a)/3.;
    }else{
      c=a+(b-a)/3.;
    }
    Int_t toymc=1000;
    toymc = Int_t(10.*1./(xa-cld)/(cld-xb));
    if(toymc>nmc_max) toymc=nmc_max;
    if(toymc<nmc_min) toymc=nmc_min;
    snee=(TH1D*)s->Clone("snee");
    snee->Scale(c);
    mydatasource = new TLimitDataSource(snee,bg,d);//&errors,&errorb,names);
    myconfidence = TLimit::ComputeLimit(mydatasource,toymc);
    xc=myconfidence->CLs();
    if(verbose){
      cout<<  "Norm = "<<b<<" -- CLs = "<<myconfidence->CLs()<<" -- CLb = "<<myconfidence->CLb()<<" -- CLsb = "<<myconfidence->CLsb()<<endl;
      anaout<<"Norm = "<<b<<" -- CLs = "<<myconfidence->CLs()<<" -- CLb = "<<myconfidence->CLb()<<" -- CLsb = "<<myconfidence->CLsb()<<endl;
    }
    if(xc>cld){
      a=c;
      xa=xc;
    }else{
      b=c;
      xb=xc;
    }
  }
  delete mydatasource;
  delete snee;
  return myconfidence;
}

Double_t MCLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebin,Double_t es,Double_t eb){
  //init histogram errors
  Float_t rintegral[400];
  
  if(nlim !=1){
    std::cout<<" MCLimit can put limit on one component only. But "<<nlim<<" components were provided"<<endl;
    std::cout<<" No limit was set"<<endl;
    return 0;
  }
  TH1D *d  = (TH1D *)data->Geth1(name)->Rebin(rebin,"tlimitdata$");
  TH1D *bg = (TH1D *)data->Geth1(name)->Clone("mcl$");
  TH1D *s  = (TH1D *)lim[0]->Geth1(name)->Rebin(rebin,"tlimitsig$");
  
  bg->Reset();
  for(Int_t i=0;i<nbgr;i++){
    bg->Add(bgr[i]->Geth1(name),bgr[i]->norm);
    //    cout<<"Adding TLimit bgr "<<bgr[i]->name<<" norm="<<bgr[i]->norm<<" integral="<<bgr[i]->Geth1(name)->Integral()<<endl;
  }
  for(Int_t i=0;i<nsig;i++){
    bg->Add(sig[i]->Geth1(name),sig[i]->norm);
  }
  bg->Rebin(rebin);
  //cout<<"Tlimit, total bgr "<<bg->Integral()<<endl;
  Double_t c,xc;
  
  
  //do not use under/overflow bins!
  Int_t nbinsx = d->GetNbinsX();
  d->SetBinContent(0,0);
  bg->SetBinContent(0,0);
  s->SetBinContent(0,0);
  d->SetBinContent(nbinsx+1,0);
  bg->SetBinContent(nbinsx+1,0);
  s->SetBinContent(nbinsx+1,0);
  TConfidenceLevel *result = MCLimit_basic(d,bg,s,c,xc,es,eb);
  
  Int_t i;
  Float_t eff;
  if(lim[0]->totaltime){
    eff = s->Integral()/lim[0]->totaltime;
  }else{
    eff = s->Integral();
  }
  
  lim[0]->setnormalization(0,c);
  
  cout<<"MC TLimit Results: "<<endl;
  cout<<"\t Data Events    = "<<fixed<<setprecision(0)<<d->Integral()<<endl;
  cout<<"\t Exp Bkg Events = "<<fixed<<setprecision(2)<<bg->Integral()<<endl;
  cout<<"\t Exp 0v Events  < "<<fixed<<setprecision(2)<<s->Integral()*c<<" @ "<<setprecision(0)<<100*(1-xc)<<" % CL "<<endl;
  cout<<"\t 0v Efficiency  = "<<fixed<<setprecision(2)<<eff*100<<" % "<<endl;
  
  anaout<<"Data Events    = "<<fixed<<setprecision(0)<<d->Integral()<<endl;
  anaout<<"Exp Bkg Events = "<<fixed<<setprecision(2)<<bg->Integral()<<endl;
  anaout<<"Exp 0v Events  < "<<fixed<<setprecision(2)<<s->Integral()*c<<" @ "<<setprecision(0)<<100*(1-xc)<<" % CL "<<endl;
  anaout<<"0v Efficiency  = "<<fixed<<setprecision(2)<<eff*100<<" % "<<endl;
  
  
  /**
  TCanvas *mp= new TCanvas("ttt","MCLimit",800,800);
  mp->Divide(2,2);
  mp->cd(1);
  d->DrawCopy();
  mp->cd(2);
  bg->DrawCopy();
  mp->cd(3);
  s->DrawCopy();
  mp->cd(4);
  result->Draw();
  **/
  
  delete d;
  delete bg;
  delete s;

  return eff;
}

Int_t MCLimitch(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[],bana10 * lim[][MAXCOMP], Int_t nlim[],Int_t nchnls,Int_t rebin,Double_t es, Double_t eb){
  //init histogram errors
  Float_t rintegral[400];

  if(nlim[0] !=1){
    std::cout<<" MCLimit can put limit on one component only. But "<<nlim<<" components were provided"<<endl;
    std::cout<<" No limit was set"<<endl;
    return 1;
  }

  Int_t nbint=0;
  // calculate total number of bins
  for(Int_t i=0;i<nchnls;i++){
    nbint+=data[i]->Geth1(name)->GetNbinsX();    
  }
  nbint=nbint/rebin;
  //book main histograms
  TH1D * d = new TH1D("tlimitdata","tlimidtat",nbint,0.,1.);
  TH1D * bg = new TH1D("mcl","mcl",nbint,0.,1.);
  TH1D * s = new TH1D("tlimitsig","tlimitsig",nbint,0.,1.);
  
  Float_t timet=0;
  for(Int_t i=0;i<nchnls;i++){  
    TH1D *d1 = (TH1D *)data[i]->Geth1(name)->Rebin(rebin,"tlimitdata1");
    TH1D *bg1 = (TH1D *)data[i]->Geth1(name)->Clone("mcl1");
    TH1D *s1 = (TH1D *)lim[i][0]->Geth1(name)->Rebin(rebin,"tlimitsig1");
    //    s1->Scale(data[i]->totaltime);
    timet+=data[i]->totaltime;
    bg1->Reset();
    d1->Sumw2();
    bg1->Sumw2();
    s1->Sumw2();

    for(Int_t j=0;j<nbgr[i];j++){
      bg1->Add(bgr[i][j]->Geth1(name),bgr[i][j]->norm);
    }
    for(Int_t j=0;j<nsig[i];j++){
      bg1->Add(sig[i][j]->Geth1(name),sig[i][j]->norm);
    }
    bg1->Rebin(rebin);
    Int_t size=d1->GetNbinsX();
    for(Int_t j=1;j<size+1;j++){
      d->SetBinContent(size*i+j,d1->GetBinContent(j));
      bg->SetBinContent(size*i+j,bg1->GetBinContent(j));
      s->SetBinContent(size*i+j,s1->GetBinContent(j));
      d->SetBinError(size*i+j,d1->GetBinError(j));
      bg->SetBinError(size*i+j,bg1->GetBinError(j));
      s->SetBinError(size*i+j,s1->GetBinError(j));
    }
    delete d1;
    delete bg1;
    delete s1;
  }
  Double_t c,xc;

  TConfidenceLevel *result = MCLimit_basic(d,bg,s,c,xc,es,eb);

  Int_t i;

  anaout <<"Result of TLimit run:"<<endl;
  anaout <<"\t Events in data = "<<d->Integral()<<endl;
  anaout <<"\t Background events = "<<bg->Integral()<<endl;
  anaout <<"\t Signal events <"<<s->Integral()*c<<" @ "<<100*(1-xc)<<"% CL, efficiency "<<s->Integral()/timet<<endl;
  lim[0][0]->setnormalization(0,c);

  TCanvas *mp= new TCanvas("tttch","MCLimit ch",800,800);
  mp->Divide(2,2);
  mp->cd(1);
  d->DrawCopy();
  mp->cd(2);
  bg->DrawCopy();
  mp->cd(3);
  s->DrawCopy();
  mp->cd(4);
  result->Draw();

  delete d;
  delete bg;
  delete s;
  return 0;
}

Int_t MCLimit2dch(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[],bana10 * lim[][MAXCOMP], Int_t nlim[],Int_t nchnls,Int_t rebinx, Int_t rebiny,Double_t es, Double_t eb){

 if(nlim[0] !=1){
    std::cout<<" MCLimit can put limit on one component only. But "<<nlim<<" components were provided"<<endl;
    std::cout<<" No limit was set"<<endl;
    return 1;
  }

  TH2D *sn2ee;
  TH1D *sproj,*bproj,*dproj;
  Double_t sbuf[20000],dbuf[20000],bbuf[20000];
  Int_t nbins;
  nbins=0;
  Int_t i,j;
  Double_t timet=0;
  Bool_t fstop;

  for(Int_t k=0;k<nchnls;k++){
  //Init histograms
    TH2D *b2eeo = (TH2D*)data[k]->Geth2(name)->Clone("mcl");
    TH2D *s2eeo = (TH2D*)lim[k][0]->Geth2(name);
    TH2D *d2eeo = (TH2D*)data[k]->Geth2(name);
    b2eeo->Reset();
    b2eeo->Sumw2();
    for(i=0;i<nbgr[k];i++){
      b2eeo->Add(bgr[k][i]->Geth2(name),bgr[k][i]->norm);
    }
    for(i=0;i<nsig[k];i++){
      b2eeo->Add(sig[k][i]->Geth2(name),sig[k][i]->norm);
    }


    Int_t multx=rebinx;
    Int_t multy=rebiny;

    Double_t binx,biny;
    Double_t bmin=1;

    TH2D *b2ee=(TH2D*)b2eeo->Clone("b2ee");
    TH2D *s2ee=(TH2D*)s2eeo->Clone("s2ee");
    TH2D *d2ee=(TH2D*)d2eeo->Clone("d2ee");

    b2ee->Rebin2D(multx,multy);
    s2ee->Rebin2D(multx,multy);
    d2ee->Rebin2D(multx,multy);
    //    s2ee->Scale(data[k]->totaltime);
    timet+=data[k]->totaltime;

    for(i=1;i<=s2ee->GetNbinsX();i++){
      for(j=1;j<=s2ee->GetNbinsY();j++){
	sbuf[nbins]=s2ee->GetBinContent(i,j);
	dbuf[nbins]=d2ee->GetBinContent(i,j);
	bbuf[nbins]=b2ee->GetBinContent(i,j);
	if(bbuf[nbins]!=0 && bbuf[nbins]<bmin) bmin=bbuf[nbins];
	//      if(dbuf[nbins]==0 && bbuf[nbins]==0) sbuf[nbins]=0;
	nbins++;
      }
    }
  }
  sproj= new TH1D("sproj","",nbins,0.,1.);
  dproj= new TH1D("dproj","",nbins,0.,1.);
  bproj= new TH1D("bproj","",nbins,0.,1.);

  for(i=0;i<nbins;i++){
   fstop=false;
   if(bbuf[i]==0 && dbuf[i]>0 &&!fstop){
     bbuf[i]=dbuf[i];
     cout<<"BIN SIZE TOO SMALL NON-ZERO DATA FOR ZERO BGR PREDICTED!!!! "<<i<<endl;
     // cout << "STOP CALCULATION"<<endl;
     //    fstop=true;
     //continue;
   }
   if(fstop) continue;
   sproj->SetBinContent(i+1,sbuf[i]);
   dproj->SetBinContent(i+1,dbuf[i]);
   bproj->SetBinContent(i+1,bbuf[i]);
 }
  Double_t c,xc;
  MCLimit_basic(dproj,bproj,sproj,c,xc,es,eb);

  anaout <<"Result of 2D TLimit run:"<<endl;
  anaout <<"\t Events in data = "<<dproj->Integral()<<endl;
  anaout <<"\t Background events = "<<bproj->Integral()<<endl;
  anaout <<"\t Signal events <"<<sproj->Integral()*c<<" @ "<<100*(1-xc)<<"% CL, efficiency "<<sproj->Integral()/timet<<endl;
  
  lim[0][0]->setnormalization(0,c);
  delete sproj;
  delete dproj;
  delete bproj;

  return 0;

}
Int_t MCLimit2d(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebinx,Int_t rebiny,Double_t es, Double_t eb){
  //init histogram errors
  Float_t rintegral[400];

  if(nlim !=1){
    std::cout<<" MCLimit can put limit on one component only. But "<<nlim<<" components were provided"<<endl;
    std::cout<<" No limit was set"<<endl;
    return 1;
  }
  //Init histograms
  TH2D *b2eeo = (TH2D*)data->Geth2(name)->Clone("mcl");
  TH2D *s2eeo = (TH2D*)lim[0]->Geth2(name);
  TH2D *d2eeo = (TH2D*)data->Geth2(name);
  b2eeo->Reset();
  b2eeo->Sumw2();
  for(Int_t i=0;i<nbgr;i++){
    b2eeo->Add(bgr[i]->Geth2(name),bgr[i]->norm);
  }
  for(Int_t i=0;i<nsig;i++){
    b2eeo->Add(sig[i]->Geth2(name),sig[i]->norm);
  }

  Int_t i,j;

  Bool_t fstop;
  Int_t multx=rebinx;
  Int_t multy=rebiny;

  TH2D *sn2ee;
  TH1D *sproj,*bproj,*dproj;
  Double_t sbuf[10000],dbuf[10000],bbuf[10000];
  Int_t nbins;
  Double_t binx,biny;
  Double_t bmin=1;

  TH2D *b2ee=(TH2D*)b2eeo->Clone("b2ee");
  TH2D *s2ee=(TH2D*)s2eeo->Clone("s2ee");
  TH2D *d2ee=(TH2D*)d2eeo->Clone("d2ee");

  b2ee->Rebin2D(multx,multy);
  s2ee->Rebin2D(multx,multy);
  d2ee->Rebin2D(multx,multy);

  nbins=0;
  for(i=1;i<=s2ee->GetNbinsX();i++){
    for(j=1;j<=s2ee->GetNbinsY();j++){
      sbuf[nbins]=s2ee->GetBinContent(i,j);
      dbuf[nbins]=d2ee->GetBinContent(i,j);
      bbuf[nbins]=b2ee->GetBinContent(i,j);
      if(bbuf[nbins]!=0 && bbuf[nbins]<bmin) bmin=bbuf[nbins];
      //      if(dbuf[nbins]==0 && bbuf[nbins]==0) sbuf[nbins]=0;
      nbins++;
    }
  }

 sproj= new TH1D("sproj","",nbins,0.,1.);
 dproj= new TH1D("dproj","",nbins,0.,1.);
 bproj= new TH1D("bproj","",nbins,0.,1.);

 for(i=0;i<nbins;i++){
   fstop=false;
   if(bbuf[i]==0 && dbuf[i]>0 &&!fstop){
     bbuf[i]=dbuf[i];
     cout<<"BIN SIZE TOO SMALL NON-ZERO DATA FOR ZERO BGR PREDICTED!!!! "<<i<<endl;
     // cout << "STOP CALCULATION"<<endl;
     //    fstop=true;
     //continue;
  }
  if(fstop) continue;
  sproj->SetBinContent(i+1,sbuf[i]);
  dproj->SetBinContent(i+1,dbuf[i]);
  bproj->SetBinContent(i+1,bbuf[i]);
 }


  Double_t c,xc;
  MCLimit_basic(dproj,bproj,sproj,c,xc,es,eb);

  anaout <<"Result of 2D TLimit run:"<<endl;
  anaout <<"\t Events in data = "<<d2ee->Integral()<<endl;
  anaout <<"\t Background events = "<<b2ee->Integral()<<endl;
  anaout <<"\t Signal events <"<<s2ee->Integral()*c<<" @ "<<100*(1-xc)<<"% CL, efficiency "<<s2ee->Integral()/lim[0]->totaltime<<endl;
  
  lim[0]->setnormalization(0,c);
  delete sproj;
  delete dproj;
  delete bproj;

  return 0;
}
