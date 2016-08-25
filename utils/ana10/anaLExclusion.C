#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include "../../utils/h10.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TGApplication.h"
#include "TLine.h"
#include "TPaveLabel.h"


#include "dollimit.hpp"

Int_t LExclusion(const char * name, bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig, bana10 ** lim, Int_t nlim, Float_t blow, Float_t bup, Int_t rebin){
  //
  //   Estimates sensetivity for a given background level and signal efficiency
  //   Use ratio of binned likelihood for the histogram
  //

  // Can do some preparation work on the histograms
  // blow, bup -- boundaries of the region to fit
  // rebin -- if one wants to rebin the content of histogram before the fit (helps to improve TFractionFitter convergence)  


  //init histogram errors
  Float_t rintegral[400];
  Int_t ux,lx;

  TH1 *bh,*sh;// bgr and signal histograms for likelihhod test

  TH1D * sum =(TH1D *) bgr[0]->Geth1(name)->Clone("lesum");  
  Double_t errors[10000];
  Int_t i1;
  sum->Reset();
  //  sum->Sumw2();

  //subtract known background
  TH1D *ref;

  // reset error bars for background
  for(i1=0;i1<sum->GetNbinsX();i1++){
    errors[i1]=0;
  }

  for(Int_t i=0;i<nbgr;i++){
    //std::cout<<"LExclusion, process bgr "<<bgr[i]->Getname()<<std::endl;
    ref=bgr[i]->Geth1(name);
    sum->Add(ref,bgr[i]->norm);
    for(i1=0;i1<sum->GetNbinsX();i1++){
      errors[i1]+=pow(ref->GetBinError(i1),2)*pow(bgr[i]->norm,2)+pow(ref->GetBinContent(i1)*bgr[i]->enorm,2);
    } 
  }
  for(Int_t i=0;i<nsig;i++){
    ref=sig[i]->Geth1(name);
    sum->Add(ref,sig[i]->norm);
    for(i1=0;i1<sum->GetNbinsX();i1++){
      errors[i1]+=pow(ref->GetBinError(i1),2)*pow(sig[i]->norm,2)+pow(ref->GetBinContent(i1)*sig[i]->enorm,2);
    } 
  }
   
   // set errors
   for(Int_t i=1;i<=sum->GetNbinsX();i++){
     errors[i]=sqrt(errors[i]);
     //std::cout<<errors[i]<<" ";
   }
   sum->SetError(errors);
   sum->Rebin(rebin);
   bh=sum;


   if(blow!=0){
      Double_t ledges[10000];
      sum->GetLowEdge(ledges);
      lx=0;
      ux=0;
      for(Int_t k=0;k<sum->GetNbinsX()+1;k++){
	if(blow<ledges[k] && lx==0) lx = k-1;
	if(bup<ledges[k] && ux==0) ux = k-1;
      }
      if(ux==0) ux=sum->GetNbinsX()+1;
      
    }else{
      lx=0;
      ux=sum->GetNbinsX()+1;
    }

   //std::cout<<"LExclusion. Boundaries "<<lx<<":"<<ux<<"  tot bgr "<<bh->Integral(lx,ux)<<std::endl;

   TH1D *sum2;
   sum2 =(TH1D *) lim[0]->Geth1(name)->Clone("lesum2");
   sum2->Reset();
   //sum2->Sumw2();
   sum2->Add(lim[0]->Geth1(name),1);
   sum2->Rebin(rebin);

   sh = sum2;
   //std::cout<<"LExclusion. tot sig "<<sh->Integral(lx,ux)<<std::endl;

   dollimit delimiter;
   delimiter.Init(bh,sh,lx,ux);
   delimiter.fit();
   Float_t ev=delimiter.get_value();//100% proportion belongs to this MC
   lim[0]->setnormalization(0,ev);
   delete sum;
   delete sum2;

}

Int_t LExclusion_pessimistic(const char * name, bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig, bana10 ** lim, Int_t nlim, Float_t blow, Float_t bup, Int_t rebin, Int_t SN){
  //Similar as above, but using pessimistic approach, add 2.3 * n expeted/ N MC
  // for every 0 bin in the background histogram;


  //
  //   Estimates sensetivity for a given background level and signal efficiency
  //   Use ratio of binned likelihood for the histogram
  //

  // Can do some preparation work on the histograms
  // blow, bup -- boundaries of the region to fit
  // rebin -- if one wants to rebin the content of histogram before the fit (helps to improve TFractionFitter convergence)  





  //init histogram errors
  Float_t rintegral[400];
  Int_t ux,lx;

  TH1 *bh,*sh;// bgr and signal histograms for likelihhod test

  TH1D * sum =(TH1D *) bgr[0]->Geth1(name)->Clone("sum");  
  Double_t errors[10000];
  Int_t i1;
  sum->Reset();
  //  sum->Sumw2();

  //subtract known background
  TH1D *ref;

  // reset error bars for background
  for(i1=0;i1<sum->GetNbinsX();i1++){
    errors[i1]=0;
  }

  sum->Rebin(rebin);

  for(Int_t i=0;i<nbgr;i++){
    ref=new TH1D(*bgr[i]->Geth1(name));
    ref->Rebin(rebin);
    for(int j=0;j<ref->GetNbinsX();j++){
      Double_t content = ref->GetBinContent(j);
      if(content==0){

	if(SN ==1){
	  ref->SetBinContent(j,2.3/bgr[i]->mcevga);
	}else{
	  ref->SetBinContent(j,2.3*bgr[i]->totaltime/bgr[i]->mcevga/1000);
	}
      }
    }
    sum->Add(ref,bgr[i]->norm);
    for(i1=0;i1<sum->GetNbinsX();i1++){
      errors[i1]+=pow(ref->GetBinError(i1),2)*pow(bgr[i]->norm,2)+pow(ref->GetBinContent(i1)*bgr[i]->enorm,2);
    } 
    delete ref; 
  }
  for(Int_t i=0;i<nsig;i++){
    ref=new TH1D(*sig[i]->Geth1(name));
    for(int j=0;j<ref->GetNbinsX();j++){
      Double_t content = ref->GetBinContent(j);
      if(content==0){
	if(SN == 1){
	  ref->SetBinContent(j,2.3/bgr[i]->mcevga);
	}else{
	  ref->SetBinContent(j,2.3*bgr[i]->totaltime/bgr[i]->mcevga);
	}
      }
    }
    sum->Add(ref,sig[i]->norm);
    for(i1=0;i1<sum->GetNbinsX();i1++){
      errors[i1]+=pow(ref->GetBinError(i1),2)*pow(sig[i]->norm,2)+pow(ref->GetBinContent(i1)*sig[i]->enorm,2);
    } 
    delete ref;
  }
   // set errors
   for(Int_t i=1;i<=sum->GetNbinsX();i++){
     errors[i]=sqrt(errors[i]);
     //std::cout<<errors[i]<<" ";
   }
   sum->SetError(errors);

   bh=sum;


   if(blow!=0){
      Double_t ledges[10000];
      sum->GetLowEdge(ledges);
      lx=0;
      ux=0;
      for(Int_t k=0;k<sum->GetNbinsX()+1;k++){
	if(blow<ledges[k] && lx==0) lx = k-1;
	if(bup<ledges[k] && ux==0) ux = k-1;
      }
      if(ux==0) ux=sum->GetNbinsX()+1;
      
    }else{
      lx=0;
      ux=sum->GetNbinsX()+1;
    }

   TH1D *sum2;
   sum2 =(TH1D *) lim[0]->Geth1(name)->Clone("sum2");
   sum2->Reset();
   sum2->Sumw2();
   sum2->Add(lim[0]->Geth1(name),1);
   sum2->Rebin(rebin);

   sh = sum2;


   TCanvas *d;
   d = new TCanvas("LEX","pessimistic PDF",800,800);
   d->Divide(2,2);
   d->cd(1);
   bh->DrawCopy();
   d->cd(2);
   sh->DrawCopy();
   d->Update();


   dollimit delimiter;
   delimiter.Init(bh,sh,lx,ux);
   delimiter.fit();
   Float_t ev=delimiter.get_value();//100% proportion belongs to this MC
   lim[0]->setnormalization(0,ev);
   delete sum;
   delete sum2;
}

