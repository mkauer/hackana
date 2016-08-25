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

#include "dolfit.hpp"

// Fit histograms. Can do simple fit (signal = data-bgr) and multiple distributions chi2 fit with ROOT FractionFitter package 

// Fit functions
Int_t FitSignal(const char *name, bana10 *data, bana10 **bgr, Int_t nbgr, bana10 **sig, Int_t nsig, Int_t lx, Int_t ux){
  
  std::cout<<"\n\nUSING REGULAR SINGLE-CHANNEL FIT!!\n\n";
  
  //init histogram errors
  Float_t rintegral[400];
  
  data->Geth1(name)->Sumw2();
  TH1D * sum =(TH1D *) data->Geth1(name)->Clone("sum");  
  Double_t errors[1000];
  Int_t i1;
  sum->Reset();
  sum->Sumw2();
  //sum2->Sumw2();
  //subtract known background
  TH1D *ref=data->Geth1(name);
  sum->Add(ref,1);
  //  std::cout<<endl<<"data ";
  for(i1=0;i1<ref->GetNbinsX();i1++){
    errors[i1]=pow(ref->GetBinError(i1),2);
    // std::cout<<errors[i1]<<" ";
  }
   for(Int_t i=0;i<nbgr;i++){
     //std::cout<<endl<<bgr[i]->Getname();
     ref=bgr[i]->Geth1(name);
     sum->Add(ref,-bgr[i]->norm);
     for(i1=0;i1<ref->GetNbinsX();i1++){
       errors[i1]+=pow(ref->GetBinError(i1),2)*pow(bgr[i]->norm,2)+pow(ref->GetBinContent(i1)*bgr[i]->enorm,2);
       //std::cout<<errors[i1]<<" ";
     } 
   }
   
   // delete negative bins
   if( nsig > 1) {
     for(Int_t i=1;i<=sum->GetNbinsX();i++){
       //std::cout<<errors[i]<<" ";
       if(sum->GetBinContent(i)<0) sum->SetBinContent(i,0.);
     }
   }
   // sqrt from errors
   for(Int_t i=1;i<=sum->GetNbinsX();i++){
     errors[i]=sqrt(errors[i]);
   }
   sum->SetError(errors);
   // cout<<" Data - BGR = "<< sum->Integral(lx,ux);
  //Init MC signal massive
  TObjArray * mc = new TObjArray(nsig);
  for(Int_t i=0;i<nsig;i++){
    if(!sig[i]->dependent){
      TH1D * sum2 =(TH1D *) sig[i]->Geth1(name)->Clone("sum2");
      sum2->Reset();
      sum2->Sumw2();
      sum2->Add(sig[i]->Geth1(name),1);
      for(Int_t j=0;j<sig[i]->nrelative;j++){
	bana10 * child=sig[i]->relative[j];
	sum2->Add(child->Geth1(name),1.*sig[i]->relfactor[j]);
      }
      mc->Add(sum2);
      rintegral[i]=sum2->Integral(lx,ux);
    }
  }
  if(mc->GetEntries()==1){
    //just adjust MC if there is only one MC sample to fit
    Int_t result=0;
    for(Int_t i=0;i<nsig;i++){
      if(!sig[i]->dependent){
	Double_t v,ev;
	v=1.;//100% proportion belongs to this MC
	result++;
	Double_t factor=sum->Integral(lx,ux)/rintegral[i];
	v=v*factor;
	Double_t sumerror=0;
	for(i1=lx;i1<=ux;i1++) sumerror+=pow(sum->GetBinError(i1),2);
	//std::cout<<" Signal "<<rintegral[i]<<endl;
	ev=v*sqrt(sumerror/pow(sum->Integral(lx,ux),2) + 1./rintegral[i]);
	sig[i]->setnormalization(v,ev);
	//cout<<sig[i]->Getname()<<" set v="<<v<<" ev="<<ev<< endl;
	for(Int_t j=0;j<sig[i]->nrelative;j++){
	  Float_t ftot=sig[i]->relfactor[j];
	  sig[i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	  cout<<sig[i]->relative[j]->Getname()<<" set v="<<v*ftot<<endl;
	}
      }
    }
  }else{ 
    TFractionFitter* fit = new TFractionFitter(sum, mc); // initialise 
    fit->SetRangeX(lx,ux);
    Int_t status = fit->Fit();               // perform the fit
    cout<<"Fit done with status "<<status<<" Chi2="<<fit->GetChisquare()<<"/"<<fit->GetNDF()<<endl;
    Int_t result=0;
    for(Int_t i=0;i<nsig;i++){
      if(!sig[i]->dependent){
	Double_t v,ev;
	fit->GetResult(result,v,ev);
	result++;
	Double_t factor=sum->Integral(lx,ux)/rintegral[i];
	v=v*factor;
	ev=ev*factor;
	sig[i]->setnormalization(v,ev);
	cout<<sig[i]->Getname()<<" set v="<<v<<endl;
	for(Int_t j=0;j<sig[i]->nrelative;j++){
	  Float_t ftot=sig[i]->relfactor[j];
	  sig[i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	  cout<<sig[i]->relative[j]->Getname()<<" set v="<<v*ftot;
	}
      }
    }
  }
}


Int_t FitSignal(const char *name, bana10 *data, bana10 **bgr, Int_t nbgr, bana10 **sig, Int_t nsig, Float_t blow, Float_t bup, Int_t rebin){
  // Improved vertion of fit function
  // Can do some preparation work on the histograms before fitting
  // blow, bup -- boundaries of the region to fit
  // rebin -- if one wants to rebin the content of histogram before the fit (helps to improve TFractionFitter convergence)  
 
  std::cout<<"\n\nUSING SINGLE-CHANNEL IMPROVED FIT!!\n\n";
  
  //init histogram errors
  Float_t rintegral[400];
  Int_t ux,lx;

  data->Geth1(name)->Sumw2();
  TH1D * sum =(TH1D *) data->Geth1(name)->Clone("sum");  
  Double_t errors[10000];
  Int_t i1;
  sum->Reset();
  sum->Sumw2();
  //sum2->Sumw2();
  //subtract known background
  TH1D *ref=data->Geth1(name);
  sum->Add(ref,1);
  //  std::cout<<endl<<"data ";
  for(i1=0;i1<ref->GetNbinsX();i1++){
    errors[i1]=pow(ref->GetBinError(i1),2);
    // std::cout<<errors[i1]<<" ";
  }
   for(Int_t i=0;i<nbgr;i++){
     //std::cout<<endl<<bgr[i]->Getname();
     ref=bgr[i]->Geth1(name);
     sum->Add(ref,-bgr[i]->norm);
     for(i1=0;i1<ref->GetNbinsX();i1++){
       errors[i1]+=pow(ref->GetBinError(i1),2)*pow(bgr[i]->norm,2)+pow(ref->GetBinContent(i1)*bgr[i]->enorm,2);
       //std::cout<<errors[i1]<<" ";
     } 
   }
   
   // set errors
   for(Int_t i=1;i<=sum->GetNbinsX();i++){
     errors[i]=sqrt(errors[i]);
     //std::cout<<errors[i]<<" ";
   }
   sum->SetError(errors);

   sum->Rebin(rebin);
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


   // delete negative bins
   for(Int_t i=1;i<=sum->GetNbinsX();i++){
     if(sum->GetBinContent(i)<0) sum->SetBinContent(i,0.);
   }
   cout<<" Data - BGR = "<< sum->Integral(lx,ux);
  //Init MC signal massive

  TObjArray * mc = new TObjArray(nsig);
  for(Int_t i=0;i<nsig;i++){
    if(!sig[i]->dependent){
      TH1D * sum2 =(TH1D *) sig[i]->Geth1(name)->Clone("sum2");
      sum2->Reset();
      sum2->Sumw2();
      sum2->Add(sig[i]->Geth1(name),1);
      for(Int_t j=0;j<sig[i]->nrelative;j++){
	bana10 * child=sig[i]->relative[j];
	sum2->Add(child->Geth1(name),1.*sig[i]->relfactor[j]);
      }
      sum2->Rebin(rebin);
      mc->Add(sum2);
      rintegral[i]=sum2->Integral(lx,ux);
    }
  }
  if(mc->GetEntries()==1){
      //just adjust MC if there is only one MC sample to fit
      Int_t result=0;
      for(Int_t i=0;i<nsig;i++){
	if(!sig[i]->dependent){
	  Double_t v,ev;
	  v=1.;//100% proportion belongs to this MC
	  result++;
	  Double_t factor=sum->Integral(lx,ux)/rintegral[i];
	  v=v*factor;
	  Double_t sumerror=0;
	  for(i1=lx;i1<=ux;i1++) sumerror+=pow(sum->GetBinError(i1),2);
	  std::cout<<" Signal "<<rintegral[i]<<endl;
	  ev=v*sqrt(sumerror/pow(sum->Integral(lx,ux),2) + 1./rintegral[i]);
	  sig[i]->setnormalization(v,ev);
	  cout<<sig[i]->Getname()<<" set v="<<v<<" ev="<<ev<< endl;
	  for(Int_t j=0;j<sig[i]->nrelative;j++){
	    Float_t ftot=sig[i]->relfactor[j];
	    sig[i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	    cout<<sig[i]->relative[j]->Getname()<<" set v="<<v*ftot<<endl;
	  }
	}
      }
  }else{
    TFractionFitter* fit = new TFractionFitter(sum, mc); // initialise 
    fit->SetRangeX(lx,ux);
    Int_t status = fit->Fit();               // perform the fit
    cout<<"Fit done with status "<<status<<" Chi2 = "<<fit->GetChisquare()<<" / "<<fit->GetNDF()<<endl;
    Int_t result=0;
    for(Int_t i=0;i<nsig;i++){
      if(!sig[i]->dependent){
	Double_t v,ev;
	fit->GetResult(result,v,ev);
	result++;
	Double_t factor=sum->Integral(lx,ux)/rintegral[i];
	v=v*factor;
	ev=ev*factor;
	sig[i]->setnormalization(v,ev);
	cout<<sig[i]->Getname()<<" set v="<<v<<endl;
	for(Int_t j=0;j<sig[i]->nrelative;j++){
	  Float_t ftot=sig[i]->relfactor[j];
	  sig[i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	  cout<<sig[i]->relative[j]->Getname()<<" set v="<<v*ftot;
	}
      }
    }
  }
}


// Fit hystograms in several channels simultaneously
Int_t FitSignalCh(const char *name, bana10 *data[], bana10 *bgr[][MAXCOMP], Int_t nbgr[], bana10 *sig[][MAXCOMP], Int_t nsig[], Int_t nchnls, Int_t lx, Int_t ux){
  
  std::cout<<"\n\nUSING MULTI-CHANNEL FIT!!\n\n";
  
  //init histogram errors
  Float_t rintegral[400];
  
  //template histogram;
  Int_t nbin0=data[0]->Geth1(name)->GetNbinsX();  
  TH1D * tdata= new TH1D("tdata","tdata",nbin0*nchnls,0.,1.);
  
  for(Int_t k=0;k<nchnls;k++){
    for(Int_t ib=1;ib<=nbin0;ib++){
      Float_t bcont=0;
      Float_t berror=0;
      //sum2->Sumw2();
      //subtract known background
      bcont+=data[k]->Geth1(name)->GetBinContent(ib);
      berror+=data[k]->Geth1(name)->GetBinError(ib);
      for(Int_t i=0;i<nbgr[k];i++){
	bcont+=-bgr[k][i]->norm*bgr[k][i]->Geth1(name)->GetBinContent(ib);
	Float_t herror=-bgr[k][i]->norm*bgr[k][i]->Geth1(name)->GetBinError(ib);
	berror=sqrt(berror*berror+herror*herror);
      }
      if(bcont<0) bcont=0;
      tdata->SetBinContent(nbin0*k+ib,bcont);
      tdata->SetBinError(nbin0*k+ib,berror);
    }
  }
  //Init MC signal massive
  TObjArray * mc = new TObjArray(nsig[0]);
  for(Int_t i=0;i<nsig[0];i++){
    if(!sig[0][i]->dependent){
      TH1D * sum2 = new TH1D((string("t")+sig[0][i]->Getname()).c_str(),"tdata",nbin0*nchnls,0.,1.);
      for(Int_t k=0;k<nchnls;k++){
	for(Int_t ib=1;ib<=nbin0;ib++){
	  Float_t bcont=0,berror=0;
	  bcont=sig[k][i]->Geth1(name)->GetBinContent(ib);
	  berror=sig[k][i]->Geth1(name)->GetBinError(ib);
	  for(Int_t j=0;j<sig[k][i]->nrelative;j++){
	    bana10 * child=sig[k][i]->relative[j];
	    // sum2->Add(child->Geth1(name),sig[k][i]->mcevga*sig[k][i]->relfactor[j]/child->mcevga);
	    bcont+=sig[k][i]->mcevga*sig[k][i]->relfactor[j]/child->mcevga*child->Geth1(name)->GetBinContent(ib);
	    Float_t herror=sig[k][i]->mcevga*sig[k][i]->relfactor[j]/child->mcevga*child->Geth1(name)->GetBinError(ib);
	    berror=sqrt(berror*berror+herror*herror);
	  }
	  sum2->SetBinContent(nbin0*k+ib,bcont);
	  sum2->SetBinError(nbin0*k+ib,berror);
	}
      }
      mc->Add(sum2);
      rintegral[i]=sum2->Integral();
    }
  }
  
  // std::cout<<"FitSignalCh:: call to TFractionFitter"<<std::endl;
  TFractionFitter* fit = new TFractionFitter(tdata, mc); // initialise 
  fit->Constrain(0,0.,1.);
  //  fit->SetRangeX(lx,ux);
  Int_t status = fit->Fit();               // perform the fit
  cout<<"Fit done with status "<<status<<" Chi2="<<fit->GetChisquare()<<"/"<<fit->GetNDF()<<endl;
  Int_t result=0;
  for(Int_t i=0;i<nsig[0];i++){
    Double_t v,ev;
    fit->GetResult(result,v,ev);
    result++;
    Double_t factor=tdata->Integral()/rintegral[i];
    v=v*factor;
    ev=ev*factor;
    for(Int_t k=0;k<nchnls;k++){
      if(!sig[k][i]->dependent){
	sig[k][i]->setnormalization(v,ev);
	cout<<sig[k][i]->Getname()<<" set v="<<v<<endl;
	for(Int_t j=0;j<sig[k][i]->nrelative;j++){
	  Float_t ftot=sig[k][i]->relfactor[j]/sig[k][i]->relative[j]->mcevga*sig[k][i]->mcevga;
	  sig[k][i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	  cout<<sig[k][i]->relative[j]->Getname()<<" set v="<<v*ftot;
	}
      }
    }
  }
}


Int_t LFitSignal(const char *name, bana10 *data, bana10 **bgr, Int_t nbgr, bana10 **sig, Int_t nsig, Float_t blow, Float_t bup, Int_t rebin, Double_t scalebkgs){
  // Improved vertion of fit function, use likelihood fit
  // Can do some preparation work on the histograms before fitting
  // blow, bup -- boundaries of the region to fit
  // rebin -- if one wants to rebin the content of histogram before the fit (helps to improve TFractionFitter convergence)  
  // scalebkgs -- you can scale the backgrounds to check systematics of fit in a quick and dirty way!
  
  std::cout<<"\n\nUSING SINGLE-CHANNEL LIKELIHOOD FIT!!\n\n";
  
  //do rebining in the fitter, after the energy cut
  Int_t fitter_rebin = rebin;
  rebin = 1;

  //init histogram errors
  Float_t rintegral[400];
  Int_t ux,lx;

  TH1 *dh,*bh,*sh;// data,bgr and signal histograms for likelihhod test

  data->Geth1(name)->Sumw2();
  dh =(TH1D *) data->Geth1(name)->Clone("sum");  
  //
  // do some gain shifting thing here eventually!!!
  //
  TH1D * sum =(TH1D *) data->Geth1(name)->Clone("sum");  
  Double_t errors[10000];
  Int_t i1;
  sum->Reset();
  sum->Sumw2();
  //sum2->Sumw2();
  //subtract known background
  TH1D *ref=data->Geth1(name);
  //  dh->Add(ref,1);
  // reset error bars for background
  for(i1=0;i1<ref->GetNbinsX();i1++){
    errors[i1]=0;
  }
  for(Int_t i=0;i<nbgr;i++){
    ref=bgr[i]->Geth1(name);
    sum->Add(ref,bgr[i]->norm);
    //bgr[i]->norm;
    //sum->Add(ref,bgr[i]->Scale(1.1));
    for(i1=0;i1<ref->GetNbinsX();i1++){
      errors[i1]+=pow(ref->GetBinError(i1),2)*pow(bgr[i]->norm,2)+pow(ref->GetBinContent(i1)*bgr[i]->enorm,2);
    } 
  }
  
  if(scalebkgs!=1.0){
    cout<<"\n\n\t !CAUTION! ==> The backgrounds are being scaled to = "<<scalebkgs<<"\n\n\n";
    sum->Scale(scalebkgs);
  }
  
  // set errors
  for(Int_t i=1;i<=sum->GetNbinsX();i++){
    errors[i]=sqrt(errors[i]);
    //std::cout<<errors[i]<<" ";
  }
  sum->SetError(errors);
  sum->Rebin(rebin);
  bh=sum;
  
  dh->Rebin(rebin);
  
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
   TObjArray * mc = new TObjArray(nsig);
   for(Int_t i=0;i<nsig;i++){
     if(!sig[i]->dependent){
       sum2 =(TH1D *) sig[i]->Geth1(name)->Clone("sum2");
       sum2->Reset();
       sum2->Sumw2();
       sum2->Add(sig[i]->Geth1(name),1);
       for(Int_t j=0;j<sig[i]->nrelative;j++){
	 bana10 * child=sig[i]->relative[j];
	 sum2->Add(child->Geth1(name),1.*sig[i]->relfactor[j]);
       }
       sum2->Rebin(rebin);
       mc->Add(sum2);
       rintegral[i]=sum2->Integral(lx,ux);
     }
   }
   if(mc->GetEntries()==1){
     sh=sum2;
     Float_t signal;
     Float_t sigma;
     
     dolfit fitter;
     fitter.Init(dh,bh,sh,lx,ux,fitter_rebin);
     fitter.fit();
     
     //just adjust MC if there is only one MC sample to fit
     Int_t result=0;
     for(Int_t i=0;i<nsig;i++){
       if(!sig[i]->dependent){
	 Double_t v,ev;
	 v=fitter.get_value();//100% proportion belongs to this MC
	 result++;
	 ev=fitter.get_righterror();
	 sig[i]->setnormalization(v,ev);
	 cout<<sig[i]->Getname()<<" set v="<<v<<" ev="<<ev<< endl;
	 cout<<"Left error="<<fitter.get_lefterror()<<" Right error="<<ev<< endl;
	 for(Int_t j=0;j<sig[i]->nrelative;j++){
	   Float_t ftot=sig[i]->relfactor[j];
	   sig[i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	   cout<<sig[i]->relative[j]->Getname()<<" set v="<<v*ftot<<endl;
	 }
       }
     }
   }
}

Int_t LFitSignalCh(const char * name, bana10 *data[],bana10 *bgr[][MAXCOMP],Int_t nbgr[],bana10 * sig[][MAXCOMP],Int_t nsig[], Int_t nchnls, Float_t blow, Float_t bup, Int_t rebin){
  // Improved vertion of fit function, use likelihood fit
  // Can do some preparation work on the histograms before fitting
  // blow, bup -- boundaries of the region to fit
  // rebin -- if one wants to rebin the content of histogram before the fit (helps to improve TFractionFitter convergence)  
  
  std::cout<<"\n\nUSING MULTI-CHANNEL IMPROVED LIKELIHOOD FIT!!\n\n";


  //do rebinning inside the fitter
  Int_t fitter_rebin = rebin;
  rebin=1;
  //init histogram errors
  Float_t rintegral[400];
  Int_t ux,lx;
  //storage for stakced histigramms for data, signal anf background
  Float_t darr[4000];
  Float_t edarr[4000];
  Float_t barr[4000];
  Float_t ebarr[4000];
  Float_t sarr[4000];
  Float_t esarr[4000];
  //number of bins in global stacked histogram 
  Int_t nbins_global=0;
  
  for(Int_t nch=0;nch<nchnls;nch++){
    cout<<"LFitCH: start to scan channel "<<nch<<endl;
    
    TH1 *dh,*bh,*sh;// data,bgr and signal histograms for likelihhod test
    
    data[nch]->Geth1(name)->Sumw2();
    dh =(TH1D *) data[nch]->Geth1(name)->Clone("sum");  
    TH1D * sum =(TH1D *) data[nch]->Geth1(name)->Clone("sum");  
    Double_t errors[10000];
    Int_t i1;
    sum->Reset();
    sum->Sumw2();
    //sum2->Sumw2();
    //subtract known background
    TH1D *ref=data[nch]->Geth1(name);
    //  dh->Add(ref,1);
    // reset error bars for background
    for(i1=0;i1<ref->GetNbinsX();i1++){
      errors[i1]=0;
    }
    for(Int_t i=0;i<nbgr[nch];i++){
      ref=bgr[nch][i]->Geth1(name);
      sum->Add(ref,bgr[nch][i]->norm);
      for(i1=0;i1<ref->GetNbinsX();i1++){
	errors[i1]+=pow(ref->GetBinError(i1),2)*pow(bgr[nch][i]->norm,2)+pow(ref->GetBinContent(i1)*bgr[nch][i]->enorm,2);
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
    
    dh->Rebin(rebin);
    cout<<"LFitCH: DATA is done"<<endl;
    
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

    cout<<"LFitCH: Boundaries done"<<endl;

    TH1D *sum2;
    TObjArray * mc = new TObjArray(nsig[nch]);
    for(Int_t i=0;i<nsig[nch];i++){
      if(!sig[nch][i]->dependent){
	sum2 =(TH1D *) sig[nch][i]->Geth1(name)->Clone("sum2");
	sum2->Reset();
	sum2->Sumw2();
	sum2->Add(sig[nch][i]->Geth1(name),1);
	for(Int_t j=0;j<sig[nch][i]->nrelative;j++){
	  bana10 * child=sig[nch][i]->relative[j];
	  sum2->Add(child->Geth1(name),1.*sig[nch][i]->relfactor[j]);
	}
	sum2->Rebin(rebin);
	mc->Add(sum2);
	rintegral[i]=sum2->Integral(lx,ux);
      }
    }
    sh=sum2;
    cout<<"LFitCH: Signals done"<<endl;

    //copy historgams from channel nch to global array
    for(Int_t counter=lx;counter<=ux;counter++){
      darr[nbins_global] = dh->GetBinContent(counter);
      edarr[nbins_global] = dh->GetBinError(counter);
      barr[nbins_global] = bh->GetBinContent(counter);
      ebarr[nbins_global] = bh->GetBinError(counter);
      sarr[nbins_global] = sh->GetBinContent(counter);
      esarr[nbins_global] = sh->GetBinError(counter);
      nbins_global++;
    }
    delete dh;
    delete bh;
    delete sh;
    delete mc;
    
  }
  cout<<"LFitCH: All channels scanned"<<endl;

  //create dummy histos
  TH1* dh = new TH1F("dh","dh",nbins_global,0.,1.);
  TH1* bh = new TH1F("bh","bh",nbins_global,0.,1.);
  TH1* sh = new TH1F("sh","sh",nbins_global,0.,1.);
  
  //copy global content
  for(Int_t count=0;count<nbins_global;count++){
    dh->SetBinContent(count+1,darr[count]);
    dh->SetBinError(count+1,edarr[count]);
    bh->SetBinContent(count+1,barr[count]);
    bh->SetBinError(count+1,ebarr[count]);
    sh->SetBinContent(count+1,sarr[count]);
    sh->SetBinError(count+1,esarr[count]);
  }
  
  TCanvas *ctest = new TCanvas("ctest","ctest",600,600);
  ctest->Divide(2,2);
  ctest->cd(1);
  dh->Draw();
  ctest->cd(2);
  bh->Draw();
  ctest->cd(3);
  sh->Draw();
  cout<<"LFitCH: Histograms ready"<<endl;

    Float_t signal;
    Float_t sigma;
    
    dolfit fitter;
    fitter.Init(dh,bh,sh,1,nbins_global,fitter_rebin);
    fitter.fit();
    
    //just adjust MC if there is only one MC sample to fit
    for(Int_t nch=0;nch<nchnls;nch++){
      Int_t result=0;
      for(Int_t i=0;i<nsig[nch];i++){
	if(!sig[nch][i]->dependent){
	  Double_t v,ev;
	  v=fitter.get_value();//100% proportion belongs to this MC
	  result++;
	  ev=fitter.get_righterror();
	  sig[nch][i]->setnormalization(v,ev);
	  cout<<sig[nch][i]->Getname()<<" set v="<<v<<" ev="<<ev<< endl;
	  cout<<"Left error="<<fitter.get_lefterror()<<" Right error="<<ev<< endl;
	  for(Int_t j=0;j<sig[nch][i]->nrelative;j++){
	    Float_t ftot=sig[nch][i]->relfactor[j];
	    sig[nch][i]->relative[j]->setnormalization(v*ftot,ev*ftot);
	    cout<<sig[nch][i]->relative[j]->Getname()<<" set v="<<v*ftot<<endl;
	  }
	}
      }
    }

    //    delete dh;
    // delete bh;
    //delete sh;
}
