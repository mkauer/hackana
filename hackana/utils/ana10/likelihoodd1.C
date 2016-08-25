#include "ana.hpp"

#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "TGApplication.h"
#include "TLine.h"
#include "TPaveLabel.h"
#include "TF1.h"
#include "TRandom.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>

TH1D *lexp,*lchi,*l2b,*lbgr;
Double_t nchi;
Double_t L0;

Double_t lm(Double_t * n2b, Double_t *par){
  // Standard distributions provided:
  //        lexp,lchi,l2b,lbgr

  Int_t nbins=lexp->GetNbinsX();
  Double_t L=0;
  Double_t mean;
  for(Int_t i=1;i<=nbins;i++){
    mean=n2b[0]*l2b->GetBinContent(i)+nchi*lchi->GetBinContent(i)+lbgr->GetBinContent(i);
    Int_t exp=(Int_t)lexp->GetBinContent(i);
    if(mean!=0){
      L-=mean;
      if(exp <10){
	L-=TMath::Log(TMath::Factorial(exp));
      }else{
	L-=(exp*(TMath::Log(exp)-1)+0.5*TMath::Log(2*3.1415926*exp));//Sterling formula
      }
      L+=exp*TMath::Log(mean);
    }
  }
  L+=L0;
  //  cout<<"\t lm output 2b,chi="<<n2b[0]<<" "<<nchi<<" L="<<L<<" mean "<<mean<<endl;
  L=TMath::Exp(L);
  return L;
}

Int_t MajoronLimit(const char * name, bana10 *data,bana10 **bgr,Int_t nbgr,bana10 ** sig,Int_t nsig,bana10 ** lim, Int_t nlim,Int_t rebin){

    lexp=data->Geth1(name);
    lchi=(TH1D *)lim[0]->Geth1(name)->Clone("lchi");
    lchi->Reset();
    lchi->Add(lim[0]->Geth1(name),1./lim[0]->mcevga);
    l2b=(TH1D *)sig[0]->Geth1(name)->Clone("l2b");
    l2b->Reset();
    l2b->Add(sig[0]->Geth1(name),1./sig[0]->mcevga);
    lbgr=(TH1D *)data->Geth1(name)->Clone("lbgr");
    lbgr->Reset();
    for(Int_t i=0;i<nbgr;i++){
      lbgr->Add(bgr[i]->Geth1(name),bgr[i]->norm);
    }

    lexp->Rebin(rebin);
    lchi->Rebin(rebin);
    l2b->Rebin(rebin);
    lbgr->Rebin(rebin);
    Double_t aprox = sig[0]->mcevga*sig[0]->norm;
    Double_t lb=aprox*0.9;
    Double_t ub= aprox*1.1;
    L0 = -TMath::Log(lm(&aprox,&aprox));
    TF1 *lf= new TF1("likelyhood",lm,0,ub,0);
    const Int_t nchimax=10000;
    Double_t Ichi[nchimax];
    Double_t step=1;
    Double_t TIchi=0;
    for(Int_t i=0;i<nchimax;i++){
      nchi=step*i;
      Ichi[i]=lf->Integral(lb-nchi,ub-nchi);
      TIchi+=Ichi[i];
      if(Ichi[i]<TIchi/10/nchimax) break;
            if(i%10==0) cout<<i*step<<" "<<Ichi[i]<<" " <<TIchi<<endl;
    }
    Double_t runner=0;
    for(Int_t i=0;i<nchimax;i++){
      runner+=Ichi[i];
      if(runner > 0.9 * TIchi){
	anaout << "Chi likelyhood limit  found "<<i*step<<" eff="<<lim[0]->Geth1(name)->Integral()/lim[0]->mcevga<<endl;
	lim[0]->eactivity = i*step/data->totaltime;
	break;
      }
    }    
return 0;
}


