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
#include "TPad.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iomanip>

////////////////////////////////
// 1e version of urun.cpp
//
// version : 09.09.17
//
// + <iomanip> in 2016-08-24
//
////////////////////////////////

Int_t Urun(bana10 *data[],bana10 *bgr[][MAXCOMP],bana10 *sig[][MAXCOMP],
	   bana10 *lim[][MAXCOMP],Int_t nbgr[],Int_t nsig[],Int_t nlim[],
	   Int_t nchnls,string namechnl[]){
  
  Int_t pre=2;
  Int_t w=9;
  cout<<setprecision(pre);
  anaout<<setprecision(pre);
  
  Bool_t itemize      = 1;
  Bool_t small_plots  = 1;
  Bool_t conf_plots   = 1;
  Bool_t big_plots    = 0;
  
  // fitting parameters
  //Double_t fit_min=2.1,fit_max=4.0;  // bi214 only fit
  //Double_t fit_min=1.5,fit_max=4.0;  // pa234m only fit
  //Double_t fit_min=1.0,fit_max=4.0;  // k40 only fit
  //Double_t fit_min=0.4,fit_max=1.5;  // bi210 and k40 fit
  Double_t fit_min=0.4,fit_max=4.0;  // bi210 only fit
  Int_t fit_rebin=2;
  
  // plotting parameters
  string ext=".pdf";
  Int_t xlen=520;
  Int_t ylen=500;
  Int_t plot_rebin=2;
  Bool_t use_fill=false;
  Int_t fillstyle=3003;
  Int_t leg=2;
  
  Int_t hdraw[6]={1,1,1,1,1,1};
  // hdraw[0] = draw data
  // hdraw[1] = draw signals
  // hdraw[2] = draw my foil internals
  // hdraw[3] = draw other foil bkgs
  // hdraw[4] = draw external bkgs
  // hdraw[5] = draw total MC
    
  vector<string> foilbkgs;
  foilbkgs.push_back("ca48");
  foilbkgs.push_back("nd150");
  
  vector<string> vetobkgs;
  vetobkgs.push_back("2b2n");
  vetobkgs.push_back("pb211");
  vetobkgs.push_back("tl207");
  vetobkgs.push_back("zr96.zr96");
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  Das Fit  ==============================================
  //////////////////////////////////////////////////////////////////////////
  Double_t nexdata,enexdata,nexbgr,enexbgr,nexsig,enexsig,eff,runtime;
  if(nsig[0]>0){
    //FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
  }
  
  nexdata=0;
  enexdata=0;
  nexbgr=0;
  enexbgr=0;
  nexsig=0;
  enexsig=0;
  eff=0;
  
  //---- BACKGROUNDS
  ////////////////////////////////////////////////////
  if(nbgr[0]>0){
    anaout<<"\n BACKGROUNDS \n";
    anaout<<"-----------------------------------------------------------------\n";
    for(Int_t i=0;i<nbgr[0];i++){
      eff=bgr[0][i]->GetExpected("etot")/bgr[0][i]->totaltime/bgr[0][i]->activity;
      anaout<<setprecision(pre);
      anaout<<scientific<<setw(25)<<left<<bgr[0][i]->Getname()<<" A = "<<setw(w)<<bgr[0][i]->norm<<" Eff = "<<setw(w)<<eff<<" Nexp = "<<setw(w)<<bgr[0][i]->GetExpected("etot")<<" +/- "<<setw(w)<<bgr[0][i]->GetErrorExpected("etot")<<endl;
      nexbgr += bgr[0][i]->GetExpected("etot");
      enexbgr += (bgr[0][i]->GetErrorExpected("etot")) * (bgr[0][i]->GetErrorExpected("etot"));
    }
  }
  enexbgr=sqrt(enexbgr);
  
  //---- DATA
  ////////////////////////////////////////////////////
  if(data[0]){
    runtime=(data[0]->totaltime)/86400.0;
    nexdata=data[0]->GetExpected("etot");
    enexdata=data[0]->GetErrorExpected("etot");
    anaout<<"\n\n "<<data[0]->Getname()<<"\n";
    anaout<<"-----------------------------------------------------------------\n";
    anaout<<fixed<<"Runtime\t= "<<runtime<<" days"<<endl;
    anaout<<fixed<<"Events \t= "<<nexdata<<endl;
  }
  
  //---- SIGNALS
  ////////////////////////////////////////////////////
  if(nsig[0]>0){
    anaout<<"\n\n SIGNALS \n";
    anaout<<"-----------------------------------------------------------------\n";
    for(Int_t i=0; i<nsig[0]; i++){
      eff=sig[0][i]->GetExpected("etot")/sig[0][i]->totaltime/sig[0][i]->activity;
      anaout<<setprecision(pre);
      anaout<<scientific<<setw(25)<<left<<sig[0][i]->Getname()<<" A = "<<setw(w)<<sig[0][i]->norm<<" +/- "<<setw(w)<<sig[0][i]->enorm<<" Eff = "<<setw(w)<<eff<<" Nexp = "<<setw(w)<<sig[0][i]->GetExpected("etot")<<" +/- "<<setw(w)<<sig[0][i]->GetErrorExpected("etot")<<endl;
      nexsig += sig[0][i]->GetExpected("etot");
      enexsig += (sig[0][i]->GetErrorExpected("etot")) * (sig[0][i]->GetErrorExpected("etot"));
    }
    anaout<<setprecision(1);
    anaout<<fixed<<"\nFit window = "<<fit_min<<" - "<<fit_max<<" MeV "<<endl;
  }  
  enexsig=sqrt(enexsig);
  
  
  //---- calculate and print out some info
  ////////////////////////////////////////////////////////////
  runtime=(data[0]->totaltime)/86400.0;
  cout<<"\n\n\n============================================================\n";
  cout<<"\n ------- "<<data[0]->Getname()<<" ------- "<<endl;
  cout<<fixed<<"Runtime\t= "<<runtime<<" days \n";
  cout<<fixed<<"Events \t= "<<nexdata<<" \n\n";
  if(nsig[0]>0){
    for(Int_t i=0; i<nsig[0]; i++){
      cout<<setprecision(pre);
      cout<<scientific<<setw(25)<<left<<sig[0][i]->Getname()<<" A = "<<sig[0][i]->norm<<"  0  BR E "<<sig[0][i]->enorm<<endl;
    }
    cout<<setprecision(1);
    cout<<fixed<<"\nFit window = "<<fit_min<<" - "<<fit_max<<" MeV "<<endl;
  }
  cout<<"\n============================================================\n\n\n";
  
  
  //---- plot the E-total histogram!
  ////////////////////////////////////////////////////////////
  
  if(itemize){  
    hdraw[3]=0;hdraw[4]=0;
    TCanvas *c990 = new TCanvas("c990","bkg-breakdown",0,0,xlen,ylen);
    DrawBreakdown("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.3,3.5,plot_rebin,"e",5);
    c990->SetLogy();
    c990->Update();
    c990->Print(("Conf-internals"+ext).c_str());
    
    hdraw[2]=0;hdraw[3]=1;hdraw[4]=1;
    TCanvas *c991 = new TCanvas("c991","bkg-breakdown",xlen,0,xlen,ylen);
    DrawBreakdown("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.3,3.5,plot_rebin,"e",5);
    c991->SetLogy();
    c991->Update();
    c991->Print(("Conf-foil-n-externals"+ext).c_str());
  }
  
  TCanvas *c501,*c502,*c801,*c802,*c803,*c804,*c805;
  if(small_plots){
    c501 = new TCanvas("c501","Electron Energy",0,0,xlen,ylen);
    Drawall_d1("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.4,1.4,plot_rebin,fillstyle,8,1.2);
    c501->Update();
    c501->Print(("Etotal"+ext).c_str());
    
    c502 = new TCanvas("c502","Electron Energy - log",xlen,0,xlen,ylen);
    Drawall_d1("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.4,3.5,plot_rebin,fillstyle,8,8);
    c502->SetLogy();
    c502->Update();
    c502->Print(("Etotal-log"+ext).c_str());
  }
  
  if(conf_plots){
    c803 = new TCanvas("c803","Conf E - stack",100,100,xlen,ylen);
    DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.4,1.4,plot_rebin,"e",1,3003,1);
    c803->Update();
    c803->Print(("Conf-E-stack"+ext).c_str());
    
    c804 = new TCanvas("c804","Conf E - log stack",xlen+100,100,xlen,ylen);
    DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.4,5.0,plot_rebin,"e",5,3003,1);
    c804->SetLogy();
    c804->Update();
    c804->Print(("Conf-E-log-stack"+ext).c_str());
    
    c805 = new TCanvas("c805","Conf E - tail chi2",xlen+200,100,xlen,ylen);
    DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],1.5,3.5,plot_rebin,"e",5,3003,1);
    c805->SetLogy();
    c805->Update();
    c805->Print(("Conf-E-tail-chi2"+ext).c_str());
    
  }
  
  if(big_plots && conf_plots){
    TCanvas *c555 = new TCanvas("c555","Analysis Results",0,0,xlen,ylen*2);
    c555->Divide(1,2,0.002,0.002,0);
    c555->cd(1);
    c801->DrawClonePad();
    c555->cd(2);
    c802->DrawClonePad();
    c555->cd(2)->SetLogy();
    c555->Update();
    c555->Print(("E-RESULTS"+ext).c_str());
  }
  
  cout<<"\n\n\n\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n";
  cout<<"\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n";
  cout<<"\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n\n"; 
  return 0;
}

