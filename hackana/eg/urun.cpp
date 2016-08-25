#include "ana.hpp"

#include "TH1.h"
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
#include "TFractionFitter.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TText.h"
#include "TLatex.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <iomanip>


////////////////////////////////
// eg version of urun.cpp     //
// version : 09.09.17         //
////////////////////////////////


Int_t Urun(bana10 *data[],bana10 *bgr[][MAXCOMP],bana10 * sig[][MAXCOMP],bana10 * lim[][MAXCOMP],Int_t nbgr[],Int_t nsig[],Int_t nlim[],Int_t nchnls, string namechnl[]){
  
  Int_t pre=2;
  Int_t w=9;
  cout<<setprecision(pre);
  anaout<<setprecision(pre);
  
  // which stuff you want done
  Bool_t itemize      = 1;
  Bool_t basic        = 1;
  Bool_t small_plots  = 1;
  Bool_t conf_plots   = 1;
  Bool_t big_plots    = 0;
  Bool_t conf_big     = 0;
  
  // fitting parameters
  //Double_t fit_min=2.7,fit_max=4.0;  // tl208 fit
  //Double_t fit_min=1.5,fit_max=4.0;  // bi214 low fit
  //Double_t fit_min=1.0,fit_max=4.0;  // bi214 high fit
  Double_t fit_min=0.4,fit_max=4.0;
  Int_t fit_rebin=4;
  
  // plotting parameters
  string ext=".pdf";
  Int_t xlen=520;
  Int_t ylen=500;
  Int_t plot_rebin=4;
  Bool_t use_fill=true;
  Int_t leg=2;
  Int_t fillstyle=3003;
  
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
  ////  ============ 2b2v Fit ==============================================
  //////////////////////////////////////////////////////////////////////////
  Double_t nexdata=0,enexdata=0,nexbgr=0,enexbgr=0,nexsig=0,enexsig=0,eff=0,runtime=0;
  if(nsig[0]>0){
    //FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
  }
  
  
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
      //cout<<setw(25)<<left<<bgr[0][i]->Getname()<<"total MC = "<<bgr[0][i]->mcevga<<endl;
    }
  }
  enexbgr=sqrt(enexbgr);
  //return 1;
  
  //---- DATA
  ////////////////////////////////////////////////////
  if(data[0]){
    runtime=(data[0]->totaltime)/86400.0;
    nexdata=data[0]->GetExpected("etot");
    enexdata=data[0]->GetErrorExpected("etot");
    anaout<<"\n\n "<<data[0]->Getname()<<endl;
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
  
  
  //---- PLOT HISTOGRAMS!!!
  /////////////////////////////////////////////////////////////
  TCanvas *c1,*c2,*c3,*c4;
  TCanvas *c801,*c802,*c803,*c804,*c805,*c806,*c807,*c808,*c809,*c810;
  
  if(itemize){
    hdraw[3]=0;hdraw[4]=0;
    TCanvas *c990 = new TCanvas("c990","bkg-breakdown",0,0,xlen,ylen);
    DrawBreakdown("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.3,5.,plot_rebin,"eg",2);
    c990->SetLogy();
    c990->Update();
    c990->Print(("Conf-internals"+ext).c_str());
    
    hdraw[2]=0;hdraw[3]=1;hdraw[4]=1;
    TCanvas *c991 = new TCanvas("c991","bkg-breakdown",xlen,0,xlen,ylen);
    DrawBreakdown("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.3,5.,plot_rebin,"eg",2);
    c991->SetLogy();
    c991->Update();
    c991->Print(("Conf-foil-n-externals"+ext).c_str());
    
    //hdraw[2]=1;hdraw[3]=0;hdraw[4]=0;
    //TCanvas *c992 = new TCanvas("c992","bkg-breakdown",0,0,xlen,ylen);
    //DrawBreakdown("cosEG",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,-1.,1.,plot_rebin,"cos",2);
    //c992->SetLogy();
    //c992->Update();
    //c992->Print(("Conf-internal-cos"+ext).c_str());
  }
  
  if(small_plots || basic){
    c1 = new TCanvas("c1","Total Energy",0,0,xlen,ylen);
    Drawall_d1("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.1,5.,plot_rebin,fillstyle,8);
    c1->Update();
    c1->Print("Etotal.png");
    
    TCanvas *c11 = new TCanvas("c11","Total Energy - log",xlen+10,0,xlen,ylen);
    Drawall_d1("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,false,false,0.1,5.,plot_rebin,fillstyle,5,9);
    c11->SetLogy();
    c11->Update();
    c11->Print("Etotal-log.png");
  }
  
  if(small_plots){
    c2 = new TCanvas("c2","Electron Energy",100,100,xlen,ylen);
    Drawall_d1("electronE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.1,3.5,plot_rebin,fillstyle,8);
    c2->Update();
    c2->Print("Eelectron.png");
    
    c3 = new TCanvas("c3","Gamma Energy",200,200,xlen,ylen);
    Drawall_d1("gammaE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.1,3.5,plot_rebin,fillstyle,8);
    c3->Update();
    c3->Print("Egamma.png");
    
    c4 = new TCanvas("c4","EG Cosine",300,300,xlen,ylen);
    Drawall_d1("cosEG",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,-1.0,1.0,plot_rebin,fillstyle,8,2);
    c4->Update();
    c4->Print("EGcosine.png");
  }// end of if(small_plots)
  
  if(conf_plots){
    c805 = new TCanvas("c805","Conf EG Total - stack",300,0,xlen,ylen);
    DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.3,5.0,plot_rebin,"eg",1,3003,1);
    c805->Update();
    c805->Print(("Conf-totE-stack"+ext).c_str());

    c810 = new TCanvas("c810","Conf EG Total - log - stack",300,0,xlen,ylen);
    DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0], 0.3 , 5.0 ,plot_rebin,"eg",8,3003,1);
    c810->SetLogy();
    c810->Update();
    c810->Print(("Conf-totE-log-stack"+ext).c_str());

    c809 = new TCanvas("c809","Conf EG - tail chi2",300,0,xlen,ylen);
    DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0], 1.5 , 4.0 ,plot_rebin,"eg",8,3003,1); // testing chi^2
    c809->SetLogy();
    c809->Update();
    c809->Print(("Conf-E-tail-chi2"+ext).c_str());

    c806 = new TCanvas("c806","Conf Electron E - stack",400,100,xlen,ylen);
    DrawBkgConf("electronE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.5,plot_rebin,"e",1,3003,1);
    c806->Update();
    c806->Print(("Conf-eleE-stack"+ext).c_str());
    
    c807 = new TCanvas("c807","Conf Gamma E - stack",500,200,xlen,ylen);
    DrawBkgConf("gammaE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.5,plot_rebin,"g",1,3003,1);
    c807->Update();
    c807->Print(("Conf-gamE-stack"+ext).c_str());
    
    c808 = new TCanvas("c808","Conf EG Cos - stack",600,300,xlen,ylen);
    DrawBkgConf("cosEG",data[0],bgr[0],nbgr[0],sig[0],nsig[0],-1.0,1.0,plot_rebin,"cos",1,3003,1);
    c808->Update();
    c808->Print(("Conf-cos-stack"+ext).c_str());
  
}// end of if(conf_plots)
  
  if(big_plots && small_plots){
    TCanvas *c77 = new TCanvas("c77","Big Canvas",0,0,xlen*2,ylen*2);
    c77->Divide(2,2,0.002,0.002,0);
    c77->cd(1);
    c1->DrawClonePad();
    c77->cd(2);
    c2->DrawClonePad();
    c77->cd(3);
    c3->DrawClonePad();
    c77->cd(4);
    c4->DrawClonePad();
    c77->Update();
    c77->Print("EG-Results.png");
  }// end of if(big_plots)
  
  if(conf_big && conf_plots){
    TCanvas *c889 = new TCanvas("c889","Big Conf Plot - stack",0,0,xlen*2,ylen*2);
    c889->Divide(2,2,0.002,0.002,0);
    c889->cd(1);
    c805->DrawClonePad();
    c889->cd(2);
    c806->DrawClonePad();
    c889->cd(3);
    c807->DrawClonePad();
    c889->cd(4);
    c808->DrawClonePad();
    c889->Update();
    c889->Print("Conf-BIG-stack.png");
    c889->Print("Conf-BIG-stack.C");
    
  }// end of if(conf_big)
  
  
  cout<<"\n\n\n\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n";
  cout<<"\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n";
  cout<<"\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n\n";
  
  return 0;
}

