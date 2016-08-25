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


//////////////////////////////////////
// ee COMBINED version of urun.cpp  //
// version : 09.09.17               //
//////////////////////////////////////


Int_t Urun(bana10 *data[],bana10 *bgr[][MAXCOMP],bana10 * sig[][MAXCOMP],
	   bana10 * lim[][MAXCOMP],Int_t nbgr[],Int_t nsig[],Int_t nlim[],
	   Int_t nchnls, string namechnl[]){
  
  string isotope="^{96}Zr";
  string decaymode="m1";
  
  // plotting parameters
  string ext=".pdf";
  Int_t xlen=520;
  Int_t ylen=500;
  Int_t plot_rebin=5;
  Bool_t use_fill=true;
  Int_t fillstyle=3003;
  Int_t leg=2;
  Int_t grey=0;
  
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
  vetobkgs.push_back("2b0n");
  vetobkgs.push_back("2b2n");
  vetobkgs.push_back("pb211");
  vetobkgs.push_back("tl207");
  vetobkgs.push_back("zr96.zr96");
  
  Int_t pre = 2;
  Int_t w   = 9;
  cout<<setprecision(pre);
  anaout<<setprecision(pre);
  
  // which source are you analyzing
  Int_t source=3;  // itep=1  inr=2  combined=3
  
  // which stuff you want done
  //////////////////////////////////////////////////////////
  Bool_t basic_anal    = 1;  // does the basic analysis
  Bool_t basic_plots   = 0;  // just show the ee lin-log plots
  Bool_t itemize       = 1;  // itemize the data, ints, foils, exts
  Bool_t gen_collie    = 0;  // creates histograms for COLLIE
  Bool_t print_syst    = 0;  // print the systematic breakdown
  Bool_t show_results  = 0;  // writes up the results on a pad
  Bool_t norm_plots    = 0;  // plots show main 10 bkgs
  Bool_t add_plots     = 1;  // blessed plots show bkgs added
  Bool_t sub_plots     = 1;  // blessed plots show bkgs subtracted
  Bool_t conf_plots    = 0;  // plots show int and ext bkgs
  Bool_t norm_big      = 0;  // same as above but all combined
  Bool_t add_big       = 0;  // same as above but all combined
  Bool_t sub_big       = 0;  // same as above but all combined
  Bool_t conf_big      = 0;  // same as above but all combined
  Bool_t bkg_syst      = 0;  // scales backgrounds to observe systematics
  Bool_t ewin_syst     = 0;  // change energy window to observe systematics
  Bool_t mclimit       = 1;  // finds the 0v limit with TLimit and Helene
  Bool_t multi_tlimit  = 0;  // changes binning for TLimit method
  Bool_t multi_helene  = 0;  // changes E window for Helene method
  Bool_t multi_LFit    = 0;  // changes binning and E window for likelihood
  Bool_t multi_Fit     = 0;  // changes binning and E window for binned fit
  
  if(0){
    basic_anal    = 1;
    basic_plots   = 1;
    itemize       = 1;
    gen_collie    = 1;
    print_syst    = 1;
    show_results  = 1;
    norm_plots    = 1;
    add_plots     = 1;
    sub_plots     = 1;
    conf_plots    = 1;
    norm_big      = 1;
    add_big       = 1;
    sub_big       = 1;
    conf_big      = 1;
    bkg_syst      = 1;
    ewin_syst     = 1;
    mclimit       = 1;
  }
  
  if(gen_collie) mclimit=1;
  
  // fitting parameters
  Double_t fit_min=0.6;
  Double_t fit_max=4.0;
  Int_t fit_rebin=4;
  
  // background systematics
  const Double_t scalebkgs=0.10;
  
  // multi_fit parameters
  Int_t rebin_max=6;
  Int_t ewin_max=20;
  
  // some 2vBB calculation constants
  const Double_t Anum     = 96.0;
  const Double_t Avogadro = 6.0221415e23;
  const Double_t lntwo    = TMath::Log(2);
  const Double_t time     = (60.*60.*24.*365.242);
  Double_t mass=0;
  if(source==1) mass      = 4.10;  //---for ITEP source
  if(source==2) mass      = 5.30;  //---for INR source
  if(source==3) mass      = 9.40;  //---COMBINED
  const Double_t Smass    = mass;
  const Double_t N_0      = (Smass/Anum)*Avogadro;
  const Double_t Nconst   = (N_0*lntwo)/time;
  //const Double_t deadtime = 0.015; // for phase 1+2
  const Double_t deadtime = 0.0201; // for phase 1+2+3
  const Double_t G_2v     = 1.8e-17;
  Double_t nexdata,enexdata,nexbgr,enexbgr,nexsig,enexsig,eff,runtime,runtime_corr;
  Double_t activity,errActiv,halflife,halflife_corr,herror,hsystp,hsystm,sig2bkg,sigmas;
  Double_t T_2v,T_2v_stat,T_2v_syst,T_2v_err,M_2v,M_2v_err;
  
  // some 0vBB calculation constants
  Double_t T_0v=0,G_0v=0;
  Double_t Me=510998.9;
  
  Double_t G_0v_Sim  = 7.36e-14;
  Double_t G_0v_Suh  = 5.70e-14;
  
  Double_t G_0vM_Sim = 1.59e-15;
  Double_t G_0vM_Suh = 1.24e-15;
  
  Double_t RQRPA_ga125,RQRPA_ga100,QRPA_ga125,QRPA_ga100;
  Double_t RQRPA_Jast_ga125,RQRPA_Jast_ga100,RQRPA_UCOM_ga125,RQRPA_UCOM_ga100;
  Double_t QRPA_Jast,QRPA_UCOM,QRPA_UCOM_ga125,QRPA_UCOM_ga100;
  
  Double_t T_helene,T_helene_corr,T_helene_eff,T_mclim,T_mclim_corr,T_mclim_eff;
  Double_t limEvents=-1;
  
  // calculate the systematics
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  const Int_t NUM=7;
  Double_t syst_p[NUM],syst_m[NUM];
  string syst_name[NUM];
  Int_t s=0;
  
  syst_name[s]="Aparatus";  
  syst_p[s]=0.0500;  // aparatus (+)
  syst_m[s]=0.0500;  // aparatus (-)
  s++;
  
  syst_name[s]="Calibration";  
  syst_p[s]=0.0293;  // calib +1% (+)
  syst_m[s]=0.0216;  // calib -1% (-)
  s++;
  
  syst_name[s]="Mass";  
  syst_p[s]=0.0200;  // mass (+)
  syst_m[s]=0.0200;  // mass (-)
  s++;
  
  syst_name[s]="Ext Bkgs +/-10%";
  syst_p[s]=0.003;  // bkgs -10% (+)
  syst_m[s]=0.003;  // bkgs +10% (-)
  s++;
  
  syst_name[s]="Nd150 +/-10%";
  syst_p[s]=0.007;  // nd150 -10% (+)
  syst_m[s]=0.007;  // nd150 +10% (-)
  s++;
  
  syst_name[s]="Int Bkgs +/-5%";
  syst_p[s]=0.019;  // bkgs -5% (+)
  syst_m[s]=0.019;  // bkgs +5% (-)
  s++;
  /*
  syst_name[s]="Pa234m";
  syst_p[s]=0.12;
  syst_m[s]=0.12;
  s++;
  
  syst_name[s]="K40";
  syst_p[s]=0.034;
  syst_m[s]=0.034;
  s++;
  
  syst_name[s]="U235";
  syst_p[s]=0.002;
  syst_m[s]=0.002;
  s++;
  
  syst_name[s]="U238";
  syst_p[s]=0.007;
  syst_m[s]=0.007;
  s++;
  
  syst_name[s]="Th232";
  syst_p[s]=0.0001;
  syst_m[s]=0.0001;
  s++;
  */
  syst_name[s]="E window";  
  syst_p[s]=0.0155;  // likelihood E window (+)
  syst_m[s]=0.0023;  // likelihood E window (-)
  s++;
  
  Double_t total_syst_p=0,total_syst_m=0;
  for(Int_t s=0;s<NUM;s++){
    total_syst_p+=TMath::Power(syst_p[s],2);
    total_syst_m+=TMath::Power(syst_m[s],2);
  }
  total_syst_p=TMath::Sqrt(total_syst_p);
  total_syst_m=TMath::Sqrt(total_syst_m);
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  
  TCanvas *c401,*c402,*c403,*c404,*c405,*c406,*c444;
  TCanvas *c501,*c502,*c503,*c504,*c505,*c506,*c555;
  TCanvas *c601,*c602,*c603,*c604,*c605,*c606,*c666;
  TCanvas *c701,*c702,*c703,*c704,*c705,*c706,*c777;
  TCanvas *c801,*c802,*c803,*c804,*c805,*c806,*c888;
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  2vBB FITTING  =========================================
  //////////////////////////////////////////////////////////////////////////
  if(basic_anal){
    if(nsig[0]>0){
      //FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
      LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
      //LFitSignal("cosa",data[0],bgr[0],nbgr[0],sig[0],nsig[0],-1,1,fit_rebin);
    }
    
    nexdata=enexdata=nexbgr=enexbgr=nexsig=enexsig=eff=runtime=runtime_corr=0;
    
    //---- BACKGROUNDS
    ////////////////////////////////////////////////////
    if(nbgr[0]>0){
      anaout<<"\n BACKGROUNDS \n";
      anaout<<"=================================================================="<<endl;
      for(Int_t i=0;i<nbgr[0];i++){
	eff=bgr[0][i]->GetExpected("etot")/bgr[0][i]->totaltime/bgr[0][i]->activity;
	anaout<<setprecision(pre);
	anaout<<scientific<<setw(25)<<left<<bgr[0][i]->Getname()<<" A = "<<setw(w)<<bgr[0][i]->norm<<" Eff = "<<setw(w)<<eff<<" Nexp = "<<setw(w)<<bgr[0][i]->GetExpected("etot")<<"+/- "<<setw(w)<<bgr[0][i]->GetErrorExpected("etot")<<endl;
	// PRINT OUT ACTIVITIES TO RE-ENTER INTO CONTROL FILE
	//anaout<<"vt_externals-p_.root    "<<setw(24)<<left<<bgr[0][i]->Getname()<<scientific<<setw(8)<<bgr[0][i]->norm<<"  0  BR "<<endl;
	nexbgr += bgr[0][i]->GetExpected("etot");
	enexbgr += (bgr[0][i]->GetErrorExpected("etot")) * (bgr[0][i]->GetErrorExpected("etot"));
      }
    }
    enexbgr=TMath::Sqrt(enexbgr);
    //enexbgr=TMath::Sqrt(nexbgr);
    
    //---- DATA
    ////////////////////////////////////////////////////
    if(data[0]){
      runtime=(data[0]->totaltime)/86400.0;
      runtime_corr=runtime*(1.0-deadtime);
      nexdata=data[0]->GetExpected("etot");
      enexdata=data[0]->GetErrorExpected("etot");
      //enexdata=TMath::Sqrt(nexdata);
      anaout<<"\n\n "<<data[0]->Getname()<<endl;
      anaout<<"=================================================================="<<endl;
      anaout<<fixed<<"Runtime     = "<<runtime<<" days"<<endl;
      anaout<<fixed<<"Corr. Time  = "<<runtime_corr<<" days"<<endl;
      anaout<<fixed<<"Events      = "<<setprecision(0)<<nexdata<<endl;
    }
    
    //---- SIGNALS
    ////////////////////////////////////////////////////
    if(nsig[0]>0){
      anaout<<"\n\n SIGNALS \n";
      anaout<<"=================================================================="<<endl;
      for(Int_t i=0; i<nsig[0]; i++){
	//eff=sig[0][i]->GetExpected("etot")/(sig[0][i]->totaltime*sig[0][i]->activity);
	eff=sig[0][i]->GetExpected("etot")/(data[0]->totaltime*sig[0][i]->activity);
	anaout<<setprecision(pre);
	anaout<<scientific<<setw(25)<<left<<sig[0][i]->Getname()<<" A = "<<setw(w)<<sig[0][i]->norm<<"+/- "<<setw(w)<<sig[0][i]->enorm<<" Eff = "<<setw(w)<<eff<<" Nexp = "<<setw(w)<<sig[0][i]->GetExpected("etot")<<"+/- "<<setw(w)<<sig[0][i]->GetErrorExpected("etot")<<endl;
	nexsig += sig[0][i]->GetExpected("etot");
	enexsig += (sig[0][i]->GetErrorExpected("etot")) * (sig[0][i]->GetErrorExpected("etot"));
      }
    }  
    enexsig=TMath::Sqrt(enexsig);
    //enexsig=TMath::Sqrt(nexsig);
    
    
    //---- calculate and print out some info
    ////////////////////////////////////////////////////////////
    //runtime=(data[0]->totaltime)/86400.0;
    //runtime_corr=runtime*(1.0-deadtime);
    cout<<"\n\n\n==================================================================\n";
    cout<<"\n ------- "<<data[0]->Getname()<<" ------- "<<endl;
    cout<<fixed<<"Runtime     = "<<setprecision(2)<<runtime<<" days \n";
    cout<<fixed<<"Events      = "<<setprecision(0)<<nexdata<<" \n\n";
    if(nsig[0]>0){
      for(Int_t i=0; i<nsig[0]; i++){
	cout<<setprecision(pre);
	cout<<scientific<<setw(25)<<left<<sig[0][i]->Getname()<<" A = "<<sig[0][i]->norm<<"  0  BR E "<<sig[0][i]->enorm<<endl;
      }
    }
    cout<<"\n==================================================================\n\n\n";
    
    if(nsig[0]>0){
      activity=sig[0][0]->norm;
      errActiv=sig[0][0]->enorm;
      halflife=Nconst/activity;
      //halflife=lntwo*(runtime/365)*(1/(TMath::Log(N_0-nexsig)-TMath::Log(N_0)));
      //halflife=lntwo*(runtime/365)*N_0*(1/(N_0-nexsig));
      halflife_corr=halflife*(1.0-deadtime);
      herror=(Nconst/(activity*activity))*(errActiv)*(1.0-deadtime);
      
      hsystp=halflife_corr*total_syst_p;
      hsystm=halflife_corr*total_syst_m;
      sig2bkg=nexsig/nexbgr;
      sigmas=TMath::Abs(nexdata-nexbgr)/(TMath::Sqrt(enexdata*enexdata+enexbgr*enexbgr));
      eff=sig[0][0]->GetExpected("etot")/sig[0][0]->totaltime/sig[0][0]->activity;
      //eff_corr=eff*(1.0+deadtime);
      
      // calculate NME for 2vBB decay
      T_2v=halflife_corr;
      T_2v_stat=herror;
      T_2v_syst=hsystp;
      T_2v_err=TMath::Sqrt((T_2v_stat*T_2v_stat)+(T_2v_syst*T_2v_syst));
      M_2v=TMath::Sqrt(1./(T_2v*G_2v));
      M_2v_err=TMath::Abs(M_2v-(TMath::Sqrt(1./((T_2v-T_2v_err)*G_2v))));
      
      
      cout<<endl;
      cout<<"Runtime      = "<<fixed<<setprecision(1)<<runtime<<" days \n";
      cout<<"Deadtime     = "<<fixed<<setprecision(2)<<deadtime*100<<" % \n";
      cout<<"Corr. Time   = "<<fixed<<setprecision(1)<<runtime_corr<<" days"<<endl;
      cout<<endl;
      cout<<"Data         = "<<fixed<<setprecision(0)<<nexdata<<endl;
      cout<<"Bgr          = "<<fixed<<setprecision(1)<<nexbgr<<" +/- "<<setprecision(1)<<enexbgr<<endl;
      cout<<"Sig          = "<<fixed<<setprecision(1)<<nexsig<<" +/- "<<setprecision(1)<<enexsig<<endl;
      cout<<"S/B          = "<<fixed<<setprecision(2)<<sig2bkg<<endl;
      cout<<"Eff          = "<<fixed<<setprecision(1)<<eff*100<<" % "<<endl;
      //cout<<"Corr. Eff    = "<<fixed<<setprecision(1)<<eff_corr*100<<" % "<<endl;
      cout<<endl;
      cout<<"96Zr Mass    = "<<fixed<<setprecision(2)<<Smass<<" gr. "<<endl;
      cout<<"Masstime     = "<<fixed<<setprecision(3)<<(Smass*runtime_corr)/1000/365<<" kg.yr "<<endl;
      cout<<"Activity     = "<<scientific<<setprecision(2)<<activity<<" +/- "<<setprecision(2)<<errActiv<<" Bq. "<<endl;
      cout<<"Syst Error   = (+) "<<fixed<<setprecision(1)<<(total_syst_p*100)<<"%  (-) "<<(total_syst_m*100)<<"% "<<endl;
      cout<<"Half-life    = "<<scientific<<setprecision(4)<<halflife<<" yr. "<<endl;
      cout<<"Corr. H-L    = "<<scientific<<setprecision(4)<<halflife_corr<<" yr. "<<endl;
      cout<<endl;
      cout<<"Final Result = "<<scientific<<setprecision(2)<<halflife_corr<<" +/- "<<setprecision(2)<<herror<<"(stat) + "
	  <<setprecision(2)<<hsystp<<"(syst) - "<<setprecision(2)<<hsystm<<"(syst) yr. "<<endl;
      cout<<"Sigmas       = "<<fixed<<setprecision(1)<<sigmas<<endl;
      cout<<"2vBB NME     = "<<fixed<<setprecision(4)<<M_2v<<" +/- "<<M_2v_err<<endl;
      cout<<endl;
      cout<<"Fit rebinning  ==>  "<<fixed<<setprecision(0)<<fit_rebin<<endl;
      cout<<"Energy Window  ==>  "<<fixed<<setprecision(1)<<fit_min<<" - "<<fit_max<<" MeV "<<endl;
      cout<<"\n==================================================================\n\n";
      
      //---- write to the ana10.dat file!!
      anaout<<endl;
      anaout<<"Runtime      = "<<fixed<<setprecision(1)<<runtime<<" days \n";
      anaout<<"Deadtime     = "<<fixed<<setprecision(2)<<deadtime*100<<" % \n";
      anaout<<"Corr. Time   = "<<fixed<<setprecision(1)<<runtime_corr<<" days"<<endl;
      anaout<<endl;
      anaout<<"Data         = "<<fixed<<setprecision(0)<<nexdata<<endl;
      anaout<<"Bgr          = "<<fixed<<setprecision(1)<<nexbgr<<" +/- "<<setprecision(1)<<enexbgr<<endl;
      anaout<<"Sig          = "<<fixed<<setprecision(1)<<nexsig<<" +/- "<<setprecision(1)<<enexsig<<endl;
      anaout<<"S/B          = "<<fixed<<setprecision(2)<<sig2bkg<<endl;
      anaout<<"Eff          = "<<fixed<<setprecision(1)<<eff*100<<" % "<<endl;
      //anaout<<"Corr. Eff    = "<<fixed<<setprecision(1)<<eff_corr*100<<" % "<<endl;
      anaout<<endl;
      anaout<<"96Zr Mass    = "<<fixed<<setprecision(2)<<Smass<<" gr. "<<endl;
      anaout<<"Masstime     = "<<fixed<<setprecision(3)<<(Smass*runtime_corr)/1000/365<<" kg.yr "<<endl;
      anaout<<"Activity     = "<<scientific<<setprecision(2)<<activity<<" +/- "<<setprecision(2)<<errActiv<<" Bq. "<<endl;
      anaout<<"Syst Error   = (+) "<<fixed<<setprecision(1)<<(total_syst_p*100)<<"%  (-) "<<(total_syst_m*100)<<"% "<<endl;
      anaout<<"Half-life    = "<<scientific<<setprecision(4)<<halflife<<" yr. "<<endl;
      anaout<<"Corr. H-L    = "<<scientific<<setprecision(4)<<halflife_corr<<" yr. "<<endl;
      anaout<<endl;
      anaout<<"Final Result = "<<scientific<<setprecision(2)<<halflife_corr<<" +/- "<<setprecision(2)<<herror<<"(stat) + "
	  <<setprecision(2)<<hsystp<<"(syst) - "<<setprecision(2)<<hsystm<<"(syst) yr. "<<endl;
      anaout<<"Sigmas       = "<<fixed<<setprecision(1)<<sigmas<<endl;
      anaout<<"2vBB NME     = "<<fixed<<setprecision(4)<<M_2v<<" +/- "<<M_2v_err<<endl;
      anaout<<endl;
      anaout<<"Fit rebinning  ==>  "<<fixed<<setprecision(0)<<fit_rebin<<endl;
      anaout<<"Energy Window  ==>  "<<fixed<<setprecision(1)<<fit_min<<" - "<<fit_max<<" MeV "<<endl;
      
    }// end of if(nsig[0]>0)
    
    
    if(0){
      hdraw[3]=0;hdraw[4]=0;
      
      int Dbreak=0;
      int Dconf=1;
      int fill=0; //3003
      
      c401 = new TCanvas("c401","Delta R",0,0,xlen,ylen);
      if(Dbreak) DrawBreakdown("deltaR",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.,5.,2,"cm");
      if(Dconf)  DrawBkgConf("deltaR",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.,5.,2,"|#DeltaR| (cm)",2.,fill,0);
      c401->SetLogy();
      c401->Update();
      c401->Print(("Delta-R"+ext).c_str());
      
      c402 = new TCanvas("c402","Delta Z",0,0,xlen,ylen);
      if(Dbreak) DrawBreakdown("deltaZ",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.,10.,2,"cm");
      if(Dconf)  DrawBkgConf("deltaZ",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.,10.,2,"|#DeltaZ| (cm)",2.,fill,0);
      c402->SetLogy();
      c402->Update();
      c402->Print(("Delta-Z"+ext).c_str());
      
      c403 = new TCanvas("c403","Track Length",0,0,xlen,ylen);
      if(Dbreak) DrawBreakdown("trkLen",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.,200.,2,"cm");
      if(Dconf)  DrawBkgConf("trkLen",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.,200.,2,"Track Length (cm)",20.,fill,0);
      c403->SetLogy();
      c403->Update();
      c403->Print(("Trk-Length"+ext).c_str());
      
      c404 = new TCanvas("c404","Int Hyp Prob",0,0,xlen,ylen);
      if(Dbreak) DrawBreakdown("pint",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.,1.,2);
      if(Dconf)  DrawBkgConf("pint",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.,1.,2,"Prob(int)",2.5,fill,0);
      c404->SetLogy();
      c404->Update();
      c404->Print(("Prob-Int"+ext).c_str());
      
      c405 = new TCanvas("c405","Ext Hyp Prob",0,0,xlen,ylen);
      if(Dbreak) DrawBreakdown("pext",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.,1.,2);
      if(Dconf)  DrawBkgConf("pext",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.,1.,2,"Prob(ext)",2.,fill,0);
      c405->SetLogy();
      c405->Update();
      c405->Print(("Prob-Ext"+ext).c_str());
      
      
    }

    // tons of plots!!
    if(norm_plots || basic_plots){
      c501 = new TCanvas("c501","EE Total",0,0,xlen,ylen);
      Drawall_d1("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.1,4.0,plot_rebin,fillstyle,5,1.2);
      c501->Update();
      c501->Print(("EE-total"+ext).c_str());
      
      c505 = new TCanvas("c505","EE Total - log",xlen+10,0,xlen,ylen);
      Drawall_d1("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,false,false,0.1,4.0,plot_rebin,fillstyle,5,8);
      c505->SetLogy();
      c505->Update();
      c505->Print(("EE-total-log"+ext).c_str());
    }
    if(norm_plots){
      c502 = new TCanvas("c502","E Min",100,100,xlen,ylen);
      Drawall_d1("minE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.1,2.0,plot_rebin-1,fillstyle,5);
      c502->Update();
      c502->Print(("EE-min"+ext).c_str());
      
      c503 = new TCanvas("c503","E Single",200,200,xlen,ylen);
      Drawall_d1("singleE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,0.1,3.0,plot_rebin,fillstyle,5);
      c503->Update();
      c503->Print(("EE-single"+ext).c_str());
      
      c504 = new TCanvas("c504","EE Cos",300,300,xlen,ylen);
      Drawall_d1("cosa",data[0],bgr[0],nbgr[0],sig[0],nsig[0],leg,use_fill,false,-1.0,1.0,plot_rebin,fillstyle,5,1.5);
      c504->Update();
      c504->Print(("EE-cos"+ext).c_str());
    }// end of if(norm_plots)
    
    if(itemize){
      hdraw[3]=0;hdraw[4]=0;
      TCanvas *c990 = new TCanvas("c990","bkg-breakdown",0,0,xlen,ylen);
      DrawBreakdown("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.0,5.,plot_rebin,"ee",2);
      c990->SetLogy();
      c990->Update();
      c990->Print(("Conf-internals"+ext).c_str());
      
      hdraw[2]=0;hdraw[3]=1;hdraw[4]=1;
      TCanvas *c991 = new TCanvas("c991","bkg-breakdown",xlen,0,xlen,ylen);
      DrawBreakdown("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],hdraw,foilbkgs,vetobkgs,0.0,5.,plot_rebin,"ee",2);
      c991->SetLogy();
      c991->Update();
      c991->Print(("Conf-foil-n-externals"+ext).c_str());
    }
    
    if(sub_plots){
      c601 = new TCanvas("c601","totE bkg-sub",100,0,xlen,ylen);
      DrawBkgSub("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,4.0,plot_rebin,"ee",isotope);
      c601->Update();
      c601->Print(("BKG-sub-totE"+ext).c_str());
      
      c602 = new TCanvas("c602","minE bkg-sub",200,100,xlen,ylen);
      DrawBkgSub("minE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",isotope);
      c602->Update();
      c602->Print(("BKG-sub-minE"+ext).c_str());
      
      c603 = new TCanvas("c603","singE bkg-sub",300,200,xlen,ylen);
      DrawBkgSub("singleE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",isotope);
      c603->Update();
      c603->Print(("BKG-sub-singE"+ext).c_str());
      
      c604 = new TCanvas("c604","cos bkg-sub",400,300,xlen,ylen);
      DrawBkgSub("cosa",data[0],bgr[0],nbgr[0],sig[0],nsig[0],-1.0,1.0,plot_rebin,"cos",isotope);
      c604->Update();
      c604->Print(("BKG-sub-cos"+ext).c_str());
    }// end of if(sub_plots)
    
    if(add_plots){
      c701 = new TCanvas("c701","totE bkg-add",200,0,xlen,ylen);
      DrawBkgAdd("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,4.0,plot_rebin,"ee",isotope);
      c701->Update();
      c701->Print(("BKG-add-totE"+ext).c_str());
      
      c702 = new TCanvas("c702","minE bkg-add",300,100,xlen,ylen);
      DrawBkgAdd("minE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",isotope);
      c702->Update();
      c702->Print(("BKG-add-minE"+ext).c_str());
      
      c703 = new TCanvas("c703","singE bkg-add",400,200,xlen,ylen);
      DrawBkgAdd("singleE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",isotope);
      c703->Update();
      c703->Print(("BKG-add-singE"+ext).c_str());
      
      c704 = new TCanvas("c704","cos bkg-add",500,300,xlen,ylen);
      DrawBkgAdd("cosa",data[0],bgr[0],nbgr[0],sig[0],nsig[0],-1.0,1.0,plot_rebin,"cos",isotope);
      c704->Update();
      c704->Print(("BKG-add-cos"+ext).c_str());
      
      c705 = new TCanvas("c705","maxE bkg-add",300,100,xlen,ylen);
      DrawBkgAdd("maxE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",isotope);
      c705->Update();
      c705->Print(("BKG-add-maxE"+ext).c_str());
      
      //c706 = new TCanvas("c706","tof bkg-add",300,100,xlen,ylen);
      //DrawBkgAdd("tof2",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,1.0,plot_rebin,"cos",isotope);
      //c706->Update();
      //c706->Print("BKG-add-tof"+ext.c_str());
      
      
    }// end of if(add_plots)
    
    if(conf_plots){
      c801 = new TCanvas("c801","Conf EE Total",300,0,xlen,ylen);
      DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,4.0,plot_rebin,"ee",1,3003,1,isotope);
      c801->Update();
      c801->Print(("Conf-totE"+ext).c_str());
      
      c805 = new TCanvas("c805","Conf EE Total - log",xlen+310,0,xlen,ylen);
      DrawBkgConf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,4.0,plot_rebin,"ee",6,3003,1,isotope);
      c805->SetLogy();
      c805->Update();
      c805->Print(("Conf-totE-log"+ext).c_str());
      
      c802 = new TCanvas("c802","Conf E Min",400,100,xlen,ylen);
      DrawBkgConf("minE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",1,3003,1,isotope);
      c802->Update();
      c802->Print(("Conf-minE"+ext).c_str());
      
      c803 = new TCanvas("c803","Conf E Single",500,200,xlen,ylen);
      DrawBkgConf("singleE",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.0,3.0,plot_rebin,"e",1,3003,1,isotope);
      c803->Update();
      c803->Print(("Conf-singE"+ext).c_str());
      
      c804 = new TCanvas("c804","Conf EE Cos",600,300,xlen,ylen);
      DrawBkgConf("cosa",data[0],bgr[0],nbgr[0],sig[0],nsig[0],-1.0,1.0,plot_rebin,"cos",1,3003,1,isotope);
      c804->Update();
      c804->Print(("Conf-cos"+ext).c_str());
      
    }// end of if(conf_plots)
    
    if(norm_big && norm_plots){
      c555 = new TCanvas("c555","Anal Results",0,0,xlen*2,ylen*2);
      c555->Divide(2,2,0.002,0.002,0);
      c555->cd(1);
      c501->DrawClonePad();
      c555->cd(2);
      c502->DrawClonePad();
      c555->cd(3);
      c503->DrawClonePad();
      c555->cd(4);
      c504->DrawClonePad();
      c555->Update();
      c555->Print(("EE-RESULTS"+ext).c_str());
    }// end of if(norm_big)
    
    if(sub_big && sub_plots){
      c666 = new TCanvas("c666","BKG-Subtracted Results",0,0,xlen*2,ylen*2);
      c666->Divide(2,2,0.002,0.002,0);
      c666->cd(1);
      c601->DrawClonePad();
      c666->cd(2);
      c602->DrawClonePad();
      c666->cd(3);
      c603->DrawClonePad();
      c666->cd(4);
      c604->DrawClonePad();
      c666->Update();
      c666->Print(("BKG-Sub-RESULTS"+ext).c_str());
    }// end of if(sub_big)
    
    if(add_big && add_plots){
      c777 = new TCanvas("c777","BKG-Added Results",0,0,xlen*2,ylen*2);
      c777->Divide(2,2,0.002,0.002,0);
      c777->cd(1);
      c701->DrawClonePad();
      c777->cd(2);
      c702->DrawClonePad();
      c777->cd(3);
      c703->DrawClonePad();
      c777->cd(4);
      c704->DrawClonePad();
      c777->Update();
      c777->Print(("BKG-Add-RESULTS"+ext).c_str());
    }// end of if(add_big)
    
    if(conf_big && conf_plots){
      c888 = new TCanvas("c888","Conf Results",0,0,xlen*2,ylen*2);
      c888->Divide(2,2,0.002,0.002,0);
      c888->cd(1);
      c801->DrawClonePad();
      c888->cd(2);
      c802->DrawClonePad();
      c888->cd(3);
      c803->DrawClonePad();
      c888->cd(4);
      c804->DrawClonePad();
      c888->Update();
      c888->Print(("Conf-RESULTS"+ext).c_str());
    }// end of if(conf_big)
    
  }// end of if(basic_anal)
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  0vBB LIMITS  ==========================================
  //////////////////////////////////////////////////////////////////////////
  TCanvas *c488;
  if(mclimit && nlim[0]){
    cout<<"\n\n ======  STARTING ANALYSIS OF 2b0v LIMITS  ====== \n\n";
    cout<<"\n\n 0vBB HELENE LIMIT"<<endl;
    cout<<"=================================================================="<<endl;
    anaout<<"\n\n 0vBB HELENE LIMIT"<<endl;
    anaout<<"=================================================================="<<endl;
    
    Double_t lx=2.70,ux=4.00;
    T_helene_eff=HelenLimitf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],lim[0],nlim[0],lx,ux);
    T_helene=Nconst/lim[0][0]->eactivity;
    cout<<"\t T_0v Helene    > "<<scientific<<setprecision(2)<<T_helene<<" yr "<<endl;
    cout<<"\t Corr. T_0v     > "<<scientific<<setprecision(2)<<T_helene*(1.0-deadtime)<<" yr "<<endl<<endl;
    anaout<< "T_0v Helene    > "<<scientific<<setprecision(2)<<T_helene<<" yr "<<endl;
    anaout<< "Corr. T_0v     > "<<scientific<<setprecision(2)<<T_helene*(1.0-deadtime)<<" yr "<<endl<<endl;
    
    if(multi_helene){
      for(lx=2.20;lx<ux;lx+=0.10){
	HelenLimitf("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],lim[0],nlim[0],lx,ux);
	T_helene=Nconst/lim[0][0]->eactivity;
	cout<<"\t T_0v Helene    > "<<scientific<<setprecision(2)<<T_helene<<" yr "<<endl;
	cout<<"\t Corr. T_0v     > "<<scientific<<setprecision(2)<<T_helene*(1.0-deadtime)<<" yr "<<endl<<endl;
	anaout<< "T_0v Helene    > "<<scientific<<setprecision(2)<<T_helene<<" yr "<<endl;
	anaout<< "Corr. T_0v     > "<<scientific<<setprecision(2)<<T_helene*(1.0-deadtime)<<" yr "<<endl<<endl;
      }
    }
    
    cout<<"\n\n 0vBB MC TLIMIT"<<endl;
    cout<<"=================================================================="<<endl;
    anaout<<"\n\n 0vBB MC TLIMIT"<<endl;
    anaout<<"=================================================================="<<endl;
    
    T_mclim_eff=MCLimit("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],lim[0],nlim[0],fit_rebin);
    T_mclim=Nconst/lim[0][0]->eactivity;
    cout<<"\t Rebin          = "<<fixed<<setprecision(0)<<fit_rebin<<endl;
    cout<<"\t T_0v MClimit   > "<<scientific<<setprecision(2)<<T_mclim<<" yr "<<endl;
    cout<<"\t Corr. T_0v     > "<<scientific<<setprecision(2)<<T_mclim*(1.0-deadtime)<<" yr "<<endl<<endl;
    anaout<< "Rebin          = "<<fixed<<setprecision(0)<<fit_rebin<<endl;
    anaout<< "T_0v MClimit   > "<<scientific<<setprecision(2)<<T_mclim<<" yr "<<endl;
    anaout<< "Corr. T_0v     > "<<scientific<<setprecision(2)<<T_mclim*(1.0-deadtime)<<" yr "<<endl<<endl;
    
    //limEvents=lim[0][0]->eactivity*data[0]->totaltime*T_mclim_eff*(1.0-deadtime);
    limEvents=lim[0][0]->eactivity*data[0]->totaltime*T_mclim_eff;
    
    // plot the limit
    //////////////////////////////////////////
    c488 = new TCanvas("c488","0v Limit",0,0,xlen,ylen);
    DrawLim("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],lim[0],nlim[0],0.0,4.0,plot_rebin,limEvents,"ee",isotope,decaymode);
    c488->Update();
    c488->Print("Conf-0v-limit.C");
    if(decaymode!="") c488->Print((decaymode+"_conf-0v-limit.C").c_str());
    c488->SetLogy();
    c488->Update();
    c488->Print(("Conf-0v-limit"+ext).c_str());
    if(decaymode!="") c488->Print((decaymode+"_conf-0v-limit"+ext).c_str());
    
    // generate collie histograms
    //////////////////////////////////////////
    if(nsig[0]==1 && nlim[0]==1 && gen_collie){
      cout<<"\nCOLLIE NORMALIZATION = "<<fixed<<setprecision(2)<<limEvents<<endl;
      createCollieHists("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],lim[0],nlim[0],fit_rebin,limEvents,decaymode);
    }
    
    if(multi_tlimit){
      for(Int_t m=1;m<7;m++){
	MCLimit("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],lim[0],nlim[0],m);
	T_mclim=Nconst/lim[0][0]->eactivity;
	cout<<"\t Rebin          = "<<fixed<<setprecision(0)<<fit_rebin<<endl;
	cout<<"\t T_0v MClimit   > "<<scientific<<setprecision(2)<<T_mclim<<" yr "<<endl;
	cout<<"\t Corr. T_0v     > "<<scientific<<setprecision(2)<<T_mclim*(1.0-deadtime)<<" yr "<<endl<<endl;
	anaout<< "Rebin          = "<<fixed<<setprecision(0)<<fit_rebin<<endl;
	anaout<< "T_0v MClimit   > "<<scientific<<setprecision(2)<<T_mclim<<" yr "<<endl;
	anaout<< "Corr. T_0v     > "<<scientific<<setprecision(2)<<T_mclim*(1.0-deadtime)<<" yr "<<endl<<endl;
      } 
    }
    
    // calculate and print out mass results
    //////////////////////////////////////////
    cout<<"\n\n LEPTON NUMBER VIOLATION "<<endl;
    cout<<"=================================================================="<<endl;
    anaout<<"\n\n LEPTON NUMBER VIOLATION "<<endl;
    anaout<<"=================================================================="<<endl;
    
    T_helene_corr     = T_helene*(1.0-deadtime);
    T_mclim_corr      = T_mclim*(1.0-deadtime);
        
    if(decaymode=="m1"){ 
      T_0v=T_mclim_corr;
      
      // Rodin 2007
      RQRPA_ga125 = TMath::Sqrt((0.050*0.050)*(0.98e+27/T_0v));
      RQRPA_ga100 = TMath::Sqrt((0.050*0.050)*(1.12e+27/T_0v));
      
      // Simkovic 2008
      RQRPA_Jast_ga125  = TMath::Sqrt((0.050*0.050)*(7.90e+26/T_0v));
      RQRPA_Jast_ga100  = TMath::Sqrt((0.050*0.050)*(13.9e+26/T_0v));
      RQRPA_UCOM_ga125  = TMath::Sqrt((0.050*0.050)*(4.43e+26/T_0v));
      RQRPA_UCOM_ga100  = TMath::Sqrt((0.050*0.050)*(8.27e+26/T_0v));
      
      // Suhonen 2007
      QRPA_UCOM_ga125 = TMath::Sqrt((1.0*1.0)*(4.70e+23/T_0v));
      QRPA_UCOM_ga100 = TMath::Sqrt((1.0*1.0)*(6.10e+23/T_0v));
      
      cout<<endl<<endl;
      cout<<"Determine mass by ratioing the half-life from publications (eV) "<<endl;
      cout<<"------------------------------------------------------------------"<<endl;
      cout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga125"<<" < "<<fixed<<setprecision(1)<<RQRPA_ga125<<endl;
      cout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga100"<<" < "<<fixed<<setprecision(1)<<RQRPA_ga100<<endl;
      cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga125"<<" < "<<fixed<<setprecision(1)<<RQRPA_Jast_ga125<<endl;
      cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga100"<<" < "<<fixed<<setprecision(1)<<RQRPA_Jast_ga100<<endl;
      cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga125"<<" < "<<fixed<<setprecision(1)<<RQRPA_UCOM_ga125<<endl;
      cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga100"<<" < "<<fixed<<setprecision(1)<<RQRPA_UCOM_ga100<<endl;
      cout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga125"<<" < "<<fixed<<setprecision(1)<<QRPA_UCOM_ga125<<endl;
      cout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga100"<<" < "<<fixed<<setprecision(1)<<QRPA_UCOM_ga100<<endl;
      
      anaout<<"Determine mass by ratioing the half-life from publications "<<endl;
      anaout<<"------------------------------------------------------------------"<<endl;
      anaout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga125"<<" < "<<fixed<<setprecision(1)<<RQRPA_ga125<<endl;
      anaout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga100"<<" < "<<fixed<<setprecision(1)<<RQRPA_ga100<<endl;
      anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga125"<<" < "<<fixed<<setprecision(1)<<RQRPA_Jast_ga125<<endl;
      anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga100"<<" < "<<fixed<<setprecision(1)<<RQRPA_Jast_ga100<<endl;
      anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga125"<<" < "<<fixed<<setprecision(1)<<RQRPA_UCOM_ga125<<endl;
      anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga100"<<" < "<<fixed<<setprecision(1)<<RQRPA_UCOM_ga100<<endl;
      anaout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga125"<<" < "<<fixed<<setprecision(1)<<QRPA_UCOM_ga125<<endl;
      anaout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga100"<<" < "<<fixed<<setprecision(1)<<QRPA_UCOM_ga100<<endl;
      
    }
    
    if(decaymode=="m5"){ 
      T_0v=T_mclim_corr;
      
      // Simkovic 2008
      RQRPA_Jast_ga125  = 1/TMath::Sqrt(T_0v*G_0vM_Sim*TMath::Power(1.01,2));
      RQRPA_UCOM_ga125  = 1/TMath::Sqrt(T_0v*G_0vM_Sim*TMath::Power(1.31,2));
      
      // Suhonen 2007
      QRPA_UCOM_ga125 = 1/TMath::Sqrt(T_0v*G_0vM_Suh*TMath::Power(3.12,2));
      
      cout<<endl<<endl;
      cout<<"Determine majoron-neutrino coupling "<<endl;
      cout<<"------------------------------------------------------------------"<<endl;
      //cout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga125"<<" < "<<scientific<<setprecision(1)<<RQRPA_ga125<<endl;
      //cout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga100"<<" < "<<scientific<<setprecision(1)<<RQRPA_ga100<<endl;
      cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga125"<<" < "<<scientific<<setprecision(1)<<RQRPA_Jast_ga125<<endl;
      //cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga100"<<" < "<<scientific<<setprecision(1)<<RQRPA_Jast_ga100<<endl;
      cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga125"<<" < "<<scientific<<setprecision(1)<<RQRPA_UCOM_ga125<<endl;
      //cout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga100"<<" < "<<scientific<<setprecision(1)<<RQRPA_UCOM_ga100<<endl;
      cout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga125"<<" < "<<scientific<<setprecision(1)<<QRPA_UCOM_ga125<<endl;
      //cout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga100"<<" < "<<scientific<<setprecision(1)<<QRPA_UCOM_ga100<<endl;
      
      anaout<<"Determine majoron-neutrino coupling "<<endl;
      anaout<<"------------------------------------------------------------------"<<endl;
      //anaout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga125"<<" < "<<scientific<<setprecision(1)<<RQRPA_ga125<<endl;
      //anaout<<"  "<<setw(30)<<left<<"Rodin - RQRPA_ga100"<<" < "<<scientific<<setprecision(1)<<RQRPA_ga100<<endl;
      anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga125"<<" < "<<scientific<<setprecision(1)<<RQRPA_Jast_ga125<<endl;
      //anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_Jast_ga100"<<" < "<<scientific<<setprecision(1)<<RQRPA_Jast_ga100<<endl;
      anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga125"<<" < "<<scientific<<setprecision(1)<<RQRPA_UCOM_ga125<<endl;
      //anaout<<"  "<<setw(30)<<left<<"Simkovic - RQRPA_UCOM_ga100"<<" < "<<scientific<<setprecision(1)<<RQRPA_UCOM_ga100<<endl;
      anaout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga125"<<" < "<<scientific<<setprecision(1)<<QRPA_UCOM_ga125<<endl;
      //anaout<<"  "<<setw(30)<<left<<"Suhonen - QRPA_UCOM_ga100"<<" < "<<scientific<<setprecision(1)<<QRPA_UCOM_ga100<<endl;
      
    }
    
    cout<<"\n\n ======  DONE WITH 2b0v LIMITS  ====== \n\n";
  }
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  FINAL RESULTS PLOT  ===================================
  //////////////////////////////////////////////////////////////////////////
  TCanvas *c998;
  if(show_results && nsig[0]){
    c998 = new TCanvas("c998","ana10 Results",xlen,ylen);
    
    Int_t precision=2;
    Int_t base=(Int_t)TMath::Log10(halflife);
    Int_t statpre=1+precision-(base-(Int_t)TMath::Log10(herror));
    Int_t systpre=1+precision-(base-(Int_t)TMath::Log10(hsystp));
    
    ostringstream o1,o2,o3,o4,o5,o10,o11;
    ostringstream o21,o22,o23,o24,o25,o30;
    ostringstream o31,o32,o33,o34,o35;
    
    if(isotope!="") o1<<isotope<<" Mass = "<<fixed<<setprecision(2)<<Smass<<" grams ";
    if(isotope=="") o1<<"Isotope Mass = "<<fixed<<setprecision(2)<<Smass<<" grams ";
    o2<<"Runtime = "<<fixed<<setprecision(1)<<runtime_corr<<" days ";
    o3<<"2#nu#beta#beta Efficiency = "<<fixed<<setprecision(1)<<eff*100<<" % ";
    o4<<"S/B = "<<fixed<<setprecision(2)<<sig2bkg<<endl;
    o10<<"#Tau^{2#nu#beta#beta}_{1/2}  =  "<<fixed<<setprecision(precision)<<halflife_corr/pow(10.,base)<<"  #pm  "
       <<setprecision(statpre)<<herror/pow(10.,base)<<"(stat) ^{+ "
       <<setprecision(systpre)<<hsystp/pow(10.,base)<<"}_{ - "
       <<setprecision(systpre)<<hsystm/pow(10.,base)<<"}(syst) #times 10^{"<<base<<"} yr ";
    o11<<"M^{2#nu#beta#beta} = "<<fixed<<setprecision(4)<<M_2v<<" #pm "<<M_2v_err<<endl;
    
    Double_t yps=0.99;
    Double_t spc=0.06;
    Double_t xps=0.05;
    
    TLatex *t1  = new TLatex(xps,yps-=spc,o1.str().c_str());
    TLatex *t2  = new TLatex(xps,yps-=spc,o2.str().c_str());
    TLatex *t3  = new TLatex(xps,yps-=spc,o3.str().c_str());
    TLatex *t4  = new TLatex(xps,yps-=spc,o4.str().c_str());
    TLatex *t10 = new TLatex(xps,yps-=.08,o10.str().c_str());
    TLatex *t11 = new TLatex(xps,yps-=.08,o11.str().c_str());
    
    t1->Draw();
    t2->Draw();
    t3->Draw();
    t4->Draw();
    t10->Draw();
    t11->Draw();
    
    
    // 0vbb stuff here!!
    /////////////////////////////////////////////////
    Int_t Tbase;
    TLatex *t21,*t22,*t23,*t24,*t25,*t30,*t31,*t32,*t33,*t34,*t35;
    if(mclimit && nlim[0]){
      Tbase=(Int_t)TMath::Log10(T_mclim_corr);
      
      if(decaymode!="") o21<<"decay mode = "<<decaymode<<" ";
      if(decaymode=="") o21<<"decay mode = (not specified)";
      o22<<"0#nu Efficiency = "<<fixed<<setprecision(1)<<T_mclim_eff*100<<" % ";
      o23<<"Confidence Level = 90 % ";
      o30<<"#Tau^{0#nu}_{1/2}  >  "<<fixed<<setprecision(2)<<T_mclim_corr/pow(10.,Tbase)<<"  #times 10^{"<<Tbase<<"} yr ";
      t21 = new TLatex(xps,yps-=.08,o21.str().c_str());
      t22 = new TLatex(xps,yps-=spc,o22.str().c_str());
      t23 = new TLatex(xps,yps-=spc,o23.str().c_str());
      t30 = new TLatex(xps,yps-=spc,o30.str().c_str());
      t21->Draw();
      t22->Draw();
      t23->Draw();
      t30->Draw();
      
      if(decaymode=="m1"){
	o32<<"Rodin - RQRPA  <  "<<fixed<<setprecision(1)<<RQRPA_ga125<<" - "<<RQRPA_ga100<<" eV ";
	o33<<"Simkovic - RQRPA(Jast)  <  "<<fixed<<setprecision(1)<<RQRPA_Jast_ga125<<" - "<<RQRPA_Jast_ga100<<" eV ";
	o34<<"Simkovic - RQRPA(UCOM)  <  "<<fixed<<setprecision(1)<<RQRPA_UCOM_ga125<<" - "<<RQRPA_UCOM_ga100<<" eV ";
	o35<<"Suhonen - QRPA(UCOM)  <  "<<fixed<<setprecision(1)<<QRPA_UCOM_ga125<<" - "<<QRPA_UCOM_ga100<<" eV ";
	t32 = new TLatex(xps,yps-=.08,o32.str().c_str());
	t33 = new TLatex(xps,yps-=spc,o33.str().c_str());
	t34 = new TLatex(xps,yps-=spc,o34.str().c_str());
	t35 = new TLatex(xps,yps-=spc,o35.str().c_str());
	t32->Draw();
	t33->Draw();
	t34->Draw();
	t35->Draw();
      }
      
      if(decaymode=="m5"){
	o32<<"Simkovic - RQRPA(Jast)  <  "<<scientific<<setprecision(1)<<RQRPA_Jast_ga125<<" ";
        o33<<"Simkovic - RQRPA(UCOM)  <  "<<scientific<<setprecision(1)<<RQRPA_UCOM_ga125<<" ";
	o34<<"Suhonen - QRPA(UCOM)  <  "<<scientific<<setprecision(1)<<QRPA_UCOM_ga125<<" ";
	t32 = new TLatex(xps,yps-=.10,o32.str().c_str());
        t33 = new TLatex(xps,yps-=.07,o33.str().c_str());
        t34 = new TLatex(xps,yps-=.07,o34.str().c_str());
	t32->Draw();
        t33->Draw();
	t34->Draw();
      }
      
    }// end of if(mclimit && nlim[0])
    
    c998->Update();
    c998->Print(("ana10-results"+ext).c_str());
    if(decaymode!="") c998->Print((decaymode+"_results"+ext).c_str());
    
  }// end of if(show_results && nsig[0])
  
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  SUMMARY PLOTS  ========================================
  //////////////////////////////////////////////////////////////////////////
  TCanvas *c999;
  if(decaymode!="" && show_results && nsig[0] && nlim[0] && add_plots && basic_anal){
    c999 = new TCanvas("c999","Results Summary",xlen*2,ylen*2);
    c999->Divide(2,2,0.002,0.002,0);
    c999->cd(4);
    c998->DrawClonePad();
    c999->cd(2);
    c488->DrawClonePad();
    c999->cd(1);
    c701->DrawClonePad();
    c999->cd(3);
    c704->DrawClonePad();
    c999->Update();
    c999->Print((decaymode+"_summary"+ext).c_str());
  }    
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  PRINT OUT SYSTEMATICS  ================================
  //////////////////////////////////////////////////////////////////////////
  if(print_syst){
    cout<<"\n\n CURRENT SYSTEMATICS"<<endl;
    cout<<"=================================================================="<<endl; 
    anaout<<"\n\n CURRENT SYSTEMATICS"<<endl;
    anaout<<"=================================================================="<<endl;
    Int_t sysW=18;
    for(Int_t s=0;s<NUM;s++){
      cout<<setw(sysW)<<syst_name[s]<<" (+) "<<fixed<<setprecision(2)<<syst_p[s]*100<<",  (-) "<<syst_m[s]*100<<" % "<<endl;
      anaout<<setw(sysW)<<syst_name[s]<<" (+) "<<fixed<<setprecision(2)<<syst_p[s]*100<<",  (-) "<<syst_m[s]*100<<" % "<<endl;
    }
    cout<<setw(sysW)<<"TOTAL"<<" (+) "<<fixed<<setprecision(2)<<total_syst_p*100<<",  (-) "<<total_syst_m*100<<" % "<<endl;
    anaout<<setw(sysW)<<"TOTAL"<<" (+) "<<fixed<<setprecision(2)<<total_syst_p*100<<",  (-) "<<total_syst_m*100<<" % "<<endl;
  }
  

  //////////////////////////////////////////////////////////////////////////
  ////  ===========  BACKGROUND SYSTEMATICS  ===============================
  //////////////////////////////////////////////////////////////////////////
  Double_t bkg_realhalf,bkg_minhalf,bkg_maxhalf,bkg_minsyst,bkg_maxsyst;
  if(bkg_syst && nsig[0]){
    cout<<"\n\n ======  CHECKING THE SYSTEMATIC ERROR OF BKGS  ====== \n\n";
    
    cout<<fixed<<setprecision(2);
    anaout<<fixed<<setprecision(2);
        
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
    bkg_realhalf=(1.0-deadtime)*Nconst/sig[0][0]->norm;
    
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin,1.0+scalebkgs);
    bkg_maxhalf=(1.0-deadtime)*Nconst/sig[0][0]->norm;
    
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin,1.0-scalebkgs);
    bkg_minhalf=(1.0-deadtime)*Nconst/sig[0][0]->norm;
    
    bkg_maxsyst=TMath::Abs((bkg_realhalf-bkg_maxhalf)/bkg_realhalf*100);
    bkg_minsyst=TMath::Abs((bkg_realhalf-bkg_minhalf)/bkg_realhalf*100);
    
    cout<<"\n\n NEW SYSTEMATICS"<<endl;
    cout<<"=================================================================="<<endl;
    cout<<"Bkgds  ==>  + "<<scalebkgs*100<<"%  ==>  + "<<bkg_maxsyst<<" % (syst) "<<endl;
    cout<<"Bkgds  ==>  - "<<scalebkgs*100<<"%  ==>  - "<<bkg_minsyst<<" % (syst) "<<endl;
    
    anaout<<"\n\n NEW SYSTEMATICS"<<endl;
    anaout<<"=================================================================="<<endl;
    anaout<<"Bkgds  ==>  + "<<scalebkgs*100<<"%  ==>  + "<<bkg_maxsyst<<" % (syst) "<<endl;
    anaout<<"Bkgds  ==>  - "<<scalebkgs*100<<"%  ==>  - "<<bkg_minsyst<<" % (syst) "<<endl;
    
    cout<<"\n\n ======  DONE WITH SYSTEMATICS OF BKGS  ====== \n\n";
  }// end of if(bkg_syst)
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  ENERGY WINDOW SYSTEMATICS  ============================
  //////////////////////////////////////////////////////////////////////////
  Double_t ewin_realhalf,ewin_minhalf,ewin_maxhalf,ewin_minsyst,ewin_maxsyst,newmin,newmax;
  if(ewin_syst && nsig[0]){
    cout<<"\n\n ======  CHECKING THE ENERGY WINDOW SYSTEMATICS  ====== \n\n";
    
    cout<<fixed<<setprecision(2);
    anaout<<fixed<<setprecision(2);
    
    newmin=0.4;
    newmax=1.1;
    
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
    ewin_realhalf=(1.0-deadtime)*Nconst/sig[0][0]->norm;
        
    Double_t temp=0;
    ewin_minhalf=ewin_realhalf;
    ewin_maxhalf=ewin_realhalf;
    for(Double_t i=newmin; i<=newmax; i+=0.1){
      LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],i,fit_max,fit_rebin);
      temp=(1.0-deadtime)*Nconst/sig[0][0]->norm;
      if(temp>ewin_maxhalf) ewin_maxhalf=temp;
      if(temp<ewin_minhalf) ewin_minhalf=temp;
    }
    
    ewin_maxsyst=TMath::Abs((ewin_realhalf-ewin_maxhalf)/ewin_realhalf*100);
    ewin_minsyst=TMath::Abs((ewin_realhalf-ewin_minhalf)/ewin_realhalf*100);
    
    cout<<fixed<<setprecision(2);
    anaout<<fixed<<setprecision(2);
    
    cout<<"\n\n SYSTEMATICS"<<endl;
    cout<<"=================================================================="<<endl;
    cout<<"E-win  ==>  "<<newmax<<"-"<<fit_max<<"  ==>  + "<<ewin_maxsyst<<" % (syst) "<<endl;
    cout<<"E-win  ==>  "<<newmin<<"-"<<fit_max<<"  ==>  - "<<ewin_minsyst<<" % (syst) "<<endl;
    
    if(!bkg_syst) anaout<<"\n\n SYSTEMATICS"<<endl;
    if(!bkg_syst) anaout<<"=================================================================="<<endl;
    anaout<<endl;
    anaout<<"E-win  ==>  "<<newmax<<"-"<<fit_max<<"  ==>  + "<<ewin_maxsyst<<" % (syst) "<<endl;
    anaout<<"E-win  ==>  "<<newmin<<"-"<<fit_max<<"  ==>  - "<<ewin_minsyst<<" % (syst) "<<endl;
    
    cout<<"\n\n ======  DONE WITH ENERGY WINDOW SYSTEMATICS  ====== \n\n";
  }// end of if(ewin_syst)
  
  
  //////////////////////////////////////////////////////////////////////////
  ////  ===========  MULTI FIT COMPARISSON  ================================
  //////////////////////////////////////////////////////////////////////////
  Double_t new_fit_min,new_fit_max,l_syst,l_real,f_syst,f_real;
  Double_t l_activity=0,l_errActiv=0,l_halflife=0,l_herror=0;
  Double_t f_activity=0,f_errActiv=0,f_halflife=0,f_herror=0;
  
  if(multi_LFit && nsig[0]){
    cout<<"\n\n ======  DOING THE MULTI-FIT  ====== \n\n";
    
    new_fit_min=0.4;
    new_fit_max=4.0;
    
    erlog<<setprecision(4);
    erlog<<"\n\n\t--- Likelyhood Fit ---\n\n";
    erlog<<setw(12)<<"Half-Life"<<setw(12)<<"+/- Error"<<setw(6)<<"min"<<setw(6)<<"max"<<setw(7)<<"Rebin"<<setw(7)<<" Syst"<<endl;
    erlog<<setw(12)<<"========="<<setw(12)<<"========="<<setw(6)<<"==="<<setw(6)<<"==="<<setw(7)<<"====="<<setw(7)<<"====="<<endl;
    
    // get the real half-life from fit
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
    l_real=(1.0-deadtime)*(Nconst/sig[0][0]->norm);
    
    // get a model half-life to scale histrograms by
    LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.1,3,1);
    l_activity=sig[0][0]->norm;
    l_errActiv=sig[0][0]->enorm;
    l_halflife=(1.0-deadtime)*(Nconst/l_activity);
    l_herror=(1.0-deadtime)*(Nconst/(l_activity*l_activity))*(l_errActiv);
    
    TH2D *gr1 = new TH2D("gr1","Likelyhood - HalfLife",ewin_max+1,0,ewin_max,rebin_max,1,rebin_max);
    gr1->SetStats(0);
    gr1->SetXTitle("Min. Energy (x0.1 MeV)");
    gr1->SetYTitle("Rebin");
    gr1->SetZTitle("Half-Life");
    gr1->SetTitleOffset(2,"X");
    gr1->SetTitleOffset(2,"Y");
    gr1->SetAxisRange(l_halflife*0.9,l_halflife*1.1,"Z");
    
    TH2D *gr2 = new TH2D("gr2","Likelyhood - Error"   ,ewin_max+1,0,ewin_max,rebin_max,1,rebin_max);
    gr2->SetStats(0);
    gr2->SetXTitle("Min. Energy (x0.1 MeV)");
    gr2->SetYTitle("Rebin");
    gr2->SetZTitle("Error");
    gr2->SetTitleOffset(2,"X");
    gr2->SetTitleOffset(2,"Y");
    gr2->SetTitleOffset(1.3,"Z");
    gr2->SetAxisRange(l_herror*0.6,l_herror*1.3,"Z");
    
    TH2D *gr3;
    TH2D *gr4;
    if(multi_Fit){
      gr3 = new TH2D("gr3","Binned Fit - Halflife",ewin_max+1,0,ewin_max,rebin_max,1,rebin_max);
      gr4 = new TH2D("gr4","Binned Fit - Error"   ,ewin_max+1,0,ewin_max,rebin_max,1,rebin_max);
      
      // get the real half-life from fit
      FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],fit_min,fit_max,fit_rebin);
      f_real=(1.0-deadtime)*(Nconst/sig[0][0]->norm);
      
      // get a model half-life to scale histrograms by
      FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],0.1,3,1);
      f_activity=sig[0][0]->norm;
      f_errActiv=sig[0][0]->enorm;
      f_halflife=(1.0-deadtime)*(Nconst/f_activity);
      f_herror=(1.0-deadtime)*(Nconst/(f_activity*f_activity))*(f_errActiv);
      
      gr3->SetStats(0);
      gr3->SetXTitle("Min. Energy (x0.1 MeV)");
      gr3->SetYTitle("Rebin");
      gr3->SetZTitle("Half-Life");
      gr3->SetTitleOffset(2,"X");
      gr3->SetTitleOffset(2,"Y");
      gr3->SetAxisRange(f_halflife*0.9,f_halflife*1.2,"Z");
      
      gr4->SetStats(0);
      gr4->SetXTitle("Min. Energy (x0.1 MeV)");
      gr4->SetYTitle("Rebin");
      gr4->SetZTitle("Error");
      gr4->SetTitleOffset(2,"X");
      gr4->SetTitleOffset(2,"Y");
      gr4->SetTitleOffset(1.3,"Z");
      gr4->SetAxisRange(f_herror*0.9,f_herror*1.2,"Z");
    }
    
    
    for(Int_t k=1;k<=rebin_max;k++){
      for(Int_t l=0;l<=ewin_max;l++){
	new_fit_min=0.1*l;
	LFitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],new_fit_min,new_fit_max,k);
	l_activity=sig[0][0]->norm;
	l_errActiv=sig[0][0]->enorm;
	l_halflife=(1.0-deadtime)*(Nconst/l_activity);
	l_herror=(1.0-deadtime)*(Nconst/(l_activity*l_activity))*(l_errActiv);
	if(isnan(l_herror)) l_herror=0;
	
	l_syst=TMath::Abs((l_real-l_halflife)/l_real*100);
	
	gr1->SetBinContent(l,k,l_halflife);
	gr2->SetBinContent(l,k,l_herror);
	
	erlog<<setw(12)<<scientific<<setprecision(2)<<l_halflife<<setw(12)<<setprecision(2)<<l_herror<<setw(6)<<fixed<<setprecision(1)<<new_fit_min<<setw(6)<<new_fit_max<<setw(5)<<k<<setw(9)<<fixed<<setprecision(2)<<l_syst<<endl;
	
      }//end of energy range loop
    }//end of rebin loop
    
    if(multi_Fit){
      erlog<<"\n\n\t--- Binned Fit ---\n\n";
      erlog<<setw(12)<<"Half-Life"<<setw(12)<<"+/- Error"<<setw(6)<<"min"<<setw(6)<<"max"<<setw(7)<<"Rebin"<<setw(7)<<" Syst"<<endl;
      erlog<<setw(12)<<"========="<<setw(12)<<"========="<<setw(6)<<"==="<<setw(6)<<"==="<<setw(7)<<"====="<<setw(7)<<"====="<<endl;
      
      for(Int_t k=1;k<=rebin_max;k++){
	for(Int_t l=0;l<=ewin_max;l++){
	  new_fit_min=0.1*l;
	  FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],new_fit_min,new_fit_max,k);
	  f_activity=sig[0][0]->norm;
	  f_errActiv=sig[0][0]->enorm;
	  f_halflife=(1.0-deadtime)*(Nconst/f_activity);
	  f_herror=(1.0-deadtime)*(Nconst/(f_activity*f_activity))*(f_errActiv);
	  if(isnan(f_herror)) f_herror=0;
	  
	  f_syst=TMath::Abs((f_real-f_halflife)/f_real*100);
	  
	  gr3->SetBinContent(l,k,f_halflife);
	  gr4->SetBinContent(l,k,f_herror);
	  
	  erlog<<setw(12)<<scientific<<setprecision(2)<<f_halflife<<setw(12)<<setprecision(2)<<f_herror<<setw(6)<<fixed<<setprecision(1)<<new_fit_min<<setw(6)<<new_fit_max<<setw(5)<<k<<setw(9)<<fixed<<setprecision(2)<<f_syst<<endl;
	  
	}//end of energy range loop
      }//end of rebin loop
    }//end of if(multi_Fit)
    
    TCanvas *c988 = new TCanvas("c988","Multi-Fit Results",xlen*2,ylen*2);
    c988->Divide(2,2,0.002,0.002,0);
    c988->cd(1);
    gr1->Draw("lego2");
    c988->Update();
    c988->cd(2);
    gr2->Draw("lego2");
    c988->Update();
    if(multi_Fit){
      c988->cd(3);
      gr3->Draw("lego2");
      c988->Update();
      c988->cd(4);
      gr4->Draw("lego2");
      c988->Update();
    }
    c988->Update();
    c988->Print(("EE-MULTIFIT"+ext).c_str());
    
    cout<<"\n\n ======  DONE WITH MULTI-FIT  ====== \n\n";
  }//end of if(multi_LFit)  
  
  
  cout<<"\n\n\n\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n";
  cout<<"\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n";
  cout<<"\n\t ======  !!! DONE WITH EVERYTHING !!!  ====== \n\n";
}//end of urun()


