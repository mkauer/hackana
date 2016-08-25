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

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>

////////////////////////////////////////////////////////////
//  Default uloop.cpp used for COMBINED-EG analysis.
//  
//  version: 09.03.03
//  
//  by Matt Kauer
////////////////////////////////////////////////////////////

void ana10::Loop(int eventnumber){
  
  Int_t i=0,j=0,k=0;
  
  //---INITIALISATIONS---
  ////////////////////////////////////////////////////////
  Int_t results[200]={0};  // 1-20 preselect, 21-99 slim, 0,100+ ana10
  
  Int_t runnum[1000000];
  Int_t evennum[1000000];
  
  Int_t nbytes=0;
  Int_t counter=0;
  Int_t nentries=Int_t(fChain->GetEntries());
  Int_t nb=-1;
  Int_t treenum=-1;
  
  //---ENERGY AND ANGLE CUTS----
  ////////////////////////////////////////////////////////
  Float_t const Ethr       = 0.400; // total energy minimum (0.4 default)
  Float_t const EeThr      = 0.200; // electron energy minimum (0.2 default)
  Float_t const EgThr      = 0.200; // gamma energy minimum (0.2 default)
  Float_t const Cos_max    = 0.9;   // cosine maximum (1.0 default)
  Float_t const minLength  = 50;    // minimum length of track
  
  //---INR DIMENSIONS----
  ////////////////////////////////////////////////////////
  Float_t const inr_sect_max  = 5.99;
  Float_t const inr_sect_min  = 5.87;
  Float_t const inr_z_max     = 68.0;
  Float_t const inr_z_min     = 0.0;
   
  //---ITEP DIMENSIONS----
  ////////////////////////////////////////////////////////
  Float_t const itep_sect_max = 5.99;
  Float_t const itep_sect_min = 5.87;
  Float_t const itep_z_max    = 114.0;
  Float_t const itep_z_min    = 68.0;
  
  //---SET RUN PERIOD---
  ////////////////////////////////////////////////////////
  SetPeriod(ana::P1,true);
  SetPeriod(ana::P2a,true);
  SetPeriod(ana::P2b,true);
  SetPeriod(ana::P3,true);
  
  //---SET RUN STATUS---
  ////////////////////////////////////////////////////////
  SetRuns(ana::STANDARD_RUN,true); //set others to false if you use this!
  SetRuns(ana::GOOD_RUN,false);
  SetRuns(ana::OK_RUN,false);
  SetRuns(ana::NOAIR_RUN,false);
  SetRuns(ana::AFTER_EC_RUN,false);
  SetRuns(ana::HV_RUN,false);
  SetRuns(ana::ALL_BUT_BAD_RUN,false);
  
  
  //-----CREATE YOUR HISTOGRAMS!!!  ana uses --30 keV per bin--
  ////////////////////////////////////////////////////////
  Addh1("etot", "Total Energy (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30keV per bin
  Geth1("etot")->SetXTitle("MeV");
  Addh1("cosEG", "eg Cosine", 100, -1.0, 1.0);  // 0.025 per bin
  Geth1("cosEG")->SetXTitle("Cosine");
  Addh1("electronE", "Electron Energy (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30 keV per bin
  Geth1("electronE")->SetXTitle("MeV");
  Addh1("gammaE", "Gamma Energy (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30 keV per bin
  Geth1("gammaE")->SetXTitle("MeV");
  Addh1("probInt", "Internal Prob of Event", 100, 0.0, 1.0);
  Addh1("probExt", "External Prob of Event", 100, 0.0, 1.0);
  Addh1("trackL","Straight Length of Track (cm)",100,0,300);
  Geth1("trackL")->SetXTitle("cm");
  Addh1("gammaL","Straight Length of Gamma (cm)",100,0,300);
  Geth1("gammaL")->SetXTitle("cm");
  Addh1("curveyL","Curvey Length of Track (cm)",100,0,300);
  Geth1("curveyL")->SetXTitle("cm");
  
  
  // bool true or false if dealing with data or not
  Bool_t isdata=(name=="data" || name=="DATA" || name=="Phase1" || name=="Phase2");
  //Bool_t isdata=(run>0);
  
  
  if (fChain == 0) return;
  std::cout<<"Number of entries in the chain "<<nentries<<std::endl;
  
  // Main selection LOOP, fills booked hystos
  for (Int_t jentry=eventnumber; jentry<nentries; jentry++) {
    
    Int_t resnum=100;  //results[i] --> pick what [i] to start at for ana cuts
    
    Bool_t preselect = false;
    Int_t ientry = LoadTree(jentry);
    Int_t nb = fChain->GetEntry(jentry);
    
    if(ientry<0) break;
    results[0]++;
    
    if(jentry%10000==0) std::cout<<jentry<<" events passed "<<std::endl;
    
    //---rerun slim2 cuts to make sure all cuts are in place!
    // Check the makefile to see which cuts are being called!
    if(!Cut(results)) continue;

    //---check that run is in the runlist 
    // save selected runs duration for normalization of MC later.	
    Int_t iitime;
    Bool_t goodrun=RunList(abs(run),iitime);
    if(!goodrun) continue;
    AddRunTime(abs(run),iitime);
    results[resnum]++;
    resnum++;
    /*
    //---check PMT status. Ignore events w/o LTC correction
    Bool_t goodpms=true;
    for(k=0;k<Nsc;k++){
      Int_t pmstat=pmstatus->get(Int_t(run),Int_t(Sc[k][1]),Int_t(Sc[k][2]),Int_t(Sc[k][3]),Int_t(Sc[k][4]));
      goodpms = goodpms && GoodPMStatus(pmstat);
    }
    if(!goodpms) continue;
    results[resnum]++;
    resnum++;
    */
    
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ///////  EXTRA CUTS! WOOOOOO!!!  ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    //---1 track, 2 scints, negative charge
    if(Nbr_tks!=1) continue;
    if(Nsc!=2) continue;
    if(Qc[0]>=0) continue;
    results[resnum]++;
    resnum++;
    
    //---setup correct scint indexes for the track and gamma!!!
    Int_t ISCtrack; // scint number associated with the track
    Int_t ISCgamma; // scint number associated with the gamma
    ISCtrack=Ind_scintil[0]-1;
    if(ISCtrack==0) ISCgamma=1;
    if(ISCtrack==1) ISCgamma=0;
    
    //---check for track associated with scint
    if(ISCtrack!=0 && ISCtrack!=1) continue;
    results[resnum]++;
    resnum++;
    
    //---check energy thresholds
    Float_t Eelec=Sc[ISCtrack][8]*1000;
    Float_t Egamma=Sc[ISCgamma][8]*1000;
    if(Eelec<EeThr) continue;
    if(Egamma<EgThr) continue;
    Float_t Etot=Eelec+Egamma;
    if(Etot<Ethr) continue;
    results[resnum]++;
    resnum++;
    
    //---checks that the electron hit the front of the scint
    if(Myimpact[0]==66) continue;
    results[resnum]++;
    resnum++;
    
    //---check that track strarts in first 2 layers
    Bool_t fclose=false;
    for(j=0;j<2;j++){
      fclose=fclose||checklayer(0,j);
    }
    if(!fclose) continue;
    results[resnum]++;
    resnum++;
    
    
    //---use more strict requirements for foil dimensions
    Bool_t inr_sector=false;
    Bool_t itep_sector=false;
    Float_t xsecv=Xvert/TMath::Sqrt(Xvert*Xvert+Yvert*Yvert);
    xsecv=TMath::ACos(xsecv);
    if(Yvert<0) xsecv=2.*TMath::Pi()-xsecv;
    xsecv=20.*(xsecv/2./TMath::Pi());
    inr_sector=(xsecv>=inr_sect_min && xsecv<=inr_sect_max);
    inr_sector=(inr_sector && Zvert<=inr_z_max && Zvert>=inr_z_min);
    itep_sector=(xsecv>=itep_sect_min && xsecv<=itep_sect_max);
    itep_sector=(itep_sector && Zvert<=itep_z_max && Zvert>=itep_z_min);
    if(!(inr_sector || itep_sector)) continue;
    results[resnum]++;
    resnum++;
    
    //---check that if both tracks on same side, then NO fast gg hits <15cm
    //---check that <=2 gg hits, not associated with tracks are >15cm away
    if(NAfasthits()!=0) continue;
    results[resnum]++;
    resnum++;
    
    //---check that total track length is >50 cm
    Float_t dlen;
    Float_t Ltrack=TMath::Sqrt(pow(Sc[ISCtrack][9]-Xvert,2)+pow(Sc[ISCtrack][10]-Yvert,2)+pow(Sc[ISCtrack][11]-Zvert,2));
    Float_t Lgamma=TMath::Sqrt(pow(Sc[ISCgamma][9]-Xvert,2)+pow(Sc[ISCgamma][10]-Yvert,2)+pow(Sc[ISCgamma][11]-Zvert,2));
    Float_t Lcurvey=getcurveylength(0,ISCtrack,dlen);
    if(Ltrack<minLength) continue;
    if(Lgamma<minLength) continue;
    if(Lcurvey<minLength) continue;
    results[resnum]++;
    resnum++;
    
    //---find cosine between electron and gamma
    Float_t dif[3];
    dif[0]=-X_foil[0]+Sc[ISCgamma][9];
    dif[1]=-Y_foil[0]+Sc[ISCgamma][10];
    dif[2]=-Z_foil[0]+Sc[ISCgamma][11];
    Float_t rdif=TMath::Sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
    Float_t prod=Cos_dir[0][0]*dif[0]+Cos_dir[0][1]*dif[1]+Cos_dir[0][2]*dif[2];
    Float_t cosineEG=prod/rdif;
    if(cosineEG>Cos_max) continue;
    results[resnum]++;
    resnum++;
    
    //---check for radon events
    checkdelayed();
    Float_t atime1;
    Float_t atime2;
    Int_t alpha1=checksingledhit(atime1);
    Int_t alpha2=checkgroupdhit(atime2);
    Bool_t alphas=(alpha1>0 || alpha2>0);
    if(alphas) continue;
    results[resnum]++;
    resnum++;
    
    //---probability of gammas being internal or external
    Float_t tfoil,dfoil;
    Float_t intProb=TOFbg_int(ISCgamma,tfoil,dfoil);
    Float_t extProb=TOFbg_ext(ISCgamma,tfoil,dfoil);
    if(intProb<0.04 || extProb>0.01) continue;
    results[resnum]++;
    resnum++;
    
    //---check that scint is alone
    Bool_t alone=amIalone(ISCtrack,ISCgamma,1);
    if(!alone) continue;
    results[resnum]++;
    resnum++;
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////
    //  DONE WITH CUTS!!  /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
    
    //---attach weight to each MC event. 
    Float_t w;
    if(run<0){
      w=GetEventWeight(); // MC
    }
    else{
      w=1.; // REAL
    }
    
    
    //--FILL YOUR HISTROGRAMS WITH DATA!!!
    Geth1("etot")->Fill(Etot,w);
    Geth1("electronE")->Fill(Eelec,w);
    Geth1("cosEG")->Fill(cosineEG,w);
    Geth1("gammaE")->Fill(Egamma,w);
    Geth1("probInt")->Fill(intProb,w); 
    Geth1("probExt")->Fill(extProb,w); 
    Geth1("trackL")->Fill(Ltrack,w);
    Geth1("gammaL")->Fill(Lgamma,w);
    Geth1("curveyL")->Fill(Lcurvey,w);
    
    
    //--SAVE THE RUN AND EVENT NUMBER FOR REAL DATA
    if(isdata){
      std::cout<<std::endl<<"@@@@@@@@@@  Runnumber = "<<run<<" and "<<"Event number = "<<Myievent<<"  @@@@@@@@@@    "<<counter<<std::endl;
      runnum[counter]=run;
      evennum[counter]=Myievent;
    }
    
    counter++;
    
  }
  
  //--WRITE OUT SOME DETAILS TO FILE
  std::cout<<std::endl;
  for(i=0;i<200;i++){ 
    if(results[i]) std::cout<<" Events kept after cut "<<i<<" = "<<results[i]<<std::endl;
  } 
  std::cout<<"\n\t Number of kept events = "<<counter<<"\n";
  
  anaout<<"Kept events for -->  "<<setw(20)<<left<<name<<" = "<<counter<<endl;
  
}

