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
#include <iomanip>

////////////////////////////////////////////////////////////
// Default uloop.cpp used for all 1E analysis.
//  
// version: 09.03.03
//  
// + <iomanip> in 2016-08-24
//
// by Matt Kauer
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
  Float_t const Ethr=0.400;   // total energy minimum (0.4 default)
  Float_t const minLength=50; // minimum length of track
  
  //---INR DIMENSIONS----
  ////////////////////////////////////////////////////////
  Float_t const inr_sect_max=5.99;
  Float_t const inr_sect_min=5.87;
  Float_t const inr_z_max=68.0;
  Float_t const inr_z_min=0.0;
   
  //---ITEP DIMENSIONS----
  ////////////////////////////////////////////////////////
  Float_t const itep_sect_max=5.99;
  Float_t const itep_sect_min=5.87;
  Float_t const itep_z_max=114.0;
  Float_t const itep_z_min=68.0;
  
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
  Addh1("etot", "Single Electron Energy (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30keV per bin
  Geth1("etot")->SetXTitle("MeV");
  Addh1("trackL","Straight Length of Track (cm)",100,0,300);
  Geth1("trackL")->SetXTitle("cm");
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
    
    ///////////////////////////////////////////////////////////////////////////
    //  EXTRA CUTS! WOOOOOO!!!  ///////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
    //---1 track, 1 scint, negative charge
    if(Nbr_tks!=1) continue;
    if(Nsc!=1) continue;
    if(Qc[0]>=0) continue;
    results[resnum]++;
    resnum++;
    
    //---setup correct scint indexes for the track and gamma!!!
    Int_t ISCtrack; // scint number associated with the track
    //Int_t ISCgamma; // scint number associated with the gamma
    ISCtrack=Ind_scintil[0]-1;
    //if(ISCtrack==0) ISCgamma=1;
    //if(ISCtrack==1) ISCgamma=0;
    
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
    
    //---tracks have associated scint deposit and E>threshold
    if(Ind_scintil[0]==0) continue;
    Float_t Etot=Sc[ISCtrack][8]*1000;
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
    
    //---check that if both tracks on same side, then NO fast gg hits <15cm
    //---check that <=2 gg hits, not associated with tracks are >15cm away
    if(NAfasthits()!=0) continue;
    results[resnum]++;
    resnum++;
    
    //---check that total track length is >50 cm
    Float_t dlen;
    Float_t Ltrack=TMath::Sqrt(pow(Sc[ISCtrack][9]-Xvert,2)+pow(Sc[ISCtrack][10]-Yvert,2)+pow(Sc[ISCtrack][11]-Zvert,2));
    Float_t Lcurvey=getcurveylength(0,ISCtrack,dlen);
    if(Ltrack<minLength) continue;
    if(Lcurvey<minLength) continue;
    results[resnum]++;
    resnum++;
    
    
    ///////////////////////////////////////////////////////////////////////////
    //  DONE WITH CUTS!!  /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    
    
    //---ATTACH WEIGHT OF EVENT 
    Float_t w;
    if(run<0){
      w=GetEventWeight(); // MC
    }
    else{
      w=1.; // REAL
    }
    
    
    //---FILL YOUR HISTROGRAMS
    Geth1("etot")->Fill(Etot,w);
    Geth1("trackL")->Fill(Ltrack,w);
    Geth1("curveyL")->Fill(Lcurvey,w);
    
    
    //---SAVE THE RUN AND EVENT NUMBER FOR REAL DATA
    if(isdata){
      std::cout<<std::endl<<"@@@@@@@@@@  Runnumber = "<<run<<" and "<<"Event number = "<<Myievent<<"  @@@@@@@@@@    "<<counter<<std::endl;
      runnum[counter]=run;
      evennum[counter]=Myievent;
    }
  
    counter++;
    
  }
  
  
  //---WRITE OUT SOME DETAILS TO FILE
  std::cout<<std::endl;
  for(i=0;i<200;i++){ 
    if(results[i]) std::cout<<" Events kept after cut "<<i<<" = "<<results[i]<<std::endl;
  } 
  std::cout<<"\n\t Number of kept events = "<<counter<<"\n";
  anaout<<"Kept events for -->  "<<setw(20)<<left<<name<<" = "<<counter<<endl;
  
}

