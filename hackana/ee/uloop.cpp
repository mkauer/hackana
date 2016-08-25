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

/*///////////////////////////////////////////////////////////
Default uloop.cpp used for all EE analysis.

version: 09.05.14
 
change log:
----------------------------------------------------------
05.14 ~ add weight to tof plot
05.14 ~ EvE plot to go to 4.2 in both axis 

by Matt Kauer
///////////////////////////////////////////////////////////*/

void ana10::Loop(int eventnumber){
  
  Int_t i=0,j=0,k=0;
  
  //---INITIALISATIONS---
  ////////////////////////////////////////////////////////
  Int_t results[200]={0};  // 1-20 preselect, 21-99 slim, 0,100+ ana10
  
  Int_t runnum[10000];
  Int_t evennum[10000];
  
  Int_t nbytes = 0;
  Int_t counter = 0;
  Int_t nentries = Int_t(fChain->GetEntries());
  Int_t nb = -1;
  Int_t treenum=-1;
  
  
  //---ENERGY AND ANGLE CUTS----
  ////////////////////////////////////////////////////////
  const Float_t Etot_min=0.400; // minimum total energy
  const Float_t Emin=0.200;     // minimum single energy
  const Float_t Cos_max=1.0;    // moller scattering cut
  const Float_t minLength=30;   // minimum length of track
  
  //---INR DIMENSIONS----
  ////////////////////////////////////////////////////////
  const Float_t inr_sect_max=5.99;
  const Float_t inr_sect_min=5.87;
  const Float_t inr_z_max=68.0;
  const Float_t inr_z_min=0.0;
   
  //---ITEP DIMENSIONS----
  ////////////////////////////////////////////////////////
  const Float_t itep_sect_max=5.99;
  const Float_t itep_sect_min=5.87;
  const Float_t itep_z_max=114.0;
  const Float_t itep_z_min=68.0;
  
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
  Addh1("etot", "2e Total Energy (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30keV per bin
  Geth1("etot")->SetXTitle("MeV");
  Addh1("cosa", "2e Cosine", 100, -1.0, 1.0);  // 0.025 per bin
  Geth1("cosa")->SetXTitle("Cosine");
  Addh1("singleE", "Single Electron (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30 keV per bin
  Geth1("singleE")->SetXTitle("MeV");
  Addh1("minE", "Min Energy Electron (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30 keV per bin
  Geth1("minE")->SetXTitle("MeV");
  Addh1("maxE", "Max Energy Electron (MeV)", 140, 0.0, 4.2);  // 0.03 per bin = 30 keV per bin
  Geth1("maxE")->SetXTitle("MeV");
  Addh1("sect", "Sector Position", 60, 5.8, 6.1);  // 0.005 per bin
  Addh2("foil", "Event Position on Foil", 200, 5.85, 6.01, 200, -5, 119);
  Addh2("EvE",  "Emin vs Emax", 140, 0.0, 4.2, 140, 0.0, 4.2);
  Addh1("tof1", "TOF Before Cuts", 100, 0.0, 1.0);
  Addh1("tof2", "TOF After Cuts", 100, 0.0, 1.0);
  
  
  // bool true or false if dealing with data or not
  Bool_t isdata=(name=="data" || name=="DATA" || name=="Phase1" || name=="Phase2");
  //Bool_t isdata=(run>0);
  Int_t wi=9;  //width of erlog setw() event number
  //erlog<<run<<" "<<Myievent<<endl;
  
  
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
    for(Int_t k=0;k<Nsc;k++){
      Int_t pmstat=pmstatus->get(Int_t(run),Int_t(Sc[k][1]),Int_t(Sc[k][2]),Int_t(Sc[k][3]),Int_t(Sc[k][4]));
      goodpms = goodpms && GoodPMStatus(pmstat);
    }
    if(!goodpms) continue;
    results[resnum]++;
    resnum++;
    */

    //---time of flight pre-check
    Float_t t1, t2;
    Float_t pin1=TOF2b_int(t1,t2);
    Geth1("tof1")->Fill(pin1);
    
    
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ///////  EXTRA CUTS! WOOOOOO!!!  ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    
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
    if(!(inr_sector || itep_sector)){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Event vertex out of sector range!  Sector = "<<xsecv<<"  Zvert = "<<Zvert<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;
    
    //---check that no other scints had energy deposited
    if(Nsc>2){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Too many scint hits = "<<Nsc<<endl;
      continue;
    }
    results[resnum]++;
    resnum++;

    //---check that if both tracks on same side, then NO fast gg hits <15cm
    //---check that <=2 gg hits, not associated with tracks are >15cm away
    if(NAfasthits()!=0){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Too many fast GG hits = "<<NAfasthits()<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;

    //---check that both tracks >30 cm
    Float_t dlen;
    Bool_t trk_lngth=(gethelixl(0,dlen)>minLength && gethelixl(1,dlen)>minLength);
    if(!trk_lngth){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Track length too small = "<<gethelixl(0,dlen)+gethelixl(1,dlen)<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;

    //---has a gg hit in one of the 2 closest layers for both tracks
    Bool_t fclose=false;
    for(j=0;j<2;j++){
      fclose=fclose||checklayer(0,j);  //check track 1
    }
    if(!fclose){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"No GG hit in first two layers for Track 1"<<endl;
      continue;
    }
    fclose=false;
    for(j=0;j<2;j++){
      fclose=fclose||checklayer(1,j);  //cheack track 2
    }
    if(!fclose){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"No GG hit in first two layers for Track 2"<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;

    //---distance between track vertexes is less than 2cm in XY and 4cm in Z
    Float_t dvertxy, dvertz;
    Bool_t xyvert;
    dvertxy=sqrt(pow(X_foil[0]-X_foil[1],2)+pow(Y_foil[0]-Y_foil[1],2));
    dvertz=fabs(Z_foil[0]-Z_foil[1]);
    xyvert=(dvertxy<2.0 && dvertz<4.0);
    if(!xyvert){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"dXY = "<<dvertxy<<" dZ = "<<dvertz<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;

    //---checks TOF probability for internal decay
    Float_t tfoil, dtfoil;
    Float_t Pint=TOF2b_int(tfoil,dtfoil);
    Float_t Pext=TOF2b_ext();
    if(Pint>1.000 || Pint<0.040){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Prob internal = "<<Pint<<endl;
      continue;
    }    
    if(Pext<0.000 || Pext>0.010){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Prob external = "<<Pext<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;

    //---checks that the electron hit the front of the scint
    if(Myimpact[0]==66 || Myimpact[1]==66){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Electron doesn't hit the face of scintilator"<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;

    //---check for radon events
    //Int_t alparray[1000];
    //Int_t alp = getalphas(alparray);
    checkdelayed();
    Float_t atime1,atime2;
    Int_t alpha1=checksingledhit(atime1);
    Int_t alpha2=checkgroupdhit(atime2);
    Bool_t alphas=(alpha1>0 || alpha2>0);
    if(alphas){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Tagged as an Alpha hit!"<<endl;
      continue; 
    }    
    results[resnum]++;
    resnum++;

    //---checks for moller scattering
    Float_t numerator=Cos_dir[0][0]*Cos_dir[1][0]+Cos_dir[0][1]*Cos_dir[1][1]+Cos_dir[0][2]*Cos_dir[1][2];
    Float_t denominator=sqrt(pow(Cos_dir[0][0],2)+pow(Cos_dir[0][1],2)+pow(Cos_dir[0][2],2))*sqrt(pow(Cos_dir[1][0],2)+pow(Cos_dir[1][1],2)+pow(Cos_dir[1][2],2));
    Float_t cosalpha=numerator/denominator;
    if(cosalpha>Cos_max){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Cosine too large = "<<cosalpha<<endl;
      continue;
    }    
    results[resnum]++;
    resnum++;
    
    
    //---total E is greater than Etot_min
    Int_t isc1=Ind_scintil[0]-1;
    Int_t isc2=Ind_scintil[1]-1;
    Float_t E1=Sc[isc1][8]*1000;
    Float_t E2=Sc[isc2][8]*1000;
    Float_t Etot=E1+E2;
    /*
      for(Int_t i=0;i<Nsc;i++){
      Etot+=Sc[i][8]*1000;
      }
    */
    if(Etot<Etot_min || E1<Emin || E2<Emin){
      if(isdata) erlog<<run<<" "<<setw(wi)<<left<<Myievent<<"Electron E too small = "<<E1<<" and "<<E2<<endl; 
      continue;
    }    
    results[resnum]++;
    resnum++;
    
    
    //---FIND MINIMUM ENERGY ELECTRON!!
    Float_t minenergy,maxenergy;    
    if(Sc[isc1][8]<=Sc[isc2][8]){
      minenergy = Sc[isc1][8]*1000;
      maxenergy = Sc[isc2][8]*1000;
    }
    else{
      minenergy = Sc[isc2][8]*1000;
      maxenergy = Sc[isc1][8]*1000;
    }
    
    
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ///////  DONE WITH CUTS!  //////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    
    //---attach weight to each MC event. 
    Float_t w;
    if(run<0){
      w=GetEventWeight(); // MC
    }
    else{
      w=1.; // REAL
    }
    
    
    //---FILL YOUR HISTROGRAMS WITH DATA!!!
    Geth1("etot")->Fill(Etot,w);
    Geth1("minE")->Fill(minenergy,w);
    Geth1("cosa")->Fill(cosalpha,w);
    Geth1("singleE")->Fill(minenergy,w);
    Geth1("singleE")->Fill(maxenergy,w);
    Geth1("maxE")->Fill(maxenergy,w);
    Geth1("sect")->Fill(xsecv); 
    Geth2("foil")->Fill(xsecv,Zvert); 
    Geth2("EvE")->Fill(maxenergy,minenergy);
    Geth1("tof2")->Fill(Pint,w);
    
    
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
  for(Int_t i=0;i<200;i++){ 
    if(results[i]) std::cout<<" Events kept after cut "<<i<<" = "<<results[i]<<std::endl;
  } 
  std::cout<<std::endl<<" Number of kept events = "<<counter<<std::endl;
  std::cout<<std::endl;
  
  if(isdata){
    anaout<<"Event Selection for "<<name<<endl;
    for(Int_t i=0;i<counter;i++){
      anaout<<runnum[i]<<" "<<evennum[i]<<endl;
    }
    anaout<<endl<<endl;
  }
  anaout<<"Kept events for -->  "<<setw(20)<<left<<name<<" = "<<counter<<endl;
  
}

