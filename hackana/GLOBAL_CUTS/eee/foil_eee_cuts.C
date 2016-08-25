#include "h10.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>


Int_t h10::Cut(Int_t * results){
  
  //---FOIL DIMENSIONS
  Float_t const sect_min = 5.8500;
  Float_t const sect_max = 6.0100;
  Float_t const z_min =   -120.0;
  Float_t const z_max =    120.0;
  /////////////////////////////////////////////////////////
 
  Int_t resnum=20;
  Bool_t keep_event=true;
  Float_t sector_pos;
  
  //--- 2 tracks and 2 scints or 3 tracks and 3 scints
  //if(!((Nbr_tks==2 && Nsc==2)||(Nbr_tks==3 && Nsc==3))) return 0;
  if(!(Nbr_tks==3 && Nsc==3)) return 0;
  for(Int_t i=0;i<Nsc;i++){
    if(Ind_scintil[i]==0) keep_event=false;
    if(Qc[i]>=0) keep_event=false;
  }
  if(!keep_event) return 0;
  results[resnum]++;
  resnum++;
  
  //--- if the recontruction fails, the track gets set to 0,0
  for(Int_t i=0;i<Nsc;i++){
    if(fabs(X_foil[i])<0.1 && fabs(Y_foil[i])<0.1) keep_event=false;
  }
  if(!keep_event) return 0;
  results[resnum]++;
  resnum++;
  
  //--- check the sector position of all tracks
  for(Int_t i=0;i<Nsc;i++){
    sector_pos=X_foil[i]/TMath::Sqrt(X_foil[i]*X_foil[i] + Y_foil[i]*Y_foil[i]);
    sector_pos=TMath::ACos(sector_pos);
    if(Y_foil[i]<0) sector_pos=2.*TMath::Pi()-sector_pos;
    sector_pos=20.*(sector_pos/2./TMath::Pi());
    if(sector_pos<sect_min || sector_pos>sect_max) keep_event=false;
    if(Z_foil[i]<z_min || Z_foil[i]>z_max) keep_event=false;
  }
  if(!keep_event) return 0;
  results[resnum]++;
  resnum++;
  
  
  return 1;
}

