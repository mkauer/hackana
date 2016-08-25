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
  Float_t const sect_min=5.8500;
  Float_t const sect_max=6.0100;
  Float_t const z_min=-120.0;
  Float_t const z_max=120.0;
  /////////////////////////////////////////////////////////
 
  Int_t resnum=20;
  
  //---only 2 tracks and 2 scint hits
  if(!(Nbr_tks==2 && Nsc==2)) return 0;
  if(Ind_scintil[0]==0 || Ind_scintil[1]==0) return 0;
  if(Qc[0]>=0 || Qc[1]>=0) return 0;
  results[resnum]++;
  resnum++;
  
  //---if the recontruction fails, the track gets set to 0
  if(fabs(X_foil[0])<0.1 && fabs(Y_foil[0])<0.1) return 0;
  if(fabs(X_foil[1])<0.1 && fabs(Y_foil[1])<0.1) return 0;
  results[resnum]++;
  resnum++;
  
  //---check the sector position of track[0]
  Float_t vertx=(X_foil[0]);
  Float_t verty=(Y_foil[0]);
  Float_t vertz=(Z_foil[0]);
  Bool_t insector=false;
  Float_t xsecv=vertx/TMath::Sqrt(vertx*vertx + verty*verty);
  xsecv=TMath::ACos(xsecv);
  if(verty<0) xsecv=2.*TMath::Pi()-xsecv;
  xsecv=20.*(xsecv/2./TMath::Pi());
  insector=(xsecv>=sect_min && xsecv<=sect_max);
  if(!insector) return 0;
  results[resnum]++;
  resnum++;
  
  //---check Z position of track[0]
  insector=(vertz<=z_max && vertz>=z_min);
  if(!insector) return 0;
  results[resnum]++;
  resnum++;
  
  //---check the sector position of track[1]
  vertx=(X_foil[1]);
  verty=(Y_foil[1]);
  vertz=(Z_foil[1]);
  insector=false;
  xsecv=vertx/TMath::Sqrt(vertx*vertx + verty*verty);
  xsecv=TMath::ACos(xsecv);
  if(verty<0) xsecv=2.*TMath::Pi()-xsecv;
  xsecv=20.*(xsecv/2./TMath::Pi());
  insector=(xsecv>=sect_min && xsecv<=sect_max);
  if(!insector) return 0;
  results[resnum]++;
  resnum++;
  
  //---check Z position of track[1]
  insector=(vertz<=z_max && vertz>=z_min);
  if(!insector) return 0;
  results[resnum]++;
  resnum++;
  
  
  
  return 1;
}

