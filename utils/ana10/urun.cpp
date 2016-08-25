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

Int_t Urun(bana10 *data[],bana10 *bgr[][MAXCOMP],bana10 * sig[][MAXCOMP],bana10 * lim[][MAXCOMP],Int_t nbgr[],Int_t nsig[],Int_t nlim[],Int_t nchnls, string namechnl[]){

// Fit signal value
  if(nsig[0]){
   FitSignal("etot",data[0],bgr[0],nbgr[0],sig[0],nsig[0],1,99);
  }

// Print the details
  cout << "analysis results"<<endl;
  if(data[0]) data[0]->Print();
  for(Int_t i=0;i<nbgr[0];i++){
    bgr[0][i]->Print();
  }
  cout << "Results of Fit"<<endl;
  for(Int_t i =0 ;i<nsig[0];i++){
    sig[0][i]->Print();
  }


  TCanvas *c[MAXCHNLS];
  cout<<"Total channels "<<nchnls<<endl;
  for(Int_t k=0;k<nchnls;k++){
    c[k] = new TCanvas((string("c")+namechnl[k]).c_str(), "Analysis results", 600, 600); 
    if(data[k] && nbgr[k]){
      Drawall_d1("etot",data[k],bgr[k],nbgr[k],sig[k],nsig[k]);
    }
    c[k]->Update();
  }

}
