#include "ana.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include "../../utils/h10.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TMath.h"
#include "TLimit.h"
#include "TLimitDataSource.h"
#include "TConfidenceLevel.h"
#include "TPaveText.h"
#include "TLatex.h"

/**************************************************************
This function is to generate the histograms needed for the 
COLLIE package to set limits on the Majoron emissions.
**************************************************************/

void createCollieHists(const char *hname,ana10 *data,bana10 **bgr,Int_t nbgr,bana10 **sig,Int_t nsig,bana10 **lim,Int_t nlim,Int_t rebin,Double_t limEvents,string decaymode){
  std::cout<<"CREATING HISTOGRAMS FOR COLLIE"<<std::endl;
  if(nsig!=1){
    std::cout<<"\t NEED ONLY 1 SIGNAL MC "<<std::endl;
    return;  
  }
  
  if(nlim!=1){
    std::cout<<"\t NEED ONLY 1 LIMIT MC "<<std::endl;
    return;  
  }
  
  string outfile;
  if(decaymode!=""){
    outfile=string((decaymode+"_collie-hists.root").c_str());
  }else{
    outfile="collie-hists.root";
  }
  std::cout<<"CREATING ROOTFILE "<<outfile<<std::endl;
  TFile* fout = new TFile(outfile.c_str(),"RECREATE");
  
  TH1D *mydat=(TH1D*)data->Geth1(hname)->Clone("data"); 
  TH1D *mybak=(TH1D*)bgr[0]->Geth1(hname)->Clone("bkgs_mc"); 
  TH1D *mysig=(TH1D*)sig[0]->Geth1(hname)->Clone("signal_mc"); 
  TH1D *mytot=(TH1D*)bgr[0]->Geth1(hname)->Clone("total_mc");
  TH1D *mylim=(TH1D*)lim[0]->Geth1(hname)->Clone("limit_mc"); 
  
  mybak->Reset();
  mysig->Reset();
  mytot->Reset();
  
  for(Int_t i=0;i<nbgr;i++){
    mybak->Add(mybak,bgr[i]->Geth1(hname),1.,bgr[i]->norm);
    mytot->Add(mytot,bgr[i]->Geth1(hname),1.,bgr[i]->norm);
  }
  mysig->Add(mysig,sig[0]->Geth1(hname),1.,sig[0]->norm);
  mytot->Add(mytot,sig[0]->Geth1(hname),1.,sig[0]->norm);
  
  if(limEvents!=-1) mylim->Scale(limEvents/mylim->Integral());
  
  mydat->Rebin(rebin);
  mybak->Rebin(rebin);
  mysig->Rebin(rebin);
  mytot->Rebin(rebin);
  mylim->Rebin(rebin);
  
  mydat->Write();
  mybak->Write();
  mysig->Write();
  mytot->Write();
  mylim->Write();
  
  fout->Close();
  std::cout<<"FINISHED CREATING HISTOGRAMS FOR COLLIE"<<std::endl;
}

