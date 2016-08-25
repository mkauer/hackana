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

Int_t CopyShape1D(TH1 *src, TH1 *tgt){
  //--------------------------------------------
  //  Copy shape of hist. src to tgt, preserves target's 
  //  normalization. error(i)=src(i)* Sqrt(sum(error_tgt(j)^2) )/(sum tgt(j))
  Double_t norm = tgt->Integral();
  Int_t nbin = tgt->GetNbinsX();
  Double_t errort=0;
  Double_t errors=0;
  if(tgt->GetNbinsX() != src->GetNbinsX()){
    std::cout<<"CopyShape1D error:: can't copy shape! src and tgt have different number of bins"<<std::endl;
    return 1;
  }

  for(int i=0;i<nbin;i++){
    errort += pow(tgt->GetBinError(i),2);
    errors += pow(src->GetBinError(i),2);
  }
  errort = sqrt(errort)/norm;
  errors = sqrt(errors)/src->Integral();
  tgt->Add(src,tgt,norm/src->Integral(),0);
  for(int i=0;i<nbin;i++){
    tgt->SetBinError(i,tgt->GetBinContent(i)*errort/errors);
  }
  return 0;
}

int compare_doubles(const void *a,const void * b){
  Float_t temp = *(Float_t *)a - *(Float_t *)b;
  if(temp >0) return 1;
  if(temp<0) return -1;
  return 0;
}
