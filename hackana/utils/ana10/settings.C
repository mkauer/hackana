#include "ana.hpp"
#include "TStyle.h"

void ROOT_settings(){
  // gROOT and gStyle settings
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadTopMargin   (0.02);
  gStyle->SetPadRightMargin (0.04);
  //gStyle->SetLegendBorderSize(0);

  //TGaxis::SetMaxDigits(2);
  
  
}

