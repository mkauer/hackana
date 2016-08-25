#ifndef DOLLIMIT_H
#define DOLLIMIT_H

#include <TH1.h>
#include <TMath.h>
#include <Math/RootFinderAlgorithms.h>
#include <vector.h>
#include <iostream>


double Eval2(double x, void* params);

class dollimit {

public:
  vector<double> s;
  vector<double> b;

  dollimit();
  ~dollimit(){};
  void fit();
  void Init(TH1*h2, TH1*h3,Int_t lx, Int_t ux);
  Double_t get_value(){return N;};

  double signal;
  double rightside;
  TH1 *sh;
  TH1 *bh;

  Double_t N;
  Double_t NL,NR;
  Bool_t initialised;
  
  ROOT::Math::Roots::Bisection RF;

};

#endif
