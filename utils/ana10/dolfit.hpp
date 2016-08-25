#ifndef DOLFIT_H
#define DOLFIT_H

#include <TH1.h>
#include <TMath.h>
#include <Math/RootFinderAlgorithms.h>
#include <vector.h>
//#include <vector>
#include <iostream>


double Eval(double x, void* params);
double LEval(double x, void* params);

class dolfit {

public:
  vector <double> s;
  vector <double> b;
  vector <double> d;
  //vector s;
  //vector b;
  //vector d;
  
  dolfit();
  ~dolfit(){};
  void fit();
  void Init(TH1*h1, TH1*h2, TH1*h3,Int_t lx, Int_t ux, Int_t rebin =1);
  Double_t get_value(){return N;};
  Double_t get_lefterror(){return N-NL;};
  Double_t get_righterror(){return NR-N;};

  double signal;
  double data;
  double rightside;
  TH1 *sh;
  TH1 *bh;
  TH1 *dh;

  Double_t N;
  Double_t NL,NR;
  Bool_t initialised;
  
  ROOT::Math::Roots::Bisection RF;

};

#endif
