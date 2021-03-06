/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOPROFILELL
#define ROOPROFILELL

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include <map>
#include <string>

class RooMinuit ;
 
class RooProfileLL : public RooAbsReal {
public:

  RooProfileLL(const char *name, const char *title, RooAbsReal& nll, const RooArgSet& observables);
  RooProfileLL(const RooProfileLL& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooProfileLL(*this,newname); }
  virtual ~RooProfileLL() ;

  void setAlwaysStartFromMin(Bool_t flag) { _startFromMin = flag ; }
  Bool_t alwaysStartFromMin() const { return _startFromMin ; }

  RooAbsReal& nll() { return const_cast<RooAbsReal&>(_nll.arg()) ; }
  const RooArgSet& bestFitParams() const ;

  virtual RooAbsReal* createProfile(const RooArgSet& paramsOfInterest) ;


protected:

  void validateAbsMin() const ;

  RooRealProxy _nll ;    // Input -log(L) function
  RooSetProxy _obs ;     // Parameters of profile likelihood
  RooSetProxy _par ;     // Marginialized parameters of likelihood
  Bool_t _startFromMin ; // Always start minimization for global minimum?

  TIterator* _piter ; //! Iterator over profile likelihood parameters to be minimized 
  TIterator* _oiter ; //! Iterator of profile likelihood output parameter(s)

  mutable RooMinuit* _minuit ; //! Internal minuit instance

  mutable Bool_t _absMinValid ; // flag if absmin is up-to-date
  mutable Double_t _absMin ; // absolute minimum of -log(L)
  mutable RooArgSet _paramAbsMin ; // Parameter values at absolute minimum
  mutable std::map<std::string,bool> _paramFixed ; // Parameter constant status at last time of use
  
  Double_t evaluate() const ;


private:

  ClassDef(RooProfileLL,0) // Real-valued function representing profile likelihood of external (likelihood) function
};
 
#endif
