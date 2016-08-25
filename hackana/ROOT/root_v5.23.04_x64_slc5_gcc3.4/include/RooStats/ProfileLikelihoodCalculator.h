// @(#)root/roostats:$Id: ProfileLikelihoodCalculator.h 26964 2008-12-16 16:30:01Z moneta $
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_ProfileLikelihoodCalculator
#define ROOSTATS_ProfileLikelihoodCalculator

#ifndef ROOSTATS_CombinedCalculator
#include "RooStats/CombinedCalculator.h"
#endif


namespace RooStats {

   class ProfileLikelihoodCalculator : public CombinedCalculator {
   public:
      ProfileLikelihoodCalculator();

      ProfileLikelihoodCalculator(RooWorkspace& ws, RooAbsData& data, RooAbsPdf& pdf, RooArgSet& paramsOfInterest, 
                                  Double_t size = 0.05, RooArgSet* nullParams = 0, RooArgSet* altParams = 0);

      ProfileLikelihoodCalculator(RooAbsData& data, RooAbsPdf& pdf, RooArgSet& paramsOfInterest, 
                                  Double_t size = 0.05, RooArgSet* nullParams = 0, RooArgSet* altParams = 0);



      virtual ~ProfileLikelihoodCalculator();
    
      // main interface, implemented
      virtual ConfInterval* GetInterval() const ; 
      // main interface, implemented
      virtual HypoTestResult* GetHypoTest() const;   
    

   protected:
      ClassDef(ProfileLikelihoodCalculator,1) // A concrete implementation of CombinedCalculator that uses the ProfileLikelihood ratio.
   };
}
#endif
