// @(#)root/tmva $Id: TNeuronInputSum.h 21630 2008-01-10 19:40:44Z brun $
// Author: Matt Jachowski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::TNeuronInputSum                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      TNeuron input calculator -- calculates the weighted sum of inputs.        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Matt Jachowski  <jachowski@stanford.edu> - Stanford University, USA       *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/
 

#ifndef ROOT_TMVA_TNeuronInputSum
#define ROOT_TMVA_TNeuronInputSum

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TNeuronInputSum                                                      //
//                                                                      //
// TNeuron input calculator -- calculates the weighted sum of inputs    //
//
// Begin_Latex
// NeuronInput = #sum_{i} y_{i}
// End_Latex
//
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"

#ifndef ROOT_TMVA_TNeuronInput
#include "TMVA/TNeuronInput.h"
#endif
#ifndef ROOT_TMVA_TNeuron
#include "TMVA/TNeuron.h"
#endif

namespace TMVA {

   class TNeuronInputSum : public TNeuronInput {
    
   public:

      TNeuronInputSum() {}
      virtual ~TNeuronInputSum() {}

      // calculate input value for neuron
      Double_t GetInput( const TNeuron* neuron ) const {
         if (neuron->IsInputNeuron()) return 0;
         Double_t result = 0;
         for (Int_t i=0; i < neuron->NumPreLinks(); i++) {
            result += neuron->PreLinkAt(i)->GetWeightedValue();
         }
         return result;
      }

      // name of class
      TString GetName() { return "Sum of weighted activations"; }

      ClassDef(TNeuronInputSum,0) // Calculates weighted sum of neuron inputs
   };

} // namespace TMVA

#endif
