// @(#)root/tmva $Id: RootFinder.h 20882 2007-11-19 11:31:26Z rdm $    
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RootFinder                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Root finding using Brents algorithm                                       *
 *      (translated from CERNLIB function RZERO)                                  *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Xavier Prudent  <prudent@lapp.in2p3.fr>  - LAPP, France                   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_RootFinder
#define ROOT_TMVA_RootFinder

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RootFinder                                                           //
//                                                                      //
// Root finding using Brents algorithm                                  //
// (translated from CERNLIB function RZERO)                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#ifndef ROOT_TMVA_MsgLogger
#include "TMVA/MsgLogger.h"
#endif

namespace TMVA {

   class RootFinder : public TObject {

   public:

      RootFinder( Double_t (*rootVal)( Double_t ),
                  Double_t rootMin, Double_t rootMax,
                  Int_t    maxIterations = 100, 
                  Double_t absTolerance  = 0.0 );
      virtual ~RootFinder( void );
      
      // returns the root of the function
      Double_t Root( Double_t refValue );

   private:

      Double_t fRootMin;  // minimum root value
      Double_t fRootMax;  // maximum root value
      Int_t    fMaxIter;  // maximum number of iterations
      Double_t fAbsTol;   // absolute tolerance deviation

      // function pointer
      Double_t (*fGetRootVal)( Double_t );

      mutable MsgLogger fLogger; // message logger

      ClassDef(RootFinder,0) // Root finding using Brents algorithm
   };

} // namespace TMVA

#endif