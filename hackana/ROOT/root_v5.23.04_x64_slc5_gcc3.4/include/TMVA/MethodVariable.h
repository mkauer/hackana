// @(#)root/tmva $Id: MethodVariable.h 21630 2008-01-10 19:40:44Z brun $
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodVariable                                                        *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Wrapper class for a single variable "MVA"; this is required for           *
 *      the evaluation of the single variable discrimination performance          *
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

#ifndef ROOT_TMVA_MethodVariable
#define ROOT_TMVA_MethodVariable

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodVariable                                                       //
//                                                                      //
// Wrapper class for a single variable "MVA"; this is required for      //
// the evaluation of the single variable discrimination performance     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMVA_MethodBase
#include "TMVA/MethodBase.h"
#endif

namespace TMVA {

   class MethodVariable : public MethodBase {

   public:

      MethodVariable( const TString& jobName, 
                      const TString& methodTitle, 
                      DataSet& theData,
                      const TString& theOption = "Variable",
                      TDirectory* theTargetDir = 0 );
  
      virtual ~MethodVariable( void );
    
      // training method
      void Train( void );

      using MethodBase::WriteWeightsToStream;
      using MethodBase::ReadWeightsFromStream;

      // write weights to file
      void WriteWeightsToStream( ostream& o ) const;

      // read weights from file
      void ReadWeightsFromStream( istream& istr );

      // calculate the MVA value
      Double_t GetMvaValue();

      // ranking of input variables
      const Ranking* CreateRanking() { return 0; }

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

   private:

      ClassDef(MethodVariable,0) // Wrapper class for a single variable "MVA" 
   };

} // namespace TMVA

#endif
