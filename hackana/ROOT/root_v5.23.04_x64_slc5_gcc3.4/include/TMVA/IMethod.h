// @(#)root/tmva $Id: IMethod.h 26050 2008-11-01 09:18:41Z brun $   
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : IMethod                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Interface for all concrete MVA method implementations                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
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

#ifndef ROOT_TMVA_IMethod
#define ROOT_TMVA_IMethod

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// IMethod                                                              //
//                                                                      //
// Interface for all concrete MVA method implementations                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

#ifndef ROOT_TMVA_Types
#include "TMVA/Types.h"
#endif

namespace TMVA {

   class Ranking;
   
   class IMethod : public virtual TObject {
      
   public:
      
      // default constructur
      IMethod() : TObject() {}
      
      // default destructur
      virtual ~IMethod() {}

      // ------- virtual member functions to be implemented by each MVA method

      // calculate the MVA value
      virtual Double_t GetMvaValue() = 0;

      // training method
      virtual void Train( void ) = 0;

      // write weights to output stream
      virtual void WriteWeightsToStream( std::ostream& ) const = 0;

      // read weights from output stream
      virtual void ReadWeightsFromStream( std::istream& ) = 0;
      
      // write method specific monitoring histograms to target file
      virtual void WriteMonitoringHistosToFile( void ) const = 0;

      // make ROOT-independent C++ class for classifier response
      virtual void MakeClass( const TString& classFileName = "" ) const = 0;

      // create ranking
      virtual const Ranking* CreateRanking() = 0;

      // print help message
      virtual void PrintHelpMessage() const = 0;

      // get specific help message from classifer
      virtual void GetHelpMessage() const = 0;

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      virtual void MakeClassSpecific( std::ostream&, const TString& ) const = 0;

      ClassDef(IMethod,0) // Method Interface

   };
} // namespace TMVA

#endif
