// @(#)root/tmva $Id: VariableDecorrTransform.h 27320 2009-02-02 06:40:36Z brun $
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableDecorrTransform                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Decorrelation of input variables                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_VariableDecorrTransform
#define ROOT_TMVA_VariableDecorrTransform

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Decorrelation transformation of input variables                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMatrixD.h"

#ifndef ROOT_TMVA_VariableTransformBase
#include "TMVA/VariableTransformBase.h"
#endif

namespace TMVA {

   class VariableDecorrTransform : public VariableTransformBase {

   public:
  
      VariableDecorrTransform( std::vector<TMVA::VariableInfo>& );
      virtual ~VariableDecorrTransform( void );

      void   ApplyTransformation( Types::ESBType type = Types::kMaxSBType ) const;
      Bool_t PrepareTransformation( TTree* inputTree );

      void WriteTransformationToStream ( std::ostream& ) const;
      void ReadTransformationFromStream( std::istream& );

      virtual void PrintTransformation( ostream & o );

      // provides string vector describing explicit transformation
      std::vector<TString>* GetTransformationStrings( Types::ESBType type ) const;

      // writer of function code
      virtual void MakeFunction( std::ostream& fout, const TString& fncName, Int_t part );

   private:

      TMatrixD* fDecorrMatrix[2];     //! Decorrelation matrix [signal/background]

      void GetSQRMats( TTree* tr );
      void GetCovarianceMatrix( TTree* tr, Bool_t isSignal, TMatrixDBase* mat );

      ClassDef(VariableDecorrTransform,0) // Variable transformation: decorrelation
   };

} // namespace TMVA

#endif 


