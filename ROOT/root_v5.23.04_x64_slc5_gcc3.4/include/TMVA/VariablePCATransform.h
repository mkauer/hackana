// @(#)root/tmva $Id: VariablePCATransform.h 21630 2008-01-10 19:40:44Z brun $
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariablePCATransform                                                  *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Principal value composition of input variables                            *
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

#ifndef ROOT_TMVA_VariablePCATransform
#define ROOT_TMVA_VariablePCATransform

#include "TPrincipal.h"

#ifndef ROOT_TMVA_VariableTransformBase
#include "TMVA/VariableTransformBase.h"
#endif

namespace TMVA {

   class VariablePCATransform : public VariableTransformBase {

   public:
  
      VariablePCATransform( std::vector<TMVA::VariableInfo>& );
      virtual ~VariablePCATransform( void );

      void   ApplyTransformation( Types::ESBType type = Types::kMaxSBType ) const;
      Bool_t PrepareTransformation( TTree* inputTree );

      void WriteTransformationToStream ( std::ostream& ) const;
      void ReadTransformationFromStream( std::istream& );

      // provides string vector describing explicit transformation
      std::vector<TString>* GetTransformationStrings( Types::ESBType type = Types::kMaxSBType ) const;

      // writer of function code
      virtual void MakeFunction( std::ostream& fout, const TString& fncName, Int_t part );

   private:

      void CalculatePrincipalComponents( TTree* originalTree );
      void X2P( const Double_t*, Double_t*, Int_t index ) const;

      TPrincipal* fPCA[2];        //! PCA [signal/background]
      
      // store relevant parts of PCA locally
      TVectorD* fMeanValues[2];   // mean values
      TMatrixD* fEigenVectors[2]; // eigenvectors

      ClassDef(VariablePCATransform,0) // Variable transformation: Principal Value Composition
   };

} // namespace TMVA

#endif 


