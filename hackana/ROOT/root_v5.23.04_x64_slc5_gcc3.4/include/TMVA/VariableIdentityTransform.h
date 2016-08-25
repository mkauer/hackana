// @(#)root/tmva $Id: VariableIdentityTransform.h 21630 2008-01-10 19:40:44Z brun $
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : VariableIdentityTransform                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Identity transform                                                        *
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

#ifndef ROOT_TMVA_VariableIdentityTransform
#define ROOT_TMVA_VariableIdentityTransform

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// VariableIdentityTransform                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMVA_VariableTransformBase
#include "TMVA/VariableTransformBase.h"
#endif

namespace TMVA {

   class VariableIdentityTransform : public VariableTransformBase {

   public:
  
      VariableIdentityTransform( std::vector<TMVA::VariableInfo>& );
      virtual ~VariableIdentityTransform( void ) {}

      void   ApplyTransformation( Types::ESBType type = Types::kMaxSBType ) const;
      Bool_t PrepareTransformation( TTree* inputTree );

      void WriteTransformationToStream ( std::ostream& ) const {}
      void ReadTransformationFromStream( std::istream& ) { SetCreated(); }

      virtual TMVA::Event& GetEvent()  const { return GetEventRaw(); }

      // provides string vector describing explicit transformation
      std::vector<TString>* GetTransformationStrings( Types::ESBType ) const;

      // writer of function code
      virtual void MakeFunction(std::ostream& fout, const TString& fncName, Int_t part);

      ClassDef(VariableIdentityTransform,0) // Variable transformation: identity
   };

} // namespace TMVA

#endif 


