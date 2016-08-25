// @(#)root/tmva $Id: MethodKNN.h 21630 2008-01-10 19:40:44Z brun $
// Author: Rustem Ospanov

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodKNN                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Analysis of k-nearest neighbor                                            *
 *                                                                                *
 * Author:                                                                        *
 *      Rustem Ospanov <rustem@fnal.gov> - U. of Texas at Austin, USA             *
 *                                                                                *
 * Copyright (c) 2007:                                                            *
 *      CERN, Switzerland                                                         * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      U. of Texas at Austin, USA                                                *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodKNN
#define ROOT_TMVA_MethodKNN

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodKNN                                                            //
//                                                                      //
// Analysis of k-nearest neighbor                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <map>

#ifndef ROOT_TMVA_MethodBase
#include "TMVA/MethodBase.h"
#endif
#ifndef ROOT_TMVA_KNN_ModulekNN
#include "TMVA/ModulekNN.h"
#endif

namespace TMVA
{   
   namespace kNN
   {
      class ModulekNN;
   }

   class MethodKNN : public MethodBase
   {
   public:

      MethodKNN(const TString& jobName, 
                const TString& methodTitle, 
                DataSet& theData,
                const TString& theOption = "KNN",
                TDirectory* theTargetDir = NULL);

      MethodKNN(DataSet& theData, 
                const TString& theWeightFile,  
                TDirectory* theTargetDir = NULL);
      
      virtual ~MethodKNN( void );
    
      void Train( void );

      Double_t GetMvaValue();

      using MethodBase::WriteWeightsToStream;
      using MethodBase::ReadWeightsFromStream;

      void WriteWeightsToStream(std::ostream& o) const;
      void WriteWeightsToStream(TFile& rf) const;

      void ReadWeightsFromStream(std::istream& istr);
      void ReadWeightsFromStream(TFile &rf);

      const Ranking* CreateRanking();

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

   private:

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // default initialisation called by all constructors
      void InitKNN( void );

      // create kd-tree (binary tree) structure
      void MakeKNN( void );

      // polynomial kernel weight function
      double PolKernel(double value) const;

   private:

      // number of events (sumOfWeights)
      Double_t fSumOfWeightsS;        // sum-of-weights for signal training events
      Double_t fSumOfWeightsB;        // sum-of-weights for background training events      

      kNN::ModulekNN *fModule;        //! module where all work is done

      Int_t fnkNN;            // number of k-nearest neighbors 
      Int_t fTreeOptDepth;    // number of binary tree levels used for optimization

      Float_t fScaleFrac;     // fraction of events used for scaling

      Bool_t fUseKernel;      // use polynomial kernel weight function
      Bool_t fTrim;           // set equal number of signal and background events

      kNN::EventVec fEvent;   // (untouched) events used for learning

      ClassDef(MethodKNN,0) // k Nearest Neighbour classifier
   };

} // namespace TMVA

#endif // MethodKNN
