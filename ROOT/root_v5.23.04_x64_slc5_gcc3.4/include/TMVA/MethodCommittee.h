// @(#)root/tmva $Id: MethodCommittee.h 21630 2008-01-10 19:40:44Z brun $ 
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCommittee                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Boosting                                                                  *
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
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodCommittee
#define ROOT_TMVA_MethodCommittee

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodCommittee                                                      //
//                                                                      //
// Committee method                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#ifndef ROOT_TMVA_MethodBase
#include "TMVA/MethodBase.h"
#endif

class TTree;
class TH2;

namespace TMVA {

   class MethodCommittee : public MethodBase {

   public:

      // constructor for training and reading
      MethodCommittee( const TString& jobName, 
                       const TString& committeeTitle,                        
                       DataSet& theData,
                       const TString& committeeOptions,
                       Types::EMVA method,
                       const TString& methodOptions,
                       TDirectory* theTargetDir = 0 );

      // constructor for calculating Committee-MVA using previously generatad members
      MethodCommittee( DataSet& theData, 
                       const TString& theWeightFile,  
                       TDirectory* theTargetDir = NULL );
  
      virtual ~MethodCommittee( void );
    
      // overloaded members from MethodBase
      void WriteStateToFile() const;

      // the training
      void Train();

      using MethodBase::WriteWeightsToStream;
      using MethodBase::ReadWeightsFromStream;

      // write weights to file
      void WriteWeightsToStream( ostream& o ) const;

      // read weights from file
      void ReadWeightsFromStream( istream& istr );

      // write method specific histos to target file
      void WriteMonitoringHistosToFile( void ) const;

      // calculate the MVA value
      Double_t GetMvaValue();

      // apply the boost algorithm to a member in the committee
      Double_t Boost(  TMVA::IMethod*, UInt_t imember );

      // ranking of input variables
      const Ranking* CreateRanking();

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // accessors
      const std::vector<TMVA::IMethod*>& GetCommittee()    const { return fCommittee; }
      const std::vector<Double_t>&       GetBoostWeights() const { return fBoostWeights; }

      //return the individual relative variable importance 
      std::vector<Double_t> GetVariableImportance();
      Double_t GetVariableImportance( UInt_t ivar );

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

   private:

      // accessors
      std::vector<IMethod*>& GetCommittee()    { return fCommittee; }
      std::vector<Double_t>& GetBoostWeights() { return fBoostWeights; }

      // boosting algorithm (adaptive boosting)
      Double_t AdaBoost( IMethod* );
 
      // boosting as a random re-weighting
      Double_t Bagging( UInt_t imember);
  
      UInt_t                          fNMembers;        // number of members requested
      std::vector<IMethod*>           fCommittee;       // the collection of members
      std::vector<Double_t>           fBoostWeights;    // the weights applied in the individual boosts
      TString                         fBoostType;       // string specifying the boost type

      // options for the MVA method
      Types::EMVA                     fMemberType;      // the MVA method to be boosted
      TString                         fMemberOption;    // the options for that method

      Bool_t                          fUseMemberDecision;  // use binary information from IsSignal
      // use average classification from the members, or have the individual members 
      
      Bool_t                          fUseWeightedMembers; // in the committee weighted from AdaBoost
    

      // Init used in the various constructors
      void InitCommittee( void );

      //some histograms for monitoring
      TH1F*                           fBoostFactorHist; // weights applied in boosting
      TH2F*                           fErrFractHist;    // error fraction vs member number
      TTree*                          fMonitorNtuple;   // monitoring ntuple
      Int_t                           fITree      ;     // ntuple var: ith member
      Double_t                        fBoostFactor;     // ntuple var: boost weight
      Double_t                        fErrorFraction;   // ntuple var: misclassification error fraction 
      Int_t                           fNnodes;          // ntuple var: nNodes

      std::vector< Double_t >         fVariableImportance; // the relative importance of the different variables 

      ClassDef(MethodCommittee,0)  // Analysis of Boosted MVA methods
   };

} // namespace TMVA

#endif
