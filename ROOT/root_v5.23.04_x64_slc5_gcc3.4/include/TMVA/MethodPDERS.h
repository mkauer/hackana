// @(#)root/tmva $Id: MethodPDERS.h 21630 2008-01-10 19:40:44Z brun $
// Author: Andreas Hoecker, Yair Mahalalel, Joerg Stelzer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodPDERS                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Multidimensional Likelihood using the "Probability density estimator      *
 *      range search" (PDERS) method suggested in                                 *
 *      T. Carli and B. Koblitz, NIM A 501, 576 (2003)                            *
 *                                                                                *
 *      The multidimensional PDFs for signal and background are modeled           *
 *      by counting the events in the "vicinity" of a test point. The volume      *
 *      that describes "vicinity" is user-defined through the option string.      *
 *      A search method based on binary-trees is used to improve the selection    *
 *      efficiency of the volume search.                                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Yair Mahalalel  <Yair.Mahalalel@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
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

#ifndef ROOT_TMVA_MethodPDERS
#define ROOT_TMVA_MethodPDERS

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodPDERS                                                          //
//                                                                      //
// Multidimensional Likelihood using the "Probability density           //
// estimator range search" (PDERS) method                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMVA_MethodBase
#include "TMVA/MethodBase.h"
#endif
#ifndef ROOT_TMVA_BinarySearchTree
#include "TMVA/BinarySearchTree.h"
#endif
#ifndef ROOT_TMVA_TVector
#include "TVector.h"
#endif

namespace TMVA {

   class Volume;
   class Event;

   class MethodPDERS : public MethodBase {

   public:

      MethodPDERS( const TString& jobName,
                   const TString& methodTitle, 
                   DataSet& theData,
                   const TString& theOption,
                   TDirectory* theTargetDir = 0 );

      MethodPDERS( DataSet& theData,
                   const TString& theWeightFile,
                   TDirectory* theTargetDir = NULL );

      virtual ~MethodPDERS( void );

      // training method
      void Train( void );

      // write weights to file
      void WriteWeightsToStream( ostream& o ) const;
      void WriteWeightsToStream( TFile& rf ) const;

      // read weights from file
      void ReadWeightsFromStream( istream& istr );
      void ReadWeightsFromStream( TFile& istr );

      // calculate the MVA value
      Double_t GetMvaValue();

   public:

      // for root finder
      static Double_t IGetVolumeContentForRoot( Double_t );
      Double_t         GetVolumeContentForRoot( Double_t );

      // static pointer to this object
      static MethodPDERS* ThisPDERS( void ) { return fgThisPDERS; }

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

      Volume*      fHelpVolume; // auxiliary variable
      Int_t        fFcnCall;    // number of external function calls (RootFinder)

      // accessors
      BinarySearchTree* GetBinaryTreeSig( void ) const { return fBinaryTreeS; }
      BinarySearchTree* GetBinaryTreeBkg( void ) const { return fBinaryTreeB; }

      Double_t KernelEstimate( const Event&, std::vector<const BinarySearchTreeNode*>&, Volume& );
      Double_t ApplyKernelFunction( Double_t normalized_distance );
      Double_t KernelNormalization( Double_t pdf );
      Double_t GetNormalizedDistance( const TMVA::Event &base_event, 
                                      const BinarySearchTreeNode &sample_event, 
                                      Double_t *dim_normalization);
      Double_t NormSinc( Double_t x );
      Double_t LanczosFilter( Int_t level, Double_t x );

      // ranking of input variables
      const Ranking* CreateRanking() { return 0; }

   private:

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // calculate the averages of the input variables needed for adaptive training
      void CalcAverages();

      // create binary search trees for signal and background
      void CreateBinarySearchTrees( TTree* tree );

      // option
      TString fVolumeRange;   // option volume range
      TString fKernelString;  // option kernel estimator

      enum EVolumeRangeMode {
         kUnsupported = 0,
         kMinMax,
         kRMS,
         kAdaptive,
         kUnscaled,
         kkNN
      } fVRangeMode;

      enum EKernelEstimator {
         kBox = 0,
         kSphere,
         kTeepee,
         kGauss,
         kSinc3,     // the sinc enumerators must be consecutive and in order!
         kSinc5,
         kSinc7,
         kSinc9,
         kSinc11,
         kLanczos2,
         kLanczos3,
         kLanczos5,
         kLanczos8,
   kTrim
      } fKernelEstimator;

      BinarySearchTree*  fBinaryTreeS;   // binary tree for signal
      BinarySearchTree*  fBinaryTreeB;   // binary tree for background

      vector<Float_t>*   fDelta;         // size of volume
      vector<Float_t>*   fShift;         // volume center
      vector<Float_t>    fAverageRMS;    // average RMS of signal and background

      Float_t            fScaleS;        // weight for signal events
      Float_t            fScaleB;        // weight for background events
      Float_t            fDeltaFrac;     // fraction of RMS
      Double_t           fGaussSigma;    // size of Gauss in adaptive volume 
      Double_t           fGaussSigmaNorm;// size of Gauss in adaptive volume (normalised to dimensions)


      // input for adaptive volume adjustment
      Float_t            fNEventsMin;    // minimum number of events in adaptive volume
      Float_t            fNEventsMax;    // maximum number of events in adaptive volume
      Float_t            fMaxVIterations;// maximum number of iterations to adapt volume size
      Float_t            fInitialScale;  // initial scale for adaptive volume

      Bool_t             fInitializedVolumeEle; // is volume element initialized ?
      
      Int_t              fkNNMin;        // min number of events in kNN tree
      Int_t              fkNNMax;       // max number of events in kNN tree
      Int_t               fkNNTests;      // maximum number of iterations to adapt volume size
      
      Double_t           fMax_distance;  // maximum distance
      Bool_t            fPrinted;       // print
      Bool_t             fNormTree;      // binary-search tree is normalised

      void    SetVolumeElement ( void );
      Float_t RScalc           ( const Event& );
      Float_t GetError         ( Float_t countS, Float_t countB,
                                 Float_t sumW2S, Float_t sumW2B ) const;

      // this carrier
      static MethodPDERS* fgThisPDERS; // this pointer (required by root finder)
      void UpdateThis() { fgThisPDERS = this; }

      void InitPDERS( void );

      ClassDef(MethodPDERS,0) // Multi-dimensional probability density estimator range search (PDERS) method
   };

} // namespace TMVA

#endif // MethodPDERS_H
