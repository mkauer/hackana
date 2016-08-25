// @(#)root/tmva $Id: Tools.h 27320 2009-02-02 06:40:36Z brun $   
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Tools                                                                 *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Global auxiliary applications and data treatment routines                 *
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

#ifndef ROOT_TMVA_Tools
#define ROOT_TMVA_Tools

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Tools (namespace)                                                    //
//                                                                      //
// Global auxiliary applications and data treatment routines            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#ifndef ROOT_TMVA_TMatrixDSymEigen
#include "TMatrixDSymEigen.h"
#endif
#ifndef ROOT_TMVA_MsgLogger
#include "TMVA/MsgLogger.h"
#endif

class TList;
class TTree;
class TString;
class TH1;
class TSpline;

namespace TMVA {

   class Event;
   class PDF;
   
   class Tools {
   private:
      Tools();

   public:

      // destructor
      ~Tools();

      // accessor to single instance
      static Tools& Instance() { return fgTools?*(fgTools): *(fgTools = new Tools()); }

      // simple statistics operations on tree entries
      void  ComputeStat( TTree* theTree, const TString theVarName,
                                Double_t&, Double_t&, Double_t&, 
                                Double_t&, Double_t&, Double_t&, Bool_t norm = kFALSE );

      // compute variance from sums
      inline Double_t ComputeVariance( Double_t sumx2, Double_t sumx, Int_t nx );

      // creates histograms normalized to one
      TH1* ProjNormTH1F( TTree* theTree, TString theVarName,
                                TString name, Int_t nbins, 
                                Double_t xmin, Double_t xmax, TString cut );

      // normalize histogram by its integral
      Double_t NormHist( TH1* theHist, Double_t norm = 1.0 );

      // parser for TString phrase with items separated by a character
      TList* ParseFormatLine( TString theString, const char * sep = ":" );

      // parse option string for ANN methods
      std::vector<Int_t>* ParseANNOptionString( TString theOptions, Int_t nvar,
                                                std::vector<Int_t>* nodes );
      
      // returns the square-root of a symmetric matrix: symMat = sqrtMat*sqrtMat
      TMatrixD* GetSQRootMatrix( TMatrixDSym* symMat );

      // turns covariance into correlation matrix
      const TMatrixD* GetCorrelationMatrix( const TMatrixD* covMat );

      // check spline quality by comparison with initial histogram
      Bool_t CheckSplines( const TH1*, const TSpline* );

      // normalization of variable output
      Double_t NormVariable( Double_t x, Double_t xmin, Double_t xmax );

      // return separation of two histograms or PDFs
      Double_t GetSeparation( const TH1& S, const TH1& B ) const;
      Double_t GetSeparation( const PDF& pdfS, const PDF& pdfB ) const;

      // vector rescaling
      std::vector<Double_t> MVADiff( std::vector<Double_t>&, std::vector<Double_t>& );
      void Scale( std::vector<Double_t>&, Double_t );
      void Scale( std::vector<Float_t>&,  Float_t  );
  
      // re-arrange a vector of arrays (vectors) in a way such that the first array
      // is ordered, and the other arrays reshuffeld accordingly
      void UsefulSortDescending( std::vector< std::vector<Double_t> >&, std::vector<TString>* vs = 0 );
      void UsefulSortAscending ( std::vector< std::vector<Double_t> >& );
    
      void UsefulSortDescending( std::vector<Double_t>& );
      void UsefulSortAscending ( std::vector<Double_t>& );

      Int_t GetIndexMaxElement ( std::vector<Double_t>& );
      Int_t GetIndexMinElement ( std::vector<Double_t>& );

      // check if input string contains regular expression
      Bool_t  ContainsRegularExpression( const TString& s );
      TString ReplaceRegularExpressions( const TString& s, const TString& replace = "+" );

      // routines for formatted output -----------------
      void FormattedOutput( const std::vector<Double_t>&, const std::vector<TString>&, 
                            const TString titleVars, const TString titleValues, MsgLogger& logger,
                            TString format = "%+1.3f" );
      void FormattedOutput( const TMatrixD&, const std::vector<TString>&, MsgLogger& logger );

      void WriteFloatArbitraryPrecision( Float_t  val, ostream& os );
      void ReadFloatArbitraryPrecision ( Float_t& val, istream& is );

      // check variable range and set var to lower or upper if out of range
      template<typename T>
      inline Bool_t VerifyRange( MsgLogger& mlog, const char *varstr, T& var, T vmin, T vmax );

      template<typename T>
      inline Bool_t VerifyRange( MsgLogger& mlog, const char *varstr, T& var, T vmin, T vmax, T vdef );

      template<typename T>
      inline Int_t VerifyRange( T& var, T vmin, T vmax );

      // output logger
      MsgLogger& Logger() const;

      const TString& Color( const TString& );

      // print welcome message (to be called from, eg, .TMVAlogon)
      enum EWelcomeMessage { kStandardWelcomeMsg = 1,
                             kIsometricWelcomeMsg,
                             kBlockWelcomeMsg,
                             kLeanWelcomeMsg,
                             kLogoWelcomeMsg,
                             kSmall1WelcomeMsg,
                             kSmall2WelcomeMsg,
                             kOriginalWelcomeMsgColor,
                             kOriginalWelcomeMsgBW };

      void TMVAWelcomeMessage();
      void TMVAWelcomeMessage( MsgLogger& logger, EWelcomeMessage m = kStandardWelcomeMsg );
      void TMVAVersionMessage( MsgLogger& logger );
      void ROOTVersionMessage( MsgLogger& logger );

      const   TString    fRegexp;
      mutable MsgLogger* fLogger;
      static  Tools*     fgTools;

   }; // Common tools

#if !defined(__CINT__) || defined(__MAKECINT__)
   Tools& gTools(); // global accessor
#endif

} // namespace TMVA

//_______________________________________________________________________
inline Double_t TMVA::Tools::ComputeVariance( Double_t sumx2, Double_t sumx, Int_t nx )
{
   // compute variance from given sums
   if (nx<2) return 0;
   return (sumx2 - ((sumx*sumx)/static_cast<Double_t>(nx)))/static_cast<Double_t>(nx-1);
}

//_______________________________________________________________________
template<typename T>
inline Int_t TMVA::Tools::VerifyRange( T& var, T vmin, T vmax )
{
   // check range and return +1 if above, -1 if below or 0 if inside
   if (var>vmax) return  1;
   if (var<vmin) return -1;
   return 0;
}
//_______________________________________________________________________
template<typename T>
inline Bool_t TMVA::Tools::VerifyRange( TMVA::MsgLogger& mlog, const char *varstr, T& var, T vmin, T vmax )
{
   // verify range and print out message
   // if outside range, set to closest limit
   Int_t dir = TMVA::Tools::VerifyRange(var,vmin,vmax);
   Bool_t modif=kFALSE;
   if (dir==1) {
      modif = kTRUE;
      var=vmax;
   }
   if (dir==-1) {
      modif = kTRUE;
      var=vmin;
   }
   if (modif) {
      mlog << kWARNING << "Option <" << varstr << "> " << (dir==1 ? "above":"below") << " allowed range. Reset to new value = " << var << Endl;
   }
   return modif;
}

//_______________________________________________________________________
template<typename T>
inline Bool_t TMVA::Tools::VerifyRange( TMVA::MsgLogger& mlog, const char *varstr, T& var, T vmin, T vmax, T vdef )
{
   // verify range and print out message
   // if outside range, set to given default value
   Int_t dir = TMVA::Tools::VerifyRange(var,vmin,vmax);
   Bool_t modif=kFALSE;
   if (dir!=0) {
      modif = kTRUE;
      var=vdef;
   }
   if (modif) {
      mlog << kWARNING << "Option <" << varstr << "> " << (dir==1 ? "above":"below") << " allowed range. Reset to default value = " << var << Endl;
   }
   return modif;
}

#endif

