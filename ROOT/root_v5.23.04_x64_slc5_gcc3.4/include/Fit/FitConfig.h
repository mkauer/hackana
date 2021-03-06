// @(#)root/mathcore:$Id: FitConfig.h 26508 2008-11-28 14:22:29Z moneta $
// Author: L. Moneta Thu Sep 21 16:21:29 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class FitConfig

#ifndef ROOT_Fit_FitConfig
#define ROOT_Fit_FitConfig


#ifndef ROOT_Fit_ParameterSettings
#include "Fit/ParameterSettings.h"
#endif

#ifndef ROOT_Fit_MinimizerControlParams
#include "Fit/MinimizerControlParams.h"
#endif

#ifndef ROOT_Math_IParamFunctionfwd
#include "Math/IParamFunctionfwd.h"
#endif


#include <vector>

namespace ROOT { 

   namespace Math { 

      class Minimizer;
      class MinimizerOptions; 
   }

   namespace Fit { 

//___________________________________________________________________________________
/** 
   Class describing the configuration of the fit, options and parameter settings
   using the ROOT::Fit::ParameterSettings class 

   @ingroup FitMain
*/ 
class FitConfig {

public: 

   /** 
      Default constructor
   */ 
   FitConfig (unsigned int npar = 0); 

   /** 
      Destructor 
   */ 
   ~FitConfig ();    

   /**
      get the parameter settings for the i-th parameter (const method)
   */
   const ParameterSettings & ParSettings(unsigned int i) const { return fSettings.at(i); }

   /**
      get the parameter settings for the i-th parameter (non-const method)
   */
   ParameterSettings & ParSettings(unsigned int i) { return fSettings.at(i); }

   /**
      get the vector of parameter settings  (const method)
   */
   const std::vector<ROOT::Fit::ParameterSettings> & ParamsSettings() const { return fSettings; }

   /**
      get the vector of parameter settings  (non-const method)
   */
   std::vector<ROOT::Fit::ParameterSettings> & ParamsSettings() { return fSettings; }


   /**
      set the parameter settings from a model function. 
      Create always new parameter setting list from a given model function  
   */
   void CreateParamsSettings(const ROOT::Math::IParamMultiFunction & func); 

   /**
      set the parameter settings from number of parameters and a vector of values and optionally step values. If there are not existing or number of parameters does not match existing one, create a new parameter setting list. 
   */
   void SetParamsSettings(unsigned int npar, const double * params, const double * vstep = 0); 


   /**
      create a new minimizer according to chosen configuration
   */
   ROOT::Math::Minimizer * CreateMinimizer(); 


   /**
      access to the minimizer  control parameter (const method) 
   */
   const MinimizerControlParams & MinimizerOptions() const { return fMinimizerOpts; } 

   /**
      access to the minimizer  control parameter (non const method) 
   */
   MinimizerControlParams & MinimizerOptions()  { return fMinimizerOpts; } 


#ifndef __CINT__   // this method fails on Windows
   /**
      set all the minimizer options using class MinimizerOptions
    */
   void SetMinimizerOptions(const ROOT::Math::MinimizerOptions & minopt); 
#endif

   
   /**
      set minimizer type 
   */
   void SetMinimizer(const char * type, const char * algo = 0) { 
      if (type) fMinimizerType = type; 
      if (algo) fMinimAlgoType = algo; 
   } 

   /**
      return type of minimizer package
   */
   const std::string & MinimizerType() const { return fMinimizerType; } 

   /**
      return type of minimizer algorithms 
   */
   const std::string & MinimizerAlgoType() const { return fMinimAlgoType; } 


   /**
      flag to check if resulting errors are be normalized according to chi2/ndf 
   */
   bool NormalizeErrors() const { return fNormErrors; } 

   ///do analysis for parabolic errors
   bool ParabErrors() const { return fParabErrors; }

   ///do minos errros analysis on all parameters
   bool MinosErrors() const { return fMinosErrors; }

   /**
      set the option to normalize the error on the result  according to chi2/ndf
   */
   void SetNormErrors(bool on) { fNormErrors= on; }

   ///set parabolic erros
   void SetParabErrors(bool on) { fParabErrors = on; } 

   ///set Minos erros
   void SetMinosErrors(bool on) { fMinosErrors = on; } 



   /**
      static function to control default minimizer type and algorithm
   */
   static void SetDefaultMinimizer(const char * type, const char * algo = 0); 


protected: 


private: 

   bool fNormErrors;    // flag for error normalization
   bool fParabErrors;   // get correct parabolic errors estimate (call Hesse after minimizing)  
   bool fMinosErrors;   // do full error analysis using Minos


   std::vector<ROOT::Fit::ParameterSettings> fSettings;  // vector with the parameter settings

   std::string fMinimizerType;  // minimizer type (MINUIT, MINUIT2, etc..)
   std::string fMinimAlgoType;  // algorithm type (MIGRAD, SIMPLEX, etc..)
   MinimizerControlParams fMinimizerOpts;   //minimizer control parameters

}; 

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_FitConfig */
