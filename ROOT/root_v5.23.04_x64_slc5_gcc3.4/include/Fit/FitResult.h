// @(#)root/mathcore:$Id: FitResult.h 26866 2008-12-12 10:50:07Z moneta $
// Author: L. Moneta Wed Aug 30 11:05:34 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class FitResult

#ifndef ROOT_Fit_FitResult
#define ROOT_Fit_FitResult

#ifndef ROOT_Fit_IFunctionfwd
#include "Math/IFunctionfwd.h"
#endif
#ifndef ROOT_Fit_IParamFunctionfwd
#include "Math/IParamFunctionfwd.h"
#endif

#include <vector>
#include <string>
#include <cmath>
#include <cassert>

namespace ROOT { 

   namespace Math { 
      class Minimizer; 
   }


   namespace Fit { 

      class FitConfig; 
      class BinData;

//___________________________________________________________________________________
/** 
   class containg the result of the fit and all the related information 
   (fitted parameter values, error, covariance matrix and minimizer result information)
   Contains a pointer also to the fitted (model) function, modified with the fit parameter values.  
   When the fit is valid, it is constructed from a  Minimizer and a model function pointer 

   @ingroup FitMain
*/ 
class FitResult {

public: 

   typedef  ROOT::Math::IParamMultiFunction IModelFunction; 

   /** 
      Default constructor for an empty (non valid) fit result
   */ 
   FitResult (); 

   /**
      Construct from a Minimizer instance 
    */
   FitResult(ROOT::Math::Minimizer & min, const FitConfig & fconfig, const IModelFunction * f, bool isValid, unsigned int sizeOfData = 0, bool binFit = true, const ROOT::Math::IMultiGenFunction * chi2func = 0, bool minosErr = false, unsigned int ncalls = 0);

   /** 
      Copy constructor. 
   */ 
   FitResult(const FitResult &);

   /** 
      Assignment operator
   */ 
   FitResult & operator = (const FitResult & rhs);  

   /** 
      Destructor 
   */ 
   ~FitResult (); 


public: 


   /** minimization quantities **/

   /// minimizer type 
   const std::string & MinimizerType() const { return fMinimType; } 

   /// True if fit successful, otherwise false.
   bool IsValid() const { return fValid; }

   /// True if a fit result does not exist (even invalid) with parameter values 
   bool IsEmpty() const { return (fParams.size() == 0);  }
 
   /// Return value of the objective function (chi2 or likelihood) used in the fit
   double MinFcnValue() const { return fVal; } 

   ///Number of function calls to find minimum
   unsigned int NCalls() const { return fNCalls; }
   
   ///Expected distance from minimum 
   double Edm() const { return fEdm; }

   ///   get total number of parameters 
   unsigned int NTotalParameters() const { return fParams.size(); } 
   
   /// get total number of free parameters
   unsigned int NFreeParameters() const { return fNFree; }

   /// minimizer status code 
   int Status() const { return fStatus; } 
 
   /** fitting quantities **/

   /// Return pointer to model (fit) function with fitted parameter values.
   const IModelFunction * FittedFunction() const { return fFitFunc; }

   /// Chi2 fit value
   /// in case of likelihood must be computed ? 
   double Chi2() const { return fChi2; } 

   /// Number of degree of freedom
   unsigned int Ndf() const { return fNdf; } 

   /// p value of the fit (chi2 probability)
   double Prob() const;  

   /// parameter errors
   const std::vector<double> & Errors() const { return fErrors; }

   /// parameter values
   const std::vector<double> & Parameters() const { return fParams; }

   /// parameter value by index
   double Value(unsigned int i) const { return fParams[i]; }

   /// parameter error by index
   double Error(unsigned int i) const { 
      return (i < fErrors.size() ) ? fErrors[i] : 0; 
   } 

//    /// Minos  Errors 
//    const std::vector<std::pair<double, double> > MinosErrors() const; 

   /// lower Minos error
   double LowerError(unsigned int i) const { 
      return (i < fMinosErrors.size() ) ? fMinosErrors[i].first : fErrors[i]; 
   } 

   /// upper Minos error
   double UpperError(unsigned int i) const { 
      return (i < fMinosErrors.size() ) ? fMinosErrors[i].second : fErrors[i]; 
   }

   /// parameter global correlation coefficient 
   double GlobalCC(unsigned int i) const { 
      return (i < fGlobalCC.size() ) ? fGlobalCC[i] : -1; 
   } 


   /// retrieve covariance matrix element 
   double CovMatrix (unsigned int i, unsigned int j) const { 
      if ( i >= fErrors.size() || j >= fErrors.size() ) return 0; 
      if (fCovMatrix.size() == 0) return 0; // no matrix is available in case of non-valid fits
      if ( j < i ) 
         return fCovMatrix[j + i* (i+1) / 2];
      else 
         return fCovMatrix[i + j* (j+1) / 2];
   }

   /// retrieve correlation elements 
   double Correlation(unsigned int i, unsigned int j ) const { 
      if ( i >= fErrors.size() || j >= fErrors.size() ) return 0; 
      if (fCovMatrix.size() == 0) return 0; // no matrix is available in case of non-valid fits
      double tmp = CovMatrix(i,i)*CovMatrix(j,j); 
      return ( tmp > 0) ? CovMatrix(i,j)/ std::sqrt(tmp) : 0; 
   }
   
   /// fill covariance matrix elements using a generic symmetric matrix class implementing operator(i,j)
   /// the matrix must be previously allocates with right size (npar * npar) 
   template<class Matrix> 
   void GetCovarianceMatrix(Matrix & mat) { 
      int npar = fErrors.size();
      assert(fCovMatrix.size() == npar*(npar+1)/2);
      for (int i = 0; i< npar; ++i) { 
         for (int j = 0; j<=i; ++i) { 
            mat(i,j) = fCovMatrix[j + i*(i+1)/2 ];
         }
      }
   }

   /// fill a correlation matrix elements using a generic symmetric matrix class implementing operator(i,j)
   /// the matrix must be previously allocates with right size (npar * npar) 
   template<class Matrix> 
   void GetCorrelationMatrix(Matrix & mat) { 
      int npar = fErrors.size(); 
      assert(fCovMatrix.size() == npar*(npar+1)/2);
      for (int i = 0; i< npar; ++i) { 
         for (int j = 0; j<=i; ++i) { 
            double tmp = fCovMatrix[i * (i +3)/2 ] * fCovMatrix[ j * (j+3)/2 ]; 
            if (tmp < 0) 
               mat(i,j) = 0; 
            else 
               mat(i,j) = fCovMatrix[j + i*(i+1)/2 ] / std::sqrt(tmp); 
         }
      }
   }

   /**
      get confidence intervals for an array of n points x. 
      stride1 indicates the stride in the coordinate space while stride2 the stride in dimension space. 
      For 1-dim points : stride1=1, stride2=1
      for multi-dim points arranged as (x0,x1,...,xN,y0,....yN)          stride1=1      stride2=n
      for multi-dim points arraged  as (x0,y0,..,x1,y1,...,xN,yN,..)     stride1=ndim,  stride2=1
      
      the confidence interval are returned in the array ci
      cl is the desired confidedence interval value
      
    */
   void GetConfidenceIntervals(unsigned int n, unsigned int stride1, unsigned int stride2, const double * x,  double * ci, double cl=0.95 ) const;     

   /**
      evaluate confidence interval for the point specified in the passed data sets
      the confidence interval are returned in the array ci
      cl is the desired confidence interval value
    */
   void GetConfidenceIntervals(const BinData & data, double * ci, double cl=0.95 ) const;


   /// get index for parameter name (return -1 if not found)
   int Index(const std::string & name) const; 


   ///normalize errors using chi2/ndf for chi2 fits
   void NormalizeErrors();

   /// flag to chek if errors are normalized
   bool NormalizedErrors() { return fNormalized; }

   /// get confidence level given an array of x data points


   /// print the result and optionaly covariance matrix and correlations
   void Print(std::ostream & os, bool covmat = false) const;

   ///print error matrix and correlations
   void PrintCovMatrix(std::ostream & os) const; 

   /// query if a parameter is bound 
   bool IsParameterBound(unsigned int ipar) const; 

   /// query if a parameter is fixed 
   bool IsParameterFixed(unsigned int ipar) const; 

   /// get name of parameter 
   std::string GetParameterName(unsigned int ipar) const;

protected: 


private: 


   /// Return pointer non const pointer to model (fit) function with fitted parameter values.
   /// used by Fitter class 
   IModelFunction * ModelFunction()  { return fFitFunc; }
   void SetModelFunction(IModelFunction * func) { fFitFunc = func; }

   friend class Fitter; 


   bool fValid;             // flag for indicating valid fit
   bool fNormalized;        // flag for indicating is errors are normalized
   unsigned int fNFree;     // number of fit free parameters (total parameters are in size of parameter vector)  
   unsigned int fNdf;       // number of degree of freedom
   unsigned int fNCalls;    // number of function calls
   int fStatus;             // minimizer status code
   double fVal;             // minimum function value
   double fEdm;             // expected distance from mimimum
   double fChi2;            // fit chi2 value (different than fval in case of chi2 fits)
   IModelFunction * fFitFunc; // model function resulting  from the fit. It is given by Fitter but it is managed by FitResult
   std::vector<unsigned int>   fFixedParams; // list of fixed parameters
   std::vector<unsigned int>   fBoundParams; // list of limited parameters
   std::vector<double>         fParams;  // parameter values. Size is total number of parameters
   std::vector<double>         fErrors;  // errors 
   std::vector<double>         fCovMatrix;  // covariance matrix (size is npar*(npar+1)/2) where npar is total parameters
   std::vector<double>         fGlobalCC;   // global Correlation coefficient
   std::vector<std::pair<double,double> > fMinosErrors;   // vector contains the two Minos errors 
   std::string fMinimType;              // string indicating type of minimizer
   std::vector<std::string> fParNames;  // parameter names (only with FCN only fites, when fFitFunc=0)

}; 

   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_FitResult */
