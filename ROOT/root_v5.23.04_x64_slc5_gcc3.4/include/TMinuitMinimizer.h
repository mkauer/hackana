// @(#)root/minuit:$Id: TMinuitMinimizer.h 26866 2008-12-12 10:50:07Z moneta $
// Author: L. Moneta Wed Oct 25 16:28:55 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class TMinuitMinimizer

#ifndef ROOT_TMinuitMinimizer
#define ROOT_TMinuitMinimizer

#ifndef ROOT_Math_Minimizer
#include "Math/Minimizer.h"
#endif

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

class TMinuit; 



namespace ROOT { 

   namespace Minuit { 


      // enumeration specifying the type of TMinuit minimizers
      enum EMinimizerType { 
         kMigrad, 
         kSimplex, 
         kCombined, 
         kMigradImproved, 
         kScan
      };

   }
}

 

/** 
   TMinuitMinimizer class: minimizer implementation based on TMinuit.
*/ 
class TMinuitMinimizer  : public ROOT::Math::Minimizer {

public: 

   /** 
      Default constructor
   */ 
   TMinuitMinimizer ( ROOT::Minuit::EMinimizerType type = ROOT::Minuit::kMigrad); 

   /** 
      Constructor from a char * (used by PM)
   */ 
   TMinuitMinimizer ( const char * type ); 

   /** 
      Destructor (no operations)
   */ 
   ~TMinuitMinimizer (); 

private:
   // usually copying is non trivial, so we make this unaccessible

   /** 
      Copy constructor
   */ 
   TMinuitMinimizer(const TMinuitMinimizer &); 

   /** 
      Assignment operator
   */ 
   TMinuitMinimizer & operator = (const TMinuitMinimizer & rhs); 

public: 

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func); 

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func); 

   /// set free variable 
   virtual bool SetVariable(unsigned int ivar, const std::string & name, double val, double step); 

   /// set upper/lower limited variable (override if minimizer supports them )
   virtual bool SetLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double /* lower */, double /* upper */); 

#ifdef LATER
   /// set lower limit variable  (override if minimizer supports them )
   virtual bool SetLowerLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double lower );
   /// set upper limit variable (override if minimizer supports them )
   virtual bool SetUpperLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double upper );
#endif

   /// set fixed variable (override if minimizer supports them )
   virtual bool SetFixedVariable(unsigned int /* ivar */, const std::string & /* name */, double /* val */);  

   /// set the value of an existing variable 
   virtual bool SetVariableValue(unsigned int , double );
 
   /// method to perform the minimization
   virtual  bool Minimize(); 

   /// return minimum function value
   virtual double MinValue() const { return fMinVal; } 

   /// return expected distance reached from the minimum
   virtual double Edm() const { return fEdm; }

   /// return  pointer to X values at the minimum 
   virtual const double *  X() const { return &fParams.front(); }

   /// return pointer to gradient values at the minimum 
   virtual const double *  MinGradient() const { return 0; } // not available in Minuit2 

   /// number of function calls to reach the minimum 
   virtual unsigned int NCalls() const { return 0; } 

   /// this is <= Function().NDim() which is the total 
   /// number of variables (free+ constrained ones) 
   virtual unsigned int NDim() const { return fDim; }   

   /// number of free variables (real dimension of the problem) 
   /// this is <= Function().NDim() which is the total 
   virtual unsigned int NFree() const { return fNFree; }  

   /// minimizer provides error and error matrix
   virtual bool ProvidesError() const { return true; } 

   /// return errors at the minimum 
   virtual const double * Errors() const { return  &fErrors.front(); }

   /** return covariance matrices elements 
       if the variable is fixed the matrix is zero
       The ordering of the variables is the same as in errors
   */ 
   virtual double CovMatrix(unsigned int i, unsigned int j) const { 
      return fCovar[i + fDim* j]; 
   }

   /// minos error for variable i, return false if Minos failed
   virtual bool GetMinosError(unsigned int i, double & errLow, double & errUp); 

   /**
      scan a parameter i around the minimum. A minimization must have been done before, 
      return false if it is not the case
    */
   virtual bool Scan(unsigned int i, unsigned int &nstep, double * x, double * y, double xmin = 0, double xmax = 0); 

   /**
      find the contour points (xi,xj) of the function for parameter i and j around the minimum
      The contour will be find for value of the function = Min + ErrorUp();
    */
   virtual bool Contour(unsigned int i, unsigned int j, unsigned int & npoints, double *xi, double *xj); 


   virtual void PrintResults();

   /// return reference to the objective function
   ///virtual const ROOT::Math::IGenFunction & Function() const; 
   

protected: 

   /// implementation of FCN for Minuit
   static void Fcn( int &, double * , double & f, double * , int);
   /// implementation of FCN for Minuit when user provided gradient is used
   static void FcnGrad( int &, double * g, double & f, double * , int);

   /// reset 
   void DoClear(); 

private: 

   bool fUsed;
   bool fMinosRun; 
   unsigned int fDim; 
   unsigned int fNFree;
   unsigned int fStrategy;
   double fMinVal;
   double fEdm; 
   std::vector<double> fParams;
   std::vector<double> fErrors;
   std::vector<double> fCovar; 

   ROOT::Minuit::EMinimizerType fType; 
   TMinuit * fMinuit; 

   // need to have a static copy of the function 
   //NOTE: This is NOT thread safe.
   static ROOT::Math::IMultiGenFunction * fgFunc;

   static TMinuit * fgMinuit; 

   static bool fgUsed;  // flag to control if static instance has done minimization

   ClassDef(TMinuitMinimizer,1)  //Implementation of Minimizer interface using TMinuit 

}; 



#endif /* ROOT_TMinuitMinimizer */
