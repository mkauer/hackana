// @(#)root/mathmore:$Id: GSLSimAnMinimizer.h 26866 2008-12-12 10:50:07Z moneta $
// Author: L. Moneta Wed Dec 20 17:16:32 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 * This library is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU General Public License        *
 * as published by the Free Software Foundation; either version 2     *
 * of the License, or (at your option) any later version.             *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   *
 * General Public License for more details.                           *
 *                                                                    *
 * You should have received a copy of the GNU General Public License  *
 * along with this library (see file COPYING); if not, write          *
 * to the Free Software Foundation, Inc., 59 Temple Place, Suite      *
 * 330, Boston, MA 02111-1307 USA, or contact the author.             *
 *                                                                    *
 **********************************************************************/

// Header file for class GSLSimAnMinimizer

#ifndef ROOT_Math_GSLSimAnMinimizer
#define ROOT_Math_GSLSimAnMinimizer



#ifndef ROOT_Math_Minimizer
#include "Math/Minimizer.h"
#endif


#ifndef ROOT_Math_IFunctionfwd
#include "Math/IFunctionfwd.h"
#endif

#ifndef ROOT_Math_IParamFunctionfwd
#include "Math/IParamFunctionfwd.h"
#endif



#ifndef ROOT_Math_GSLSimAnnealing
#include "Math/GSLSimAnnealing.h"
#endif

namespace ROOT { 

   namespace Math { 


//_____________________________________________________________________________________
/** 
   GSLSimAnMinimizer class for minimization using simulated annealing
   using the algorithm from 
   <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Simulated-Annealing.html">
   GSL</A>.
   It implements the ROOT::Minimizer interface and 
   a plug-in (name "GSLSimAn") exists to instantiate this class via the plug-in manager

   @ingroup Min1D
*/ 
class GSLSimAnMinimizer : public  ROOT::Math::Minimizer {

public: 

   /** 
      Default constructor
   */ 
   GSLSimAnMinimizer (int type = 0); 

   /** 
      Destructor (no operations)
   */ 
   ~GSLSimAnMinimizer ();  

private:
   // usually copying is non trivial, so we make this unaccessible

   /** 
      Copy constructor
   */ 
   GSLSimAnMinimizer(const GSLSimAnMinimizer &) : ROOT::Math::Minimizer() {} 

   /** 
      Assignment operator
   */ 
   GSLSimAnMinimizer & operator = (const GSLSimAnMinimizer & rhs)  {
      if (this == &rhs) return *this;  // time saving self-test
      return *this;
   }

public: 

   /// set the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGenFunction & func); 

   /// set gradient the function to minimize
   virtual void SetFunction(const ROOT::Math::IMultiGradFunction & func); 

   /// set free variable 
   virtual bool SetVariable(unsigned int ivar, const std::string & name, double val, double step); 

   /// set fixed variable (override if minimizer supports them )
   virtual bool SetFixedVariable(unsigned int /* ivar */, const std::string & /* name */, double /* val */);  

#ifdef LATER
   /// set lower limit variable  (override if minimizer supports them )
   virtual bool SetLowerLimitedVariable(unsigned int  ivar , const std::string & name , double val , double step , double lower );
   /// set upper limit variable (override if minimizer supports them )
   virtual bool SetUpperLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double upper ); 
   /// set upper/lower limited variable (override if minimizer supports them )
   virtual bool SetLimitedVariable(unsigned int ivar , const std::string & name , double val , double step , double /* lower */, double /* upper */); 
#endif

   /// set the value of an existing variable 
   virtual bool SetVariableValue(unsigned int ivar, double val );
   /// set the values of all existing variables (array must be dimensioned to the size of the existing parameters)
   virtual bool SetVariableValues(const double * x);


   /// method to perform the minimization
   virtual  bool Minimize(); 

   /// return minimum function value
   virtual double MinValue() const { return fMinVal; } 

   /// return expected distance reached from the minimum
   virtual double Edm() const { return 0; } // not impl. }

   /// return  pointer to X values at the minimum 
   virtual const double *  X() const { return &fValues.front(); } 

   /// return pointer to gradient values at the minimum 
   virtual const double *  MinGradient() const { return 0; } // not impl.  

   /// number of function calls to reach the minimum 
   virtual unsigned int NCalls() const { return 0; } // not yet ipl.  

   /// this is <= Function().NDim() which is the total 
   /// number of variables (free+ constrained ones) 
   virtual unsigned int NDim() const { return fDim; }   

   /// number of free variables (real dimension of the problem) 
   /// this is <= Function().NDim() which is the total 
   virtual unsigned int NFree() const { return fDim; }  

   /// minimizer provides error and error matrix
   virtual bool ProvidesError() const { return false; } 

   /// return errors at the minimum 
   virtual const double * Errors() const { return 0; }

   /** return covariance matrices elements 
       if the variable is fixed the matrix is zero
       The ordering of the variables is the same as in errors
   */ 
   virtual double CovMatrix(unsigned int , unsigned int ) const { return 0; }

   /// minos error for variable i, return false if Minos failed
   virtual bool GetMinosError(unsigned int , double & /* errLow */ , double & /* errUp */ ) { return false; }

   /// return reference to the objective function
   ///virtual const ROOT::Math::IGenFunction & Function() const; 


protected: 

private: 
   
   // dimension of the function to be minimized 
   unsigned int fDim; 
   // number of fixed variables 
   unsigned int fNFix; 

   ROOT::Math::GSLSimAnnealing  fSolver; 
   const ROOT::Math::IMultiGenFunction * fObjFunc; 
   
   double fMinVal; 

   mutable std::vector<double> fValues;

   std::vector<double> fSteps;
   std::vector<std::string> fNames;
   std::vector<bool> fVarFix; 

}; 

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_GSLSimAnMinimizer */
