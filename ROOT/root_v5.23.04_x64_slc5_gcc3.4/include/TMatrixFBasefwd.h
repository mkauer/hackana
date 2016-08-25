// @(#)root/matrix:$Id: TMatrixFBasefwd.h 20882 2007-11-19 11:31:26Z rdm $
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixFBasefwd
#define ROOT_TMatrixFBasefwd

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrixFBase                                                         //
//                                                                      //
//  Forward declaration of TMatrixTBase<Float_t>                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

template<class Element> class TMatrixTBase;
typedef TMatrixTBase<Float_t> TMatrixFBase;

#endif
