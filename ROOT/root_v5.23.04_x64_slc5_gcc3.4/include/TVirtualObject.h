// @(#)root/io:$Id: TVirtualObject.h 25450 2008-09-18 21:13:42Z pcanal $
// Author: Philippe Canal July, 2008

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
  *************************************************************************/

#ifndef ROOT_TVirtualObject
#define ROOT_TVirtualObject


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualObject                                                       //
//                                                                      //
// Wrapper around an object and giving indirect access to its content    //
// even if the object is not of a class in the Cint/Reflex dictionary.  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClassRef.h"

class TVirtualObject {
public:
   TClassRef  fClass;
   void      *fObject;

   TVirtualObject( TClass *cl ) : fClass(cl), fObject( cl ? cl->New() : 0 ) {};
   ~TVirtualObject() { if (fClass) fClass->Destructor( fObject ); }

   TClass *GetClass() { return fClass; }
   void   *GetObject() { return fObject; }

};

#endif // ROOT_TVirtualObject
