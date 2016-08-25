// @(#)root/eve:$Id: TEveCompound.h 27157 2009-01-15 14:05:12Z brun $
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveCompound
#define ROOT_TEveCompound

#include "TEveElement.h"
#include "TEveProjectionBases.h"


//==============================================================================
// TEveCompound
//==============================================================================

class TEveCompound : public TEveElementList,
                     public TEveProjectable
{
private:
   TEveCompound(const TEveCompound&);            // Not implemented
   TEveCompound& operator=(const TEveCompound&); // Not implemented

protected:
   Short_t  fCompoundOpen; // If more than zero, tag new children as compound members.

public:
   TEveCompound(const char* n="TEveCompound", const char* t="",
                Bool_t doColor=kTRUE);
   virtual ~TEveCompound() {}

   void   OpenCompound()         { ++fCompoundOpen;  }
   void   CloseCompound()        { --fCompoundOpen; }
   Bool_t IsCompoundOpen() const { return fCompoundOpen > 0; }

   virtual void SetMainColor(Color_t color);

   virtual void AddElement(TEveElement* el);
   virtual void RemoveElementLocal(TEveElement* el);
   virtual void RemoveElementsLocal();

   virtual void FillImpliedSelectedSet(Set_t& impSelSet);

   virtual TClass* ProjectedClass() const;

   ClassDef(TEveCompound, 0); // Container for managing compounds of TEveElements.
};


//==============================================================================
// TEveCompoundProjected
//==============================================================================

class TEveCompoundProjected : public TEveCompound,
                              public TEveProjected
{
private:
   TEveCompoundProjected(const TEveCompoundProjected&);            // Not implemented
   TEveCompoundProjected& operator=(const TEveCompoundProjected&); // Not implemented

protected:

public:
   TEveCompoundProjected();
   virtual ~TEveCompoundProjected() {}

   virtual void SetMainColor(Color_t color);

   // Abstract from TEveProjected, we seem not to care.
   virtual void SetDepth(Float_t /*d*/) {}
   virtual void UpdateProjection()      {}

   ClassDef(TEveCompoundProjected, 0); // Projected TEveCompund container.
};

#endif
