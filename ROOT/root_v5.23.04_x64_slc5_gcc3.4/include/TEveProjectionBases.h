// @(#)root/eve:$Id: TEveProjectionBases.h 24004 2008-05-24 20:08:56Z matevz $
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveProjectionBases
#define ROOT_TEveProjectionBases

#include "TEveUtil.h"

#include <list>

class TBuffer3D;

class TEveElement;

class TEveProjected;
class TEveProjectionManager;

////////////////////////////////////////////////////////////////
//                                                            //
// TEveProjectable                                            //
//                                                            //
// Abstract base class for non-linear projectable objects.    //
//                                                            //
////////////////////////////////////////////////////////////////

class TEveProjectable
{
private:
   TEveProjectable(const TEveProjectable&);            // Not implemented
   TEveProjectable& operator=(const TEveProjectable&); // Not implemented

protected:
   typedef std::list<TEveProjected*>            ProjList_t;
   typedef std::list<TEveProjected*>::iterator  ProjList_i;

   ProjList_t       fProjectedList; // references to projected instances.

public:
   TEveProjectable();
   virtual ~TEveProjectable();

   virtual TClass* ProjectedClass() const = 0;

   virtual Bool_t HasProjecteds() const { return ! fProjectedList.empty(); }

   virtual void AddProjected(TEveProjected* p)    { fProjectedList.push_back(p); }
   virtual void RemoveProjected(TEveProjected* p) { fProjectedList.remove(p);    }

   virtual void AddProjectedsToSet(std::set<TEveElement*>& set);

   virtual void PropagateVizParams(TEveElement* el=0);
   virtual void PropagateRenderState(Bool_t rnr_self, Bool_t rnr_children);
   virtual void PropagateMainColor(Color_t color, Color_t old_color);

   ClassDef(TEveProjectable, 0); // Abstract base class for classes that can be transformed with non-linear projections.
};


////////////////////////////////////////////////////////////////
//                                                            //
// TEveProjected                                              //
//                                                            //
// Abstract base class for non-linear projected objects.      //
//                                                            //
////////////////////////////////////////////////////////////////

class TEveProjected
{
private:
   TEveProjected(const TEveProjected&);            // Not implemented
   TEveProjected& operator=(const TEveProjected&); // Not implemented

protected:
   TEveProjectionManager *fManager;       // manager
   TEveProjectable       *fProjectable;   // link to original object
   Float_t                fDepth;         // z coordinate

public:
   TEveProjected();
   virtual ~TEveProjected();

   TEveProjectable* GetProjectable() const { return fProjectable; }

   virtual void SetProjection(TEveProjectionManager* mng, TEveProjectable* model);
   virtual void UnRefProjectable(TEveProjectable* assumed_parent);

   virtual void SetDepth(Float_t d) = 0;

   virtual void UpdateProjection() = 0;

   void SetDepthCommon(Float_t d, TEveElement* el, Float_t* bbox);

   ClassDef(TEveProjected, 0); // Abstract base class for classes that hold results of a non-linear projection transformation.
};

#endif
