// @(#)root/eve:$Id: TEveTrackProjected.h 25422 2008-09-16 20:50:49Z matevz $
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveTrackProjected
#define ROOT_TEveTrackProjected

#include "TEveTrack.h"
#include "TEveProjectionBases.h"

class TEveProjection;

class TEveTrackProjected : public TEveTrack,
                           public TEveProjected
{
   friend class TEveTrackProjectedGL;

private:
   TEveTrackProjected(const TEveTrackProjected&);            // Not implemented
   TEveTrackProjected& operator=(const TEveTrackProjected&); // Not implemented

   Int_t GetBreakPointIdx(Int_t start);
   void  GetBreakPoint(Int_t N, Bool_t back, Float_t& x, Float_t& y, Float_t& z);

   TEveVector*          fOrigPnts;     // original track points

protected:
   std::vector<Int_t>   fBreakPoints; // indices of track break-points
   TEveProjection      *fProjection;  // projection

public:
   TEveTrackProjected();
   virtual ~TEveTrackProjected() {}

   virtual void SetProjection(TEveProjectionManager* mng, TEveProjectable* model);

   virtual void SetDepth(Float_t d);

   virtual void UpdateProjection();
   virtual void MakeTrack(Bool_t recurse=kTRUE);

   void         PrintLineSegments();

   virtual void SecSelected(TEveTrack*); // marked as signal in TEveTrack

   ClassDef(TEveTrackProjected, 1); // Projected copy of a TEveTrack.
};


/******************************************************************************/
// TEveTrackListProjected
/******************************************************************************/

class TEveTrackListProjected : public TEveTrackList,
                               public TEveProjected
{
private:
   TEveTrackListProjected(const TEveTrackListProjected&);            // Not implemented
   TEveTrackListProjected& operator=(const TEveTrackListProjected&); // Not implemented

public:
   TEveTrackListProjected();
   virtual ~TEveTrackListProjected() {}

   virtual void SetProjection(TEveProjectionManager* proj, TEveProjectable* model);
   virtual void SetDepth(Float_t d);
   virtual void UpdateProjection()  {}

   virtual void SetDepth(Float_t d, TEveElement* el);

   ClassDef(TEveTrackListProjected, 1); // Specialization of TEveTrackList for holding TEveTrackProjected objects.
};

#endif
