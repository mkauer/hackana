// @(#)root/eve:$Id: TGLCameraOverlay.h 27643 2009-02-27 16:13:24Z matevz $
// Author: Alja Mrak-Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLCameraOverlay
#define ROOT_TGLCameraOverlay

#include "TAttAxis.h"
#include "TGLOverlay.h"
#include "TGLUtil.h"

class TGLAxisPainter;
class TGLFont;

class TAttAxis;
class TAxis;

class TGLCameraOverlay : public TGLOverlayElement
{
public:
   enum EMode { kPlaneIntersect, kBar, kAxis };

private:
   TGLCameraOverlay(const TGLCameraOverlay&);            // Not implemented
   TGLCameraOverlay& operator=(const TGLCameraOverlay&); // Not implemented

   Double_t       fFrustum[4];

protected:
   Bool_t         fShowOrthographic;
   Bool_t         fShowPerspective;

   EMode          fOrthographicMode;
   EMode          fPerspectiveMode;

   TGLAxisPainter *fAxisPainter;
   TAxis          *fAxis;
   Float_t         fAxisExtend;

   TGLPlane       fExternalRefPlane;
   Bool_t         fUseExternalRefPlane;

   void    RenderPlaneIntersect(TGLRnrCtx& rnrCtx);
   void    RenderAxis(TGLRnrCtx& rnrCtx);
   void    RenderBar(TGLRnrCtx& rnrCtx);

public:
   TGLCameraOverlay(Bool_t showOrtho=kTRUE, Bool_t showPersp=kFALSE);
   virtual ~TGLCameraOverlay();

   virtual  void   Render(TGLRnrCtx& rnrCtx);

   TGLPlane& RefExternalRefPlane() { return fExternalRefPlane; }
   void      UseExternalRefPlane(Bool_t x) { fUseExternalRefPlane=x; }
   Bool_t    GetUseExternalRefPlane() const { return fUseExternalRefPlane; }

   Int_t    GetPerspectiveMode() const { return fPerspectiveMode;}
   void     SetPerspectiveMode(EMode m) {fPerspectiveMode = m;}
   Int_t    GetOrthographicMode() const { return fOrthographicMode;}
   void     SetOrthographicMode(EMode m) {fOrthographicMode = m;}

   Bool_t   GetShowOrthographic() const { return fShowOrthographic; }
   void     SetShowOrthographic(Bool_t x) {fShowOrthographic =x;}
   Bool_t   GetShowPerspective() const { return fShowPerspective; }
   void     SetShowPerspective(Bool_t x) {fShowPerspective =x;}

   TAttAxis* GetAttAxis();

   ClassDef(TGLCameraOverlay, 1); // Show coorinates of current camera frustum.
};

#endif