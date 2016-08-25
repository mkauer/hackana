// @(#)root/gl:$Id: TGLParametricEquationGL.h 21252 2007-12-07 01:39:32Z matevz $
// Author:  Matevz Tadel, Jun 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLParametricEquationGL
#define ROOT_TGLParametricEquationGL

#ifndef ROOT_TGLObject
#include "TGLObject.h"
#endif
#ifndef ROOT_TGLPlotPainter
#include "TGLPlotPainter.h"
#endif

class TGLRnrCtx;
class TGLParametricEquation;
class TH2;


class TGLParametricEquationGL : public TGLObject
{
private:
   TGLParametricEquationGL(const TGLParametricEquationGL&);            // Not implemented
   TGLParametricEquationGL& operator=(const TGLParametricEquationGL&); // Not implemented

protected:
   TGLParametricEquation  *fM;
   TGLPlotPainter         *fPlotPainter;

public:
   TGLParametricEquationGL();
   virtual ~TGLParametricEquationGL();

   virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
   virtual void   SetBBox();
   virtual void   DirectDraw(TGLRnrCtx & rnrCtx) const;

   virtual Bool_t KeepDuringSmartRefresh() const { return kFALSE; }

   // To support two-level selection
   // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
   // virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

   ClassDef(TGLParametricEquationGL, 0) // GL renderer for TGLParametricEquation
};

#endif
