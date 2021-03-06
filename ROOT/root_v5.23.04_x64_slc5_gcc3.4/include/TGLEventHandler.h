// @(#)root/gl:$Id: TGLEventHandler.h 25245 2008-08-25 21:44:09Z matevz $
// Author: Bertrand Bellenot   29/01/2008

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLEventHandler
#define ROOT_TGLEventHandler

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLEventHandler                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGEventHandler
#include "TGEventHandler.h"
#endif
#ifndef ROOT_TPoint
#include "TPoint.h"
#endif

class TGToolTip;

class TGLViewer;
class TGLPhysicalShape;

class TGLEventHandler : public TGEventHandler
{
protected:
   TGLViewer        *fGLViewer;
   TTimer           *fMouseTimer;     // mouse delay timer
   TPoint            fLastPos;
   TPoint            fLastMouseOverPos;
   TGLPhysicalShape *fLastMouseOverShape;
   TGToolTip        *fTooltip;        // tooltip for highlight
   TPoint            fLastGlobalPos;
   TPoint            fTooltipPos;
   UInt_t            fActiveButtonID;
   UInt_t            fLastEventState;
   Bool_t            fInPointerGrab;
   Bool_t            fMouseTimerRunning;
   Bool_t            fTooltipShown;
   Int_t             fTooltipPixelTolerance;

   virtual Bool_t Rotate(Int_t xDelta, Int_t yDelta, Bool_t mod1, Bool_t mod2);

public:
   TGLEventHandler(const char *name, TGWindow *w, TObject *obj, const char *title="");
   virtual ~TGLEventHandler();

   virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual Bool_t HandleEvent(Event_t *event);
   virtual Bool_t HandleExpose(Event_t * event);
   virtual Bool_t HandleFocusChange(Event_t *event);
   virtual Bool_t HandleCrossing(Event_t *event);
   virtual Bool_t HandleButton(Event_t * event);
   virtual Bool_t HandleDoubleClick(Event_t *event);
   virtual Bool_t HandleConfigureNotify(Event_t *event);
   virtual Bool_t HandleKey(Event_t *event);
   virtual Bool_t HandleMotion(Event_t * event);
   virtual Bool_t HandleTimer(TTimer *t);
   virtual void   StartMouseTimer();
   virtual void   StopMouseTimer();
   virtual void   ClearMouseOver();
   virtual void   Repaint();

   virtual void   TriggerTooltip(const char* text);
   virtual void   RemoveTooltip();

   void SetMouseOverSelectDelay(Int_t ms);
   void SetMouseOverTooltipDelay(Int_t ms);

   Int_t GetTooltipPixelTolerance()  const { return fTooltipPixelTolerance; }
   void  SetTooltipPixelTolerance(Int_t t) { fTooltipPixelTolerance = t; }

   ClassDef(TGLEventHandler, 0); // Base-class and default implementation of event-handler for TGLViewer.
};

#endif

