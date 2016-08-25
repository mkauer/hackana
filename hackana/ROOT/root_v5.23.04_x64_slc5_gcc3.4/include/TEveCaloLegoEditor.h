// @(#)root/eve:$Id: TEveCaloLegoEditor.h 25245 2008-08-25 21:44:09Z matevz $
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveCaloLegoEditor
#define ROOT_TEveCaloLegoEditor

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGComboBox;
class TEveGValuator;

class TEveCaloLego;

class TEveCaloLegoEditor : public TGedFrame
{
private:
   TEveCaloLegoEditor(const TEveCaloLegoEditor&);            // Not implemented
   TEveCaloLegoEditor& operator=(const TEveCaloLegoEditor&); // Not implemented
   TGComboBox*  MakeLabeledCombo(const char* name, Int_t off);

protected:
   TEveCaloLego      *fM; // Model object.

   TGCheckButton     *fTopViewUseMaxColor;
   TGColorSelect     *fTopViewTowerColor;
   TGColorSelect     *fGridColor;
   TGColorSelect     *fFontColor;
   TGColorSelect     *fPlaneColor;
   TGNumberEntry     *fTransparency;

   TEveGValuator     *fNZSteps;

   TGComboBox        *fProjection;
   TGComboBox        *f2DMode;
   TGComboBox        *fBoxMode;

   TGVerticalFrame   *fRebinFrame;
   TGCheckButton     *fAutoRebin;
   TEveGValuator     *fPixelsPerBin;
   TGCheckButton     *fNormalizeRebin;

   void               MakeRebinFrame();
public:
   TEveCaloLegoEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
         UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
   virtual ~TEveCaloLegoEditor() {}

   virtual void SetModel(TObject* obj);

   // Declare callback/slot methods
   void DoTopViewUseMaxColor();
   void DoTopViewTowerColor(Pixel_t color);
   void DoGridColor(Pixel_t color);
   void DoFontColor(Pixel_t color);
   void DoPlaneColor(Pixel_t color);
   void DoTransparency();

   void DoNZSteps();

   void DoProjection();
   void Do2DMode();
   void DoBoxMode();

   void DoAutoRebin();
   void DoPixelsPerBin();
   void DoNormalize();

   ClassDef(TEveCaloLegoEditor, 0); // GUI editor for TEveCaloLego.
};

#endif
