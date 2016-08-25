// @(#)root/eve:$Id: TEveRGBAPalette.h 21215 2007-12-05 17:19:23Z matevz $
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveRGBAPalette
#define ROOT_TEveRGBAPalette

#include "TEveUtil.h"

#include "TObject.h"

class TEveRGBAPalette : public TObject, public TEveRefCnt
{
   friend class TEveRGBAPaletteEditor;
   friend class TEveRGBAPaletteSubEditor;

public:
   enum ELimitAction_e { kLA_Cut, kLA_Mark, kLA_Clip, kLA_Wrap };

private:
   TEveRGBAPalette(const TEveRGBAPalette&);            // Not implemented
   TEveRGBAPalette& operator=(const TEveRGBAPalette&); // Not implemented

protected:
   Int_t     fLowLimit;  // Low  limit for Min/Max values (used by editor)
   Int_t     fHighLimit; // High limit for Min/Max values (used by editor)
   Int_t     fMinVal;
   Int_t     fMaxVal;
   Int_t     fNBins;

   Bool_t    fInterpolate;
   Bool_t    fShowDefValue;
   Int_t     fUnderflowAction;
   Int_t     fOverflowAction;

   Color_t   fDefaultColor;   // Color for when value is not specified
   UChar_t   fDefaultRGBA[4];
   Color_t   fUnderColor;     // Underflow color
   UChar_t   fUnderRGBA[4];
   Color_t   fOverColor;      // Overflow color
   UChar_t   fOverRGBA[4];

   mutable UChar_t* fColorArray; //[4*fNBins]

   void SetupColor(Int_t val, UChar_t* pix) const;

   static TEveRGBAPalette* fgDefaultPalette;

public:
   TEveRGBAPalette();
   TEveRGBAPalette(Int_t min, Int_t max, Bool_t interp=kFALSE, Bool_t showdef=kTRUE);
   virtual ~TEveRGBAPalette();

   void SetupColorArray() const;
   void ClearColorArray();

   Bool_t   WithinVisibleRange(Int_t val) const;
   const UChar_t* ColorFromValue(Int_t val) const;
   void     ColorFromValue(Int_t val, UChar_t* pix, Bool_t alpha=kTRUE) const;
   Bool_t   ColorFromValue(Int_t val, Int_t defVal, UChar_t* pix, Bool_t alpha=kTRUE) const;

   Int_t  GetMinVal() const { return fMinVal; }
   Int_t  GetMaxVal() const { return fMaxVal; }

   void   SetLimits(Int_t low, Int_t high);
   void   SetLimitsScaleMinMax(Int_t low, Int_t high);
   void   SetMinMax(Int_t min, Int_t max);
   void   SetMin(Int_t min);
   void   SetMax(Int_t max);

   Int_t  GetLowLimit()  const { return fLowLimit;  }
   Int_t  GetHighLimit() const { return fHighLimit; }

   // ================================================================

   Bool_t GetInterpolate() const { return fInterpolate; }
   void   SetInterpolate(Bool_t b);

   Bool_t GetShowDefValue() const { return fShowDefValue; }
   void   SetShowDefValue(Bool_t v) { fShowDefValue = v; }

   Int_t GetUnderflowAction() const  { return fUnderflowAction; }
   Int_t GetOverflowAction()  const  { return fOverflowAction;  }
   void  SetUnderflowAction(Int_t a) { fUnderflowAction = a;    }
   void  SetOverflowAction(Int_t a)  { fOverflowAction  = a;    }

   // ================================================================

   Color_t  GetDefaultColor() const { return fDefaultColor; }
   Color_t* PtrDefaultColor() { return &fDefaultColor; }
   UChar_t* GetDefaultRGBA()  { return fDefaultRGBA;  }
   const UChar_t* GetDefaultRGBA() const { return fDefaultRGBA;  }

   void   SetDefaultColor(Color_t ci);
   void   SetDefaultColor(Pixel_t pix);
   void   SetDefaultColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

   // ----------------------------------------------------------------

   Color_t  GetUnderColor() const { return fUnderColor; }
   Color_t* PtrUnderColor() { return &fUnderColor; }
   UChar_t* GetUnderRGBA()  { return fUnderRGBA;  }
   const UChar_t* GetUnderRGBA() const { return fUnderRGBA;  }

   void   SetUnderColor(Color_t ci);
   void   SetUnderColor(Pixel_t pix);
   void   SetUnderColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

   // ----------------------------------------------------------------

   Color_t  GetOverColor() const { return fOverColor; }
   Color_t* PtrOverColor() { return &fOverColor; }
   UChar_t* GetOverRGBA()  { return fOverRGBA;  }
   const UChar_t* GetOverRGBA() const { return fOverRGBA;  }

   void   SetOverColor(Color_t ci);
   void   SetOverColor(Pixel_t pix);
   void   SetOverColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a=255);

   // ================================================================

   // ?? Should we emit some *SIGNALS* ??
   // ?? Should we have a RendererTimeStamp ??

   ClassDef(TEveRGBAPalette, 1); // A generic, speed-optimised mapping from value to RGBA color supporting different wrapping and range truncation modes.
};


/******************************************************************************/
// Inlines for TEveRGBAPalette
/******************************************************************************/

//______________________________________________________________________________
inline Bool_t TEveRGBAPalette::WithinVisibleRange(Int_t val) const
{
   if ((val < fMinVal && fUnderflowAction == kLA_Cut) ||
       (val > fMaxVal && fOverflowAction  == kLA_Cut))
      return kFALSE;
   else
      return kTRUE;
}

//______________________________________________________________________________
inline const UChar_t* TEveRGBAPalette::ColorFromValue(Int_t val) const
{
   // Here we expect that kLA_Cut has been checked; we further check
   // for kLA_Wrap and kLA_Clip otherwise we proceed as for kLA_Mark.

   if (!fColorArray)  SetupColorArray();
   if (val < fMinVal) {
      if (fUnderflowAction == kLA_Wrap)
         val = (val+1-fMinVal)%fNBins + fMaxVal;
      else if (fUnderflowAction == kLA_Clip)
         val = fMinVal;
      else
         return fUnderRGBA;
   }
   else if(val > fMaxVal) {
      if (fOverflowAction == kLA_Wrap)
         val = (val-1-fMaxVal)%fNBins + fMinVal;
      else if (fOverflowAction == kLA_Clip)
         val = fMaxVal;
      else
         return fOverRGBA;
   }
   return fColorArray + 4 * (val - fMinVal);
}

//______________________________________________________________________________
inline void TEveRGBAPalette::ColorFromValue(Int_t val, UChar_t* pix, Bool_t alpha) const
{
   const UChar_t* c = ColorFromValue(val);
   pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2];
   if (alpha) pix[3] = c[3];
}

//______________________________________________________________________________
inline Bool_t TEveRGBAPalette::ColorFromValue(Int_t val, Int_t defVal, UChar_t* pix, Bool_t alpha) const
{
   if (val == defVal) {
      if (fShowDefValue) {
         pix[0] = fDefaultRGBA[0];
         pix[1] = fDefaultRGBA[1];
         pix[2] = fDefaultRGBA[2];
         if (alpha) pix[3] = fDefaultRGBA[3];
         return kTRUE;
      } else {
         return kFALSE;
      }
   }

   if (WithinVisibleRange(val)) {
      ColorFromValue(val, pix, alpha);
      return kTRUE;
   } else {
      return kFALSE;
   }
}

#endif
