#ifndef ROOT_TGLBoxPainter
#define ROOT_TGLBoxPainter

#ifndef ROOT_TGLPlotPainter
#include "TGLPlotPainter.h"
#endif
#ifndef ROOT_TGLQuadric
#include "TGLQuadric.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif
#ifndef ROOT_TGLUtil
#include "TGLUtil.h"
#endif

class TGLPlotCamera;
class TAxis;
class TH1;

class TGLBoxPainter : public TGLPlotPainter {
private:
   TGLTH3Slice fXOZSlice;
   TGLTH3Slice fYOZSlice;
   TGLTH3Slice fXOYSlice;

   enum EBoxType {
      kBox, //boxes, sizes are proportional to bin content
      kBox1 //spheres, not boxes
   };

   mutable EBoxType        fType;

   TString                 fPlotInfo;
   Rgl::Range_t            fMinMaxVal;

   mutable TGLQuadric      fQuadric;

   TGLBoxPainter(const TGLBoxPainter &);
   TGLBoxPainter &operator = (const TGLBoxPainter &);

public:
   TGLBoxPainter(TH1 *hist, TGLPlotCamera *camera, TGLPlotCoordinates *coord, TGLPaintDevice *dev = 0);

   char   *GetPlotInfo(Int_t px, Int_t py);
   Bool_t  InitGeometry();
   void    StartPan(Int_t px, Int_t py);
   void    Pan(Int_t px, Int_t py);
   void    AddOption(const TString &stringOption);
   void    ProcessEvent(Int_t event, Int_t px, Int_t py);

private:
   //Overriders
   void    InitGL()const;
   void    DrawPlot()const;

   void    SetPlotColor()const;

   void    DrawSectionXOZ()const;
   void    DrawSectionYOZ()const;
   void    DrawSectionXOY()const;

   void    DrawPalette()const;
   void    DrawPaletteAxis()const;

   Bool_t  HasSections()const;

   ClassDef(TGLBoxPainter, 0)//Box painter
};

#endif
