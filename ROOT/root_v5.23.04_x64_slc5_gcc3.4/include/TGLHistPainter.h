#ifndef ROOT_TGLHistPainter
#define ROOT_TGLHistPainter

#include <memory>

#ifndef ROOT_TVirtualHistPainter
#include "TVirtualHistPainter.h"
#endif
#ifndef ROOT_TGLPlotPainter
#include "TGLPlotPainter.h"
#endif
#ifndef ROOT_TGLPlotCamera
#include "TGLPlotCamera.h"
#endif
#ifndef ROOT_TGLAdapter
#include "TGLAdapter.h"
#endif

/*
   TGLHistPainter is a proxy class. It inherits TVirtualHistPainter and
   overrides its virtual functions, but all actual work is done by :
      THistPainter - I name it "default" painter, it's the member of type
                     TVirtualHistPainter * and loaded via plugin-manager;
      TGLLegoPainter - it draws different legos (lego/lego1/lego2/lego3);
      TGLSurfacePainter - supports surfaces (surf/surf1/surf2/surf3/surf4/surf5);
      TGLBoxPainter - box option for TH3;
      TGLTF3Painter - TF3.
*/

class TGLParametricEquation;
class TString;
class TList;
class TF3;
class TH1;

class TGLHistPainter : public TVirtualHistPainter {
private:
   //Dynamic type is THistPainter, no problems with simultaneous inheritance and membership
   //TGLHistPainter delegates unsupported options/calls to this object
   std::auto_ptr<TVirtualHistPainter> fDefaultPainter;
   //This member can have different dynamic types: TGLLegoPainter, etc.
   std::auto_ptr<TGLPlotPainter>      fGLPainter;

   TGLParametricEquation *fEq;
   TH1                   *fHist;
   TF3                   *fF3;
   TList                 *fStack;
   EGLPlotType            fPlotType;
   TGLPlotCamera          fCamera;
   TGLPlotCoordinates     fCoord;

   TGLAdapter             fGLDevice;

public:
   TGLHistPainter(TH1 *hist);
   TGLHistPainter(TGLParametricEquation *equation);

   //TVirtualHistPainter final overriders
   Int_t          DistancetoPrimitive(Int_t px, Int_t py);
   void           DrawPanel();
   void           ExecuteEvent(Int_t event, Int_t px, Int_t py);
   TList         *GetContourList(Double_t contour)const;
   char          *GetObjectInfo(Int_t px, Int_t py)const;
   TList         *GetStack()const;
   Bool_t         IsInside(Int_t x, Int_t y);
   Bool_t         IsInside(Double_t x, Double_t y);
   void           Paint(Option_t *option);
   void           PaintStat(Int_t dostat, TF1 *fit);
   void           ProcessMessage(const char *message, const TObject *obj);
   void           SetHistogram(TH1 *hist);
   void           SetStack(TList *stack);
   Int_t          MakeCuts(char *cutsOpt);
   void           SetShowProjection(const char *option, Int_t nbins);

private:

   struct PlotOption_t;

   PlotOption_t   ParsePaintOption(const TString &option)const;
   void           CreatePainter(const PlotOption_t &parsed,
                                const TString &option);

   TGLHistPainter(const TGLHistPainter &);
   TGLHistPainter &operator = (const TGLHistPainter &);

   ClassDef(TGLHistPainter, 0) //Proxy class for GL hist painter
};

#endif
