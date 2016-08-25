// @(#)root/gl:$Id: TGLViewer.h 28197 2009-04-14 13:59:27Z matevz $
// Author:  Richard Maunder  25/05/2005

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLViewer
#define ROOT_TGLViewer

#include "TGLViewerBase.h"
#include "TGLRnrCtx.h"
#include "TGLSelectRecord.h"

#include "TVirtualViewer3D.h"
#include "TBuffer3D.h"

#include "TGLPerspectiveCamera.h"
#include "TGLOrthoCamera.h"

#include "TTimer.h"
#include "TPoint.h"

#include "TGEventHandler.h"

#include "GuiTypes.h"
#include "TQObject.h"

#include <vector>

class TGLSceneBase;
class TGLRedrawTimer;
class TGLViewerEditor;
class TGLWidget;
class TGLLightSet;
class TGLClipSet;
class TGLManipSet;
class TGLCameraOverlay;
class TGLContextIdentity;
class TTimer;

class TContextMenu;


class TGLViewer : public TVirtualViewer3D,
                  public TGLViewerBase,
                  public TQObject

{
   friend class TGLOutput;
   friend class TGLEventHandler;
public:

   enum ECameraType { kCameraPerspXOZ,  kCameraPerspYOZ,  kCameraPerspXOY,
                      kCameraOrthoXOY,  kCameraOrthoXOZ,  kCameraOrthoZOY,
                      kCameraOrthoXnOY, kCameraOrthoXnOZ, kCameraOrthoZnOY };

private:
   TGLViewer(const TGLViewer &);             // Not implemented
   TGLViewer & operator=(const TGLViewer &); // Not implemented

   void InitSecondaryObjects();

protected:
   // External handles
   TVirtualPad   *fPad;         //! external pad - remove replace with signal

   // GUI Handles
   TContextMenu  *fContextMenu; //!

   // Cameras
   // TODO: Put in vector and allow external creation
   TGLPerspectiveCamera fPerspectiveCameraXOZ; //!
   TGLPerspectiveCamera fPerspectiveCameraYOZ; //!
   TGLPerspectiveCamera fPerspectiveCameraXOY; //!
   TGLOrthoCamera       fOrthoXOYCamera;       //!
   TGLOrthoCamera       fOrthoXOZCamera;       //!
   TGLOrthoCamera       fOrthoZOYCamera;       //!
   TGLOrthoCamera       fOrthoXnOYCamera;       //!
   TGLOrthoCamera       fOrthoXnOZCamera;       //!
   TGLOrthoCamera       fOrthoZnOYCamera;       //!
   TGLCamera           *fCurrentCamera;        //!

   // Lights
   TGLLightSet         *fLightSet;             //!
   // Clipping
   TGLClipSet          *fClipSet;              //!
   // Selected physical
   TGLSelectRecord      fCurrentSelRec;        //! select record in use as selected
   TGLSelectRecord      fSelRec;               //! select record from last select (should go to context)
   TGLSelectRecord      fSecSelRec;            //! select record from last secondary select (should go to context)
   TGLManipSet         *fSelectedPShapeRef;    //!
   // Overlay
   TGLOverlayElement   *fCurrentOvlElm;        //! current overlay element
   TGLOvlSelectRecord   fOvlSelRec;            //! select record from last overlay select

   TGEventHandler      *fEventHandler;         //! event handler

   // Mouse ineraction
public:
   enum EPushAction   { kPushStd,
                        kPushCamCenter, kPushAnnotate };
   enum EDragAction   { kDragNone,
                        kDragCameraRotate, kDragCameraTruck, kDragCameraDolly,
                        kDragOverlay };
protected:
   EPushAction          fPushAction;
   EDragAction          fDragAction;

   // Redraw timer
   TGLRedrawTimer      *fRedrawTimer;        //! timer for triggering redraws
   Float_t              fMaxSceneDrawTimeHQ; //! max time for scene rendering at high LOD (in ms)
   Float_t              fMaxSceneDrawTimeLQ; //! max time for scene rendering at high LOD (in ms)

   TGLRect        fViewport;       //! viewport - drawn area
   TGLColorSet    fDarkColorSet;   //! color-set with dark background
   TGLColorSet    fLightColorSet;  //! color-set with dark background
   Int_t          fAxesType;       //! axes type
   Bool_t         fAxesDepthTest;  //! remove guides hidden-lines
   Bool_t         fReferenceOn;    //! reference marker on?
   TGLVertex3     fReferencePos;   //! reference position
   Bool_t         fDrawCameraCenter; //! reference marker on?
   TGLCameraOverlay  *fCameraOverlay; //! markup size of viewport in scene units

   Bool_t         fInitGL;         //! has GL been initialised?
   Bool_t         fSmartRefresh;   //! cache logicals during scene rebuilds, use TAtt3D time-stamp to determine if they are still valid

   // Debug tracing (for scene rebuilds)
   Bool_t         fDebugMode;            //! debug mode (forced rebuild + draw scene/frustum/interest boxes)
   Bool_t         fIsPrinting;           //!
   TString        fPictureFileName;      //! default file-name for SavePicture()
   Float_t        fFader;                //! fade the view (0 - no fade/default, 1 - full fade/no rendering done)

   static TGLColorSet fgDefaultColorSet;                 //! a shared, default color-set
   static Bool_t      fgUseDefaultColorSetForNewViewers; //! name says it all


   ///////////////////////////////////////////////////////////////////////
   // Methods
   ///////////////////////////////////////////////////////////////////////
   // Drawing - can tidy up/remove lots when TGLManager added
   void InitGL();
   void PreDraw();
   void PostDraw();
   void FadeView(Float_t alpha);
   void MakeCurrent() const;
   void SwapBuffers() const;

   // Cameras
   void        SetViewport(Int_t x, Int_t y, Int_t width, Int_t height);
   void        SetupCameras(Bool_t reset);

protected:
   TGLWidget          *fGLWidget;
   Int_t               fGLDevice; //!for embedded gl viewer
   TGLContextIdentity *fGLCtxId;  //!for embedded gl viewer

   // Updata/camera-reset behaviour
   Bool_t           fIgnoreSizesOnUpdate;      // ignore sizes of bounding-boxes on update
   Bool_t           fResetCamerasOnUpdate;     // reposition camera on each update
   Bool_t           fResetCamerasOnNextUpdate; // reposition camera on next update
   Bool_t           fResetCameraOnDoubleClick; // reposition camera on double-click

public:
   TGLViewer(TVirtualPad* pad, Int_t x, Int_t y, Int_t width, Int_t height);
   TGLViewer(TVirtualPad* pad);
   virtual ~TGLViewer();

   void ResetInitGL();

   // TVirtualViewer3D interface ... mostly a facade

   // Forward to TGLScenePad
   virtual Bool_t CanLoopOnPrimitives() const { return kTRUE; }
   virtual void   PadPaint(TVirtualPad* pad);
   // Actually used by GL-in-pad
   virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
   virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
   // Only implemented because they're abstract ... should throw an
   // exception or assert they are not called.
   virtual Bool_t PreferLocalFrame() const { return kTRUE; }
   virtual void   BeginScene() {}
   virtual Bool_t BuildingScene() const { return kFALSE; }
   virtual void   EndScene() {}
   virtual Int_t  AddObject(const TBuffer3D&, Bool_t* = 0) { return TBuffer3D::kNone; }
   virtual Int_t  AddObject(UInt_t, const TBuffer3D&, Bool_t* = 0) { return TBuffer3D::kNone; }
   virtual Bool_t OpenComposite(const TBuffer3D&, Bool_t* = 0) { return kFALSE; }
   virtual void   CloseComposite() {}
   virtual void   AddCompositeOp(UInt_t) {}

   virtual void   PrintObjects();
   virtual void   ResetCameras()                { SetupCameras(kTRUE); }
   virtual void   ResetCamerasAfterNextUpdate() { fResetCamerasOnNextUpdate = kTRUE; }

   virtual void   RefreshPadEditor(TObject* = 0) {}

   TGLWidget* GetGLWidget() { return fGLWidget; }

   virtual void  CreateGLWidget()  {}
   virtual void  DestroyGLWidget() {}

   Int_t   GetDev()          const           { return fGLDevice; }
   Bool_t  GetSmartRefresh() const           { return fSmartRefresh; }
   void    SetSmartRefresh(Bool_t smart_ref) { fSmartRefresh = smart_ref; }

   TGLColorSet& RefDarkColorSet()  { return fDarkColorSet;  }
   TGLColorSet& RefLightColorSet() { return fLightColorSet; }
   TGLColorSet& ColorSet()         { return * fRnrCtx->GetBaseColorSet(); }
   void         UseDarkColorSet();
   void         UseLightColorSet();
   void         SwitchColorSet();

   void         UseDefaultColorSet(Bool_t x);
   Bool_t       IsUsingDefaultColorSet() const;

   void         SetClearColor(Color_t col);

   static TGLColorSet& GetDefaultColorSet();
   static void         UseDefaultColorSetForNewViewers(Bool_t x);
   static Bool_t       IsUsingDefaultColorSetForNewViewers();

   TGLLightSet* GetLightSet() const { return fLightSet; }
   TGLClipSet * GetClipSet()  const { return fClipSet; }

   // External GUI component interface
   TGLCamera & CurrentCamera() const { return *fCurrentCamera; }
   TGLCamera & RefCamera(ECameraType camera);
   void SetCurrentCamera(ECameraType camera);
   void SetOrthoCamera(ECameraType camera, Double_t zoom, Double_t dolly,
                             Double_t center[3], Double_t hRotate, Double_t vRotate);
   void SetPerspectiveCamera(ECameraType camera, Double_t fov, Double_t dolly,
                             Double_t center[3], Double_t hRotate, Double_t vRotate);
   void GetGuideState(Int_t & axesType, Bool_t & axesDepthTest, Bool_t & referenceOn, Double_t* referencePos) const;
   void SetGuideState(Int_t axesType, Bool_t axesDepthTest, Bool_t referenceOn, const Double_t* referencePos);
   void SetDrawCameraCenter(Bool_t x);
   Bool_t GetDrawCameraCenter() { return fDrawCameraCenter; }
   void   PickCameraCenter()    { fPushAction = kPushCamCenter; RefreshPadEditor(this); }
   void   PickAnnotate()        { fPushAction = kPushAnnotate;  RefreshPadEditor(this); }
   TGLCameraOverlay* GetCameraOverlay() const { return fCameraOverlay; }
   void SetCameraOverlay(TGLCameraOverlay* m) { fCameraOverlay = m; }

   EPushAction GetPushAction() const { return fPushAction; }
   EDragAction GetDragAction() const { return fDragAction; }

   const TGLPhysicalShape * GetSelected() const;


   // Draw and selection

   // Scene rendering timeouts
   Float_t GetMaxSceneDrawTimeHQ() const    { return fMaxSceneDrawTimeHQ; }
   Float_t GetMaxSceneDrawTimeLQ() const    { return fMaxSceneDrawTimeLQ; }
   void    SetMaxSceneDrawTimeHQ(Float_t t) { fMaxSceneDrawTimeHQ = t; }
   void    SetMaxSceneDrawTimeLQ(Float_t t) { fMaxSceneDrawTimeLQ = t; }

   // Request methods post cross thread request via TROOT::ProcessLineFast().
   void RequestDraw(Short_t LOD = TGLRnrCtx::kLODMed); // Cross thread draw request
   virtual void PreRender();
   void DoDraw();

   void DrawGuides();
   void DrawDebugInfo();

   Bool_t RequestSelect(Int_t x, Int_t y, Bool_t trySecSel=kFALSE); // Cross thread select request
   Bool_t DoSelect(Int_t x, Int_t y, Bool_t trySecSel=kFALSE);      // Window coords origin top left
   void   ApplySelection();

   Bool_t RequestOverlaySelect(Int_t x, Int_t y); // Cross thread select request
   Bool_t DoOverlaySelect(Int_t x, Int_t y);      // Window coords origin top left

   // Saveing of screen image
   Bool_t       SavePicture(const TString &fileName);
   Bool_t       SavePicture();
   const char*  GetPictureFileName() const { return fPictureFileName.Data(); }
   void         SetPictureFileName(const TString& f) { fPictureFileName = f; }
   Float_t      GetFader() const { return fFader; }
   void         SetFader(Float_t x) { fFader = x; }
   void         AutoFade(Float_t fade, Float_t time=1, Int_t steps=10);

   // Update/camera-reset
   void   UpdateScene();
   Bool_t GetIgnoreSizesOnUpdate() const        { return fIgnoreSizesOnUpdate; }
   void   SetIgnoreSizesOnUpdate(Bool_t v)      { fIgnoreSizesOnUpdate = v; }
   void   ResetCurrentCamera();
   Bool_t GetResetCamerasOnUpdate() const       { return fResetCamerasOnUpdate; }
   void   SetResetCamerasOnUpdate(Bool_t v)     { fResetCamerasOnUpdate = v; }
   Bool_t GetResetCameraOnDoubleClick() const   { return fResetCameraOnDoubleClick; }
   void   SetResetCameraOnDoubleClick(Bool_t v) { fResetCameraOnDoubleClick = v; }

   virtual void PostSceneBuildSetup(Bool_t resetCameras);

   virtual void SelectionChanged();    // *SIGNAL*
   virtual void OverlayDragFinished(); // *SIGNAL*
   virtual void MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t); // *SIGNAL*
   virtual void MouseOver(TGLPhysicalShape*); // *SIGNAL*
   virtual void MouseOver(TGLPhysicalShape*, UInt_t state); // *SIGNAL*
   virtual void Activated() { Emit("Activated()"); } // *SIGNAL*
   virtual void Clicked(TObject *obj); //*SIGNAL*
   virtual void Clicked(TObject *obj, UInt_t button, UInt_t state); //*SIGNAL*
   virtual void DoubleClicked() { Emit("DoubleClicked()"); } // *SIGNAL*

   TGEventHandler *GetEventHandler() const { return fEventHandler; }
   virtual void    SetEventHandler(TGEventHandler *handler);

   virtual void RemoveOverlayElement(TGLOverlayElement* el);

   TGLSelectRecord&    GetSelRec()    { return fSelRec; }
   TGLOvlSelectRecord& GetOvlSelRec() { return fOvlSelRec; }
   TGLOverlayElement*  GetCurrentOvlElm() const { return fCurrentOvlElm; }
   void                ClearCurrentOvlElm();

   ClassDef(TGLViewer,0) // Standard ROOT GL viewer.
};



// TODO: Find a better place/way to do this
class TGLRedrawTimer : public TTimer
{
private:
   TGLViewer & fViewer;
   Short_t     fRedrawLOD;
   Bool_t      fPending;
public:
   TGLRedrawTimer(TGLViewer & viewer) :
      fViewer(viewer), fRedrawLOD(TGLRnrCtx::kLODHigh), fPending(kFALSE) {}
   ~TGLRedrawTimer() {}
   void RequestDraw(Int_t milliSec, Short_t redrawLOD)
   {
      if (fPending) TurnOff(); else fPending = kTRUE;
      if (redrawLOD < fRedrawLOD) fRedrawLOD = redrawLOD;
      TTimer::Start(milliSec, kTRUE);
   }
   Bool_t IsPending() const { return fPending; }
   virtual void Stop()
   {
      if (fPending) { TurnOff(); fPending = kFALSE; }
   }
   Bool_t Notify()
   {
      TurnOff();
      fPending = kFALSE;
      fViewer.RequestDraw(fRedrawLOD);
      fRedrawLOD = TGLRnrCtx::kLODHigh;
      return kTRUE;
   }
};

class TGLFaderHelper {
private:
   TGLFaderHelper(const TGLFaderHelper&); // Not implemented
   TGLFaderHelper& operator=(const TGLFaderHelper&); // Not implemented

public:
   TGLViewer *fViewer;
   Float_t    fFadeTarget;
   Float_t    fTime;
   Int_t      fNSteps;

   TGLFaderHelper() :
      fViewer(0), fFadeTarget(0), fTime(0), fNSteps(0) {}
   TGLFaderHelper(TGLViewer* v, Float_t fade, Float_t time, Int_t steps) :
      fViewer(v),fFadeTarget(fade), fTime(time), fNSteps(steps) {}
   virtual ~TGLFaderHelper() {}

   void MakeFadeStep();

   ClassDef(TGLFaderHelper, 0);
};

#endif // ROOT_TGLViewer
