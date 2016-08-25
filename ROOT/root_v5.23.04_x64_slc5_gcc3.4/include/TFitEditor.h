// @(#)root/fitpanel:$Id: TFitEditor.h 27549 2009-02-20 14:31:41Z moneta $
// Author: Ilka Antcheva, Lorenzo Moneta, David Gonzalez Maline 10/08/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFitEditor
#define ROOT_TFitEditor

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFitEditor                                                           //
//                                                                      //
// Allows to explore and compare various fits.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif

#include "Foption.h"
#include "Math/MinimizerOptions.h"
#include "Fit/DataRange.h"

#include <vector>

//--- Object types
enum EObjectType {
   kObjectHisto,
   kObjectGraph,
   kObjectGraph2D,
   kObjectHStack,
   kObjectTree,
   kObjectMultiGraph
};


class TGTab;
class TVirtualPad;
class TCanvas;
class TGLabel;
class TGComboBox;
class TGTextEntry;
class TGNumberEntry;
class TGDoubleHSlider;
class TGNumberEntry;
class TGNumberEntryField;
class TGStatusBar;
class TAxis;
class TF1;


class TFitEditor : public TGMainFrame {

protected:
   TGTab               *fTab;              // tab widget holding the editor
   TGCompositeFrame    *fTabContainer;     // main tab container
   TGCompositeFrame    *fGeneral;          // general tab
   TGCompositeFrame    *fMinimization;     // minimization tab
   TGTextButton        *fFitButton;        // performs fitting
   TGTextButton        *fResetButton;      // resets fit parameters
   TGTextButton        *fCloseButton;      // close the fit panel
   TGLabel             *fSelLabel;         // contains selected fit function
   TGComboBox          *fDataSet;          // contains list of data set to be fitted
   TGComboBox          *fTypeFit;          // contains the types of functions to be selected
   TGComboBox          *fFuncList;         // contains function list
   TGTextEntry         *fEnteredFunc;      // contains user function file name
   TGTextButton        *fUserButton;       // opens a dialog for user-defined fit method
   TGRadioButton       *fNone;             // set no operation mode
   TGRadioButton       *fAdd;              // set addition mode
   TGRadioButton       *fConv;             // set convolution mode
   TGLayoutHints       *fLayoutNone;       // layout hints of fNone radio button
   TGLayoutHints       *fLayoutAdd;        // layout hints of fAdd radio button
   TGLayoutHints       *fLayoutConv;       // layout hints of fConv radio button
   TGTextButton        *fSetParam;         // open set parameters dialog
   TGCheckButton       *fIntegral;         // switch on/off option 'integral'
   TGCheckButton       *fBestErrors;       // switch on/off option 'improve errors'
   TGCheckButton       *fUseRange;         // switch on/off option 'use function range'
   TGCheckButton       *fAdd2FuncList;     // switch on/off option 'add to list'
   TGCheckButton       *fUseGradient ;     // switch on/off option 'use gradient'
   TGCheckButton       *fAllWeights1;      // switch on/off option 'all weights=1'
   TGCheckButton       *fImproveResults;   // switch on/off option 'improve fit results'
   TGCheckButton       *fEmptyBinsWghts1;  // switch on/off option 'include empry bins'
   TGComboBox          *fMethodList;       // contains method list
   TGCheckButton       *fLinearFit;        // switch on/off linear fit option
   TGCheckButton       *fNoChi2;           // switch on/off option 'No Chi-square'
   TGCheckButton       *fNoStoreDrawing;   // switch on/off 'no store/drwing' option
   TGCheckButton       *fNoDrawing;        // switch on/off 'no drawing' option
   TGCheckButton       *fDrawSame;         // switch on/off fit function drawing
   TGTextButton        *fDrawAdvanced;     // opens a dialog for advanced draw options
   TGDoubleHSlider     *fSliderX;          // slider to set fit range along x-axis
   TGNumberEntry       *fSliderXMax;       // entry to set the maximum in the range
   TGNumberEntry       *fSliderXMin;       // entry to set the minumum in the range
   TGDoubleHSlider     *fSliderY;          // slider to set fit range along y-axis
   TGNumberEntry       *fSliderYMax;       // entry to set the maximum in the range
   TGNumberEntry       *fSliderYMin;       // entry to set the minumum in the range
   TGDoubleHSlider     *fSliderZ;          // slider to set fit range along z-axis
   TGHorizontalFrame   *fSliderXParent;    // parent of fSliderX
   TGHorizontalFrame   *fSliderYParent;    // parent of fSliderY
   TGHorizontalFrame   *fSliderZParent;    // parent of fSliderZ
   TGNumberEntry       *fRobustValue;      // contains robust value for linear fit
   TGRadioButton       *fOptDefault;       // set default printing mode
   TGRadioButton       *fOptVerbose;       // set printing mode to 'Verbose'
   TGRadioButton       *fOptQuiet;         // set printing mode to 'Quiet'
   TVirtualPad         *fParentPad;        // pad containing the object
   TObject             *fFitObject;        // selected object to fit
   EObjectType          fType;             // object type info
   Int_t                fDim;              // object dimension
   TAxis               *fXaxis;            // x-axis
   TAxis               *fYaxis;            // y-axis
   TAxis               *fZaxis;            // z-axis

   // structure holding parameter value and limits
   struct FuncParamData_t { 
      FuncParamData_t() {
         fP[0] = 0; fP[1] = 0; fP[2] = 0; 
      }
      Double_t & operator[](UInt_t i) { return fP[i];}
      Double_t fP[3];
   };
   std::vector<FuncParamData_t>  fFuncPars;         // function parameters (value + limits)

   std::vector<TF1*> fPrevFit;             // Previous succesful fits.

   TGRadioButton       *fLibMinuit;        // set default minimization library (Minuit)
   TGRadioButton       *fLibMinuit2;       // set Minuit2 as minimization library
   TGRadioButton       *fLibFumili;        // set Fumili as minimization library
   TGRadioButton       *fMigrad;           // set default minimization method (MIGRAD)
   TGRadioButton       *fSimplex;          // set Simplex as minimization method
   TGRadioButton       *fFumili;           // set Fumili as minimization method
   TGRadioButton       *fCombination;      // set Combination as minimization method
   TGRadioButton       *fScan;             // set Scan as minimization method
   TGNumberEntryField  *fErrorScale;       // contains error scale set for minimization
   TGNumberEntryField  *fTolerance;        // contains tolerance set for minimization
   TGNumberEntryField  *fIterations;       // contains maximum number of iterations

   TGStatusBar         *fStatusBar;        // statusbar widget
   
   static TFitEditor *fgFitDialog;         // singleton fit panel

   TGComboBox* BuildDataSetList(TGFrame *parent, Int_t id);
   TGComboBox* BuildMethodList(TGFrame *parent, Int_t id);
   void        GetRanges(ROOT::Fit::DataRange&);
   TF1*        GetFitFunction();
   TList*      GetFitObjectListOfFunctions();
   void        DrawSelection(bool restore = false);
   Int_t       CheckFunctionString(const char* str);
   void        CreateFunctionGroup();
   void        CreateGeneralTab();
   void        CreateMinimizationTab();
   void        MakeTitle(TGCompositeFrame *parent, const char *title);
   TF1*        HasFitFunction();
   void        SetEditable(Bool_t);

private:
   TFitEditor(const TFitEditor&);              // not implemented
   TFitEditor& operator=(const TFitEditor&);   // not implemented

   void RetrieveOptions(Foption_t&, TString&, ROOT::Math::MinimizerOptions&, Int_t);

public:
   TFitEditor(TVirtualPad* pad, TObject *obj);
   virtual ~TFitEditor();

//   static TFitEditor *&GetFP();
   static  TFitEditor *GetInstance(TVirtualPad* pad, TObject *obj);
   virtual Option_t  *GetDrawOption() const;
   virtual void       Hide();
   virtual void       Show(TVirtualPad* pad, TObject *obj);

           void       ShowObjectName(TObject* obj);
           Bool_t     SetObjectType(TObject* obj);
   virtual void       Terminate();
           void       UpdateGUI();

   virtual void   CloseWindow();
   virtual void   ConnectSlots();
   virtual void   DisconnectSlots();
   virtual void   RecursiveRemove(TObject* obj);

protected:
   virtual void   SetCanvas(TCanvas *c);

public:
   virtual void   SetFitObject(TVirtualPad *pad, TObject *obj, Int_t event);
   virtual void   SetFunction(const char *function);

   // slot methods 'General' tab
   void           FillFunctionList(Int_t selected = -1);
   virtual void   DoAddition(Bool_t on);
   virtual void   DoAdvancedOptions();
   virtual void   DoAllWeights1();
   virtual void   DoClose();
   virtual void   DoEmptyBinsAllWeights1();
   virtual void   DoEnteredFunction();
   virtual void   DoFit();
   virtual void   DoMaxIterations();
   virtual void   DoDataSet(Int_t sel);
   virtual void   DoFunction(Int_t sel);
   virtual void   DoLinearFit();
   virtual void   DoNoChi2();
   virtual void   DoNoSelection();
   virtual void   DoNoStoreDrawing();
   virtual void   DoReset();
   virtual void   DoSetParameters();
   virtual void   DoSliderXMoved();
   virtual void   DoNumericSliderXChanged();
   virtual void   DoSliderYMoved();
   virtual void   DoNumericSliderYChanged();
   virtual void   DoSliderZMoved();
   virtual void   DoUserDialog();
   virtual void   DoUseFuncRange();

   // slot methods 'Minimization' tab
   virtual void   DoLibrary(Bool_t on);
   virtual void   DoMinMethod(Bool_t on);
   virtual void   DoPrintOpt(Bool_t on);
   
public:
   typedef std::vector<FuncParamData_t > FuncParams_t; 

   
   ClassDef(TFitEditor,0)  //Fit Panel interface
};

#endif
