// @(#)root/sessionviewer:$Id: TProofProgressDialog.h 25072 2008-08-06 09:26:41Z ganis $
// Author: Fons Rademakers   21/03/03

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProofProgressDialog
#define ROOT_TProofProgressDialog


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProofProgressDialog                                                 //
//                                                                      //
// This class provides a query progress bar.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TTime
#include "TTime.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

class TGTransientFrame;
class TGProgressBar;
class TGTextButton;
class TGCheckButton;
class TGLabel;
class TGTextBuffer;
class TGTextEntry;
class TProof;
class TProofProgressLog;
class TProofProgressMemoryPlot;
class TNtuple;
class TGraph;


class TProofProgressDialog {

   friend class TProofProgressLog;
   friend class TProofProgressMemoryPlot;

private:
   enum EQueryStatus { kRunning = 0, kDone, kStopped, kAborted, kIncomplete };

   TGTransientFrame   *fDialog;  // transient frame, main dialog window
   TGProgressBar      *fBar;     // progress bar
   TGTextButton       *fClose;
   TGTextButton       *fStop;
   TGTextButton       *fAbort;
   TGTextButton       *fLog;
   TGTextButton       *fRatePlot;
   TGTextButton       *fMemPlot;
   TGCheckButton      *fKeepToggle;
   TGCheckButton      *fLogQueryToggle;
   TGTextBuffer       *fTextQuery;
   TGTextEntry        *fEntry;
   TGLabel            *fTitleLab;
   TGLabel            *fFilesEvents;
   TGLabel            *fProcessed;
   TGLabel            *fTotal;
   TGLabel            *fRate;
   TGLabel            *fInit;
   TGLabel            *fSelector;
   TProofProgressLog  *fLogWindow;       // transient frame for logs
   TProofProgressMemoryPlot *fMemWindow;  // transient frame for memory plots
   TProof             *fProof;
   TTime               fStartTime;
   TTime               fEndTime;
   Long64_t            fPrevProcessed;
   Long64_t            fPrevTotal;
   Long64_t            fFirst;
   Long64_t            fEntries;
   Int_t               fFiles;
   EQueryStatus        fStatus;
   Bool_t              fKeep;
   Bool_t              fLogQuery;
   TNtuple            *fRatePoints;
   TGraph             *fRateGraph;
   Float_t             fProcTime;
   Double_t            fAvgRate;
   Double_t            fAvgMBRate;
   Int_t               fSVNRev;

   TString             fSessionUrl;

   static Bool_t       fgKeepDefault;
   static Bool_t       fgLogQueryDefault;
   static TString      fgTextQueryDefault;

public:
   TProofProgressDialog(TProof *proof, const char *selector,
                        Int_t files, Long64_t first, Long64_t entries);
   virtual ~TProofProgressDialog();

   void ResetProgressDialog(const char *sel, Int_t sz, Long64_t fst, Long64_t ent);
   void Progress(Long64_t total, Long64_t processed);
   void Progress(Long64_t total, Long64_t processed, Long64_t bytesread,
                 Float_t initTime, Float_t procTime,
                 Float_t evtrti, Float_t mbrti);
   void IndicateStop(Bool_t aborted);
   void LogMessage(const char *msg, Bool_t all);

   void CloseWindow();
   void DoClose();
   void DoLog();
   void DoKeep(Bool_t on);
   void DoSetLogQuery(Bool_t on);
   void DoStop();
   void DoAbort();
   void DoPlotRateGraph();
   void DoMemoryPlot();

   ClassDef(TProofProgressDialog,0)  //PROOF progress dialog
};

#endif
