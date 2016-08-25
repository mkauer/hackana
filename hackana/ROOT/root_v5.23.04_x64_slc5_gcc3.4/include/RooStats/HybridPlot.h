// @(#)root/roostats:$Id: HybridPlot.h 26434 2008-11-24 21:29:32Z moneta $

/*************************************************************************
 * Project: RooStats                                                     *
 * Package: RooFit/RooStats                                              *
 * Authors:                                                              *
 *   Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke       *
 *************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_HybridPlot
#define ROOSTATS_HybridPlot

#include <vector>
#include <iostream>

#ifndef ROOT_TNamed 
#include "TNamed.h"
#endif

// these  should be maybe forward decleared 
// by moving implementations in source file 
#include "TH1.h"
#include "TCanvas.h"


class TLine; 
class TLegend; 


namespace RooStats {

   class HybridPlot : public TNamed {

   public:

      /// Constructor
      HybridPlot(const char* name,
                 const char* title,
                 std::vector<double> sb_values,
                 std::vector<double> b_values,
                 double m2lnQ_data,
                 int n_bins,
                 bool verbosity=true);

      /// Destructor
      ~HybridPlot();

      /// Draw on canvas
      void Draw (const char* options="");

      /// All the objects are written to rootfile
      void DumpToFile (const char* RootFileName, const char* options);

      /// Get B histo mean
      double GetBmean(){return fB_histo->GetMean();};

      /// Get B histo RMS
      double GetBrms(){return fB_histo->GetRMS();};

      /// Get B histo
      TH1F * GetBhisto(){return fB_histo;}

      /// Get B histo center
      double GetBCenter(double n_sigmas=1, bool display=false)
      {return GetHistoCenter(fB_histo,n_sigmas,display);};

      /// Get B histo integration extremes to obtain the requested area fraction
      double* GetBIntExtremes(double frac)
      {return GetHistoPvals(fB_histo,frac);};

      /// Get SB histo mean
      double GetSBmean(){return fSb_histo->GetMean();};

      /// Get SB histo center
      double GetSBCenter(double n_sigmas=1, bool display=false)
      {return GetHistoCenter(fSb_histo,n_sigmas,display);};

      /// Get SB histo RMS
      double GetSBrms(){return fSb_histo->GetRMS();};

      /// Get SB histo integration extremes to obtain the requested area fraction
      double* GetSBIntExtremes(double frac)
      {return GetHistoPvals(fSb_histo,frac);};

      /// Get B histo
      TH1F* GetSBhisto(){return fSb_histo;}

      /// from Statistical plot

      /// Get the canvas
      TCanvas* GetCanvas(){return fCanvas;}

      /// Set the canvas
      void SetCanvas(TCanvas* new_canvas){fCanvas=new_canvas;}

      /// Write an image on disk
      void DumpToImage (const char* filename){fCanvas->Print(filename);}

      // moved from Rsc.h

      /// Get the center of the histo
      double GetHistoCenter(TH1* histo, double n_rms=1,bool display_result=false);

      /// Get the "effective sigmas" of the histo
      double* GetHistoPvals (TH1* histo, double percentage);

      /// Get the median of an histogram
      double GetMedian(TH1* histo);

   private:

      TH1F* fSb_histo; // The sb Histo
      TH1F* fSb_histo_shaded; // The sb Histo shaded
      TH1F* fB_histo; // The b Histo
      TH1F* fB_histo_shaded; // The b Histo shaded
      TLine* fData_m2lnQ_line; // The line for the data -2lnQ
      TLegend* fLegend; // The legend of the plot
      bool fVerbose; // verbosity flag
      TCanvas* fCanvas; // plot canvas

      ClassDef(HybridPlot,1)   // Provides the plots for an HybridResult
   };
}

#endif
