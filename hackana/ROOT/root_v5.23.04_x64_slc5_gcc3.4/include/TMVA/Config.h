// @(#)root/tmva $Id: Config.h 26050 2008-11-01 09:18:41Z brun $   
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Config                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      GLobal configuration settings (singleton)                                 *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker    <Andreas.Hocker@cern.ch>     - CERN, Switzerland       *
 *      Joerg Stelzer      <Joerg.Stelzer@cern.ch>      - CERN, Switzerland       *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *      Helge Voss         <Helge.Voss@cern.ch>         - MPI-K Heidelberg, GER   *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U., USA                                                        *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

#ifndef ROOT_TMVA_Config
#define ROOT_TMVA_Config

#include "Rtypes.h"
#include "TString.h"

#ifndef ROOT_TMVA_MsgLogger
#include "TMVA/MsgLogger.h"
#endif

namespace TMVA {

   class Config {
               
   public:

      static Config& Instance() { return fgConfigPtr ? *fgConfigPtr : *(fgConfigPtr = new Config()); }
      virtual ~Config();

      Bool_t UseColor() { return fUseColoredConsole; }
      void   SetUseColor( Bool_t uc ) { fUseColoredConsole = uc; }

      Bool_t IsSilent() { return fSilent; }
      void   SetSilent( Bool_t s ) { fSilent = s; }

      Bool_t WriteOptionsReference() { return fWriteOptionsReference; }
      void   SetWriteOptionsReference( Bool_t w ) { fWriteOptionsReference = w; }

   public:

      class VariablePlotting;
      class IONames;

      VariablePlotting& GetVariablePlotting() { return fVariablePlotting; }
      IONames&          GetIONames()          { return fIONames; }

      // publicly accessible global settings
      class VariablePlotting {
         // data collection class to configure plotting of variables
      public:
         Float_t fTimesRMS;
         Int_t   fNbins1D;
         Int_t   fNbins2D;
         Int_t   fMaxNumOfAllowedVariablesForScatterPlots;
         Int_t   fNbinsXOfROCCurve;
      } fVariablePlotting; // Customisable plotting properties

      // for file names and similar
      class IONames {
      public:
         TString fWeightFileDir;
         TString fWeightFileExtension;
         TString fOptionsReferenceFileDir;
      } fIONames; // Customisable weight file properties
         
      
   private:

      // private constructor
      Config();
      static Config* fgConfigPtr;
                  
   private:

      Bool_t fUseColoredConsole;     // coloured standard output
      Bool_t fSilent;                // no output at all
      Bool_t fWriteOptionsReference; // if set true: Configurable objects write file with option reference

      mutable MsgLogger fLogger;   // message logger
         
      ClassDef(Config,0) // Singleton class for global configuration settings
   };

   // global accessor
   Config& gConfig();
}

#endif
