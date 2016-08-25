// @(#)root/proof:$Id: TProofDataSetManagerFile.h 25547 2008-09-25 21:10:17Z rdm $
// Author: Jan Fiete Grosse-Oetringhaus, 08.08.07

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TProofDataSetManagerFile
#define ROOT_TProofDataSetManagerFile

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProofDataSetManagerFile                                             //
//                                                                      //
// Implementation of TProofDataSetManager handling datasets from root   //
// files under a specific directory path                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TProofDataSetManager
#include "TProofDataSetManager.h"
#endif


class TProofDataSetManagerFile : public TProofDataSetManager {

private:
   TString fDataSetDir;        // Location of datasets
   TString fMSSUrl;            // URL for the Mass Storage System
   TString fDataSetLockFile;   // Dataset lock file
   Int_t   fLockFileTimeLimit; // Limit in seconds after a lock automatically expires

protected:
   const char *GetDataSetPath(const char *group, const char *user, const char *dsName);
   Bool_t BrowseDataSets(const char *group, const char *user, UInt_t option, TObject *target);

   Bool_t RemoveDataSet(const char *group, const char *user, const char *dsName);
   Bool_t ExistsDataSet(const char *group, const char *user, const char *dsName);

   Int_t  ScanDataSet(const char *group, const char *user, const char *dsName, UInt_t option);

   void UpdateUsedSpace();

public:
   TProofDataSetManagerFile(const char *group = 0, const char *user = 0, const char *ins = 0);
   virtual ~TProofDataSetManagerFile() { }

   void ParseInitOpts(const char *opts);

   TFileCollection *GetDataSet(const char *uri);
   TMap *GetDataSets(const char *uri, UInt_t /*option*/ = 0);
   Bool_t ExistsDataSet(const char *uri);
   Bool_t RemoveDataSet(const char *uri);

   Int_t RegisterDataSet(const char *uri, TFileCollection *dataSet, const char *opt);
   Int_t ScanDataSet(const char *uri, UInt_t /*option*/ = 0);

   // These should / could be private but they are used directly by the external daemon
   TFileCollection *GetDataSet(const char *group, const char *user, const char *dsName,
                               UInt_t option = 0, TMD5 **checksum = 0);
   TMap *GetDataSets(const char *group, const char *user, UInt_t option = 0);
   Int_t ScanDataSet(TFileCollection *dataset, UInt_t option, Int_t filesmax = -1);
   Int_t WriteDataSet(const char *group, const char *user, const char *dsName,
                      TFileCollection *dataset, UInt_t option = 0, TMD5 *checksum = 0);

   ClassDef(TProofDataSetManagerFile, 0) // DataSet manager for files
};

#endif
