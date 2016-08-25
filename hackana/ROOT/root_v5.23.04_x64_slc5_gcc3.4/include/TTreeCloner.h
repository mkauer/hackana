// @(#)root/tree:$Id: TTreeCloner.h 26028 2008-10-30 20:09:44Z brun $
// Author: Philippe Canal 07/11/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeCloner
#define ROOT_TTreeCloner

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTreeCloner                                                          //
//                                                                      //
// Class implementing or helping  the various TTree cloning method      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

#include <vector>

#ifdef R__OLDHPACC
namespace std {
   using ::string;
   using ::vector;
}
#endif

class TBranch;
class TTree;

class TTreeCloner {
   Bool_t     fIsValid;
   TTree     *fFromTree;
   TTree     *fToTree;
   Option_t  *fMethod;
   TObjArray  fFromBranches;
   TObjArray  fToBranches;

   UInt_t     fMaxBaskets;
   UInt_t    *fBasketBranchNum;  //[fMaxBaskets] Index of the branch(es) of the basket.
   UInt_t    *fBasketNum;        //[fMaxBaskets] index of the basket within the branch.

   Long64_t  *fBasketSeek;       //[fMaxBaskets] list of basket position to be read.
   Long64_t  *fBasketEntry;      //[fMaxBaskets] list of basket start entries.
   UInt_t    *fBasketIndex;      //[fMaxBaskets] ordered list of basket indices to be written.

   UShort_t   fPidOffset;        //Offset to be added to the copied key/basket.

   UInt_t     fCloneMethod;      //Indicates which cloning method was selected.
   Long64_t   fToStartEntries;   //Number of entries in the target tree before any addition.

   enum ECloneMethod {
      kDefault             = 0,
      kSortBasketsByBranch = 1,
      kSortBasketsByOffset = 2,
      kSortBasketsByEntry  = 3
   };
   
public:
   TTreeCloner(TTree *from, TTree *to, Option_t *method);
   virtual ~TTreeCloner();

   void   CloseOutWriteBaskets();
   UInt_t CollectBranches(TBranch *from, TBranch *to);
   UInt_t CollectBranches(TObjArray *from, TObjArray *to);
   UInt_t CollectBranches();
   void   CollectBaskets();
   void   CopyMemoryBaskets();
   void   CopyStreamerInfos();
   void   CopyProcessIds();
   Bool_t Exec();
   Bool_t IsValid() { return fIsValid; }
   void   SortBaskets();
   void   WriteBaskets();

   ClassDef(TTreeCloner,0); // helper used for the fast cloning of TTrees.
};

#endif
