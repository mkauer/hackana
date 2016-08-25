// @(#)root/tmva $Id: DecisionTree.h 26050 2008-11-01 09:18:41Z brun $
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : DecisionTree                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation of a Decision Tree                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Xavier Prudent  <prudent@lapp.in2p3.fr>  - LAPP, France                   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 *                                                                                *
 **********************************************************************************/

#ifndef ROOT_TMVA_DecisionTree
#define ROOT_TMVA_DecisionTree

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DecisionTree                                                         //
//                                                                      //
// Implementation of a Decision Tree                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TH2.h"
#include "TRandom.h"
#include "TRandom2.h"

#ifndef ROOT_TMVA_DecisionTreeNode
#include "TMVA/DecisionTreeNode.h"
#endif
#ifndef ROOT_TMVA_BinarySearchTree
#include "TMVA/BinaryTree.h"
#endif
#ifndef ROOT_TMVA_BinarySearchTree
#include "TMVA/BinarySearchTree.h"
#endif
#ifndef ROOT_TMVA_SeparationBase
#include "TMVA/SeparationBase.h"
#endif

using std::vector;

namespace TMVA {
   
   class Event;
   
   class DecisionTree : public BinaryTree {
      
   public:
      
      // the constructur needed for the "reading" of the decision tree from weight files
      DecisionTree( void );

      // the constructur for a decsion tree from a root node
      DecisionTree( DecisionTreeNode *n );

      // the constructur needed for constructing the decision tree via training with events
      DecisionTree( SeparationBase *sepType,Int_t minSize, 
                    Int_t nCuts, SeparationBase *qtype=NULL,
                    Bool_t randomisedTree=kFALSE, Int_t useNvars=0, Int_t iSeed=0);

      // copy constructor
      DecisionTree (const DecisionTree &d);

      virtual ~DecisionTree( void );

      virtual Node * CreateNode() { return new DecisionTreeNode(); }
  
      // building of a tree by recursivly splitting the nodes 
      Int_t BuildTree( vector<TMVA::Event*> & eventSample, 
                       DecisionTreeNode *node = NULL );

      // determine the way how a node is split (which variable, which cut value)
      Double_t TrainNode( vector<TMVA::Event*> & eventSample,  DecisionTreeNode *node ) { return TrainNodeFast( eventSample, node ); }
      Double_t TrainNodeFast( vector<TMVA::Event*> & eventSample,  DecisionTreeNode *node );
      Double_t TrainNodeFull( vector<TMVA::Event*> & eventSample,  DecisionTreeNode *node );

      // fill at tree with a given structure already (just see how many signa/bkgr
      // events end up in each node 
      void FillTree( vector<TMVA::Event*> & eventSample);

      // fill the existing the decision tree structure by filling event
      // in from the top node and see where they happen to end up
      void FillEvent( TMVA::Event & event,  
                      TMVA::DecisionTreeNode *node  );

      // returns: 1 = Signal (right),  -1 = Bkg (left)
      Double_t CheckEvent( const TMVA::Event & , Bool_t UseYesNoLeaf = kFALSE ); 

      // return the individual relative variable importance 
      vector< Double_t > GetVariableImportance();
      Double_t GetVariableImportance(Int_t ivar);

      // clear the tree nodes (their S/N, Nevents etc), just keep the structure of the tree
      void ClearTree();

      // set pruning method
      enum EPruneMethod { kExpectedErrorPruning=0, kCostComplexityPruning, kNoPruning };
      void SetPruneMethod( EPruneMethod m = kCostComplexityPruning ) { fPruneMethod = m; }

      // recursive pruning of the tree
      void PruneTree();

      void SetPruneStrength( Double_t p ) { fPruneStrength = p; }

      void SetNodePurityLimit( Double_t p ) { fNodePurityLimit = p; }
      Double_t GetNodePurityLimit( ) const { return fNodePurityLimit; }


      void DescendTree( DecisionTreeNode *n = NULL );
      void SetParentTreeInNodes( DecisionTreeNode *n = NULL );

      DecisionTreeNode* GetLeftDaughter( DecisionTreeNode *n );

      DecisionTreeNode* GetRightDaughter( DecisionTreeNode *n );

      // retrieve node from the tree. Its position (up to a maximal tree depth of 64)
      // is coded as a sequence of left-right moves starting from the root, coded as
      // 0-1 bit patterns stored in the "long-integer" together with the depth
      DecisionTreeNode* GetNode( ULong_t sequence, UInt_t depth );

      void CleanTree(DecisionTreeNode *node=NULL);


      void PruneTreeEEP(DecisionTreeNode *node);
      void PruneTreeCC();

      void PruneNode(TMVA::DecisionTreeNode *node);

      UInt_t CountLeafNodes(TMVA::DecisionTreeNode *n = NULL);

   private:

      Double_t GetNodeError(DecisionTreeNode *node);

      Double_t GetSubTreeError(DecisionTreeNode *node);

      // utility functions

      // calculates the min and max values for each variable in event sample
      void FindMinAndMax(vector<TMVA::Event*> & eventSample,
                         vector<Double_t> & min,
                         vector<Double_t> & max);
      
      // set up the grid for the cut scan
      void SetCutPoints(vector<Double_t> & cut_points,
                        Double_t xmin,
                        Double_t xmax,
                        Int_t num_gridpoints);
            
      // calculate the Purity out of the number of sig and bkg events collected
      // from individual samples.
      
      // calculates the purity S/(S+B) of a given event sample
      Double_t SamplePurity(vector<Event*> eventSample);
      
      Int_t     fNvars;          // number of variables used to separate S and B
      Int_t     fNCuts;          // number of grid point in variable cut scans
      SeparationBase *fSepType;  // the separation crition
      
      Double_t  fMinSize;        // min number of events in node
      Double_t  fMinSepGain;     // min number of separation gain to perform node splitting
      
      Bool_t    fUseSearchTree;  // cut scan done with binary trees or simple event loop.
      Double_t  fPruneStrength;  // a parameter to set the "amount" of pruning..needs to be adjusted 
      
      EPruneMethod fPruneMethod; // method used for prunig 

      Double_t  fNodePurityLimit;// purity limit to decide whether a node is signal

      Bool_t    fRandomisedTree; // choose at each node splitting a random set of variables 
      Int_t     fUseNvars;       // the number of variables used in randomised trees;
      
      TRandom2  *fMyTrandom;     // random number generator for randomised trees

      std::vector< Double_t > fVariableImportance; // the relative importance of the different variables 
      
      SeparationBase *fQualityIndex;  // separation/quality criterio for CC-pruning

      static const Int_t  fgDebugLevel = 0;     // debug level determining some printout/control plots etc.

      ClassDef(DecisionTree,0)                  // implementation of a Decision Tree
   };

} // namespace TMVA

#endif 

