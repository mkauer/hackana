// @(#)root/tmva $Id: BinaryTree.h 21630 2008-01-10 19:40:44Z brun $    
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : BinaryTree                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      BinaryTree: A base class for BinarySearch- or Decision-Trees              *
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
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_BinaryTree
#define ROOT_TMVA_BinaryTree

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BinaryTree                                                           //
//                                                                      //
// Base class for BinarySearch and Decision Trees                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"

#ifndef ROOT_TMVA_Node
#include "TMVA/Node.h"
#endif
#ifndef ROOT_TMVA_Event
#include "TMVA/Event.h"
#endif
#ifndef ROOT_TMVA_MsgLogger
#include "TMVA/MsgLogger.h"
#endif

// -----------------------------------------------------------------------------

// the actual tree class
// Handles allocation, deallocation, and sorting of nodes.
// the Tree consists of a "root-node" wich might have  0 to 2 daughther nodes

namespace TMVA {
   
   class BinaryTree;
   ostream& operator<< ( ostream& os, const BinaryTree& tree );
   istream& operator>> ( istream& istr,     BinaryTree& tree );
   
   class BinaryTree {
      
      friend ostream& operator<< ( ostream& os, const BinaryTree& tree );
      friend istream& operator>> ( istream& istr,     BinaryTree& tree );
      
   public:
      
      // or a tree with Root node "n", any daughters of this node are automatically in the tree
      BinaryTree( void );

      virtual ~BinaryTree();

      virtual Node* CreateNode() = 0;

      // set the root node of the tree
      void SetRoot( Node* r ) { fRoot = r; }
    
      // Retrieves the address of the root node
      Node* GetRoot() const { return fRoot; }
    
      // get number of Nodes in the Tree as counted while booking the nodes;
      UInt_t GetNNodes() const { return fNNodes; }

      // count the number of Nodes in the Tree by looping through the tree and updates
      // the stored number. (e.g. useful when pruning, as the number count is updated when
      // building the tree.
      UInt_t CountNodes( Node* n = NULL );

      UInt_t GetTotalTreeDepth() const { return fDepth; }

      void SetTotalTreeDepth( Int_t depth ) { fDepth = depth;};

      void SetTotalTreeDepth( Node* n = NULL );

      Node* GetLeftDaughter ( Node* n);    
      Node* GetRightDaughter( Node* n);

      void Print( ostream& os ) const;
      void Read ( istream& istr );

   private:
  
      Node*      fRoot;                //the root node of the tree
      // the tree only has it's root node, the "daughters" are taken car 
      // of by the "node" properties of the "root"
  
   protected:

      // delete a node (and the corresponding event if owned by the tree)
      void       DeleteNode( Node* );

      Int_t      fNNodes;           // total number of nodes in the tree (counted)
      UInt_t     fDepth;            // maximal depth in tree reached

      mutable MsgLogger  fLogger;   // message loggera    

      ClassDef(BinaryTree,0) // Base class for BinarySearch and Decision Trees
   };  

} // namespace TMVA

#endif

