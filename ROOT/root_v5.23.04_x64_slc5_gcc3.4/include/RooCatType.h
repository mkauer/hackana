/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCatType.h,v 1.20 2007/05/11 09:11:30 verkerke Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_CAT_TYPE
#define ROO_CAT_TYPE

#include "Riosfwd.h"
#include "TObject.h"
#include "RooPrintable.h"

class RooCatType : public TObject, public RooPrintable {
public:
  inline RooCatType() : TObject(), RooPrintable() { 
    // Default constructor
    _value = 0 ; _label[0] = 0 ; 
  }
 
  inline RooCatType(const char* name, Int_t value) : TObject(), RooPrintable(), _value(value) { 
    // Constructor with state name and index value
    SetName(name) ; 
  } 
  inline RooCatType(const RooCatType& other) : 
    TObject(other), RooPrintable(other), _value(other._value) { 
    // Copy constructor
    strcpy(_label,other._label) ;
  } ;

  virtual ~RooCatType() {
    // Destructor
  } ;
  virtual TObject* Clone(const char*) const { return new RooCatType(*this); }

  virtual const Text_t* GetName() const { 
    // Return state name
    return _label ; 
  }
  virtual void SetName(const Text_t* name) ;

  inline RooCatType& operator=(const RooCatType& other) { 
    // Assignment operator from other RooCatType
    if (&other==this) return *this ;
    SetName(other.GetName()) ; 
    _value = other._value ; 
    return *this ; } 

  inline Bool_t operator==(const RooCatType& other) {
    // Equality operator with other RooCatType
    return ( _value==other._value && !strcmp(_label,other._label)) ;
  }

  inline Bool_t operator==(Int_t index) { 
    // Return true if index value matches integer
    return (_value==index) ; 
  }

  Bool_t operator==(const char* label) { 
    // Return true if state name matchins string
    return !strcmp(_label,label) ; 
  }

  inline Int_t getVal() const { 
    // Return index value
    return _value ; 
  }
  void setVal(Int_t newValue) { 
    // Set index value
  _value = newValue ; 
  }

  virtual void printName(ostream& os) const ;
  virtual void printTitle(ostream& os) const ;
  virtual void printClassName(ostream& os) const ;
  virtual void printValue(ostream& os) const ;
  
  inline virtual void Print(Option_t *options= 0) const {
    // Printing interface
    printStream(defaultPrintStream(),defaultPrintContents(options),defaultPrintStyle(options));
  }

protected:
  friend class RooAbsCategory ;
  Int_t _value ;     // Index value
  char _label[256] ; // State name
	
  ClassDef(RooCatType,1) // Category state, (name,index) pair
} ;


#endif

