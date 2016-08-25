/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooCintUtils.h 28259 2009-04-16 16:21:16Z wouter $
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
#ifndef ROO_CINT_UTILS
#define ROO_CINT_UTILS

#include "Rtypes.h"
#include <list>
#include <string>
namespace RooCintUtils {
  
  std::pair<std::list<std::string>,unsigned int> ctorArgs(const char* classname, UInt_t nMinArgs=0) ;
  Bool_t isEnum(const char* typeName) ;
  Bool_t isValidEnumValue(const char* typeName, const char* value) ;
  const char* functionName(void* func) ;
  Bool_t matchFuncPtrArgs(void* func, const char* args) ;
  
};

#endif
