/* @(#)root/thread:$Id: PosixThreadInc.h 20882 2007-11-19 11:31:26Z rdm $ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_PosixThreadInc
#define ROOT_PosixThreadInc


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// PosixThreadInc                                                       //
//                                                                      //
// Common includes for the Posix thread implementation.                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef __CINT__

#ifndef ROOT_TError
#include "TError.h"
#endif
#ifndef ROOT_TSystem
#include "TSystem.h"
#endif

#include <stdlib.h>
#include <time.h>

#if defined(R__LINUX) || defined(R__SOLARIS) || defined(R__ALPHA) || \
    defined(R__MACOSX) || defined(R__HPUX11)
#define PthreadDraftVersion 10
#elif defined(R__HPUX) || defined(R__AIX)
#define PthreadDraftVersion 4
#else
#if !defined(R__KCC)
#warning PthreadDraftVersion not specified
#endif
#endif


#if (PthreadDraftVersion <= 6)
#define ERRNO(x) (((x) != 0) ? (TSystem::GetErrno()) : 0)
#else
#define ERRNO(x) (x)
#endif

#endif /* __CINT__ */

#endif
