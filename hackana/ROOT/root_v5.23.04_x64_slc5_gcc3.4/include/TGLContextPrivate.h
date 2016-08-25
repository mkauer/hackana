// @(#)root/gl:$Id: TGLContextPrivate.h 24162 2008-06-06 11:33:13Z matevz $
// Author:  Timur Pocheptsov, Matevz Tadel, June 2007

#ifndef ROOT_TGLContextPrivate
#define ROOT_TGLContextPrivate

#include <map>

#ifndef ROOT_TGLIncludes
#include "TGLIncludes.h"
#endif
#ifndef ROOT_TGLContext
#include "TGLContext.h"
#endif

#ifdef WIN32

class TGLContextPrivate {
public:
   HWND        fHWND;
   HDC         fHDC;
   HGLRC       fGLContext;

   TGLContextPrivate()
      : fHWND(0),
        fHDC(0),
        fGLContext(0)
   {
   }
   static void RegisterContext(TGLContext *ctx);
   static void RemoveContext(TGLContext *ctx);
   static TGLContext *GetCurrentContext();


private:
   TGLContextPrivate(const TGLContextPrivate &);
   TGLContextPrivate &operator = (const TGLContextPrivate &);

   static std::map<HGLRC, TGLContext *> fgContexts;
};


#else

class TGLContextPrivate {
public:
   Display     *fDpy;
   XVisualInfo *fVisualInfo;
   GLXContext   fGLContext;
   Window       fWindowID;
   //GLXPbuffer   fPBDC;

   TGLContextPrivate()
      : fDpy(0),
        fVisualInfo(0),
        fGLContext(0),
        fWindowID(0)
   {
   }

   static void RegisterContext(TGLContext *ctx);
   static void RemoveContext(TGLContext *ctx);
   static TGLContext *GetCurrentContext();

private:
   TGLContextPrivate(const TGLContextPrivate &);
   TGLContextPrivate &operator = (const TGLContextPrivate &);

   static std::map<GLXContext, TGLContext *> fgContexts;
};

#endif

#endif
