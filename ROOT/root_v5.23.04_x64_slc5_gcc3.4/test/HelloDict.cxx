//
// File generated by rootcint at Thu Apr 23 11:09:26 2009

// Do NOT change. Changes will be lost next time file is generated
//

#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "HelloDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void Hello_ShowMembers(void *obj, TMemberInspector &R__insp, char *R__parent);
   static void *new_Hello(void *p = 0);
   static void *newArray_Hello(Long_t size, void *p);
   static void delete_Hello(void *p);
   static void deleteArray_Hello(void *p);
   static void destruct_Hello(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Hello*)
   {
      ::Hello *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Hello >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Hello", ::Hello::Class_Version(), "Hello.h", 43,
                  typeid(::Hello), DefineBehavior(ptr, ptr),
                  &::Hello::Dictionary, isa_proxy, 0,
                  sizeof(::Hello) );
      instance.SetNew(&new_Hello);
      instance.SetNewArray(&newArray_Hello);
      instance.SetDelete(&delete_Hello);
      instance.SetDeleteArray(&deleteArray_Hello);
      instance.SetDestructor(&destruct_Hello);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Hello*)
   {
      return GenerateInitInstanceLocal((::Hello*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Hello*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *Hello::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *Hello::Class_Name()
{
   return "Hello";
}

//______________________________________________________________________________
const char *Hello::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hello*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Hello::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hello*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void Hello::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hello*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *Hello::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hello*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void Hello::Streamer(TBuffer &R__b)
{
   // Stream an object of class Hello.

   TTimer::Streamer(R__b);
}

//______________________________________________________________________________
void Hello::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
      // Inspect the data members of an object of class Hello.
      TClass *R__cl = ::Hello::IsA();
      Int_t R__ncp = strlen(R__parent);
      if (R__ncp || R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__parent, "*fList", &fList);
      R__insp.Inspect(R__cl, R__parent, "fI", &fI);
      R__insp.Inspect(R__cl, R__parent, "*fPad", &fPad);
      TTimer::ShowMembers(R__insp, R__parent);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Hello(void *p) {
      return  p ? new(p) ::Hello : new ::Hello;
   }
   static void *newArray_Hello(Long_t nElements, void *p) {
      return p ? new(p) ::Hello[nElements] : new ::Hello[nElements];
   }
   // Wrapper around operator delete
   static void delete_Hello(void *p) {
      delete ((::Hello*)p);
   }
   static void deleteArray_Hello(void *p) {
      delete [] ((::Hello*)p);
   }
   static void destruct_Hello(void *p) {
      typedef ::Hello current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Hello

/********************************************************
* HelloDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && (__GNUC__ > 3) && (__GNUC_MINOR__ > 1)
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableHelloDict();

extern "C" void G__set_cpp_environmentHelloDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("Hello.h");
  G__cpp_reset_tagtableHelloDict();
}
#include <new>
extern "C" int G__cpp_dllrevHelloDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* Hello */
static int G__HelloDict_225_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   Hello* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 1:
     //m: 1
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new Hello((const char*) G__int(libp->para[0]));
     } else {
       p = new((void*) gvp) Hello((const char*) G__int(libp->para[0]));
     }
     break;
   case 0:
     int n = G__getaryconstruct();
     if (n) {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new Hello[n];
       } else {
         p = new((void*) gvp) Hello[n];
       }
     } else {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new Hello;
       } else {
         p = new((void*) gvp) Hello;
       }
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__HelloDictLN_Hello));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((Hello*) G__getstructoffset())->GetWidth());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((Hello*) G__getstructoffset())->GetList());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) Hello::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) Hello::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) Hello::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      Hello::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((Hello*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) Hello::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) Hello::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) Hello::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__HelloDict_225_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) Hello::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef Hello G__THello;
static int G__HelloDict_225_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (Hello*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((Hello*) (soff+(sizeof(Hello)*i)))->~G__THello();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (Hello*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((Hello*) (soff))->~G__THello();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* Hello */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncHelloDict {
 public:
  G__Sizep2memfuncHelloDict(): p(&G__Sizep2memfuncHelloDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncHelloDict::*p)();
};

size_t G__get_sizep2memfuncHelloDict()
{
  G__Sizep2memfuncHelloDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceHelloDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__HelloDictLN_Hello))) {
     Hello *G__Lderived;
     G__Lderived=(Hello*)0x1000;
     {
       TTimer *G__Lpbase=(TTimer*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello),G__get_linked_tagnum(&G__HelloDictLN_TTimer),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TSysEvtHandler *G__Lpbase=(TSysEvtHandler*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello),G__get_linked_tagnum(&G__HelloDictLN_TSysEvtHandler),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello),G__get_linked_tagnum(&G__HelloDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
     {
       TQObject *G__Lpbase=(TQObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello),G__get_linked_tagnum(&G__HelloDictLN_TQObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableHelloDict() {

   /* Setting up typedef entry */
   G__search_typename2("Float_t",102,-1,0,-1);
   G__setnewtype(-1,"Float 4 bytes (float)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<TSchemaHelper>",117,G__get_linked_tagnum(&G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__HelloDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__HelloDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* Hello */
static void G__setup_memvarHello(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello));
   { Hello *p; p=(Hello*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__HelloDictLN_TList),-1,-1,4,"fList=",0,"list of characters");
   G__memvar_setup((void*)0,104,0,0,-1,G__defined_typename("UInt_t"),-1,4,"fI=",0,"\"infinit\"  counter");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__HelloDictLN_TPad),-1,-1,4,"fPad=",0,"pad where this text is drawn");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__HelloDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarHelloDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncHello(void) {
   /* Hello */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello));
   G__memfunc_setup("Hello",500,G__HelloDict_225_0_1, 105, G__get_linked_tagnum(&G__HelloDictLN_Hello), -1, 0, 1, 1, 1, 0, "C - - 10 '\"Hello, World!\"' text", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Notify",633,(G__InterfaceMethod) NULL,103, -1, G__defined_typename("Bool_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ExecuteEvent",1237,(G__InterfaceMethod) NULL,121, -1, -1, 0, 3, 1, 1, 0, 
"i - 'Int_t' 0 - event i - 'Int_t' 0 - px "
"i - 'Int_t' 0 - py", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("DistancetoPrimitive",1991,(G__InterfaceMethod) NULL,105, -1, G__defined_typename("Int_t"), 0, 2, 1, 1, 0, 
"i - 'Int_t' 0 - - i - 'Int_t' 0 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("GetWidth",800,G__HelloDict_225_0_5, 102, -1, G__defined_typename("Float_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Paint",508,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "C - 'Option_t' 10 '\"\"' option", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Print",525,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 8, "C - 'Option_t' 10 '\"\"' -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ls",223,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 8, "C - 'Option_t' 10 '\"\"' -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("GetList",700,G__HelloDict_225_0_9, 85, G__get_linked_tagnum(&G__HelloDictLN_TList), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__HelloDict_225_0_10, 85, G__get_linked_tagnum(&G__HelloDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&Hello::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__HelloDict_225_0_11, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&Hello::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__HelloDict_225_0_12, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&Hello::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__HelloDict_225_0_13, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&Hello::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__HelloDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 2, 1, 1, 0, 
"u 'TMemberInspector' - 1 - insp C - - 0 - parent", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__HelloDict_225_0_17, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__HelloDict_225_0_18, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&Hello::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__HelloDict_225_0_19, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&Hello::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__HelloDict_225_0_20, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&Hello::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__HelloDict_225_0_21, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&Hello::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~Hello", 626, G__HelloDict_225_0_22, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncHelloDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {
}

static void G__cpp_setup_global3() {
}

static void G__cpp_setup_global4() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalHelloDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
  G__cpp_setup_global4();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcHelloDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__HelloDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TList = { "TList" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TTimer = { "TTimer" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TQObject = { "TQObject" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TSysEvtHandler = { "TSysEvtHandler" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_TPad = { "TPad" , 99 , -1 };
G__linked_taginfo G__HelloDictLN_Hello = { "Hello" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableHelloDict() {
  G__HelloDictLN_TClass.tagnum = -1 ;
  G__HelloDictLN_TBuffer.tagnum = -1 ;
  G__HelloDictLN_TMemberInspector.tagnum = -1 ;
  G__HelloDictLN_TObject.tagnum = -1 ;
  G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__HelloDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__HelloDictLN_TList.tagnum = -1 ;
  G__HelloDictLN_TTimer.tagnum = -1 ;
  G__HelloDictLN_TQObject.tagnum = -1 ;
  G__HelloDictLN_TSysEvtHandler.tagnum = -1 ;
  G__HelloDictLN_TPad.tagnum = -1 ;
  G__HelloDictLN_Hello.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableHelloDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TList);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TTimer);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TQObject);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TSysEvtHandler);
   G__get_linked_tagnum_fwd(&G__HelloDictLN_TPad);
   G__tagtable_setup(G__get_linked_tagnum(&G__HelloDictLN_Hello),sizeof(Hello),-1,62720,"animated text with cool wave effect",G__setup_memvarHello,G__setup_memfuncHello);
}
extern "C" void G__cpp_setupHelloDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupHelloDict()");
  G__set_cpp_environmentHelloDict();
  G__cpp_setup_tagtableHelloDict();

  G__cpp_setup_inheritanceHelloDict();

  G__cpp_setup_typetableHelloDict();

  G__cpp_setup_memvarHelloDict();

  G__cpp_setup_memfuncHelloDict();
  G__cpp_setup_globalHelloDict();
  G__cpp_setup_funcHelloDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncHelloDict();
  return;
}
class G__cpp_setup_initHelloDict {
  public:
    G__cpp_setup_initHelloDict() { G__add_setup_func("HelloDict",(G__incsetup)(&G__cpp_setupHelloDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initHelloDict() { G__remove_setup_func("HelloDict"); }
};
G__cpp_setup_initHelloDict G__cpp_setup_initializerHelloDict;

