# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch
-include ../MyConfig.mk

#------------------------------------------------------------------------------

EVENTO        = Event.$(ObjSuf) EventDict.$(ObjSuf)
EVENTS        = Event.$(SrcSuf) EventDict.$(SrcSuf)
EVENTSO       = libEvent.$(DllSuf)
EVENT         = Event$(ExeSuf)
ifeq ($(PLATFORM),win32)
EVENTLIB      = libEvent.lib
else
EVENTLIB      = $(shell pwd)/$(EVENTSO)
endif

MAINEVENTO    = MainEvent.$(ObjSuf)
MAINEVENTS    = MainEvent.$(SrcSuf)

HWORLDO       = hworld.$(ObjSuf)
HWORLDS       = hworld.$(SrcSuf)
HWORLD        = hworld$(ExeSuf)

CTORTUREO     = ctorture.$(ObjSuf)
CTORTURES     = ctorture.$(SrcSuf)
CTORTURE      = ctorture$(ExeSuf)

HSIMPLEO      = hsimple.$(ObjSuf)
HSIMPLES      = hsimple.$(SrcSuf)
HSIMPLE       = hsimple$(ExeSuf)

MINEXAMO      = minexam.$(ObjSuf)
MINEXAMS      = minexam.$(SrcSuf)
MINEXAM       = minexam$(ExeSuf)

TSTRINGO      = tstring.$(ObjSuf)
TSTRINGS      = tstring.$(SrcSuf)
TSTRING       = tstring$(ExeSuf)

TCOLLEXO      = tcollex.$(ObjSuf)
TCOLLEXS      = tcollex.$(SrcSuf)
TCOLLEX       = tcollex$(ExeSuf)

TCOLLBMO      = tcollbm.$(ObjSuf)
TCOLLBMS      = tcollbm.$(SrcSuf)
TCOLLBM       = tcollbm$(ExeSuf)

VVECTORO      = vvector.$(ObjSuf)
VVECTORS      = vvector.$(SrcSuf)
VVECTOR       = vvector$(ExeSuf)

VMATRIXO      = vmatrix.$(ObjSuf)
VMATRIXS      = vmatrix.$(SrcSuf)
VMATRIX       = vmatrix$(ExeSuf)

STRESSLO      = stressLinear.$(ObjSuf)
STRESSLS      = stressLinear.$(SrcSuf)
STRESSL       = stressLinear$(ExeSuf)

STRESSGO      = stressGraphics.$(ObjSuf)
STRESSGS      = stressGraphics.$(SrcSuf)
STRESSG       = stressGraphics$(ExeSuf)

ifeq ($(shell $(ROOTCONFIG) --has-genvector),yes)
STRESSVECO   = stressVector.$(ObjSuf)
STRESSVECS   = stressVector.$(SrcSuf)
STRESSVEC    = stressVector$(ExeSuf)

STRESSMATHO   = stressMathCore.$(ObjSuf) 
STRESSMATHS   = stressMathCore.$(SrcSuf) 
ifneq ($(USE_REFLEX),)
CXXFLAGS += -DUSE_REFLEX
ifeq ($(PLATFORM),win32)
STRESSMATHLIBS = '$(ROOTSYS)/lib/libGenVector.lib' '$(ROOTSYS)/lib/libReflex.lib' '$(ROOTSYS)/lib/libCintex.lib' 
else
STRESSMATHLIBS = -lGenVector -lReflex -lCintex
endif

TRACKMATHSRC = TrackMathCoreRflx.$(SrcSuf)
TRACKMATHOBJ = TrackMathCoreRflx.$(ObjSuf)
TRACKMATHLIB = libTrackMathCoreRflx.$(DllSuf)

else 

ifeq ($(PLATFORM),win32)
STRESSMATHLIBS = '$(ROOTSYS)/lib/libGenVector.lib'
else
STRESSMATHLIBS = -lGenVector
endif

TRACKMATHSRC = TrackMathCoreDict.$(SrcSuf)
TRACKMATHOBJ = TrackMathCoreDict.$(ObjSuf)
TRACKMATHLIB = libTrackMathCoreDict.$(DllSuf)
 
endif

STRESSMATH    = stressMathCore$(ExeSuf)

endif

ifeq ($(shell $(ROOTCONFIG) --has-mathmore),yes)
STRESSMATHMOREO   = stressMathMore.$(ObjSuf) 
STRESSMATHMORES   = stressMathMore.$(SrcSuf)  
ifeq ($(PLATFORM),win32)
STRESSMATHMORELIBS = '$(ROOTSYS)/lib/libMathMore.lib'  
else
STRESSMATHMORELIBS = -lMathMore
endif
STRESSMATHMORE    = stressMathMore$(ExeSuf)
endif

VLAZYO        = vlazy.$(ObjSuf)
VLAZYS        = vlazy.$(SrcSuf)
VLAZY         = vlazy$(ExeSuf)

GUITESTO      = guitest.$(ObjSuf)
GUITESTS      = guitest.$(SrcSuf)
GUITEST       = guitest$(ExeSuf)

GUIVIEWERO    = guiviewer.$(ObjSuf) guiviewerDict.$(ObjSuf)
GUIVIEWERS    = guiviewer.$(SrcSuf) guiviewerDict.$(SrcSuf)
GUIVIEWER     = guiviewer$(ExeSuf)

HELLOO        = Hello.$(ObjSuf) HelloDict.$(ObjSuf)
HELLOS        = Hello.$(SrcSuf) HelloDict.$(SrcSuf)
HELLOSO       = Hello.$(DllSuf)

ACLOCKO       = Aclock.$(ObjSuf) AclockDict.$(ObjSuf)
ACLOCKS       = Aclock.$(SrcSuf) AclockDict.$(SrcSuf)
ACLOCKSO      = Aclock.$(DllSuf)

TETRISO       = Tetris.$(ObjSuf) TetrisDict.$(ObjSuf)
TETRISS       = Tetris.$(SrcSuf) TetrisDict.$(SrcSuf)
TETRISSO      = Tetris.$(DllSuf)

STRESSO       = stress.$(ObjSuf)
STRESSS       = stress.$(SrcSuf)
STRESS        = stress$(ExeSuf)

STRESSGEOMETRYO   = stressGeometry.$(ObjSuf)
STRESSGEOMETRYS   = stressGeometry.$(SrcSuf)
STRESSGEOMETRY    = stressGeometry$(ExeSuf)

STRESSSHAPESO   = stressShapes.$(ObjSuf)
STRESSSHAPESS   = stressShapes.$(SrcSuf)
STRESSSHAPES    = stressShapes$(ExeSuf)

ifeq ($(shell $(ROOTCONFIG) --has-roofit),yes)
STRESSROOFITO  = stressRooFit.$(ObjSuf)
STRESSROOFITS  = stressRooFit.$(SrcSuf)
STRESSROOFIT   = stressRooFit$(ExeSuf)
endif

STRESSFITO    = stressFit.$(ObjSuf)
STRESSFITS    = stressFit.$(SrcSuf)
STRESSFIT     = stressFit$(ExeSuf)

STRESSTMVAO    = stressTMVA.$(ObjSuf)
STRESSTMVAS    = stressTMVA.$(SrcSuf)
STRESSTMVA     = stressTMVA$(ExeSuf)

STRESSENTRYLISTO = stressEntryList.$(ObjSuf)
STRESSENTRYLISTS = stressEntryList.$(SrcSuf)
STRESSENTRYLIST  = stressEntryList$(ExeSuf)

STRESSHEPIXO  = stressHepix.$(ObjSuf)
STRESSHEPIXS  = stressHepix.$(SrcSuf)
STRESSHEPIX   = stressHepix$(ExeSuf)

STRESSSPO     = stressSpectrum.$(ObjSuf)
STRESSSPS     = stressSpectrum.$(SrcSuf)
STRESSSP      = stressSpectrum$(ExeSuf)

STRESSPROOFO  = stressProof.$(ObjSuf)
STRESSPROOFS  = stressProof.$(SrcSuf)
STRESSPROOF   = stressProof$(ExeSuf)

TBENCHO       = TBench.$(ObjSuf) TBenchDict.$(ObjSuf)
TBENCHS       = TBench.$(SrcSuf) TBenchDict.$(SrcSuf)
TBENCHSO      = TBench.$(DllSuf)

BENCHO        = bench.$(ObjSuf)
BENCHS        = bench.$(SrcSuf)
BENCH         = bench$(ExeSuf)

TESTBITSO     = testbits.$(ObjSuf)
TESTBITSS     = testbits.$(SrcSuf)
TESTBITS      = testbits$(ExeSuf)

QPRANDOMO     = QpRandomDriver.$(ObjSuf)
QPRANDOMS     = QpRandomDriver.$(SrcSuf)
QPRANDOM      = QpRandomDriver$(ExeSuf)

THREADSO      = threads.$(ObjSuf)
THREADSS      = threads.$(SrcSuf)
THREADS       = threads$(ExeSuf)

STRESSITERO    = stressIterators.$(ObjSuf)
STRESSITERS    = stressIterators.$(SrcSuf)
STRESSITER     = stressIterators$(ExeSuf)

STRESSHISTO   = stressHistogram.$(ObjSuf)
STRESSHISTS   = stressHistogram.$(SrcSuf)
STRESSHIST    = stressHistogram$(ExeSuf)


OBJS          = $(EVENTO) $(MAINEVENTO) $(HWORLDO) $(HSIMPLEO) $(MINEXAMO) \
                $(TSTRINGO) $(TCOLLEXO) $(VVECTORO) $(VMATRIXO) $(VLAZYO) \
                $(HELLOO) $(ACLOCKO) $(STRESSO) $(TBENCHO) $(BENCHO) \
                $(STRESSSHAPESO) $(TCOLLBMO) $(STRESSGEOMETRYO) $(STRESSLO) \
                $(STRESSGO) $(STRESSSPO) $(TESTBITSO) $(STRESSTMVAO) \
                $(CTORTUREO) $(QPRANDOMO) $(THREADSO) $(STRESSVECO) \
                $(STRESSMATHO) $(STRESSFITO) $(STRESSHEPIXO) \
                $(STRESSENTRYLISTO) $(STRESSROOFITO) $(STRESSPROOFO) \
                $(STRESSMATHMOREO) $(STRESSITERO) $(STRESSHISTO)

PROGRAMS      = $(EVENT) $(HWORLD) $(HSIMPLE) $(MINEXAM) $(TSTRING) \
                $(TCOLLEX) $(TCOLLBM) $(VVECTOR) $(VMATRIX) $(VLAZY) \
                $(HELLOSO) $(ACLOCKSO) $(STRESS) $(TBENCHSO) $(BENCH) \
                $(STRESSSHAPES) $(STRESSGEOMETRY) $(STRESSL) $(STRESSG) \
                $(TESTBITS) $(CTORTURE) $(QPRANDOM) $(THREADS) $(STRESSSP) \
                $(STRESSVEC) $(STRESSFIT) $(STRESSHEPIX) $(STRESSENTRYLIST) \
                $(STRESSROOFIT) $(STRESSPROOF) $(STRESSMATH) $(STRESSTMVA) \
                $(STRESSMATHMORE) $(STRESSITER) $(STRESSHIST)


OBJS         += $(GUITESTO) $(GUIVIEWERO) $(TETRISO)
PROGRAMS     += $(GUITEST) $(GUIVIEWER) $(TETRISSO)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    Aclock Hello Tetris

all:            $(PROGRAMS)

$(EVENTSO):     $(EVENTO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
		bindexplib $* $^ > $*.def
		lib -nologo -MACHINE:IX86 $^ -def:$*.def \
		   $(OutPutOpt)$(EVENTLIB)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) \
		   $(OutPutOpt)$@
		$(MT_DLL)
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
endif
endif
endif
endif
		@echo "$@ done"

$(EVENT):       $(EVENTSO) $(MAINEVENTO)
		$(LD) $(LDFLAGS) $(MAINEVENTO) $(EVENTLIB) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(HWORLD):      $(HWORLDO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(CTORTURE):    $(CTORTUREO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(HSIMPLE):     $(HSIMPLEO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(MINEXAM):     $(MINEXAMO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(TSTRING):     $(TSTRINGO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(TCOLLEX):     $(TCOLLEXO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(TCOLLBM):     $(TCOLLBMO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(VVECTOR):     $(VVECTORO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(VMATRIX):     $(VMATRIXO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(VLAZY):       $(VLAZYO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSL):     $(STRESSLO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSG):     $(STRESSGO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSSP):    $(STRESSSPO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libSpectrum.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lSpectrum $(OutPutOpt)$@
endif
		@echo "$@ done"

$(STRESSVEC):   $(STRESSVECO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libGenVector.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lGenVector $(OutPutOpt)$@
endif
		@echo "$@ done"






$(TRACKMATHLIB):	$(TRACKMATHOBJ)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(LIBS) $(STRESSMATHLIBS) $(OutPutOpt)$@
		$(MT_DLL)
#endif
endif
endif
endif

$(STRESSMATH):   $(STRESSMATHO) $(TRACKMATHLIB)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $(STRESSMATHO) $(LIBS) $(STRESSMATHLIBS)  $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $(STRESSMATHO) $(LIBS) $(STRESSMATHLIBS)  $(OutPutOpt)$@
endif
		@echo "$@ done"


$(STRESSMATHMORE):   $(STRESSMATHMOREO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(STRESSMATHMORELIBS)  $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) $(STRESSMATHMORELIBS)  $(OutPutOpt)$@
endif
		@echo "$@ done"


$(TESTBITS):    $(TESTBITSO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(THREADS):     $(THREADSO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libThread.lib' $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"
else
ifeq ($(HASTHREAD),yes)
		$(LD) $(LDFLAGS) $^ $(LIBS) -lThread $(OutPutOpt)$@
		@echo "$@ done"
else
		@echo "This version of ROOT has no thread support, $@ not built"
endif
endif

$(QPRANDOM):    $(QPRANDOMO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libQuadp.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lQuadp $(OutPutOpt)$@
endif
		@echo "$@ done"

$(GUITEST):     $(GUITESTO)
		$(LD) $(LDFLAGS) $^ $(GLIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(GUIVIEWER):   $(GUIVIEWERO)
		$(LD) $(LDFLAGS) $^ $(GLIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESS):      $(STRESSO) $(EVENT)
		$(LD) $(LDFLAGS) $(STRESSO) $(EVENTLIB) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSGEOMETRY):  $(STRESSGEOMETRYO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libGeom.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lGeom $(OutPutOpt)$@
endif
		@echo "$@ done"

$(STRESSSHAPES):  $(STRESSSHAPESO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libGeom.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lGeom $(OutPutOpt)$@
endif
		@echo "$@ done"

$(STRESSFIT):   $(STRESSFITO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSTMVA): $(STRESSTMVAO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libTMVA.lib'   '$(ROOTSYS)/lib/libMinuit.lib' '$(ROOTSYS)/lib/libMLP.lib' '$(ROOTSYS)/lib/libTreePlayer.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lTMVA -lMinuit -lMLP -lTreePlayer $(OutPutOpt)$@
endif
		@echo "$@ done"

$(STRESSENTRYLIST):	$(STRESSENTRYLISTO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(STRESSHEPIX): $(STRESSHEPIXO) $(STRESSGEOMETRY) $(STRESSFIT) $(STRESSL) \
                $(STRESSSP) $(STRESS)
		$(LD) $(LDFLAGS) $(STRESSHEPIXO) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSROOFITO): stressRooFit_tests.cxx

$(STRESSROOFIT): $(STRESSROOFITO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libRooFit.lib' '$(ROOTSYS)/lib/libRooFitCore.lib' '$(ROOTSYS)/lib/libHtml.lib' '$(ROOTSYS)/lib/libThread.lib' '$(ROOTSYS)/lib/libMinuit.lib' '$(ROOTSYS)/lib/libFoam.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $^ $(LIBS) -lRooFit -lRooFitCore -lHtml -lThread -lMinuit -lFoam $(OutPutOpt)$@
endif
		@echo "$@ done"

$(STRESSPROOF): $(STRESSPROOFO)
ifeq ($(PLATFORM),win32)
		$(LD) $(LDFLAGS) $^ $(LIBS) '$(ROOTSYS)/lib/libProof.lib' '$(ROOTSYS)/lib/libThread.lib' $(OutPutOpt)$@
		$(MT_EXE)
else
		$(LD) $(LDFLAGS) $(STRESSPROOFO) $(LIBS) -lProof -lThread $(OutPutOpt)$@
endif
		@echo "$@ done"

$(BENCH):       $(BENCHO) $(TBENCHSO)
		$(LD) $(LDFLAGS) $(BENCHO) $(TBENCHO) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

Hello:          $(HELLOSO)
$(HELLOSO):     $(HELLOO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
endif
endif
endif

Aclock:         $(ACLOCKSO)
$(ACLOCKSO):    $(ACLOCKO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
endif
endif
endif

Tetris:         $(TETRISSO)
$(TETRISSO):    $(TETRISO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(GLIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(GLIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	 	   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
endif
endif
endif

$(TBENCHSO):    $(TBENCHO)
ifeq ($(ARCH),aix)
		/usr/ibmcxx/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(ARCH),aix5)
		/usr/vacpp/bin/makeC++SharedLib $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS) $^ -lstdc++ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
endif
endif
endif

$(STRESSITER):	$(STRESSITERO) 
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(STRESSHIST):  $(STRESSHISTO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"


clean:
		@rm -f $(OBJS) $(TRACKMATHSRC) core

distclean:      clean
		-@mv -f linearIO.root linearIO.roott
		@rm -f $(PROGRAMS) $(EVENTSO) $(EVENTLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations
		@rm -rf cxx_repository files
		-@mv -f linearIO.roott linearIO.root
		-@cd RootShower && $(MAKE) distclean
		-@cd rhtml && $(MAKE) distclean
		-@cd RootIDE && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)

###
stressIterators.$(ObjSuf): stressIterators.h
 
Event.$(ObjSuf): Event.h
MainEvent.$(ObjSuf): Event.h

EventDict.$(SrcSuf): Event.h EventLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^

Hello.$(ObjSuf): Hello.h
HelloDict.$(SrcSuf): Hello.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^

Aclock.$(ObjSuf): Aclock.h
AclockDict.$(SrcSuf): Aclock.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^

Tetris.$(ObjSuf): Tetris.h
TetrisDict.$(SrcSuf): Tetris.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^

TBench.$(ObjSuf): TBench.h
TBenchDict.$(SrcSuf): TBench.h benchLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^

guiviewer.$(ObjSuf): guiviewer.h
guiviewerDict.$(SrcSuf): guiviewer.h guiviewerLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTCINT) -f $@ -c $^

stressMathCore.$(ObjSuf): 	TrackMathCore.h
TrackMathCoreDict.$(SrcSuf): 	TrackMathCore.h TrackMathCoreLinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^	

TrackMathCoreRflx.$(SrcSuf): 	TrackMathCore.h TrackMathCoreRflx.xml
			@echo "Generating dictionary $@ using gccxml ..."
			genreflex TrackMathCore.h --selection_file=TrackMathCoreRflx.xml -o TrackMathCoreRflx.cxx  -I$(ROOTSYS)/include

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<

ifeq ($(GCC_MAJOR),4)
ifeq ($(GCC_MINOR),1)
TBenchDict.o: CXXFLAGS += -Wno-strict-aliasing
endif
endif
