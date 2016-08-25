##########################################################################
# makefile to compile user ana program. Author V. Vasiliev
# Copy it to your analysis directory and set up paths to
# valid urun.C and uloop.C and cut.C used for slimming.
#
# R10SYS global variable should be set in your login setup to rootana dir
# LD_LIBRARY_PATH should contain $R10SYS/utils/lib to use rootana shared 
# labraries and need N3Db package installed and compiled (see main nemo 
# page software section for details)
##########################################################################

############# USERS SETTINGS #########################################
URUNC       = ${shell pwd}/urun.cpp
ULOOPC	    = ${shell pwd}/uloop.cpp
SLIMC	    = ../GLOBAL_CUTS/eg/eg_cuts.C
######################################################################

GC = gcc34

UTILS       = ${R10SYS}/utils
SRCS        = $(shell find ${UTILS} -name "*\.C"  -print)
LIBRARY     = $(shell find ${UTILS} -name "*\.so"  -print)

OPT         = -O3  
#OPT	    = -Og -ggdb # for optimized debugging
CXX         = $(GC) $(OPT)
INCLUDES    = $(shell ${UTILS}/r10-config --include)
CXXFLAGS    = -fPIC ${INCLUDES} -Wno-deprecated
LD          = $(GC) $(OPT)
LDFLAGS     = 
R10LIBS     = -L${UTILS}/lib -lana10
LIBS        = $(shell ${UTILS}/r10-config --libs) ${R10LIBS}

# command to erase all files
RM          = /bin/rm -f

# list of object files
OBJS        = cut.o urun.o uloop.o

# program executable file name
PROG        = ana


# top level rule, to compile everything
all :  r10lib $(PROG)

r10lib:
	@echo -e "\n Recompiling the Rootana Libraries ... \n"
	gmake -C $(UTILS)/ 

# rule for file "cut.o"
cut.o : $(SLIMC)
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

# rule for file "urun.o"
urun.o : $(URUNC)  
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

# rule for file "urun.o"
uloop.o : $(ULOOPC)  
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

# rule to link the program
$(PROG) : $(OBJS)
	@echo -e "\n <**Linking**> " $@
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROG)

# meta-rule for compiling any "C" source file
%.o : %.C
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

# clean up this mess
clean:
	$(RM) $(PROG) $(OBJS)
	$(RM) core*

