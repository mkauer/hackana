#############################################################
#  makefile that can be copied anywhere with the slim2.C
#  and sets up correctly as long as the system variables
#  have been set up correctly!
#
#  Matt Kauer
#  26 Sept 2008
#############################################################

#*********  USERS SETTINGS  ******************************
#CUT        = ${shell pwd}/eg_cuts.C
#SLIM       = ${shell pwd}/slim2.C
#CUT	   = /unix/nemo2/users/kauer/GLOBAL_CUTS/eg/eg_cuts.C
CUT	   = /unix/nemo2/users/kauer/GLOBAL_CUTS/eg/foil_eg_cuts.C
SLIM       = /unix/nemo2/users/kauer/GLOBAL_CUTS/eg/slim2.C
#*********************************************************

UTILS      =  $(R10SYS)/utils
SRCS       =  ${shell find $(UTILS) -name "*.C"  -print}
LIBRARY    =  ${shell find $(UTILS) -name "*.so"  -print}

OPT        = -g 
CXX        = g++ -Wno-deprecated
INCLUDES   = $(shell ${R10SYS}/utils/r10-config --include)
CXXFLAGS   = $(OPT) -fPIC $(INCLUDES) $(shell $(ROOTSYS)/bin/root-config --cflags) -m32
LD         = g++
LDFLAGS    = -m32

ROOTLIBS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
R10LIBS    = $(shell $(R10SYS)/utils/r10-config --libs)
LIBS       = $(R10LIBS) $(ROOTLIBS) -lMinuit


# command to erase all files
RM         = /bin/rm -f

# list of object files
OBJS       = cut.o slim2.o

# program executable file name
PROG       = slim2

# top level rule, to compile everything
all :  $(LIBRARY) $(PROG)

$(LIBRARY): $(SRCS)
	@echo " Recompile rootana libraries"
	gmake -C $(UTILS)/ 

# rule for file "cut.o"
cut.o : $(CUT)
	@echo " <**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

# rule for file "slim2.o"
slim2.o : $(SLIM)
	@echo " <**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o $@

# rule to link the program
$(PROG) : $(OBJS)
	@echo " <**Linking**> "	
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROG)

# meta-rule for compiling any "C" source file
%.o : %.C
	@echo " <**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean that mess up
clean:
	$(RM) $(PROG) $(OBJS) 
	$(RM) core*

