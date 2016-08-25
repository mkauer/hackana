#
# Author: M. Kordosky
#
# a very simple makefile, which builds a shared library <package_name>.so
# from all the .C files in the package directory
#
# The only dependencies are on ROOT
#
# The shared library is deposited in the local (e.g. build) directory
# in principle it is better to move all libraries into one lib/ directory
# in a "make install" step but this makefile doesn't do that
#
#

GC = gcc34

package_name  = rootana_utils

OPT           = -O3
#OPT	      = -Og -ggdb # for optimized debugging
INCLUDES      = $(shell ${R10SYS}/utils/r10-config --include)
CXX           = $(GC) $(OPT)
CXXFLAGS      = -fPIC $(INCLUDES) -Wno-deprecated
LD            = $(GC) $(OPT)
LDFLAGS       = -shared
LIBS          = $(shell ${R10SYS}/utils/r10-config --libs)


# make an object file <filename>.o out of each source file <filename>.c
objects = $(shell ls *\.C | sed 's/\.C/\.o/g')

lib = lib$(package_name).so

lib2 = lib$(package_name).a

all : $(lib)
	gmake -C slim2
	gmake -C ana10

$(lib) : $(objects)

echo:
	@echo -e "\n objects = $(objects) "
	@echo -e "\n library = $(lib) "

mo100 : $(lib) mo100.o
	@echo -e "\n <**Linking**> " $<
	$(LD) mo100.o $(LIBS) $(R10LIBS) -o mo100

mo100.o : mo100.cpp
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@


%.o : %.C
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.so :
	@echo -e "\n <**Linking**> " $@
	$(LD) $(LDFLAGS) $(objects) -o $@
	@mkdir -p lib
	@mv $@ $(R10SYS)/utils/lib

%.a :
	@echo -e "\n <**Linking**> " $@
	ar cru $@ $(objects) 
	@mkdir -p lib
	@mv $@ $(R10SYS)/utils/lib


# remove the object files and the library
clean:
	@/bin/rm -rf *.o
	@/bin/rm -rf lib/*.so
	gmake clean -C slim2
	gmake clean -C ana10

