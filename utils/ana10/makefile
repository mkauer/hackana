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

package_name  = ana10

OPT           = -O3
#OPT	      = -Og -ggdb # for optimized debugging
CXX           = $(GC) $(OPT)
INCLUDES      = $(shell ${R10SYS}/utils/r10-config --include)
CXXFLAGS      = -fPIC $(INCLUDES) -Wno-deprecated
LD            = $(GC) $(OPT)
LDFLAGS       = -shared
LIBS          = $(shell $(R10SYS)/utils/r10-config --libs)

# make an object file <filename>.o out of each source file <filename>.c
objects = $(shell ls *\.C | sed 's/\.C/\.o/g') 

lib = lib$(package_name).so


all : $(lib)


echo:
	@echo -e "\n objects = $(objects) "
	@echo -e "\n library = $(lib) "

$(lib) : $(objects)


%.o : %.C
	@echo -e "\n <**Compiling**> " $<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.so :
	@echo -e "\n <**Linking**> " $@
	$(LD) $(LDFLAGS) $(objects) -o $@
	@mv $@ $(R10SYS)/utils/lib

# just remove the unnecessary object files
tidy:
	@rm -rf *.o

# remove the object files and the library
clean:
	@/bin/rm -rf *.o
	@/bin/rm -rf lib/$(lib)
