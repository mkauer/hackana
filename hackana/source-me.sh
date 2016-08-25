#!/bin/bash

### Fedora 22 version!!

source reset-root.sh

#export ROOTSYS=/home/mkauer/software/root_v5.34.36_x64_fedora20_gcc4.8
#export ROOTSYS=/home/mkauer/software/root_v5.34.36_x64_slc6_gcc5.1
#export ROOTSYS=/home/mkauer/software/root_v5.23.04_i386_slc5_gcc3.4
export ROOTSYS=/home/mkauer/software/root_v5.23.04_x64_slc5_gcc3.4

export R10SYS=/home/mkauer/COSINE/nemo/hackana
export LD_LIBRARY_PATH=${R10SYS}/utils/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export PATH=${ROOTSYS}/bin:${PATH}

### user aliases
hackana=/home/mkauer/COSINE/nemo/hackana
alias epub="emacs $hackana/utils/ana10/anaPublish.C &"
alias cd1e="cd $hackana/1e"
alias cdeg="cd $hackana/eg"
alias cdee="cd $hackana/ee"
alias cdutils="cd $hackana/utils"
alias cdana="cd $hackana/utils/ana10"

alias mana="gmake && ./ana  2>/dev/null"
alias mana6="gmake && ./ana master-model-6.dat  2>/dev/null"
alias mana7="gmake && ./ana master-model-7.dat  2>/dev/null"
alias mana6l="gmake && ./ana master-model-6-lim.dat  2>/dev/null"
alias mana7l="gmake && ./ana master-model-7-lim.dat  2>/dev/null"
alias urun="emacs ./urun.cpp &"


###  needed to do the following  #################
#
# add <iomanip> to uloop and urun
# install all i386 glib things - maybe don't need anymore
# add "-lstdc++" to r10-config libs
# remove -m32 from makefiles
# keep gcc34
# 
##################################################

