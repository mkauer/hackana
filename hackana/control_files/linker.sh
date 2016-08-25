#!/bin/bash

channels="1e eg ee"
#phases="combined phase123"

#files="master-model-5.dat master-model-6.dat master-model-7.dat"
files="pb211-tl207.dat"

here=`pwd`
for chan in $channels
  do
#  for phase in $phases
#    do
#    there="/unix/nemo2/users/kauer/$chan/zr96_anal/$phase"
    there="/home/kauer/hackana/$chan"
    if [ ! -d "$there" ];then
	echo "$there doesn't exist! "
    else
	for file in $files
	  do
	  echo "$there/$file"
	  /bin/rm -vf "$there/$file" 2>/dev/null
	  ln -sv "$here/$file" "$there/$file"
	done
    fi
#  done
done

