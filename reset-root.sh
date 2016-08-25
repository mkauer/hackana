#!/bin/bash

# you must source this script for it to take effect

# PATH
temps=`echo $PATH | sed 's/:/\n/g' | uniq | grep -ve "root.*"`
unset PATH
for temp in $temps
  do
  PATH=$temp:${PATH}
done
export PATH

# LD_LIBRARY_PATH
temps=`echo $LD_LIBRARY_PATH | sed 's/:/\n/g' | uniq | grep -ve "root.*"`
unset LD_LIBRARY_PATH
for temp in $temps
  do
  LD_LIBRARY_PATH=$temp:${LD_LIBRARY_PATH}
done
export LD_LIBRARY_PATH

