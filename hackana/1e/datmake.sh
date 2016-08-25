#!/bin/bash

infile="ana10.dat"
searches="foils.zr96 nd150 ca48"
here=`pwd`
#there="/unix/nemo2/users/kauer/calc_sumbkgs"
there=`pwd`

for search in $searches
  do
  cat $infile | grep -e foils -e 2b2n | grep -e $search | awk '{printf "%f\t%f\n",$10,$12}' > "$there/$search.dat"
done
cat $infile | grep -e exbg -e scin -e wire -e sfoil | awk '{printf "%f\t%f\n",$10,$12}' > "$there/exbg.dat"


echo ""
echo "  INTERNALS"
echo "====================================="
for search in $searches
  do
  echo " $search "
  echo "-----------------------"
  cat $infile | grep -e foils -e 2b2n | grep -e $search | awk '{printf "%25s   %f\n",$1,$10}'
  echo ""
done

echo ""
echo "  EXTERNALS"
echo "====================================="
cat $infile | grep -e exbg -e scin -e wire -e sfoil | awk '{printf "%25s   %f\n",$1,$10}'
echo ""
echo ""


cd $there
root -l -q -b sumbkgs.cxx
cd $here


