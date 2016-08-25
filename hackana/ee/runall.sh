#!/bin/bash

modes="m15 m8 m7 m6 m5 m2 m1_1p m1"
reuse="master-model-6-lim.dat"
urun="urun.cpp"

for mode in $modes
do
    echo -e "\n\n Working on ==> $mode \n"
    echo -e "==========================================\n\n"
    
    sed -i -c /2b0n/s/LR/BR/g $reuse
    find="2b0n_${mode}_zr96"
    sed -i -c /$find/s/BR/LR/ $reuse
    sed -i -c s/string\ decaymode=\".*\"\;/string\ decaymode=\"$mode\"\;/ $urun
    gmake && ./ana $reuse -b
    cp "ana10.dat" "${mode}_log.dat"
    
done


echo -e "\n\n=================================================="
echo -e     "             DONE WITH MASTER SCRIPT              "
echo -e     "==================================================\n\n"

