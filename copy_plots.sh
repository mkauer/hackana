#!/bin/bash

hackana=/home/kauer/hackana
thesis=/home/kauer/myLatex/das_thesis/n3pics

for chan in "1e" "eg" "ee"
do
    here="$hackana/$chan/plots"
    there="$thesis/anal-$chan"
    /bin/cp -fv $here/*.png $there/
done

