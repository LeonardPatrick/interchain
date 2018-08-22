#!/bin/bash


for((i=1; i<=24; i++))
do

sed "s/REP/$i/g" FM\ Green\ Main.cpp >  FM\ Green\ Main\ $i.cpp

icpc -framework Accelerate -o FM\ Green\ Main\ $i FM\ Green\ Main\ $i.cpp FMcoeff.cpp  FMDSF.cpp  FMBS.cpp

./FM\ Green\ Main\ $i > FMGF$i.txt &

done



