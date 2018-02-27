#!/bin/bash

echo "#####  starting scan..."


arr_partIntDist=(0.15 0.2   0.2   0.25  0.25  0.3   0.3   0.35  0.35  0.4   0.4   0.45  0.45  0.5   0.5   0.55  0.55  0.6   0.65  0.7   0.75  0.8 )
arr_nPart=( 20 15.5  19.5  13  16  11  13.5  9.5  11.5  8   10  7   8.5 6.5 7.5 5.5 7   6   5   4.5 4.5 4   )

COUNTER=0
for partIntDist in ${arr_partIntDist[*]}
do
    #for nPart in ${arr_nPart[*]}
    #do
        nPart=${arr_nPart[$COUNTER]}
        COUNTER=$((COUNTER+1))

        for rangeId in {0..8}
        do
            echo "#####  do $partIntDist $nPart, rangeId=$rangeId..."
            # Nov 2017:
            root -b "launchNucleiCollisionGenerator.C( 20,2,$rangeId,$partIntDist,$nPart)"

            echo $partIntDist $nPart $rangeId `cat outputs_NucleiCollision/tmpTextOutput1_0.txt`  >> OUTPUT_MEAN_N_STRINGS.txt
            echo $partIntDist $nPart $rangeId `cat outputs_NucleiCollision/tmp_cross_section.txt`  >> OUTPUT_ALL_CROSS_SECTIONS.txt
        done
    #done
done
