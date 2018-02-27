#!/bin/bash

echo "#####  starting scan..."

#arr_partIntDist=( 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 )

arr_partIntDist=( 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 )


#arr_partIntDist=( 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45  )
#arr_partIntDist=( 0.2 0.25 0.3 0.35 0.4 0.45  )
#arr_partIntDist=( 0.5 0.55 0.6 0.65 0.7 0.75 0.8 )
#arr_partIntDist=( 0.6 0.65 0.7 0.75 0.8 )

#arr_nPart=( 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 )
#arr_nPart=( 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 20.5 21 21.5 22 22.5 23 23.5 24 24.5 25 )
arr_nPart=( 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 )
#20.5 21 21.5 22 22.5 23 23.5 24 24.5 25 )


#arr_partIntDist=( 0.1 0.15  )
#arr_nPart=( 3 4  )

for partIntDist in ${arr_partIntDist[*]}
do
    for nPart in ${arr_nPart[*]}
    do
        #partIntDist=$(( 0.05 + 0.05*$counter ))
        echo "#####  do $partIntDist $nPart..."
        # Nov 2017:
        #root -l "launchNucleiCollisionGenerator.C( 5000,2,$rangeId,0.4,10)"
        #root -l "launchNucleiCollisionGenerator.C( 5000,2,$rangeId,0.4,9)"

        root -b "launchNucleiCollisionGenerator.C( 20000,0,0,$partIntDist,$nPart)"
       # root -b "launchNucleiCollisionGenerator.C( 4000,0,0,$partIntDist,$nPart)"

        echo $partIntDist $nPart `cat outputs_NucleiCollision/tmpTextOutput1_0.txt`  >> OUTPUT_MEAN_N_STRINGS.txt
        echo $partIntDist $nPart `cat outputs_NucleiCollision/tmp_cross_section.txt`  >> OUTPUT_ALL_CROSS_SECTIONS.txt
    done
done
