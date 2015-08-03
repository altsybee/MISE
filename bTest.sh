#!/bin/bash


echo "#####  start..."

for meanPartonsInNucleon in {10,} #5,6,7,8,9,10,11}
#for meanPartonsInNucleon in {17,18}
do
    partonInteractionDist=0.25
    clusterFormDist=20 #0.4 #0.3
    stringInteractionDist=3
    for b in {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14}
    do
        fileLabel=pDist${partonInteractionDist}_strDist${stringInteractionDist}_clust${clusterFormDist}_nP${meanPartonsInNucleon}
        echo "#####  impact parameter $b"
        echo points_$fileLabel.txt
    done
done
#root -l "launchNuclearStructure.C( $winId )"
