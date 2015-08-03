#!/bin/bash

echo "#####  starting generating..."

#for centrClassId in {11,12,13}
#for b in {16..0}
for nPartons in {3..20}
do
echo "#####  do b ${b}..."
#echo "#####  do b ${b}..."
#root -l "launchNuclearStructure.C( 1000,1,$b,0.2,1.,15,0.002,0)"
#root -l "launchNuclearStructure.C( 1000,1,$b,0.2,2.,15,0.0001,0)"
#root -l "launchNuclearStructure.C( 1000,1,$b,0.2,0.25,15,0.2,0)"
root -l "launchNuclearStructure.C( 1000,1,0,0.2,2.,$nPartons,0.0001,0)"
done
