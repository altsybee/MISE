#export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

g++ -c MinDistanceFinder.cpp DistanceEntry.cpp
g++ -dynamiclib MinDistanceFinder.o DistanceEntry.o -o MinDistanceFinder.dylib

g++ -c testMinDistanceFinder.cpp
g++ -v -Wall testMinDistanceFinder.o ./MinDistanceFinder.dylib -o testMinDistanceFinder
