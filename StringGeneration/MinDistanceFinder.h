#ifndef MINDISTANCEFINDER_H
#define MINDISTANCEFINDER_H

#include "DistanceEntry.h"

const int NMAX = 5000;

class MinDistanceFinder
{
    //current matrix of distances
    int fNrows;
    int fNcols;
    float fMatrix[NMAX][NMAX];

    DistanceEntry arrDist[NMAX*NMAX];
    int fSizeOfArrayWithMinimums;
//    int fNumberOfInteractions;
    float fMinDistance;

    //occupying particles flags
    bool occupiedX[NMAX];
    bool occupiedY[NMAX];

    void quickSortR( DistanceEntry* a, long N );

public:
    MinDistanceFinder(float minDistance = 100000);
    void FindMinDistancesBetweenPairs( const float *x1, const float *y1, const float *x2, const float *y2, int Nrows, int Ncols,
                                       bool flagOnlyOneInteractionPerParton );
    // not really used now? Was used to find clusters of strings, as far as I remember... (comment from May 10, 2018)
    void FindMinDistancesWithinArray( const float *x, const float *y, int arrSize );
    DistanceEntry* GetMinDistanceArray( int &arraySize ) { arraySize = fSizeOfArrayWithMinimums; return arrDist; }
//    int GetNumberOfInteractions() { return fNumberOfInteractions; }
    void PrintMinDistanceMatrix();

};

#endif // MINDISTANCEFINDER_H
