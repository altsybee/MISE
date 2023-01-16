#include <iostream>
#include "math.h"
//#include "TRandom3.h"
//#include <algorithm>    // std::sort
//#include <vector>       // std::vector
#include "DistanceEntry.h"
#include "MinDistanceFinder.h"


using namespace std;


MinDistanceFinder::MinDistanceFinder(float minDistance)
{
    fNrows = 0;
    fNcols = 0;

    fSizeOfArrayWithMinimums = 0;
    fMinDistance = minDistance;
}

void MinDistanceFinder::quickSortR( DistanceEntry* a, long N )
{
    // На входе - массив a[], a[N] - его последний элемент.

    long i = 0, j = N;        // поставить указатели на исходные места
    DistanceEntry temp, p;

    p = a[ N>>1 ];        // центральный элемент

    // процедура разделения
    do {
        while ( a[i] < p ) i++;
        while ( a[j] > p ) j--;

        if (i <= j) {
            temp = a[i]; a[i] = a[j]; a[j] = temp;
            i++; j--;
        }
    } while ( i<=j );

    // рекурсивные вызовы, если есть, что сортировать
    if ( j > 0 ) quickSortR(a, j);
    if ( N > i ) quickSortR(a+i, N-i);
}


void MinDistanceFinder::FindMinDistancesBetweenPairs(const float *x1, const float *y1
                                                      , const float *x2, const float *y2, int Nrows, int Ncols , bool flagOnlyOneInteractionPerParton )
{
    fNrows = Nrows;
    fNcols = Ncols;
    DistanceEntry *d, *d2;

    //make all rows and cols available
    for ( int i = 0; i < Nrows; i++ )
        occupiedX[i] = false;
    for ( int i = 0; i < Ncols; i++ )
        occupiedY[i] = false;

    //calculate matrix of distances, write array of distances
    int nArrDist = 0;
    for ( int i = 0; i < Nrows; i++ )
    {
        for ( int j = 0; j < Ncols; j++ )
        {
            float dx = x1[i] - x2[j];
            float dy = y1[i] - y2[j];
            if ( fabs(dx) > fMinDistance || fabs(dy) > fMinDistance )
                continue;
            fMatrix[i][j] = dx*dx + dy*dy;
            d = &arrDist[nArrDist];
            d->dist = fMatrix[i][j];
            d->x = i;
            d->y = j;
            d->inInteraction = false;
            nArrDist++;
        }
    }
    fSizeOfArrayWithMinimums = nArrDist;

    quickSortR( arrDist, nArrDist-1 ); //Nrows * Ncols-1);

    //check sorting:
    for ( int i = 0; i < nArrDist; i++ ) //Nrows * Ncols; i++ )
    {
        d = &arrDist[i];
        if (0)
        {
            cout << sqrt(d->dist) << " ";
            cout << d->x << " ";
            cout << d->y << " ";
            cout << "\n";
        }
        //check occupation
        // !!! ADDED FLAG ON MAY 2018: if true - works as before, if false - partons not "occupied", each parton can interact with many others, like in WQM (to be confirmed that it works like desired)
        if ( flagOnlyOneInteractionPerParton )
        {
            if ( occupiedX[d->x] || occupiedY[d->y] )
                continue;
        }

        //check distance (HERE DIST IS IN ^2!)
        if ( d->dist < fMinDistance*fMinDistance ) // !occupiedX[d->x] && !occupiedY[d->y] )
        {
//            cout << d->dist << endl;
            occupiedX[d->x] = true;
            occupiedY[d->y] = true;
            d->inInteraction = true;
//            fNumberOfInteractions++;

            for ( int j = 0; j < nArrDist; j++ ) //Nrows * Ncols; j++ )
            {
                if ( i == j )
                    continue;
                d2 = &arrDist[j];
                if ( d->x == d2->x )
                    occupiedX[d2->x] = true;
                if ( d->y == d2->y )
                    occupiedY[d2->y] = true;
            }
        }
    }

    if (0)
    {
        cout << endl;
        cout << "check occupation..." << endl;
        for ( int i = 0; i < Nrows; i++ )
            cout << occupiedX[i] << " ";
        cout << endl;

        for ( int i = 0; i < Nrows; i++ )
            cout << occupiedY[i] << " ";
        cout << endl;

        if(0)
        {
            for ( int i = 0; i < Nrows*Ncols; i++ )
                cout << arrDist[i].x << " ";
            cout << endl;

            for ( int i = 0; i < Nrows*Ncols; i++ )
                cout << arrDist[i].y << " ";
            cout << endl;

            for ( int i = 0; i < Nrows*Ncols; i++ )
                cout << arrDist[i].inInteraction << " ";
            cout << endl;
        }

        cout << "after taking minimums:" << endl;
    }

    // May 2018: total number of occupied (=wounded) partons in nuclA and nuclB - to compare with wounded quark model:
    nWoundedPartonsInA = 0;
    nWoundedPartonsInB = 0;
    for ( int i = 0; i < Nrows; i++ )
        if ( occupiedX[i] )
            nWoundedPartonsInA++;
    for ( int i = 0; i < Ncols; i++ )
        if ( occupiedY[i] )
            nWoundedPartonsInB++;
    if(0)
        cout << "nPartonsA = " << Nrows << ", nPartonsB = " << Ncols
             << ", nWoundedPartonsInA = " << nWoundedPartonsInA << ", nWoundedPartonsInB = " << nWoundedPartonsInB << endl;
}


void MinDistanceFinder::FindMinDistancesWithinArray( const float *x, const float *y, int arrSize )
{
    //find min 2D-distances for a bunch of particles
    //remark: use variables for rows, cols are not used!

//    cout << ">>>>> same array dist calc:" << endl;
    fNrows = arrSize;
    DistanceEntry *d, *d2;

    //make all rows available
    for ( int i = 0; i < arrSize; i++ )
        occupiedX[i] = false;

    //calculate matrix of distances, write array of distances
    int nArrDist = 0;
    for ( int i = 0; i < arrSize; i++ )
    {
        for ( int j = 0; j < arrSize; j++ )
        {
            if ( i == j ) //same particle
                continue;
            float dx = x[i] - x[j];
            float dy = y[i] - y[j];
            if ( fabs(dx) > fMinDistance || fabs(dy) > fMinDistance )
                continue;
            fMatrix[i][j] = dx*dx + dy*dy;
            d = &arrDist[nArrDist];
            d->dist = fMatrix[i][j];
            d->x = i;
            d->y = j;
            d->inInteraction = false;
            nArrDist++;
        }
    }
    fSizeOfArrayWithMinimums = nArrDist;

    quickSortR( arrDist, nArrDist-1 ); //Nrows * Ncols-1);

    //check sorting:
    for ( int i = 0; i < nArrDist; i++ ) //Nrows * Ncols; i++ )
    {
        d = &arrDist[i];
        if (0)
        {
            cout << sqrt(d->dist) << " ";
            cout << d->x << " ";
            cout << d->y << " ";
            cout << "\n";
        }
        //check occupation
        if ( occupiedX[d->x] || occupiedX[d->y] )
            continue;

        //check distance
        if ( d->dist < fMinDistance*fMinDistance ) // !occupiedX[d->x] && !occupiedY[d->y] )
        {
            occupiedX[d->x] = true;
            occupiedX[d->y] = true;
            d->inInteraction = true;
//            fNumberOfInteractions++;

            for ( int j = 0; j < nArrDist; j++ ) //Nrows * Ncols; j++ )
            {
                if ( i == j )
                    continue;
                d2 = &arrDist[j];
                if ( d->x == d2->x )
                    occupiedX[d2->x] = true;
                if ( d->y == d2->y )
                    occupiedX[d2->y] = true;
            }
        }
    }



}


void MinDistanceFinder::PrintMinDistanceMatrix()
{
    //print matrix
    cout << ">> TEST matrix of distances:" << endl;
    for ( int i = 0; i < fNrows; i++ )
    {
        for ( int j = 0; j < fNcols; j++ )
        {
            cout << fMatrix[i][j] << " ";
        }
        cout << endl;
    }
}
