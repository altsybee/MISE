#include <iostream>
#include "DistanceEntry.h"
//#include "TRandom3.h"

//#include <algorithm>    // std::sort
//#include <vector>       // std::vector

using namespace std;

const int NMAX = 2000;

void quickSortR( DistanceEntry* a, long N )
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



float matrix[NMAX][NMAX];
DistanceEntry arrDist[NMAX];

DistanceEntry* minDistancePairsFinder( const float *x1, const float *y1, const float *x2, const float *y2, int Nrows, int Ncols )
{

    //occupying particles
    bool occupiedX[NMAX];
    bool occupiedY[NMAX];

    DistanceEntry *d, *d2;

    //make all rows and cols available
    for ( int i = 0; i < Nrows; i++ )
        occupiedX[i] = false;
    for ( int i = 0; i < Ncols; i++ )
        occupiedY[i] = false;

    int iArrDist = 0;

    //calculate matrix of distances
    for ( int i = 0; i < Nrows; i++ )
    {
        for ( int j = 0; j < Ncols; j++ )
        {
            float dx = x1[i] - x2[j];
            float dy = y1[i] - y2[j];
            matrix[i][j] = dx*dx + dy*dy;
            d = &arrDist[iArrDist];
            d->dist = matrix[i][j];
            d->x = i;
            d->y = j;
            d->inInteraction = false;
            iArrDist++;
        }

    }
    for ( int i = 0; i < Nrows * Ncols; i++ )
    {
        d = &arrDist[i];
        cout << d->dist << " ";
        cout << d->x << " ";
        cout << d->y << " ";
        cout << "\n";
    }
    cout << "sorting..." << endl;
    quickSortR( arrDist, Nrows * Ncols-1);
    //check sorting:
    for ( int i = 0; i < Nrows * Ncols; i++ )
    {
        d = &arrDist[i];
        cout << d->dist << " ";
        cout << d->x << " ";
        cout << d->y << " ";
        cout << "\n";
        if ( !occupiedX[d->x] && !occupiedY[d->y] )
        {
            occupiedX[d->x] = true;
            occupiedY[d->y] = true;
            d->inInteraction = true;

            for ( int j = 0; j < Nrows * Ncols; j++ )
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
    cout << endl;

//    cout << "check x,y" << endl;

//    for ( int i = 0; i < Nrows*Ncols; i++ )
//    {
//        cout << arrDist[i].x << " ";
//    }
//    cout << endl;

//    for ( int i = 0; i < Nrows*Ncols; i++ )
//    {
//        cout << arrDist[i].y << " ";
//    }
//    cout << endl;

    cout << "check occupation..." << endl;

    for ( int i = 0; i < Nrows; i++ )
    {
        cout << occupiedX[i] << " ";
    }
    cout << endl;

    for ( int i = 0; i < Nrows; i++ )
    {
        cout << occupiedY[i] << " ";
    }
    cout << endl;

    for ( int i = 0; i < Nrows*Ncols; i++ )
    {
        cout << arrDist[i].x << " ";
    }
    cout << endl;
    for ( int i = 0; i < Nrows*Ncols; i++ )
    {
        cout << arrDist[i].y << " ";
    }
    cout << endl;
    for ( int i = 0; i < Nrows*Ncols; i++ )
    {
        cout << arrDist[i].inInteraction << " ";
    }
    cout << endl;

    cout << "after taking minimums:" << endl;
    for ( int i = 0; i < Nrows * Ncols; i++ )
    {
        d = &arrDist[i];
        if ( !d->inInteraction )
            continue;
        cout << d->dist << " ";
        cout << d->x << " ";
        cout << d->y << " ";
        cout << "\n";
    }

    //print matrix
    cout << ">> matrix of distances:" << endl;
    for ( int i = 0; i < Nrows; i++ )
    {
        for ( int j = 0; j < Ncols; j++ )
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    return arrDist;

}

