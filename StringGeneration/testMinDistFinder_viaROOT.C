void testMinDistFinder_viaROOT()
{
        gROOT->ProcessLine(".L DistanceEntry.cpp+");
        gROOT->ProcessLine(".L MinDistanceFinder.cpp+");


        cout << "Hello!" << endl;

        const int nRows = 50;
        const int nCols = 50;
        float x1[nRows]; // = {5,3,2,5,3,1,1,7};
        float y1[nCols]; // = {0,4,2,5,3,1,7,3};
        float x2[nRows]; // = {3,2,4,5,5,3,1,7};
        float y2[nCols]; // = {11,4,2,8,7,6,10};

        for ( int i = 0; i < nRows; i++ )
        {
            x1[i] = rand() % 100 * 0.01;
            x2[i] = rand() % 100 * 0.01;
        }
        for ( int i = 0; i < nCols; i++ )
        {
            y1[i] = rand() % 100 * 0.01;
            y2[i] = rand() % 100 * 0.01;
        }

        cout << ">> x1:  ";
        for ( int i = 0; i < nRows; i++ )
            cout << x1[i] << " ";
        cout << endl;
        cout << ">> y1:  ";
        for ( int i = 0; i < nCols; i++ )
            cout << y1[i] << " ";
        cout << endl;

        cout << ">> x2:  ";
        for ( int i = 0; i < nRows; i++ )
            cout << x2[i] << " ";
        cout << endl;
        cout << ">> y2:  ";
        for ( int i = 0; i < nCols; i++ )
            cout << y2[i] << " ";
        cout << endl;

        //    srand(time(NULL));
        //    for ( int i = 0; i < 1000; i++ )
        //        cout << rand() % 100 << " ";
        //    cout << endl;


        MinDistanceFinder minDistanceFinder(0.25);
        for ( int i = 0; i < 10000; i++ )
        {
    //        for ( int i = 0; i < nRows; i++ )
    //        {
    //            x1[i] = rand() % 100 * 0.01;
    //            x2[i] = rand() % 100 * 0.01;
    //        }
    //        for ( int i = 0; i < nCols; i++ )
    //        {
    //            y1[i] = rand() % 100 * 0.01;
    //            y2[i] = rand() % 100 * 0.01;
    //        }
            minDistanceFinder.FindMinDistancesBetweenPairs(x1,y1,x2,y2,nRows,nCols);
    //        minDistanceFinder.PrintMinDistanceMatrix();
        }



//        const int nRows = 4;
//        const int nCols = 4;
//        float x1[nRows] = {5,3,1,7};
//        float y1[nCols] = {0,4,2,3};
//        float x2[nRows] = {3,2,4,5};
//        float y2[nCols] = {11,4,2,8};

//        cout << ">> x1:  ";
//        for ( int i = 0; i < nRows; i++ )
//            cout << x1[i] << " ";
//        cout << endl;
//        cout << ">> y1:  ";
//        for ( int i = 0; i < nCols; i++ )
//            cout << y1[i] << " ";
//        cout << endl;

//        cout << ">> x2:  ";
//        for ( int i = 0; i < nRows; i++ )
//            cout << x2[i] << " ";
//        cout << endl;
//        cout << ">> y2:  ";
//        for ( int i = 0; i < nCols; i++ )
//            cout << y2[i] << " ";
//        cout << endl;

//        MinDistanceFinder minDistanceFinder;
//        minDistanceFinder.FindMinDistancesBetweenPairs(x1,y1,x2,y2,nRows,nCols);
}
