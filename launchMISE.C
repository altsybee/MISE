void launchMISE( int nEvents = 1, int flagMBorFixedB = 0, float bImpact=0
        , double partonInteractionDist = 0.25
        , double stringInteractionDist = 0.4
//        , double strIntDistFractionParA = 10.
        , double meanPartonsInNucleon = 6
        , double stringOverlapEnergyDensity = 0.1
//        , double coeffPtKickPerUnitMagn = 0.1
        , int bMode = 0 )
{
    // dist to form a cluster: now (from Aug 2014) we consider All strings to be in cluster
    double clusterFormDist = 100; //0.4;

    gROOT->ProcessLine(".L StringGeneration/DistanceEntry.cpp+");
    gROOT->ProcessLine(".L StringGeneration/MinDistanceFinder.cpp+");
    gROOT->ProcessLine(".L StringGeneration/NucleusStructure.cpp+");

    gROOT->ProcessLine(".L StringDecayer/ParticleDescr.cpp+");
    gROOT->ProcessLine(".L StringDecayer/DecayInTwo.cpp+");
    gROOT->ProcessLine(".L StringDecayer/StringFragmentation.cpp+");
    gROOT->ProcessLine(".L StringDecayer/StringDescr.cpp+");




    //    gROOT->ProcessLine(".L ../../newLRCan/AliSimpleEvent.cxx+");
//    gROOT->LoadMacro( "/Users/macbook/alice/aliceAnalysis/analysisTask/analysis/AliSimpleEvent.cxx+g" );
    gROOT->LoadMacro( "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AliSimpleEvent.cxx+g" );    // "../../simpleEventAnalyzer/AliSimpleEvent.cxx+g" );

    gROOT->ProcessLine( ".L /Users/macbook/alice/simpleAnalysis/commonTools/QAforWindows.cxx+" );     // ".L ../../commonTools/QAforWindows.cxx+");
    //    gROOT->ProcessLine(".L ../../analysers/AnalyserBase.cxx+");
    //    gROOT->ProcessLine(".L ../../analysers/diHadronMethod/DiHadronAnalyser.cxx+");
    //    gROOT->ProcessLine(".L ../../analysers/tileCorrelations/TileCorrelations.cxx+");

    gROOT->ProcessLine(".L EventManager.cpp+");

    gROOT->ProcessLine(".mkdir tmpOutputs");

    gRandom->SetSeed(0);

    TStopwatch timer;
    timer.Start();

    if(1)
    {
        NucleusStructure *d = new NucleusStructure;
        d->setRandomGenerator(gRandom);

        // set nuclei type
//          d->setNucleusType( nucleus_Pb );
          d->setNucleusType( nucleus_Au );
//        d->setNucleusType( nucleus_proton );

        // tune impact parameter
        d->setImpactParSpecification( flagMBorFixedB ); //0,1,2 - MB event, precise b, b in range
        d->setImpactParameterByHand( bImpact );
//        d->setImpactParameterRange( 0, 15 );
//        d->setImpactParameterRange( 0.0, 0.5 );
        d->setImpactParameterRange( 0, 13.5 );

        //tune model parameters
        d->setPartonInteractionDistance( partonInteractionDist );
        //        d->setMaxNumberOfPartons( 4000 );
        d->setClusterFormationDistance( clusterFormDist );
        d->setStringInteractionDistance( stringInteractionDist );
//        float const parAf_forStringInteractionDist = stringInteractionDist / strIntDistFractionParA;
//        d->setStringInteractionParA( parAf_forStringInteractionDist );
        d->setMeanNofPartonsInNucleon( meanPartonsInNucleon );
        d->setStringOverlapEnergyDensity(stringOverlapEnergyDensity);

        d->setComputeStringRepulsion(1);

        //string decayer
        StringDescr *strDescr = new StringDescr();
        strDescr->setRandomGenerator(gRandom);
//        strDescr->setCoeffPtKickPerUnitMagn( coeffPtKickPerUnitMagn ); //kick in GeV per unit magnitude of string "boost"


        EventManager evMan;
        evMan.setFillEventTree( true );

        //draw (yes/no) event view canvases to .eps .png
        if (bMode==0)
            d->setWriteEventViewCanvas(true);
        else
            d->setWriteEventViewCanvas(false);

        // ### JUST BLOCK EVENT CANVASES:
        //        d->setWriteEventViewCanvas(true);
        d->setWriteEventViewCanvas(false);

        // ##### generate events
        evMan.generateEvents( d, strDescr, nEvents ); //, analysersArray, 0);//nAnalysers );

        if (bMode==0)
        {
            d->drawEventStructure(); //draw only last event
        }
        else if (bMode==1) //batch mode - no histos
        {
            evMan.setDrawHistos(false);
        }
        else if (bMode==2) //spec mode - draw last event and out
        {
            d->drawEventStructure();
            gROOT->ProcessLine(".q");
        }
        //        d->drawEventStructure();

    }
//    gROOT->ProcessLine(".q");


    if(0)
    {
        gROOT->ProcessLine(".L DistanceEntry.cpp+");
        gROOT->ProcessLine(".L MinDistanceFinder.cpp+");
        const int nRows = 4;
        const int nCols = 4;
        float x1[nRows] = {5,3,1,7};
        float y1[nCols] = {0,4,2,3};
        float x2[nRows] = {3,2,4,5};
        float y2[nCols] = {11,4,2,8};

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

        MinDistanceFinder minDistanceFinder;
        minDistanceFinder.minDistancePairsFinder(x1,y1,x2,y2,nRows,nCols);
    }

    //time estimation
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

}
