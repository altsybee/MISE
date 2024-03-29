void launchNucleiCollisionGenerator( int nEvents = 100, int flagMBorFixedB = 0, float bImpact=0
        , double partonInteractionDist = 0.25
//        , double stringInteractionRadius = 0.4
//        , double strIntDistFractionParA = 10.
        , double meanPartonsInNucleon = 6
//        , double stringOverlapEnergyDensity = 0.1
//        , double coeffPtKickPerUnitMagn = 0.1
        , double nucleonGaussianRadius = 0.4
        , double meanMultFromOneStringForFictiveMultStudies = 1.1
        , bool flagUsePoissonianNpartonsFluctuations = 1
        , bool flagOnlyOneInteractionPerParton = 1
        , bool flagConsiderWoundedPartonsAsStrings = 0
        , int drawMode = 0 )
{
    // dist to form a cluster: now (from Aug 2014) we consider All strings to be in cluster
//    double clusterFormDist = 100; //0.4;

    //nuclei generation routine
    gROOT->ProcessLine(".L StringGeneration/DistanceEntry.cpp+");
    gROOT->ProcessLine(".L StringGeneration/MinDistanceFinder.cpp+");
    gROOT->ProcessLine(".L StringGeneration/NucleiCollision.cpp+");

    //string fragmentation and particle decays
//    gROOT->ProcessLine(".L StringDecayer/ParticleDescr.cpp+");
//    gROOT->ProcessLine(".L StringDecayer/DecayInTwo.cpp+");
//    gROOT->ProcessLine(".L StringDecayer/StringFragmentation.cpp+");
//    gROOT->ProcessLine(".L StringDecayer/StringDescr.cpp+");


    //simple event routine
//    gROOT->LoadMacro( "/Users/macbook/alice/aliceAnalysis/analysisTask/analysis/AliSimpleEvent.cxx+g" );
//    gROOT->LoadMacro( "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AliSimpleEvent.cxx+g" );    // "../../simpleEventAnalyzer/AliSimpleEvent.cxx+g" );

//    gROOT->ProcessLine( ".L /Users/macbook/alice/simpleAnalysis/commonTools/QAforWindows.cxx+" );     // ".L ../../commonTools/QAforWindows.cxx+");
    //    gROOT->ProcessLine(".L ../../analysers/AnalyserBase.cxx+");
    //    gROOT->ProcessLine(".L ../../analysers/diHadronMethod/DiHadronAnalyser.cxx+");
    //    gROOT->ProcessLine(".L ../../analysers/tileCorrelations/TileCorrelations.cxx+");

    gROOT->ProcessLine(".L ManagerNucleiCollisions.cpp+");

    //create directory for output files
//    gROOT->ProcessLine(".mkdir tmpOutputs");

    //prepare directories for outputs
    char *strOutputDirName_NucleiCollision = "outputs_NucleiCollision";
//    char *strOutputDirName_EventManager = "outputs_EventManager";
    gROOT->ProcessLine( Form( ".mkdir %s", strOutputDirName_NucleiCollision ) );
//    gROOT->ProcessLine( Form( ".mkdir %s", strOutputDirName_EventManager ) );

    gRandom->SetSeed(0);

    TStopwatch timer;
    timer.Start();

    if(1)
    {
        NucleiCollision *d = new NucleiCollision;
        d->setRandomGenerator(gRandom);
        d->setOutputDirectoryName( strOutputDirName_NucleiCollision );

        // set nuclei type
//        d->setNucleusType( nucleus_proton );
          d->setNucleusType( nucleus_Pb );
//          d->setNucleusType( nucleus_Au );
//        d->setNucleusType( nucleus_p_Pb );

        // tune impact parameter
        d->setImpactParSpecification( flagMBorFixedB ); //0,1,2 - MB event, precise b, b in range
        d->setImpactParameterByHand( bImpact );
//        d->setImpactParameterRange( 0, 15 );
//        d->setImpactParameterRange( 0.0, 0.1 ); //used for "central" events!
//        d->setImpactParameterRange( 0.8, 1.4 );
//        d->setImpactParameterRange( 0, 13.5 );
//        d->setImpactParameterRange( 4.5, 5 );
        // from https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentralityCodeSnippets:
        // 10 - 20 %	4.96	7.01
//        d->setImpactParameterRange( 5, 7);
        // 0 - 5 %
//        d->setImpactParameterRange( 0, 3.5 );
        // 40 - 50 %
//        d->setImpactParameterRange( 9.92, 11.1 );

        // try 30-50% range (for D meson v2)
//        d->setImpactParameterRange( 8.59, 11.1 );
        // try 20-40% range (for D meson v2)
//        d->setImpactParameterRange( 7.01, 9.92 );

//        d->setImpactParameterRange( 13, 18 );


        // Nov 2017 BY HAND:     (ALICE centrality 2.76 paper: https://arxiv.org/pdf/1301.4361.pdf)
//        d->setImpactParameterRange( 0, 1.57 ); // 0-5%
//        d->setImpactParameterRange( 9.88, 12.09 ); // 40–60%
//        d->setImpactParameterRange( 12.09, 13.97  );  // 60–80%
        if ( flagMBorFixedB == 2) // bImpPar ranges cooked using centrality paper: https://arxiv.org/pdf/1301.4361.pdf
        {
            if (int(bImpact)==0)         d->setImpactParameterRange( 0, 3.50  );  // 0-5%
            else if (int(bImpact)==1)    d->setImpactParameterRange( 3.5, 4.94  );  // 5-10%
            else if (int(bImpact)==2)    d->setImpactParameterRange( 4.94, 6.98  );  // 10-20%
            else if (int(bImpact)==3)    d->setImpactParameterRange( 6.98, 8.55  );  // 20-30%
            else if (int(bImpact)==4)    d->setImpactParameterRange( 8.55, 9.88  );  // 30-40%
            else if (int(bImpact)==5)    d->setImpactParameterRange( 9.88, 11.04  );  // 40-50%
            else if (int(bImpact)==6)    d->setImpactParameterRange( 11.04, 12.09  );  // 50-60%
            else if (int(bImpact)==7)    d->setImpactParameterRange( 12.09, 13.05  );  // 60-70%
            else if (int(bImpact)==8)    d->setImpactParameterRange( 13.05, 13.97  );  // 70-80%

        }


//        d->setImpactParameterRange( 0, 0.01 );
//        d->setImpactParameterRange( 7, 7.2 );
//        d->setImpactParameterRange( 10, 11 );

        //used for "central" PP events!
//        d->setImpactParameterRange( 0.0, 0.7 );

        //tune model parameters
        d->setPartonInteractionDistance( partonInteractionDist );
        //        d->setMaxNumberOfPartons( 4000 );
//        d->setClusterFormationDistance( clusterFormDist );
//        d->setStringInteractionRadius( stringInteractionRadius );  //setStringInteractionDistance( stringInteractionDist );
//        float const parAf_forStringInteractionDist = stringInteractionDist / strIntDistFractionParA;
//        d->setStringInteractionParA( parAf_forStringInteractionDist );
        d->setMeanNofPartonsInNucleon( meanPartonsInNucleon );
//        d->setStringOverlapEnergyDensity(stringOverlapEnergyDensity);

//        d->setHardScatteringProbability(0);//0.03);

        d->setFlagUsePoissonianNpartonsFluctuations( flagUsePoissonianNpartonsFluctuations );
        d->setFlagOnlyOneInteractionPerParton( flagOnlyOneInteractionPerParton );
        d->setFlagConsiderWoundedPartonsAsStrings( flagConsiderWoundedPartonsAsStrings );

        // Nov 2017: EXPERIMENTAL FEATURE!
        // FROM Ulrich Heinz, J. Scott Moreland (2011) - https://arxiv.org/pdf/1108.5379.pdf
        // TABLE 1: sqrt(B) ---> r_transv = sqrt(2B)
        // ASSUME PROPORTIONALITY B/N parameter B and avNpart
//        d->setNucleonGaussianRadius( 0.22 * sqrt(2*meanPartonsInNucleon) ); // 0.1038 is a conversion factor = 0.935/9.  , where avNpart=9.
//        d->setNucleonGaussianRadius( 0.4 ); // 0.1038 is a conversion factor = 0.935/9.  , where avNpart=9.
        d->setNucleonGaussianRadius( nucleonGaussianRadius );


        // Nov 2017 - Feb 2018: av mult from one string = 1.1 particle per unit of rapidity (FOR FICTIVE MULT NucleiCollision-only studies)
//        d->setMeanMultFromOneStringForFictiveMultStudies( 1.1 );
        d->setMeanMultFromOneStringForFictiveMultStudies( meanMultFromOneStringForFictiveMultStudies );

        //d->setComputeStringRepulsion(1);

        //string decayer
//        StringDescr *strDescr = new StringDescr();
//        strDescr->setRandomGenerator(gRandom);
//        strDescr->setCoeffPtKickPerUnitMagn( coeffPtKickPerUnitMagn ); //kick in GeV per unit magnitude of string "boost"

        // ########## Event Manager";

        ManagerNucleiCollisions evMan;
//        evMan.setNumberOfCentralityBins(10);
        evMan.initOutputObjects();

//        evMan.setFillEventTree( true );
        evMan.setOutputDirectoryName( strOutputDirName_NucleiCollision ); //strOutputDirName_EventManager );

        //draw (yes/no) event view canvases to .eps .png
        if (drawMode==0)
            d->setWriteEventViewCanvas(true);
        else
            d->setWriteEventViewCanvas(false);

        // ### JUST BLOCK EVENT CANVASES:
        //        d->setWriteEventViewCanvas(true);
        d->setWriteEventViewCanvas(false);

        // ##### generate events
        evMan.generateEvents( d, nEvents ); //, analysersArray, 0);//nAnalysers );

        if (drawMode==0)
        {
            d->drawEventStructure(); //draw only last event
        }
        else if (drawMode==1) //batch mode - no histos
        {
            evMan.setDrawHistos(false);
        }
        else if (drawMode==2) //spec mode - draw last event and out
        {
            d->drawEventStructure();
            gROOT->ProcessLine(".q");
        }
        //        d->drawEventStructure();

    }

    //time estimation
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

//        gROOT->ProcessLine(".q");

}
