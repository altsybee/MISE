#include "StringGeneration/DistanceEntry.h"
#include "StringGeneration/MinDistanceFinder.h"
#include "StringGeneration/StringBoosting.h"
#include "ManagerStringBoostsCalculator.h"


void launchStringBoostsCalculator( TString inputFileName = "nuclCollisionsTree_nEv1000.root"
//        , double partonInteractionDist = 0.25
        , double stringInteractionRadius =  0.4 // was 0.03 for QM2018 (and maybe also in 2017 for D-meson)
//        , double strIntDistFractionParA = 10.
//        , double meanPartonsInNucleon = 6
        , double stringOverlapEnergyDensity =  0.001 // was 0.03 for QM2018 (and maybe also in 2017 for D-meson)
//        , double coeffPtKickPerUnitMagn = 0.1
//        , int bMode = 0
        , int nEvents = -1
        )
{
    // dist to form a cluster: now (from Aug 2014) we consider All strings to be in cluster
    double clusterFormDist = 100; //0.4;

    //nuclei generation routine
    gROOT->ProcessLine(".L StringGeneration/DistanceEntry.cpp+");
    gROOT->ProcessLine(".L StringGeneration/MinDistanceFinder.cpp+");
    gROOT->ProcessLine(".L StringGeneration/StringBoosting.cpp+");
    gROOT->ProcessLine(".L ManagerStringBoostsCalculator.cpp+");

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

//    gROOT->ProcessLine(".L ManagerStringBoostsCalculator.cpp+");

    //create directory for output files
//    gROOT->ProcessLine(".mkdir tmpOutputs");

    //prepare directories for outputs
    TString strOutputDirName_StringBoosting = "outputs_NucleiCollision_REWAKE_2023";
    gROOT->ProcessLine( Form( ".! mkdir %s", strOutputDirName_StringBoosting.Data() ) );

    gRandom->SetSeed(0);

    TStopwatch timer;
    timer.Start();

    if(1)
    {
        StringBoosting *booster = new StringBoosting;
        booster->setRandomGenerator(gRandom);
        booster->setOutputDirectoryName( strOutputDirName_StringBoosting.Data() );

        //tune model parameters
        booster->setClusterFormationDistance( clusterFormDist );
        booster->setStringInteractionRadius( stringInteractionRadius );  //setStringInteractionDistance( stringInteractionDist );
        booster->setStringOverlapEnergyDensity( stringOverlapEnergyDensity );

//        booster->setHardScatteringProbability(0);
        booster->setHardScatteringProbability(0);//0.03); // was 0.03 for QM2018 (and maybe also in 2017 for D-meson)

        booster->setComputeStringRepulsion(1);


        // ########## Event Manager";

        ManagerStringBoostsCalculator boostManager;
        boostManager.initOutputObjects();
        boostManager.setInputFileName( inputFileName );
        boostManager.setOutputDirectoryName( strOutputDirName_StringBoosting.Data() ); //strOutputDirName_EventManager );

        booster->initDataMembers();
        booster->setWriteEventViewCanvas(false);

        // ##### generate events
        boostManager.generateBoosts( booster, nEvents ); //, analysersArray, 0);//nAnalysers );

//        booster->drawEventStructure(); //draw only last event
        //        d->drawEventStructure();

    }
//    gROOT->ProcessLine(".q");

    //time estimation
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

}
