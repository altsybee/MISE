void launchStringFragmentationGenerator(
//        TString inputFileName_NucleiCollisions = "nuclCollisionsTree_nEv20000.root"
        TString inputFileName_NucleiCollisions = "nuclCollisionsTree_nEv60000.root"
//        TString inputFileName_NucleiCollisions = "nuclCollisionsTree_nEv10000000.root"
//        TString inputFileName_NucleiCollisions = "nuclCollisionsTree_nEv8000000.root"
//        TString inputFileName_NucleiCollisions = "nuclCollisionsTree_nEv1100000.root"
//        TString inputFileName_NucleiCollisions = "nuclCollisionsTree_nEv50000.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv50000.root_StringBoosts_Edensity_0_0001_pp.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv50000.root_StringBoosts_Edensity_0_0001_pp_WITH_3perc_HARD_SCATTERING.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv20000.root_StringBoosts_density_0_00002.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv10000000.root_StringBoosts_density_0_00002.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv8000000.root_StringBoosts_density_0_0001.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv1100000.root_StringBoosts.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv20000.root_StringBoosts_r_3_density_0_00004.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv20000.root_StringBoosts_r_1_density_0_0004.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv20000.root_StringBoosts_r_2_density_0_0001_WITH_D.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv8000000.root_StringBoosts_r_2_density_0_0001_WITH_D.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv8000000.root_StringBoosts_r_2_density_0_00002_WITH_D.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv8000000.root_StringBoosts_r_2_density_0_0005_WITH_D.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv20000.root_StringBoosts_r_2_density_0_00002_WITH_D.root"
//        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv60000.root_StringBoosts_r_2_density_0_0001_WITH_D.root"
        , TString inputFileName_StringBoosts = "nuclCollisionsTree_nEv60000.root_StringBoosts_r_2_density_0_00002_WITH_D.root"

        , int nEvents = -1//100 //)
//        , int nEvents = 1000 //1100000
//        , int nEvents = 20000 //1100000
//, int flagMBorFixedB = 0, float bImpact=0
        , int drawMode = 0 )
{
    //string fragmentation and particle decays
    gROOT->ProcessLine(".L StringDecayer/ParticleDescr.cpp+");
    gROOT->ProcessLine(".L StringDecayer/DecayInTwo.cpp+");
    gROOT->ProcessLine(".L StringDecayer/StringFragmentation.cpp+");
    gROOT->ProcessLine(".L StringDecayer/StringDescr.cpp+");

    gROOT->ProcessLine(".L ManagerStringFragmentation.cpp+");

    //create directory for output files
//    gROOT->ProcessLine(".mkdir tmpOutputs");

    //prepare directories for outputs
    char *strOutputDirName_ManagerStringFragmentation = "outputs_ManagerStringFragmentation";
    gROOT->ProcessLine( Form( ".mkdir %s", strOutputDirName_ManagerStringFragmentation ) );

    gRandom->SetSeed(0);

    TStopwatch timer;
    timer.Start();

    if(1)
    {
        //string decayer
        StringDescr *strDescr = new StringDescr();
        strDescr->setRandomGenerator(gRandom);
//        strDescr->setCoeffPtKickPerUnitMagn( coeffPtKickPerUnitMagn ); //kick in GeV per unit magnitude of string "boost"

        // ########## Event Manager";

        ManagerStringFragmentation evMan;
        evMan.setNumberOfCentralityBins(10);
        evMan.initOutputObjects();

        evMan.setFillEventTree( true );
        evMan.setOutputDirectoryName( strOutputDirName_ManagerStringFragmentation );

        // ##### take nuclei collisions and generate tracks
//        evMan.setInputFileName( inputFileName_NucleiCollisions );
        evMan.setInputFileNameNucleiCollisions( inputFileName_NucleiCollisions );
        evMan.setInputFileNameStringBoosts( inputFileName_StringBoosts );

        evMan.setCutMinNumberOfParticles(2); //20 was used for 2016 studies of Pb-Pb!
        // if have "hard scatterings":
        evMan.setWhatToDoWithHardScattering(1); // 0 - makeTwoJets, 1 - particle pair with random pt from Power law

        evMan.applyFragmentationToEvents( strDescr, nEvents ); //, analysersArray, 0);//nAnalysers );

        if (drawMode==0)
        {
            //d->drawEventStructure(); //draw only last event
        }
        else if (drawMode==1) //batch mode - no histos
        {
            evMan.setDrawHistos(false);
        }
        else if (drawMode==2) //spec mode - draw last event and out
        {
            //d->drawEventStructure();
            gROOT->ProcessLine(".q");
        }
        //        d->drawEventStructure();

    }
//    gROOT->ProcessLine(".q");

    //time estimation
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

}
