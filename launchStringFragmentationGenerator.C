void launchStringFragmentationGenerator( int nEvents = 100, int flagMBorFixedB = 0, float bImpact=0
        , int bMode = 0 )
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
        evMan.setCutMinNumberOfParticles(20);
        evMan.applyFragmentationToEvents( strDescr, nEvents ); //, analysersArray, 0);//nAnalysers );

        if (bMode==0)
        {
            //d->drawEventStructure(); //draw only last event
        }
        else if (bMode==1) //batch mode - no histos
        {
            evMan.setDrawHistos(false);
        }
        else if (bMode==2) //spec mode - draw last event and out
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
