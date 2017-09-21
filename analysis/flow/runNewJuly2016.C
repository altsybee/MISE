void runNewJuly2016(int pid = -1)
{
//    gStyle->SetOptStat(kFALSE);

    // FROM FLOW PACKAGE:
    gSystem->AddIncludePath("-Wno-deprecated");
    gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

    // Load the needed libraries for root (in AliRoot already loaded)
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
    gSystem->Load("libPWGflowBase");

//    TString fileDir("../../outputs_ManagerStringFragmentation");
    TString fileDir("../flow_toy");

//    int nEventsToProcess = 100000;//500000;//1000000;//20001;
//    TString fileName("../../outputs_EventManager/eventTree_nEv100000.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_lowPtFromString.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_GaussianMeanPtFromString035.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_GaussianMeanPtFromString_NEW_PIDS_try6_ExpPt.root");
//    TString fileName("../flow_toy/output_toy_events.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_try9_GF_pT_CUT_BOOST_MAG_05.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv30000_try12_TSALLIS.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_try13_TSALLIS_r2fm.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv1000000_PP_TRY3_TSALLIS_r2fm.root");

//    TString fileName("nuclCollisionsTree_nEv20000_cImpPar5_7.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    TString fileName("nuclCollisionsTree_nEv20000_cImpPar5_7_r3fm.root_StringFragm_nEv20000_TSALLIS_try1.root");
//    TString fileName( "nuclCollisionsTree_nEv10000_cImpPar0_3.5_r3fm.root_StringFragm_nEv10000_TSALLIS_try1.root" );

//    TString fileName( "nuclCollisionsTree_nEv40000_cImpPar9.92_11.1_r2fm.root_StringFragm_nEv40000_TSALLIS_try1.root" );
//    TString fileName( "nuclCollisionsTree_nEv20000_cImpPar5_7_r2fm.root_StringFragm_nEv20000_TSALLIS_try1.root" );
//    TString fileName( "nuclCollisionsTree_nEv10000_cImpPar0_3.5_r2fm.root_StringFragm_nEv10000_TSALLIS_try1.root" );

    TString fileName("output_toy_100_particles_10k_events.root");


    // ##### prepare analysers and run
    gROOT->ProcessLine(".L flowAnalyserNewJuly2016.cxx+");

//    int pid = 4;
    flowAnalyserNewJuly2016( fileDir, fileName, pid );


    gROOT->ProcessLine(".q");
}


