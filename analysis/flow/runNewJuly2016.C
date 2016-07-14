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

//    int nEventsToProcess = 100000;//500000;//1000000;//20001;
//    TString fileName("../../outputs_EventManager/eventTree_nEv100000.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_lowPtFromString.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_GaussianMeanPtFromString035.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_GaussianMeanPtFromString_NEW_PIDS_try6_ExpPt.root");
//    TString fileName("../flow_toy/output_toy_events.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_try9_GF_pT_CUT_BOOST_MAG_05.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv30000_try12_TSALLIS.root");
//    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv10000_try13_TSALLIS_r2fm.root");
    TString fileName("../../outputs_ManagerStringFragmentation/eventTree_nEv1000000_PP_TRY3_TSALLIS_r2fm.root");



    // ##### prepare analysers and run
    gROOT->ProcessLine(".L flowAnalyserNewJuly2016.cxx+");

//    int pid = 4;
    flowAnalyserNewJuly2016( fileName, pid );


//    gROOT->ProcessLine(".q");
}


