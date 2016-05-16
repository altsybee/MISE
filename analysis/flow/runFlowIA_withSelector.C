void runFlowIA_withSelector(int centralityClassId = 0)
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



    gROOT->LoadMacro( "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AliSimpleEvent.cxx+g" );

    gROOT->ProcessLine(".L /Users/macbook/alice/simpleAnalysis/commonTools/QAforWindows.cxx+");
    gROOT->ProcessLine(".L /Users/macbook/alice/simpleAnalysis/analysers/AnalyserBase.cxx+");
    gROOT->ProcessLine(".L /Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AnalysisSelector.C+g");


    AnalysisSelector *selector = TSelector::GetSelector("/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AnalysisSelector.C+g");

     /*double*/ //centralityClassId = 16;
    //    selector->SetCentralityRange( centralityClassId, centralityClassId);
    //    int centralityClassId = 12;
    //    selector->SetCentralityRange( centralityClassId, centralityClassId);
//    selector->SetCentralityRange( 4, 10);//19 );
//    selector->SetCentralityRange( 2, 2);
//        selector->SetCentralityRange( centralityClassId, centralityClassId);
//        selector->SetCentralityRange( 9, 9);
        selector->SetCentralityRange( -1, 10);

//        selector->SetNchCuts(5500,100000);
//            selector->SetNchCuts(12*4,100000);
    int nEventsToProcess = 100000;//500000;//1000000;//20001;
//        TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv500000_r2_longStrings.root");
//    TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv400000_r2.root");
//    TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv500000_r1.root");
//    TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv50014.root");
//    TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv20006.root");
//    TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv100000.root");
//    TString fileName("/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv400200.root");
//    TString fileName("../../outputs_EventManager/eventTree_nEv100000.root");
    TString fileName("../../outputs_EventManager/eventTree_nEv9000.root");


    // ##### prepare analysers and run
    gROOT->ProcessLine(".L AnalyserForFlowIA.cxx+");



    gROOT->ProcessLine(".L PrepareAndRun_AnalyseFlow.cxx+");
    AnalyserForFlowIA **analysersArray = PrepareAnalyzers();

    int nAnalysers = 1;
    selector->SetAnalysers(analysersArray, nAnalysers);

    RunAnalysis( selector, analysersArray, nEventsToProcess, fileName
                 , centralityClassId );

//    gROOT->ProcessLine(".q");
}


