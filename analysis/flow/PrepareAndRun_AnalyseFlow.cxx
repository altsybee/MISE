#include "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AnalysisSelector.h"
#include "/Users/macbook/alice/simpleAnalysis/analysers/AnalyserBase.h"
//#include "/Users/macbook/alice/simpleAnalysis/analysers/diHadronMethod/DiHadronAnalyser.cxx"
#include "AnalyserForFlowIA.h"
#include "/Users/macbook/alice/simpleAnalysis/commonTools/Tools.cxx"

#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

//TH2D* getRidge2D( TList *outList, int histId );

//AnalysisSelector *selector;
const int nAnalysers = 1;



AnalyserForFlowIA ** PrepareAnalyzers()
{
    AnalyserForFlowIA **analysersArray;
    analysersArray = new AnalyserForFlowIA*[nAnalysers];
    for ( int anId = 0; anId < nAnalysers; anId++ )
    {
        AnalyserForFlowIA *analyser = new AnalyserForFlowIA;
        analysersArray[anId] = analyser;

//        analyser->SetMaxDeltaEtaForHistBounds(4.0);
////        analyser->SetEtaRanges( -2.0, 2.0, -2.0, 2.0 );
//        analyser->SetEtaRanges( -2.0, 2.0, -2.0, 2.0 );
////        analyser->SetEtaRanges( -1.0, 1.0, -1.0, 1.0 );
//        //        analyser->SetEtaRanges( -2.5, -0.5, 0.5, 2.5 );
//        //        analyser->SetEtaRanges( etaEdges[anId][0], etaEdges[anId][1], etaEdges[anId][2], etaEdges[anId][3] );
//        analyser->SetPtRanges( ptEdges[anId][0], ptEdges[anId][1], ptEdges[anId][2], ptEdges[anId][3] );

        analyser->SetShortDef( Form( "list%d", anId ));
        // !!! PID!
        analyser->SetPID(-1);
        analyser->InitDataMembers();
//        analyser->SetMeanPt(0.438); //tmp for checks?

    }
    return analysersArray;
}


void RunAnalysis( AnalysisSelector *selector, AnalyserForFlowIA **analysersArray, int nEventsToProcess, TString strFileData, int centralityClassId = 0 )
{
    //    selector->SetCentralityRange(10,40);
//    selector->SetCentralityRange( 0.9, 1.1 );


    // ###### new!!!! filter bits!!!
    selector->SetFilterMap(768);

    // ##### extract multiplicity boundaries
//    TString strFileData(  );
    TFile *fileData = new TFile( strFileData );
    TH1D *fHistMultClassBoundaries = (TH1D*)fileData->Get( "fHistMultiplicityBoundaries" );

    int nCentralityBins = fHistMultClassBoundaries->GetNbinsX();
    Float_t *centrMultBounds = new Float_t[nCentralityBins];
    for ( int iBin = 0; iBin < fHistMultClassBoundaries->GetNbinsX(); iBin++ )
    {
        centrMultBounds[iBin] = fHistMultClassBoundaries->GetBinContent( iBin+1 );
    }

    selector->SetMultiplicityClasses( centrMultBounds, nCentralityBins );


    // ###### output file name
    //    selector->SetOuputFileName( Form("output_phi_%d", phiShiftId ) );
    selector->SetOuputFileName( Form( "output_nAn%d_nEv%d", nAnalysers, nEventsToProcess ) );

    // ##### run analysis
    TChain *chain = new TChain( "my_data", "My Chain for Example N-Tuple" );
    chain->Add( Form( "%s/EventTree", strFileData.Data() ) );
    chain->Process ( selector , "", nEventsToProcess );
}






