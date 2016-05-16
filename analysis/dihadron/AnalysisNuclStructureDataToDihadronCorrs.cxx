#include "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AnalysisSelector.h"
#include "/Users/macbook/alice/simpleAnalysis/analysers/AnalyserBase.h"
#include "/Users/macbook/alice/simpleAnalysis/analysers/diHadronMethod/DiHadronAnalyser.cxx"
#include "/Users/macbook/alice/simpleAnalysis/commonTools/Tools.cxx"

#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

TH2D* getRidge2D( TList *outList, int histId );

//AnalysisSelector *selector;
const int nAnalysers = 2;



DiHadronAnalyser ** PrepareAnalyzers()
{
    DiHadronAnalyser **analysersArray;

    // ########### tune analysers
    double ptEdges[nAnalysers][4] = {
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },

//        { 0.0, 100.0, 0.0, 100.0 },
//        { 0.2, 0.8, 3, 100.0 },
//        { 0.1, 0.7, 0.1, 0.7 },
//        { 4, 10, 4, 10 },
//        { 3, 6, 3, 6 },
//        { 2, 4, 2, 4 },
//        { 6, 10, 6, 10 },
//        { 8, 10, 8, 10 },
//        { 1, 100, 1, 100 },
//        { 0.15, 100, 0.15, 100 },
        { 0.2, 2, 0.2, 2 }, //a-la balance functions
        { 0.2, 2, 0.2, 2 }, //a-la balance functions
//        { 1, 100, 1, 100 },
//        { 1, 100, 1, 100 },
//        { 0.15, 100, 0.15, 100 },
//                { 3.5, 6, 3.5, 6 },
//        { 3., 5, 3., 5 },
//        { 1.5, 2., 2., 2.5 },
//        { 0.7, 100, 0.7, 100 },

//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
//        { 0.0, 10.0 },
        //            { 0.2, 2.0 },
        //            { 0.0, 10.0 },
        //            { 0.2, 2.0 },
    };


    analysersArray = new DiHadronAnalyser*[nAnalysers];
    for ( int anId = 0; anId < nAnalysers; anId++ )
    {
        cout << "test " << anId << endl;
        DiHadronAnalyser *analyser = new DiHadronAnalyser;
        analysersArray[anId] = analyser;

        analyser->SetMaxDeltaEtaForHistBounds(2.0);
//        analyser->SetMaxDeltaEtaForHistBounds(1.6);
//        analyser->SetMaxDeltaEtaForHistBounds(4.0);
//        analyser->SetEtaRanges( -2.0, 2.0, -2.0, 2.0 );
//        analyser->SetEtaRanges( -2.0, 2.0, -2.0, 2.0 );
        analyser->SetEtaRanges( -1.0, 1.0, -1.0, 1.0 );
//        analyser->SetEtaRanges( -0.8, 0.8, -0.8, 0.8 );
        //        analyser->SetEtaRanges( -2.5, -0.5, 0.5, 2.5 );
        //        analyser->SetEtaRanges( etaEdges[anId][0], etaEdges[anId][1], etaEdges[anId][2], etaEdges[anId][3] );
        analyser->SetPtRanges( ptEdges[anId][0], ptEdges[anId][1], ptEdges[anId][2], ptEdges[anId][3] );
        if ( anId == 0 )
            analyser->SetChargesFB(-1,1);
        else if ( anId == 1 )
            analyser->SetChargesFB(1,1);


        analyser->SetShortDef( Form( "list%d", anId ));
        analyser->InitDataMembers();
//        analyser->SetMeanPt(0.438); //tmp for checks?

    }
    return analysersArray;
}

void RunAnalysis( AnalysisSelector *selector, DiHadronAnalyser **analysersArray, int nEventsToProcess, TString strFileData, int centralityClassId = 0 )
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

//    TTree *tree = new TTree( "/Users/macbook/alice/simpleAnalysis/toyAnalysis/NuclearStructure/tmpOutputs/eventTree_nEv200.root/EventTree" );
//    tree->Process ( selector , "", nEventsToProcess );





    TCanvas *canvRidge = new TCanvas("Ridge_plot","ridge plot",10,10,600,600 ); //550,600 );
    canvRidge->Divide(2,2);
    TH2D *hist2D_ridges[nAnalysers];
    for ( int iAn = 0; iAn < nAnalysers; iAn++ )
    {
        //        cout << iAn << endl;
        TVirtualPad * pad = canvRidge->cd(iAn+1);
        pad->SetTheta(45);
        pad->SetPhi(45);//45);

        TList *outList = analysersArray[iAn]->GetOutputList();
        outList->ls();
        //        TH1D* hist = ExtractEtaIntral(outList);
        TH2D* histR = getRidge2D(outList, iAn);
        tuneHistogram2D(histR);
        histR->SetTitle("");
        if ( histR )
//            histR->DrawCopy("lego2");
        histR->DrawCopy("surf1");

        //        delete histR;
        //        cout << hist << endl;
        //        hist->SetMarkerColor(kMyColors[iAn]);
        //        cout << hist << endl;
        if (1)
        {
            //                hist->Draw("P same");
            //                double min = hist->GetMinimum();
            //                double max = hist->GetMaximum();
            //                hist->GetYaxis()->SetRangeUser(0.2,1.8);
        }
        histR->SaveAs( Form("plots/hist2D_an%d_centrClass%d.root", iAn, centralityClassId) );
        hist2D_ridges[iAn] = histR;
    }
//    canvRidge->SaveAs( Form("phiProfile_eventObjects_%.2f.png", rhoPt) );
    canvRidge->SaveAs( Form("corrPlot_rhoPt_%.2f.gif", 111) );

    //subtracted ridges (a-la balance function)
    int nAnToPlot = 0;
    addTH2D_withCoeff(hist2D_ridges[nAnToPlot], hist2D_ridges[nAnToPlot+1], -1);

    TVirtualPad * pad = canvRidge->cd(nAnalysers+1);
    pad->SetTheta(45);
    pad->SetPhi(45);//45);
    hist2D_ridges[nAnToPlot]->DrawCopy("surf1");



}




TH2D* getRidge2D(TList *outList , int histId)
{
    TH2D *histSn = (TH2D*)outList->FindObject("histSn_mult0-100000");
    TH2D *histBn = (TH2D*)outList->FindObject("histBn_mult0-100000");
    //    return outList;

    cout << "histSn->GetEntries()=" << histSn->GetEntries() << endl;
    cout << "histBn->GetEntries()=" << histBn->GetEntries() << endl;

    if ( !histSn->GetEntries() || !histBn->GetEntries() )
        return 0x0;
    smoothBn( histBn ); // !!!!

    TH2D *histR = (TH2D*)histSn->Clone(Form("histR_%d",histId));
//    recalcBinsRelative( histBn, histR );
    calcDihadronRatioLikeInSTAR( histR, histBn );
    const int nBinsPhiShift = 6;
    shiftPhiIn2Dplot( histR, nBinsPhiShift );
    histR = mirrorPlotInEta( histR, nBinsPhiShift );
    //    histR->Draw("surf1");

    return histR;
}



