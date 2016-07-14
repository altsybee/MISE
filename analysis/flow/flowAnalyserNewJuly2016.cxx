//#include "AnalyserForFlowIA.h"
//#include "/Users/macbook/alice/simpleAnalysis/commonTools/SimpleTrack.cxx"
//#include "/Users/macbook/alice/simpleAnalysis/commonTools/MiniEvent.cxx"


// aliroot includes
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackSimpleCuts.h"

#include "AliFlowAnalysisWithMCEventPlane.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "AliFlowAnalysisWithLeeYangZeros.h"
#include "AliFlowAnalysisWithScalarProduct.h"

#include "AliFlowCommonHistResults.h"



#include "Riostream.h"
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
//#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include "TRandom3.h"


bool calc_qc_v2 = 1;
bool calc_qc_v3 = 0;//true;
bool calc_lyz = 0;
bool calc_sp = 1;

// k) Define the ranges for two subevents separated with eta gap (needed only for SP method):
Double_t etaMinA = -3; // minimum eta of subevent A
Double_t etaMaxA = -2; // maximum eta of subevent A
Double_t etaMinB = 2; // minimum eta of subevent B
Double_t etaMaxB = 3; // maximum eta of subevent B



TRandom3 pRandom;



const int NMaxTrack = 15000;

Int_t fNumberOfTracks;
Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack],fTrackPhi[NMaxTrack]; //
Int_t fTrackCharge[NMaxTrack];
Int_t fTrackPID[NMaxTrack];

inline void fixDeltaPhi( Double_t &dphi)
{
    if ( dphi < 0 )
        dphi = dphi + 2*TMath::Pi();
}

void ConnectTreeToVars(TTree* tree)
{
    if(!tree) return;

    tree->SetBranchAddress("fNumberOfTracks",&fNumberOfTracks);
    tree->SetBranchAddress("fTrackPt",fTrackPt);
    tree->SetBranchAddress("fTrackPhi",fTrackPhi);
    tree->SetBranchAddress("fTrackEta",fTrackEta);
    tree->SetBranchAddress("fTrackCharge",fTrackCharge);
    tree->SetBranchAddress("fTrackPID",fTrackPID);
}

void flowAnalyserNewJuly2016( TString fileName, Int_t fPID = -1 )
{
    TFile* inputFile = 0x0;
    inputFile = new TFile( fileName );

    if ( !inputFile )
    {
        cout << "NO INPUT FILE!" << endl;
        return;
    }

    TTree *fTree = (TTree*)inputFile->Get("EventTree");
    ConnectTreeToVars( fTree );


//    Int_t fNtracks;

//    Int_t fNeventsQA;

    //centrality class
//    Float_t fMinCentrality;    // min bound on centrality percentile
//    Float_t fMaxCentrality;    // max bound on centrality percentile
//    Float_t fEventCentrality; //event centrality

//    TH1D    *fHistNparticlesInTile; //!
//    MiniEvent *fMiniEvent;


    AliFlowTrackSimpleCuts *cutsRP;
    AliFlowTrackSimpleCuts *cutsPOI;

    AliFlowAnalysisWithQCumulants* qc_v2;
    AliFlowAnalysisWithQCumulants* qc_v3;
    AliFlowAnalysisWithLeeYangZeros *lyz;

    AliFlowTrackSimpleCuts *cutsPOI_forSP;
    AliFlowAnalysisWithScalarProduct *sp;







//    fMiniEvent = new MiniEvent;

    cutsRP = new AliFlowTrackSimpleCuts();
    cutsRP->SetPtMin(0.1);
    cutsRP->SetPtMax(4.0);
    cutsRP->SetEtaMin(-3.0);
    cutsRP->SetEtaMax(3.0);//-0.5);

    cutsPOI = new AliFlowTrackSimpleCuts();
    cutsPOI->SetPtMin(0.1);
    cutsPOI->SetPtMax(4.0);
    //    cutsPOI->SetEtaMin(-2.0);//0.5);
    //    cutsPOI->SetEtaMax(2.0);
//    cutsPOI->SetEtaMin(-1);
//    cutsPOI->SetEtaMax(1);

    cutsPOI->SetEtaMin(-3);
    cutsPOI->SetEtaMax(3);


    //    cutsPOI_forSP = new AliFlowTrackSimpleCuts();
    //    cutsPOI_forSP->SetPtMin(0.2);
    //    cutsPOI_forSP->SetPtMax(4.0);
    //    cutsPOI_forSP->SetEtaMin(-0.5);
    //    cutsPOI_forSP->SetEtaMax(0.5);

    // Initialize the flow methods for default analysis:
    //    AliFlowAnalysisWithMCEventPlane *mcep = new AliFlowAnalysisWithMCEventPlane();
    //mcep->SetHarmonic(2); // default is v2
    //    mcep->Init();

    if ( calc_qc_v2 )
    {
        qc_v2 = new AliFlowAnalysisWithQCumulants();
        qc_v2->SetHarmonic(2); // default is v2
        qc_v2->Init();
    }

    if ( calc_qc_v3 )
    {
        qc_v3 = new AliFlowAnalysisWithQCumulants();
        qc_v3->SetHarmonic(3); // default is v2
        qc_v3->Init();
    }

    if ( calc_lyz )
    {
        lyz = new AliFlowAnalysisWithLeeYangZeros();
        lyz->Init();
    }

    if ( calc_sp )
    {
        sp = new AliFlowAnalysisWithScalarProduct();
        sp->SetUsePhiWeights(kFALSE);
        sp->SetApplyCorrectionForNUA(kFALSE);
        sp->SetHarmonic(2); // default is v2
        sp->Init();
    }




    //    fEventCount = 0;
    //    fFlagHavePrevEvent = kFALSE;
    TList *fOutList = new TList();
    fOutList->SetName("flowData");


    int nEvents = fTree->GetEntries();
    cout << "nEvents in tree: " << nEvents << endl;
//    return;

    for ( int ev = 0; ev < nEvents; ev++)
    {
        if ( ev % 10 == 0 )
            cout << "analysing " << (int)ev << "\r"; cout.flush();

        fTree->GetEntry(ev);

        // ####### CUT ON NUMBER OF TRACKS IN EVENT!
        if ( fNumberOfTracks < 20 )
            continue;

        AliFlowEventSimple* flowEvent = new AliFlowEventSimple(fNumberOfTracks);
        for( int tr = 0; tr < fNumberOfTracks; tr++ )
        {
            double pt = fTrackPt[tr];
            double eta = fTrackEta[tr];
            double phi = fTrackPhi[tr];
            int charge = fTrackCharge[tr];
            int pid = fTrackPID[tr];


            AliFlowTrackSimple* flowtrack = new AliFlowTrackSimple();
            flowtrack->SetPhi(phi);//tracks[i].phi);
            flowtrack->SetEta(eta);
            flowtrack->SetPt(pt);
            flowtrack->SetCharge(charge);

            if (cutsRP)
                flowtrack->TagRP(cutsRP ->PassesCuts(flowtrack));
            if (cutsPOI)
            {
                bool isPOI = cutsPOI->PassesCuts(flowtrack);
                bool isPID_we_want =  ( fPID == -1 || fPID == pid );

                flowtrack->TagPOI( isPOI && isPID_we_want); //tag POIs
            }

            //        if (calc_sp)
            //        {
            //            if( flowtrack->Eta() >= etaMinA && flowtrack->Eta() < etaMaxA )
            //                flowtrack->SetForSubevent(0);
            //            if( flowtrack->Eta() >= etaMinB && flowtrack->Eta() < etaMaxB )
            //                flowtrack->SetForSubevent(1);
            //        }
            flowEvent->AddTrack(flowtrack);

            //for SP
            //        if (cutsPOI_forSP)
            //            flowtrack->TagPOI(cutsPOI_forSP->PassesCuts(flowtrack)); //tag POIs
            //        flowEventSP->AddTrack(flowtrack);

        }


        if ( calc_qc_v2 )   qc_v2->Make(flowEvent);
        if ( calc_qc_v3 )   qc_v3->Make(flowEvent);
        if ( calc_lyz )     lyz->Make(flowEvent);

        if ( calc_sp )
        {
            flowEvent->TagSubeventsInEta(etaMinA,etaMaxA,etaMinB,etaMaxB);
            sp->Make(flowEvent);
        }
        delete flowEvent;


    } // end of event loop



    // ##### terminate

    cout << "FillAliFlowEvent: nEvents = " << fTree->GetEntries() << endl;

    // calculate the final results
    //    fOutputFileId = -1000;
    //    mcep->Finish();
    if ( calc_qc_v2 )
    {
        qc_v2->Finish();
//        qc_v2->WriteHistograms("outputCumulants_v2.root");
        qc_v2->WriteHistograms("outputQC_100k_pp.root");
    }

    if ( calc_qc_v3 )
    {
        qc_v3->Finish();
        qc_v3->WriteHistograms("outputCumulants_v3.root");
    }

    if ( calc_lyz )
    {
        lyz->Finish();
        lyz->WriteHistograms("outputLYZ_100k_kaon_boltzmanPt_Schwinger_WITH_SMEARING.root");
    }

    if ( calc_sp )
    {
        sp->Finish();
//        TFile* outFileSP = new TFile( "outputSP_12k_kaons_boltzmanPt_strIntRad1.root", "RECREATE" );
//        TFile* outFileSP = new TFile( "outputSP_11k_kaons_boltzmanPt_tryEP.root", "RECREATE" );
//        TFile* outFileSP = new TFile( "outputSP_10k_protons_boltzmanPt_LOW_T.root", "RECREATE" );
//        TFile* outFileSP = new TFile( "outputSP_6k_protons_boltzmanPt_Schwinger.root", "RECREATE" );
//        TFile* outFileSP = new TFile( "outputSP_100k_kaon_boltzmanPt_Schwinger_WITH_SMEARING.root", "RECREATE" );
//        TFile* outFileSP = new TFile( "outputSP_TEST_10k_kaon_boltzmanPt_Schwinger_WITH_SMEARING.root", "RECREATE" );
//        TFile* outFileSP = new TFile( Form( "outputSP_TEST_10k_pid%d_lowPtFromString.root", fPID) , "RECREATE" );
//        TFile* outFileSP = new TFile( Form( "outputSP_TEST_10k_pid%d_GaussianMeanPtFromString035.root", fPID) , "RECREATE" );
//        TFile* outFileSP = new TFile( Form( "outputSP_TEST_10k_pid%d_GaussianMeanPtFromString_NEW_PIDS_try6_ExpPt.root", fPID) , "RECREATE" );
//        TFile* outFileSP = new TFile( "outputSP_TEST_10k_TOY.root", "RECREATE" );
//        TFile* outFileSP = new TFile( Form( "outputSP_TEST_11k_pid%d_try13_TSALLIS_r2fm.root", fPID) , "RECREATE" );
        TFile* outFileSP = new TFile( Form( "outputSP_100k_pid%d_PP_TRY1_TSALLIS_r2fm.root", fPID) , "RECREATE" );

        sp->WriteHistograms(outFileSP);


        //QA CHECK:
        TList *list = (TList*)outFileSP->Get("cobjSP");
        AliFlowCommonHistResults *commonHistRes =
                        dynamic_cast<AliFlowCommonHistResults*> (list->FindObject("AliFlowCommonHistResults_SP") );
        TH1D *histSP = commonHistRes->GetHistDiffFlowPtPOI(); //GetHistDiffFlowPtRP();
        histSP->DrawCopy();



        outFileSP->Close();
    }



}


