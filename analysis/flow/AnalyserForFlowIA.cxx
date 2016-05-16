#include "AnalyserForFlowIA.h"
//#include "/Users/macbook/alice/simpleAnalysis/commonTools/SimpleTrack.cxx"
#include "/Users/macbook/alice/simpleAnalysis/commonTools/MiniEvent.cxx"


// aliroot includes
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackSimpleCuts.h"

#include "AliFlowAnalysisWithMCEventPlane.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "AliFlowAnalysisWithLeeYangZeros.h"
#include "AliFlowAnalysisWithScalarProduct.h"



#include "Riostream.h"
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
//#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"

#include "TRandom3.h"


bool calc_qc_v2 = false;
bool calc_qc_v3 = false;//true;
bool calc_lyz = false;
bool calc_sp = true;//true;

// k) Define the ranges for two subevents separated with eta gap (needed only for SP method):
Double_t etaMinA = -3; // minimum eta of subevent A
Double_t etaMaxA = -2; // maximum eta of subevent A
Double_t etaMinB = 2; // minimum eta of subevent B
Double_t etaMaxB = 3; // maximum eta of subevent B



TRandom3 pRandom;


ClassImp(AnalyserForFlowIA)

inline void fixDeltaPhi( Double_t &dphi)
{
    if ( dphi < 0 )
        dphi = dphi + 2*TMath::Pi();
}

AnalyserForFlowIA::AnalyserForFlowIA()
    :
      //    fIsEventOpened( kFALSE )
      //  , fIsOnline( kFALSE )
      //  , fEventCount( 0 )
      fNtracks (0),
      fPID(-1)
{
    fMiniEvent = new MiniEvent;

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
    cutsPOI->SetEtaMin(-1);
    cutsPOI->SetEtaMax(1);

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


}

Bool_t AnalyserForFlowIA::InitDataMembers()
{
    //    fShortDef
    cout << " # Init for  AnalyserForFlowIA \n";
    //    if( fIsOnline )
    //    {
    //        Printf("Can't init data members more then one time! \n");
    //        return kFALSE;
    //    }


    //    fEventCount = 0;
    //    fFlagHavePrevEvent = kFALSE;
    fOutList = new TList();
    fOutList->SetName("flowData");

    fHistNparticlesInTile  = new TH1D( "fHistNparticlesInTile", "N particles in a tile;n;dN/dn", 501, -0.5, 500.5);
    fOutList->Add( fHistNparticlesInTile );

    return kTRUE;
}

AnalyserForFlowIA::~AnalyserForFlowIA()
{
    //Destructor

}

void AnalyserForFlowIA::StartEvent()
{
    // Open new Event for track by track event import
    //    if( fIsEventOpened )                     // Check if trying to open event more than once !
    //    {
    //        FinishEvent();
    //    }

    //    if( !fIsOnline )  // Autocreating histos if InitDataMembers was not called by hand
    //    {
    //        Printf("InitDataMembers was not called by hand ! Autocreating histos...\n");
    //        InitDataMembers();
    //    }

    //    fIsEventOpened = kTRUE;
}


void AnalyserForFlowIA::AddTrack(const SimpleTrack *track) //AddTrack( Double_t Eta , Double_t Phi, Double_t Pt )
{
//    if ( fPID == -1 || fPID == track->pid )
    {
        Double_t Eta = track->eta;
        Double_t Phi = track->phi;
        Double_t Pt = track->pt;
        Double_t Charge = track->charge;
        Double_t Pid = track->pid;

        fHistNparticlesInTile->Fill(5);

        UpdateEvent ( fMiniEvent, Eta, Phi, Pt, Charge, Pid );

        fNtracks++;
    }
}

// !!!!!!!! to use UpdateEvent in MiniEvent!!! (?)
void AnalyserForFlowIA::UpdateEvent(MiniEvent *event, Double_t Eta , Double_t Phi, Double_t Pt , Int_t Charge, Int_t Pid )
{
    SimpleTrack *track = &event->tracks[event->nTracks]; //get current track to update
    track->id = fNtracks;
    track->eta = Eta;
    track->phi = Phi;
    track->pt = Pt;
    track->charge = Charge;
    track->pid = Pid;
    event->nTracks++;
}

void AnalyserForFlowIA::FillAliFlowEvent()
{
    int nTracks = fMiniEvent->nTracks;
    AliFlowEventSimple* flowEvent = new AliFlowEventSimple(nTracks);//,AliFlowEventSimple::kGenerate);
    //    AliFlowEventSimple* flowEventSP = new AliFlowEventSimple(nTracks);//,AliFlowEventSimple::kGenerate);
    //    flowEventSP->TagSubeventsInEta(-2,-1,1,2);

    SimpleTrack *tracks = fMiniEvent->tracks;
    for ( Int_t i = 0; i < nTracks; i++ )
    {
        if ( tracks[i].pid == 2 ) //PID!
            if(0)  cout << "pt = " << tracks[i].pt
                        << ", eta = " << tracks[i].eta
                        << ", phi=" << tracks[i].phi
                        << ", charge=" << tracks[i].charge
                        << ", pid=" << tracks[i].pid
                        << endl;

        AliFlowTrackSimple* flowtrack = new AliFlowTrackSimple();
        flowtrack->SetPhi(tracks[i].phi);
        flowtrack->SetEta(tracks[i].eta);
        flowtrack->SetPt(tracks[i].pt);
        flowtrack->SetCharge(tracks[i].charge);

        if (cutsRP)
            flowtrack->TagRP(cutsRP ->PassesCuts(flowtrack));
        if (cutsPOI)
        {
            bool isPOI = cutsPOI->PassesCuts(flowtrack);
            bool isPID_we_want =  ( fPID == -1 || fPID == tracks[i].pid );

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
    flowEvent->TagSubeventsInEta(etaMinA,etaMaxA,etaMinB,etaMaxB);

    //    mcep->Make(flowEvent);
    if ( calc_qc_v2 )   qc_v2->Make(flowEvent);
    if ( calc_qc_v3 )   qc_v3->Make(flowEvent);
    if ( calc_lyz )     lyz->Make(flowEvent);

    if ( calc_sp )     sp->Make(flowEvent);
    //        cout <<"Event: " << i+1 << "\r"; cout.flush();
    delete flowEvent;
    //    delete flowEventSP;

}



void AnalyserForFlowIA::FinishEvent()
{
    // Track by track event import : Close opened event and fill event summary histos
    //    if( !fIsEventOpened )
    //    {
    //        Printf("Event is not opened!\n");
    //        return;
    //    }
    FillAliFlowEvent();


    //    //refresh info
    fMiniEvent->nTracks = 0;
    fNtracks = 0;
    //    fIsEventOpened = kFALSE;
}

void AnalyserForFlowIA::Terminate()
{
    // calculate the final results
    //    fOutputFileId = -1000;
    //    mcep->Finish();
    if ( calc_qc_v2 )
    {
        qc_v2->Finish();
        qc_v2->WriteHistograms("outputCumulants_v2.root");
    }

    if ( calc_qc_v3 )
    {
        qc_v3->Finish();
        qc_v3->WriteHistograms("outputCumulants_v3.root");
    }

    if ( calc_lyz )
    {
        lyz->Finish();
        lyz->WriteHistograms("outputLYZ.root");
    }

    if ( calc_sp )
    {
        sp->Finish();
        TFile* outFileSP = new TFile( "outputSP.root", "RECREATE" );
        sp->WriteHistograms(outFileSP);
        outFileSP->Close();
    }

}



