#include "/Users/macbook/alice/simpleAnalysis/commonTools/Tools.cxx"  //  "../../commonTools/Tools.cxx"
//#include "../../commonTools/SimpleTrack.cxx"


#include "ManagerStringFragmentation.h"
#include "StringGeneration/NucleiCollision.h"
#include "StringDecayer/StringDescr.h"
#include "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AliSimpleEvent.h"    // "../../simpleEventAnalyzer/AliSimpleEvent.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"

#include "TRandom3.h"

#include <fstream>
#include <iostream>
using namespace std;

const float kEtaBoundForCountingParticles = 0.8;


const int nMaxPartForHist = 4000;
const int nBinsForNparticlesHist = 4000;//500;//4000;//500;
//const int nMaxPartForHist = 400;
//const int nBinsForNparticlesHist = 400;

const double kPtCutsForEtaDistr[5] =
{
    0.1,
    0.5,
    1.0,
    2.0,
    3.0
};


inline void FixAngleInTwoPi( float &lPhi )
{
    if ( lPhi > 2 * TMath::Pi() )
        lPhi -= 2 * TMath::Pi();
    else if ( lPhi < 0 )
        lPhi += 2 * TMath::Pi();
}

//int whichCentralityClass( int mult, const int nq, double *yq)
//{
//    for ( Int_t i = 0; i < nq; i++ )
//        if ( mult < yq[i] )
//        {
//            return i;
//        }
//    return -1;
//}

void getQuantiles(TH1D *h, const int nq, double *yq)
{
    //    const Int_t nq = 10;
    //    const Int_t nshots = 10;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    //    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++)
        xq[i] = Float_t(i+1)/nq;

    //    TGraph *gr70 = new TGraph(nshots);
    //   TH1F *h = new TH1F("h","demo quantiles",50,-3,3);

    h->GetQuantiles(nq,yq,xq);
    //    for (Int_t i=0;i<nq;i++)
    //        cout << yq[i] << ", ";
    //    cout << endl;

    //    TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,600,900);
    //    TGraph *gr = new TGraph(nq,xq,yq);
    //    gr->SetTitle("final quantiles");
    //    gr->SetMarkerStyle(21);
    //    gr->SetMarkerColor(kRed);
    //    gr->SetMarkerSize(0.3);
    //    gr->Draw("ap");
}


ManagerStringFragmentation::ManagerStringFragmentation():
    //    fPtrNuclStruct(0x0)
    fFillEventTree(false)
  , fNumberOfCentralityBins(10)
  , fCutMinNumberOfParticles(1)
  , whatToDoWithHardScattering(0)
  //    fFlagGenerateCentralEventByHand(false)
  //  , fFlagGenerateSemicentralEventByHand(false)
  //  , fImpactParameterByHand(-1)
{
    fOutputDirName = "outputs_ManagerStringFragmentation";

    fInputFileName_NucleiCollisions = "input.root";
    fInputFileName_StringBoosts = "input.root";
    fOutputFileName = "testOutput.root";
    fOutputFile = 0x0;
    fDrawHistos = true;
}

void ManagerStringFragmentation::initOutputObjects()
{
    int fLowMultHor = 0;
    int fHiMultHor = 150;//2000;
    int fMultBinsHor = fHiMultHor - fLowMultHor;
    double lowMultHor = fLowMultHor - 0.5;
    double hiMultHor = fHiMultHor - 0.5;

    fHistSources = new TH1D("fHistSources", "N sources", fMultBinsHor, lowMultHor, hiMultHor);
    fHistSources->GetXaxis()->SetTitle("N_{sources}");
    fHistSources->GetYaxis()->SetTitle("N_{events}");
    fHistSources->SetMarkerStyle(kFullCircle);

    fHistParticlesInSource = new TH1D("fHistParticlesInSource", "N particles in sources", fMultBinsHor, lowMultHor, hiMultHor);
    fHistParticlesInSource->GetXaxis()->SetTitle("N_{particles}");
    fHistParticlesInSource->GetYaxis()->SetTitle("dN_{str}/dN_{part}" );
    fHistParticlesInSource->SetMarkerStyle(kFullCircle);
    
    fHistParticlesInSourceIncludingJets = new TH1D("fHistParticlesInSourceIncludingJets", "N particles in sources including decays and jets", fMultBinsHor, lowMultHor, hiMultHor);
    fHistParticlesInSourceIncludingJets->GetXaxis()->SetTitle("N_{particles}");
    fHistParticlesInSourceIncludingJets->GetYaxis()->SetTitle("dN_{str}/dN_{part}" );
    fHistParticlesInSourceIncludingJets->SetMarkerStyle(kFullCircle);

    fHistParticlesInJets = new TH1D("fHistParticlesInJets", "N particles in jets and decays", 20, -0.5, 19.5);
    fHistParticlesInJets->GetXaxis()->SetTitle("N_{particles}");
    fHistParticlesInJets->GetYaxis()->SetTitle("dN/dN_{part}" );
    fHistParticlesInJets->SetMarkerStyle(kFullCircle);

    fHistParticlesInEvent = new TH1D("fHistParticlesInEvent", "N particles in Event", nBinsForNparticlesHist, -0.5, nMaxPartForHist-0.5); //fMultBinsHor, lowMultHor, hiMultHor);
    fHistParticlesInEvent->GetXaxis()->SetTitle("N_{particles}");
    fHistParticlesInEvent->GetYaxis()->SetTitle("dN/dN_{part}" );
    fHistParticlesInEvent->SetMarkerStyle(kFullCircle);
    
    fHistParticlesInCutConditionInEvent = new TH1D("fHistParticlesInCutConditionInEvent", "N particles in cut conditions in Event", nBinsForNparticlesHist, -0.5, nMaxPartForHist-0.5);
    fHistParticlesInCutConditionInEvent->GetXaxis()->SetTitle("N_{tracks}");
    fHistParticlesInCutConditionInEvent->GetYaxis()->SetTitle("dN/dN_{tracks}" );
    fHistParticlesInCutConditionInEvent->SetMarkerStyle(kFullCircle);

    //    fHistParticlesInEventInEta = new TH1D("fHistParticlesInEventInEta", Form( "N particles in Event in |#eta|<%.2f", kEtaBoundForCountingParticles), fMultBinsHor, lowMultHor, hiMultHor);
    //    fHistParticlesInEventInEta->GetXaxis()->SetTitle("N_{particles}");
    //    fHistParticlesInEventInEta->GetYaxis()->SetTitle("dN/dN_{part}" );
    //    fHistParticlesInEventInEta->SetMarkerStyle(kFullCircle);

    //    fHistParticlesInCutConditionInEventInEta = new TH1D("fHistParticlesInCutConditionInEventInEta", Form("N particles in cut conditions in Event in |#eta|<%.2f", kEtaBoundForCountingParticles), fMultBinsHor, lowMultHor, hiMultHor);
    //    fHistParticlesInCutConditionInEventInEta->GetXaxis()->SetTitle("N_{tracks}");
    //    fHistParticlesInCutConditionInEventInEta->GetYaxis()->SetTitle("dN/dN_{tracks}" );
    //    fHistParticlesInCutConditionInEventInEta->SetMarkerStyle(kFullCircle);


    fHistParticlesInCutConditionVsNu = new TH2D("fHistParticlesInCutConditionVsNu", "N ch Vs #nu"
                                                , nBinsForNparticlesHist, -0.5, nMaxPartForHist-0.5, 80, 1 , 7 );
    fHistParticlesInCutConditionVsNu->GetXaxis()->SetTitle( "N_{tracks}");
    fHistParticlesInCutConditionVsNu->GetYaxis()->SetTitle( "#nu" );
    //    fHistParticlesInCutConditionVsNu->SetMarkerStyle(kFullCircle);



    fHistPt = new TH1D("fHistPt", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);
    fHistPtAfterCuts = new TH1D("fHistPtAfterCuts", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);
    fHistPtBeforeKick = new TH1D("fHistPtBeforeKick", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);

    fHistPtAfterCutsPID[0] = new TH1D("fHistPtAfterCutsPID_pions", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);
    fHistPtAfterCutsPID[1] = new TH1D("fHistPtAfterCutsPID_kaons", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);
    fHistPtAfterCutsPID[2] = new TH1D("fHistPtAfterCutsPID_protons", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);
    fHistPtAfterCutsPID[3] = new TH1D("fHistPtAfterCutsPID_phis", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);
    fHistPtAfterCutsPID[4] = new TH1D("fHistPtAfterCutsPID_lambdas", "p_{T} distribution;p_{T} (GeV/c);1/2#pi p_{T} dN/dp_{T} (GeV/c)^{-2}", 1000, 0.0, 20);


    fHistNeventsInCentralityClasses = new TH1D("fHistNeventsInCentralityClasses", "N particles in class;class;n events", fNumberOfCentralityBins, -0.5, fNumberOfCentralityBins-0.5 );
    //    fHistNeventsInCentralityClasses->SetMarkerStyle(kFullCircle);

    fHistEta = new TH1D("fEta", "#eta distribution;#eta;dN/d#eta", 100, -4, 4);
    for ( int iPtBin = 0; iPtBin < 5; iPtBin++)
        fHistEtaInPtCuts[iPtBin] = new TH1D( Form("fHistEtaInPtCuts%d", iPtBin), "#eta distribution;#eta;dN/d#eta", 100, -4, 4);

    fHistPhi = new TH1D("fHistPhi", "#phi distribution", 200, -0.5, 7 );
    fHistPhi->GetXaxis()->SetTitle("#phi");
    fHistPhi->GetYaxis()->SetTitle("dN/d#phi");
    fHistPhi->SetMarkerStyle(kFullCircle);
    
    fHist2DEtaPhi = new TH2D("fHist2DEtaPhi", "#Delta#eta#Delta#phi b/n child and parent for "
                             , 31, -2, 2, 31, 0, 2*TMath::Pi() );
    fHist2DEtaPhi->GetXaxis()->SetTitle("#eta");
    fHist2DEtaPhi->GetYaxis()->SetTitle("#phi");
    fHist2DEtaPhi->GetZaxis()->SetTitle("dN/(d#eta,d#phi");


    //    fPtCutMin = 0.;
    //    fPtCutMax = 100.;
}

ManagerStringFragmentation::ManagerStringFragmentation(const ManagerStringFragmentation& ) {}

ManagerStringFragmentation::~ManagerStringFragmentation() {}


void ManagerStringFragmentation::applyFragmentationToEvents(StringDescr *strDescr, int nEvents )
{
    cout << "fragmentation of the strings: start looping over events..." << endl;
    
    //reset all hists
    fHistSources->Reset();
    fHistParticlesInSource->Reset();
    fHistParticlesInEvent->Reset();
    fHistPt->Reset();
    fHistPtAfterCuts->Reset();
    for ( int i = 0; i < 3; i++)
        fHistPtAfterCutsPID[i]->Reset();
    fHistPtBeforeKick->Reset();
    fHistEta->Reset();
    for ( int iPtBin = 0; iPtBin < 5; iPtBin++)
        fHistEtaInPtCuts[iPtBin]->Reset();
    fHistPhi->Reset();

    //    fPtrNuclStruct = nuclStruct;
    
    //    SimpleTrack *fSimpleTrack = new SimpleTrack;

    //Form(
    //                  "%s/%s_%devents.root"
    //                , dirOutName.Data(), fileOutName.Data(), nev), "RECREATE" );



    // ##### 10.07.2016 - attempt to save NucleiCollision info into the tree
    TFile* inputFileNuclColTree = 0x0;
    TString strInputFile = Form( "outputs_NucleiCollision/%s", fInputFileName_NucleiCollisions.Data() );
//    TString strInputFile = fInputFileName_NucleiCollisions;//Form( "%s/%s", "outputs_NucleiCollision", fInputFileName.Data() );
    inputFileNuclColTree = new TFile( strInputFile );

    if ( !inputFileNuclColTree )
    {
        cout << "NO INPUT FILE!" << endl;
        return;
    }

    cout << "#### Input files:" << endl;
    cout << fInputFileName_NucleiCollisions << endl;
    cout << fInputFileName_StringBoosts << endl;

    TTree *fNucleiCollisionsTree = (TTree*) inputFileNuclColTree->Get("NucleiCollisionsTree");
    const int NMaxStrings = 5000;

    // check number of events!
    if ( nEvents < 0 || nEvents > fNucleiCollisionsTree->GetEntries() )
        nEvents = fNucleiCollisionsTree->GetEntries();
    cout << "nEvents=" << nEvents << endl;

    Int_t fNuclTreeNumberOfStrings = 0;
    Float_t fNuclTreeStringBoostAngle[NMaxStrings];
    Float_t fNuclTreeStringBoostMagn[NMaxStrings];
    Float_t fNuclTreeStringRadiusVectorAngle[NMaxStrings];
    Bool_t fNuclTreeIsHardInteractionString[NMaxStrings];

    Float_t fNuclTreeRandomEventPlanePhi = 0;
    Float_t fNuclTreeNu = 0;

    fNucleiCollisionsTree->SetBranchAddress( "numberOfStrings", &fNuclTreeNumberOfStrings );

    // ##### Friend tree
//    fNucleiCollisionsTree->AddFriend( "StringBoostsTree", Form( "%s_StringBoosts.root", strInputFile.Data() ) );
    fNucleiCollisionsTree->AddFriend( "StringBoostsTree", Form( "outputs_NucleiCollision/%s", fInputFileName_StringBoosts.Data() ) );

    fNucleiCollisionsTree->SetBranchAddress( "stringBoostAngle", fNuclTreeStringBoostAngle );
    fNucleiCollisionsTree->SetBranchAddress( "stringBoostMagn", fNuclTreeStringBoostMagn );
    fNucleiCollisionsTree->SetBranchAddress( "stringRadiusVectorAngle", fNuclTreeStringRadiusVectorAngle );
    fNucleiCollisionsTree->SetBranchAddress( "isHardInteraction",fNuclTreeIsHardInteractionString );

    fNucleiCollisionsTree->SetBranchAddress( "randomEventPlanePhi", &fNuclTreeRandomEventPlanePhi );
    fNucleiCollisionsTree->SetBranchAddress( "nu", &fNuclTreeNu );

//    return;

    //output file for the events
//    fOutputFileName = Form( "%s/eventTree_nEv%d_try13_TSALLIS_r2fm.root", fOutputDirName.Data(), nEvents );
    fOutputFileName = Form( "%s/%s_StringFragm_nEv%d_TSALLIS_try1.root", fOutputDirName.Data(), fInputFileName_StringBoosts.Data(), nEvents );
            //Form( "%s/eventTree_nEv%d_PP_TRY3_TSALLIS_r2fm.root", fOutputDirName.Data(), nEvents );
    TFile* outFile = new TFile( fOutputFileName, "RECREATE" );


    // ##### 10.07.2016 - attempt to use even-more-simple root tree
    TTree *fEventTrackTree = new TTree("EventTree","EventTree");
    const int NMaxTrack = 15000;
    Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack],fTrackPhi[NMaxTrack]; //
    Int_t fTrackCharge[NMaxTrack];
    Int_t fTrackPID[NMaxTrack];
    Int_t fNumberOfTracks = 0;
    fEventTrackTree->Branch("fNumberOfTracks",&fNumberOfTracks,"fNumberOfTracks/I");
    fEventTrackTree->Branch("fTrackPt",fTrackPt,"fTrackPt[fNumberOfTracks]/F");
    fEventTrackTree->Branch("fTrackPhi",fTrackPhi,"fTrackPhi[fNumberOfTracks]/F");
    fEventTrackTree->Branch("fTrackEta",fTrackEta,"fTrackEta[fNumberOfTracks]/F");
    fEventTrackTree->Branch("fTrackCharge",fTrackCharge,"fTrackCharge[fNumberOfTracks]/I");
    fEventTrackTree->Branch("fTrackPID",fTrackPID,"fTrackPID[fNumberOfTracks]/I");



    const int nPIDs = 5;
    int nParticlesPID[nPIDs];
    for ( int i = 0; i < nPIDs; i++)
        nParticlesPID[i] = 0;

    cout << ">>> nEvents in input tree: " << fNucleiCollisionsTree->GetEntries() << endl;

    // ##### event loop
    for ( int iEvent = 0; iEvent < nEvents; iEvent++)
    {
        if ( iEvent % 1 == 0 )
            cout <<"generating " << (int)iEvent << "\r"; cout.flush();
        //            printf("generating %d event...\n",(int)iEvent );

        fNucleiCollisionsTree->GetEntry(iEvent);

        //need some event checks/triggers here?..

        //do analysis, store events with tracks in a tree
        if (1)
        {
            // ##### get random EP for the event!
            float phiEP = fNuclTreeRandomEventPlanePhi;

            // ##### loop over strings: hadronized and count multiplicity in VZERO
            const double cutEtaWithinDetectorAcceptance = 3;
            int nParticles = 0;
            fNumberOfTracks = 0; // for a new event tree
            int nParticlesWithinEtaPtCuts = 0;

            if ( fNuclTreeNumberOfStrings >= NMaxStrings )
            {
                cout << "AHTUNG: fNuclTreeNumberOfStrings >= NMaxStrings!!!" << endl;
                return;
            }

            for ( int iString = 0; iString < fNuclTreeNumberOfStrings; iString++)
            {
                // TMP CUT!!! July 2016 tests
//                if ( fNuclTreeStringBoostMagn[iString] < 0.5 )
//                    continue;

                if (iString==0) //QA plotting
                    cout << " fNuclTreeStringBoostMagn[iString] = " << fNuclTreeStringBoostMagn[iString] << endl;

                if ( !fNuclTreeIsHardInteractionString[iString] ) // this string is a soft interaction
                    strDescr->hadronizeString( fNuclTreeStringBoostMagn[iString], fNuclTreeStringBoostAngle[iString] );
                else //hard interaction
                {
                    if (whatToDoWithHardScattering==0) // create two jets
                        strDescr->makeTwoJets();
                    else if (whatToDoWithHardScattering==1) // particle pair with random pt from Power law
                        strDescr->makeTwoParticlesWithRandomPtEtaPhi();
                }

                int nParticlesInString = strDescr->getNparticles();

                //for "rescattering by hand" below
                double stringRadiusVectorAngle = fNuclTreeStringRadiusVectorAngle[iString];

                for ( int iP = 0; iP < nParticlesInString; iP++ )
                {
                    ParticleDescr *p = strDescr->getParticle(iP);
                    double phi = p->phi;

                    // 20.05.2016 - try to do "rescattering" of particles which are emitted "inside the medium" (to probably get right mass ordering!)
                    if(0)
                    {
                        //double phi_wrt_plane = phi - stringRadiusVectorAngle;
                        //FixAngleInTwoPi(phi_wrt_plane);

                        float phi_wrt_plane = stringRadiusVectorAngle-phi;
                        FixAngleInTwoPi(phi_wrt_plane);
                        const float extra = TMath::Pi()/8; //spec angle to constain even more
                        if ( phi_wrt_plane > TMath::PiOver2()-extra && phi_wrt_plane < 3*TMath::PiOver2()+extra )
                            phi = gRandom->Uniform( 0, 2*TMath::Pi() );
                        // continue;
                    }

                    //rotate all tracks in phi by event plane phi-angle
                    float phiMod = phiEP + phi;
                    FixAngleInTwoPi(phiMod);

                    //spec counting to compare with publications
                    //                    if ( fabs(p->eta) < 0.3 && ( p->pt > 0.15 && p->pt < 10 ) )
                    //                    if ( p->eta < -3.7 || p->eta > 1. )//fabs(p->eta) > 3.7 ) // out of our intresting range
                    //                        continue;
                    double weightForPt = 1;
                    //                    double weightForPt = 1./ ( 2*TMath::Pi()*p->pt );
                    fHistPt->Fill(p->pt, weightForPt );

                    if ( fabs(p->eta) > cutEtaWithinDetectorAcceptance )//fabs(p->eta) > 3.7 ) // out of our intresting range
                        continue;

                    //                    if ( p->eta > -3.7 && p->eta < -1.7 )  //VZERO-C
                    // if ( fabs(p->eta) < 1. )  //STAR TPC

                    //                    if ( fabs(p->eta) < 0.5 )


                    if ( fabs(p->eta) < cutEtaWithinDetectorAcceptance )//0.5 )//&& ( p->pt > 0. && p->pt < 1000 ) )
                    {
                        nParticlesWithinEtaPtCuts++;

                        //                        if ( p->pt > 0.05 )
                        fHistPtAfterCuts->Fill(p->pt, weightForPt );
                        if ( p->pid == 0 ) { fHistPtAfterCutsPID[0]->Fill(p->pt, weightForPt ); nParticlesPID[0]++; }
                        if ( p->pid == 1 ) { fHistPtAfterCutsPID[1]->Fill(p->pt, weightForPt ); nParticlesPID[1]++; }
                        if ( p->pid == 2 ) { fHistPtAfterCutsPID[2]->Fill(p->pt, weightForPt ); nParticlesPID[2]++; }
                        if ( p->pid == 3 ) { fHistPtAfterCutsPID[3]->Fill(p->pt, weightForPt ); nParticlesPID[3]++; }
                        if ( p->pid == 4 ) { fHistPtAfterCutsPID[4]->Fill(p->pt, weightForPt ); nParticlesPID[4]++; }

                        if ( 0 )//not used... ( p->ptBeforeKick > 0.0 )
                        {
                            double weightForPtBeforeKick = 1./ ( 2*TMath::Pi()*p->ptBeforeKick );
                            fHistPtBeforeKick->Fill(p->ptBeforeKick, weightForPtBeforeKick );
                        }
                    }

                    fHistEta->Fill(p->eta);
                    if ( p->pt > kPtCutsForEtaDistr[0] ) fHistEtaInPtCuts[0]->Fill(p->eta);
                    if ( p->pt > kPtCutsForEtaDistr[1] ) fHistEtaInPtCuts[1]->Fill(p->eta);
                    if ( p->pt > kPtCutsForEtaDistr[2] ) fHistEtaInPtCuts[2]->Fill(p->eta);
                    if ( p->pt > kPtCutsForEtaDistr[3] ) fHistEtaInPtCuts[3]->Fill(p->eta);
                    if ( p->pt > kPtCutsForEtaDistr[4] ) fHistEtaInPtCuts[4]->Fill(p->eta);
                    fHistPhi->Fill(phiMod);
                    fHist2DEtaPhi->Fill( p->eta, phiMod );

                    // ##### add particle to analysers
                    //                    fSimpleTrack->id = nParticles;
                    //                    fSimpleTrack->eta = p->eta;
                    //                    fSimpleTrack->phi = phiMod;
                    //                    fSimpleTrack->pt = p->pt;
                    //                    fSimpleTrack->pt = p->pid;
                    //                    fSimpleTrack->charge = 0;//charge;

                    //                    cout << "pid = " << p->pid << endl;

                    nParticles++;

                    if ( fFillEventTree )
                    {
                        //fill arrays for NEW tree
                        if ( fNumberOfTracks >= NMaxTrack )
                        {
                            cout << "AHTUNG: fNumberOfTracks >= NMaxTrack for tree arrays!!!" << endl;
                            return;
                        }
                        fTrackPt[fNumberOfTracks]  = p->pt;
                        fTrackPhi[fNumberOfTracks] = phiMod;
                        fTrackEta[fNumberOfTracks] = p->eta;
                        fTrackCharge[fNumberOfTracks] = p->charge;
                        fTrackPID[fNumberOfTracks] = p->pid;

                        fNumberOfTracks++;
                    }
                } // end of loop over particles from string
            }
            fHistParticlesInEvent->Fill(nParticles);

            fHistParticlesInCutConditionInEvent->Fill( (double)nParticlesWithinEtaPtCuts/(2*cutEtaWithinDetectorAcceptance ) ); // !
            fHistParticlesInCutConditionVsNu->Fill( (double)nParticlesWithinEtaPtCuts/(2*cutEtaWithinDetectorAcceptance ),
                                                    fNuclTreeNu );

            // cut on nTracks!!! July 2016
            if ( nParticles < fCutMinNumberOfParticles )
                continue;

            if ( fFillEventTree )
            {
                //Aug15 (IA): why was /2 ?..
                //                fSimpleEvent->GetHeader()->SetCentrality(nParticlesWithinEtaPtCuts/2/cutEtaWithinDetectorAcceptance);
                //                fSimpleEvent->GetHeader()->SetCentrality(nParticlesWithinEtaPtCuts/cutEtaWithinDetectorAcceptance);
                fEventTrackTree->Fill();
            }
            //fNucleiCollisionsTree->Fill();
        }

    } // end of the event loop

    // ##### print particle ratios
    int nParticlesPID_All = 0;
    for ( int i = 0; i < nPIDs; i++)
        nParticlesPID_All += nParticlesPID[i];
    cout << ">> pions/all = " << (double)nParticlesPID[0]/nParticlesPID_All<< endl;
    cout << ">> kaons/all = " << (double)nParticlesPID[1]/nParticlesPID_All<< endl;
    cout << ">> protons/all = " << (double)nParticlesPID[2]/nParticlesPID_All<< endl;
    cout << ">> phi-s/all = " << (double)nParticlesPID[3]/nParticlesPID_All<< endl;
    cout << ">> Lambda-s/all = " << (double)nParticlesPID[4]/nParticlesPID_All<< endl;


    //fill centralities
    //    Float_t         fEvtHdr_fCentrality;
    //    TBranch        *b_event_fEvtHdr_fCentrality;   //!
    //    fEventTree->SetBranchAddress("fEvtHdr.fCentrality", &fEvtHdr_fCentrality, &b_event_fEvtHdr_fCentrality);
    //    for ( int iEv = 0; iEv < fEventTree->GetEntries(); iEv++ )
    //    {
    //        fEventTree->GetEntry ( iEv ) ;
    //        cout << fSimpleEvent->fVertexX << endl;
    //        fSimpleEvent->SetCentrality( iEv );
    //        fEventTree->Fill();  //fill the tree
    //    }

    //draw the last generated event
    if (fDrawHistos)
    {
        //        fPtrNuclStruct->drawEventStructure();

        //        fPtrNuclStruct->drawStatisticHists();
        drawStatHists();
    }
    //    fPtrNuclStruct->finalActions();


    //    delete d;
    //    delete strDescr;

    //    fQAhistos.fHistNch->SetName("hist70-80");
    //    fQAhistos.fHistNch->Write();
    //    f->Close();


    // ##### multiplicity quantiles
    //prepare quantiles
    Double_t *centrMultBounds = new Double_t[fNumberOfCentralityBins]; // array to contain the quantiles

    //    bool flagHaveQuantiles = false; //flag that we have calculated multiplicity quantiles

    TCanvas *fCanvMultClasses = new TCanvas("canvMultClasses","multiplicity classes",50,50,800,600);
    fCanvMultClasses->cd()->SetLogy();

    fHistParticlesInCutConditionInEvent->SetLineColor(kRed);
    fHistParticlesInCutConditionInEvent->DrawCopy();
    fCanvMultClasses->SaveAs(Form( "%s/histParticlesInCutCondition.root", fOutputDirName.Data() ) );
    fCanvMultClasses->SaveAs(Form( "%s/histParticlesInCutCondition.png", fOutputDirName.Data() ) );
    fCanvMultClasses->SaveAs(Form( "%s/histParticlesInCutCondition.eps", fOutputDirName.Data() ) );
    fCanvMultClasses->SaveAs(Form( "%s/histParticlesInCutCondition.C", fOutputDirName.Data() ) );
    //    fHistParticlesInCutConditionInEvent->GetQuantiles();

    cout << ">>> Mean Nch overall = " << fHistParticlesInCutConditionInEvent->GetMean() << endl;

    //remove zero bin:
    if ( fHistParticlesInCutConditionInEvent->GetBinContent(1) > 0 )
        fHistParticlesInCutConditionInEvent->SetBinContent(1, 0);
    // get quantiles
    getQuantiles( fHistParticlesInCutConditionInEvent, fNumberOfCentralityBins, centrMultBounds );
    //    flagHaveQuantiles = true;

    //fill multiplicity boundaries for each centrality class
    TH1D *fHistMultClassBoundaries = new TH1D( "fHistMultiplicityBoundaries", "fHistMultClassBoundaries"
                                               , fNumberOfCentralityBins, -0.5, fNumberOfCentralityBins-0.5 );
    TH1D *fHistMultClassMeanNch = new TH1D( "fHistMultClassMeanNch", "fHistMultClassMeanNch"
                                            , fNumberOfCentralityBins, -0.5, fNumberOfCentralityBins-0.5 );
    for ( int iBin = 0; iBin < fNumberOfCentralityBins; iBin++ )
    {
        fHistMultClassBoundaries->SetBinContent( iBin+1, centrMultBounds[iBin] );
    }

    //draw centrality classes on mult hist
    //    int nBinsForClassHists = fHistParticlesInCutConditionInEvent->GetNbinsX();
    TH1D **fHistCentrClass = new TH1D*[fNumberOfCentralityBins];
    for ( int iCentrClass = 0; iCentrClass < fNumberOfCentralityBins; iCentrClass++ )
        fHistCentrClass[iCentrClass] = new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass), nBinsForNparticlesHist, -0.5, nMaxPartForHist-0.5);

    int centrClassId = 0;
    for ( int iBin = 0; iBin < fHistParticlesInCutConditionInEvent->GetNbinsX(); iBin++ )
    {
        //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
        if ( fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
            centrClassId++;
        if ( centrClassId >= fNumberOfCentralityBins )
            break;
        double binContent = fHistParticlesInCutConditionInEvent->GetBinContent( iBin+1 );
        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    }
    for ( int iCentrClass = 0; iCentrClass < fNumberOfCentralityBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
        fHistCentrClass[iCentrClass]->DrawCopy( "same" );

        double meanNchInCentrBin = fHistCentrClass[iCentrClass]->GetMean();
        fHistMultClassMeanNch->SetBinContent( iCentrClass+1, meanNchInCentrBin );
        cout << "meanNch in centrality bin " << iCentrClass << ": " << meanNchInCentrBin << endl;
    }


    //write objects
    outFile->cd();
    fEventTrackTree->Write();
    fHistMultClassBoundaries->Write();
    fHistMultClassMeanNch->Write();
    fHistParticlesInCutConditionVsNu->Write();
    outFile->Close();

    inputFileNuclColTree->Close();


    //    delete outFile;

}

void ManagerStringFragmentation::cleanup()
{
    fOutputFile->Close();
    fOutputFile = 0x0;
}

void ManagerStringFragmentation::drawStatHists()
{
    fCanv = new TCanvas("EvManagerStats","Event Manager Statistics",280,250,800,600);
    fCanv->Divide(3,2);
    //    fCanv->cd(1);
    //    fHistSources->DrawCopy();
    //    fCanv->cd(2);
    //    fHistParticlesInSource->DrawCopy();
    fCanv->cd(1)->SetLogy();
    fHistParticlesInEvent->DrawCopy();
    //    fHistParticlesInCutConditionInEvent->SetLineColor(kRed);
    //    fHistParticlesInCutConditionInEvent->DrawCopy("same");
    //    //    fHistParticlesInCutConditionInEvent->GetQuantiles();
    //    //quantiles
    //    const Int_t fNumberOfCentralityBins = 10;
    //    Double_t centrMultBounds[fNumberOfCentralityBins];  // array to contain the quantiles
    //    getQuantiles( fHistParticlesInCutConditionInEvent, fNumberOfCentralityBins, centrMultBounds );
    //    //draw centrality classes on mult hist
    //    TH1D *fHistCentrClass[fNumberOfCentralityBins];
    //    for ( int iCentrClass = 0; iCentrClass < fNumberOfCentralityBins; iCentrClass++ )
    //        fHistCentrClass[iCentrClass] = new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass), fHistParticlesInCutConditionInEvent->GetNbinsX(), 0, 8000);

    //    int centrClassId = 0;
    //    for ( int iBin = 0; iBin < fHistParticlesInCutConditionInEvent->GetNbinsX(); iBin++ )
    //    {
    //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
    //        if ( fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
    //            centrClassId++;
    //        if ( centrClassId >= fNumberOfCentralityBins )
    //            break;
    //        double binContent = fHistParticlesInCutConditionInEvent->GetBinContent( iBin+1 );
    //        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    //    }
    //    for ( int iCentrClass = 0; iCentrClass < fNumberOfCentralityBins; iCentrClass++ )
    //    {
    //        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
    //        fHistCentrClass[iCentrClass]->DrawCopy( "same" );
    //    }



    fCanv->cd(2)->SetLogy();
    fHistPtAfterCuts->SetLineColor(kRed);
    fHistPtBeforeKick->SetLineColor(kGreen);
    //    fHistPt->DrawNormalized();
    //    fHistPtAfterCuts->DrawNormalized("same");
    //    fHistPtBeforeKick->DrawNormalized("same");
    fHistPt->DrawNormalized();
    fHistPtAfterCuts->DrawNormalized("same");
    //    fHistPtBeforeKick->Draw("same");

    fHistPtAfterCutsPID[0]->SetLineColor(kGreen);
    fHistPtAfterCutsPID[1]->SetLineColor(kOrange+1);
    fHistPtAfterCutsPID[2]->SetLineColor(kMagenta);
    fHistPtAfterCutsPID[0]->DrawNormalized("same");
    fHistPtAfterCutsPID[1]->DrawNormalized("same");
    fHistPtAfterCutsPID[2]->DrawNormalized("same");

    cout << ">> <pT> pions = " <<       fHistPtAfterCutsPID[0]->GetMean()<< endl;
    cout << ">> <pT> kaons = " <<       fHistPtAfterCutsPID[1]->GetMean()<< endl;
    cout << ">> <pT> protons = " <<     fHistPtAfterCutsPID[2]->GetMean()<< endl;
    cout << ">> <pT> phi-s = " <<     fHistPtAfterCutsPID[3]->GetMean()<< endl;
    cout << ">> <pT> Lambda-s = " <<     fHistPtAfterCutsPID[4]->GetMean()<< endl;


    TLegend *legPt = new TLegend(0.65,0.5,0.95,0.8);
    legPt->SetFillStyle(0);
    legPt->SetBorderSize(0);
    legPt->AddEntry(fHistPt,  "no #eta cut",  "l");
    legPt->AddEntry(fHistPtAfterCuts,  "after #eta cut",  "l");
    //    legPt->AddEntry(fHistPtBeforeKick,  "after #eta cut, before kick",  "l");
    legPt->AddEntry(fHistPtAfterCutsPID[0],  "pions in #eta cut",  "l");
    legPt->AddEntry(fHistPtAfterCutsPID[1],  "kaons in #eta cut",  "l");
    legPt->AddEntry(fHistPtAfterCutsPID[2],  "protons in #eta cut",  "l");
    legPt->Draw();


    fCanv->cd(3);
    TLegend *legEta = new TLegend(0.65,0.5,0.95,0.8);
    legEta->SetFillStyle(0);
    legEta->SetBorderSize(0);
    legEta->AddEntry(fHistEta,  "all p_{T}",  "l");

    fHistEta->DrawNormalized();
    fHistEtaInPtCuts[0]->SetLineColor(kRed);
    fHistEtaInPtCuts[1]->SetLineColor(kGreen);
    fHistEtaInPtCuts[2]->SetLineColor(kMagenta);
    fHistEtaInPtCuts[3]->SetLineColor(kOrange);
    fHistEtaInPtCuts[4]->SetLineColor(kBlack);
    for ( int iPtBin = 0; iPtBin < 3/*5*/; iPtBin++)
    {
        fHistEtaInPtCuts[iPtBin]->DrawNormalized("same");
        legEta->AddEntry(fHistEtaInPtCuts[iPtBin],  Form("p_{T}<%.2f", kPtCutsForEtaDistr[iPtBin]),  "l");
    }
    legEta->Draw();

    fCanv->cd(4);
    fHistPhi->DrawCopy();
    //fCanv->cd(7);
    //fEvent->getHistStringLength()->Draw();
    //    fCanv->cd(8);
    //    if ( fEvent->getStringCreator()->getDealWithJets() )
    //        fHist2DEtaPhi->Draw("surf1");

    //    fCanv->cd(9);
    //    fEvent->getStringCreator()->getHistStringLength()->Draw();

    fCanv->cd(5);
    fHistNeventsInCentralityClasses->DrawCopy();

    fCanv->cd(6);
    fHistParticlesInCutConditionVsNu->DrawCopy("colz");

    fCanv->SaveAs(Form( "%s/canvEvManagerStats.root", fOutputDirName.Data() ) );
    fCanv->SaveAs(Form( "%s/canvEvManagerStats.eps", fOutputDirName.Data() ));
    fCanv->SaveAs(Form( "%s/canvEvManagerStats.png", fOutputDirName.Data() ));
    fCanv->SaveAs(Form( "%s/canvEvManagerStats.C", fOutputDirName.Data() ));


    //write stats to file
    TFile *fileEvManagerStats = new TFile(Form( "%s/stats_EvManager.root", fOutputDirName.Data() ),"RECREATE");
    fHistParticlesInEvent->Write();
    fHistParticlesInCutConditionInEvent->Write();
    fHistPt->Write();
    fHistPtAfterCuts->Write();
    fHistPtAfterCutsPID[0]->Write();
    fHistPtAfterCutsPID[1]->Write();
    fHistPtAfterCutsPID[2]->Write();
    //    fHistPtBeforeKick->Write();
    fHistEta->Write();
    fHistPhi->Write();
    fileEvManagerStats->Close();

    // ###### spec output txt files
    ofstream fout( Form( "%s/tmpTextOutput_ManagerStringFragmentation_ptVsNch.txt", fOutputDirName.Data() ), ios::out | ios::binary);
    fout << fHistParticlesInCutConditionInEvent->GetMean() //fPtrNuclStruct->getImpactParameterByHand() //_0_100
         << " " << fHistPtAfterCuts->GetMean()
         << " " << fHistPtAfterCuts->GetMeanError()
         << endl;
    fout.close();


}
