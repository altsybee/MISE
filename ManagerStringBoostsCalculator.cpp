#include "ManagerStringBoostsCalculator.h"
#include "StringGeneration/StringBoosting.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TMath.h"

#include "TRandom3.h"

#include <fstream>
#include <iostream>
using namespace std;

inline void FixAngleInTwoPi( float &lPhi )
{
    if ( lPhi > 2 * TMath::Pi() )
        lPhi -= 2 * TMath::Pi();
    else if ( lPhi < 0 )
        lPhi += 2 * TMath::Pi();
}

ManagerStringBoostsCalculator::ManagerStringBoostsCalculator() //:
//    fPtrNuclStruct(0x0)
  //  , fFillEventTree(false)
  //  , fNumberOfCentralityBins(10)
  //    fFlagGenerateCentralEventByHand(false)
  //  , fFlagGenerateSemicentralEventByHand(false)
  //  , fImpactParameterByHand(-1)
{
    fOutputDirName = "outputs_ManagerStringBoostsCalculator";

//    fOutputFileName = "testOutput.root";
//    fOutputFile = 0x0;
    fDrawHistos = true;
}

void ManagerStringBoostsCalculator::initOutputObjects()
{
//    int fLowMultHor = 0;
//    int fHiMultHor = 150;//2000;
//    int fMultBinsHor = fHiMultHor - fLowMultHor;
//    double lowMultHor = fLowMultHor - 0.5;
//    double hiMultHor = fHiMultHor - 0.5;

//    fHistSources = new TH1D("fHistSources", "N sources", fMultBinsHor, lowMultHor, hiMultHor);
//    fHistSources->GetXaxis()->SetTitle("N_{sources}");
//    fHistSources->GetYaxis()->SetTitle("N_{events}");
//    fHistSources->SetMarkerStyle(kFullCircle);
}

ManagerStringBoostsCalculator::ManagerStringBoostsCalculator(const ManagerStringBoostsCalculator& ) {}

ManagerStringBoostsCalculator::~ManagerStringBoostsCalculator() {}


void ManagerStringBoostsCalculator::generateBoosts( StringBoosting *booster, int nEvents )
{
    cout << "generating boosts for strings..." << endl;

    const int NMaxStrings = 5000;

    // ##### get NucleiCollision info from file
    TFile* inputFileNuclColTree = new TFile( Form( "%s/%s"
                   , fOutputDirName.Data() //"outputs_NucleiCollision"
                   , fInputFileName.Data() ) );

    if ( !inputFileNuclColTree )
    {
        cout << "NO INPUT FILE!" << endl;
        return;
    }

    TTree *fNucleiCollisionsTree = (TTree*) inputFileNuclColTree->Get("NucleiCollisionsTree");
    Int_t fNuclTreeNumberOfStrings = 0;
    Float_t fNuclTreeStringRadiusVectorAngle[NMaxStrings];
    Float_t fNuclTreeStringX[NMaxStrings];
    Float_t fNuclTreeStringY[NMaxStrings];

    Float_t fDistanceBetweenPartonsForString[NMaxStrings];
    short fStringOrigin[NMaxStrings];

    fNucleiCollisionsTree->Print();


    fNucleiCollisionsTree->SetBranchAddress("numberOfStrings",&fNuclTreeNumberOfStrings );

    fNucleiCollisionsTree->SetBranchAddress("stringRadiusVectorAngle", fNuclTreeStringRadiusVectorAngle );
    fNucleiCollisionsTree->SetBranchAddress("stringX", fNuclTreeStringX );
    fNucleiCollisionsTree->SetBranchAddress("stringY", fNuclTreeStringY );

    fNucleiCollisionsTree->SetBranchAddress( "distanceBetweenPartonsForString", fDistanceBetweenPartonsForString );
    fNucleiCollisionsTree->SetBranchAddress( "stringOrigin", fStringOrigin );


    fOutputFileName = Form( "%s/%s_StringBoosts.root", fOutputDirName.Data(), fInputFileName.Data() );
    TFile* outFile = new TFile( fOutputFileName, "RECREATE" );

    // ##### string boosts tree
    TTree *fStringBoostsTree = new TTree("StringBoostsTree", "StringBoostsTree");


    Float_t fNuclTreeStringBoostAngle[NMaxStrings];
    Float_t fNuclTreeStringBoostMagn[NMaxStrings];
    short fNuclTreeIsHardInteractionString[NMaxStrings]; // was Bool_t

    fStringBoostsTree->Branch("numberOfStrings",    &fNuclTreeNumberOfStrings,"fNuclTreeNumberOfStrings/I");
    fStringBoostsTree->Branch("stringBoostAngle",   fNuclTreeStringBoostAngle,"fNuclTreeStringBoostAngle[fNuclTreeNumberOfStrings]/F");
    fStringBoostsTree->Branch("stringBoostMagn",    fNuclTreeStringBoostMagn,"fNuclTreeStringBoostMagn[fNuclTreeNumberOfStrings]/F");


    fStringBoostsTree->Branch("isHardInteraction",fNuclTreeIsHardInteractionString,"fNuclTreeIsHardInteractionString[fNuclTreeNumberOfStrings]/S"); //O");


    // QA counters
    int countNstr = 0;
    int countNstr_boostZero = 0;
    int countNstr_angleZero = 0;
    int countNstr_isHardInt = 0;

    // check number of events!
    if ( nEvents < 0 || nEvents > fNucleiCollisionsTree->GetEntries() )
        nEvents = fNucleiCollisionsTree->GetEntries();
    cout << "nEvents=" << nEvents << endl;

    // ##### event loop
    for ( int iEvent = 0; iEvent < nEvents; iEvent++)
    {
        if ( iEvent % 100 == 0 )
            cout <<"generating boosts for event " << (int)iEvent << "\r"; cout.flush();

        fNucleiCollisionsTree->GetEntry(iEvent);

//        cout << "check: fNuclTreeNumberOfStrings = " << fNuclTreeNumberOfStrings << endl;
//        cout << "check: fNuclTreeStringX[0] = " << fNuclTreeStringX[0] << endl;

        //process event
        booster->setCoordinatesForStrings( fNuclTreeNumberOfStrings, fNuclTreeStringX, fNuclTreeStringY
                                           , fDistanceBetweenPartonsForString, fStringOrigin );
        booster->processEvent();

        // ##### loop over strings to write their boosts
        for ( int iString = 0; iString < fNuclTreeNumberOfStrings; iString++)
        {
            countNstr++;
            if ( booster->getStringBoostMagn(iString) < 0.00001)
                countNstr_boostZero++;
            if ( booster->getStringBoostAngle(iString) < 0.00001)
                countNstr_angleZero++;
            if ( booster->isHardInteractionString(iString) == 1 )
                countNstr_isHardInt++;

            if(0)if ( booster->getStringBoostMagn(iString) < 0.00001)
                    cout << "   check: iString = " << iString << ", boost mag = " << booster->getStringBoostMagn(iString)
                     << " " << booster->isHardInteractionString(iString) << endl;
            fNuclTreeStringBoostAngle[iString] = booster->getStringBoostAngle(iString);
            fNuclTreeStringBoostMagn[iString] = booster->getStringBoostMagn(iString);
            fNuclTreeIsHardInteractionString[iString] = booster->isHardInteractionString(iString);
//            if (fNuclTreeIsHardInteractionString[iString]==3)cout << "CHECK: " << fNuclTreeIsHardInteractionString[iString] << endl;
        }

        //need some event checks/triggers here?..


        if ( fNuclTreeNumberOfStrings >= NMaxStrings )
        {
            cout << "AHTUNG: fNuclTreeNumberOfStrings >= NMaxStrings!!!" << endl;
            return;
        }

        fStringBoostsTree->Fill();
    } // end of the event loop

    // QA print:
    cout << ">>> Nstr_boostZero=" << countNstr_boostZero << ", angleZero=" << countNstr_angleZero
         << ", isHardInt=" << countNstr_isHardInt << endl;
    cout << ">>> fraction of hardInt=" << (double)countNstr_isHardInt/countNstr << endl;


    //draw the last generated event
    if (fDrawHistos)
    {
        booster->drawStatisticHists();
        booster->drawEventStructure();
//        drawStatHists();
    }
//    booster->finalActions();

//    TCanvas *canv_impParVSnStrings = new TCanvas("canv_impParVSnStrings","canvas impPar VS nStrings",380,280,800,600);
//    fHist2D_impParVSnStrings->DrawCopy("colz");

    //write tree
    outFile->cd();
    fStringBoostsTree->Write();
    outFile->Close();

    inputFileNuclColTree->Close();

}


