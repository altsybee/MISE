#include "ManagerNucleiCollisions.h"
#include "StringGeneration/NucleiCollision.h"

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

ManagerNucleiCollisions::ManagerNucleiCollisions() //:
//    fPtrNuclStruct(0x0)
  //  , fFillEventTree(false)
  //  , fNumberOfCentralityBins(10)
  //    fFlagGenerateCentralEventByHand(false)
  //  , fFlagGenerateSemicentralEventByHand(false)
  //  , fImpactParameterByHand(-1)
{
    fOutputDirName = "outputs_ManagerNucleiCollisions";

//    fOutputFileName = "testOutput.root";
//    fOutputFile = 0x0;
    fDrawHistos = true;
}

void ManagerNucleiCollisions::initOutputObjects()
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

ManagerNucleiCollisions::ManagerNucleiCollisions(const ManagerNucleiCollisions& ) {}

ManagerNucleiCollisions::~ManagerNucleiCollisions() {}


void ManagerNucleiCollisions::generateEvents( NucleiCollision *fPtrNuclStruct, int nEvents )
{
    cout << "generating events..." << endl;
    
    //reset all hists
//    fHistSources->Reset();

    TH2D *fHist2D_impParVSnStrings = new TH2D("fHist2D_impParVSnStrings", ";impact parameter, fm;n strings"
                                                , 400, 0, 20, 2001, -0.5 , 2000.5 );



    // ##### 10.07.2016 - save NucleiCollision info into the tree
    TTree *fNucleiCollisionsTree = new TTree("NucleiCollisionsTree","NucleiCollisionsTree");
    const int NMaxStrings = 35000;//8000;

    Float_t fImpactParameter = 0;
    Float_t fNuclTreeRandomEventPlanePhi = 0;
    Float_t fNuclTreeNu = 0;
    Int_t fNuclTreeNumberOfStrings = 0;
    Int_t fMultFromAllStringsFictiveV0 = 0;
    Int_t fMultFromAllStringsFictiveMidEta = 0;
//    Float_t fNuclTreeStringBoostAngle[NMaxStrings];
//    Float_t fNuclTreeStringBoostMagn[NMaxStrings];
    Float_t fDistanceBetweenPartonsForString[NMaxStrings];
    Float_t fNuclTreeStringRadiusVectorAngle[NMaxStrings];
    Float_t fNuclTreeStringX[NMaxStrings];
    Float_t fNuclTreeStringY[NMaxStrings];
    short fStringOrigin[NMaxStrings];
//    Bool_t fNuclTreeIsHardInteractionString[NMaxStrings];

    Int_t fNumberOfWoundedPartonsA = 0;
    Int_t fNumberOfWoundedPartonsB = 0;
    Int_t fNumberOfWoundedOnlyOncePartonsA = 0;
    Int_t fNumberOfWoundedOnlyOncePartonsB = 0;

    Int_t fNumberOfWoundedNucleonsA = 0;
    Int_t fNumberOfWoundedNucleonsB = 0;
    Int_t fNumberOfWoundedOnlyOnceNucleonsA = 0;
    Int_t fNumberOfWoundedOnlyOnceNucleonsB = 0;

    // May 2018: for core-corona:
    const int maxBinsCoreCorona = 23;
    int fStringsInCells[maxBinsCoreCorona][maxBinsCoreCorona];

    fNucleiCollisionsTree->Branch("impactParameter",&fImpactParameter,"fImpactParameter/F");
    fNucleiCollisionsTree->Branch("randomEventPlanePhi",&fNuclTreeRandomEventPlanePhi,"fNuclTreeRandomEventPlanePhi/F");
    fNucleiCollisionsTree->Branch("nu",&fNuclTreeNu,"fNuclTreeNu/F");
    fNucleiCollisionsTree->Branch("numberOfStrings",&fNuclTreeNumberOfStrings,"fNuclTreeNumberOfStrings/I");

    fNucleiCollisionsTree->Branch("nWoundedPartonsA",&fNumberOfWoundedPartonsA,"fNumberOfWoundedPartonsA/I");
    fNucleiCollisionsTree->Branch("nWoundedPartonsB",&fNumberOfWoundedPartonsB,"fNumberOfWoundedPartonsB/I");
    fNucleiCollisionsTree->Branch("nWoundedOnlyOncePartonsA",&fNumberOfWoundedOnlyOncePartonsA,"fNumberOfWoundedOnlyOncePartonsA/I");
    fNucleiCollisionsTree->Branch("nWoundedOnlyOncePartonsB",&fNumberOfWoundedOnlyOncePartonsB,"fNumberOfWoundedOnlyOncePartonsB/I");

    fNucleiCollisionsTree->Branch("nWoundedNucleonsA",&fNumberOfWoundedNucleonsA,"fNumberOfWoundedNucleonsA/I");
    fNucleiCollisionsTree->Branch("nWoundedNucleonsB",&fNumberOfWoundedNucleonsB,"fNumberOfWoundedNucleonsB/I");
    fNucleiCollisionsTree->Branch("nWoundedOnlyOnceNucleonsA",&fNumberOfWoundedOnlyOnceNucleonsA,"fNumberOfWoundedOnlyOnceNucleonsA/I");
    fNucleiCollisionsTree->Branch("nWoundedOnlyOnceNucleonsB",&fNumberOfWoundedOnlyOnceNucleonsB,"fNumberOfWoundedOnlyOnceNucleonsB/I");

    if(1)
        fNucleiCollisionsTree->Branch("nStringsInCells", fStringsInCells,"fStringsInCells[23][23]/I");

    if(0)
    {
        //    fNucleiCollisionsTree->Branch("stringBoostAngle", fNuclTreeStringBoostAngle,"fNuclTreeStringBoostAngle[fNuclTreeNumberOfStrings]/F");
        //    fNucleiCollisionsTree->Branch("stringBoostMagn", fNuclTreeStringBoostMagn,"fNuclTreeStringBoostMagn[fNuclTreeNumberOfStrings]/F");
            fNucleiCollisionsTree->Branch("stringX", fNuclTreeStringX,"fNuclTreeStringX[fNuclTreeNumberOfStrings]/F");
            fNucleiCollisionsTree->Branch("stringY", fNuclTreeStringY,"fNuclTreeStringY[fNuclTreeNumberOfStrings]/F");
            fNucleiCollisionsTree->Branch("distanceBetweenPartonsForString", fDistanceBetweenPartonsForString,"fDistanceBetweenPartonsForString[fNuclTreeNumberOfStrings]/F");
            fNucleiCollisionsTree->Branch("stringRadiusVectorAngle", fNuclTreeStringRadiusVectorAngle,"fNuclTreeStringRadiusVectorAngle[fNuclTreeNumberOfStrings]/F");
            fNucleiCollisionsTree->Branch("stringOrigin", fStringOrigin,"fStringOrigin[fNuclTreeNumberOfStrings]/S");
        //    fNucleiCollisionsTree->Branch("isHardInteraction",fNuclTreeIsHardInteractionString,"fNuclTreeIsHardInteractionString[fNuclTreeNumberOfStrings]/O");
    }

    fNucleiCollisionsTree->Branch( "multFromAllStringsFictiveV0", &fMultFromAllStringsFictiveV0,"fMultFromAllStringsFictiveV0/I");
    fNucleiCollisionsTree->Branch( "multFromAllStringsFictiveMidEta", &fMultFromAllStringsFictiveMidEta,"fMultFromAllStringsFictiveMidEta/I");


    // ##### event loop
    for ( int iEvent = 0; iEvent < nEvents; iEvent++)
    {
        if ( iEvent % 10 == 0 )
            cout <<"generating " << (int)iEvent << "\r"; cout.flush();

        //building event, spec impact parameter if requested:
        fPtrNuclStruct->setEventId( iEvent ); //to be used in fPtrNuclStruct for output files (tmp?)
        fPtrNuclStruct->buildEvent();
        //        fPtrNuclStruct->drawEventStructure();

        //need some event checks/triggers here?..
        if ( !fPtrNuclStruct->isMBcollision() )
        {
            iEvent--;
            continue;
        }

        fImpactParameter = fPtrNuclStruct->getImpactParameter();
        fNuclTreeNumberOfStrings = fPtrNuclStruct->getNstrings();
        fNuclTreeRandomEventPlanePhi = fPtrNuclStruct->getRandomEventPlanePhi();
        fNuclTreeNu = fPtrNuclStruct->getNu();
        fMultFromAllStringsFictiveV0 = fPtrNuclStruct->getMultFromAllStringsFictiveV0();
        fMultFromAllStringsFictiveMidEta = fPtrNuclStruct->getMultFromAllStringsFictiveMidEta();
        //cout << fMultFromAllStringsFictive << endl;

        fNumberOfWoundedPartonsA = fPtrNuclStruct->getNumberOfWoundedPartonsA();
        fNumberOfWoundedPartonsB = fPtrNuclStruct->getNumberOfWoundedPartonsB();
        fNumberOfWoundedOnlyOncePartonsA = fPtrNuclStruct->getNumberOfWoundedOnlyOncePartonsA();
        fNumberOfWoundedOnlyOncePartonsB = fPtrNuclStruct->getNumberOfWoundedOnlyOncePartonsB();

        fNumberOfWoundedNucleonsA = fPtrNuclStruct->getNumberOfWoundedNucleonsA();
        fNumberOfWoundedNucleonsB = fPtrNuclStruct->getNumberOfWoundedNucleonsB();
        fNumberOfWoundedOnlyOnceNucleonsA = fPtrNuclStruct->getNumberOfWoundedOnlyOnceNucleonsA();
        fNumberOfWoundedOnlyOnceNucleonsB = fPtrNuclStruct->getNumberOfWoundedOnlyOnceNucleonsB();

        if ( fNuclTreeNumberOfStrings >= NMaxStrings )
        {
            cout << "AHTUNG: fNuclTreeNumberOfStrings >= NMaxStrings!!!" << endl;
            return;
        }

        // ##### loop over strings
        for ( int iString = 0; iString < fNuclTreeNumberOfStrings; iString++)
        {
//            fNuclTreeStringBoostAngle[iString] = fPtrNuclStruct->getStringBoostAngle(iString);
//            fNuclTreeStringBoostMagn[iString] = fPtrNuclStruct->getStringBoostMagn(iString);
            fNuclTreeStringRadiusVectorAngle[iString] = fPtrNuclStruct->getStringRadiusVectorAngle(iString);
            fNuclTreeStringX[iString] = fPtrNuclStruct->getStringXminusBover2(iString);
            fNuclTreeStringY[iString] = fPtrNuclStruct->getStringY(iString);
            fStringOrigin[iString] = fPtrNuclStruct->getStringOrigin(iString);
            fDistanceBetweenPartonsForString[iString] = fPtrNuclStruct->getDistanceBetweenPartonsForString(iString);
//            fNuclTreeIsHardInteractionString[iString] = fPtrNuclStruct->isHardInteractionString(iString);
        }

        fHist2D_impParVSnStrings->Fill( fImpactParameter, fNuclTreeNumberOfStrings );

        // ##### May 2018: core-corona:
        TH2I *fHist2DNstringsXY_thisEv_CoreCorona = fPtrNuclStruct->getHist2DNstringsXY_thisEv_CoreCorona();
        for ( int i = 0; i < fHist2DNstringsXY_thisEv_CoreCorona->GetNbinsX(); i++) //x
        {
            for ( int j = 0; j < fHist2DNstringsXY_thisEv_CoreCorona->GetNbinsY(); j++) //y
            {
                fStringsInCells[i][j] = fHist2DNstringsXY_thisEv_CoreCorona->GetBinContent(i+1,j+1);
//                cout << fStringsInCells[i][j] << " "; //endl;
            }
//            cout << endl;
        }


        fNucleiCollisionsTree->Fill();

    } // end of the event loop

    //draw the last generated event
    if (fDrawHistos)
    {
        fPtrNuclStruct->drawStatisticHists();
//        drawStatHists();
    }
    fPtrNuclStruct->finalActions();


    // ##### additional histogram with info (cross-section, etc) - March 2018
    const int nBinLabels = 2;
    TH1D *fAdditionalInfo = new TH1D( "fAdditionalInfo", "additional information from simulations;bin;value", nBinLabels,-0.5,nBinLabels-0.5);
    TString strBinLabels[nBinLabels] =
    {
//        "cross_seciton", // TUT OSHIBKA V NAPISANII!!!
        "cross_section",
        "mean_n_strings",
    };
    for( Int_t i = 1; i <= nBinLabels; i++ )
        fAdditionalInfo->GetXaxis()->SetBinLabel( i, strBinLabels[i-1].Data() );

    fAdditionalInfo->SetBinContent( 1, fPtrNuclStruct->getCrossSection() );
    fAdditionalInfo->SetBinContent( 2, fPtrNuclStruct->getMeanNstrings() );


    //
    TCanvas *canv_impParVSnStrings = new TCanvas("canv_impParVSnStrings","canvas impPar VS nStrings",380,280,800,600);
    fHist2D_impParVSnStrings->DrawCopy("colz");

    //write tree
    TFile* outFileNuclColTree = new TFile( Form( "%s/nuclCollisionsTree_nEv%d.root", fOutputDirName.Data(), nEvents ), "RECREATE" );
    outFileNuclColTree->cd();
    fNucleiCollisionsTree->Write();
    fHist2D_impParVSnStrings->Write();
    fAdditionalInfo->Write();
    outFileNuclColTree->Close();
}


