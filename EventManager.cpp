#include "/Users/macbook/alice/simpleAnalysis/commonTools/Tools.cxx"  //  "../../commonTools/Tools.cxx"
//#include "../../commonTools/SimpleTrack.cxx"


#include "EventManager.h"
#include "StringGeneration/NucleusStructure.h"
#include "StringDecayer/StringDescr.h"
#include "/Users/macbook/alice/simpleAnalysis/simpleEventAnalyzer/AliSimpleEvent.h"    // "../../simpleEventAnalyzer/AliSimpleEvent.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
using namespace std;

const float kEtaBoundForCountingParticles = 0.8;
const Int_t nCentralityBins = 20;//10;


const int nMaxPartForHist = 4000;
const int nBinsForNparticlesHist = 500;
//const int nMaxPartForHist = 400;
//const int nBinsForNparticlesHist = 400;



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


EventManager::EventManager():
    fPtrNuclStruct(0x0)
  , fFillEventTree(false)
  //    fFlagGenerateCentralEventByHand(false)
  //  , fFlagGenerateSemicentralEventByHand(false)
  //  , fImpactParameterByHand(-1)
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

    fHistNeventsInCentralityClasses = new TH1D("fHistNeventsInCentralityClasses", "N particles in class;class;n events", nCentralityBins, -0.5, nCentralityBins-0.5 );
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



    
    fOutputFileName = "testOutput.root";
    fOutputFile = 0x0;
    fDrawHistos = true;

    //    fPtCutMin = 0.;
    //    fPtCutMax = 100.;
}

EventManager::EventManager(const EventManager& ) {}

EventManager::~EventManager() {}


void EventManager::generateEvents( NucleusStructure *nuclStruct, StringDescr *strDescr, int nEvents )
{
    cout << "generating events..." << endl;
    
    //reset all hists
    fHistSources->Reset();
    fHistParticlesInSource->Reset();
    fHistParticlesInEvent->Reset();
    fHistPt->Reset();
    fHistPtAfterCuts->Reset();
    fHistPtBeforeKick->Reset();
    fHistEta->Reset();
    for ( int iPtBin = 0; iPtBin < 5; iPtBin++)
        fHistEtaInPtCuts[iPtBin]->Reset();
    fHistPhi->Reset();

    fPtrNuclStruct = nuclStruct;
    
//    SimpleTrack *fSimpleTrack = new SimpleTrack;

    //output file for the events
    fOutputFileName = Form( "tmpOutputs/eventTree_nEv%d.root", nEvents );
    TFile* outFile = new TFile( fOutputFileName, "RECREATE" );
    //Form(
    //                  "%s/%s_%devents.root"
    //                , dirOutName.Data(), fileOutName.Data(), nev), "RECREATE" );


    //create simple event object
    AliSimpleEvent *fSimpleEvent = new AliSimpleEvent();// By setting the value, we own the pointer and must delete it.

    //organize a tree
    TTree *fEventTree = new TTree("EventTree","An example of a ROOT tree");
    fEventTree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
    fEventTree->SetCacheSize(10000000);  // set a 10 MBytes cache (useless when writing local files)
    Int_t bufsize  = 64000;
    Int_t split  = 2;
    if ( split )
        bufsize /= 4;

    TBranch *branch = fEventTree->Branch("event", &fSimpleEvent, bufsize, split );
    branch->SetAutoDelete(kFALSE);

    if( split >= 0 )
        fEventTree->BranchRef();


    // ##### event loop
    for ( int iEvent = 0; iEvent < nEvents; iEvent++)
    {
        if ( iEvent % 100 == 0 )
            printf("generating %d event...\n",(int)iEvent );

        //building event, spec impact parameter if requested:
        fPtrNuclStruct->setEventId( iEvent ); //to be used in d for output files (tmp?)
        fPtrNuclStruct->buildEvent();
        //        fPtrNuclStruct->drawEventStructure();

        //need some event checks/triggers here?..
        if ( !fPtrNuclStruct->isMBcollision() )
        {
            iEvent--;
            continue;
        }

        //do analysis using NuclearStructure data
        if (1)
        {
            // ##### get random EP for the event!
            float phiEP = fPtrNuclStruct->getRandomEventPlanePhi();

            //prepare simple event
            if ( fFillEventTree )
            {
                fSimpleEvent->SetHeader( iEvent, -1, -1, 2.0, 3.0 );
                fSimpleEvent->SetVertexPos( iEvent*2,0,0 );//vertex->GetXv(), vertex->GetYv(), vertex->GetZv() );
                fSimpleEvent->StartEventFilling();
            }

            // ##### loop over strings: hadronized and count multiplicity in VZERO
            const double cutEtaCMS = 2;
            int nParticles = 0;
            int nParticlesWithinEtaPtCuts = 0;
            for ( int iString = 0; iString < fPtrNuclStruct->getNstrings(); iString++)
            {
                int nParticlesInString = 0;
                if ( fPtrNuclStruct->isHardInteractionString(iString) == 0 ) //soft interaction - string
                {
                    strDescr->hadronizeString(fPtrNuclStruct->getStringBoostMagn(iString), fPtrNuclStruct->getStringBoostAngle(iString) );
                    nParticlesInString = strDescr->getNparticles();
                }
                else //hard interaction - two jets
                {
                    strDescr->makeTwoJets();
                    nParticlesInString = strDescr->getNparticles();
                }

                for ( int iP = 0; iP < nParticlesInString; iP++ )
                {
                    ParticleDescr *p = strDescr->getParticle(iP);
                    //rotate all tracks in phi
                    float phiMod = phiEP + p->phi;
                    FixAngleInTwoPi(phiMod);

                    //spec counting to compare with publications
                    //                    if ( fabs(p->eta) < 0.3 && ( p->pt > 0.15 && p->pt < 10 ) )
//                    if ( p->eta < -3.7 || p->eta > 1. )//fabs(p->eta) > 3.7 ) // out of our intresting range
//                        continue;
                    if ( fabs(p->eta) > cutEtaCMS )//fabs(p->eta) > 3.7 ) // out of our intresting range
                        continue;

//                    if ( p->eta > -3.7 && p->eta < -1.7 )  //VZERO-C
                   // if ( fabs(p->eta) < 1. )  //STAR TPC
                    if ( fabs(p->eta) < cutEtaCMS )  //CMS
                    {
                        nParticlesWithinEtaPtCuts++;

                    }


                    if ( fabs(p->eta) < cutEtaCMS )//0.5 )//&& ( p->pt > 0. && p->pt < 1000 ) )
                    {
                        if ( p->pt > 0.05 )
                        {
                            double weight = 1./ ( 2*TMath::Pi()*p->pt );
                            fHistPtAfterCuts->Fill(p->pt, weight );
                        }
//                        fHistPtAfterCuts->Fill( p->pt, 1./ ( 2*TMath::Pi()*p->pt)  );
                        //nParticlesWithinEtaPtCuts++;
                    }
                    fHistEta->Fill(p->eta);
                    if ( p->pt > 0.1 ) fHistEtaInPtCuts[0]->Fill(p->eta);
                    if ( p->pt > 0.5 ) fHistEtaInPtCuts[1]->Fill(p->eta);
                    if ( p->pt > 1.0 ) fHistEtaInPtCuts[2]->Fill(p->eta);
                    if ( p->pt > 2.0 ) fHistEtaInPtCuts[3]->Fill(p->eta);
                    if ( p->pt > 3.0 ) fHistEtaInPtCuts[4]->Fill(p->eta);
                    fHistPhi->Fill(phiMod);
                    if ( fabs(p->eta) < 0.5 )
                    {
                        if ( p->pt > 0.0 )
                        {
                            double weight = 1./ ( 2*TMath::Pi()*p->pt );
                            fHistPt->Fill(p->pt, weight );
                        }
                        if ( p->ptBeforeKick > 0.0 )
                        {
                            double weight = 1./ ( 2*TMath::Pi()*p->ptBeforeKick );
                            fHistPtBeforeKick->Fill(p->ptBeforeKick, weight );
                        }

                    }
                    fHist2DEtaPhi->Fill( p->eta, phiMod );

                    // ##### add particle to analysers
//                    fSimpleTrack->id = nParticles;
//                    fSimpleTrack->eta = p->eta;
//                    fSimpleTrack->phi = phiMod;
//                    fSimpleTrack->pt = p->pt;
//                    fSimpleTrack->charge = 0;//charge;

                    nParticles++;

                    if ( fFillEventTree )
                        fSimpleEvent->AddTrack(  p->pt,  p->eta, phiMod, p->charge, -1. );//random,ptmin);
                }
            }
            fHistParticlesInEvent->Fill(nParticles);
            fHistParticlesInCutConditionInEvent->Fill(nParticlesWithinEtaPtCuts/2/cutEtaCMS); // !
            fHistParticlesInCutConditionVsNu->Fill(nParticlesWithinEtaPtCuts/2/cutEtaCMS, fPtrNuclStruct->getNu() );

            if ( fFillEventTree )
            {
                fSimpleEvent->GetHeader()->SetCentrality(nParticlesWithinEtaPtCuts/2/cutEtaCMS);
                fSimpleEvent->FinishEventFilling();
            }
            fEventTree->Fill();  //fill the tree
        }

    } // end of the event loop


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

        fPtrNuclStruct->drawStatisticHists();
        drawStatHists();
    }
    fPtrNuclStruct->finalActions();


    //    delete d;
    //    delete strDescr;

    //    fQAhistos.fHistNch->SetName("hist70-80");
    //    fQAhistos.fHistNch->Write();
    //    f->Close();


    // ##### multiplicity quantiles
    //prepare quantiles
    Double_t centrMultBounds[nCentralityBins];  // array to contain the quantiles
    //    bool flagHaveQuantiles = false; //flag that we have calculated multiplicity quantiles

    TCanvas *fCanvMultClasses = new TCanvas("canvMultClasses","multiplicity classes",50,50,800,600);
    fCanvMultClasses->cd()->SetLogy();

    fHistParticlesInCutConditionInEvent->SetLineColor(kRed);
    fHistParticlesInCutConditionInEvent->DrawCopy();
    fHistParticlesInCutConditionInEvent->SaveAs("tmpOutputs/histParticlesInCutCondition.root");
    //    fHistParticlesInCutConditionInEvent->GetQuantiles();

    cout << ">>> Mean Nch overall = " << fHistParticlesInCutConditionInEvent->GetMean() << endl;

    //remove zero bin:
    if ( fHistParticlesInCutConditionInEvent->GetBinContent(1) > 0 )
        fHistParticlesInCutConditionInEvent->SetBinContent(1, 0);
    // get quantiles
    getQuantiles( fHistParticlesInCutConditionInEvent, nCentralityBins, centrMultBounds );
    //    flagHaveQuantiles = true;

    //fill multiplicity boundaries for each centrality class
    TH1D *fHistMultClassBoundaries = new TH1D( "fHistMultiplicityBoundaries", "fHistMultClassBoundaries"
                                               , nCentralityBins, -0.5, nCentralityBins-0.5 );
    TH1D *fHistMultClassMeanNch = new TH1D( "fHistMultClassMeanNch", "fHistMultClassMeanNch"
                                               , nCentralityBins, -0.5, nCentralityBins-0.5 );
    for ( int iBin = 0; iBin < nCentralityBins; iBin++ )
    {
        fHistMultClassBoundaries->SetBinContent( iBin+1, centrMultBounds[iBin] );
    }

    //draw centrality classes on mult hist
//    int nBinsForClassHists = fHistParticlesInCutConditionInEvent->GetNbinsX();
    TH1D *fHistCentrClass[nCentralityBins];
    for ( int iCentrClass = 0; iCentrClass < nCentralityBins; iCentrClass++ )
        fHistCentrClass[iCentrClass] = new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass), nBinsForNparticlesHist, -0.5, nMaxPartForHist-0.5);

    int centrClassId = 0;
    for ( int iBin = 0; iBin < fHistParticlesInCutConditionInEvent->GetNbinsX(); iBin++ )
    {
        //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
        if ( fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
            centrClassId++;
        if ( centrClassId >= nCentralityBins )
            break;
        double binContent = fHistParticlesInCutConditionInEvent->GetBinContent( iBin+1 );
        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    }
    for ( int iCentrClass = 0; iCentrClass < nCentralityBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
        fHistCentrClass[iCentrClass]->DrawCopy( "same" );

        double meanNchInCentrBin = fHistCentrClass[iCentrClass]->GetMean();
        fHistMultClassMeanNch->SetBinContent( iCentrClass+1, meanNchInCentrBin );
        cout << "meanNch in centrality bin " << iCentrClass << ": " << meanNchInCentrBin << endl;
    }


    //write objects
    outFile->cd();
    fEventTree->Write();
    fHistMultClassBoundaries->Write();
    fHistMultClassMeanNch->Write();
    fHistParticlesInCutConditionVsNu->Write();
    outFile->Close();
    //    delete outFile;

}

void EventManager::cleanup()
{
    fOutputFile->Close();
    fOutputFile = 0x0;
}

void EventManager::drawStatHists()
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
    //    const Int_t nCentralityBins = 10;
    //    Double_t centrMultBounds[nCentralityBins];  // array to contain the quantiles
    //    getQuantiles( fHistParticlesInCutConditionInEvent, nCentralityBins, centrMultBounds );
    //    //draw centrality classes on mult hist
    //    TH1D *fHistCentrClass[nCentralityBins];
    //    for ( int iCentrClass = 0; iCentrClass < nCentralityBins; iCentrClass++ )
    //        fHistCentrClass[iCentrClass] = new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass), fHistParticlesInCutConditionInEvent->GetNbinsX(), 0, 8000);

    //    int centrClassId = 0;
    //    for ( int iBin = 0; iBin < fHistParticlesInCutConditionInEvent->GetNbinsX(); iBin++ )
    //    {
    //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
    //        if ( fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
    //            centrClassId++;
    //        if ( centrClassId >= nCentralityBins )
    //            break;
    //        double binContent = fHistParticlesInCutConditionInEvent->GetBinContent( iBin+1 );
    //        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    //    }
    //    for ( int iCentrClass = 0; iCentrClass < nCentralityBins; iCentrClass++ )
    //    {
    //        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
    //        fHistCentrClass[iCentrClass]->DrawCopy( "same" );
    //    }



    fCanv->cd(2)->SetLogy();
    fHistPtAfterCuts->SetLineColor(kRed);
    fHistPtBeforeKick->SetLineColor(kGreen);
    fHistPt->DrawNormalized();
    fHistPtAfterCuts->DrawNormalized("same");
//    fHistPtBeforeKick->DrawNormalized("same");
//    fHistPt->Draw();
//    fHistPtAfterCuts->Draw("same");
//    fHistPtBeforeKick->Draw("same");

    fCanv->cd(3);
    fHistEta->DrawNormalized();
    fHistEtaInPtCuts[0]->SetLineColor(kRed);
    fHistEtaInPtCuts[1]->SetLineColor(kGreen);
    fHistEtaInPtCuts[2]->SetLineColor(kMagenta);
    fHistEtaInPtCuts[3]->SetLineColor(kOrange);
    fHistEtaInPtCuts[4]->SetLineColor(kBlack);
    for ( int iPtBin = 0; iPtBin < 3/*5*/; iPtBin++)
        fHistEtaInPtCuts[iPtBin]->DrawNormalized("same");
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

    fCanv->SaveAs("tmpOutputs/canvEvManagerStats.root");
    fCanv->SaveAs("tmpOutputs/canvEvManagerStats.eps");
    fCanv->SaveAs("tmpOutputs/canvEvManagerStats.png");


    //write stats to file
    TFile *fileEvManagerStats = new TFile("tmpOutputs/stats_EvManager.root","RECREATE");
    fHistParticlesInEvent->Write();
    fHistParticlesInCutConditionInEvent->Write();
    fHistPt->Write();
    fHistPtAfterCuts->Write();
    fHistEta->Write();
    fHistPhi->Write();
    fileEvManagerStats->Close();

    // ###### spec output txt files
    ofstream fout( "tmpOutputs/tmpTextOutput_EventManager_ptVsNch.txt", ios::out | ios::binary);
    fout << fHistParticlesInCutConditionInEvent->GetMean() //fPtrNuclStruct->getImpactParameterByHand() //_0_100
         << " " << fHistPtAfterCuts->GetMean()
         << " " << fHistPtAfterCuts->GetMeanError()
         << endl;
    fout.close();


}
