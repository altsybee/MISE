#include "Riostream.h"
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
//#include "TProfile.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include "TArrow.h"
#include "TCanvas.h"


#include "TRandom3.h"


inline void FixAngleInTwoPi( double &lPhi )
{
    if ( lPhi > 2 * TMath::Pi() )
        lPhi -= 2 * TMath::Pi();
    else if ( lPhi < 0 )
        lPhi += 2 * TMath::Pi();
}


Double_t funcPt(Double_t *x, Double_t *par)
{
    return par[0]*x[0]*TMath::Exp(-TMath::Pi()*(par[1]*par[1]+x[0]*x[0])/par[2]);
}

void flowToy()
{

    double R = 6.7;
    double d = 5; // d = b/2



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



    TH1D *fHistPt = new TH1D("fHistPt", "p_{T} distribution;p_{T} (GeV/c);n tracks", 1000, 0.0, 20);
    TH1D *fHistAlpha = new TH1D("fHistAlpha", "#alpha distribution", 200, -0.5, 7 );

    TCanvas *fCanvEventView = new TCanvas("Event Canvas","Event View",600,70,700,700);
    TArrow **arrowsR = new TArrow* [1000];
    TArrow **arrowsP = new TArrow* [1000];
    TArrow **arrowsSpecAngle = new TArrow* [1000];

    gRandom->SetSeed(0);

    TF1 *funcPartDistr = new TF1( "funcPartDistr", "1+0.5*TMath::Cos(2*(x-TMath::PiOver2()))", -TMath::PiOver2(), TMath::PiOver2() ); //0, TMath::TwoPi() );
//    TF1 *funcPartDistr = new TF1( "funcPartDistr", "1+0.5*TMath::Cos(2*(x-0))", -TMath::PiOver2(), TMath::PiOver2() ); //0, TMath::TwoPi() );


    const double mRho = 0.775;// fRand->Gaus(0.77,0.05);
    const double mRhoWidth = 0.16;
    const double mPion = 0.1395;
    const double mKaon = 0.494;
    const double mProton = 0.938;

    TF1 *funcPtBoltzmanLikePion = new TF1( "funcPtBoltzmanLikePion", "[0]*x*TMath::Exp(-TMath::Pi()*([1]*[1]+x*x)/[2])", 0, 10 );
//    TF1 *funcPtBoltzmanLikePion = new TF1( "funcPtBoltzmanLikePion", funcPt, 0, 5, 3 );
    funcPtBoltzmanLikePion->SetNpx(1000);
    funcPtBoltzmanLikePion->SetParameter( 0, 1 );
    funcPtBoltzmanLikePion->SetParameter( 1, mPion );
    funcPtBoltzmanLikePion->SetParameter( 2, 0.07 ); // t = 0.568 ± 0.001 GeV2

    TF1 *funcPtBoltzmanLikeKaon = new TF1( "funcPtBoltzmanLikeKaon", "[0]*x*TMath::Exp(-TMath::Pi()*([1]*[1]+x*x)/[2])", 0, 10 );
//    TF1 *funcPtBoltzmanLikeKaon = new TF1( "funcPtBoltzmanLikeKaon", funcPt, 0, 5, 3 );
    funcPtBoltzmanLikeKaon->SetParameter( 0, 5 );
    funcPtBoltzmanLikeKaon->SetParameter( 1, mKaon );
    funcPtBoltzmanLikeKaon->SetParameter( 2, 0.25 );//0.568 );

    TF1 *funcPtBoltzmanLikeProton = new TF1( "funcPtBoltzmanLikeProton", "[0]*x*TMath::Exp(-TMath::Pi()*([1]*[1]+x*x)/[2])", 0, 10 );
//    TF1 *funcPtBoltzmanLikeProton = new TF1( "funcPtBoltzmanLikeProton", funcPt, 0, 5, 3 );
    funcPtBoltzmanLikeProton->SetParameter( 0, 18 );
    funcPtBoltzmanLikeProton->SetParameter( 1, mProton );
    funcPtBoltzmanLikeProton->SetParameter( 2, 0.568 );


    // Tsallis:
    // (from http://www.nikhef.nl/pub/services/biblio/theses_pdf/thesis_M_Chojnacki.pdf
    // and http://arxiv.org/pdf/1504.00024v2.pdf)
    funcPtBoltzmanLikePion = new TF1( "funcPtBoltzmanLikePion", "[0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])", 0, 5 );
    funcPtBoltzmanLikePion->SetParameter(0,5.7);
    funcPtBoltzmanLikePion->SetParameter(1,0.12);
    funcPtBoltzmanLikePion->SetParameter(2,mPion);

    funcPtBoltzmanLikeKaon = new TF1( "funcPtBoltzmanLikeKaon", "[0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])", 0, 5 );
    funcPtBoltzmanLikeKaon->SetParameter( 0, 6.7 );
    funcPtBoltzmanLikeKaon->SetParameter( 1, 0.195 );
    funcPtBoltzmanLikeKaon->SetParameter( 2, mKaon );//0.568 );

    funcPtBoltzmanLikeProton = new TF1( "funcPtBoltzmanLikeProton", "[0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])", 0, 5 );
    funcPtBoltzmanLikeProton->SetParameter( 0, 6.3 );
    funcPtBoltzmanLikeProton->SetParameter( 1, 0.212 );
    funcPtBoltzmanLikeProton->SetParameter( 2, mProton );



    if(1)
    {
        funcPtBoltzmanLikePion->SetLineColor( kBlue );
        funcPtBoltzmanLikePion->Draw();
        funcPtBoltzmanLikeKaon->SetLineColor( kGreen );
        funcPtBoltzmanLikeKaon->Draw("same");
        funcPtBoltzmanLikeProton->SetLineColor( kRed );
        funcPtBoltzmanLikeProton->Draw("same");

        return;
    }




    int nEvents = 1;
    for ( int ev = 0; ev < nEvents; ev++)
    {
        if ( ev % 10 == 0 )
            cout << "generating " << (int)ev << "\r"; cout.flush();

        //random "reaction plane"
        double phiRP = gRandom->Uniform(0, TMath::TwoPi() );

        fNumberOfTracks = 50;//500;
//        double fluctPt = gRandom->Uniform( 0., 1 );
        double fluctPt = gRandom->Gaus(1,0.1); //1.;//gRandom->Exp( 0.8 );
        for( int tr = 0; tr < fNumberOfTracks; tr++ )
        {
//            double alpha = gRandom->Uniform( -TMath::PiOver2(), TMath::PiOver2() );
            double alpha = funcPartDistr->GetRandom();

            double beta = TMath::ASin( d/(R/TMath::Sin(alpha)) );
            double gamma = TMath::Pi() - (TMath::Pi()-alpha) - beta;
            double r = TMath::Sin(gamma)*R/TMath::Sin(alpha);

            if ( gRandom->Uniform( 0, 1 ) < 0.5 )
                alpha = TMath::Pi() - alpha;
            FixAngleInTwoPi( alpha );

            fHistAlpha->Fill( alpha );

//            cout << "alpha1 = " << alpha << endl;
            alpha += phiRP;
            FixAngleInTwoPi( alpha );
//            cout << "alpha2 = " << alpha << endl;



//            cout << r << endl;
            double rx = r*cos(alpha)/10;
            double ry = r*sin(alpha)/10;

            //spec angle for p:
            double angleP = alpha + TMath::PiOver2();
            FixAngleInTwoPi( angleP );

//            double p = TMath::Exp(1/r);// /10;//(r-1)/3;//1/r;
//            double p = gRandom->Exp(r/3);// /10;//(r-1)/3;//1/r;
//            double p = 0.3+gRandom->Exp(r/3);// /10;//(r-1)/3;//1/r;
//            double p = funcPtBoltzmanLikeProton->GetRandom()*1.2*r;// /10;//(r-1)/3;//1/r;
            double p = funcPtBoltzmanLikePion->GetRandom()*1.2*r;// /10;//(r-1)/3;//1/r;

            p *= fluctPt; // fluctuating pt!

            double px = p*cos(angleP)/2;
            double py = p*sin(angleP)/2;

            fHistPt->Fill(p);


//            if ( gRandom->Uniform( 0, 1 ) < 0.5 )
//            {
//                rx *= -1;
//                px *= -1;
//            }

            fTrackPhi[tr] = alpha;
            fTrackPt[tr] = p;
            fTrackEta[tr] = gRandom->Uniform( -3, 3 );

            fTrackCharge[tr] = 0;
            fTrackPID[tr] = 0;

            if ( ev == 0 ) //draw only first event
            {
                // draw positions
                arrowsR[tr] = new TArrow( 0.5, 0.5,
                                          0.5+rx, 0.5+ry,
                                          0.025);
                arrowsR[tr]->SetLineColor( kGreen );//kMagenta+2 );
                arrowsR[tr]->SetLineWidth( 2 );
                arrowsR[tr]->Draw();

                // draw momentum
                arrowsP[tr] = new TArrow( 0.5, 0.5,
                                          0.5+px, 0.5+py,
                                          0.025);
                arrowsP[tr]->SetLineColor( kRed );//kMagenta+2 );
                arrowsP[tr]->SetLineWidth( 2 );
//                arrowsP[tr]->Draw();


                if(1)
                {
                    for( int i = 0; i < 20; i++ )
                    {
                        double specAngle = gRandom->Uniform(0, TMath::TwoPi() );
                        double diff = alpha-specAngle;
                        FixAngleInTwoPi(diff);
                        const double extra = TMath::Pi()/8; //spec angle to constain even more
                        if ( diff > TMath::PiOver2()-extra && diff < 3*TMath::PiOver2()+extra )
                            continue;

                        double specx = 0.1*cos(specAngle)/2;
                        double specy = 0.1*sin(specAngle)/2;


                        arrowsSpecAngle[i] = new TArrow( 0.5+rx, 0.5+ry,
                                                  0.5+rx+specx, 0.5+ry+specy,
                                                  0.025);
                        arrowsSpecAngle[i]->SetLineColor( kBlue );
                        arrowsSpecAngle[i]->SetLineWidth( 2 );
                        arrowsSpecAngle[i]->Draw();
                    }
                }

            } // end of draw

        } // end of track loop

        fEventTrackTree->Fill();

    } // end of event loop

    //output file for the events
    TFile* outFile = new TFile( "output_toy_events.root", "RECREATE" );
    fEventTrackTree->Write();
    outFile->Close();

    fCanvEventView->SaveAs( "evView.eps" );


    TCanvas *canvQA = new TCanvas("canvQA","canvQA",600,70,800,600);
    canvQA->Divide(2,2);
    canvQA->cd(1);
    fHistPt->DrawCopy();
    canvQA->cd(2);
    fHistAlpha->DrawCopy();

    cout << "finished." << endl;



}


