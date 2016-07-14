#include <Math/Vector3D.h>
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
//#include "TLorentzVector.h"
//#include "TVector3.h"
#include "TF1.h"
#include "TString.h"

#include "DecayInTwo.h"
#include "StringFragmentation.h"
#include "StringDescr.h"

#include <iostream>
using namespace std;

using namespace ROOT::Math;

//const double kMassPion = 0.14;
//const double particleMass = 0.14; //0.938; //0.14;
//const double particleMass = 0.938;

//const double mRho = 0.77;// fRand->Gaus(0.77,0.05);
//const double mPion = 0.14;

// gaussian boost profile
Double_t stringBoostGaus(Double_t *x, Double_t *par)
{
    return TMath::Exp( -0.5 * ( x[0]*x[0]/par[0]/par[0] ) );
}


//const double kMinEta = -3;
//const double kMaxEta = 3;
StringDescr::StringDescr() //: //ParticleArr(fNparticles)
{
    fRand = 0x0; //new TRandom3;
    //    fPtKickPerUnitMagn = 0.1;
    //    fLorVect = new TLorentzVector();
    decayInTwo = new DecayInTwo;
    strFr = new StringFragmentation;
    //    strFr->setStringEndPoinsY(-1,1,3);
//    strFr->setStringEndPoinsY(kMinEta,kMaxEta,2);
//    strFr->setStringEndPoinsY(-0.75,0.75,2);
//    strFr->setStringEndPoinsY(-3,3,2);
    strFr->setStringEndPoinsY(-4,4,2);

    if(0)
    {
        //    funcPt = new TF1( "funcPt", "[0]*TMath::Exp(-sqrt([1]*[1]+x*x)/[2])", 0., 5 );
        //    funcPt = new TF1( "funcPt", "[0]*TMath::Exp(-([1]*[1]+x*x)/[2])", 0., 50 );
        funcPt = new TF1( "funcPt", "[0]*TMath::Exp(-x*x/[2])", 0., 20 );
        funcPt->SetParameter(0,1);
        //    funcPt->SetParameter(1,0.14 ); //938);//14);
        //    funcPt->SetParameter(2,0.6);

        funcPt->SetParameter(1,0.6);

    }
    fFuncGausYBoost = new TF1( "funcGausYBoost", stringBoostGaus, -10, 10, 1 );
    //    fFuncGausYBoost->SetParameter(0,2.); //gauss sigma is of the order of mean string size (iz obshih soobrajeniy)
    fFuncGausYBoost->SetParameter(0,2.); //gauss sigma is of the order of mean string size (iz obshih soobrajeniy)
    //    double parVafabs( fRand->Gaus(2,0.2) ) + fabs( fRand->Gaus(0,2) );

    //    funcPt = new TF1( "funcPt", "x*[0] * x/sqrt(0.14*0.14+x*x) * pow( 1 + x/[1], -[2])", 0., 20 ); //2*TMath::Pi() );

    //Hagedorn formula, from ALICE pt spectra 900 GeV paper
    //    TString strFormula = Form("x*[0] * x/sqrt(%.2f*%.2f+x*x) * pow( 1 + x/[1], -[2])", particleMass, particleMass );

    //    funcPt = new TF1( "funcPt", strFormula.Data(), 0., 20 ); //2*TMath::Pi() );
    //    funcPt->SetParameter(0,1);
    //    funcPt->SetParameter(1,1.05);
    //    funcPt->SetParameter(2,7.92);





    //    TLorentzVector *rhoVector = new TLorentzVector;
    //    TLorentzVector vMother;
    vMothers = new TLorentzVector[N_LORENTZ_VECTORS_FROM_STRING];


}

StringDescr::~StringDescr() {}

void StringDescr::setRandomGenerator(TRandom* rand)
{
    fRand = rand;
    decayInTwo->setRandomGenerator(fRand);
    strFr->setRandomGenerator(fRand);
}

void StringDescr::hadronizeString( double stringBeta, double boostPhiDir )
{
    if ( !fRand )
    {
        cout << "StringDescr: fRand IS NOT INITIALIZED!!! exiting..." << endl;
        return;
    }

    //    int nParticlesInString = TMath::Nint( fRand->Gaus(5,0.5) );
    //    calcStringFragmentation();
    int nParticlesInString = strFr->decayStringIntoParticles( vMothers, /*fictionRhoPt*/0.6 );

    //    double boostX = fRand->Uniform( 0.6 , 0.90 );
    //    fFuncGausYBoost->SetParameter( 0, 1.8 + stringBeta*0.5 );


    //TMP!!!? make mother eta=0
    double motherPt, motherPhi, motherEta;
    if(1)for ( int iMother = 0; iMother < nParticlesInString; iMother++ )
    {
        vMother = &vMothers[iMother];
        motherPt = vMother->Pt();
        motherPhi = vMother->Phi();
        etaArrTmp[iMother] = vMother->Rapidity(); //keep values to set back (see below)
        motherEta = 0;
        vMother->SetPtEtaPhiM( motherPt, motherEta, motherPhi, vMother->M() );
    }

    //boost mothers in string
    stringBoostT.SetXYZ(0,0,0);
    if(1)for ( int iMother = 0; iMother < nParticlesInString; iMother++ )
    {
        vMother = &vMothers[iMother];
        double stringBetaAtY = stringBeta;// * fFuncGausYBoost->Eval( vMother->Rapidity() );
        stringBoostT.SetX( stringBetaAtY * cos(boostPhiDir) );//gRandom->Uniform( 0.0, 0.6 ) );
        stringBoostT.SetY( stringBetaAtY * sin(boostPhiDir) );//gRandom->Uniform( 0.0, 0.6 ) );

        if ( stringBoostT.Mag() > 0 )
            vMother->Boost( stringBoostT );
    }

    //boost mothers in string ALONG Z (TMP? redo?..) 26.12.2014
    // MOTIVATION: to smooth y distribution of rho and products!
    double boostBetaZ = stringBeta;
    //    boostBetaZ = fRand->Gaus( 0, stringBeta/2 );
    //    while ( fabs(boostBetaZ) > 1 )
    //        boostBetaZ = fRand->Gaus( 0, stringBeta/2 );

    stringBoostZ.SetZ( fRand->Uniform( -boostBetaZ, boostBetaZ ) );

    if(0)for ( int iMother = 0; iMother < nParticlesInString; iMother++ )
    {
        vMother = &vMothers[iMother];
        //                 prepareMother( vMother );
        if ( stringBoostZ.Mag() > 0 )
            vMother->Boost( stringBoostZ );
    }

    // TMP!!!? randomize eta
    if(1)for ( int iMother = 0; iMother < nParticlesInString; iMother++ )
    {
        vMother = &vMothers[iMother];
        //                 prepareMother( vMother );
        //if ( fabs(vMother->Eta()) < 2 )
        motherPt = vMother->Pt();
        motherPhi = vMother->Phi();
        //                    motherEta = vMother->Eta();
        //                    motherEta = fRand->Uniform(-3,3);
        //                    vMother->SetPtEtaPhiM( motherPt, motherEta, motherPhi, vMother->M() );

//        double y = etaArrTmp[iMother];//fRand->Uniform(-4,4);
        double y = fRand->Uniform(-3, 3);
        double mT = sqrt( motherPt*motherPt + vMother->M()*vMother->M() );
        double pX = motherPt * cos(motherPhi);
        double pY = motherPt * sin(motherPhi);
        double pZ = mT*sinh(y);
        vMother->SetXYZM( pX, pY, pZ, vMother->M() );
    }



    // ##### particles from string decay
    //    cout << "nParticlesInString=" << nParticlesInString << endl;
    int nP = 0;
    for ( int iMother = 0; iMother < nParticlesInString; iMother++ )
    {
        vMother = &vMothers[iMother];

        //            cout << "mother after boost: eta, pt: " << vMother->Eta() << " " << vMother->Pt() <<  endl;
        //        cout << "vMother->M(): " << vMother->M() <<  endl;

//        cout << "StringDescr::hadronizeString : vMother->M() = " << vMother->M() << endl;

        if ( fabs( vMother->M() - mRho ) < 2*mRhoWidth/2 ) //rho meson
        {

            //perform decay to pions if it is rho-meson
            decayInTwo->generateDecayProducts( *vMother, v1, v2 );

            // take into account rho type - charged or neutral, assume equal probabilities
            int rhoCharge = strFr->probabilityChargePlusMinusNeutral();
            int charge[2];
            //            cout << "rhoCharge=" << rhoCharge << endl;
            if ( rhoCharge == 0 ) //rho neutral -> two charged pions
            {
                charge[0] = +1;
                charge[1] = -1;
            }
            else //rho charged -> charged pion and neutral pion
            {
                charge[0]  = 0;
                charge[1]  = rhoCharge;
                //                        cout << products[1]->charge << endl;
            }
            // fill array with particles
            if ( charge[0] != 0 )
            {
                fParticles[nP].eta = v1.Eta();
                fParticles[nP].phi = v1.Phi();
                fParticles[nP].pt  = v1.Pt();
                fParticles[nP].charge  = charge[0];
                fParticles[nP].pid  = kPid_pion;
                FixAngleInTwoPi(fParticles[nP].phi);
                nP++;
            }
            if ( charge[1] != 0 )
            {
                fParticles[nP].eta = v2.Eta();
                fParticles[nP].phi = v2.Phi();
                fParticles[nP].pt  = v2.Pt();
                fParticles[nP].charge  = charge[1];
                fParticles[nP].pid  = kPid_pion;
                FixAngleInTwoPi(fParticles[nP].phi);
                nP++;
            }
        }
        else // not a resonanse
        {
            int charge  = 0;
            int pid  = -1;

            //charge and pid depending on mass:
            if ( fabs( vMother->M() - mPion ) < 0.01 ) // pions
            {
                charge  = strFr->probabilityChargePlusMinusNeutral();
                pid  = kPid_pion;
            }
            else if ( fabs( vMother->M() - mKaon ) < 0.01 ) // kaons
            {
                charge  = strFr->probabilityChargePlusMinusNeutral();
                pid  = kPid_kaon;
            }
            else if ( fabs( vMother->M() - mProton ) < 0.01 ) // protons
            {
                charge  = strFr->probabilityChargePlusMinus();
                pid  = kPid_proton;
            }
            else if ( fabs( vMother->M() - mPhi ) < 0.01 ) // phi
            {
                charge  = 0;
                pid  = kPid_phi;
            }
            else if ( fabs( vMother->M() - mLambda ) < 0.01 ) // lambda
            {
                charge  = 0;//strFr->probabilityChargePlusMinus();
                pid  = kPid_Lambda;
            }

            // fill array with particle
            if ( charge != 0 || pid  == kPid_Lambda || pid  == kPid_phi )
            {
                fParticles[nP].eta = vMother->Eta();
                fParticles[nP].phi = vMother->Phi();
                fParticles[nP].pt  = vMother->Pt();
                fParticles[nP].charge  = charge;
                fParticles[nP].pid  = pid;
                FixAngleInTwoPi(fParticles[nP].phi);
                nP++;
            }
        }


    }

    //total number of charged particles
    fNparticles = nP;
    //    cout << "fNparticles=" << fNparticles << endl;



}



//int StringDescr::probabilityChargePlusMinusNeutral()
//{
//    //all +,-,0 have probability 1/3
//    int chargeSign = 0;

//    double randForCharge = fRand->Uniform(-2,2);

//    if ( randForCharge < 0 )
//        chargeSign = 1; // neutral
//    else if ( randForCharge > 1 )
//        chargeSign = 1; // positive
//    else
//        chargeSign = -1; // negative
//    return chargeSign;
//}


void StringDescr::makeTwoJets()
{
    int nP = 0;

    double phiJet = fRand->Uniform( 0, 2*TMath::Pi() );
    int nJets = 1;//( fRand->Uniform() > 0.5 ? 2 : 1 ); // quenching of some fraction of jets

    if(1)for ( int iJet = 0; iJet < nJets; iJet++ )
    {
        int nParticlesInJet = TMath::Max(2,TMath::Nint( fRand->Gaus(4.,0.5) ));
        if (iJet == 1 )
        {
            phiJet += TMath::Pi() + fRand->Gaus(0,0.05);
            FixAngleInTwoPi(phiJet);
            FixAngleInTwoPi(phiJet);
            nParticlesInJet = TMath::Max(1,TMath::Nint( fRand->Gaus(1.5,0.5) )); //fewer particles due to quenching
        }
        double etaJet = fRand->Uniform( -2.5,2.5); //kMinEta,kMaxEta);

        for ( int iP = 0; iP < nParticlesInJet; iP++ )
        {
//            vMother->SetPtEtaPhiM( motherPt, motherEta, motherPhi, mPions );
            //        fParticles[nP].eta = vMother->Eta();
            //        fParticles[nP].phi = vMother->Phi();
            //        fParticles[nP].pt  = vMother->Pt();
            fParticles[nP].eta = etaJet + fRand->Gaus(0,0.15);//0.2);
            double phi = phiJet + fRand->Gaus(0,0.12);//15);
            FixAngleInTwoPi(phi);
            fParticles[nP].phi = phi;
            fParticles[nP].pt  = fRand->Exp(2.);
            fParticles[nP].charge  = ( fRand->Uniform() > 0.5 ? 1 : -1 );
//            cout << "phi=" << phi << ", eta=" << fParticles[nP].eta
//                 << ", pt=" << fParticles[nP].pt << endl;
            nP++;

        }
    }

    //total number of charged particles
    fNparticles = nP;
//        cout << "fNparticles=" << fNparticles << endl;


}


