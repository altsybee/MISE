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

    fFuncPowerLawPt = new TF1( "funcPowerLawPt", "[0]*TMath::Power(x,[1])", 1.0, 50);
//    fFuncPowerLawPt->SetParameters( 1, -4.9 ); // from CMS paper https://arxiv.org/pdf/1104.3547.pdf
    fFuncPowerLawPt->SetParameters( 1, -2.838 ); // from 0.2/distImpPar - /opt/mygit/MISE/analysis/RAA/play_with_partonIntDistHistos


    // spec fit function constructed for pt in [3,20]
    // to have smooth rise of the HARD PART (tuned using LHC10h Pb-Pb minBias spectra)
//    fFunc_HARD_COMBINED = new TF1( "fFunc_HARD_COMBINED", "x<3.6 ? 7100 : [0]+[1]*TMath::Power(x,[2]) + [3]*(TMath::Exp([4]*(x+[5])))", 0.2/*3.0*/, 20);
//    fFunc_HARD_COMBINED = new TF1( "fFunc_HARD_COMBINED", "x<4 ? 2e4-x*12800./4.0 :[0]+[1]*TMath::Power(x,[2]) + [3]*(TMath::Exp([4]*(x+[5])))", 0.5/*3.0*/, 20);
    fFunc_HARD_COMBINED = new TF1( "fFunc_HARD_COMBINED", "x<4 ? 0.92*1e5*TMath::Exp(-0.65*x) :[0]+[1]*TMath::Power(x,[2]) + [3]*(TMath::Exp([4]*(x+[5])))", 0.5/*3.0*/, 20);
//    fFunc_HARD_COMBINED = new TF1( "fFunc_HARD_COMBINED", "[0]+[1]*TMath::Power(x,[2]) + [3]*(TMath::Exp([4]*(x+[5])))", 3.0, 20);
//    fFunc_HARD_COMBINED->SetParameters( 1e6, 2, 1, 0, -2 );
    fFunc_HARD_COMBINED->SetParameter( 0, 1.54936e+00 );
    fFunc_HARD_COMBINED->SetParameter( 1, 1.65448e+07 );
    fFunc_HARD_COMBINED->SetParameter( 2, -5.13362e+00 );
    fFunc_HARD_COMBINED->SetParameter( 3, -7.46413e+02 );
    fFunc_HARD_COMBINED->SetParameter( 4, -2.14778e+00 );
    fFunc_HARD_COMBINED->SetParameter( 5, -5.02438e+00 );


    // November 9, 2017: try "hard components" separately for PIDs
    fFunc_HARD_PION = new TF1( "fFunc_HARD_PION", "x<[5] ? 0 : [0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])    -   [3]*x*TMath::Exp([4]*x)", 0.4, 20 );
    fFunc_HARD_PION->SetParameter(0, 6.09662   );
    fFunc_HARD_PION->SetParameter(1, 0.12469   );
    fFunc_HARD_PION->SetParameter(2, -0.172349 );
    fFunc_HARD_PION->SetParameter(3, 18.042    );
    fFunc_HARD_PION->SetParameter(4, -5.60966  );
    fFunc_HARD_PION->SetParameter(5, 0.57  );
            // integral from fFuncExp = 0.573339
            // integral from funcPtDiff = 0.0446483


    fFunc_HARD_KAON = new TF1( "fFunc_HARD_KAON", "x<[5] ? 0 : [0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])    -   [3]*x*TMath::Exp([4]*x)", 0.4, 20 );
    fFunc_HARD_KAON->SetParameter(0, 6.34756   );
    fFunc_HARD_KAON->SetParameter(1, 0.120234  );
    fFunc_HARD_KAON->SetParameter(2, -0.462141 );
    fFunc_HARD_KAON->SetParameter(3, 0.732649  );
    fFunc_HARD_KAON->SetParameter(4, -3.14015  );
    fFunc_HARD_KAON->SetParameter(5, 1.45  );
    //        integral from fFuncExp = 0.0743011
    //        integral from funcPtDiff = 0.00184029

    fFunc_HARD_PROTONS = new TF1( "fFunc_HARD_PROTONS", "x<[5] ? 0 : [0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])    -   [3]*x*TMath::Exp([4]*x)", 0.4, 20 );
    fFunc_HARD_PROTONS->SetParameter(0, 9.38768   );
    fFunc_HARD_PROTONS->SetParameter(1, 0.166251  );
    fFunc_HARD_PROTONS->SetParameter(2, -0.713964 );
    fFunc_HARD_PROTONS->SetParameter(3, 0.231222  );
    fFunc_HARD_PROTONS->SetParameter(4, -2.5661   );
    fFunc_HARD_PROTONS->SetParameter(5, 2.28  );
    //        integral from fFuncExp = 0.0351141
    //        integral from funcPtDiff = 0.000208108


    fFunc_HARD_D0 = new TF1( "fFunc_HARD_D0", "x<[5] ? 0 : [0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])    -   [3]*x*TMath::Exp([4]*x)", 0.4, 20 );
    fFunc_HARD_D0->SetParameter(0, 2.78777    );
    fFunc_HARD_D0->SetParameter(1, 0.00184544 );
    fFunc_HARD_D0->SetParameter(2, 251.561    );
    fFunc_HARD_D0->SetParameter(3, 2.679      );
    fFunc_HARD_D0->SetParameter(4, -1.28389   );
    fFunc_HARD_D0->SetParameter(5, 4.0  );
    //        integral from fFuncExp = 1.62524
    //        integral from funcPtDiff = 0.0231375


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
            else if ( fabs( vMother->M() - mD0 ) < 0.01 ) // D meson
            {
                charge  = 0;
                pid  = kPid_D0;
            }

            // fill array with particle
            if ( charge != 0 || pid  == kPid_Lambda || pid  == kPid_phi || pid  == kPid_D0 )
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

            // put PID of hard particle BY HAND:
            double randPid = fRand->Uniform();
            int pid = -1;
            if ( randPid < 0.748542 ) pid = kPid_pion;
            else if ( randPid < 0.748542 + 0.134895 ) pid = kPid_kaon;
            else if ( randPid < 0.748542 + 0.134895 + 0.0876396 ) pid = kPid_proton;
            else if ( randPid < 0.748542 + 0.134895 + 0.0876396 + 0.00913536 ) pid = kPid_phi;
            else if ( randPid < 0.748542 + 0.134895 + 0.0876396 + 0.00913536 + 0.00818593 ) pid = kPid_Lambda;
            else pid = kPid_D0;

            fParticles[nP].pid  = pid;
            nP++;

        }
    }

    //total number of charged particles
    fNparticles = nP;
//        cout << "fNparticles=" << fNparticles << endl;


}

void StringDescr::makeTwoParticlesWithRandomPtEtaPhi(int _pid)
{
    int nP = 0;

    double phiJet = fRand->Uniform( 0, 2*TMath::Pi() );
    int nJets = 2;//( fRand->Uniform() > 0.5 ? 2 : 1 ); // quenching of some fraction of jets

    for ( int iJet = 0; iJet < nJets; iJet++ )
    {
        int nParticlesInJet = 1; //TMath::Max(2,TMath::Nint( fRand->Gaus(4.,0.5) ));
        if (iJet == 1 )
        {
            phiJet += TMath::Pi() + fRand->Gaus(0,0.05);
            FixAngleInTwoPi(phiJet);
            FixAngleInTwoPi(phiJet);
            nParticlesInJet = 1;//TMath::Max(1,TMath::Nint( fRand->Gaus(1.5,0.5) )); //fewer particles due to quenching
        }
        double etaJet = fRand->Uniform( -2.5,2.5); //kMinEta,kMaxEta);

        for ( int iP = 0; iP < nParticlesInJet; iP++ )
        {
//            vMother->SetPtEtaPhiM( motherPt, motherEta, motherPhi, mPions );
            //        fParticles[nP].eta = vMother->Eta();
            //        fParticles[nP].phi = vMother->Phi();
            //        fParticles[nP].pt  = vMother->Pt();
            fParticles[nP].eta = etaJet + fRand->Gaus(0,0.15);//0.2);
            double phi = phiJet;// + fRand->Gaus(0,0.12);//15);
            //FixAngleInTwoPi(phi);

            double pt = 1.; // TMP!
            if ( _pid == kPid_pion )        pt = fFunc_HARD_PION->GetRandom();
            else if ( _pid == kPid_kaon )   pt = fFunc_HARD_KAON->GetRandom();
            else if ( _pid == kPid_proton)  pt = fFunc_HARD_PROTONS->GetRandom();
            else if ( _pid == kPid_D0 )      pt = fFunc_HARD_D0->GetRandom();

            // put PID of hard particle BY HAND:
//            double randPid = fRand->Uniform();
//            int pid = -1;
//            if ( randPid < 0.748542 ) pid = kPid_pion;
//            else if ( randPid < 0.748542 + 0.134895 ) pid = kPid_kaon;
//            else if ( randPid < 0.748542 + 0.134895 + 0.0876396 ) pid = kPid_proton;
//            else if ( randPid < 0.748542 + 0.134895 + 0.0876396 + 0.00913536 ) pid = kPid_phi;
//            else if ( randPid < 0.748542 + 0.134895 + 0.0876396 + 0.00913536 + 0.00818593 ) pid = kPid_Lambda;
//            else pid = kPid_D0;

            fParticles[nP].pid  = _pid; //pid;




//            double pt = fFuncPowerLawPt->GetRandom(); //fRand->Exp(2.);
            //double pt = fFunc_HARD_COMBINED->GetRandom(); //fRand->Exp(2.);
//            cout << "pt = " << pt << endl;
            fParticles[nP].phi = phi;
            fParticles[nP].pt  = pt;
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
