#include "TRandom.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "DecayInTwo.h"

//const double mMother = 0.77;
const double mProd = 0.14;

DecayInTwo::DecayInTwo()
{
}

void DecayInTwo::generateDecayProductsInRestFrame(TLorentzVector &v1, TLorentzVector &v2, const double &motherMass )
{
    double momProd = sqrt( motherMass*motherMass/2/2 - mProd*mProd );

    double thetaProd = TMath::ACos( fRand->Uniform( -1, 1 ) ); //fRand->Uniform(0,2*TMath::Pi());
    double phiProd = fRand->Uniform( 0, 2*TMath::Pi() );

    v1.SetPtEtaPhiM(momProd,0,0, mProd);
    v1.SetPhi(phiProd);
    v1.SetTheta(thetaProd);

    v2.SetPtEtaPhiM(momProd,0,0, mProd);
    v2.SetPhi( phiProd + TMath::Pi() );
    v2.SetTheta( TMath::Pi() - thetaProd );

    //    cout << "energy of pions: " << v1.Energy() << " " << v2.Energy() << endl;
}

void DecayInTwo::generateDecayProducts( const TLorentzVector &vMother, TLorentzVector &v1, TLorentzVector &v2 )
{
    //##### pair of particles from mother-particle decay in their rest frame
    generateDecayProductsInRestFrame( v1, v2, vMother.M() );

    //##### boost decay products!
    TVector3 vBoost = vMother.BoostVector();
//    cout << vBoost.X() << " " << vBoost.Y() << " " << vBoost.Z() << endl;
    if ( vBoost.Mag() > 0 )
    {
        v1.Boost( vBoost );
        v2.Boost( vBoost );
    }
}

