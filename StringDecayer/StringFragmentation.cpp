//#include "../../commonTools/Tools.cxx"
#include "TRandom.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "DecayInTwo.h"
#include "TLatex.h"

#include "StringFragmentation.h"


const int nAnalysers = 1;



StringFragmentation::StringFragmentation() :
    fRand(0x0)
  , fYmin(-2)
  , fYmax(2)
  , fStringShiftSigma(0)
{
    funcPt = new TF1( "myFuncPt", "x * exp(-x/1.)", 0.0, 20.); //2*TMath::Pi() );
    //    funcPt = new TF1( "myFuncPt", "exp(-x/0.6)", 0., 20.); //2*TMath::Pi() );
    //funcPt->SetParameter( 0, 1. );
    //        funcPt->SetParameter(1,1.05);
    //        funcPt->SetParameter(2,7.92);


}

double stringDecay(Double_t *x, Double_t *par)
{
    //return 1;
    double value = 0;
    double sigma = par[1];
    if ( x[0] > par[0]/2 )
        value = TMath::Exp( -0.5 * ( (x[0]-0)*(x[0]-0)/sigma/sigma ) );
    else if ( x[0] < -par[0]/2 )
        value = TMath::Exp( -0.5 * ( (x[0]+0)*(x[0]+0)/sigma/sigma ) );
    else
        value = 1;
    return value;
}


void StringFragmentation::setStringEndPoinsY(double yMin, double yMax, double yStringShiftSigma )
{
    fYmin = yMin;
    fYmax = yMax;
    fStringShiftSigma = yStringShiftSigma;

    funcStringDecay = new TF1( "funcStringDecay", stringDecay, -8, 8, 2 );
    funcStringDecay->SetParameter(0,1.0); //string plateau length
    funcStringDecay->SetParameter(1,0.1); //gauss ends sigma

}

int StringFragmentation::decayStringIntoParticles( TLorentzVector *vArr, double fictionRhoPt )
{
    if ( !fRand )
    {
        printf("StringFragmentation: fRand is not set!!!\n" );
        return 0;
    }
    //    double yStringSize = fabs( fRand->Gaus(1,0.2) ) + fabs( fRand->Gaus(0,fYmax-fYmin) );
    //    double yStringShift = fRand->Gaus(0,fStringShiftSigma);

    //    double yStringSize = fabs( fRand->Gaus(4,0.4) ) + 2*fabs( fRand->Gaus(0,2) );
    //    double yStringSize = fabs( fRand->Gaus(4,0.2) ) + 2*fabs( fRand->Gaus(0,4) );
    //WAS USED FOR PROCEEDING...    double yStringSize = 10;//fRand->Uniform(2,10);
    double yStringSize = fYmax-fYmin;//fRand->Uniform(2,10);
    //    double yStringSize = fabs( fRand->Gaus(fYmax-fYmin,0.2) ) + 2*fabs( fRand->Gaus(0,0.5) );
    //    double yStringSize = fYmax-fYmin;//fabs( fRand->Gaus(3.6,0.2) ) + 2*fabs( fRand->Gaus(0,5) );
    //    double yStringShift = fRand->Uniform(-fStringShiftSigma,fStringShiftSigma) + fRand->Gaus(0,fStringShiftSigma);
    double yStringShift = (fYmax+fYmin)/2;//fRand->Gaus(0,fStringShiftSigma);
    //    double yStringShift = fRand->Gaus(0,2);
    //    double yStringSize = fabs( fRand->Gaus(3.6,0.2) ) + 2*fabs( fRand->Gaus(0,5) );
    //    double yStringSize = fabs( fRand->Gaus(3.6,0.2) ) + 2*fabs( fRand->Gaus(0,5) );
    //        double yStringSize = fabs( fRand->Gaus(3.6,0.2) ) + 2*fabs( fRand->Gaus(0,5) );
    //    double yStringSize = fRand->Uniform(0., fYmax-fYmin);
    //        double yStringSize = fabs( fRand->Gaus(2.,0.001) ) + 2*fabs( fRand->Gaus(0,0.2) );



    //    funcStringDecay->SetParameter( 0, yStringSize );
    //    funcStringDecay->SetParameter( 1, 1 );

    //    cout << nParticlesInString << endl;
    if(0)if ( fRand->Uniform() > 0.25) //ministrings
    {
        //        yStringSize = TMath::Max(4., fabs(fRand->Gaus(0,2)) );
        yStringSize = fabs( fRand->Gaus(2,0.2) ) + 2*fabs( fRand->Gaus(0,1) );
        //        yStringSize = TMath::Max(0.5, fabs(fRand->Gaus(0,1)) );
        //TMath::Max( 0.5, fabs( fRand->Gaus(1.5,0.3) ) + 2*fabs( fRand->Gaus(0,0.2) ) );
    }
    int nParticlesInString = 0.72 /*coeffToTuneMult*/ * TMath::Max(1,TMath::Nint( fRand->Gaus(yStringSize,yStringSize/10) )); // /2270*1550;
    nParticlesInString *= 1.14; //for energy-dependence! (30.01.2015, tuning basing on STAR)
    //nParticlesInString; // to take into account two pions from single rho!!!!!!!
    //    int nParticlesInString = TMath::Max(1,TMath::Nint( fRand->Gaus(yStringSize,yStringSize/5) ));
    //WAS USED FOR PROCEEDING...    int nParticlesInString = TMath::Max(2,TMath::Nint( fRand->Gaus(1.2,0.5) ));

    int nCutPoints = nParticlesInString + 1;
    for ( int iCut = 0; iCut < nCutPoints; iCut++ )
    {
        //        double y = fRand->Uniform( fRand->Gaus(fYmin,0.2)+yStringShift, fRand->Gaus(fYmax,0.2)+yStringShift );
        //        double y = fRand->Uniform( fYmin+yStringShift, fYmax+yStringShift );
        //        double y = funcStringDecay->GetRandom();
        double y = fRand->Uniform( -yStringSize/2 + yStringShift, yStringSize/2+yStringShift );
        //        double y = fRand->Gaus( 0, fRand->Gaus(fYmin,0.2) ) + yStringShift;
        //y += yStringShift; //fRand->Gaus(0,2);
        yCutPoints[iCut] = y;
    }

    //sort cut points
    TMath::Sort<double, int>( nCutPoints, yCutPoints, indecesCutsSorted, kFALSE );

    //fill array with sorted y of the string cuts
    for ( int iCut = 0; iCut < nCutPoints; iCut++ )
        yCutPointsSorted[iCut] = yCutPoints[ indecesCutsSorted[iCut] ];

    //    double particleMass = mRho;
    //        double particleMass = mPion;
    // LOGICAL ERROR HERE!!! :
    double particleMass = fRand->Gaus(mRho,mRhoWidth/2);//( fRand->Uniform(0,1) > 0.5 ? fRand->Gaus(mRho,mRhoWidth/2) : mPion );
    //double factorToPtExp = ( particleMass == mPion ? 2.5 : 1.57 );
    //make pTs at string cuts
    for ( int iCut = 0; iCut < nCutPoints; iCut++ )
    {
        // !!! use some numerical factor to MATCH MEAN PT when later merge two quarks of string fragments
        //        if ( particleMass == mPion )
        //            cutPointPt[iCut] = fRand->Exp( 0.25 ); //0.3); //funcPt->GetRandom();//fRand->Exp(pTtau);
        //        else //rho
        cutPointPt[iCut] = fRand->Exp( 0.75 );//fictionRhoPt ); //0.3); //funcPt->GetRandom();//fRand->Exp(pTtau);
        //        cutPointPt[iCut] = funcPt->GetRandom( /*fictionRhoPt*/ ); //0.3); //funcPt->GetRandom();//fRand->Exp(pTtau);


        //        if ( fictionRhoPt == 0 )
        //            cutPointPt[iCut] = 0;
        //        else
        //            cutPointPt[iCut] = fRand->Gaus(fictionRhoPt,0.01); //funcPt->GetRandom();//fRand->Exp(pTtau);

        //        cutPointPt[iCut] = fRand->Uniform( 0.1, 0.5 );
        cutPointPhi[iCut] = fRand->Uniform( 0, 2*TMath::Pi() );
    }

    //calc kinematic params of "particles"=string fragments
    for ( int iCut = 0; iCut < nParticlesInString; iCut++ )
    {
        double pT1 = cutPointPt[iCut];
        double pT2 = cutPointPt[iCut+1];
        double phi1 = cutPointPhi[iCut];
        double phi2 = cutPointPhi[iCut+1];
        phi2 += TMath::Pi();
        FixAngleInTwoPi(phi2);

        double ptX = pT1*cos( phi1 ) + pT2*cos( phi2 );
        double ptY = pT1*sin( phi1 ) + pT2*sin( phi2 );
        double ptParticle = sqrt( ptX*ptX+ptY*ptY );
        //        double ptParticle = fRand->Exp(0.9);
        if ( ptParticle < 0.01 )
            ptParticle = 0.01;


        // FOR TESTS!!!
        if(0)
        {
            if ( fictionRhoPt == 0 )
                ptParticle = 0.01;
            else
                ptParticle = fRand->Exp(fictionRhoPt);
            //            ptParticle = fRand->Gaus(fictionRhoPt,0.01); //funcPt->GetRandom();//fRand->Exp(pTtau);
        }

        // !!! RANDOMIZE PT for generated particles:
        //        fParticles[iCut].pt  = funcPt->GetRandom();
        //        cout << fParticles[iP].pt << endl;

        double phiVectorSum = asin( ptY/ptParticle );
        if ( ptX < 0 )
            phiVectorSum = TMath::Pi()-phiVectorSum;
        FixAngleInTwoPi(phiVectorSum);
        FixAngleInTwoPi(phiVectorSum);

        //        histPhi->Fill(phiVectorSum);
        //        histPtRho->Fill(ptParticle );
        //        histPtRhoWithWeight->Fill(ptParticle, 1./2/TMath::Pi()/ptParticle );

        double phiParticle = phiVectorSum;
        if(0)
            phiParticle = fRand->Uniform( 0, 2*TMath::Pi() );//phiVectorSum;

        double yParticle = (yCutPointsSorted[iCut] + yCutPointsSorted[iCut+1])/2;

        if(0)
            yParticle = fRand->Uniform( -yStringSize/2, yStringSize/2 );
        //            yParticle = fRand->Uniform( fYmin+yStringShift, fYmax+yStringShift ); //use it to "shuffle" particles in y!
        //        yParticle += fRand->Gaus(0,1);
//        yParticle = fRand->Uniform( -yStringSize/2, yStringSize/2 );
        yParticle = fRand->Uniform( -yStringSize/2+yStringShift, yStringSize/2+yStringShift );

        // prepare lorentz vector
        double mT = sqrt( ptParticle*ptParticle + particleMass*particleMass );
        double pX = ptParticle * cos(phiParticle);
        double pY = ptParticle * sin(phiParticle);
        double pZ = mT*sinh(yParticle);
        vArr[iCut].SetXYZM( pX, pY, pZ, particleMass );
    }

    return nParticlesInString;

}

// TMP but good function for random mother (rho)!
//void prepareMother( TLorentzVector *vMother )
//{
//    double motherPt = fRand->Exp( 0.4 );
//    //        double motherPt = fRand->Gaus( 0.7, 0.1 );
//    double motherPhi = fRand->Uniform( 0, 2*TMath::Pi() );

//    //    double y = fRand->Uniform( -2, 2 );
//    double y = fRand->Uniform( -2, 2 );
//    y += fRand->Gaus(0,2);

//    double mT = sqrt(motherPt*motherPt+mRho*mRho);
//    double pX = motherPt*cos(motherPhi);
//    double pY = motherPt*sin(motherPhi);
//    double pZ = mT*sinh(y);
//    vMother->SetXYZM( pX, pY, pZ, mRho );

//}



int StringFragmentation::probabilityChargePlusMinusNeutral()
{
    //all +,-,0 have probability 1/3
    int chargeSign = 0;

    double randForCharge = fRand->Uniform(-2,2);

    if ( randForCharge < 0 )
        chargeSign = 0; // neutral (probability is same as for charged! quark combinatorics...)
    else if ( randForCharge > 1 )
        chargeSign = 1; // positive
    else
        chargeSign = -1; // negative
    return chargeSign;
}



