#ifndef StringFragmentation_H
#define	StringFragmentation_H

#define N_STRING_CUTS 1000


class TRandom;
class TLorentzVector;
//class TH1D;
class TF1;
//class TH2D;

const double mRho = 0.775;// fRand->Gaus(0.77,0.05);
const double mRhoWidth = 0.16;
const double mPion = 0.1395;
const double mKaon = 0.494;
const double mProton = 0.938;

const double mPhi = 1.020;
const double mLambda = 1.115;


class StringFragmentation
{
public:
    StringFragmentation();
    void setRandomGenerator(TRandom* rand) {  fRand = rand; }
    void setStringEndPoinsY( double yMin, double yMax, double yStringEndSigma );

    int decayStringIntoParticles( TLorentzVector *vArr, double fictionRhoPt );
    int probabilityChargePlusMinusNeutral();
    int probabilityChargePlusMinus();

private:
    TF1 *funcStringDecay;
    TF1 *funcPt;
    TF1 *funcPtBoltzmanLikePion; // from OnTheFlyDoc.pdf
    TF1 *funcPtBoltzmanLikeKaon; // from OnTheFlyDoc.pdf
    TF1 *funcPtBoltzmanLikeProton; // from OnTheFlyDoc.pdf
//    double stringDecay(Double_t *x, Double_t *par);
    TRandom *fRand;

    //string ends
    double fYmin;
    double fYmax;
    double fStringShiftSigma;

    //string cuts arrays
    double yBreakPoints[N_STRING_CUTS];
    double yBreakPointsSorted[N_STRING_CUTS];
    int indecesCutsSorted[N_STRING_CUTS];

    double breakPointPt[N_STRING_CUTS];
    double breakPointPhi[N_STRING_CUTS];
    double breakPointType[N_STRING_CUTS]; //u,d,s,diquark?

    inline void FixAngleInTwoPi( double &lPhi )
    {
        if ( lPhi > 2 * TMath::Pi() )
            lPhi -= 2 * TMath::Pi();
        else if ( lPhi < 0 )
            lPhi += 2 * TMath::Pi();
    }

};


#endif	/* StringFragmentation_H */
