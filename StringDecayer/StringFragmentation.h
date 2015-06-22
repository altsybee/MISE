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
const double mPion = 0.14;

class StringFragmentation
{
public:
    StringFragmentation();
    void setRandomGenerator(TRandom* rand) {  fRand = rand; }
    void setStringEndPoinsY( double yMin, double yMax, double yStringEndSigma );

    int decayStringIntoParticles( TLorentzVector *vArr, double fictionRhoPt );
    int probabilityChargePlusMinusNeutral();

private:
    TF1 *funcStringDecay;
    TF1 *funcPt;
//    double stringDecay(Double_t *x, Double_t *par);
    TRandom *fRand;

    //string ends
    double fYmin;
    double fYmax;
    double fStringShiftSigma;

    //string cuts arrays
    double yCutPoints[N_STRING_CUTS];
    double yCutPointsSorted[N_STRING_CUTS];
    int indecesCutsSorted[N_STRING_CUTS];

    double cutPointPt[N_STRING_CUTS];
    double cutPointPhi[N_STRING_CUTS];

    inline void FixAngleInTwoPi( double &lPhi )
    {
        if ( lPhi > 2 * TMath::Pi() )
            lPhi -= 2 * TMath::Pi();
        else if ( lPhi < 0 )
            lPhi += 2 * TMath::Pi();
    }

};


#endif	/* StringFragmentation_H */
