#ifndef STRING_H
#define	STRING_H

//#define MAX_PARTICLES   10000

#define MAX_NJETS   1000
#define N_STRING_CUTS 1000

#include "TLorentzVector.h"
#include "TVector3.h"
#include "ParticleDescr.h"
//#include "Jet.h"
#include "DecayInTwo.h"


class TH1D;
class TRandom;
class TLorentzVector;
class TF1;
class TVector3;

class DecayInTwo;
class StringFragmentation;


class StringDescr : public ParticleArr
{
public:
    StringDescr();
    void setRandomGenerator(TRandom* rand);
//    void setCoeffPtKickPerUnitMagn(float ptKickPerUnitMagn) { fPtKickPerUnitMagn = ptKickPerUnitMagn; }

    StringFragmentation *getStringFragmentationPointer() { return strFr; }

    //String(const String& orig);
    virtual ~StringDescr();
    void hadronizeString(double stringBeta = 0, double boostPhi = 0 );
    void makeTwoJets();
//    int getNparticles() { return fNparticles; }

//    void setState( short state ) { fStringState = state; }
//    short getState() { return fStringState; }
//    ParticleArr * getJetArray( int jetId ) { return &fJetArray[jetId]; }

//    double etaStep();
//    void calcStringFragmentation();
private:
//    void decayStringIntoParticles( int nParticlesInString, TLorentzVector *vArr, double fictionRhoPt );
//    int probabilityChargePlusMinusNeutral();
    inline void FixAngleInTwoPi( double &lPhi )
    {
        if ( lPhi > 2 * TMath::Pi() )
            lPhi -= 2 * TMath::Pi();
        else if ( lPhi < 0 )
            lPhi += 2 * TMath::Pi();
    }


//    int fNparticles;
    TRandom *fRand;
//    float fPtKickPerUnitMagn; // GeV //apply this kick per unit of magnitude taken from the NuclearStructure model
//    TLorentzVector *fLorVect;
    TF1 *funcPt;
    TF1 *fFuncGausYBoost; //gaussian profile of the boost (set by hand)

    // ##### string cut generation
//    double yCutPoints[N_STRING_CUTS];
//    double yCutPointsSorted[N_STRING_CUTS];
//    int indecesCutsSorted[N_STRING_CUTS];

//    double cutPointPt[N_STRING_CUTS];
//    double cutPointPhi[N_STRING_CUTS];

//    double cutPointEta[N_STRING_CUTS];
//    double cutPointEtaTmp[N_STRING_CUTS];

//    double cutPointPt[N_STRING_CUTS];
//    double cutPointPhi[N_STRING_CUTS];

    //    ParticleArr fJetArray[MAX_NJETS];



    TLorentzVector *vMother;
    TLorentzVector *vMothers;
    TLorentzVector v1, v2;
    double etaArrTmp[N_STRING_CUTS];

    TVector3 stringBoostT;
    TVector3 stringBoostZ;

    //string fragmentation
    StringFragmentation *strFr;
    DecayInTwo *decayInTwo;

};





#endif	/* STRING_H */

