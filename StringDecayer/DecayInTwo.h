#ifndef RhoDecay_H
#define	RhoDecay_H

//struct SimpleTrack;
class TLorentzVector;
class TRandom;
class TVector3;

class DecayInTwo
{
public:
    DecayInTwo();
    void setRandomGenerator(TRandom* rand) { fRand = rand; }
//    void generateDecayProducts( const SimpleTrack *trMother, SimpleTrack **trProducts, const TVector3 *motherBoost = 0x0 );
    void generateDecayProducts( const TLorentzVector &trMother, TLorentzVector &v1, TLorentzVector &v2 );

private:
    void generateDecayProductsInRestFrame(TLorentzVector &v1, TLorentzVector &v2, const double &motherMass );

    TRandom *fRand;

//    TLorentzVector *lLorV_mother;
//    TLorentzVector *lLorV_products[2];
};


#endif	/* RhoDecay_H */
