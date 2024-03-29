/*
 * File:   NucleiCollision.h
 * Author: altsybee
 *
 * Created on 25 Февраль 2012 г., 11:15
 */

#ifndef NucleiCollision_H
#define	NucleiCollision_H

#include "TString.h"

//#define MAX_PARTONS 10000
//#define MAX_PARTONS 208 // for Pb
//#define MAX_PARTONS 10 // partons in toy proton

class TF1;
class TH1D;
class TH2D;
class TH2I;
class TCanvas;
class MinDistanceFinder;
class TRandom;
//class TRandom3;

enum en_nucleus_type
{
    nucleus_proton = 0,
    nucleus_Pb = 1,
    nucleus_Au = 2,
    nucleus_p_Pb = 3,
};

struct Nucleus
{
    int nNucleons; //n partons in nucleus
    int nPartons; //n partons in nucleus
    int maxNofPartons; //max n partons in nucleus, to allocate arrays

    float bx;
    float by;

    float *nX; // x array nucleons in nucleus
    float *nY; // y array nucleons in nucleus
//    float *nZ; // z array nucleons in nucleus
//    float *nTheta; // theta angle of nucleons in nucleus
    int *nuclWounded; // array of flags for wounded nucleons

    float *pX; // x array partons in nucleus
    float *pY; // y array partons in nucleus
    int *pNid; // which nucleon each parton is assigned to
    int *pBusy; // (TMP?) 1.12.2014: array of flags for partons which are in interactioт with other nucleus

    //May 2018: for WQM mode:
    int *isWoundedParton;

    Nucleus(int nNucleons_, int maxNofPartons_)
    {
        bx = 0;
        by = 0;
        nNucleons = nNucleons_;
        nPartons = 0;
        nX = new float[nNucleons];
        nY = new float[nNucleons];
//        nZ = new float[nNucleons];
//        nTheta = new float[nNucleons];
        nuclWounded = new int[nNucleons];

        maxNofPartons = maxNofPartons_;
        pX = new float[maxNofPartons];
        pY = new float[maxNofPartons];
        pNid = new int[maxNofPartons];
        pBusy = new int[maxNofPartons];

        isWoundedParton = new int[maxNofPartons];
    }
    virtual ~Nucleus()
    {
        delete [] nX;
        delete [] nY;
//        delete [] nZ;
//        delete [] nTheta;
        delete [] nuclWounded;

        delete [] pX;
        delete [] pY;
        delete [] pNid;
    }
    void reset()
    {
        for ( int i = 0; i < nPartons; i++ )
            isWoundedParton[i] = 0;
        nPartons = 0;
        for ( int i = 0; i < nNucleons; i++ )
            nuclWounded[i] = 0;
    }
    int getNumberOfWoundedNucleons()
    {
        int nWounded = 0;
        for ( int i = 0; i < nNucleons; i++ )
            if ( nuclWounded[i] > 0 )
                nWounded++;
        return nWounded;
    }
    int getNumberOfWoundedOnlyOnceNucleons()
    {
        int nWoundedOnlyOnce = 0;
        for ( int i = 0; i < nNucleons; i++ )
            if ( nuclWounded[i] == 1 )
                nWoundedOnlyOnce++;
        return nWoundedOnlyOnce;
    }

};

class NucleiCollision {
public:
    NucleiCollision();
    NucleiCollision(const NucleiCollision&);
    virtual ~NucleiCollision();
    
    void setRandomGenerator(TRandom* rand) { fRand = rand; }
    void setOutputDirectoryName(TString strDirName) { fOutDirName = strDirName; }

    void initDataMembers();
    void setEventId(int eventId) { fEventId = eventId; }
//    void changePartonState(int &state, bool hardScatteringFlag = false );
    void changePartonPairState( int &stateA, int &stateB );

    //##### event construction
    void setNucleusType(en_nucleus_type nucleusType) { fNucleusType = nucleusType; }
//    void setNumberOfNucleons(int nNucleons) { fNumberOfNucleons = nNucleons; }
    void setMaxNumberOfPartons(int nPartons) { fMaxPartons = nPartons; }
    void setMeanNofPartonsInNucleon(float n) { fMeanNofPartonsInNucleon = n; }
    void setNucleonGaussianRadius(float r) { fNucleonGaussianRadius = r; }

    void setMeanMultFromOneStringForFictiveMultStudies(float avMult) { fAvMultFromOneStringPerEtaUnitFictive = avMult; }

    float getRandomEventPlanePhi() { return fRandomEventPlane; }

    void setPartonInteractionDistance(float dist) { fPartonInteractionDistance = dist; }
//    void setClusterFormationDistance(float dist) { fClusterFormationDist = dist; }
//    void setStringInteractionRadius(float dist) { fStringInteractionRadius = dist; }
//    void setStringInteractionParA(float parA) { fStringInteractionParA = parA; }
//    void setStringOverlapEnergyDensity(float energy) { fStringOverlapEnergyDensity = energy; }

//    void setHardScatteringProbability(float prob) { fHardScatteringProbability = prob; }


//    void setImpactParameterByHand_0_100(float bImpByHand) { fImpactParameterByHand_0_100 = bImpByHand; }
//    void setImpactParameterRange_0_100(float min, float max) { fImpactParameter_rangeMin_0_100 = min; fImpactParameter_rangeMax_0_100 = max; }
    void setImpactParSpecification(int par) { fImpactParSpecification = par; }
    void setImpactParameterByHand(float bImpByHand) { fImpactParameterByHand = bImpByHand; }
    void setImpactParameterRange(float min, float max) { fImpactParameter_rangeMin = min; fImpactParameter_rangeMax = max; }

    void setWriteEventViewCanvas(int flag) { fFlagWriteEventViewCanvases = flag; }
//    void setComputeStringRepulsion(int flag) { fFlagComputeStringRepulsion = flag; }

    void setFlagUsePoissonianNpartonsFluctuations(int flag) { fFlagUsePoissonianNpartonsFluctuations = flag; }
    void setFlagOnlyOneInteractionPerParton(int flag) { fFlagOnlyOneInteractionPerParton = flag; }
    void setFlagConsiderWoundedPartonsAsStrings(int flag) { fFlagConsiderWoundedPartonsAsStrings = flag; }

    void buildEvent(); //int flagSpecImpactPar=0);
    void finalActions();

    //##### drawings
    void drawStatisticHists();
    void drawEventStructure();

    //Getters
    float getImpactParameterByHand() { return fImpactParameterByHand; }
    float getImpactParameter() { return fImpactParameter; }
    int getMaxNpartons() { return fMaxPartons; }

    int getNstrings() { return fNumberOfStrings; }

    int getNumberOfWoundedPartonsA() { return fNumberOfWoundedPartonsA; }
    int getNumberOfWoundedPartonsB() { return fNumberOfWoundedPartonsB; }
    int getNumberOfWoundedOnlyOncePartonsA() { return fNumberOfWoundedOnlyOncePartonsA; }
    int getNumberOfWoundedOnlyOncePartonsB() { return fNumberOfWoundedOnlyOncePartonsB; }

    int getNumberOfWoundedNucleonsA() { return fA->getNumberOfWoundedNucleons(); }
    int getNumberOfWoundedNucleonsB() { return fB->getNumberOfWoundedNucleons(); }
    int getNumberOfWoundedOnlyOnceNucleonsA() { return fA->getNumberOfWoundedOnlyOnceNucleons(); }
    int getNumberOfWoundedOnlyOnceNucleonsB() { return fB->getNumberOfWoundedOnlyOnceNucleons(); }


    double getNu();
    bool isMBcollision() { return fFlagHaveMBcollision; }
    void calcEccentricity();
    void fillNumberOfStringsInXY();
    void getIngredientsForEccentricity(Nucleus *nucl, Nucleus *nucl2, int &nWounded, float &num, float &denom );

    int getMultFromAllStringsFictiveV0() { return fMultFromAllStringsFictiveV0; }
    int getMultFromAllStringsFictiveMidEta() { return fMultFromAllStringsFictiveMidEta; }


//    float*  getArrStringX() const { return fXstring; }
//    float*  getArrStringY() const { return fYstring; }
//    float*  getArrStringBoostMagn() const { return fStringBoostMagn;   }
//    float*  getArrStringBoostAngle() const { return fStringBoostAngle; }
//    float  getStringBoostMagn(int stringId) const { return fStringBoostMagn[stringId];   }
//    float  getStringBoostAngle(int stringId) const { return fStringBoostAngle[stringId]; }

    float  getStringXminusBover2(int stringId) const { return fXstring[stringId]-fImpactParameter/2;   }
    float  getStringY(int stringId) const { return fYstring[stringId];   }
    float  getDistanceBetweenPartonsForString(int stringId) const { return fDistanceBetweenPartonsForString[stringId];   }
    float  getStringRadiusVectorAngle(int stringId) const { return fStringRadiusVectorAngle[stringId];   }
    float  getStringOrigin(int stringId) const { return fStringOrigin[stringId];   }

    //    bool*   isStringInInteractionArr() const { return fFlagStringInInteraction; }
//    bool isHardInteractionString(int stringId) const { return fFlagStringIsHardInteraction[stringId]; }

    float  getCrossSection() const { return fCrossSectionFinal;   }
    float  getMeanNstrings() const { return fMeanNstringsFinal;   }

    TH2I * getHist2DNstringsXY_thisEv_CoreCorona() { return fHist2DNstringsXY_thisEv_CoreCorona; }

private:
    //##### event construction
    void createNucleus(int nucleusId, Nucleus *nucl, float bx, float by ); //(float *nuclX, float *nuclY, float bx=0, float by=0, float *x, float *y ); //x, y - tmp, potom ubrat'!!!
    void createPartonsOld(float *x, float *y, float bx=0, float by=0 ); //bx, by - used only for backward capability!
    void createPartons( Nucleus *nucl, int nId ); //(float nucleonX, float nucleonY, float *x=0, float *y=0 ); //bx, by - used only for backward capability!
    void createNucleiPair();
    void createStrings();
//    void startClusterSearch();
//    void findStringClusters(int, int& nLinked, int clusterId);
//    void createForcesInsideCluster(int clusterId);
//    void createStringRepulsion();
//    void createStringRepulsionOld();

    //##### drawings
    void drawPartons();
    void drawStrings();
//    void drawStringInteractions();
//    void drawForcesInsideClusters();
//    void drawStringBoosts();

    //##### event construction helpers
    MinDistanceFinder *fPartonInteractionsFinder;
    MinDistanceFinder *fStringInteractionsFinder;
    TRandom *fRand;
    TF1 *funcNucleonDensity;
    TF1 *funcPartonDensity;
//    TF1 *funcStringInteractionDist;
    TF1 *funcImpactPar1D;
    bool fFlagDataMembersInitialized;
    bool fFlagHaveMBcollision;
    bool fFlagWriteEventViewCanvases;
//    bool fFlagComputeStringRepulsion; //if true - make clusters and repulsion
    int fEventId; // number of event, used for output file names
    float fRandomEventPlane; // to be used outsice

    //##### model parameters
    en_nucleus_type fNucleusType; //type of nucleus
    int fNumberOfNucleons; //n nucleons in nucleus
    int fMaxPartons; //max partons in the structure (nucleus or proton or whatever)
    double fNucleusRadius; //nucleus radius, fm
    double fNucleusWSa; //Woods-Saxon parameter a, fm
//    double fRadiusParton; //proton radius
    //    double fRstring; //string radius
    double fPartonInteractionDistance; //distanse when partons interact, fm
//    double fStringInteractionRadius; //distanse when strings interact, fm
//    double fStringOverlapEnergyDensity; //overlap energy density per fm2 of area
//    double fStringInteractionParA;  //parameterization a-la Woods-Saxon, fm
//    double fClusterFormationDist; //distance b/n strings to form a cluster, fm
    bool fFlagUsePoissonianNpartonsFluctuations; // by default - YES // ADDED ON MAY 9, 2018!
    bool fFlagOnlyOneInteractionPerParton; // by default - YES, in NO - should be like in WQM! // ADDED ON MAY 10, 2018
    bool fFlagConsiderWoundedPartonsAsStrings; // by default - NO

    double fMeanNofPartonsInNucleon; //mean n of partons in nucleon
    double fNucleonGaussianRadius; // nucleon radius, fm

//    float fHardScatteringProbability; //probability for a string to be a hard partonic collision

    float fImpactParameter; //impact parameter, fm
//    float fImpactParameterByHand_0_100; //impact parameter, 0-100
//    float fImpactParameter_rangeMin_0_100; //impact parameter, 0-100
//    float fImpactParameter_rangeMax_0_100; //impact parameter, 0-100
    float fImpactParameterByHand; //impact parameter, fm
    float fImpactParameter_rangeMin; //impact parameter, fm
    float fImpactParameter_rangeMax; //impact parameter, fm
    int fImpactParSpecification; //impact parameter tuning (MB, set presicely, or in range)

    Nucleus *fA; //nucleus A
    Nucleus *fB; //nucleus B

//    StringCluster **fStringClusters; // array of StringCluster struct

    bool **matrixNcoll; //matrix of binary nucleon collisions

//    int fNpartonsInA; //n partons in nucleus A
//    int fNpartonsInB; //n partons in nucleus B

//    float *fNuclX1; //x array nucleons in nucleus 1
//    float *fNuclY1; //y array nucleons in nucleus 1
//    float *fNuclX2; //x array nucleons in nucleus 2
//    float *fNuclY2; //y array nucleons in nucleus 2

//    float *fX1; //x array partons in nucleus 1
//    float *fY1; //y array partons in nucleus 1
//    float *fX2; //x array partons in nucleus 2
//    float *fY2; //y array partons in nucleus 2

    //##### calculated event properties
    float *fXstring; //strings position, x coord array
    float *fYstring; //strings position, y coord array
    float *fDistanceBetweenPartonsForString; //what was the distance between partons which give this string
    float *fStringRadiusVectorAngle; //remember strings radius vector angles wrt (0,0)
    short *fStringOrigin; //from valence quarks? sea quarks/gluons?..

//    float *fXstringInteraction; //strings interaction point, x coord array
//    float *fYstringInteraction; //strings interaction point, y coord array
//    float *fStringIntDist;      //distances btwn pairs of interacting strings
//    float *fStringIntAngles;    //angles of string pairs (phi), in radians
//    bool *fFlagStringInInteraction; //strings in interaction flag
//    int *fNStringsInCluster; //number of strings in cluster, to which current string is linked
//    int *fClusterIdForString; //if string in cluster -> bind cluster id
//    bool *fFlagStringIsHardInteraction; //strings is a hard interaction (for jet simulation)

//    int **fClusterIdsWithLinkedStringIds; // 2D array: 1-cluster id, 2- string ids in this cluster

//    float *fStringBoostMagn;
//    float *fStringBoostAngle;

    int fNumberOfStrings;   //n of strings in the event
//    int fNumberOfClusters;   //n of clusters in the event
    int fNumberOfStringInteractions;    //n of interacting string pairs in the event

    int fNumberOfWoundedPartonsA;   //
    int fNumberOfWoundedPartonsB;   //
    int fNumberOfWoundedOnlyOncePartonsA;   //
    int fNumberOfWoundedOnlyOncePartonsB;   //

    float fEccentricity;

    int fEvTrialsSuccess;
    int fEvTrialsFailed;

    float fCrossSectionFinal;  // cross section after calculation of all events
    float fMeanNstringsFinal;  // mean n stringsafter calculation of all events

    int fNparticipants; // for nu calculation
    int fNcollisions; // for nu calculation

    float fAvMultFromOneStringPerEtaUnitFictive;  // Nov 2017 - Feb 2018: fictive multiplicities from strings (to control mult distr and for centrality determination)
    int fMultFromAllStringsFictiveV0;  // Nov 2017: fictive multiplicities from strings (to control mult distr and for centrality determination)
    int fMultFromAllStringsFictiveMidEta;  // Feb 2018

    // ##### output directory name
    TString fOutDirName;

    //##### histos
    TCanvas *fCanvEventView;
    TCanvas *fCanvEventStatistics;

    TH1D *fHistNstrings; // n of strings distribution
    TH1D *fHistStringInteractions; // n of string intersections distribution
    TH1D *fHistNucleonsR; // nucleon r pos distribution
    TH1D *fHistPartonsR; // partons r pos distribution
    TH1D *fHistBusyPartonsR; // test how many partons become busy after interaction as function of radius
    TH1D *fHistPartonsValenceR; // "valence partons" radius
    TH1D *fHistBusyPartonsValenceR; // test how many valence partons become busy after interaction as function of radius

    TH1D *fHistMinDistBetweenInteractingPartons; // (min) distances b/n interecting partons // ADDED: Sept 2017
    TH1D *fHistMinDistBetweenInteractingPartonsSoft; // (min) distances b/n interecting partons // ADDED: Sept 2017
    TH1D *fHistMinDistBetweenInteractingPartonsValence; // (min) distances b/n interecting partons // ADDED: Sept 2017

    TH1D *fHistStringPositionRadius; // position of the string: radius
    TH1D *fHistImpactParameter; // impact parameter distribution
    TH1D *fHistNucleonEccentricity; // hist with eccentricities
    TH1D *fHistImpParNucleonsFromAB_central    ; // to test impact parameter distibution for nucleon from A wrt nucleons from B
    TH1D *fHistImpParNucleonsFromAB_semicentral; // to test impact parameter distibution for nucleon from A wrt nucleons from B
    TH1D *fHistImpParNucleonsFromAB_peripheral ; // to test impact parameter distibution for nucleon from A wrt nucleons from B

    TH2D *fHist2DNstringsXY_thisEvent; //2d plot n strings in xy (this event)
    TH2D *fHist2DNstringsXY; //2d plot n strings in xy
    TH2D *fHist2DN2stringsXY; //2d plot n2 strings in xy
    TH2D *fHist2DSigmaNstringsXY; //2d plot sigma n strings in xy

    TH2I *fHist2DNstringsXY_thisEv_CoreCorona; //2d plot n strings in xy (this event) - May 2018 - for Core-corona study

    TH1D *fHistStringRadiusVectorPhi; //radius-vectors phi-s of strings wrt (0,0)

    //string interactions QA plots
//    TH1D *fHistStringInteractionsPhi; //phi directions of the interacting string pairs
//    TH1D *fHistStringInteractionsDistance; //distances b/n string pairs
//    TH1D *fHistStringInteractionsMagnitude; //distances b/n string pairs

    TH2D *fHist2DNInteractionsVsNstrings; //2d plot n strings vs n interactions b/n them
    TH2D *fHist2DNInteractionsVsImpactPar; //2d plot impact par vs n interactions

    //force on strings
//    TH1D *fHistForcesInClustersX; // forces acting on string, x
//    TH1D *fHistForcesInClustersY; // forces acting on string, y
//    TH1D *fHistForcesInClustersMagn; // forces acting on string, magn

//    TH1D *fHistForcesInClustersAbsX; // forces acting on string, abs x
//    TH1D *fHistForcesInClustersAbsY; // forces acting on string, abs y
//    TH1D *fHistForcesInClustersRatioAbsXY; // forces acting on string, ratio absX/absY

    //string clusters
//    TH1D *fHistStringsInClusters; //number of strings in clusters
//    TH2D *fHist2DStringsInClustersVsB; //>=some number of strings in clusters vs impact par


    TH2D *fHist2DStringsSVsImpactS; //Nstrings*S0/S_impact
    TH1D *fHistNStringsInLargestCluster; //number of strings in largest cluster
    TH1D *fHistNucleonDensityR; //nucleon r density hist
    TH1D *fHistNumberWoundedNucleons; //number of wounded nucleons
    TH1D *fHistNcoll; //number of binary nucleon collisions

    // Nov 2017: fictive multiplicities from strings (to control mult distr and for centrality determination)
    TH1D *fHistFictiveMultDistrV0; // mult distribution in V0
    TH1D *fHistFictiveMultDistrMidEta; // mult distribution at mid rap

    //##### visualisation parameters
    double fVisNucleusRadiusNucleus; //visual size nucleus
    //    double rVisParton; //visual nucleon

};

#endif	/* NucleiCollision_H */

