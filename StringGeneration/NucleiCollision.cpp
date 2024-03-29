#include "NucleiCollision.h"
#include "DistanceEntry.h"

#include "TFile.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "TEllipse.h"
#include "TLine.h"
#include "TArrow.h"
#include "TMath.h"
#include "TRandom.h"
//#include "TRandom3.h"
#include "MinDistanceFinder.h"

#include <fstream>
#include <iostream>
using namespace std;

#define PI 3.1415926
//const float kNucleonSize = 1.12;  // fm
const float kNucleonSize = 0.85; // fm
const float kNucleonParA = 0.02; // fm

const float Mstring = 1.0; // string energy density 1 GeV (from Abromovsky) //0.6; // stringMassOverRapUnit  //GeV
//const float coeffOverlapEnergy = 0.01; //0.125;//0.005; //GeV

//const
double kCoeffRadiusForMB = 1.75;  //tuned to be Ok for MB! (for NN)

const double sizeForXY_for_CoreCorona = 11.5; // fm // Added in May 2018
const int nBinsXY_for_CoreCorona = 23; // May 2018, after discussion with EA


const int nClusterColors = 14;
const int kClusterColors[nClusterColors] = {
    kRed, kRed+2,
    kMagenta, kMagenta+2,
    kOrange+7, kOrange+9,
    kPink+1, kPink+6,
    kViolet, kViolet-1,
    //    kViolet-2, kViolet-3,
    kViolet-3, kViolet-5,
    //    kViolet-6, kViolet-7,
    kViolet-7, kViolet-9,
};

const int kNstringsWhenClustersIsConsideredLarge = 3;

inline void FixAngleInTwoPi( double &lPhi )
{
    if ( lPhi > 2 * TMath::Pi() )
        lPhi -= 2 * TMath::Pi();
    else if ( lPhi < 0 )
        lPhi += 2 * TMath::Pi();
}

int maxElement( int *arr, int size )
{
    //    cout << "size=" << size << endl;
    //    cout << "arr[0]=" << arr[0] << endl;
    //    cout << "arr[1]=" << arr[1] << endl;
    if (size<=0)
        return -1;
    
    int idMax = -1;
    int maxVal = -1;//arr[0];
    for ( int i = 0; i < size; i++ )
    {
        if ( arr[i] > maxVal )
        {
            maxVal = arr[i];
            //            cout << maxVal << endl;
            idMax = i;
        }
    }
    return idMax;
}


double segmentArea(double R, double x)
{
    //    return R*R*TMath::ACos(x/R)-x*sqrt(R*R-x*x);
    double alpha = 2*TMath::ACos(x/R);
    return R*R/2 * (alpha-TMath::Sin(alpha));

}

double circlesIntersectionArea( double d, double r, double R )
{
    d = fabs(d);
    if ( d/R < 0.0000001 ) // if impact par is 0, return area of circle r
        return TMath::Pi() * r * r;
    double x = (d*d-r*r+R*R) / (2*d); //proj of the intersection points on d
    //    double y = sqrt(R*R-x*x); //chord half-length
    return segmentArea(r,x) + segmentArea(R,d-x);
}

//int minIndex( float *arr, int size, int startIndex=0 )
//{
//    if ( size <= 0 || startIndex < 0 || startIndex >= size )
//        return -1;
//    int minId = -1;
//    float minValue = arr[startIndex];
//    for ( int i = startIndex; i < size; i++ )
//    {
//        if ( arr[i] < minValue )
//        {
//            minValue = arr[i];
//            minId = i;
//        }
//    }
//    return minId;
//}


//double rhoNucleon( double r, double R, double a )
//{
//    // (omega is omitted in this formula)
//    const double rho0 = 1.;
//    return rho / (1 + TMath::Exp( (r-R) / a ));
//}


//void NucleiCollision::changePartonState( int &state, bool hardScatteringFlag )
//{
//    //change parton state if the parton participates in string formation

//    //2.12.2014 flag is not used (always FALSE), but may be changed in some future:
//    if ( hardScatteringFlag )
//    {
//        state = 3;
//        return;
//    }

//    if ( state == 0 ) //it's "sea quark" or "gluon"
//        state = 1;
//    if ( state == -1 ) //it's "valence quark"
//        state = 2;
//}
void NucleiCollision::changePartonPairState( int &stateA, int &stateB )
{
    //change parton state if the parton participates in string formation
    if ( stateA == -1 && stateB == -1 ) //both partons are "valence quarks"
    {
        stateA = 2;
        stateB = 2;
    }
    else //one of the quarks is "sea quark" or "gluon"
    {
        stateA = 1;
        stateB = 1;
    }
}


NucleiCollision::NucleiCollision(/*int seed*/) :
    fRand(0x0)
  , funcNucleonDensity(0x0)
  , funcPartonDensity(0x0)
  //  , funcStringInteractionDist(0x0)
  , funcImpactPar1D(0)
  , fFlagDataMembersInitialized(false)
  , fFlagHaveMBcollision(false)
  , fFlagWriteEventViewCanvases(false)
  //  , fFlagComputeStringRepulsion(true)
  , fNucleusType(nucleus_Pb)
  , fNumberOfNucleons(208) //Pb nucleus
  , fMaxPartons(8000)
  , fNucleusRadius(1.) //6.413)//5.5) //0.877) //radius proton, fm
  , fNucleusWSa(1.)
  //  , fRadiusParton(0.1)    //radius "parton", fm
  , fPartonInteractionDistance(0.2)//0.2)
  , fFlagUsePoissonianNpartonsFluctuations(true)
  , fFlagOnlyOneInteractionPerParton(true)
  , fFlagConsiderWoundedPartonsAsStrings(false)
  , fMeanNofPartonsInNucleon(5)
  , fNucleonGaussianRadius(0.4)
  //  , fHardScatteringProbability(0.03)

  , fAvMultFromOneStringPerEtaUnitFictive(1.1) // 1.1 particle per unit of rapidity
  , fMultFromAllStringsFictiveV0(0)
  , fMultFromAllStringsFictiveMidEta(0)

  , fImpactParameter(0.) //fm (is MC-played)
  //  , fImpactParameterByHand_0_100(0.)
  //  , fImpactParameter_rangeMin_0_100(0)
  //  , fImpactParameter_rangeMax_0_100(100)
  , fImpactParameterByHand(0.)
  , fImpactParameter_rangeMin(0.)
  , fImpactParameter_rangeMax(1.)
  , fImpactParSpecification(0) //MB
  , fCanvEventView(0x0)
  , fCanvEventStatistics(0x0)
  , fOutDirName("outputs_NucleiCollision")
{
}

void NucleiCollision::initDataMembers()
{
    fFlagDataMembersInitialized = true;
    fEventId = 0;
    fRandomEventPlane = 0;
    
    switch (fNucleusType)
    {
    case nucleus_proton:
        fNumberOfNucleons = 1;
        fNucleusRadius = 0.877; //6min 08 sec in https://www.youtube.com/watch?v=IN9METK5nYM#t=368
        fNucleusWSa = 0.02; // from ASer
        //        cout << "Pb" << endl;
        kCoeffRadiusForMB = 2.1;
        break;
    case nucleus_Pb:
        fNumberOfNucleons = 208; // !! Centrality p-Pb ALICE: arXiv:1412.6828       //207;//8;
        fNucleusRadius = 6.62;// 6.62 ± 0.06 // Loizides
        fNucleusWSa = 0.546; // 0.546 ± 0.010 // Loizides
        //        cout << "Pb" << endl;
        break;
    case nucleus_Au:
        fNumberOfNucleons = 197;
        fNucleusRadius = 6.38; // 6.38 ± 0.13 // Loizides
        fNucleusWSa = 0.535; // 0.535 ± 0.053 // Loizides
        break;
    case nucleus_p_Pb:
        fNumberOfNucleons = 197;
        fNucleusRadius = 6.38; // 6.38 ± 0.13 // Loizides
        fNucleusWSa = 0.535; // 0.535 ± 0.053 // Loizides
        break;
    default:
        break;
    }
    
    //    fRand = new TRandom3;
    //    fRand->SetSeed(0);
    
    funcNucleonDensity = new TF1( "funcNucleonDensity"
                                  , "4.*3.1415926*x*x / (1 + TMath::Exp( (x-[0]) / [1] ))", 0, fNucleusRadius*2 );
    //    , "1. / (1 + TMath::Exp( (x-[0]) / [1] ))", 0, fNucleusRadius*2 );
    funcNucleonDensity->SetParameter(0, fNucleusRadius);//parameter R
    funcNucleonDensity->SetParameter(1, fNucleusWSa); //parameter a
    //    funcNucleonDensity->SetNpx(500);
    
    funcPartonDensity = new TF1( "funcPartonDensity"
                                 , "4.*3.1415926*x*x / (1 + TMath::Exp( (x-[0]) / [1] ))", 0, kNucleonSize*2 );
    funcPartonDensity->SetParameter(0, kNucleonSize);//parameter R
    funcPartonDensity->SetParameter(1, kNucleonParA); //parameter a
    
    funcImpactPar1D = new TF1( "funcImpactPar1D", "x", 0, 2*fNucleusRadius*kCoeffRadiusForMB ); // !!! tmp!!! make wider!

    fA = new Nucleus( fNumberOfNucleons, fMaxPartons);
    fB = new Nucleus( fNumberOfNucleons, fMaxPartons);
    
    //init 2D array for clusterId->stringsIds
    //    fStringClusters = new StringCluster*[fMaxPartons]; // CAN'T BE MORE THEN fMaxPartons/2 clusters
    //    for( int i = 0; i < fMaxPartons; ++i )
    //        fStringClusters[i] = new StringCluster( fMaxPartons ); // CAN'T BE MORE THEN fMaxPartons/2 strings in one cluster

    //init N coll matrix
    matrixNcoll = new bool*[fNumberOfNucleons];
    for( int i = 0; i < fNumberOfNucleons; ++i )
        matrixNcoll[i] = new bool[fNumberOfNucleons];

    
    fXstring = new float[fMaxPartons];
    fYstring = new float[fMaxPartons];
    fDistanceBetweenPartonsForString = new float[fMaxPartons];
    fStringRadiusVectorAngle = new float[fMaxPartons];
    fStringOrigin = new short[fMaxPartons];

    //    fXstringInteraction = new float[fMaxPartons];
    //    fYstringInteraction = new float[fMaxPartons];
    //    fStringIntDist = new float[fMaxPartons];
    //    fStringIntAngles = new float[fMaxPartons];
    
    //    fFlagStringInInteraction = new bool[fMaxPartons];
    //    fNStringsInCluster = new int[fMaxPartons];
    //    fClusterIdForString = new int[fMaxPartons];
    //    fFlagStringIsHardInteraction = new bool[fMaxPartons];

    //    fStringBoostMagn = new float[fMaxPartons];
    //    fStringBoostAngle = new float[fMaxPartons];
    
    fNumberOfStrings = 0;
    fNumberOfWoundedPartonsA = 0;
    fNumberOfWoundedPartonsB = 0;
    fNumberOfWoundedOnlyOncePartonsA = 0;
    fNumberOfWoundedOnlyOncePartonsB = 0;
    //    fNumberOfClusters = 0;
    fNumberOfStringInteractions = 0;
    fEccentricity = -1;

    fEvTrialsSuccess = 0;
    fEvTrialsFailed = 0;

    fCrossSectionFinal = 0;
    fMeanNstringsFinal = 0;

    fNparticipants = 0;
    fNcollisions = 0;
    
    fHistNucleonsR = new TH1D("fHistNucleonsR","nucleon r position;r, fm;entries",100,0, fNucleusRadius*3 );
    fHistPartonsR = new TH1D("fHistPartonsR","parton r position;r, fm;entries",100,0, fNucleusRadius*3 );
    fHistBusyPartonsR = new TH1D("fHistBusyPartonsR","busy partons r position;r, fm;entries",100,0, fNucleusRadius*3 );
    fHistPartonsValenceR = new TH1D("fHistPartonsValenceR","parton r position;r, fm;entries",100,0, fNucleusRadius*3 );
    fHistBusyPartonsValenceR = new TH1D("fHistBusyPartonsValenceR","busy partons r position;r, fm;entries",100,0, fNucleusRadius*3 );

    fHistMinDistBetweenInteractingPartons = new TH1D("fHistMinDistBetweenInteractingPartons","min distances between interacting partons;r, fm;n entries",400,0, fPartonInteractionDistance*1.2 );
    fHistMinDistBetweenInteractingPartonsSoft = new TH1D("fHistMinDistBetweenInteractingPartonsSoft","min distances b/n interacting partons - SOFT INTERACTION;r, fm;n entries",400,0, fPartonInteractionDistance*1.2 );
    fHistMinDistBetweenInteractingPartonsValence = new TH1D("fHistMinDistBetweenInteractingPartonsValence","min distances b/n interacting VALENCE QUARKS;r, fm;n entries",400,0, fPartonInteractionDistance*1.2 );

    fHistNstrings = new TH1D("fHistNstrings", "n strings;N;entries", 4*2*2*5*2*fNumberOfNucleons+1+20, -0.5, 4*2*2*5*2*fNumberOfNucleons+0.5+20 );
    fHistStringInteractions = new TH1D("fHistStringInteractions", "n string intersections", fMaxPartons+1, -0.5, fMaxPartons+0.5 );
    fHistStringPositionRadius = new TH1D("fHistStringPositionRadius", "string r position/R", 50, 0, 2.5 );
    fHistImpactParameter = new TH1D("fHistImpactParameter", "impact parameter", 50, 0, fNucleusRadius*4 );
    fHistNucleonEccentricity = new TH1D("fHistNucleonEccentricity","nucleon eccentricity;#epsilon;entries", 100, 0, 1 );

    //1.12.2014 //to test impact parameter distibution for nucleon from A wrt nucleons from B
    fHistImpParNucleonsFromAB_central     = new TH1D("fHistImpParNucleonsFromAB_central","imp par for nucleon from A wrt nucleon from b;r, fm;entries",200,0, kNucleonSize*2 );
    fHistImpParNucleonsFromAB_semicentral = new TH1D("fHistImpParNucleonsFromAB_semicentral","imp par for nucleon from A wrt nucleon from b;r, fm;entries",200,0, kNucleonSize*2 );
    fHistImpParNucleonsFromAB_peripheral  = new TH1D("fHistImpParNucleonsFromAB_peripheral","imp par for nucleon from A wrt nucleon from b;r, fm;entries",200,0, kNucleonSize*2 );

    const int nBinsXY = 41;
    const double sizeForXY = 2*fNucleusRadius;
    fHist2DNstringsXY_thisEvent = new TH2D("fHist2DNstringsXY_thisEvent", "nStr_xy_this", nBinsXY, -sizeForXY, sizeForXY , nBinsXY, -sizeForXY, sizeForXY );
    fHist2DNstringsXY = new TH2D("fHist2DNstringsXY", "n strings in xy"                      , nBinsXY, -sizeForXY, sizeForXY , nBinsXY, -sizeForXY, sizeForXY );
    fHist2DN2stringsXY = new TH2D("fHist2DN2stringsXY", "nStr in xy"                    , nBinsXY, -sizeForXY, sizeForXY , nBinsXY, -sizeForXY, sizeForXY );
    fHist2DSigmaNstringsXY = new TH2D("fHist2DSigmaNstringsXY", "#sigma n strings in xy", nBinsXY, -sizeForXY, sizeForXY , nBinsXY, -sizeForXY, sizeForXY );

    fHist2DNstringsXY_thisEv_CoreCorona = new TH2I("fHist2DNstringsXY_thisEv_CoreCorona", "nStr_xy_this_for_CoreCorona"
                                                   , nBinsXY_for_CoreCorona, -sizeForXY_for_CoreCorona, sizeForXY_for_CoreCorona , nBinsXY_for_CoreCorona, -sizeForXY_for_CoreCorona, sizeForXY_for_CoreCorona );


    fHistStringRadiusVectorPhi = new TH1D("fHistStringRadiusVectorPhi", "string #phi angle wrt (0,0);#phi;n strings", 100, 0, 2*TMath::Pi() );

    fHist2DNInteractionsVsNstrings = new TH2D("fHist2DNInteractionsVsNstrings", "n string interactions vs nstrings;n strings;n interactions", fMaxPartons+1, -0.5, fMaxPartons+0.5 , fMaxPartons+1, -0.5, fMaxPartons+0.5 );
    fHist2DNInteractionsVsImpactPar = new TH2D("fHist2DNInteractionsVsImpactPar", "n string interactions;b,fm;n interactions", 50, 0, 2*fNucleusRadius , fMaxPartons+1, -0.5, fMaxPartons+0.5 );

    //    fHist2DStringsSVsImpactS= new TH2D("fHist2DStringsSVsImpactS", "n*s0/S;b/R;n*s0/S", 50, 0, 100 , 100, 0, 2 );
    fHist2DStringsSVsImpactS= new TH2D("fHist2DStringsSVsImpactS", "n*s0/S;b;n*s0/S", 50, 0, 2*fNucleusRadius , 100, 0, 10 );
    
    fHistNStringsInLargestCluster = new TH1D("fHistNStringsInLargestCluster", "N strings in largest cluster;n strings", fMaxPartons, -0.5, fMaxPartons-0.5 );
    
    fHistNucleonDensityR = new TH1D("fHistNucleonDensityR", "Nucleon density;r, fm;#rho (r)", 200, 0, fNucleusRadius*2 );
    
    fHistNumberWoundedNucleons = new TH1D("fHistNumberWoundedNucleons", ";number of wounded nucleons;n events", 5*2*fNumberOfNucleons+1, -0.5, 5*2*fNumberOfNucleons+0.5 );
    fHistNcoll = new TH1D("fHistNcoll", ";N_{coll};n events", 5*2*fNumberOfNucleons+1, -0.5, 5*2*fNumberOfNucleons+0.5 );
    
    // Nov 2017: fictive multiplicities from strings (to control mult distr and for centrality determination)
    fHistFictiveMultDistrV0 = new TH1D("fHistFictiveMultDistrV0", ";n particles;entries", 12001, -0.5, 12000.5 );
    // Feb 2018:
    fHistFictiveMultDistrMidEta = new TH1D("fHistFictiveMultDistrMidEta", ";n particles;entries", 12001, -0.5, 12000.5 );

    //
    fVisNucleusRadiusNucleus = 0.35; //visual size nucleus
    //    rVisParton = fRadiusParton / fNucleusRadius * fVisNucleusRadiusNucleus; //visual nucleon
    
    //event construction helpers
    fPartonInteractionsFinder = new MinDistanceFinder( fPartonInteractionDistance );
    //    fStringInteractionsFinder = new MinDistanceFinder( fStringInteractionRadius ); //NOT USED NOW!!! it's for old implementation of strings interactions
}

NucleiCollision::NucleiCollision(const NucleiCollision& ) {
}

NucleiCollision::~NucleiCollision() {
    //    delete [] fX1;
    //    delete [] fY1;
    //    delete [] fX2;
    //    delete [] fY2;
}

//######### build event
void NucleiCollision::buildEvent()
{
    if ( !fRand )
    {
        cout << "NucleiCollision: fRand IS NOT INITIALIZED!!! exiting..." << endl;
        return;
    }

    if ( !fFlagDataMembersInitialized )
        initDataMembers();

    fA->reset();
    fB->reset();
    //    for( int i = 0; i < fNumberOfClusters; ++i )
    //        fStringClusters[i]->reset();

    //reset N coll matrix
    for( int i = 0; i < fNumberOfNucleons; ++i )
        for( int j = 0; j < fNumberOfNucleons; ++j )
            matrixNcoll[i][j] = 0;
    
    fFlagHaveMBcollision = false;
    fEccentricity = -1;

    fMultFromAllStringsFictiveV0 = 0;
    fMultFromAllStringsFictiveMidEta = 0;

    // generate random event plane (EP) for later use outside
    fRandomEventPlane = fRand->Uniform( 0, 2*TMath::Pi() );
    
    createNucleiPair();
    createStrings();
    
    //count successes/failures for cross-section calculation
    fFlagHaveMBcollision ? fEvTrialsSuccess++ : fEvTrialsFailed++;

    if ( fFlagHaveMBcollision )
    {
        fHistImpactParameter->Fill( fImpactParameter );
        fHistNstrings->Fill( fNumberOfStrings );
        fillNumberOfStringsInXY();


        //fill histos with partons from A
        for ( int iP = 0; iP < fA->nPartons; iP++ )
        {
            double rParton = sqrt(fA->pX[iP]*fA->pX[iP] + fA->pY[iP]*fA->pY[iP] );
            fHistPartonsR->Fill( rParton );
            //1.12.2014: test how many partons become busy after interaction as function of radius
            if(1)
            {
                int partState = fA->pBusy[iP];
                if ( partState == 1  ) fHistBusyPartonsR->Fill(rParton);
                if ( partState == -1  || partState == 2 ) fHistPartonsValenceR->Fill(rParton); //"valence quark"
                if ( partState == 2 ) fHistBusyPartonsValenceR->Fill(rParton); //"valence quark" is BUSY
            }
        }

        //calc nucleon eccentricity
        calcEccentricity();
        fHistNucleonEccentricity->Fill( fEccentricity );
        
        fHist2DNInteractionsVsNstrings->Fill( fNumberOfStrings, fNumberOfStringInteractions );
        fHist2DNInteractionsVsImpactPar->Fill( fImpactParameter, fNumberOfStringInteractions );
        
        //wounded nucleons
        int nWoundedA = fA->getNumberOfWoundedNucleons();
        int nWoundedB = fB->getNumberOfWoundedNucleons();
        fNparticipants = nWoundedA + nWoundedB;
        fHistNumberWoundedNucleons->Fill( fNparticipants );
        //    cout << "nW_A=" << nWoundedA << ", nW_B=" << nWoundedB << endl;
        
        //calc Ncoll
        int nColl = 0;
        for( int i = 0; i < fNumberOfNucleons; ++i )
            for( int j = 0; j < fNumberOfNucleons; ++j )
                if ( matrixNcoll[i][j] )
                    nColl++;
        fNcollisions = nColl;
        fHistNcoll->Fill( fNcollisions );

        //spec draw for saving
        if ( fFlagWriteEventViewCanvases && fEventId < 10 )
        {
            drawEventStructure();
            fCanvEventView->SaveAs( Form( "%s/canv_eventView_%d.eps", fOutDirName.Data(), fEventId));
            fCanvEventView->SaveAs( Form( "%s/canv_eventView_%d.png", fOutDirName.Data(), fEventId));
        }

        // Nov 2017: fictive multiplicities from strings (to control mult distr and for centrality determination)
        for( int i = 0; i < fNumberOfStrings; ++i )
        {
            fMultFromAllStringsFictiveV0 += TMath::Nint( fRand->Poisson( 4.3*fAvMultFromOneStringPerEtaUnitFictive ) ); // 4.3=acceptance of V0, 1.1 particle per unit of rapidity
            fMultFromAllStringsFictiveMidEta += TMath::Nint( fRand->Poisson( 1.0*fAvMultFromOneStringPerEtaUnitFictive ) ); // "|dN/dEta|<0.5"
        }
        fHistFictiveMultDistrV0->Fill( fMultFromAllStringsFictiveV0 );
        fHistFictiveMultDistrMidEta->Fill( fMultFromAllStringsFictiveMidEta );
    } // endl of if fFlagHaveMBcollision
}

void NucleiCollision::createNucleiPair()
{
    bool flagGoodConfiguration = false;
    //nucleus 1 is in (0,0), make pos of the 2nd
    double b2 = 10000; //impact parameter ^2
    double bx = 0;
    double by = 0;
    //    int nTrials = 0;
    while ( !flagGoodConfiguration )
    {
        //        nTrials++;
        if ( fImpactParSpecification == 0) // MB event
        {
            bx = funcImpactPar1D->GetRandom();
            by = 0; // !!!! b always along x!!! //fRand->Uniform( -2*fNucleusRadius*kCoeffRadiusForMB, 2*fNucleusRadius*kCoeffRadiusForMB );
            //            bx = fRand->Uniform( -fNucleusRadius*kCoeffRadiusForMB, fNucleusRadius*kCoeffRadiusForMB );
            //            by = fRand->Uniform( -fNucleusRadius*kCoeffRadiusForMB, fNucleusRadius*kCoeffRadiusForMB );
        }
        else if ( fImpactParSpecification == 1) // precise b
        {
            bx = fImpactParameterByHand; //fImpactParameterByHand_0_100/100. * 2*fNucleusRadius;
            by = 0; //fImpactParameterByHand_0_100/2;
        }
        else if ( fImpactParSpecification == 2) // b in range
        {
            //            bx = fRand->Uniform(fImpactParameter_rangeMin_0_100, fImpactParameter_rangeMax_0_100) /100. * 2*fNucleusRadius;
            bx = funcImpactPar1D->GetRandom( fImpactParameter_rangeMin, fImpactParameter_rangeMax );
            //bx = fRand->Uniform(fImpactParameter_rangeMin, fImpactParameter_rangeMax);
            by = 0; //fImpactParameterByHand_0_100/2;
        }
        
        b2 = bx*bx + by*by;
        if ( sqrt(b2) < 2*fNucleusRadius*kCoeffRadiusForMB )
        {
            fImpactParameter = sqrt(b2);
            //            fHistImpactParameter->Fill( fImpactParameter ); //moved after "have MB event"!
            flagGoodConfiguration = true;
        }
    }
    //    cout << "nTrials=" << nTrials << endl;
    
    //create nuclei
    //    createPartonsOld( fX1, fY1, 0, 0 );
    //    createPartonsOld( fX2, fY2, bx, by );
    createNucleus( 0, fA, 0, 0 );  //fNuclX1, fNuclY1, 0, 0, fX1, fY1 );
    createNucleus( 1, fB, bx, by ); //fNuclX2, fNuclY2, bx, by, fX2, fY2 );
    
    //fill histos with nucleons from A
    for ( int iP = 0; iP < fA->nNucleons; iP++ )
        fHistNucleonsR->Fill( sqrt(fA->nX[iP]*fA->nX[iP] + fA->nY[iP]*fA->nY[iP] ) );

    // 1.12.2014: test impact parameter distibution for nucleon from A wrt nucleons from B
    // can be commented! it doesn't affect the dynamics of the model!
    double xA,yA,xB,yB,dx,dy,nuclR_A,nuclDeltaR;
    if(0)
        for ( int iP = 0; iP < fA->nNucleons; iP++ )
        {
            xA = fA->nX[iP];
            yA = fA->nY[iP];
            nuclR_A = sqrt( xA*xA+yA*yA );
            for ( int iP2 = 0; iP2 < fB->nNucleons; iP2++ )
            {
                xB = fB->nX[iP2];
                yB = fB->nY[iP2];
                dx = xA-xB;
                dy = yA-yB;
                nuclDeltaR = sqrt( dx*dx + dy*dy );
                if ( nuclDeltaR < kNucleonSize )
                {
                    //now fill histograms in dependence of were nucleon from A is situated:
                    if ( nuclR_A < 1 ) //take only "central" nucleons
                        fHistImpParNucleonsFromAB_central->Fill( nuclDeltaR );
                    if ( nuclR_A > 2 && nuclR_A < 4 ) //semicentral
                        fHistImpParNucleonsFromAB_semicentral->Fill( nuclDeltaR );
                    if ( nuclR_A > 6 )
                        fHistImpParNucleonsFromAB_peripheral->Fill( nuclDeltaR );
                }
            }
        }


}

void NucleiCollision::createNucleus( int nucleusId, Nucleus *nucl, float bx, float by )
{
    //remember impact parameter (bx,by)
    nucl->bx = bx;
    nucl->by = by;
    double rRand = 0;
    double phiRand = 0;
    //    double thetaRand = 0;
    if (fNucleusType == nucleus_proton || (fNucleusType == nucleus_p_Pb && nucleusId == 0 ) )
    {
        nucl->nX[0] = bx;
        nucl->nY[0] = by;
        createPartons( nucl, 0 );
    }
    else //all nuclei
    {
        for ( int iP = 0; iP < nucl->nNucleons; iP++ )
        {
            rRand = funcNucleonDensity->GetRandom();
            //        thetaRand = fRand->Uniform( 0, PI );
            phiRand = fRand->Uniform( 0, 2*PI );
            //        double radialPart = 4*PI * rRand * rRand;
            double randCosTheta = fRand->Uniform( -1, 1 );
            double randSinTheta = sqrt(1 - randCosTheta*randCosTheta );
            nucl->nX[iP] = rRand * randSinTheta/*TMath::Sin(thetaRand) */ * TMath::Cos( phiRand ) + bx;
            nucl->nY[iP] = rRand * randSinTheta/*TMath::Sin(thetaRand) */ * TMath::Sin( phiRand ) + by;
            //        nucl->nZ[iP] = rRand * TMath::Cos(thetaRand);
            //        nucl->nTheta[iP] = thetaRand;
            //        fHistNucleonDensityR->Fill( nucl->nX[iP] ); //rRand );

            createPartons( nucl, iP ); //nuclX[iP], nuclY[iP], x, y );
        }
    }
}

void NucleiCollision::createPartons( Nucleus *nucl, int nId ) //float nucleonX, float nucleonY, float *x, float *y )
{
    //nId - id of a nucleon inside nucleus
    
    //generate n of partons in this nucleon
    int curN = nucl->nPartons; //overall number of partons in this NUCLEUS
    // MAY 9, 2018: add possibility to fix nPartons in eachNucleon - like in WQM (?) )
    int nPartonsInThisNucleon = ( fFlagUsePoissonianNpartonsFluctuations ? TMath::Nint( fRand->Poisson(fMeanNofPartonsInNucleon) )
                                                                         : fMeanNofPartonsInNucleon ); // ADDED ON MAY 9, 2018!
    if ( curN + nPartonsInThisNucleon > nucl->maxNofPartons )
    {
        cout << "MAX N OF PARTONS EXCEEDED!!! => don't generate partons..." << endl;
        return;
    }
    
    //distribute partons
    double rRand = 0;
    double phiRand = 0;
    //    double thetaRand = 0;
    float a,b;
    for ( int iP = 0; iP < nPartonsInThisNucleon; iP++ )
    {
        int id = curN+iP;
        nucl->pNid[id] = nId;

        // (TMP?) 1.12.2014 - for some spec studies:
        // make first three partons to be "valence quarks" (pBusy=-1):
        nucl->pBusy[id] = (iP < 3 ? -1 : 0);

        if(1)//variant with Gaus2D position of partons
        {
            fRand->Rannor( a, b ); //mean=0 and sigma=1
            const float kNucleonGausSigma = fNucleonGaussianRadius; //0.4 - WAS A DEFAULT BEFORE 27 NOV 2017; //0.6
            nucl->pX[id] = /*nucl->nX[nId] +*/ a*kNucleonGausSigma;//kNucleonSize; //sigma 0.6?..
            nucl->pY[id] = /*nucl->nY[nId] +*/ b*kNucleonGausSigma;//kNucleonSize;
        }
        else if(0)//variant with density func
        {
            rRand = funcPartonDensity->GetRandom();
            //            thetaRand = fRand->Uniform( 0, PI );
            phiRand = fRand->Uniform( 0, 2*PI );
            double randCosTheta = fRand->Uniform( -1, 1 );
            double randSinTheta = sqrt(1 - randCosTheta*randCosTheta );
            
            //            fHistNucleonDensityR->Fill( rRand );
            //            nucl->pX[id] = nucl->nX[nId] + rRand * sin(thetaRand) * cos( phiRand );
            //            nucl->pY[id] = nucl->nY[nId] + rRand * sin(thetaRand) * sin( phiRand );
            nucl->pX[id] = /*nucl->nX[nId] +*/ rRand * randSinTheta/*TMath::Sin(thetaRand) */ * cos( phiRand );
            nucl->pY[id] = /*nucl->nY[nId] +*/ rRand * randSinTheta/*TMath::Sin(thetaRand) */ * sin( phiRand );
        }
        else //variant with EXPONENTIAL density
        {
            rRand = fRand->Exp(kNucleonSize);
            //            thetaRand = fRand->Uniform( 0, PI );
            phiRand = fRand->Uniform( 0, 2*PI );
            // 3D proton
            //            double randCosTheta = fRand->Uniform( -1, 1 );
            //            double randSinTheta = sqrt(1 - randCosTheta*randCosTheta );
            //            nucl->pX[id] = nucl->nX[nId] + rRand * randSinTheta/*TMath::Sin(thetaRand) */ * cos( phiRand );
            //            nucl->pY[id] = nucl->nY[nId] + rRand * randSinTheta/*TMath::Sin(thetaRand) */ * sin( phiRand );
            // 2D proton
            nucl->pX[id] = /*nucl->nX[nId] +*/ rRand * cos( phiRand );
            nucl->pY[id] = /*nucl->nY[nId] +*/ rRand * sin( phiRand );
        }
    } // end of parton loop for nucleon

    // RECENTER PARTONS IN NUCLEON:
    // move partons to have center of nucleon = a center of mass for 3 valence quarks (Sept 2017)
    double xCM = 0, yCM = 0;
    const int nPartonsForCM = (nPartonsInThisNucleon > 3 ? 3 : nPartonsInThisNucleon); // usually: take 3 valence quark for that
    for (int vQ=0; vQ<nPartonsForCM; vQ++)
    {
        int id = curN+vQ;
        xCM += nucl->pX[id];
        yCM += nucl->pY[id];
    }
    xCM /= nPartonsForCM;
    yCM /= nPartonsForCM;
    //    cout << "xCM = " << xCM << ", yCM = " << yCM << endl;
    for ( int iP = 0; iP < nPartonsInThisNucleon; iP++ )
    {
        int id = curN+iP;
        //        cout << "id = " << id << endl;
        nucl->pX[id] -= xCM;
        nucl->pY[id] -= yCM;
    }

    // finally: place partons according to nucleon position
    for ( int iP = 0; iP < nPartonsInThisNucleon; iP++ )
    {
        int id = curN+iP;

        nucl->pX[id] += nucl->nX[nId];
        nucl->pY[id] += nucl->nY[nId];
    }

    nucl->nPartons += nPartonsInThisNucleon;
}

void NucleiCollision::createPartonsOld(float *x, float *y, float bx, float by )
{
    double rRand = 0;
    double phiRand = 0;
    
    //distribute partons in hard sphere
    double m = 0;
    double mnoj = 4./3.* PI;
    double R3 = fNucleusRadius*fNucleusRadius*fNucleusRadius;
    for ( int iP = 0; iP < fMaxPartons; iP++ )
    {
        //r.Sphere(x,y);
        m = fRand->Uniform(0,mnoj*R3);
        //        if (iP==0) cout << m << endl;
        rRand = sqrt( fNucleusRadius*fNucleusRadius - pow( R3 - m / mnoj, 2./3. ) );
        phiRand = fRand->Uniform( 0, 2*PI );
        x[iP] = rRand * cos( phiRand ) + bx;
        y[iP] = rRand * sin( phiRand ) + by;
    }
}

void NucleiCollision::createStrings()
{

    fNumberOfStrings = 0;
    fNumberOfWoundedPartonsA = 0;
    fNumberOfWoundedPartonsB = 0;
    fNumberOfWoundedOnlyOncePartonsA = 0;
    fNumberOfWoundedOnlyOncePartonsB = 0;

    //    fNumberOfClusters = 0;
    fHist2DNstringsXY_thisEvent->Reset();
    fHist2DNstringsXY_thisEv_CoreCorona->Reset();


//    int nWoundedPartons = nWoundedPartonsA + nWoundedPartonsB;
//    int nSources = ( !fFlagConsiderWoundedPartonsAsStrings ? distArraySize /*default*/
//                                                       : nWoundedPartons ); // May 2018: TRY WQM

    if ( !fFlagConsiderWoundedPartonsAsStrings ) // default! "burning-out" partons or "wrong" mode of WQM used for QM2018 poster
    {
        //find the closest parton pairs
        fPartonInteractionsFinder->FindMinDistancesBetweenPairs(fA->pX,fA->pY,fB->pX,fB->pY,fA->nPartons,fB->nPartons
                                                                , fFlagOnlyOneInteractionPerParton ); //fX1,fY1,fX2,fY2,fMaxPartons,fMaxPartons);

        //define "strings" as the middles between interacting pairs
        int distArraySize = 0;
        DistanceEntry* arrDist = fPartonInteractionsFinder->GetMinDistanceArray( distArraySize );

        for ( int iP = 0; iP < distArraySize; iP++ )
        {
            if ( !arrDist[iP].inInteraction )
                continue;
            int id1 = arrDist[iP].x;
            int id2 = arrDist[iP].y;

            //mark nucleons as wounded
            int idNucleonInA = fA->pNid[id1];
            int idNucleonInB = fB->pNid[id2];
            fA->nuclWounded[idNucleonInA]++;
            fB->nuclWounded[idNucleonInB]++;

            matrixNcoll[idNucleonInA][idNucleonInB] = 1;

            //mark partons to be in interaction with each other // (TMP?) 1.12.2014 - some spec studies
            if (0)//1) // commented in Sept 2017 - to use stringOrigin (below)
            {
                //            double rParton = sqrt(fA->pX[id1]*fA->pX[id1] + fA->pY[id2]*fA->pY[id2] );
                //            bool hardScatteringFlag = false; // 2.12.14: may be useful for the future: ( rParton < fPartonInteractionDistance/2 ? true: false );
                //            changePartonState( fA->pBusy[id1] );//, hardScatteringFlag );
                //            changePartonState( fB->pBusy[id2] );//, hardScatteringFlag );
                changePartonPairState( fA->pBusy[id1], fB->pBusy[id2] );
            }

            fXstring[fNumberOfStrings] = ( fA->pX[id1] + fB->pX[id2] ) / 2.;
            fYstring[fNumberOfStrings] = ( fA->pY[id1] + fB->pY[id2] ) / 2.;
            //cout << "id1=" << id1 << ", id2=" << id2 << endl;

            fDistanceBetweenPartonsForString[fNumberOfStrings] = sqrt( arrDist[iP].dist );

            // string origin (Sept 2017)
            short stringOrigin = 0; // sea quarks/gluons
            if(0)cout << "status of partons for string id " << fNumberOfStrings << ": "
                      << fA->pBusy[id1] << " " << fB->pBusy[id2] << endl;
            if ( fA->pBusy[id1] == -1 && fB->pBusy[id2] == -1 ) // both partons are valence quarks
                stringOrigin = 2;
            else if ( fA->pBusy[id1] == -1 || fB->pBusy[id2] == -1 ) // one valence, one sea/gluon
                stringOrigin = 1;
            fStringOrigin[fNumberOfStrings] = stringOrigin;

            fHistMinDistBetweenInteractingPartons->Fill( sqrt( arrDist[iP].dist ) );
            if ( stringOrigin == 2)  // interaction b/n valence quarks
                fHistMinDistBetweenInteractingPartonsValence->Fill( sqrt( arrDist[iP].dist ) );
            else // all other cases
                fHistMinDistBetweenInteractingPartonsSoft->Fill( sqrt( arrDist[iP].dist ) );

            fNumberOfStrings++;
        }
    }
    else // May 2018: TRY WQM
    {
        //
        for ( int i = 0; i < fA->nPartons; i++ )
        {
            for ( int j = 0; j < fB->nPartons; j++ )
            {
                float dx = fA->pX[i] - fB->pX[j];
                float dy = fA->pY[i] - fB->pY[j];
                if ( fabs(dx) > fPartonInteractionDistance || fabs(dy) > fPartonInteractionDistance )
                    continue;
                if( dx*dx + dy*dy < fPartonInteractionDistance*fPartonInteractionDistance )
                {
                    fA->isWoundedParton[i]++;
                    fB->isWoundedParton[j]++;
                }
            }

        }

        // count wounded partons in nucleus A, mark nucleons as wounded, add strings:
        for ( int i = 0; i < fA->nPartons; i++ )
        {
            if ( fA->isWoundedParton[i] == 0 )
                continue;

            fA->nuclWounded[ fA->pNid[i] ]++;
            fXstring[fNumberOfStrings] = fA->pX[i];
            fYstring[fNumberOfStrings] = fA->pY[i];

            fNumberOfWoundedPartonsA++;
            fNumberOfStrings++;

            if ( fA->isWoundedParton[i] == 1 ) // count if wounded only once
                fNumberOfWoundedOnlyOncePartonsA++;
        }

        // count wounded partons in nucleus B, mark nucleons as wounded, add strings:
        for ( int i = 0; i < fB->nPartons; i++ )
        {
            if ( fB->isWoundedParton[i] == 0 )
                continue;

            fB->nuclWounded[ fB->pNid[i] ]++;
            fXstring[fNumberOfStrings] = fB->pX[i];
            fYstring[fNumberOfStrings] = fB->pY[i];

            fNumberOfWoundedPartonsB++;
            fNumberOfStrings++;

            if ( fB->isWoundedParton[i] == 1 ) // count if wounded only once
                fNumberOfWoundedOnlyOncePartonsB++;
        }

    } // end of if-else - burning-off or WQM

    // some actions on strings:
    for ( int iP = 0; iP < fNumberOfStrings; iP++ )
    {
        //        fStringBoostMagn[iP] = 0;
        //        fStringBoostAngle[iP] = 0;
        //        fFlagStringInInteraction[iP] = 0;//false;
        //        fNStringsInCluster[iP] = 0;//false;
        //        fClusterIdForString[iP] = -1;
        //        fFlagStringIsHardInteraction[iP] = 0;//false;

        //        fFlagStringIsHardInteraction[iP] = ( fRand->Uniform() < fHardScatteringProbability ? 1 : 0 );

        // calc angle of the radius-vector of the string
        float middleX_shifted = fXstring[iP] - fImpactParameter/2;
        double stringRadiusVectorAngle = atan( fYstring[iP]/middleX_shifted ); //   sqrt(middleX*middleX+middleY*middleY) );
        if ( middleX_shifted < 0 )
            stringRadiusVectorAngle = TMath::Pi()-stringRadiusVectorAngle;
        FixAngleInTwoPi(stringRadiusVectorAngle);
        FixAngleInTwoPi(stringRadiusVectorAngle);
        fStringRadiusVectorAngle[iP] = stringRadiusVectorAngle;
        fHistStringRadiusVectorPhi->Fill( stringRadiusVectorAngle );


        fHistStringPositionRadius->Fill( sqrt(fXstring[iP]*fXstring[iP]+fYstring[iP]*fYstring[iP])/fNucleusRadius );
        fHist2DNstringsXY_thisEvent->Fill( middleX_shifted, fYstring[iP] );
        // May 2018: Core-corona:
        fHist2DNstringsXY_thisEv_CoreCorona->Fill( middleX_shifted, fYstring[iP] );
    }
    
    //set flag MB collision if have at least one string
    if ( fNumberOfStrings > 0 )
        fFlagHaveMBcollision = true;
}



void NucleiCollision::calcEccentricity()
{
    float num = 0;
    float denom = 0;
    int nWounded = 0;

    // accumulate values from nucleus A and B
    getIngredientsForEccentricity( fA, fB, nWounded, num, denom );
    getIngredientsForEccentricity( fB, fA, nWounded, num, denom );

    //calc eccentricity
    if ( nWounded>0 )
    {
        num /= nWounded;
        denom /= nWounded;
        fEccentricity = num/denom;
        //        cout << "eccentricity: " << num << " " << denom << " " << fEccentricity << endl;
    }
    else
        fEccentricity = -1;
}

void NucleiCollision::fillNumberOfStringsInXY()
{
    for ( int i = 0; i < fHist2DNstringsXY_thisEvent->GetNbinsX(); i++) //x
    {
        for ( int j = 0; j < fHist2DNstringsXY_thisEvent->GetNbinsY(); j++) //y
        {
            double binValue = fHist2DNstringsXY_thisEvent->GetBinContent(i+1,j+1);
            //            cout << binValue << endl;

            double nStr = fHist2DNstringsXY->GetBinContent( i+1, j+1 );
            fHist2DNstringsXY->SetBinContent( i+1, j+1, nStr + binValue );

            double nStr2 = fHist2DN2stringsXY->GetBinContent( i+1, j+1 );
            fHist2DN2stringsXY->SetBinContent( i+1, j+1, nStr2 + binValue*binValue );

            //            cout << fHist2DNstringsXY->GetBinContent(i+1,j+1) << " ";
        }
    }
    //    cout << endl;


}

void NucleiCollision::getIngredientsForEccentricity( Nucleus *nucl, Nucleus *nucl2, int &nWounded, float &num, float &denom )
{
    // nucl2 is ONLY for calculation of the middle b/n nuclei A&B centers
    float midX = (nucl->bx + nucl2->bx) / 2;
    float midY = (nucl->by + nucl2->by) / 2;

    double x, y, x2, y2;
    for ( int i = 0; i < nucl->nNucleons; i++ )
    {
        if ( nucl->nuclWounded[i] > 0 )
        {
            nWounded++;
            x = nucl->nX[i] - midX;
            y = nucl->nY[i] - midY;
            x2 = x*x;
            y2 = y*y;
            num += y2 - x2;
            denom += y2 + x2;
        }
    }
}


//######### drawings
void NucleiCollision::drawStatisticHists()
{
    //save analyzer lists to file
    TFile *fileNuclStructStats = new TFile( Form( "%s/stats_NuclearStructure.root", fOutDirName.Data() ),"RECREATE");

    if (!fCanvEventStatistics)
    {
        fCanvEventStatistics = new TCanvas("Stats Canvas","Nuclei Generation Statistics",180,50,1200,800);
        fCanvEventStatistics->Divide(3,3);
    }
    if (1)
    {
        int padId = 1;

        //pad
        fCanvEventStatistics->cd(padId++);
        //        fHistForcesInClustersAbsX->SetLineColor(kGreen);
        //        fHistForcesInClustersAbsY->SetLineColor(kRed);
        //        fHistForcesInClustersAbsX->DrawCopy();
        //        fHistForcesInClustersAbsY->DrawCopy("same");

        //        fHistForcesInClustersAbsX->Write();
        //        fHistForcesInClustersAbsY->Write();

        fHistImpParNucleonsFromAB_semicentral->SetLineColor(kRed);
        fHistImpParNucleonsFromAB_peripheral->SetLineColor(kBlack);
        fHistImpParNucleonsFromAB_central->DrawNormalized();
        fHistImpParNucleonsFromAB_semicentral->DrawNormalized("same");
        fHistImpParNucleonsFromAB_peripheral->DrawNormalized("same");
        //        fHistImpParNucleonsFromAB->Write();

        //pad
        fCanvEventStatistics->cd(padId++);

        fHistMinDistBetweenInteractingPartons->DrawNormalized();
        fHistMinDistBetweenInteractingPartonsSoft->SetLineColor(kGreen);
        fHistMinDistBetweenInteractingPartonsSoft->DrawNormalized("same");
        fHistMinDistBetweenInteractingPartonsValence->SetLineColor(kRed);
        fHistMinDistBetweenInteractingPartonsValence->DrawNormalized("same");

        //pad 1
        fCanvEventStatistics->cd(padId++);
        fHistNucleonsR->DrawNormalized();
        fHistPartonsR->SetLineColor(kRed);
        fHistPartonsR->DrawNormalized("same");
        //        fHistBusyPartonsR->SetLineColor(kGreen);
        //        fHistBusyPartonsR->DrawNormalized("same");

        fHistNucleonsR->Write();
        fHistPartonsR->Write();
        //pad 2
        //        fCanvEventStatistics->cd(padId++);
        //        fHistNstrings->DrawCopy();
        //pad 3
        //        fCanvEventStatistics->cd(padId++);
        //        fHistStringInteractions->DrawCopy();
        //pad 4
        fCanvEventStatistics->cd(padId++);
        fHistStringPositionRadius->DrawCopy();
        fHistStringPositionRadius->Write();
        //pad 5
        fCanvEventStatistics->cd(padId++);
        //        fHistStringInteractionsPhi->DrawCopy();
        //        fHistStringInteractionsPhi->Write();

        fHistStringRadiusVectorPhi->SetLineColor(kRed);
        fHistStringRadiusVectorPhi->DrawCopy( "same" );
        //pad 6
        //        fCanvEventStatistics->cd(6);
        //        fHistStringInteractionsDistance->DrawCopy();
        //pad 7
        fCanvEventStatistics->cd(padId++);//->SetLogy();
        //        fHistStringInteractionsMagnitude->DrawCopy();
        //        fHistStringInteractionsMagnitude->Write();
        //        cout << ">> mean string boost (<beta>) = " << fHistStringInteractionsMagnitude->GetMean() << endl;
        //pad 8
        fCanvEventStatistics->cd(padId++);
        fHistImpactParameter->DrawCopy();
        fHistImpactParameter->Write();
        //        //pad 9
        //        fCanvEventStatistics->cd(padId++);
        //        fHist2DNInteractionsVsNstrings->DrawCopy( "colz" );
        //        //pad 10
        //        fCanvEventStatistics->cd(10);
        //        fHistStringsInClusters->DrawCopy();
        //        cout << "mean cluster size: " << fHistStringsInClusters->GetMean() << endl;
        //        //pad 11
        //        fCanvEventStatistics->cd(11);
        //        fHist2DStringsInClustersVsB->ProfileX()->DrawCopy("colz");
        //        //pad 12
        //        fCanvEventStatistics->cd(12);
        //        (fHist2DNInteractionsVsNstrings->ProfileX())->DrawCopy();
        //        //pad 13
        //        fCanvEventStatistics->cd(13);
        //        (fHist2DNInteractionsVsImpactPar->ProfileX())->DrawCopy();
        //        //pad 14
        //        fCanvEventStatistics->cd(14);
        //        fHist2DStringsSVsImpactS->DrawCopy("colz");
        //        //    (fHist2DStringsSVsImpactS->ProfileX())->DrawCopy();
        //        fCanvEventStatistics->cd(15);
        //        (fHist2DStringsSVsImpactS->ProfileY())->DrawCopy();
        //        //pad 15
        //        fCanvEventStatistics->cd(16);
        //        fHistNStringsInLargestCluster->DrawCopy();

        //pad
        fCanvEventStatistics->cd(padId++);
        fHistNstrings->SetLineColor(kRed);
        fHistNstrings->DrawNormalized();
        fHistNcoll->SetLineColor(kGreen);
        fHistNcoll->DrawNormalized("same");
        fHistNumberWoundedNucleons->SetLineColor(kBlue);
        fHistNumberWoundedNucleons->DrawNormalized("same");

        fHistNstrings->Write();
        fHistNcoll->Write();
        fHistNumberWoundedNucleons->Write();

        fHistFictiveMultDistrV0->Write();
        fHistFictiveMultDistrMidEta->Write();

        //pad
        fCanvEventStatistics->cd(padId++);
        //        fHistNucleonEccentricity->DrawCopy();
        //        fHistNucleonEccentricity->Write();
        fHistBusyPartonsR->Divide(fHistPartonsR);
        fHistBusyPartonsR->DrawCopy();
        fHistBusyPartonsR->Write();

        fHistBusyPartonsValenceR->Divide(fHistPartonsValenceR);
        fHistBusyPartonsValenceR->SetLineColor(kRed);
        fHistBusyPartonsValenceR->DrawCopy("same");
        fHistBusyPartonsValenceR->Write();

        //ratio of ratio:
        fHistBusyPartonsValenceR->Divide(fHistBusyPartonsR);
        fHistBusyPartonsValenceR->SetLineColor(kBlack);
        fHistBusyPartonsValenceR->DrawCopy("same");
        fHistBusyPartonsValenceR->Write();

        fileNuclStructStats->Close();

        fCanvEventStatistics->SaveAs( Form( "%s/canvNuclearStructureStats.root", fOutDirName.Data()) );
        fCanvEventStatistics->SaveAs( Form( "%s/canvNuclearStructureStats.eps", fOutDirName.Data()));
        fCanvEventStatistics->SaveAs( Form( "%s/canvNuclearStructureStats.png", fOutDirName.Data()));
    }
    
    if (0)
    {
        //nucleon density QA hist
        //        TCanvas *lCanvNuclDensity = new TCanvas("lCanvNuclDensity","Nucleon Density",250,80,600,600);
        //        lCanvNuclDensity->cd(); //to avoid warning
        //        fHistNucleonDensityR->DrawCopy();
        
        //wounded nucleons QA hist
        TCanvas *lCanvNwoundedNucl = new TCanvas("lCanvNwoundedNucl","N part nucleons, N coll",280,80,600,600);
        lCanvNwoundedNucl->cd(); //to avoid warning
        fHistNcoll->DrawCopy();
        fHistNumberWoundedNucleons->SetLineColor(kRed);
        fHistNumberWoundedNucleons->DrawCopy("same");
    }
    
}

void NucleiCollision::finalActions()
{    
    cout << "########## Begin NucleiCollision::finalActions()" << endl;
    if (1)
    {
        //percolation parameter vs bImpact
        double S_impact = circlesIntersectionArea(fImpactParameter,fNucleusRadius,fNucleusRadius);
        const double S_string = TMath::Pi() * 0.25*0.25;//fPartonInteractionDistance*fPartonInteractionDistance;
        //        //    fHist2DStringsSVsImpactS->Fill( S_impact/*fImpactParameter/fNucleusRadius*/, fNumberOfStrings * S_string / S_impact );
        //        fHist2DStringsSVsImpactS->Fill( fImpactParameter, fNumberOfStrings * S_string / S_impact );
        double percPar = fHistNstrings->GetMean() * S_string / S_impact;
        double percParErr = fHistNstrings->GetRMS() * S_string / S_impact;
        cout << "percolation par for this b: " << fImpactParameter << " " <<  fNumberOfStrings * S_string / S_impact << endl;

        ofstream fout( Form( "%s/tmpTextOutput0.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout << fImpactParameterByHand
             << " " << percPar
             << " " << percParErr
             << endl;
        fout.close();

    }

    if (1)
    {
        //string statistics
        float nStringsMean = fHistNstrings->GetMean();
        float ratioToNstrings = fHistNStringsInLargestCluster->GetMean() / nStringsMean;
        cout << "b vs mean number of strings in largest cluster: "
             << fImpactParameterByHand  << " "
                //         << fHist2DStringsSVsImpactS->ProfileY()->GetMean() << " "
             << ratioToNstrings << endl;
        cout << "mean N strings = " << nStringsMean << endl;

        fMeanNstringsFinal = nStringsMean; // to take out into spec hist in ManagerNucleiCollisions

        double err1 = fHistNStringsInLargestCluster->GetMeanError();
        double err2 = fHistNstrings->GetMeanError();
        double ratioToNstringsTotal = sqrt(err1*err1 + err2*err2) / nStringsMean;

        //print fraction of strings in cluster to n of all strings
        ofstream fout1( Form( "%s/tmpTextOutput1.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout1 << fImpactParameterByHand //_0_100
              << " " << ratioToNstrings
              << " " << ratioToNstringsTotal
              << endl;
        fout1.close();

        //print mean number of strings in the EVENT vs impact par
        float nStringsRMS= fHistNstrings->GetRMS();
        ofstream fout1_0( Form( "%s/tmpTextOutput1_0.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout1_0 << fImpactParameterByHand //_0_100
                << " " << nStringsMean
                << " " << nStringsRMS
                << endl;
        fout1_0.close();
    }
    
    //print nWounded
    cout << "mean N wounded = " << fHistNumberWoundedNucleons->GetMean() << endl;
    if (1)
    {
        ofstream fout2( Form( "%s/tmpTextOutput2.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout2 << fImpactParameterByHand //_0_100
              << " " << fHistNumberWoundedNucleons->GetMean()
              << " " << fHistNumberWoundedNucleons->GetMeanError()
              << endl;
        fout2.close();
    }
    
    //print Ncoll
    if (1)
    {
        cout << "mean Ncoll = " << fHistNcoll->GetMean() << endl;
        ofstream fout3( Form( "%s/tmpTextOutput3.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout3 << fImpactParameterByHand //_0_100
              << " " << fHistNcoll->GetMean()
              << " " << fHistNcoll->GetMeanError()
              << endl;
        fout3.close();
    }
    
    // print Fx, Fy results
    //    float meanFx = fHistForcesInClustersAbsX->GetMean();
    //    float meanFy = fHistForcesInClustersAbsY->GetMean();
    //    double errFx = fHistForcesInClustersAbsX->GetMeanError();
    //    double errFy = fHistForcesInClustersAbsY->GetMeanError();
    //    double errNumDenom = errFx + errFy;
    //    double errTotal = sqrt(errNumDenom*errNumDenom + errNumDenom*errNumDenom);

    //    float forceEccentricity = (meanFx-meanFy) / (meanFx+meanFy);
    //    if (1)
    //    {
    //        cout << "forceEccentricity = " << forceEccentricity << endl;
    //        ofstream fout( Form( "%s/tmpTextOutput4.txt", fOutDirName.Data() ), ios::out | ios::binary);
    //        fout << fImpactParameterByHand
    //             << " " << forceEccentricity
    //             << " " << errTotal
    //             << endl;
    //        fout.close();
    //    }

    if (1)
    {
        //        ofstream fout( Form( "%s/tmpTextOutput5.txt", fOutDirName.Data() ), ios::out | ios::binary);
        //        fout << fImpactParameterByHand
        //             << " " << fHistStringInteractionsMagnitude->GetMean()
        //             << " " << fHistStringInteractionsMagnitude->GetMeanError()
        //             << endl;
        //        fout.close();
    }

    // transverse density fluctuations
    if (1)
    {
        for ( int i = 0; i < fHist2DNstringsXY->GetNbinsX(); i++) //x
        {
            for ( int j = 0; j < fHist2DNstringsXY->GetNbinsY(); j++) //y
            {
                double binValue     = fHist2DNstringsXY->GetBinContent(i+1,j+1);
                //                cout << binValue << endl;
                double binValue2    = fHist2DN2stringsXY->GetBinContent(i+1,j+1);
                binValue /= fEvTrialsSuccess;
                binValue2 /= fEvTrialsSuccess;
                fHist2DNstringsXY->SetBinContent( i+1,j+1, binValue);
                fHist2DN2stringsXY->SetBinContent( i+1,j+1, binValue2);

                double sigma =sqrt( binValue2 - binValue*binValue );
                double omega = ( binValue > 0 ? sigma/binValue : 0 );
                //                fHist2DSigmaNstringsXY->SetBinContent( i+1,j+1, sigma );
                fHist2DSigmaNstringsXY->SetBinContent( i+1,j+1, omega );
            }
        }
        TCanvas *lCanvTransvStringDensity = new TCanvas("lCanvTransvStringDensity","string density in xy",280,10,800,800);
        lCanvTransvStringDensity->Divide(2,2);
        lCanvTransvStringDensity->cd(1); //to avoid warning
        fHist2DSigmaNstringsXY->DrawCopy( "cont1z" );
        lCanvTransvStringDensity->cd(2);
        fHist2DNstringsXY->DrawCopy( "colz" );

        // ### May 2018: add profiles across centers (along x-axis)
        // ### projection for omega
        int nBinsTmp = fHist2DSigmaNstringsXY->GetNbinsY();
        TH1D *profXomega = fHist2DSigmaNstringsXY->ProjectionX( "profNstringsFluctOmega", (nBinsTmp-1)/2-3, (nBinsTmp-1)/2+3 );
        lCanvTransvStringDensity->cd(3);
        profXomega->DrawCopy();

        // ### projection for omega
        nBinsTmp = fHist2DNstringsXY->GetNbinsY();
        TH1D *profXmeanNstrings = fHist2DNstringsXY->ProjectionX( "profNstrings", (nBinsTmp-1)/2-3, (nBinsTmp-1)/2+3 );
        lCanvTransvStringDensity->cd(4);
        profXmeanNstrings->DrawCopy();

    }


    //print mean string boost
    //    cout << "(b, mean string boost, err) = " << fImpactParameterByHand << " " << fHistStringInteractionsMagnitude->GetMean()
    //    << " " << fHistStringInteractionsMagnitude->GetMeanError() << endl;

    //print cross section
    cout << " >>>>EvTrialsSuccess=" << fEvTrialsSuccess << endl;
    cout << " >>>>EvTrialsFailed=" << fEvTrialsFailed << endl;
    int nTrials = fEvTrialsSuccess + fEvTrialsFailed;
    if ( nTrials > 0 )
    {
        double radiusOfTrials = 2*fNucleusRadius*kCoeffRadiusForMB;
        double areaOfTrials = TMath::Pi() * radiusOfTrials * radiusOfTrials;
        cout << " >>>>> areaOfTrials = " << areaOfTrials << endl;
        double crossSection = (double)fEvTrialsSuccess / nTrials * areaOfTrials / 100; //100 fm2 in 1 barn
        cout << " >>>>> crossSection (in barns?..)= " << crossSection << endl;

        ofstream fout_crossSection( Form( "%s/tmp_cross_section.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout_crossSection << crossSection << endl;
        fout_crossSection.close();

        fCrossSectionFinal = crossSection;
    }

    cout << "########## End NucleiCollision::finalActions()" << endl;
}

double NucleiCollision::getNu()
{
    return (2.*fNcollisions ) / fNparticipants;
}


void NucleiCollision::drawEventStructure()
{
    if (!fCanvEventView)
    {
        fCanvEventView = new TCanvas("Event Canvas","Event View",600,70,700,700);
        //fCanvEventView->Divide(2,2);
    }
    fCanvEventView->Clear();

    drawPartons();
    drawStrings();
    if ( fFlagHaveMBcollision )
    {
        //        drawStringInteractions();
        //        drawForcesInsideClusters();
        //            drawStringBoosts();
    }
    //    fCanvEventView->SaveAs( Form( "eventViews/eventView_%d.eps", fEventId) );
}

void NucleiCollision::drawPartons()
{
    fCanvEventView->cd();
    //TF1 *f1 = new TF1( "ball", "4./3. * 3.1415926 * ( 1 - ( 1 - x*x )^(3./2.) )", 0, 1 );
    //f1->Draw();
    
    //    double fVisNucleusRadiusNucleus = 0.4; //visual size nucleus
    TEllipse *el = new TEllipse( 0.5, 0.5, fVisNucleusRadiusNucleus, fVisNucleusRadiusNucleus );
    el->Draw();
    
    //    TH2D *hist2DNuclear = new TH2D("histNuclear2D","histNuclear2D"
    //                                   ,10,-fNucleusRadius,fNucleusRadius,10,-fNucleusRadius,fNucleusRadius );
    
    //visualize "nucleons"
    double rVisNucleon = fVisNucleusRadiusNucleus * kNucleonSize/fNucleusRadius;
    TEllipse **elNucl;
    elNucl = new TEllipse* [fNumberOfNucleons];
    TEllipse **elNucl2;
    elNucl2 = new TEllipse* [fNumberOfNucleons];
    
    for ( int iP = 0; iP < fNumberOfNucleons; iP++ )
    {
        elNucl[iP] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * fA->nX[iP] / fNucleusRadius, 0.5 + fVisNucleusRadiusNucleus * fA->nY[iP] / fNucleusRadius
                                   , rVisNucleon, rVisNucleon );
        elNucl[iP]->SetLineColor( kGray );
        elNucl[iP]->Draw();
        //        hist2DNuclear->Fill( fX1[iP], fY1[iP] );
        
        elNucl2[iP] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * fB->nX[iP] / fNucleusRadius, 0.5 + fVisNucleusRadiusNucleus * fB->nY[iP] / fNucleusRadius
                                    , rVisNucleon, rVisNucleon );
        elNucl2[iP]->SetLineColor( kRed-7 );
        elNucl2[iP]->Draw();
    }
    
    //visualize "partons"
    double rVisParton = fVisNucleusRadiusNucleus * fPartonInteractionDistance/2/fNucleusRadius; //0.01; //fRadiusParton / fNucleusRadius * fVisNucleusRadiusNucleus; //visual nucleon
    TEllipse **elPoints;
    elPoints = new TEllipse* [fA->nPartons];
    for ( int iP = 0; iP < fA->nPartons; iP++ )
    {
        elPoints[iP] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * fA->pX[iP] / fNucleusRadius, 0.5 + fVisNucleusRadiusNucleus * fA->pY[iP] / fNucleusRadius
                                     , rVisParton, rVisParton );
        if ( fA->pBusy[iP] == -1 ) // valence quark
            elPoints[iP]->SetFillColorAlpha( kBlue-10, 0.6 );
        elPoints[iP]->Draw();
        
    }
    
    TEllipse **elPoints2;
    elPoints2 = new TEllipse* [fB->nPartons];
    for ( int iP = 0; iP < fB->nPartons; iP++ )
    {
        
        elPoints2[iP] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * fB->pX[iP] / fNucleusRadius, 0.5 + fVisNucleusRadiusNucleus * fB->pY[iP] / fNucleusRadius
                                      , rVisParton, rVisParton );
        elPoints2[iP]->SetLineColor( kRed );
        if ( fB->pBusy[iP] == -1 ) // valence quark
            elPoints2[iP]->SetFillColorAlpha( kRed-10, 0.6 );
        elPoints2[iP]->Draw();
    }
    
    
    //    fCanvEventView->cd(3);
    //    hist2DNuclear->DrawCopy();
    //    fCanvEventView->cd(2);
    delete [] elNucl;
    delete [] elNucl2;
    delete [] elPoints;
    delete [] elPoints2;
}

void NucleiCollision::drawStrings()
{
    //visualize "strings"
    TEllipse **elPointsStrings = new TEllipse* [fNumberOfStrings];
    //    float rVisString = 0.008;//rVisParton * 0.2;
    float rVisString = 0.006;//rVisParton * 0.2;
    //    float rVisString = 0.002;//rVisParton * 0.2;
    //    int nUsedClusterColors = 0;
    for ( int iP = 0; iP < fNumberOfStrings; iP++ )
    {
        elPointsStrings[iP] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * fXstring[iP] / fNucleusRadius, 0.5 + fVisNucleusRadiusNucleus * fYstring[iP] / fNucleusRadius
                                            , rVisString, rVisString );
        //        elPointsStrings[iP]->SetLineColor( fNStringsInCluster[iP] ? kGreen+fNStringsInCluster[iP] : kBlue );
        //        elPointsStrings[iP]->SetFillColor( fNStringsInCluster[iP] ? kGreen+fNStringsInCluster[iP] : kBlue );
        //        int stringIdInCluster = fNStringsInCluster[iP];
        //        int clusterId = fClusterIdForString[iP];
        
        //        cout << "draw string: " << stringIdInCluster << endl;
        //        int clusterColor = ( stringIdInCluster < 4 ? kGreen+stringIdInCluster : kRed+1 );
        
        //        int clusterColor = kOrange;
        //        if ( fNStringsInCluster[iP] > kNstringsWhenClustersIsConsideredLarge )
        //        {
        //            clusterColor = kBlue;//( clusterId < nClusterColors ? kClusterColors[clusterId] : kAzure+10 );
        //        }

        int stringColor = kBlue;
        //mark hard interactions
        //        if ( fFlagStringIsHardInteraction[iP] )
        //        {
        ////            cout << "string is hard: " << iP << endl;
        //            stringColor = kRed;
        //        }

        elPointsStrings[iP]->SetLineColor( stringColor );
        elPointsStrings[iP]->SetFillColor( stringColor );
        elPointsStrings[iP]->Draw();
    }
    delete [] elPointsStrings;

    // draw string interaction radius (July 2016)
    //    double visStringIntRad = fVisNucleusRadiusNucleus * fStringInteractionRadius / fNucleusRadius;
    //    TEllipse *elStringInteractionRadius = new TEllipse( 0.12, 0.88
    //                                                        , visStringIntRad, visStringIntRad );
    //    elStringInteractionRadius->SetLineColor( kOrange+1 );
    //    elStringInteractionRadius->SetFillColor( kOrange+1 );
    //    elStringInteractionRadius->SetFillStyle( 3013 );
    //    elStringInteractionRadius->SetLineStyle( 9 );
    //    elStringInteractionRadius->SetLineWidth( 2 );
    //    elStringInteractionRadius->Draw();
    //    delete elStringInteractionRadius;
}





