#include "StringBoosting.h"
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

#include "../StringDecayer/ParticleDescr.h"

#include <fstream>
#include <iostream>
using namespace std;

#define PI 3.1415926
//const float kNucleonSize = 1.12;  // fm
//const float kNucleonSize = 0.85; // fm
//const float kNucleonParA = 0.02; // fm

const float Mstring = 1.0; // string energy density 1 GeV (from Abromovsky) //0.6; // stringMassOverRapUnit  //GeV
//const float coeffOverlapEnergy = 0.01; //0.125;//0.005; //GeV

//const
//double kCoeffRadiusForMB = 1.75;  //tuned to be Ok for MB! (for NN)


const float fNucleusRadius = 6.7; // tmp!

const int fMaxStrings = 8000;

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

const int kNstringsWhenClustersIsConsideredLarge = 0;//3;

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



StringBoosting::StringBoosting(/*int seed*/) :
    fRand(0x0)
  , fFlagDataMembersInitialized(false)
  , fFlagWriteEventViewCanvases(false)
  , fFlagComputeStringRepulsion(true)

  , fStringInteractionRadius(0.2)//2) //fm
  , fStringOverlapEnergyDensity(0.1)
  //  , fStringInteractionParA(0.05)
  , fClusterFormationDist(0.4)
  , fHardScatteringProbability(0.03)
  
  , fCanvEventView(0x0)
  , fCanvEventStatistics(0x0)
  , fOutDirName("outputs_StringBoosting")
{
}

void StringBoosting::initDataMembers()
{
    fFlagDataMembersInitialized = true;

    //init 2D array for clusterId->stringsIds
    fStringClusters = new StringCluster*[fMaxStrings]; // CAN'T BE MORE THEN fMaxStrings/2 clusters
    for( int i = 0; i < fMaxStrings; ++i )
        fStringClusters[i] = new StringCluster( fMaxStrings ); // CAN'T BE MORE THEN fMaxStrings/2 strings in one cluster


    
    //    fX1 = new float[fMaxStrings];
    //    fY1 = new float[fMaxStrings];
    //    //    fParton1Participated = new bool[fMaxStrings];
    
    //    fX2 = new float[fMaxStrings];
    //    fY2 = new float[fMaxStrings];
    //    fParton2Participated = new bool[fMaxStrings];
    
    fXstring = new float[fMaxStrings];
    fYstring = new float[fMaxStrings];
    fStringRadiusVectorAngle = new float[fMaxStrings];

    fXstring_FOR_TEST = new float[fMaxStrings];
    fYstring_FOR_TEST = new float[fMaxStrings];

    fXstringInteraction = new float[fMaxStrings];
    fYstringInteraction = new float[fMaxStrings];
    fStringIntDist = new float[fMaxStrings];
    fStringIntAngles = new float[fMaxStrings];
    
    fFlagStringInInteraction = new bool[fMaxStrings];
    fNStringsInCluster = new int[fMaxStrings];
    fClusterIdForString = new int[fMaxStrings];
    fFlagStringIsHardInteraction = new short[fMaxStrings];

    fStringBoostMagn = new float[fMaxStrings];
    fStringBoostAngle = new float[fMaxStrings];
    
//    fNumberOfStrings = 0; // shouldn't set to 0 since now this number is set from outside!
    fNumberOfClusters = 0;
    fNumberOfStringInteractions = 0;
    
    fHistStringInteractions = new TH1D("fHistStringInteractions", "n string intersections", fMaxStrings+1, -0.5, fMaxStrings+0.5 );


    fHistStringInteractionsPhi = new TH1D("fHistStringInteractionsPhi", "string kick #phi;#phi;n strings", 100, 0, 2*TMath::Pi() );
    fHistStringInteractionsDistance = new TH1D("fHistStringInteractionsDistance", "string pairs dist", 100, 0, fNucleusRadius );
    // fHistStringInteractionsMagnitude = new TH1D("fHistStringInteractionsMagnitude", "string interaction magnitude; arb. units;n strings", 500, 0, 100 );
    fHistStringInteractionsMagnitude = new TH1D("fHistStringInteractionsMagnitude", "beta of the strings; #beta;n strings", 220, 0, 1.1 );
    
    fHistForcesInClustersX = new TH1D("fHistForcesInClustersX", "forces acting on string, x; force x, arb units; entries", 200, -10, 10 );
    fHistForcesInClustersY = new TH1D("fHistForcesInClustersY", "forces acting on string, y; force y, arb units; entries", 200, -10, 10 );
    fHistForcesInClustersMagn = new TH1D("fHistForcesInClustersMagn", "forces acting on string, magnitude; force magn, arb units; entries", 200, -10, 10 );

    fHistForcesInClustersAbsX = new TH1D("fHistForcesInClustersAbsX", "force abs(x); force abs(x), arb units; entries", 1000, 0, 100 );
    fHistForcesInClustersAbsY = new TH1D("fHistForcesInClustersAbsY", "force abs(y); force abs(y), arb units; entries", 1000, 0, 100 );
    //    fHistForcesInClustersRatioAbsXY = new TH1D("fHistForcesInClustersRatioAbsXY", "force abs(x)/abs(y); force abs(x)/abs(y), arb units; entries", 50, 0, 5 );


    fHistStringsInClusters = new TH1D("fHistStringsInClusters", "N strings in clusters", fMaxStrings, -0.5, fMaxStrings-0.5 );
    fHist2DStringsInClustersVsB = new TH2D("fHist2DStringsInClustersVsB", "n 'large' clusters vs b;b,fm;n clusters", 50, 0, 2*fNucleusRadius , 41, -0.5, 40.5 );
    

    //    fHist2DStringsSVsImpactS= new TH2D("fHist2DStringsSVsImpactS", "n*s0/S;b/R;n*s0/S", 50, 0, 100 , 100, 0, 2 );
    
    fHistNStringsInLargestCluster = new TH1D("fHistNStringsInLargestCluster", "N strings in largest cluster;n strings", fMaxStrings, -0.5, fMaxStrings-0.5 );
    
    
    fVisNucleusRadiusNucleus = 0.35; //visual size nucleus
    //    rVisParton = fRadiusParton / fNucleusRadius * fVisNucleusRadiusNucleus; //visual nucleon
    
}

StringBoosting::StringBoosting(const StringBoosting& ) {
}

StringBoosting::~StringBoosting() {
    //    delete [] fX1;
    //    delete [] fY1;
    //    delete [] fX2;
    //    delete [] fY2;
}

//######### build event
void StringBoosting::processEvent()
{
    if ( !fRand )
    {
        cout << "StringBoosting: fRand IS NOT INITIALIZED!!! exiting..." << endl;
        return;
    }

    if ( !fFlagDataMembersInitialized )
        initDataMembers();

    for( int i = 0; i < fNumberOfClusters; ++i )
        fStringClusters[i]->reset();

    // TMP! Copy string coords to the test arrays
    for ( int i = 0; i < fNumberOfStrings; i++ )
    {
        fXstring_FOR_TEST[i] = fXstring[i];
        fYstring_FOR_TEST[i] = fYstring[i];
    }


    for ( int i = 0; i < fNumberOfStrings; i++ )
    {
//        cout << "string id " << i << ": fXstring[i]" << fXstring[i] << ", fYstring[i]" << fYstring[i] << endl;
        fStringBoostMagn[i] = 0;
        fStringBoostAngle[i] = 0;
        fFlagStringInInteraction[i] = 0;//false;
        fNStringsInCluster[i] = 0;//false;
        fClusterIdForString[i] = -1;
//        fFlagStringIsHardInteraction[fNumberOfStrings] = 0;//false;

        // istead of "hand-made" probability, try new appoach
        // - with partonic dist for valence quarks! (Sept 2017)
        if(1)
            fFlagStringIsHardInteraction[i] = ( fRand->Uniform() < fHardScatteringProbability ? 1 : -1 );
        else if(0) // use info about string origin from which quarks:
            fFlagStringIsHardInteraction[i] = ( fStringOrigin[i] == 2 && fDistanceBetweenPartonsForString[i]<0.3/*0.2*//*0.1*/ ? 1 : 0 );
        else // USE COOKED BY HAND PIDs*probHardScat! (Nov 2017)
        {
            // prob pids:   (from StringDescr.cpp)
            double prob_pid_pi = 0.748542;
            double prob_pid_K = 0.134895;
            double prob_pid_p = 0.0876396;
            double prob_pid_D0 = 0.0289233; // TMP!!! Nov2017 analysis. Just 1-(pi+K+p)

            //double prob_pid_phi = 0.00913536;
            //double prob_pid_Lambda = 0.00818593;

            // prob hard scatterings for pids (inside a fraction for this pid in whole pool of particles):
            double prob_hard_pi = 1-0.927752;     // PIONS, 0.573339/(0.573339+0.0446483)
            double prob_hard_K  = 1-0.9758306;    // KAONS, 0.0743011/(0.0743011+0.00184029)
            double prob_hard_p  = 1-0.9941083;    // PROTONS, 0.0351141/(0.0351141+0.000208108)
            double prob_hard_D0  = 1-0.985963;    // D0, 1.62524/(1.62524+0.0231375)

            // prob hard scatterings for pids (a fraction from the whole pool of particles):
            double Prob_hard_pi = prob_pid_pi*prob_hard_pi;
            double Prob_hard_K  = prob_pid_K*prob_hard_K;
            double Prob_hard_p  = prob_pid_p*prob_hard_p;
            double Prob_hard_D0  = prob_pid_D0*prob_hard_D0;

            double prob_HARD = Prob_hard_pi + Prob_hard_K + Prob_hard_p + Prob_hard_D0;

//            cout << Prob_hard_pi << " " << Prob_hard_K << " " << Prob_hard_p << endl;

            fFlagStringIsHardInteraction[i] = ( gRandom->Uniform()<prob_HARD ? 1 : -1 );
            if ( fFlagStringIsHardInteraction[i] >= 0 ) // is hard interaction -> now try to set PID:
            {
                double randPid = fRand->Uniform(0,prob_HARD);
                int pid = -1;
                if ( randPid < Prob_hard_pi ) pid = kPid_pion;
                else if ( randPid < Prob_hard_pi + Prob_hard_K ) pid = kPid_kaon;
                else if ( randPid < Prob_hard_pi + Prob_hard_K + Prob_hard_p )
                {
                    pid = kPid_proton;
//                    cout << randPid << endl;
                }
                else if ( randPid < Prob_hard_pi + Prob_hard_K + Prob_hard_p + Prob_hard_D0 )
                    pid = kPid_D0;
                //else if ( randPid < Prob_hard_pi + Prob_hard_K + Prob_hard_p + 0.00913536 ) pid = kPid_phi;
                //else if ( randPid < Prob_hard_pi + Prob_hard_K + Prob_hard_p + 0.00913536 + 0.00818593 ) pid = kPid_Lambda;
                else pid = kPid_Lambda; // tmp! just to assign something
                fFlagStringIsHardInteraction[i] = pid;
//                if (pid==2)
//                    cout << fFlagStringIsHardInteraction[i] << endl;
//                fFlagStringIsHardInteraction[i] = 3;
            }

        }
    }


    if ( fFlagComputeStringRepulsion )
    {
        //find string clusters
        startClusterSearch();

        //string repulsion
        createStringRepulsion();
    }

    //print forces (for debug)
    if(1)for( int i = 0; i < fNumberOfClusters; ++i )
    {
        StringCluster *cluster = fStringClusters[i];
        int nStringsInCluster = cluster->nStrings;
        for ( int iP = 0; iP < nStringsInCluster; iP++ )
        {
            double Fx = cluster->Fx[iP];
            double Fy = cluster->Fy[iP];
            fHistForcesInClustersX->Fill(Fx);
            fHistForcesInClustersY->Fill(Fy);
            fHistForcesInClustersMagn->Fill( sqrt(Fx*Fx+Fy*Fy) );

            fHistForcesInClustersAbsX->Fill(fabs(Fx));
            fHistForcesInClustersAbsY->Fill(fabs(Fy));
            //                if ( fabs(Fy) > 0 ) //tmp?
            //                    fHistForcesInClustersRatioAbsXY->Fill( fabs(Fx)/fabs(Fy) );
        }
        //            cout << ">>>>>>>>>>>>> cluster id " << i << ": " << endl;
        //            fStringClusters[i]->printForces();
    }

    //spec draw for saving
//    if ( fFlagWriteEventViewCanvases && fEventId < 10 )
//    {
//        drawEventStructure();
//        fCanvEventView->SaveAs( Form( "%s/canv_eventView_%d.eps", fOutDirName.Data(), fEventId));
//        fCanvEventView->SaveAs( Form( "%s/canv_eventView_%d.png", fOutDirName.Data(), fEventId));
//    }

    // QA CHECK:
//    cout << "##### QA CHECK FOR Strings X in processEvent()" << endl;
    if(0)for ( int i = 0; i < fNumberOfStrings; i++ )
    {
        if ( fabs(fXstring_FOR_TEST[i] - fXstring[i]) > 0.00000001 )
        cout << "string " << i << ": dX = " << fXstring_FOR_TEST[i] - fXstring[i]
                << ", dY = " << fYstring_FOR_TEST[i] - fYstring[i]
                   << ", fXstring_FOR_TEST is " << fXstring_FOR_TEST[i]
                   << endl;
    }

}



void StringBoosting::startClusterSearch()
{
//    cout << "fNumberOfStrings = " << fNumberOfStrings << endl;

    int nOfLargeClusters = 0;
    int clusterId = 0;
    for ( int i = 0; i < fNumberOfStrings; i++ )
    {
        if ( fClusterIdForString[i]>=0 ) //string is already in some cluster
            continue;
        if ( fFlagStringIsHardInteraction[i] >= 0 ) // this "string" is a hard interaction
            continue;
        int nLinkedStrings = 0;
        findStringClusters( i, nLinkedStrings, clusterId );
        if ( nLinkedStrings > 0 )
        {
            fHistStringsInClusters->Fill( nLinkedStrings );
            //loop again: bind number of strings in cluster to the special array
            for ( int i2 = 0; i2 < fNumberOfStrings; i2++ )
            {
                if ( fClusterIdForString[i2] == clusterId )
                    fNStringsInCluster[i2] = nLinkedStrings;
            }
            //calc forces inside cluster
            createForcesInsideCluster(clusterId);
            clusterId++;
            //            cout << "string " << i << ": " << nLinkedStrings << endl;
        }
        if ( nLinkedStrings > kNstringsWhenClustersIsConsideredLarge )
            nOfLargeClusters++;
    }
    fNumberOfClusters = clusterId;
//    cout << "fNumberOfClusters = " << fNumberOfClusters << endl;
    //GOOD TO KEEP IN FUTURE? !
//    fHist2DStringsInClustersVsB->Fill( fImpactParameter, nOfLargeClusters );
    //    cout << "nOfLargeClusters=" << nOfLargeClusters << endl;
    
    int maxClusterSize = fNStringsInCluster[ maxElement(fNStringsInCluster,fNumberOfStrings) ];
    //    cout << "max cluster size=" << maxClusterSize << endl;
    fHistNStringsInLargestCluster->Fill(maxClusterSize);
}

void StringBoosting::findStringClusters(int iString, int &nLinked, int clusterId)
{
    //    cout << "deep: string id=" << iString << endl;
    const float rCluster2 = fClusterFormationDist*fClusterFormationDist;
    
    float x = fXstring[iString];
    float y = fYstring[iString];
    for ( int i1 = 0; i1 < fNumberOfStrings; i1++ )
    {
        if ( i1 == iString || fClusterIdForString[i1]>=0 )
            continue;
        if ( fFlagStringIsHardInteraction[i1] >= 0 ) // this "string" is a hard interaction
            continue;
        float x1 = fXstring[i1];
        float y1 = fYstring[i1];
        float dx = x1-x;
        float dy = y1-y;
        float dr2 = dx*dx+dy*dy;
        if ( dr2 < rCluster2 )
        {
            if ( nLinked == 0 ) //this is the first pair of linked strings
                fStringClusters[clusterId]->addStringId(iString);

            nLinked++;
            //bind current stringId into the 2D array
            fStringClusters[clusterId]->addStringId(i1);

            //            cout << "  found neighbour with id=" << i1 << endl;
            //            fNStringsInCluster[iString] = nLinked;//++;
            //            fNStringsInCluster[i1] = nLinked;//++;
            fClusterIdForString[iString] = clusterId;
            fClusterIdForString[i1] = clusterId;

            findStringClusters(i1, nLinked, clusterId);
        }
    }
}


void StringBoosting::createForcesInsideCluster(int clusterId)
{
    float maxStringOverlapArea = TMath::Pi()*fStringInteractionRadius*fStringInteractionRadius; //circlesIntersectionArea(0, fStringInteractionRadius, fStringInteractionRadius);

    // !!!!! playing with params!
    // take DOUBLED parameter (distance) for cuts:
    const float fInteractionDist = 2*fStringInteractionRadius; //3;//1.8; //4.*fClusterFormationDist;//3*fNucleusRadius; //10.*fClusterFormationDist;
    const float rInteraction2 = fInteractionDist*fInteractionDist;

    StringCluster *thisCluster = fStringClusters[clusterId];
    for ( int i = 0; i < thisCluster->nStrings-1; i++ )
    {
        int strId = thisCluster->linkedStringIds[i];
        float x = fXstring[strId];
        float y = fYstring[strId];
        for ( int i2 = i+1; i2 < thisCluster->nStrings; i2++ )
        {
            int strId2 = thisCluster->linkedStringIds[i2];
            //            if ( strId2 == strId ) //it's the same string
            //                continue;
            // ##### create forces from nearby strings
            //check x
            float x2 = fXstring[strId2];
            if ( fabs(x-x2) > fInteractionDist )
                continue;
            //check y
            float y2 = fYstring[strId2];
            if ( fabs(y-y2) > fInteractionDist )
                continue;
            //check r2
            float dx = x2-x;
            float dy = y2-y;
            float dr2 = dx*dx+dy*dy;
            if ( dr2 < rInteraction2 ) //strings are in interaction
            {
                float dr = sqrt(dr2);
//                cout << "dr = " << dr << endl;
                // 29.11.2014 - old version: "non-relativistic" forces
                if (0)
                {
                    //                    float forceCoeff = 0.05 * funcStringInteractionDist->Eval( dr );
                    //                    //                if ( dr > fClusterFormationDist )
                    //                    //                    forceCoeff = 1/dr;
                    //                    thisCluster->Fx[i] += -dx/dr * forceCoeff;
                    //                    thisCluster->Fy[i] += -dy/dr * forceCoeff;
                    //                    thisCluster->Fx[i2] += dx/dr * forceCoeff;
                    //                    thisCluster->Fy[i2] += dy/dr * forceCoeff;
                }
                float stringOverlapArea = circlesIntersectionArea(dr, fStringInteractionRadius, fStringInteractionRadius);
                float U = stringOverlapArea/maxStringOverlapArea * fStringOverlapEnergyDensity;
                //                cout << circlesIntersectionArea(0, fStringInteractionRadius, fStringInteractionRadius)/maxStringOverlapArea << endl;
                float momAbs = sqrt( (Mstring+U)*(Mstring+U) - Mstring*Mstring ); //acquired momentum, see Abramovsky paper
                //                cout << "momString = " << momAbs << endl;
                thisCluster->Fx[i] += -dx/dr * momAbs;
                thisCluster->Fy[i] += -dy/dr * momAbs;
                thisCluster->Fx[i2] += dx/dr * momAbs;
                thisCluster->Fy[i2] += dy/dr * momAbs;

//                cout << "thisCluster->Fx[i] = " << thisCluster->Fx[i] << endl;
            }
        }
    }
}


void StringBoosting::createStringRepulsion()
{
    for ( int iCluster = 0; iCluster < fNumberOfClusters; iCluster++ )
    {
//        cout << "test... iCluster = " << iCluster << endl;
        StringCluster *cluster = fStringClusters[iCluster];
        for ( int i = 0; i < cluster->nStrings; i++ )
        {
//            cout << "test... iString = " << i << endl;
            int strId = cluster->linkedStringIds[i];

            //            float Fx = cluster->Fx[i];
            //            float Fy = cluster->Fy[i];
            //            float absF = sqrt(Fx*Fx+Fy*Fy);
            float Px = cluster->Fx[i];
            float Py = cluster->Fy[i];
//            cout << "Px = " << Px << ", Py = " << Py << endl;
            float absMom = sqrt(Px*Px+Py*Py);
            float beta = absMom/sqrt( absMom*absMom + Mstring*Mstring );
            //            cout << "absMom = " << absMom << endl;
//                        cout << "beta = " << beta << endl;

            double interactionAngle = asin( Py/absMom );
            if ( Px < 0 )
                interactionAngle = TMath::Pi()-interactionAngle;
            FixAngleInTwoPi(interactionAngle);
            FixAngleInTwoPi(interactionAngle);

            fStringBoostMagn[strId] = beta;//absF;
            fStringBoostAngle[strId] = interactionAngle;
            fHistStringInteractionsMagnitude->Fill(beta); //absF);
            fHistStringInteractionsPhi->Fill(interactionAngle);

        }
    }
}


void StringBoosting::createStringRepulsionOld()
{
    //string interactions - find the closest string pairs
    fStringInteractionsFinder->FindMinDistancesWithinArray(fXstring,fYstring,fNumberOfStrings);

    //define "string interactions" as pairs of strings
    //    cout << "fNumberOfStrings=" << fNumberOfStrings << endl;
    int distBetweenStringsArraySize = 0;
    DistanceEntry* arrDistBetweenStrings = fStringInteractionsFinder->GetMinDistanceArray( distBetweenStringsArraySize );
    fNumberOfStringInteractions = 0;
    for ( int iP = 0; iP < distBetweenStringsArraySize; iP++ )
    {
        if ( !arrDistBetweenStrings[iP].inInteraction )
            continue;
        int id1 = arrDistBetweenStrings[iP].x;
        int id2 = arrDistBetweenStrings[iP].y;

        double x1 = fXstring[id1];
        double y1 = fYstring[id1];
        double x2 = fXstring[id2];
        double y2 = fYstring[id2];

        fXstringInteraction[fNumberOfStringInteractions] = (x1+x2)/2;
        fYstringInteraction[fNumberOfStringInteractions] = (y1+y2)/2;
        double dx=x1-x2;
        double dy=y1-y2;
        double dist = sqrt(dx*dx+dy*dy);
        fStringIntDist[fNumberOfStringInteractions] = dist;
        double interactionAngle = asin(dy/dist);
        if (dx<0)
            interactionAngle = TMath::Pi()-interactionAngle;
        FixAngleInTwoPi(interactionAngle);
        FixAngleInTwoPi(interactionAngle);
        fStringIntAngles[fNumberOfStringInteractions] = interactionAngle;

        double interactionAngle2 = interactionAngle + TMath::Pi();
        FixAngleInTwoPi(interactionAngle2);

        //for each string write boost information: the magnitude and the angle
        if(0) //OLD BOOSTS, new are in createStringRepulsion
        {
            const double coeffInteraction = 1.;
            const double maxMagnitude = 0.1;
            double interacitonMagnitude = coeffInteraction*( dist>0 ? 1./dist : maxMagnitude );  // dist;
            fStringBoostMagn[id1] = interacitonMagnitude;
            fStringBoostAngle[id1] = interactionAngle;
            fStringBoostMagn[id2] = interacitonMagnitude;
            fStringBoostAngle[id2] = interactionAngle2;
            fHistStringInteractionsDistance->Fill(dist);
            fHistStringInteractionsMagnitude->Fill(interacitonMagnitude);
            //        fHistStringInteractionsPhi->Fill(interactionAngle);
            fHistStringInteractionsPhi->Fill(interactionAngle2);
        }

        fFlagStringInInteraction[id1] = true;
        fFlagStringInInteraction[id2] = true;

        fNumberOfStringInteractions++;
    }
    fHistStringInteractions->Fill( fNumberOfStringInteractions );
}


//######### drawings
void StringBoosting::drawStatisticHists()
{
    cout << "test1" << endl;

    //save to file
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

        //        fHistImpParNucleonsFromAB->Write();

        //pad
        fCanvEventStatistics->cd(padId++);
        fHistForcesInClustersX->SetLineColor(kGreen);
        fHistForcesInClustersY->SetLineColor(kRed);
        fHistForcesInClustersMagn->DrawCopy();
        fHistForcesInClustersX->DrawCopy("same");
        fHistForcesInClustersY->DrawCopy("same");

        fHistForcesInClustersMagn->Write();
        fHistForcesInClustersX->Write();
        fHistForcesInClustersY->Write();

        cout << "test2" << endl;
        //pad 5
        fCanvEventStatistics->cd(padId++);
        fHistStringInteractionsPhi->DrawCopy();
        fHistStringInteractionsPhi->Write();

//        fHistStringRadiusVectorPhi->SetLineColor(kRed);
//        fHistStringRadiusVectorPhi->DrawCopy( "same" );

        cout << "test3" << endl;

        //pad 6
        //        fCanvEventStatistics->cd(6);
        //        fHistStringInteractionsDistance->DrawCopy();
        //pad 7
        fCanvEventStatistics->cd(padId++);//->SetLogy();
        fHistStringInteractionsMagnitude->DrawCopy();
        fHistStringInteractionsMagnitude->Write();
        cout << ">> mean string boost (<beta>) = " << fHistStringInteractionsMagnitude->GetMean() << endl;

        cout << "test4" << endl;

        fileNuclStructStats->Close();

        fCanvEventStatistics->SaveAs( Form( "%s/canvNuclearStructureStats.root", fOutDirName.Data()) );
        fCanvEventStatistics->SaveAs( Form( "%s/canvNuclearStructureStats.eps", fOutDirName.Data()));
        fCanvEventStatistics->SaveAs( Form( "%s/canvNuclearStructureStats.png", fOutDirName.Data()));
        cout << "test5" << endl;
    }

}

void StringBoosting::finalActions()
{    
    cout << "########## Begin StringBoosting::finalActions()" << endl;
    
    // print Fx, Fy results
    float meanFx = fHistForcesInClustersAbsX->GetMean();
    float meanFy = fHistForcesInClustersAbsY->GetMean();
    double errFx = fHistForcesInClustersAbsX->GetMeanError();
    double errFy = fHistForcesInClustersAbsY->GetMeanError();
    double errNumDenom = errFx + errFy;
    double errTotal = sqrt(errNumDenom*errNumDenom + errNumDenom*errNumDenom);

    float forceEccentricity = (meanFx-meanFy) / (meanFx+meanFy);
    if (1)
    {
        cout << "forceEccentricity = " << forceEccentricity << endl;
        ofstream fout( Form( "%s/tmpTextOutput4.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout //<< fImpactParameterByHand
             << " " << forceEccentricity
             << " " << errTotal
             << endl;
        fout.close();
    }

    if (1)
    {
        ofstream fout( Form( "%s/tmpTextOutput5.txt", fOutDirName.Data() ), ios::out | ios::binary);
        fout //<< fImpactParameterByHand
             << " " << fHistStringInteractionsMagnitude->GetMean()
             << " " << fHistStringInteractionsMagnitude->GetMeanError()
             << endl;
        fout.close();
    }

    //print mean string boost
    cout << "(mean string boost, err) = " << //fImpactParameterByHand << " " <<
            fHistStringInteractionsMagnitude->GetMean()
         << " " << fHistStringInteractionsMagnitude->GetMeanError() << endl;

    cout << "########## End StringBoosting::finalActions()" << endl;
}


void StringBoosting::drawEventStructure()
{
    if (!fCanvEventView)
    {
        fCanvEventView = new TCanvas("Event Canvas","Event View",600,70,700,700);
        //fCanvEventView->Divide(2,2);
    }
//    fCanvEventView->Clear();

//    drawPartons();
    drawStrings();
//    if ( fFlagHaveMBcollision )
    {
        //        drawStringInteractions();
        drawForcesInsideClusters();
        //            drawStringBoosts();
    }
    //    fCanvEventView->SaveAs( Form( "eventViews/eventView_%d.eps", fEventId) );
//        fCanvEventView->SaveAs( Form( "eventView_%d.eps", 0) );
}

void StringBoosting::drawStringInteractions()
{
    //visualize "string interactions"
    TLine **linesStringInteracions = new TLine* [fNumberOfStringInteractions];
    for ( int iP = 0; iP < fNumberOfStringInteractions; iP++ )
    {
        double x1 = fXstringInteraction[iP] - fStringIntDist[iP]/2 * cos(fStringIntAngles[iP]);
        double y1 = fYstringInteraction[iP] - fStringIntDist[iP]/2 * sin(fStringIntAngles[iP]);
        double x2 = fXstringInteraction[iP] + fStringIntDist[iP]/2 * cos(fStringIntAngles[iP]);
        double y2 = fYstringInteraction[iP] + fStringIntDist[iP]/2 * sin(fStringIntAngles[iP]);
        //lines
        linesStringInteracions[iP] = new TLine( 0.5+fVisNucleusRadiusNucleus * x1/fNucleusRadius, 0.5+fVisNucleusRadiusNucleus * y1/fNucleusRadius,
                                                0.5+fVisNucleusRadiusNucleus * x2/fNucleusRadius, 0.5+fVisNucleusRadiusNucleus * y2/fNucleusRadius );
        linesStringInteracions[iP]->SetLineColor( kGreen+2 );
        linesStringInteracions[iP]->SetLineWidth( 2 );
        linesStringInteracions[iP]->Draw();
    }
    delete [] linesStringInteracions;
}

void StringBoosting::drawForcesInsideClusters()
{
    //find min and max forces
    float min_dr = 0;
    float max_dr = 0;
    for ( int i = 0; i < fNumberOfClusters; i++ )
    {
        StringCluster *cluster = fStringClusters[i];
        int nStringsInCluster = cluster->nStrings;
        //visualize "string interactions"
        for ( int iP = 0; iP < nStringsInCluster; iP++ )
        {
            float Fx = cluster->Fx[iP];
            float Fy = cluster->Fy[iP];
            float dr = sqrt(Fx*Fx+Fy*Fy);

            if( dr > max_dr )
                max_dr = dr;
            if( dr < min_dr )
                min_dr = dr;
        }
    }
    cout << "min,max dr: " << min_dr << " " << max_dr << endl;

    //visualize "string interactions"
    const float kVizForce = 1;
    for ( int i = 0; i < fNumberOfClusters; i++ )
    {
        StringCluster *cluster = fStringClusters[i];
        int nStringsInCluster = cluster->nStrings;
        TArrow **arrowsStringBoosts = new TArrow* [nStringsInCluster];
        cout << "nStringsInCluster=" << nStringsInCluster << endl;
        for ( int iP = 0; iP < nStringsInCluster; iP++ )
        {
            float xString = fXstring[ cluster->linkedStringIds[iP] ];
            float yString = fYstring[ cluster->linkedStringIds[iP] ];
            float Fx = cluster->Fx[iP];
            float Fy = cluster->Fy[iP];

            float xStringBoost = xString + kVizForce*Fx/max_dr;
            float yStringBoost = yString + kVizForce*Fy/max_dr;
            arrowsStringBoosts[iP] = new TArrow( 0.5+fVisNucleusRadiusNucleus * xString/fNucleusRadius, 0.5+fVisNucleusRadiusNucleus * yString/fNucleusRadius,
                                                 0.5+fVisNucleusRadiusNucleus * xStringBoost/fNucleusRadius, 0.5+fVisNucleusRadiusNucleus * yStringBoost/fNucleusRadius,
                                                 0.025);
            arrowsStringBoosts[iP]->SetLineColor( kGreen + (int)fRand->Uniform(-9,4) );//kMagenta+2 );
            arrowsStringBoosts[iP]->SetLineWidth( 2 );
            if( sqrt(Fx*Fx+Fy*Fy) > /*0.33*/0.55*(max_dr-min_dr) ) //draw only long arrows
                arrowsStringBoosts[iP]->Draw();
        }
        delete [] arrowsStringBoosts;
    }
}


void StringBoosting::drawStringBoosts()
{
    //arrows
    TArrow **arrowsStringBoosts = new TArrow* [fNumberOfStrings];
    for ( int iP = 0; iP < fNumberOfStrings; iP++ )
    {
        if ( !fFlagStringInInteraction[iP] )
            continue;
        double xString = fXstring[iP];
        double yString = fYstring[iP];
        
        double xStringBoost = xString + 0.01*fStringBoostMagn[iP] * cos(fStringBoostAngle[iP]);
        double yStringBoost = yString + 0.01*fStringBoostMagn[iP] * sin(fStringBoostAngle[iP]);
        arrowsStringBoosts[iP] = new TArrow( 0.5+fVisNucleusRadiusNucleus * xString/fNucleusRadius, 0.5+fVisNucleusRadiusNucleus * yString/fNucleusRadius,
                                             0.5+fVisNucleusRadiusNucleus * xStringBoost/fNucleusRadius, 0.5+fVisNucleusRadiusNucleus * yStringBoost/fNucleusRadius,
                                             0.025);
        arrowsStringBoosts[iP]->SetLineColor( kMagenta+2 );
        arrowsStringBoosts[iP]->SetLineWidth( 2 );
        arrowsStringBoosts[iP]->Draw();
    }
    delete [] arrowsStringBoosts;
    
}


void StringBoosting::drawStrings()
{
    //visualize "strings"
    TEllipse **elPointsStrings = new TEllipse* [fNumberOfStrings];
    float rVisString = 0.008;//rVisParton * 0.2;
    //    int nUsedClusterColors = 0;
//    cout << "fNumberOfStrings in drawStrings() =" << fNumberOfStrings << endl;

    // QA CHECK:
//   cout << "##### QA CHECK FOR Strings X in drawStrings()" << endl;
   if(0)for ( int i = 0; i < fNumberOfStrings; i++ )
    {
        if ( fabs(fXstring_FOR_TEST[i] - fXstring[i]) > 0.00000001 )
        cout << "string " << i << ": dX = " << fXstring_FOR_TEST[i] - fXstring[i]
                << ", dY = " << fYstring_FOR_TEST[i] - fYstring[i]
                   << ", fXstring_FOR_TEST is " << fXstring_FOR_TEST[i]
                   << endl;
    }


    for ( int iP = 0; iP < fNumberOfStrings; iP++ )
    {
        elPointsStrings[iP] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * fXstring[iP] / fNucleusRadius, 0.5 + fVisNucleusRadiusNucleus * fYstring[iP] / fNucleusRadius
                                            , rVisString, rVisString );
        //        elPointsStrings[iP]->SetLineColor( fNStringsInCluster[iP] ? kGreen+fNStringsInCluster[iP] : kBlue );
        //        elPointsStrings[iP]->SetFillColor( fNStringsInCluster[iP] ? kGreen+fNStringsInCluster[iP] : kBlue );
        //        int stringIdInCluster = fNStringsInCluster[iP];
        int clusterId = fClusterIdForString[iP];

        //        cout << "draw string: " << stringIdInCluster << endl;
        //        int clusterColor = ( stringIdInCluster < 4 ? kGreen+stringIdInCluster : kRed+1 );

        int clusterColor = kOrange;
        if ( fNStringsInCluster[iP] > kNstringsWhenClustersIsConsideredLarge )
        {
            clusterColor = kBlue;//+2;//( clusterId < nClusterColors ? kClusterColors[clusterId] : kAzure+10 );
        }

        int colorWhenNotInCluster = kBlue;
        //mark hard interactions
        if ( fFlagStringIsHardInteraction[iP]>=0 )
        {
//            cout << "string is hard: " << iP << endl;
            colorWhenNotInCluster = kRed;
        }

        elPointsStrings[iP]->SetLineColor( clusterId>=0 ? clusterColor : colorWhenNotInCluster );
        elPointsStrings[iP]->SetFillColor( clusterId>=0 ? clusterColor : colorWhenNotInCluster );
//        if( fabs(fXstring[iP]) < 0.000001 )
        {
//            cout << "something wrong with the string: " << iP << endl;
            elPointsStrings[iP]->Draw();
        }
//        cout << "TEST" << endl;
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

