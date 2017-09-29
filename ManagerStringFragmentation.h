#ifndef ManagerStringFragmentation_H
#define	ManagerStringFragmentation_H

//#include "Event.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"

//#include "ManagerStringFragmentationBase.h"

class TH1D;
class TH2D;
class TCanvas;
//class Event;
//class NucleiCollision;
class StringDescr;
//class AnalyserBase;

class ManagerStringFragmentation {
public:
    ManagerStringFragmentation();
    ManagerStringFragmentation(const ManagerStringFragmentation&);
    virtual ~ManagerStringFragmentation();
    
    // file name with NuclColl data
    void setInputFileNameNucleiCollisions(TString name) { fInputFileName_NucleiCollisions = name; }
    void setInputFileNameStringBoosts(TString name) { fInputFileName_StringBoosts = name; }

    void setOutputDirectoryName(TString strDirName) { fOutputDirName = strDirName; }
    void setOutputFileName(TString name) { fOutputFileName = name; }
    void setDrawHistos(bool flag) { fDrawHistos = flag; }

    void setFillEventTree(bool flag) { fFillEventTree = flag; }

    void setNumberOfCentralityBins(int nCentralityBins) { fNumberOfCentralityBins = nCentralityBins; }
    void setCutMinNumberOfParticles(int n) { fCutMinNumberOfParticles = n; }

    void setWhatToDoWithHardScattering(int what) { whatToDoWithHardScattering = what; }

    //    Event *getEvent() /*const*/
//    {
//        //if (!fEvent)
//        //    fEvent = new Event();
//        return fEvent;
//    }

    void initOutputObjects();
    void applyFragmentationToEvents( StringDescr *strDescr, int nEvents );
    void drawStatHists();

//    void SetFlagGenerateCentralEvent(bool flag) { fFlagGenerateCentralEventByHand = flag; }
//    void SetFlagGenerateSemicentralEvent(bool flag) { fFlagGenerateSemicentralEventByHand = flag; }
//    void SetImpactParameterByHand(float par) { fImpactParameterByHand = par; }

    void cleanup();
    //TList * getProcessorsList() { return &fLRCproc; }
private:
    //TString fStrSpecTitle;

//    NucleiCollision *fPtrNuclStruct;

    bool fPrintInfo;
    bool fDrawHistos;
    bool fFillEventTree;

    int whatToDoWithHardScattering; // 0 - makeTwoJets, 1 - particle pair with random pt from Power law

    TString fInputFileName_NucleiCollisions;
    TString fInputFileName_StringBoosts;

    TString fOutputDirName;
    TString fOutputFileName;
    TString fOutputListName;
//    Event *fEvent;
//    int fNevents;
    Int_t fNumberOfCentralityBins;
    //Distribution *fNsourcesDistr;
    //Distribution *fNparticlesDistr;
    //Distribution *fPtDistr;
    Int_t fCutMinNumberOfParticles;

//    bool fFlagGenerateCentralEventByHand;
//    bool fFlagGenerateSemicentralEventByHand;
//    float fImpactParameterByHand;

    TH1D *fHistSources;
    TH1D *fHistParticlesInSource;
    TH1D *fHistParticlesInSourceIncludingJets;

    TH1D *fHistParticlesInJets;

    TH1D *fHistParticlesInEvent;
    TH1D *fHistParticlesInCutConditionInEvent;
    TH2D *fHistParticlesInCutConditionVsNu;

//    TH1D *fHistParticlesInEventInEta;
//    TH1D *fHistParticlesInCutConditionInEventInEta;

    TH1D* fHistPt;
    TH1D* fHistEta; 
    TH1D* fHistEtaInPtCuts[5];
    TH1D* fHistPhi;
    TH1D* fHistPtAfterCuts;
    TH1D* fHistPtAfterCutsPID[5];
    TH1D* fHistPtBeforeKick;

    TH1D* fHistNeventsInCentralityClasses;

    TH2D* fHist2DEtaPhi;

    TCanvas *fCanv;
    TFile *fOutputFile;
    //TList fLRCproc;       //  AliLRCProcess objects list
    //int fLrcNum;

//    //window sets
//    TList * fWindowSets[100];
//    int fNumberOfWindowSets;
    
//    double fPtCutMin;
//    double fPtCutMax;

    //TRandom3 fRand;
};

#endif	/* ManagerStringFragmentation_H */

