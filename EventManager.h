#ifndef EVENTMANAGER_H
#define	EVENTMANAGER_H

//#include "Event.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"

//#include "EventManagerBase.h"

class TH1D;
class TH2D;
class TCanvas;
//class Event;
class NucleusStructure;
class StringDescr;
class AnalyserBase;

class EventManager {
public:
    EventManager();
    EventManager(const EventManager&);
    virtual ~EventManager();
    
    void setOutputFileName(TString name) {    fOutputFileName = name; }
    void setDrawHistos(bool flag) { fDrawHistos = flag; }

    void setFillEventTree(bool flag) { fFillEventTree = flag; }

    //    Event *getEvent() /*const*/
//    {
//        //if (!fEvent)
//        //    fEvent = new Event();
//        return fEvent;
//    }

    void generateEvents(NucleusStructure *d, StringDescr *strDescr, int nEvents );
    void drawStatHists();

//    void SetFlagGenerateCentralEvent(bool flag) { fFlagGenerateCentralEventByHand = flag; }
//    void SetFlagGenerateSemicentralEvent(bool flag) { fFlagGenerateSemicentralEventByHand = flag; }
//    void SetImpactParameterByHand(float par) { fImpactParameterByHand = par; }

    void cleanup();
    //TList * getProcessorsList() { return &fLRCproc; }
private:
    //TString fStrSpecTitle;

    NucleusStructure *fPtrNuclStruct;

    bool fPrintInfo;
    bool fDrawHistos;
    bool fFillEventTree;
    TString fOutputFileName;
    TString fOutputListName;
//    Event *fEvent;
    int fNevents;
    //Distribution *fNsourcesDistr;
    //Distribution *fNparticlesDistr;
    //Distribution *fPtDistr;

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

#endif	/* EVENTMANAGER_H */

