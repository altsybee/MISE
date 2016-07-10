#ifndef AnalyserForFlowIA_H
#define AnalyserForFlowIA_H

//#include "/Users/macbook/alice/simpleAnalysis/commonTools/QAforWindows.h"
#include "/Users/macbook/alice/simpleAnalysis/analysers/AnalyserBase.h"

#include "TString.h"


class TH1D;
struct MiniEvent;

class AliFlowTrackSimpleCuts;
class AliFlowTrackSimpleCuts;
class AliFlowAnalysisWithQCumulants;
class AliFlowAnalysisWithLeeYangZeros;
class AliFlowAnalysisWithScalarProduct;

class AnalyserForFlowIA : public AnalyserBase {
public:
    // Constructors
    AnalyserForFlowIA();

    // Destructor
    virtual ~AnalyserForFlowIA();

    virtual Bool_t InitDataMembers(); //Is to be called in CreateOutputObjects method

    virtual void StartEvent();  // Open new Event for track by track event import
//    virtual void AddTrack(Double_t Eta, Double_t Phi, Double_t Pt);
    virtual void AddTrack(const SimpleTrack *track);
    virtual void FinishEvent(); // Close opened event
    virtual void Terminate(); // Close opened event

    void SetPID( int pid ) { fPID = pid; } //pid under flow study

protected:
    virtual void UpdateEvent(MiniEvent *event, Double_t Eta , Double_t Phi , Double_t Pt, Int_t Charge = 0, Int_t Pid = -1 );
    virtual void FillAliFlowEvent();
    //Output data
    Int_t fNtracks;
    Int_t fPID;

    Int_t fNeventsQA;

    //centrality class
    Float_t fMinCentrality;    // min bound on centrality percentile
    Float_t fMaxCentrality;    // max bound on centrality percentile
    Float_t fEventCentrality; //event centrality

    TH1D    *fHistNparticlesInTile; //!
    MiniEvent *fMiniEvent;


    AliFlowTrackSimpleCuts *cutsRP;
    AliFlowTrackSimpleCuts *cutsPOI;

    AliFlowAnalysisWithQCumulants* qc_v2;
    AliFlowAnalysisWithQCumulants* qc_v3;
    AliFlowAnalysisWithLeeYangZeros *lyz;

    AliFlowTrackSimpleCuts *cutsPOI_forSP;
    AliFlowAnalysisWithScalarProduct *sp;

    AnalyserForFlowIA(const AnalyserForFlowIA&); // not implemented
    AnalyserForFlowIA& operator=(const AnalyserForFlowIA&); // not implemented

    ClassDef(AnalyserForFlowIA, 1);
};

#endif

