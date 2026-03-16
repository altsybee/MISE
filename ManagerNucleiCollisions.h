#ifndef ManagerNucleiCollisions_H
#define ManagerNucleiCollisions_H

#include "TList.h"
#include "TFile.h"
#include "TString.h"

class TH1D;
class TH2D;
class TCanvas;
class NucleiCollision;

class ManagerNucleiCollisions
{
public:
    ManagerNucleiCollisions();
    // ManagerNucleiCollisions(const ManagerNucleiCollisions &);
    ~ManagerNucleiCollisions();

    void setOutputDirectoryName(TString strDirName) { fOutputDirName = strDirName; }
    //    void setOutputFileName(TString name) { fOutputFileName = name; }
    void setDrawHistos(bool flag) { fDrawHistos = flag; }
    void setFillTree(bool flag) { fFillTree = flag; }
    // void initOutputObjects();
    void generateEvents(NucleiCollision *d, int nEvents);

private:
    //    NucleiCollision *fPtrNuclStruct;
    TString fOutputDirName;
    bool fPrintInfo;
    bool fDrawHistos;

    bool fFillTree;
};

#endif /* ManagerNucleiCollisions_H */
