#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TMath.h"

#include "/Users/igor/ALICE_Run2_analyses/utils.C"

void Theta_with_AlphaS_via_Q()
{

  const int colors[] = {
      kRed,
      kRed + 2,
      kMagenta,
      kMagenta + 2,
      kOrange + 7,
      kOrange + 9,
      kPink + 1,
      kPink + 6,
      kViolet,
      kViolet - 1,
      //    kViolet-2, kViolet-3,
      kViolet - 3,
      kViolet - 5,
      //    kViolet-6, kViolet-7,
      kViolet - 7,
      kViolet - 9,
  };

  const double kFactorQcdVsQedByHand = 3.;
  const double hc = 0.197;                     // GeV*fm
  const double H = kFactorQcdVsQedByHand * hc; // GeV*fm

  const double Lambda0 = 0.3; // GeV - to avoid too high alphaS for low Q (important for low Ecm!)
  const double A = 0.118 / 0.236 * 12 * TMath::Pi() / 23;

  // prepare for loops and drawing
  TCanvas *canv_Theta = new TCanvas("canv_Theta", "canv_Theta", 90, 90, 800, 600);

  const int nBinsImp = 10;
  TGraphErrors *grAtB[nBinsImp];

  const int nPzSlices = 5;
  TGraphErrors *grThetaVsImpAtPz[nPzSlices];
  TGraphErrors *grPtVsImpAtPz[nPzSlices];
  TGraphErrors *grPtVsImpAtPzOldCalc[nPzSlices];

  // loop over b and pz, get Theta
  for (int i = 0; i < nPzSlices; i++)
  {
    grThetaVsImpAtPz[i] = new TGraphErrors;
    grPtVsImpAtPz[i] = new TGraphErrors;
    grPtVsImpAtPzOldCalc[i] = new TGraphErrors;
  }

  for (int i = 0; i < nBinsImp; i++)
  {
    float b = 0.02 + 0.02 * i; // fm

    grAtB[i] = new TGraphErrors;

    const int nPzBins = 200;
    for (int j = 0; j < nPzBins; j++)
    {
      // float pz = 0.5 + 0.2 * j; // GeV
      // float pz = 1.0 + 0.5 * j; // GeV
      float pz = 1.0 + 0.5 * j; // GeV

      //
      double B = H * A / (Lambda0 * b);
      double W = ROOT::Math::lambert_W0(B);

      double Theta = (Lambda0 / pz) * B / W;

      // old calc
      double QforAlphaS = 2 * pz;
      double alphaS = 0.2;//A / log(QforAlphaS * QforAlphaS / Lambda0 / Lambda0);
      double ThetaOld = 2 * atan(alphaS * H / (2 * b * pz)); // alphaS * 0.197 GeV*fm / b[fm] / p[GeV]

      // cout << "W(" << x << ") = " << w << endl;
      // cout << "Theta = " << Theta << endl;
      // cout << "b=" << b << ", pz=" << pz << ", Theta = " << Theta << endl;

      grAtB[i]->SetPoint(j, pz, Theta);

      if (j == 0)
      {
        grThetaVsImpAtPz[0]->SetPoint(grThetaVsImpAtPz[0]->GetN(), b, Theta);
        grPtVsImpAtPz[0]->SetPoint(grPtVsImpAtPz[0]->GetN(), b, sin(Theta) * pz);
        grPtVsImpAtPzOldCalc[0]->SetPoint(grPtVsImpAtPzOldCalc[0]->GetN(), b, sin(ThetaOld) * pz);
      }
      if (j == 2)
      {
        grThetaVsImpAtPz[1]->SetPoint(grThetaVsImpAtPz[1]->GetN(), b, Theta);
        grPtVsImpAtPz[1]->SetPoint(grPtVsImpAtPz[1]->GetN(), b, sin(Theta) * pz);
        grPtVsImpAtPzOldCalc[1]->SetPoint(grPtVsImpAtPzOldCalc[1]->GetN(), b, sin(ThetaOld) * pz);
      }
      if (j == nPzBins / 2)
      {
        grThetaVsImpAtPz[2]->SetPoint(grThetaVsImpAtPz[2]->GetN(), b, Theta);
        grPtVsImpAtPz[2]->SetPoint(grPtVsImpAtPz[2]->GetN(), b, sin(Theta) * pz);
        grPtVsImpAtPzOldCalc[2]->SetPoint(grPtVsImpAtPzOldCalc[2]->GetN(), b, sin(ThetaOld) * pz);
      }
      if (j == nPzBins - 1)
      {
        grThetaVsImpAtPz[3]->SetPoint(grThetaVsImpAtPz[3]->GetN(), b, Theta);
        grPtVsImpAtPz[3]->SetPoint(grPtVsImpAtPz[3]->GetN(), b, sin(Theta) * pz);
        grPtVsImpAtPzOldCalc[3]->SetPoint(grPtVsImpAtPzOldCalc[3]->GetN(), b, sin(ThetaOld) * pz);
      }
    }

    grAtB[i]->SetTitle(";p_{z} (GeV);#theta (rad)");
    drawGraph(grAtB[i], 1, colors[i], i == 0 ? "APLz" : "PLz", 1, 1);
  }

  gPad->SetGrid();

  //
  TCanvas *canv_ThetaVsImpAtPz = new TCanvas("canv_ThetaVsImpAtPz", "canv_ThetaVsImpAtPz", 20, 20, 800, 600);
  grThetaVsImpAtPz[0]->SetTitle(";b (fm);#theta (rad)");

  drawGraph(grThetaVsImpAtPz[0], 20, colors[0], "APLz");
  drawGraph(grThetaVsImpAtPz[1], 20, colors[1], "PLz");
  drawGraph(grThetaVsImpAtPz[2], 20, colors[2], "PLz");
  drawGraph(grThetaVsImpAtPz[3], 20, colors[3], "PLz");
  gPad->SetGrid();

  //
  TCanvas *canv_PtVsImpAtPz = new TCanvas("canv_PtVsImpAtPz", "canv_PtVsImpAtPz", 20, 20, 800, 600);
  grPtVsImpAtPz[0]->SetTitle(";b (fm);p_{T} (GeV/c)");

  drawGraph(grPtVsImpAtPz[0], 20, colors[0], "APLz");
  drawGraph(grPtVsImpAtPz[1], 20, colors[1], "PLz");
  drawGraph(grPtVsImpAtPz[2], 20, colors[2], "PLz");
  drawGraph(grPtVsImpAtPz[3], 20, colors[3], "PLz");

  grPtVsImpAtPzOldCalc[0]->SetLineStyle(2);
  grPtVsImpAtPzOldCalc[1]->SetLineStyle(2);
  grPtVsImpAtPzOldCalc[2]->SetLineStyle(2);
  grPtVsImpAtPzOldCalc[3]->SetLineStyle(2);

  drawGraph(grPtVsImpAtPzOldCalc[0], 24, colors[0], "PLz");
  drawGraph(grPtVsImpAtPzOldCalc[1], 24, colors[1], "PLz");
  drawGraph(grPtVsImpAtPzOldCalc[2], 24, colors[2], "PLz");
  drawGraph(grPtVsImpAtPzOldCalc[3], 24, colors[3], "PLz");

  gPad->SetGrid();
}