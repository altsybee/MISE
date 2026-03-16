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

void draw_ratio_pT()
{

  // TFile *f_pp = new TFile("../outputs_NucleiCollision_REWAKE_2024_pp/stats_NuclearStructure.root");
  // TFile *f_AA = new TFile("../outputs_NucleiCollision_REWAKE_2024_0_5/stats_NuclearStructure.root");

  // TFile *f_pp = new TFile("../outputs_NucleiCollision_REWAKE_2026_pp_MB_1MEv_NoHardCoreInProton/stats_NuclearStructure.root");
  // TFile *f_AA = new TFile("../outputs_NucleiCollision_REWAKE_2026_PbPb_0_5_1000ev_NoHardCoreInProton/stats_NuclearStructure.root");

  // TFile *f_pp = new TFile("../outputs_NucleiCollision_REWAKE_2026_pp_MB_10MEv_NoHardCoreInProton/stats_NuclearStructure.root");
  // TFile *f_AA = new TFile("../outputs_NucleiCollision_REWAKE_2026_PbPb_0_5_5000ev_NoHardCoreInProton/stats_NuclearStructure.root");
  // TFile *f_pp = new TFile("../outputs_NucleiCollision_REWAKE_2026_pp_MB_25MEv_NoHardCoreInProton_gPDFrangeAdjusted/stats_NuclearStructure.root");
  // TFile *f_AA = new TFile("../outputs_NucleiCollision_REWAKE_2026_PbPb_0_5_10kEv_NoHardCoreInProton_gPDFrangeAdjusted/stats_NuclearStructure.root");
  TFile *f_pp = new TFile("../outputs_NucleiCollision_REWAKE_2026_pp_MB_5MEv_NoHardCoreInProton_alphaS_const_0.2_ShrinkedValenceGauss/stats_NuclearStructure.root");
  TFile *f_AA = new TFile("../outputs_NucleiCollision_REWAKE_2026_PbPb_0_5_5kEv_NoHardCoreInProton_alphaS_const_0.2_ShrinkedValenceGauss/stats_NuclearStructure.root");

  // TFile *f_pp = new TFile("../outputs_NucleiCollision_REWAKE_2026_pp_MB_1MEv_NoHardCoreInProton_changedPDFs/stats_NuclearStructure.root");
  // TFile *f_AA = new TFile("../outputs_NucleiCollision_REWAKE_2026_PbPb_0_5_200ev_NoHardCoreInProton_changedPDFs/stats_NuclearStructure.root");

  // fTheta->SetNpx(1000);
  // fTheta->SetParameter(1, 0.2);           // alpha_s

  TH2 *h2_pp = (TH2 *)f_pp->Get("fHist_partons_pZ_vs_pT_QA");
  TH2 *h2_AA = (TH2 *)f_AA->Get("fHist_partons_pZ_vs_pT_QA");
  TH2 *h2_Ncoll_pp = (TH2 *)f_pp->Get("fHistNcoll");
  TH2 *h2_Ncoll_AA = (TH2 *)f_AA->Get("fHistNcoll");

  TProfile *prof_pp = h2_pp->ProfileX("prof_pp");
  TProfile *prof_AA = h2_AA->ProfileX("prof_AA");

  TCanvas *canv = new TCanvas("canv", "canv", 10, 10, 800, 600);

  prof_pp->DrawCopy();
  prof_AA->SetLineColor(kRed);
  prof_AA->DrawCopy("same");

  TCanvas *canv_ratio = new TCanvas("canv_ratio", "canv_ratio", 40, 40, 800, 600);

  prof_AA->Divide(prof_pp);
  prof_AA->DrawCopy();

  gPad->SetLogx();
  gPad->SetGrid();

  // ###### ratio pT
  TCanvas *canv_partons_pT_QA = new TCanvas("canv_partons_pT_QA", "canv_partons_pT_QA", 90, 90, 800, 600);

  TH1 *h2_pt_pp = (TH1 *)f_pp->Get("fHist_partons_pT_QA");
  TH1 *h2_pt_AA = (TH1 *)f_AA->Get("fHist_partons_pT_QA");

  h2_pt_pp->Rebin(1);
  h2_pt_AA->Rebin(1);

  h2_pt_AA->Divide(h2_pt_pp);
  h2_pt_AA->Scale(h2_Ncoll_pp->GetMean());
  h2_pt_AA->Scale(1. / h2_Ncoll_AA->GetMean());

  h2_pt_AA->DrawCopy();

  gPad->SetGrid();

  // ###### ratio x
  TCanvas *canv_ratio_pdf = new TCanvas("canv_ratio_pdf", "canv_ratio_pdf", 50, 50, 800, 600);

  TH1 *h2_pdf_pp = (TH1 *)f_pp->Get("fHist_partons_x");
  TH1 *h2_pdf_AA = (TH1 *)f_AA->Get("fHist_partons_x");

  h2_pdf_AA->Divide(h2_pdf_pp);
  h2_pdf_AA->DrawCopy();

  gPad->SetGrid();
}