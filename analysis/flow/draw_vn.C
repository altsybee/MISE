void tuneHist1D( TH1 *hist )
{
    hist->GetYaxis()->SetTitleOffset( 1.45 );
    hist->GetYaxis()->SetTitleSize( 0.048 );
    hist->GetYaxis()->SetLabelSize( 0.042 );

    hist->GetXaxis()->SetTitleOffset( 0.95 );
    hist->GetXaxis()->SetTitleSize( 0.048 );
    hist->GetXaxis()->SetLabelSize( 0.042 );
}


void tuneCanvas(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.11);
//    canvas->SetGridy();
}
void tuneLegend( TLegend *leg )
{
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}


void draw_vn()
{
    int kColors[5] = { kBlue, kGreen, kRed, kMagenta, kYellow};
    int kMarkers[3] = { 24, 20, 21 };

    TString strFile[4];
//    strFile[0] = Form("outputCumulants_v2_5k_pions.root");
//    strFile[1] = Form("outputCumulants_v2_5k_kaons.root");
//    strFile[2] = Form("outputCumulants_v2_5k_protons.root");
//    strFile[3] = Form("outputCumulants_v2_5k_allSpecies.root");

//    strFile[0] = Form("outputSP_15k_pions_boltzmanPt.root");
//    strFile[1] = Form("outputSP_15k_kaons_boltzmanPt.root");
//    strFile[2] = Form("outputSP_15k_protons_boltzmanPt.root");

//    strFile[0] = Form("outputSP_10k_pions_boltzmanPt_LOW_T.root");
//    strFile[1] = Form("outputSP_10k_kaons_boltzmanPt_LOW_T.root");
//    strFile[2] = Form("outputSP_10k_protons_boltzmanPt_LOW_T.root");

//    strFile[0] = Form("outputSP_7k_pions_boltzmanPt_T_014.root");
//    strFile[1] = Form("outputSP_7k_kaons_boltzmanPt_T_014.root");
//    strFile[2] = Form("outputSP_7k_protons_boltzmanPt_T_014.root");

//    strFile[0] = Form("outputSP_8k_pions_boltzmanPt_T_014_highRepulsion.root");
//    strFile[1] = Form("outputSP_8k_kaons_boltzmanPt_T_014_highRepulsion.root");
//    strFile[2] = Form("outputSP_8k_protons_boltzmanPt_T_014_highRepulsion.root");

//    strFile[0] = Form("outputSP_15k_pions_boltzmanPt_Schwinger.root");
//    strFile[1] = Form("outputSP_15k_kaons_boltzmanPt_Schwinger.root");
//    strFile[2] = Form("outputSP_15k_protons_boltzmanPt_Schwinger.root");

//    strFile[0] = Form("outputSP_5k_pions_boltzmanPt_Schwinger_WITH_SMEARING.root");
//    strFile[1] = Form("outputSP_5k_kaons_boltzmanPt_Schwinger_WITH_SMEARING.root");
//    strFile[2] = Form("outputSP_5k_protons_boltzmanPt_Schwinger_WITH_SMEARING.root");

    strFile[0] = Form("outputSP_100k_pions_boltzmanPt_Schwinger_WITH_SMEARING.root");
    strFile[1] = Form("outputSP_100k_kaons_boltzmanPt_Schwinger_WITH_SMEARING.root");
    strFile[2] = Form("outputSP_100k_protons_boltzmanPt_Schwinger_WITH_SMEARING.root");



//    strFile[0] = Form("outputSP_10k_pions_boltzmanPt.root");
//    strFile[1] = Form("outputSP_10k_kaons_boltzmanPt.root");
//    strFile[2] = Form("outputSP_10k_protons_boltzmanPt.root");

//    strFile[0] = Form("outputSP_11k_pions_boltzmanPt_tryEP.root");
//    strFile[1] = Form("outputSP_11k_kaons_boltzmanPt_tryEP.root");
//    strFile[2] = Form("outputSP_11k_protons_boltzmanPt_tryEP.root");

//    strFile[0] = Form("outputSP_12k_pions_boltzmanPt_strIntRad1.root");
//    strFile[1] = Form("outputSP_12k_kaons_boltzmanPt_strIntRad1.root");
//    strFile[2] = Form("outputSP_12k_protons_boltzmanPt_strIntRad1.root");

    TCanvas *canv_v = new TCanvas( "canv_v", "canv_v", 100,150,700,600 );
    tuneCanvas(canv_v);

    gStyle->SetOptStat(0);

    TH1D *histSP[10];
    for ( int i = 0; i < 3; i++ )
    {
//        TFile *file = new TFile( Form( "c%d/outputCumulants.root", i) );
        TFile *file = new TFile( strFile[i] );

//        TList *list = (TList*)file->Get("cobjQC");

//        AliFlowCommonHistResults *commonHistRes =
//                        dynamic_cast<AliFlowCommonHistResults*> list->FindObject("AliFlowCommonHistResults2ndOrderQC");
//        TH1D *hist2 = commonHistRes->GetHistDiffFlowPtRP();
//        hist2->SetLineColor( kColors[histId] );
////        if (i!=7)
//            hist2->Draw("same");

//        AliFlowCommonHistResults *commonHistRes =
//                dynamic_cast<AliFlowCommonHistResults*> list->FindObject("AliFlowCommonHistResults4thOrderQC");
//        TH1D *hist4 = commonHistRes->GetHistDiffFlowPtRP();
//        hist4->SetLineColor( kColors[histId] );
//        hist4->SetMarkerStyle( 24 );
//        hist4->SetMarkerSize( 1.1 );
//        hist4->SetMarkerColor( kColors[histId] );
//        hist4->Draw("same");

        TList *list = (TList*)file->Get("cobjSP");

        //from
        // /opt/mygit/MISE/analysis/flow/c6/compareFlowResults.C
        AliFlowCommonHistResults *commonHistRes =
                        dynamic_cast<AliFlowCommonHistResults*> list->FindObject("AliFlowCommonHistResults_SP");
//        TH1D *histSP = commonHistRes->GetHistDiffFlowPtRP();
        histSP[i] = commonHistRes->GetHistDiffFlowPtPOI();


        histSP[i]->SetMarkerColor( kColors[i] );
        histSP[i]->SetLineColor( kColors[i] );
        histSP[i]->SetMarkerStyle( kMarkers[i] );

        tuneHist1D( histSP[i] );
//        if (i!=7)
        histSP[i]->SetTitle(";p_{T}, GeV/#it{c};v_{2} {SP}");
        histSP[i]->DrawCopy( i==0 ? "P" : "P same");

    }


    TLegend *leg = new TLegend(0.6,0.20,0.95,0.55);
    tuneLegend(leg);

    leg->AddEntry( histSP[0], "#pi", "p");
    leg->AddEntry( histSP[1], "K", "p");
    leg->AddEntry( histSP[2], "p", "p");

    leg->Draw();


}


